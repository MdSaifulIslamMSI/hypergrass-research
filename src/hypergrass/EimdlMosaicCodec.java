package hypergrass;

import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

/**
 * Evolutionary-Information Minimum Description Length (EIMDL) Mosaic Parser, Pricer, Encoder, and
 * Decoder implementing Laws I–IV of {@code docs/COMPRESSION-LAW.md}:
 *
 * <ul>
 *   <li><b>Law I:</b> Klein four-group ({@code V_4}) strand-invariant Kimura + CpG trinucleotide
 *       substitution coding ({@code delta = tgt ^ ref in {1, 2, 3}}).</li>
 *   <li><b>Law II:</b> Bijective codon-wobble phase ({@code (gap + 1) mod 3}) and coalescent
 *       divergence-regime gap coding with exact boundary clamping.</li>
 *   <li><b>Law III:</b> Sticky twin/parent donor priority, per-donor collinear coordinate memory,
 *       MDL-optimal tail trimming, and fast SWAR syntenic micro-indel splicing.</li>
 *   <li><b>Law IV:</b> 5-regime competition across panel homology, bidirectional ({@code +/-})
 *       approximate self-homology, overlapping tandem periodicity ({@code P in [1..64]}), and
 *       extracted novelty DNA segments.</li>
 * </ul>
 */
public final class EimdlMosaicCodec {

    public static final int MIN_MATCH_LEN = 20;
    public static final int MIN_SPLICE_ANCHOR = 12;
    public static final int MIN_TANDEM_LEN = 18;
    public static final int MAX_TANDEM_PERIOD = 64;

    private static final int PROB_BITS = 12;
    private static final int PROB_SCALE = 1 << PROB_BITS; // 4096
    private static final int PROB_HALF = 1 << (PROB_BITS - 1); // 2048
    private static final int TOP = 1 << 24;

    private static final int MAX_BUCKET_PROBES = 48;
    private static final int MAX_SPLICE_CHAIN = 1024;

    private EimdlMosaicCodec() {}

    /** Ablation &amp; mode configuration serialized into a 1-byte mosaic header. */
    public record MosaicConfig(
            boolean allowApprox,
            boolean allowIndelSplice,
            boolean allowKimuraV4,
            boolean allowWobbleMod3,
            boolean allowSelfMatch,
            boolean allowTandem) {

        public static final MosaicConfig DEFAULT =
                new MosaicConfig(true, true, true, true, true, true);

        int toHeaderByte(boolean hasPanel) {
            int f = 0;
            if (hasPanel) f |= 1;
            if (allowApprox) f |= 2;
            if (allowIndelSplice) f |= 4;
            if (allowKimuraV4) f |= 8;
            if (allowWobbleMod3) f |= 16;
            if (allowSelfMatch) f |= 32;
            if (allowTandem) f |= 64;
            return f;
        }

        static MosaicConfig fromHeaderByte(int f) {
            return new MosaicConfig(
                    (f & 2) != 0,
                    (f & 4) != 0,
                    (f & 8) != 0,
                    (f & 16) != 0,
                    (f & 32) != 0,
                    (f & 64) != 0);
        }
    }

    /** Summary telemetry emitted by mosaic encoding. */
    public record MosaicStats(
            int targetBases,
            int panelMatchBases,
            int selfMatchBases,
            int tandemMatchBases,
            int noveltyBases,
            int panelTokens,
            int selfTokens,
            int tandemTokens,
            int literalRuns,
            int indelSplices,
            int substitutions) {}

    /** Output of encoding one sequence's mosaic structure + extracted novelty bases. */
    public record EncodedMosaic(
            byte[] mosaicBytes, byte[] noveltyBases, int noveltyLen, MosaicStats stats) {}

    /** Evolutionary telemetry decoded directly from Mosaic tokens. */
    public record MosaicInspectStats(
            long targetLen,
            long noveltyLen,
            long panelMatchBases,
            long selfMatchBases,
            long tandemMatchBases,
            long literalBases,
            long panelTokens,
            long selfTokens,
            long tandemTokens,
            long literalRuns,
            long snps,
            long transitions,
            long transversions,
            long wobbleMod3Hits,
            long splices,
            long insertedBases,
            long deletedBases,
            long recombinationSwitches) {}

    /** Exact variant call extracted directly from compressed Mosaic tokens. */
    public record VariantRecord(
            int kind,
            int strand,
            int isTransition,
            int prevMod3,
            int donorId,
            long refPos,
            long tgtPos,
            int refLen,
            int altLen,
            byte[] refBases,
            byte[] altBases) {}

    /** Disjoint target span descriptor for O(tokens) compressed-domain k-mer queries. */
    public record SpanEvent(
            int kind,
            int strand,
            boolean hasEdits,
            int donorId,
            long refPos,
            long tgtPos,
            long tgtLen) {}

    public record MosaicInspectResult(
            MosaicInspectStats stats, List<VariantRecord> variants, List<SpanEvent> spans) {}

    private static final class MosaicObserver {
        final boolean collectVariants;
        final boolean collectSpans;
        final List<VariantRecord> variants;
        final List<SpanEvent> spans;
        long targetLen;
        long noveltyLen;
        long panelMatchBases;
        long selfMatchBases;
        long tandemMatchBases;
        long literalBases;
        long panelTokens;
        long selfTokens;
        long tandemTokens;
        long literalRuns;
        long snps;
        long transitions;
        long transversions;
        long wobbleMod3Hits;
        long splices;
        long insertedBases;
        long deletedBases;
        long recombinationSwitches;
        int lastPanelDonor = -1;

        MosaicObserver(boolean collectVariants, boolean collectSpans) {
            this.collectVariants = collectVariants;
            this.collectSpans = collectSpans;
            this.variants = collectVariants ? new ArrayList<>() : List.of();
            this.spans = collectSpans ? new ArrayList<>() : List.of();
        }

        MosaicInspectStats toStats() {
            return new MosaicInspectStats(
                    targetLen,
                    noveltyLen,
                    panelMatchBases,
                    selfMatchBases,
                    tandemMatchBases,
                    literalBases,
                    panelTokens,
                    selfTokens,
                    tandemTokens,
                    literalRuns,
                    snps,
                    transitions,
                    transversions,
                    wobbleMod3Hits,
                    splices,
                    insertedBases,
                    deletedBases,
                    recombinationSwitches);
        }
    }

    private record Subblock(int len, int[] snpOffsets, int snpCount, int dTgtAfter, int dRefAfter) {}

    private static final class CandidateChain {
        int kind; // 1 = PANEL, 2 = SELF, 3 = TANDEM
        int donorId;
        int strand;
        int startRefPos;
        int totalTgtSpan;
        int totalRefSpan;
        int totalSnps;
        int totalNovelSpliceBases;
        double mdlNetBits;
        Subblock[] subblocks;
        int numSubblocks;
    }

    /**
     * Emits a minimal zero-cut single-literal-run mosaic when a panel-less sequence has sparse
     * macro-repeats ({@code < 6.25%}), preserving 100% unbroken codon reading-frame phase and
     * in-place repeat templates for {@link FastNoveltyCoder}.
     */
    public static EncodedMosaic encodeLiteralOnlyMosaic(
            Packed2BitSequence target, MosaicConfig cfg) {
        final int n = target.length;
        ByteOut headerAndBody = new ByteOut(32);
        headerAndBody.writeByte(cfg.toHeaderByte(false));
        headerAndBody.writeVarInt(n);
        if (n == 0) {
            headerAndBody.writeVarInt(0);
            return new EncodedMosaic(
                    headerAndBody.toByteArray(),
                    new byte[0],
                    0,
                    new MosaicStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
        }
        headerAndBody.writeVarInt(n);
        BitRangeEnc enc = new BitRangeEnc(headerAndBody);
        MosaicModels models = new MosaicModels();
        enc.encodeBit(models.isLiteral, 1, 1);
        models.writeCount(enc, models.litLenSlots, 0, 0, n - 1);
        enc.finish();
        byte[] novelty = target.extractTwoBitSlice(0, n);
        return new EncodedMosaic(
                headerAndBody.toByteArray(),
                novelty,
                n,
                new MosaicStats(n, 0, 0, 0, n, 0, 0, 0, 1, 0, 0));
    }

    private static final class ParsedRegion {
        final ArrayList<ParsedToken> tokens;
        byte[] noveltyBuf;
        final int noveltyLen;
        final int statPanelBases;
        final int statSelfBases;
        final int statTandemBases;
        final int statPanelTok;
        final int statSelfTok;
        final int statTandemTok;
        final int statLitRuns;
        final int statSplices;
        final int statSnps;

        ParsedRegion(
                ArrayList<ParsedToken> tokens,
                byte[] noveltyBuf,
                int noveltyLen,
                int statPanelBases,
                int statSelfBases,
                int statTandemBases,
                int statPanelTok,
                int statSelfTok,
                int statTandemTok,
                int statLitRuns,
                int statSplices,
                int statSnps) {
            this.tokens = tokens;
            this.noveltyBuf = noveltyBuf;
            this.noveltyLen = noveltyLen;
            this.statPanelBases = statPanelBases;
            this.statSelfBases = statSelfBases;
            this.statTandemBases = statTandemBases;
            this.statPanelTok = statPanelTok;
            this.statSelfTok = statSelfTok;
            this.statTandemTok = statTandemTok;
            this.statLitRuns = statLitRuns;
            this.statSplices = statSplices;
            this.statSnps = statSnps;
        }
    }

    /**
     * Parses and encodes {@code target} against optional {@code panelIdx} (restricted to donors
     * {@code [0, maxDonorExclusive)}) plus bidirectional self-prefix and tandem periodicity.
     */
    public static EncodedMosaic encodeMosaic(
            Packed2BitSequence target,
            EimdlSyncmerIndex panelIdx,
            int maxDonorExclusive,
            int preferredTwinDonor,
            int[] donorOrderBySimilarity,
            MosaicConfig cfg) {
        return encodeMosaic(
                target, panelIdx, maxDonorExclusive, preferredTwinDonor, donorOrderBySimilarity, cfg, 1);
    }

    /** Parallelizable Stage-1 mosaic encoder preserving identical single-stream {@code .hgs2} wire format. */
    public static EncodedMosaic encodeMosaic(
            Packed2BitSequence target,
            EimdlSyncmerIndex panelIdx,
            int maxDonorExclusive,
            int preferredTwinDonor,
            int[] donorOrderBySimilarity,
            MosaicConfig cfg,
            int threads) {
        final int n = target.length;
        final boolean hasPanel = panelIdx != null && maxDonorExclusive > 0;

        ByteOut headerAndBody = new ByteOut(Math.max(64, n / 64));
        headerAndBody.writeByte(cfg.toHeaderByte(hasPanel));
        headerAndBody.writeVarInt(n);
        if (n == 0) {
            headerAndBody.writeVarInt(0);
            return new EncodedMosaic(
                    headerAndBody.toByteArray(),
                    new byte[0],
                    0,
                    new MosaicStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
        }

        // Precompute target primary (k=16) and secondary (k=13) syncmers in one fast linear pass
        TargetSyncmers tgtSync16 = extractTargetSyncmers(target, true);
        TargetSyncmers tgtSync13 =
                (hasPanel && panelIdx.hasSecondary()) ? extractTargetSyncmers(target, false) : null;

        // Build donor priority rank lookup so twin donors are prioritized during bucket probing
        int[] donorRank = null;
        if (hasPanel) {
            donorRank = new int[maxDonorExclusive];
            Arrays.fill(donorRank, maxDonorExclusive);
            if (donorOrderBySimilarity != null) {
                int r = 0;
                for (int d : donorOrderBySimilarity) {
                    if (d >= 0 && d < maxDonorExclusive && donorRank[d] == maxDonorExclusive) {
                        donorRank[d] = r++;
                    }
                }
            }
            if (preferredTwinDonor >= 0 && preferredTwinDonor < maxDonorExclusive) {
                donorRank[preferredTwinDonor] = -1;
            }
        }

        int numChunks =
                (threads > 1 && n >= 16_000_000)
                        ? Math.min(threads, Math.min(16, n / 4_000_000))
                        : 1;

        ArrayList<ParsedToken> tokens;
        byte[] noveltyBuf;
        int noveltyLen;
        int statPanelBases = 0, statSelfBases = 0, statTandemBases = 0;
        int statPanelTok = 0, statSelfTok = 0, statTandemTok = 0;
        int statLitRuns = 0, statSplices = 0, statSnps = 0;

        if (numChunks <= 1) {
            int selfSetBits = (!hasPanel && n >= 1_000_000) ? 17 : 15;
            ParsedRegion reg =
                    parseRegion(
                            target,
                            0,
                            n,
                            panelIdx,
                            maxDonorExclusive,
                            preferredTwinDonor,
                            donorRank,
                            tgtSync16,
                            tgtSync13,
                            cfg,
                            selfSetBits);
            tokens = reg.tokens;
            noveltyBuf = reg.noveltyBuf;
            noveltyLen = reg.noveltyLen;
            statPanelBases = reg.statPanelBases;
            statSelfBases = reg.statSelfBases;
            statTandemBases = reg.statTandemBases;
            statPanelTok = reg.statPanelTok;
            statSelfTok = reg.statSelfTok;
            statTandemTok = reg.statTandemTok;
            statLitRuns = reg.statLitRuns;
            statSplices = reg.statSplices;
            statSnps = reg.statSnps;
        } else {
            int selfSetBits = !hasPanel ? 16 : 15;
            ParsedRegion[] regions = new ParsedRegion[numChunks];
            final int[] finalDonorRank = donorRank;
            final TargetSyncmers finalSync16 = tgtSync16;
            final TargetSyncmers finalSync13 = tgtSync13;
            ForkJoinPool pool = new ForkJoinPool(numChunks);
            try {
                pool.submit(
                                () ->
                                        IntStream.range(0, numChunks)
                                                .parallel()
                                                .forEach(
                                                        c -> {
                                                            int cStart =
                                                                    (int) (((long) n * c) / numChunks);
                                                            int cEnd =
                                                                    (int) (((long) n * (c + 1L)) / numChunks);
                                                            regions[c] =
                                                                    parseRegion(
                                                                            target,
                                                                            cStart,
                                                                            cEnd,
                                                                            panelIdx,
                                                                            maxDonorExclusive,
                                                                            preferredTwinDonor,
                                                                            finalDonorRank,
                                                                            finalSync16,
                                                                            finalSync13,
                                                                            cfg,
                                                                            selfSetBits);
                                                        }))
                        .join();
            } finally {
                pool.shutdown();
            }

            int totalNov = 0;
            int totalTokCap = 16;
            for (ParsedRegion r : regions) {
                totalNov += r.noveltyLen;
                totalTokCap += r.tokens.size();
                statPanelBases += r.statPanelBases;
                statSelfBases += r.statSelfBases;
                statTandemBases += r.statTandemBases;
                statPanelTok += r.statPanelTok;
                statSelfTok += r.statSelfTok;
                statTandemTok += r.statTandemTok;
                statLitRuns += r.statLitRuns;
                statSplices += r.statSplices;
                statSnps += r.statSnps;
            }
            noveltyLen = totalNov;
            noveltyBuf = new byte[totalNov];
            tokens = new ArrayList<>(totalTokCap);
            int novOutOff = 0;
            ParsedToken prevTok = null;
            for (ParsedRegion r : regions) {
                if (r.noveltyLen > 0) {
                    System.arraycopy(r.noveltyBuf, 0, noveltyBuf, novOutOff, r.noveltyLen);
                    novOutOff += r.noveltyLen;
                }
                r.noveltyBuf = null;
                for (ParsedToken tok : r.tokens) {
                    if (prevTok != null && prevTok.kind == 0 && tok.kind == 0) {
                        ParsedToken merged = ParsedToken.literal(prevTok.litLen + tok.litLen);
                        tokens.set(tokens.size() - 1, merged);
                        prevTok = merged;
                        statLitRuns--;
                    } else {
                        tokens.add(tok);
                        prevTok = tok;
                    }
                }
            }
        }

        tgtSync16 = null;
        tgtSync13 = null;
        donorRank = null;

        headerAndBody.writeVarInt(noveltyLen);

        // Entropy-encode the token & edit stream symmetrically
        BitRangeEnc enc = new BitRangeEnc(headerAndBody);
        MosaicModels models = new MosaicModels();

        int[] encMtf = initMtf(maxDonorExclusive);
        int[] encLastRefEnd = hasPanel ? new int[maxDonorExclusive * 2] : null;
        int[] encLastTgtEnd = hasPanel ? new int[maxDonorExclusive * 2] : null;
        if (hasPanel) {
            Arrays.fill(encLastRefEnd, -1);
            Arrays.fill(encLastTgtEnd, -1);
        }

        boolean prevWasLiteral = false;
        boolean prevWasPanel = true;
        int cursorPos = 0;

        for (ParsedToken tok : tokens) {
            if (tok.kind == 0) {
                // Literal run: prevWasLiteral is guaranteed false here
                enc.encodeBit(models.isLiteral, prevWasPanel ? 1 : 2, 1);
                models.writeCount(enc, models.litLenSlots, 0, 0, tok.litLen - 1);
                cursorPos += tok.litLen;
                prevWasLiteral = true;
                continue;
            }

            boolean afterLit = prevWasLiteral;
            if (!prevWasLiteral) {
                enc.encodeBit(models.isLiteral, prevWasPanel ? 1 : 2, 0);
            }
            prevWasLiteral = false;

            boolean isPanel = tok.kind == 1;
            if (hasPanel) {
                if (cfg.allowSelfMatch || cfg.allowTandem) {
                    enc.encodeBit(models.isPanelMatch, prevWasPanel ? 1 : 0, isPanel ? 1 : 0);
                }
            }
            prevWasPanel = isPanel;

            if (isPanel) {
                int ds = (tok.donorId << 1) | tok.strand;
                int isMtf0 = (ds == encMtf[0]) ? 1 : 0;
                encodeMtfDonorStrand(enc, models, encMtf, maxDonorExclusive, ds, afterLit ? 1 : 0);

                int donorStart = panelIdx.donorOffsets[tok.donorId];
                int donorLen = panelIdx.donorOffsets[tok.donorId + 1] - donorStart;
                encodePanelRefPos(
                        enc,
                        models,
                        tok.startRefPos,
                        cursorPos,
                        donorLen,
                        encLastRefEnd[ds],
                        encLastTgtEnd[ds],
                        isMtf0);

                encodeSubblocks(
                        enc,
                        models,
                        cfg,
                        panelIdx.bank,
                        donorStart,
                        donorLen,
                        tok.strand,
                        tok.startRefPos,
                        target,
                        cursorPos,
                        tok.subblocks);

                encLastRefEnd[ds] = tok.startRefPos + tok.totalRefSpan;
                encLastTgtEnd[ds] = cursorPos + tok.totalTgtSpan;
                cursorPos += tok.totalTgtSpan;
            } else if (tok.kind == 3) {
                // Tandem match
                if (cfg.allowSelfMatch) {
                    enc.encodeBit(models.isTandem, 0, 1);
                }
                models.writeFixedTree(enc, models.tandemPeriodTree, tok.donorId - 1, 6);
                models.writeCount(enc, models.matchLenSlots, 0, 1, tok.totalTgtSpan - MIN_TANDEM_LEN);
                cursorPos += tok.totalTgtSpan;
            } else {
                // Bidirectional Self Match
                if (cfg.allowTandem) {
                    enc.encodeBit(models.isTandem, 0, 0);
                }
                enc.encodeBit(models.selfStrand, 0, tok.strand);
                int selfDistCode =
                        (tok.strand == 0)
                                ? (cursorPos - 1 - tok.startRefPos)
                                : tok.startRefPos;
                models.writeCount(enc, models.selfDistSlots, tok.strand, 2, selfDistCode);

                encodeSubblocks(
                        enc,
                        models,
                        cfg,
                        target,
                        0,
                        cursorPos,
                        tok.strand,
                        tok.startRefPos,
                        target,
                        cursorPos,
                        tok.subblocks);
                cursorPos += tok.totalTgtSpan;
            }
        }

        enc.finish();
        if (cursorPos != n) {
            throw new IllegalStateException("Mosaic encode cursor mismatch: " + cursorPos + " vs " + n);
        }

        MosaicStats stats =
                new MosaicStats(
                        n,
                        statPanelBases,
                        statSelfBases,
                        statTandemBases,
                        noveltyLen,
                        statPanelTok,
                        statSelfTok,
                        statTandemTok,
                        statLitRuns,
                        statSplices,
                        statSnps);

        return new EncodedMosaic(
                headerAndBody.toByteArray(),
                noveltyLen == noveltyBuf.length
                        ? noveltyBuf
                        : Arrays.copyOf(noveltyBuf, noveltyLen),
                noveltyLen,
                stats);
    }

    private static int lowerBoundPos(TargetSyncmers ts, int minPos) {
        int lo = 0;
        int hi = ts.count;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (ts.positions[mid] < minPos) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        return lo;
    }

    private static ParsedRegion parseRegion(
            Packed2BitSequence target,
            int regionStart,
            int regionEnd,
            EimdlSyncmerIndex panelIdx,
            int maxDonorExclusive,
            int preferredTwinDonor,
            int[] donorRank,
            TargetSyncmers tgtSync16,
            TargetSyncmers tgtSync13,
            MosaicConfig cfg,
            int selfSetBits) {
        final boolean hasPanel = panelIdx != null && maxDonorExclusive > 0;
        final int regionLen = regionEnd - regionStart;
        EimdlSyncmerIndex.SelfSyncmerTable selfTable =
                cfg.allowSelfMatch ? new EimdlSyncmerIndex.SelfSyncmerTable(selfSetBits) : null;

        ArrayList<ParsedToken> tokens = new ArrayList<>(Math.max(16, regionLen / 2048));
        byte[] noveltyBuf = new byte[Math.min(regionLen, Math.max(1024, Math.min(65536, regionLen / 8)))];
        int noveltyLen = 0;

        int[] mtf = initMtf(maxDonorExclusive);
        int[] lastRefEnd = hasPanel ? new int[maxDonorExclusive * 2] : null;
        int[] lastTgtEnd = hasPanel ? new int[maxDonorExclusive * 2] : null;
        if (hasPanel) {
            Arrays.fill(lastRefEnd, -1);
            Arrays.fill(lastTgtEnd, -1);
            if (preferredTwinDonor >= 0 && preferredTwinDonor < maxDonorExclusive) {
                promoteMtf(mtf, preferredTwinDonor << 1);
            }
        }

        int sync16Cursor = regionStart > 0 ? lowerBoundPos(tgtSync16, regionStart) : 0;
        int sync13Cursor =
                (tgtSync13 != null && regionStart > 0) ? lowerBoundPos(tgtSync13, regionStart) : 0;
        int selfInsertedUpTo =
                (selfTable != null && regionStart > 2_000_000)
                        ? lowerBoundPos(tgtSync16, regionStart - 2_000_000)
                        : 0;

        int litStart = regionStart;
        int pos = regionStart;

        int statPanelBases = 0, statSelfBases = 0, statTandemBases = 0;
        int statPanelTok = 0, statSelfTok = 0, statTandemTok = 0;
        int statLitRuns = 0, statSplices = 0, statSnps = 0;

        int[] scratchSnpBuf = new int[4096];
        int[] horizonCache = new int[] {-1, regionEnd};
        int[] scanCache = new int[] {-1, -1, -1}; // deadPanelPIdx, deadSecS13, deadSelfSIdx
        int[] bucketScratch = new int[MAX_BUCKET_PROBES * 4 + 1];

        while (pos < regionEnd) {
            if (selfTable != null) {
                int inserted = 0;
                while (selfInsertedUpTo < tgtSync16.count
                        && tgtSync16.positions[selfInsertedUpTo] + Packed2BitSequence.K_PRIMARY
                                <= pos) {
                    selfTable.insert(
                            tgtSync16.canonKmers[selfInsertedUpTo] & 0xFFFFFFFFL,
                            tgtSync16.positions[selfInsertedUpTo],
                            tgtSync16.strands[selfInsertedUpTo]);
                    selfInsertedUpTo++;
                    inserted++;
                }
                if (inserted > 0) {
                    scanCache[2] = -1;
                }
            }

            while (sync16Cursor < tgtSync16.count && tgtSync16.positions[sync16Cursor] < pos) {
                sync16Cursor++;
            }
            if (tgtSync13 != null) {
                while (sync13Cursor < tgtSync13.count && tgtSync13.positions[sync13Cursor] < pos) {
                    sync13Cursor++;
                }
            }

            CandidateChain best =
                    findBestChainAt(
                            pos,
                            litStart,
                            regionEnd,
                            target,
                            panelIdx,
                            maxDonorExclusive,
                            donorRank,
                            mtf,
                            lastRefEnd,
                            lastTgtEnd,
                            selfTable,
                            tgtSync16,
                            sync16Cursor,
                            tgtSync13,
                            sync13Cursor,
                            cfg,
                            scratchSnpBuf,
                            horizonCache,
                            scanCache,
                            bucketScratch);

            if (best != null && best.totalTgtSpan < 256 && pos + 8 + MIN_MATCH_LEN < regionEnd) {
                int nextPos = findNextSyncmerPos(tgtSync16, sync16Cursor, pos + 1, pos + 7);
                if (nextPos > pos) {
                    CandidateChain lazy =
                            findBestChainAt(
                                    nextPos,
                                    litStart,
                                    regionEnd,
                                    target,
                                    panelIdx,
                                    maxDonorExclusive,
                                    donorRank,
                                    mtf,
                                    lastRefEnd,
                                    lastTgtEnd,
                                    selfTable,
                                    tgtSync16,
                                    sync16Cursor,
                                    tgtSync13,
                                    sync13Cursor,
                                    cfg,
                                    scratchSnpBuf,
                                    horizonCache,
                                    scanCache,
                                    bucketScratch);
                    if (lazy != null
                            && lazy.mdlNetBits > best.mdlNetBits + 14.0
                            && lazy.totalTgtSpan > best.totalTgtSpan + 24) {
                        pos = nextPos;
                        continue;
                    }
                }
            }

            if (best == null) {
                int step = 1;
                if (pos - litStart >= 32) {
                    int nextS16 =
                            (sync16Cursor < tgtSync16.count)
                                    ? tgtSync16.positions[sync16Cursor]
                                    : regionEnd;
                    int safeJump = nextS16 - (pos + Packed2BitSequence.W_PRIMARY - 1);
                    if (safeJump > 1) {
                        step = Math.min(safeJump, cfg.allowTandem ? 4 : 8);
                    }
                }
                pos = Math.min(regionEnd, pos + step);
                continue;
            }

            if (pos > litStart) {
                int litLen = pos - litStart;
                noveltyBuf = appendNovelBases(noveltyBuf, noveltyLen, target, litStart, litLen);
                noveltyLen += litLen;
                tokens.add(ParsedToken.literal(litLen));
                statLitRuns++;
            }

            int curT = pos;
            for (int b = 0; b < best.numSubblocks; b++) {
                Subblock sb = best.subblocks[b];
                curT += sb.len;
                if (sb.dTgtAfter > 0) {
                    noveltyBuf =
                            appendNovelBases(noveltyBuf, noveltyLen, target, curT, sb.dTgtAfter);
                    noveltyLen += sb.dTgtAfter;
                    curT += sb.dTgtAfter;
                }
            }

            tokens.add(ParsedToken.fromChain(pos, best));

            if (best.kind == 1) {
                int ds = (best.donorId << 1) | best.strand;
                promoteMtf(mtf, ds);
                lastRefEnd[ds] = best.startRefPos + best.totalRefSpan;
                lastTgtEnd[ds] = pos + best.totalTgtSpan;
                statPanelBases += best.totalTgtSpan - best.totalNovelSpliceBases;
                statPanelTok++;
            } else if (best.kind == 2) {
                statSelfBases += best.totalTgtSpan - best.totalNovelSpliceBases;
                statSelfTok++;
            } else {
                statTandemBases += best.totalTgtSpan;
                statTandemTok++;
            }
            statSplices += Math.max(0, best.numSubblocks - 1);
            statSnps += best.totalSnps;

            pos += best.totalTgtSpan;
            litStart = pos;
        }

        if (litStart < regionEnd) {
            int litLen = regionEnd - litStart;
            noveltyBuf = appendNovelBases(noveltyBuf, noveltyLen, target, litStart, litLen);
            noveltyLen += litLen;
            tokens.add(ParsedToken.literal(litLen));
            statLitRuns++;
        }

        return new ParsedRegion(
                tokens,
                noveltyBuf,
                noveltyLen,
                statPanelBases,
                statSelfBases,
                statTandemBases,
                statPanelTok,
                statSelfTok,
                statTandemTok,
                statLitRuns,
                statSplices,
                statSnps);
    }

    /** Inspects just the header of {@code mosaicBytes} to read {@code (targetLen, noveltyLen)}. */
    public static int[] peekLengths(byte[] mosaicBytes) {
        ByteIn in = new ByteIn(mosaicBytes);
        in.readByte(); // configFlags
        int targetLen = in.readVarInt();
        int noveltyLen = in.readVarInt();
        return new int[] {targetLen, noveltyLen};
    }

    /**
     * Reconstructs exact 2-bit sequence from {@code mosaicBytes} + decoded {@code noveltyBases}
     * against optional {@code panelBank} ({@code donorOffsets}). Runs at memory-copy speed (hundreds
     * of MB/s, zero hash tables needed!).
     */
    public static MosaicInspectResult inspectMosaic(
            byte[] mosaicBytes,
            byte[] noveltyBases,
            Packed2BitSequence panelBank,
            int[] donorOffsets,
            int maxDonorExclusive,
            boolean collectVariants,
            boolean collectSpans) {
        MosaicObserver obs = new MosaicObserver(collectVariants, collectSpans);
        decodeMosaicImpl(
                mosaicBytes,
                noveltyBases,
                panelBank,
                donorOffsets,
                maxDonorExclusive,
                true,
                obs);
        return new MosaicInspectResult(obs.toStats(), obs.variants, obs.spans);
    }

    public static byte[] decodeMosaicToTwoBit(
            byte[] mosaicBytes,
            byte[] noveltyBases,
            Packed2BitSequence panelBank,
            int[] donorOffsets,
            int maxDonorExclusive) {
        return decodeMosaicImpl(
                mosaicBytes,
                noveltyBases,
                panelBank,
                donorOffsets,
                maxDonorExclusive,
                false,
                null);
    }

    private static byte[] decodeMosaicImpl(
            byte[] mosaicBytes,
            byte[] noveltyBases,
            Packed2BitSequence panelBank,
            int[] donorOffsets,
            int maxDonorExclusive,
            boolean inspectOnly,
            MosaicObserver obs) {
        ByteIn in = new ByteIn(mosaicBytes);
        int flags = in.readByte() & 0xFF;
        boolean hasPanel = (flags & 1) != 0;
        MosaicConfig cfg = MosaicConfig.fromHeaderByte(flags);
        int targetLen = in.readVarInt();
        int noveltyLen = in.readVarInt();

        if (targetLen < 0 || noveltyLen < 0 || noveltyLen > noveltyBases.length) {
            throw new IllegalStateException("Malformed mosaic lengths");
        }
        if (obs != null) {
            obs.targetLen = targetLen;
            obs.noveltyLen = noveltyLen;
        }
        if (targetLen == 0) {
            return inspectOnly ? null : new byte[0];
        }
        if (hasPanel && (panelBank == null || maxDonorExclusive <= 0)) {
            throw new IllegalStateException("Archive requires reference panel for mosaic decode");
        }

        boolean needOut = !inspectOnly || cfg.allowSelfMatch || cfg.allowTandem || !cfg.allowKimuraV4;
        byte[] out = needOut ? new byte[targetLen] : null;
        BitRangeDec dec = new BitRangeDec(in);
        MosaicModels models = new MosaicModels();

        int[] mtf = initMtf(maxDonorExclusive);
        int[] lastRefEnd = hasPanel ? new int[maxDonorExclusive * 2] : null;
        int[] lastTgtEnd = hasPanel ? new int[maxDonorExclusive * 2] : null;
        if (hasPanel) {
            Arrays.fill(lastRefEnd, -1);
            Arrays.fill(lastTgtEnd, -1);
        }

        boolean prevWasLiteral = false;
        boolean prevWasPanel = true;
        int pos = 0;
        int novCursor = 0;

        while (pos < targetLen) {
            boolean isLit = false;
            if (!prevWasLiteral) {
                isLit = dec.decodeBit(models.isLiteral, prevWasPanel ? 1 : 2) == 1;
            }
            if (isLit) {
                int litLen = models.readCount(dec, models.litLenSlots, 0, 0) + 1;
                if (pos + litLen > targetLen || novCursor + litLen > noveltyLen) {
                    throw new IllegalStateException("Literal run exceeds target or novelty bounds");
                }
                if (out != null) {
                    System.arraycopy(noveltyBases, novCursor, out, pos, litLen);
                }
                if (obs != null) {
                    obs.literalRuns++;
                    obs.literalBases += litLen;
                    if (obs.collectSpans) {
                        obs.spans.add(new SpanEvent(0, 0, true, -1, novCursor, pos, litLen));
                    }
                    if (obs.collectVariants) {
                        byte[] alt =
                                (litLen <= 65536)
                                        ? Arrays.copyOfRange(noveltyBases, novCursor, novCursor + litLen)
                                        : null;
                        obs.variants.add(
                                new VariantRecord(
                                        4, 0, 0, 0, -1, novCursor, pos, 0, litLen, null, alt));
                    }
                }
                novCursor += litLen;
                pos += litLen;
                prevWasLiteral = true;
                continue;
            }
            boolean afterLit = prevWasLiteral;
            prevWasLiteral = false;

            boolean isPanel = false;
            if (hasPanel) {
                if (cfg.allowSelfMatch || cfg.allowTandem) {
                    isPanel = dec.decodeBit(models.isPanelMatch, prevWasPanel ? 1 : 0) == 1;
                } else {
                    isPanel = true;
                }
            }
            prevWasPanel = isPanel;

            if (isPanel) {
                int prevMtf0 = mtf[0];
                int ds = decodeMtfDonorStrand(dec, models, mtf, maxDonorExclusive, afterLit ? 1 : 0);
                int isMtf0 = (ds == prevMtf0) ? 1 : 0;
                int donorId = ds >>> 1;
                int strand = ds & 1;
                int donorStart = donorOffsets[donorId];
                int donorLen = donorOffsets[donorId + 1] - donorStart;
                int startRefPos =
                        decodePanelRefPos(
                                dec,
                                models,
                                pos,
                                donorLen,
                                lastRefEnd[ds],
                                lastTgtEnd[ds],
                                isMtf0);

                if (obs != null) {
                    obs.panelTokens++;
                    if (obs.lastPanelDonor >= 0 && obs.lastPanelDonor != donorId) {
                        obs.recombinationSwitches++;
                    }
                    obs.lastPanelDonor = donorId;
                }
                int[] adv =
                        decodeSubblocksFromPanel(
                                dec,
                                models,
                                cfg,
                                panelBank,
                                donorStart,
                                donorLen,
                                strand,
                                startRefPos,
                                out,
                                pos,
                                targetLen,
                                noveltyBases,
                                novCursor,
                                noveltyLen,
                                donorId,
                                obs);
                lastRefEnd[ds] = adv[0];
                lastTgtEnd[ds] = adv[1];
                pos = adv[1];
                novCursor = adv[2];
            } else {
                boolean isTandem;
                if (cfg.allowSelfMatch && cfg.allowTandem) {
                    isTandem = dec.decodeBit(models.isTandem, 0) == 1;
                } else {
                    isTandem = cfg.allowTandem;
                }
                if (isTandem) {
                    int period = models.readFixedTree(dec, models.tandemPeriodTree, 6) + 1;
                    int len = models.readCount(dec, models.matchLenSlots, 0, 1) + MIN_TANDEM_LEN;
                    if (period > pos || pos + len > targetLen) {
                        throw new IllegalStateException("Malformed tandem copy geometry");
                    }
                    if (out != null) {
                        int src = pos - period;
                        for (int k = 0; k < len; k++) {
                            out[pos + k] = out[src + k];
                        }
                    }
                    if (obs != null) {
                        obs.tandemTokens++;
                        obs.tandemMatchBases += len;
                        if (obs.collectSpans) {
                            obs.spans.add(new SpanEvent(3, 0, false, -1, period, pos, len));
                        }
                    }
                    pos += len;
                } else {
                    int strand = dec.decodeBit(models.selfStrand, 0);
                    int distCode = models.readCount(dec, models.selfDistSlots, strand, 2);
                    int startRefPos = (strand == 0) ? (pos - 1 - distCode) : distCode;
                    int selfRefLen = pos;
                    if (startRefPos < 0 || startRefPos >= selfRefLen) {
                        throw new IllegalStateException("Malformed self copy start position");
                    }
                    if (obs != null) {
                        obs.selfTokens++;
                    }
                    int[] adv =
                            decodeSubblocksFromSelf(
                                    dec,
                                    models,
                                    cfg,
                                    out,
                                    selfRefLen,
                                    strand,
                                    startRefPos,
                                    pos,
                                    targetLen,
                                    noveltyBases,
                                    novCursor,
                                    noveltyLen,
                                    obs);
                    pos = adv[1];
                    novCursor = adv[2];
                }
            }
        }

        if (novCursor != noveltyLen) {
            throw new IllegalStateException(
                    "Unconsumed novelty bases after mosaic decode: " + novCursor + " vs " + noveltyLen);
        }
        return inspectOnly ? null : out;
    }

    // =========================================================================
    // Candidate Search, MDL Pricing, & Micro-Indel Splice Chaining
    // =========================================================================

    private static CandidateChain findBestChainAt(
            int pos,
            int litStart,
            int tgtLimit,
            Packed2BitSequence target,
            EimdlSyncmerIndex panelIdx,
            int maxDonorExclusive,
            int[] donorRank,
            int[] mtf,
            int[] lastRefEnd,
            int[] lastTgtEnd,
            EimdlSyncmerIndex.SelfSyncmerTable selfTable,
            TargetSyncmers tgtSync16,
            int sync16Cursor,
            TargetSyncmers tgtSync13,
            int sync13Cursor,
            MosaicConfig cfg,
            int[] scratchSnpBuf,
            int[] horizonCache,
            int[] scanCache,
            int[] bucketScratch) {
        final int n = tgtLimit;
        if (pos + MIN_TANDEM_LEN > n) {
            return null;
        }

        CandidateChain best = null;

        // 1. Fast O(1) Tandem Periodicity Check (periods 1..64 via 4-word SWAR shifts)
        if (cfg.allowTandem && pos >= 1) {
            int maxP = Math.min(MAX_TANDEM_PERIOD, pos);
            int bestP = -1;
            int bestTLen = 0;
            int tgt16 = (int) target.word32Forward(pos);
            if (maxP == MAX_TANDEM_PERIOD) {
                long w1 = target.word32Forward(pos - 13);
                for (int j = 12; j >= 0; j--) {
                    if (((int) (w1 >>> (j << 1))) == tgt16) {
                        int p = 13 - j;
                        int tLen =
                                Packed2BitSequence.countExactMatch(
                                        target, 0, n, 0, pos - p, target, pos, n - pos);
                        if (tLen >= MIN_TANDEM_LEN && tLen > bestTLen) {
                            bestTLen = tLen;
                            bestP = p;
                            if (tLen >= 4096) break;
                        }
                    }
                }
                if (bestTLen < 4096) {
                    long w2 = target.word32Forward(pos - 30);
                    for (int j = 16; j >= 0; j--) {
                        if (((int) (w2 >>> (j << 1))) == tgt16) {
                            int p = 30 - j;
                            int tLen =
                                    Packed2BitSequence.countExactMatch(
                                            target, 0, n, 0, pos - p, target, pos, n - pos);
                            if (tLen >= MIN_TANDEM_LEN && tLen > bestTLen) {
                                bestTLen = tLen;
                                bestP = p;
                                if (tLen >= 4096) break;
                            }
                        }
                    }
                }
                if (bestTLen < 4096) {
                    long w3 = target.word32Forward(pos - 47);
                    for (int j = 16; j >= 0; j--) {
                        if (((int) (w3 >>> (j << 1))) == tgt16) {
                            int p = 47 - j;
                            int tLen =
                                    Packed2BitSequence.countExactMatch(
                                            target, 0, n, 0, pos - p, target, pos, n - pos);
                            if (tLen >= MIN_TANDEM_LEN && tLen > bestTLen) {
                                bestTLen = tLen;
                                bestP = p;
                                if (tLen >= 4096) break;
                            }
                        }
                    }
                }
                if (bestTLen < 4096) {
                    long w4 = target.word32Forward(pos - 64);
                    for (int j = 16; j >= 0; j--) {
                        if (((int) (w4 >>> (j << 1))) == tgt16) {
                            int p = 64 - j;
                            int tLen =
                                    Packed2BitSequence.countExactMatch(
                                            target, 0, n, 0, pos - p, target, pos, n - pos);
                            if (tLen >= MIN_TANDEM_LEN && tLen > bestTLen) {
                                bestTLen = tLen;
                                bestP = p;
                                if (tLen >= 4096) break;
                            }
                        }
                    }
                }
            } else {
                for (int p = 1; p <= maxP; p++) {
                    if (((int) target.word32Forward(pos - p)) != tgt16) {
                        continue;
                    }
                    int tLen =
                            Packed2BitSequence.countExactMatch(
                                    target, 0, n, 0, pos - p, target, pos, n - pos);
                    if (tLen >= MIN_TANDEM_LEN && tLen > bestTLen) {
                        bestTLen = tLen;
                        bestP = p;
                        if (tLen >= 4096) break;
                    }
                }
            }
            int minTandem = (panelIdx != null && maxDonorExclusive > 0) ? MIN_TANDEM_LEN : 36;
            if (bestTLen >= minTandem) {
                double cost = 11.0 + approxCountBits(bestTLen - MIN_TANDEM_LEN);
                double gain = bestTLen * 1.865 - cost;
                if (gain > 6.5) {
                    CandidateChain tc = new CandidateChain();
                    tc.kind = 3;
                    tc.donorId = bestP;
                    tc.strand = 0;
                    tc.startRefPos = pos - bestP;
                    tc.totalTgtSpan = bestTLen;
                    tc.totalRefSpan = bestTLen;
                    tc.mdlNetBits = gain;
                    tc.subblocks = new Subblock[] {new Subblock(bestTLen, new int[0], 0, 0, 0)};
                    tc.numSubblocks = 1;
                    best = tc;
                    if (bestTLen >= 512) {
                        return best;
                    }
                }
            }
        }

        if (pos + MIN_MATCH_LEN > n) {
            return best;
        }

        // 2. Sticky Active-Donor Collinear / Indel-Shifted Probe (O(1) Twin Lock!)
        boolean multiDonor = (panelIdx != null && maxDonorExclusive > 1);
        if (panelIdx != null && maxDonorExclusive > 0) {
            for (int mSlot = 0; mSlot < 2; mSlot++) {
                int ds = mtf[mSlot];
                int dId = ds >>> 1;
                int strand = ds & 1;
                if (dId >= maxDonorExclusive || lastRefEnd[ds] < 0) continue;
                int dStart = panelIdx.donorOffsets[dId];
                int dLen = panelIdx.donorOffsets[dId + 1] - dStart;
                int predCol = lastRefEnd[ds] + (pos - lastTgtEnd[ds]);
                int predIns = lastRefEnd[ds];

                best =
                        probePanelStart(
                                pos,
                                n,
                                target,
                                panelIdx,
                                dId,
                                dStart,
                                dLen,
                                strand,
                                predCol,
                                mtf,
                                lastRefEnd[ds],
                                lastTgtEnd[ds],
                                cfg,
                                multiDonor,
                                scratchSnpBuf,
                                best);
                if (best != null
                        && best.totalTgtSpan >= 4096
                        && (multiDonor
                                ? (best.totalSnps == 0 && best.numSubblocks == 1)
                                : (best.totalSnps * 128 < best.totalTgtSpan))) {
                    return best;
                }
                if (predIns != predCol) {
                    best =
                            probePanelStart(
                                    pos,
                                    n,
                                    target,
                                    panelIdx,
                                    dId,
                                    dStart,
                                    dLen,
                                    strand,
                                    predIns,
                                    mtf,
                                    lastRefEnd[ds],
                                    lastTgtEnd[ds],
                                    cfg,
                                    multiDonor,
                                    scratchSnpBuf,
                                    best);
                    if (best != null
                            && best.totalTgtSpan >= 4096
                            && (multiDonor
                                    ? (best.totalSnps == 0 && best.numSubblocks == 1)
                                    : (best.totalSnps * 128 < best.totalTgtSpan))) {
                        return best;
                    }
                }
            }
        }

        // Locate first primary syncmer within [pos .. pos + W_PRIMARY - 1]
        int sIdx = sync16Cursor;
        while (sIdx < tgtSync16.count && tgtSync16.positions[sIdx] < pos) {
            sIdx++;
        }
        boolean hasNearSync16 =
                sIdx < tgtSync16.count
                        && tgtSync16.positions[sIdx] < pos + Packed2BitSequence.W_PRIMARY
                        && tgtSync16.positions[sIdx] + Packed2BitSequence.K_PRIMARY <= n;

        int gear2Limit = -1;

        // 3. Panel CSR Canonical Open-Syncmer Lookup (Two-Gear Coalescent Haplotype Copying!)
        if (panelIdx != null && maxDonorExclusive > 0) {
            int pIdx = sIdx;
            boolean hasNearPanelSync16;
            if (panelIdx.primarySubsampleMask == 0) {
                hasNearPanelSync16 = hasNearSync16;
            } else {
                int maxPanelWin = Packed2BitSequence.W_PRIMARY << panelIdx.primarySubsampleShift;
                while (pIdx < tgtSync16.count
                        && tgtSync16.positions[pIdx] < pos + maxPanelWin
                        && !panelIdx.acceptsPrimary(tgtSync16.canonKmers[pIdx] & 0xFFFFFFFFL)) {
                    pIdx++;
                }
                hasNearPanelSync16 =
                        pIdx < tgtSync16.count
                                && tgtSync16.positions[pIdx] < pos + maxPanelWin
                                && tgtSync16.positions[pIdx] + Packed2BitSequence.K_PRIMARY <= n
                                && panelIdx.acceptsPrimary(tgtSync16.canonKmers[pIdx] & 0xFFFFFFFFL);
            }
            if (hasNearPanelSync16 && pIdx != scanCache[0]) {
                int syncOff = tgtSync16.positions[pIdx] - pos;
                long qCanon = tgtSync16.canonKmers[pIdx] & 0xFFFFFFFFL;
                int qStrand = tgtSync16.strands[pIdx];
                int bucket = panelIdx.bucketForPrimary(qCanon);
                short qSig = panelIdx.sigForPrimary(qCanon);

                int bStart = panelIdx.bucketStart(bucket);
                int bEnd = panelIdx.bucketEnd(bucket);
                if (bStart >= bEnd) {
                    scanCache[0] = pIdx;
                } else {
                    best =
                            probePanelBucket(
                                    pos,
                                    n,
                                    syncOff,
                                    Packed2BitSequence.K_PRIMARY,
                                    qStrand,
                                    qSig,
                                    bStart,
                                    bEnd,
                                    false,
                                    target,
                                    panelIdx,
                                    maxDonorExclusive,
                                    donorRank,
                                    mtf,
                                    lastRefEnd,
                                    lastTgtEnd,
                                    cfg,
                                    multiDonor,
                                    scratchSnpBuf,
                                    bucketScratch,
                                    best);
                    if (bucketScratch[MAX_BUCKET_PROBES * 4] == 0) {
                        scanCache[0] = pIdx;
                    }
                    if (best != null
                            && best.totalTgtSpan >= 2048
                            && (multiDonor
                                    ? (best.totalSnps == 0 && best.numSubblocks == 1)
                                    : (best.totalSnps * 64 <= best.totalTgtSpan))) {
                        return best;
                    }
                    // Gear 2 (Horizon-Bounded Diverged Homolog Gear): never bulldoze downstream clonal blocks!
                    if (multiDonor && pIdx != scanCache[0] && (best == null || best.totalTgtSpan < 96)) {
                        gear2Limit =
                                findClonalReentryHorizon(
                                        pos,
                                        n,
                                        sIdx,
                                        target,
                                        tgtSync16,
                                        panelIdx,
                                        maxDonorExclusive,
                                        horizonCache);
                        if (gear2Limit - pos >= MIN_MATCH_LEN) {
                            for (int mSlot = 0; mSlot < 2; mSlot++) {
                                int ds = mtf[mSlot];
                                int dId = ds >>> 1;
                                int strand = ds & 1;
                                if (dId >= maxDonorExclusive || lastRefEnd[ds] < 0) continue;
                                int dStart = panelIdx.donorOffsets[dId];
                                int dLen = panelIdx.donorOffsets[dId + 1] - dStart;
                                int predCol = lastRefEnd[ds] + (pos - lastTgtEnd[ds]);
                                best =
                                        probePanelStart(
                                                pos,
                                                gear2Limit,
                                                target,
                                                panelIdx,
                                                dId,
                                                dStart,
                                                dLen,
                                                strand,
                                                predCol,
                                                mtf,
                                                lastRefEnd[ds],
                                                lastTgtEnd[ds],
                                                cfg,
                                                false,
                                                scratchSnpBuf,
                                                best);
                            }
                            best =
                                    probePanelBucket(
                                            pos,
                                            gear2Limit,
                                            syncOff,
                                            Packed2BitSequence.K_PRIMARY,
                                            qStrand,
                                            qSig,
                                            bStart,
                                            bEnd,
                                            false,
                                            target,
                                            panelIdx,
                                            maxDonorExclusive,
                                            donorRank,
                                            mtf,
                                            lastRefEnd,
                                            lastTgtEnd,
                                            cfg,
                                            false,
                                            scratchSnpBuf,
                                            bucketScratch,
                                            best);
                        }
                    }
                }
            }
        }

        // 4. Secondary k=13 Syncmer Diverged-Homology Rescue (horizon-bounded in multi-donor panels)
        if (panelIdx != null
                && maxDonorExclusive > 0
                && tgtSync13 != null
                && (best == null || best.mdlNetBits < 24.0)
                && (pos - litStart >= 16 || !hasNearSync16)) {
            if (multiDonor && gear2Limit < 0) {
                gear2Limit =
                        findClonalReentryHorizon(
                                pos,
                                n,
                                sIdx,
                                target,
                                tgtSync16,
                                panelIdx,
                                maxDonorExclusive,
                                horizonCache);
            }
            int secLimit = multiDonor ? gear2Limit : n;
            if (secLimit - pos >= MIN_MATCH_LEN) {
                int s13 = sync13Cursor;
                while (s13 < tgtSync13.count && tgtSync13.positions[s13] < pos) {
                    s13++;
                }
                if (s13 < tgtSync13.count
                        && s13 != scanCache[1]
                        && tgtSync13.positions[s13] < pos + Packed2BitSequence.W_SECONDARY
                        && tgtSync13.positions[s13] + Packed2BitSequence.K_SECONDARY <= secLimit) {
                    int syncOff = tgtSync13.positions[s13] - pos;
                    long qCanon13 = tgtSync13.canonKmers[s13] & 0xFFFFFFFFL;
                    int b13 = panelIdx.secBucketFor(qCanon13);
                    if (b13 < 0) {
                        scanCache[1] = s13;
                    } else {
                        short qSig13 = panelIdx.secSigFor(qCanon13);
                        int bStart = panelIdx.secBucketStart(b13);
                        int bEnd = panelIdx.secBucketEnd(b13);
                        if (bStart >= bEnd) {
                            scanCache[1] = s13;
                        } else {
                            best =
                                    probePanelBucket(
                                            pos,
                                            secLimit,
                                            syncOff,
                                            Packed2BitSequence.K_SECONDARY,
                                            tgtSync13.strands[s13],
                                            qSig13,
                                            bStart,
                                            bEnd,
                                            true,
                                            target,
                                            panelIdx,
                                            maxDonorExclusive,
                                            donorRank,
                                            mtf,
                                            lastRefEnd,
                                            lastTgtEnd,
                                            cfg,
                                            false,
                                            scratchSnpBuf,
                                            bucketScratch,
                                            best);
                            if (bucketScratch[MAX_BUCKET_PROBES * 4] == 0) {
                                scanCache[1] = s13;
                            }
                        }
                    }
                }
            }
        }

        // 5. Bidirectional (+/-) Exact & Approximate Self-Homology Lookup
        if (selfTable != null && hasNearSync16 && sIdx != scanCache[2] && pos >= MIN_MATCH_LEN) {
            if (multiDonor && gear2Limit < 0 && (best == null || best.totalTgtSpan < 96)) {
                gear2Limit =
                        findClonalReentryHorizon(
                                pos,
                                n,
                                sIdx,
                                target,
                                tgtSync16,
                                panelIdx,
                                maxDonorExclusive,
                                horizonCache);
            }
            boolean selfClonal =
                    multiDonor
                            && ((best != null && best.totalTgtSpan >= 96)
                                    || (gear2Limit >= 0 && gear2Limit - pos < MIN_MATCH_LEN));
            int selfLimit = (!multiDonor || selfClonal) ? n : gear2Limit;

            int syncOff = tgtSync16.positions[sIdx] - pos;
            long qCanon = tgtSync16.canonKmers[sIdx] & 0xFFFFFFFFL;
            int qStrand = tgtSync16.strands[sIdx];
            int setBase = selfTable.setBase(qCanon);
            short qSig = selfTable.sigFor(qCanon);

            boolean anySelfSeedHit = false;
            for (int w = 0; w < EimdlSyncmerIndex.SelfSyncmerTable.WAYS; w++) {
                int slot = setBase + w;
                int packed = selfTable.slotPackedPosStrand(slot);
                if (packed == -1 || selfTable.slotSig(slot) != qSig) continue;
                int hitFwdPos = packed & 0x7FFFFFFF;
                int hitCanonStrand = packed >>> 31;
                if (hitFwdPos + Packed2BitSequence.K_PRIMARY > pos) continue;

                int relStrand = qStrand ^ hitCanonStrand;
                int selfRefLen = pos; // [0 .. pos) is immutable already-decoded prefix!
                int refSyncPos =
                        (relStrand == 0)
                                ? hitFwdPos
                                : (selfRefLen - Packed2BitSequence.K_PRIMARY - hitFwdPos);
                int candRefStart = refSyncPos - syncOff;
                if (candRefStart < 0 || candRefStart + MIN_MATCH_LEN > selfRefLen) continue;
                if (pos + syncOff + Packed2BitSequence.K_PRIMARY > selfLimit) continue;

                // Verify exact syncmer seed match before full extension
                if (Packed2BitSequence.countExactMatch(
                                target,
                                0,
                                selfRefLen,
                                relStrand,
                                refSyncPos,
                                target,
                                pos + syncOff,
                                Packed2BitSequence.K_PRIMARY)
                        < Packed2BitSequence.K_PRIMARY) {
                    continue;
                }
                anySelfSeedHit = true;

                // Left-adjust if syncOff > 0 had leading mismatches
                int adjShift = 0;
                while (adjShift < syncOff
                        && target.baseAtStrand(0, selfRefLen, relStrand, candRefStart + adjShift)
                                != target.baseAt(pos + adjShift)) {
                    adjShift++;
                }
                if (adjShift > 0) {
                    continue; // Wait until pos reaches pos + adjShift cleanly!
                }

                CandidateChain sc =
                        buildChain(
                                2,
                                0,
                                relStrand,
                                candRefStart,
                                target,
                                0,
                                selfRefLen,
                                target,
                                pos,
                                selfLimit,
                                mtf,
                                -1,
                                -1,
                                cfg,
                                selfClonal,
                                scratchSnpBuf);
                boolean passSelfMin =
                        (panelIdx != null && maxDonorExclusive > 0)
                                || (sc != null
                                        && sc.totalTgtSpan >= 64
                                        && sc.mdlNetBits >= 64.0);
                if (sc != null && passSelfMin && (best == null || sc.mdlNetBits > best.mdlNetBits)) {
                    best = sc;
                }
            }
            if (!anySelfSeedHit) {
                scanCache[2] = sIdx;
            }
        }

        return best;
    }

    /**
     * Locates the earliest downstream coordinate {@code > pos} where any reference panel donor
     * resumes a clonal exact haplotype of length {@code >= 192 bp}. Cached in {@code horizonCache}
     * for O(1) amortized cost across novel/diverged islands.
     */
    private static int findClonalReentryHorizon(
            int pos,
            int tgtLimit,
            int sync16Cursor,
            Packed2BitSequence target,
            TargetSyncmers tgtSync16,
            EimdlSyncmerIndex panelIdx,
            int maxDonorExclusive,
            int[] horizonCache) {
        final int n = tgtLimit;
        if (pos < horizonCache[0]) {
            return horizonCache[1];
        }
        final int minClonalLen = 192;
        if (pos + minClonalLen >= n) {
            horizonCache[0] = n;
            horizonCache[1] = n;
            return n;
        }
        final int maxScanAhead = 65536;
        final int scanEnd = Math.min(n - 144, pos + maxScanAhead);
        int nextCheckPos = pos + 8;

        for (int idx = sync16Cursor; idx < tgtSync16.count; idx++) {
            int sp = tgtSync16.positions[idx];
            if (sp < nextCheckPos) continue;
            if (sp > scanEnd) break;
            long qCanon = tgtSync16.canonKmers[idx] & 0xFFFFFFFFL;
            if (!panelIdx.acceptsPrimary(qCanon)) continue;
            nextCheckPos = sp + 32;

            int qStrand = tgtSync16.strands[idx];
            int bucket = panelIdx.bucketForPrimary(qCanon);
            short qSig = panelIdx.sigForPrimary(qCanon);
            int bStart = panelIdx.bucketStart(bucket);
            int bEnd = panelIdx.bucketEnd(bucket);
            if (bStart >= bEnd) continue;

            int earliestStart = Integer.MAX_VALUE;
            int probes = 0;
            for (int b = bStart; b < bEnd && probes < 24; b++) {
                if (panelIdx.entrySig(b) != qSig) continue;
                int packed = panelIdx.entryPackedPosAndStrand(b);
                int globalPos = packed & 0x7FFFFFFF;
                int dId = panelIdx.donorOfGlobalPos(globalPos);
                if (dId >= maxDonorExclusive) break;
                probes++;

                int hitCanonStrand = packed >>> 31;
                int dStart = panelIdx.donorOffsets[dId];
                int dLen = panelIdx.donorOffsets[dId + 1] - dStart;
                int relFwdPos = globalPos - dStart;
                int strand = qStrand ^ hitCanonStrand;
                int refSyncPos =
                        (strand == 0) ? relFwdPos : (dLen - Packed2BitSequence.K_PRIMARY - relFwdPos);
                if (refSyncPos < 0 || refSyncPos + 144 > dLen) continue;

                int maxFwd = Math.min(256, Math.min(dLen - refSyncPos, n - sp));
                int fwdExact =
                        Packed2BitSequence.countExactMatch(
                                panelIdx.bank, dStart, dLen, strand, refSyncPos, target, sp, maxFwd);
                if (fwdExact < 144) continue;

                int back = 0;
                int maxBack = Math.min(refSyncPos, sp);
                while (back < maxBack
                        && panelIdx.bank.baseAtStrand(dStart, dLen, strand, refSyncPos - 1 - back)
                                == target.baseAt(sp - 1 - back)) {
                    back++;
                }
                if (back + fwdExact >= minClonalLen) {
                    int rawStart = sp - back;
                    if (rawStart < earliestStart) {
                        earliestStart = rawStart;
                    }
                }
            }
            if (earliestStart != Integer.MAX_VALUE) {
                if (earliestStart > pos) {
                    horizonCache[0] = earliestStart;
                    horizonCache[1] = earliestStart;
                    return earliestStart;
                } else {
                    horizonCache[0] = pos + 1;
                    horizonCache[1] = pos + 1;
                    return pos + 1;
                }
            }
        }

        if (pos + maxScanAhead >= n - 144) {
            horizonCache[0] = n;
            horizonCache[1] = n;
            return n;
        }
        int hor = Math.min(n, pos + maxScanAhead);
        horizonCache[0] = pos + (maxScanAhead >>> 1);
        horizonCache[1] = hor;
        return hor;
    }

    private static CandidateChain probePanelStart(
            int pos,
            int tgtLimit,
            Packed2BitSequence target,
            EimdlSyncmerIndex panelIdx,
            int dId,
            int dStart,
            int dLen,
            int strand,
            int candRefPos,
            int[] mtf,
            int lastRef,
            int lastTgt,
            MosaicConfig cfg,
            boolean clonalGear,
            int[] scratchSnpBuf,
            CandidateChain best) {
        if (candRefPos < 0 || candRefPos + MIN_MATCH_LEN > dLen || pos + MIN_MATCH_LEN > tgtLimit) {
            return best;
        }
        // Quick 8-base check so non-collinear breaks bail in 2 ns
        int exactLead =
                Packed2BitSequence.countExactMatch(
                        panelIdx.bank, dStart, dLen, strand, candRefPos, target, pos, 12);
        if (exactLead < 8) {
            return best;
        }
        CandidateChain c =
                buildChain(
                        1,
                        dId,
                        strand,
                        candRefPos,
                        panelIdx.bank,
                        dStart,
                        dLen,
                        target,
                        pos,
                        tgtLimit,
                        mtf,
                        lastRef,
                        lastTgt,
                        cfg,
                        clonalGear,
                        scratchSnpBuf);
        if (c != null && (best == null || c.mdlNetBits > best.mdlNetBits)) {
            return c;
        }
        return best;
    }

    private static CandidateChain probePanelBucket(
            int pos,
            int tgtLimit,
            int syncOff,
            int seedK,
            int qCanonStrand,
            short qSig,
            int bStart,
            int bEnd,
            boolean secondary,
            Packed2BitSequence target,
            EimdlSyncmerIndex panelIdx,
            int maxDonorExclusive,
            int[] donorRank,
            int[] mtf,
            int[] lastRefEnd,
            int[] lastTgtEnd,
            MosaicConfig cfg,
            boolean clonalGear,
            int[] scratchSnpBuf,
            int[] bucketScratch,
            CandidateChain best) {
        bucketScratch[MAX_BUCKET_PROBES * 4] = 1; // conservative default unless we prove zero seed match
        if (pos + MIN_MATCH_LEN > tgtLimit || pos + syncOff + seedK > tgtLimit) {
            return best;
        }
        // Single-pass zero-alloc collection: [0..48)=candIdx, [48..96)=candScoreRank, [96..144)=dupIdx, [144..192)=dupRank
        final int[] candIdx = bucketScratch;
        int candCount = 0;
        int dupCount = 0;
        long seenDonorBitsLo = 0L;
        long seenDonorBitsHi = 0L;

        for (int idx = bStart; idx < bEnd; idx++) {
            short eSig = secondary ? panelIdx.secEntrySig(idx) : panelIdx.entrySig(idx);
            if (eSig != qSig) continue;
            int packed =
                    secondary
                            ? panelIdx.secEntryPackedPosAndStrand(idx)
                            : panelIdx.entryPackedPosAndStrand(idx);
            int globalPos = packed & 0x7FFFFFFF;
            int dId = panelIdx.donorOfGlobalPos(globalPos);
            if (dId >= maxDonorExclusive) break; // donors are appended in ascending ID order!

            boolean alreadySeen =
                    (dId < 64)
                            ? ((seenDonorBitsLo & (1L << dId)) != 0L)
                            : (dId < 128 && (seenDonorBitsHi & (1L << (dId - 64))) != 0L);
            if (alreadySeen) {
                if (dupCount < MAX_BUCKET_PROBES) {
                    candIdx[MAX_BUCKET_PROBES * 2 + dupCount] = idx;
                    candIdx[MAX_BUCKET_PROBES * 3 + dupCount] =
                            (donorRank != null ? donorRank[dId] : dId) + 1000;
                    dupCount++;
                }
                continue;
            }
            if (dId < 64) seenDonorBitsLo |= (1L << dId);
            else if (dId < 128) seenDonorBitsHi |= (1L << (dId - 64));

            int r = donorRank != null ? donorRank[dId] : dId;
            if (((dId << 1) == mtf[0]) || (((dId << 1) | 1) == mtf[0])) {
                r = -2; // sticky active donor gets top priority
            }
            if (candCount < MAX_BUCKET_PROBES) {
                candIdx[candCount] = idx;
                candIdx[MAX_BUCKET_PROBES + candCount] = r;
                candCount++;
            } else {
                int worstSlot = 0;
                for (int u = 1; u < candCount; u++) {
                    if (candIdx[MAX_BUCKET_PROBES + u] > candIdx[MAX_BUCKET_PROBES + worstSlot]) {
                        worstSlot = u;
                    }
                }
                if (r < candIdx[MAX_BUCKET_PROBES + worstSlot]) {
                    candIdx[worstSlot] = idx;
                    candIdx[MAX_BUCKET_PROBES + worstSlot] = r;
                }
            }
        }

        // Fill remaining candidate slots from single-pass paralog backup buffer without re-scanning CSR
        for (int p = 0; p < dupCount && candCount < MAX_BUCKET_PROBES; p++) {
            candIdx[candCount] = candIdx[MAX_BUCKET_PROBES * 2 + p];
            candIdx[MAX_BUCKET_PROBES + candCount] = candIdx[MAX_BUCKET_PROBES * 3 + p];
            candCount++;
        }

        if (candCount == 0) {
            bucketScratch[MAX_BUCKET_PROBES * 4] = 0;
            return best;
        }

        // Insertion-sort small candidate list by priority rank so twin donor evaluates FIRST!
        for (int i = 1; i < candCount; i++) {
            int kIdx = candIdx[i];
            int kRank = candIdx[MAX_BUCKET_PROBES + i];
            int j = i - 1;
            while (j >= 0 && candIdx[MAX_BUCKET_PROBES + j] > kRank) {
                candIdx[j + 1] = candIdx[j];
                candIdx[MAX_BUCKET_PROBES + j + 1] = candIdx[MAX_BUCKET_PROBES + j];
                j--;
            }
            candIdx[j + 1] = kIdx;
            candIdx[MAX_BUCKET_PROBES + j + 1] = kRank;
        }

        boolean anySeedVerified = false;
        for (int u = 0; u < candCount; u++) {
            int idx = candIdx[u];
            int packed =
                    secondary
                            ? panelIdx.secEntryPackedPosAndStrand(idx)
                            : panelIdx.entryPackedPosAndStrand(idx);
            int globalPos = packed & 0x7FFFFFFF;
            int hitCanonStrand = packed >>> 31;
            int dId = panelIdx.donorOfGlobalPos(globalPos);
            int dStart = panelIdx.donorOffsets[dId];
            int dLen = panelIdx.donorOffsets[dId + 1] - dStart;
            int relFwdPos = globalPos - dStart;

            int strand = qCanonStrand ^ hitCanonStrand;
            int refSyncPos = (strand == 0) ? relFwdPos : (dLen - seedK - relFwdPos);
            int candRefStart = refSyncPos - syncOff;
            if (candRefStart < 0 || candRefStart + MIN_MATCH_LEN > dLen) continue;

            // Verify exact k-mer seed match via SWAR word32
            if (Packed2BitSequence.countExactMatch(
                            panelIdx.bank, dStart, dLen, strand, refSyncPos, target, pos + syncOff, seedK)
                    < seedK) {
                continue;
            }
            anySeedVerified = true;

            // If syncOff > 0 and base at pos itself doesn't match (and !allowApprox or >1 mismatch), skip until pos+adj
            if (syncOff > 0
                    && panelIdx.bank.baseAtStrand(dStart, dLen, strand, candRefStart)
                            != target.baseAt(pos)) {
                continue;
            }

            int ds = (dId << 1) | strand;
            CandidateChain c =
                    buildChain(
                            1,
                            dId,
                            strand,
                            candRefStart,
                            panelIdx.bank,
                            dStart,
                            dLen,
                            target,
                            pos,
                            tgtLimit,
                            mtf,
                            lastRefEnd[ds],
                            lastTgtEnd[ds],
                            cfg,
                            clonalGear,
                            scratchSnpBuf);
            if (c != null && (best == null || c.mdlNetBits > best.mdlNetBits)) {
                best = c;
                // Early exit when top-ranked twin yields a long, ultra-clean chromosomal block
                if (best.totalTgtSpan >= 8192
                        && (clonalGear
                                ? (best.totalSnps * 2048 <= best.totalTgtSpan
                                        && best.numSubblocks == 1)
                                : (best.totalSnps * 128 <= best.totalTgtSpan))) {
                    return best;
                }
            }
        }
        if (!anySeedVerified) {
            bucketScratch[MAX_BUCKET_PROBES * 4] = 0;
        }
        return best;
    }

    /**
     * Extends an ungapped subblock with MDL-optimal tail trimming and chains syntenic micro-indel
     * splices when profitable.
     */
    private static CandidateChain buildChain(
            int kind,
            int donorId,
            int strand,
            int startRefPos,
            Packed2BitSequence refBank,
            int refOffset,
            int refLen,
            Packed2BitSequence target,
            int startTgtPos,
            int tgtLimit,
            int[] mtf,
            int lastRef,
            int lastTgt,
            MosaicConfig cfg,
            boolean clonalGear,
            int[] scratchSnpBuf) {
        Subblock first =
                extendSubblockMdl(
                        refBank,
                        refOffset,
                        refLen,
                        strand,
                        startRefPos,
                        target,
                        startTgtPos,
                        tgtLimit,
                        MIN_MATCH_LEN,
                        cfg.allowApprox,
                        clonalGear,
                        scratchSnpBuf);
        if (first == null) {
            return null;
        }

        ArrayList<Subblock> subs = null;
        int curRef = startRefPos + first.len;
        int curTgt = startTgtPos + first.len;
        int totalSnps = first.snpCount;
        int totalNovelSplice = 0;

        if (cfg.allowIndelSplice && !clonalGear) {
            int chains = 0;
            Subblock prev = first;
            int prevRefStart = startRefPos;
            int prevTgtStart = startTgtPos;
            while (chains < MAX_SPLICE_CHAIN
                    && curTgt + MIN_SPLICE_ANCHOR <= tgtLimit
                    && curRef + MIN_SPLICE_ANCHOR <= refLen) {
                int[] splice =
                        tryFindSyntenicSplice(
                                refBank, refOffset, refLen, strand, curRef, target, curTgt, tgtLimit);
                Subblock trimmedPrev = null;
                int trimmedSnpsRemoved = 0;
                if (splice == null && prev.snpCount > 0) {
                    int minPrev = (subs == null) ? MIN_MATCH_LEN : MIN_SPLICE_ANCHOR;
                    int trimLen = prev.snpOffsets[prev.snpCount - 1];
                    int newSnpCnt = prev.snpCount - 1;
                    while (newSnpCnt > 0 && prev.snpOffsets[newSnpCnt - 1] >= trimLen - 1) {
                        trimLen = prev.snpOffsets[newSnpCnt - 1];
                        newSnpCnt--;
                    }
                    if (trimLen >= minPrev && (prev.len - trimLen) <= 16) {
                        int candRef = prevRefStart + trimLen;
                        int candTgt = prevTgtStart + trimLen;
                        int[] trimSplice =
                                tryFindSyntenicSplice(
                                        refBank, refOffset, refLen, strand, candRef, target, candTgt, tgtLimit);
                        if (trimSplice != null) {
                            splice = trimSplice;
                            trimmedSnpsRemoved = prev.snpCount - newSnpCnt;
                            trimmedPrev =
                                    new Subblock(
                                            trimLen,
                                            newSnpCnt == 0
                                                    ? new int[0]
                                                    : Arrays.copyOf(prev.snpOffsets, newSnpCnt),
                                            newSnpCnt,
                                            0,
                                            0);
                        }
                    }
                }
                if (splice == null) break;
                int effEndRef = (trimmedPrev != null) ? (prevRefStart + trimmedPrev.len) : curRef;
                int effEndTgt = (trimmedPrev != null) ? (prevTgtStart + trimmedPrev.len) : curTgt;
                int dTgt = splice[0];
                int dRef = splice[1];
                Subblock next =
                        extendSubblockMdl(
                                refBank,
                                refOffset,
                                refLen,
                                strand,
                                effEndRef + dRef,
                                target,
                                effEndTgt + dTgt,
                                tgtLimit,
                                MIN_SPLICE_ANCHOR,
                                cfg.allowApprox,
                                false,
                                scratchSnpBuf);
                if (next == null) break;

                // Exact MDL profitability check: splice cost vs sticky re-entry (for >=24 bp) or literals (<24 bp)
                double spliceMetaBits =
                        4.2
                                + (dTgt > 0 ? (1.0 + approxCountBits(dTgt - 1)) : 1.0)
                                + (dRef > 0 ? ((dTgt > 0 ? 1.0 : 0.0) + approxCountBits(dRef - 1)) : 0.0);
                double nextBodyBits =
                        approxCountBits(next.len - MIN_SPLICE_ANCHOR)
                                + estimateEditBits(next.len, next.snpCount, cfg);
                double altCostIfBroken =
                        (next.len < MIN_MATCH_LEN)
                                ? (next.len * 1.865)
                                : (16.5
                                        + (dTgt > 0 ? (4.5 + approxCountBits(dTgt - 1)) : 0.0)
                                        + approxCountBits(next.len - MIN_MATCH_LEN)
                                        + estimateEditBits(next.len, next.snpCount, cfg));
                double subGain = altCostIfBroken - nextBodyBits;
                if (subGain <= spliceMetaBits + 2.0) {
                    break;
                }

                if (trimmedPrev != null) {
                    prev = trimmedPrev;
                    curRef = effEndRef;
                    curTgt = effEndTgt;
                    totalSnps -= trimmedSnpsRemoved;
                }
                if (subs == null) {
                    subs = new ArrayList<>(4);
                }
                subs.add(new Subblock(prev.len, prev.snpOffsets, prev.snpCount, dTgt, dRef));
                prevRefStart = curRef + dRef;
                prevTgtStart = curTgt + dTgt;
                prev = next;
                curTgt = prevTgtStart + next.len;
                curRef = prevRefStart + next.len;
                totalSnps += next.snpCount;
                totalNovelSplice += dTgt;
                chains++;
            }
            if (subs != null) {
                subs.add(prev);
            }
        }

        Subblock[] arr =
                (subs == null) ? new Subblock[] {first} : subs.toArray(new Subblock[0]);

        int totalTgtSpan = curTgt - startTgtPos;
        int totalRefSpan = curRef - startRefPos;
        int matchedBases = totalTgtSpan - totalNovelSplice;

        double hdrBits =
                estimateHeaderBits(kind, donorId, strand, startRefPos, startTgtPos, mtf, lastRef, lastTgt);
        double bodyBits = 0.0;
        for (int b = 0; b < arr.length; b++) {
            Subblock sb = arr[b];
            int minL = (b == 0) ? MIN_MATCH_LEN : MIN_SPLICE_ANCHOR;
            bodyBits +=
                approxCountBits(sb.len - minL)
                        + estimateEditBits(sb.len, sb.snpCount, cfg)
                        + (cfg.allowIndelSplice ? (sb.snpCount > 0 ? 0.35 : 0.03) : 0.0);
            if (b + 1 < arr.length) {
                bodyBits +=
                        3.5
                                + (sb.dTgtAfter > 0 ? approxCountBits(sb.dTgtAfter - 1) : 0.0)
                                + (sb.dRefAfter > 0 ? approxCountBits(sb.dRefAfter - 1) : 0.0);
            }
        }

        double mdlNetBits = matchedBases * 1.865 - (hdrBits + bodyBits);
        if (mdlNetBits < 6.0) {
            return null;
        }

        CandidateChain chain = new CandidateChain();
        chain.kind = kind;
        chain.donorId = donorId;
        chain.strand = strand;
        chain.startRefPos = startRefPos;
        chain.totalTgtSpan = totalTgtSpan;
        chain.totalRefSpan = totalRefSpan;
        chain.totalSnps = totalSnps;
        chain.totalNovelSpliceBases = totalNovelSplice;
        chain.mdlNetBits = mdlNetBits;
        chain.subblocks = arr;
        chain.numSubblocks = arr.length;
        return chain;
    }

    private static Subblock extendSubblockMdl(
            Packed2BitSequence refBank,
            int refOffset,
            int refLen,
            int strand,
            int refStart,
            Packed2BitSequence target,
            int tgtStart,
            int tgtLimit,
            int minLen,
            boolean allowApprox,
            boolean clonalGear,
            int[] scratchSnpBuf) {
        int maxLen = Math.min(refLen - refStart, tgtLimit - tgtStart);
        if (maxLen < minLen) {
            return null;
        }
        if (!allowApprox) {
            int exact =
                    Packed2BitSequence.countExactMatch(
                            refBank, refOffset, refLen, strand, refStart, target, tgtStart, maxLen);
            return exact >= minLen ? new Subblock(exact, new int[0], 0, 0, 0) : null;
        }

        int k = 0;
        int snpCount = 0;
        int[] snpBuf = scratchSnpBuf;
        int consecMm = 0;
        int windowMm = 0; // sliding decay of local mismatch burst

        final double baseGain = clonalGear ? 0.008 : 1.865;
        final double snpPenalty = clonalGear ? 11.2 : 7.8;
        final double xDrop = clonalGear ? -16.0 : -22.0;

        double curScore = 0.0;
        double bestScore = -1e9;
        int bestLen = 0;
        int bestSnpCount = 0;

        while (k < maxLen) {
            int exactRun =
                    Packed2BitSequence.countExactMatch(
                            refBank,
                            refOffset,
                            refLen,
                            strand,
                            refStart + k,
                            target,
                            tgtStart + k,
                            maxLen - k);
            if (exactRun > 0) {
                k += exactRun;
                consecMm = 0;
                windowMm =
                        clonalGear
                                ? Math.max(0, windowMm - exactRun)
                                : Math.max(0, windowMm - (exactRun >>> 2));
                curScore += exactRun * baseGain;
                if (k >= minLen && curScore > bestScore) {
                    bestScore = curScore;
                    bestLen = k;
                    bestSnpCount = snpCount;
                }
                if (k >= maxLen) break;
            }

            // Mismatch at k
            if (clonalGear && k < 16) {
                break;
            }
            consecMm++;
            if (clonalGear) {
                if (consecMm >= 2 || windowMm > 0) {
                    break;
                }
                windowMm = 512;
            } else {
                windowMm += 4;
                if (consecMm >= 4 || windowMm >= 28 || (k >= 32 && snpCount * 7 > k)) {
                    break;
                }
            }
            if (curScore - bestScore < xDrop && bestLen >= minLen) {
                break;
            }
            if (snpCount == snpBuf.length) {
                snpBuf = Arrays.copyOf(snpBuf, snpBuf.length << 1);
            }
            snpBuf[snpCount++] = k;
            curScore -= snpPenalty;
            k++;
        }

        if (bestLen < minLen) {
            return null;
        }
        // Guarantee final base of subblock is an exact match (and >= 128 exact tail in clonalGear)
        int reqTail = clonalGear ? 128 : 1;
        while (bestSnpCount > 0 && (bestLen - 1 - snpBuf[bestSnpCount - 1]) < (reqTail - 1)) {
            bestLen = snpBuf[bestSnpCount - 1];
            bestSnpCount--;
        }
        while (bestSnpCount > 0 && snpBuf[bestSnpCount - 1] >= bestLen - 1) {
            bestLen = snpBuf[bestSnpCount - 1];
            bestSnpCount--;
        }
        if (bestLen < minLen) {
            return null;
        }
        int[] snps = bestSnpCount == 0 ? new int[0] : Arrays.copyOf(snpBuf, bestSnpCount);
        return new Subblock(bestLen, snps, bestSnpCount, 0, 0);
    }

    /**
     * Ultra-fast SWAR 16-mer anchor re-locking for syntenic micro-indels on the same {@code (donor,
     * strand)}:
     *
     * <ol>
     *   <li>Pure deletion in target ({@code dTgt = 0, 1 <= dRef <= 192})</li>
     *   <li>Pure insertion in target ({@code 1 <= dTgt <= 48, dRef = 0})</li>
     *   <li>Compound frameshift / hotspot ({@code 1 <= dTgt <= 16, 1 <= dRef <= 32})</li>
     * </ol>
     */
    private static int[] tryFindSyntenicSplice(
            Packed2BitSequence refBank,
            int refOffset,
            int refLen,
            int strand,
            int endRef,
            Packed2BitSequence target,
            int endTgt,
            int tgtLimit) {
        int remTgt = tgtLimit - endTgt;
        int remRef = refLen - endRef;
        if (remTgt < MIN_SPLICE_ANCHOR || remRef < MIN_SPLICE_ANCHOR) {
            return null;
        }

        final int anchorMask = 0xFFFFFF; // 12-base (24-bit) syntenic anchor

        // 1. Pure target deletion: dTgt = 0, dRef in [1 .. 192] via SWAR 64-bit word shifts
        int maxDel = Math.min(192, remRef - MIN_SPLICE_ANCHOR);
        long tgtWord0 = target.word32Forward(endTgt);
        int tgtAnchor12 = ((int) tgtWord0) & anchorMask;
        if (maxDel >= 1) {
            for (int baseRef = 1; baseRef <= maxDel; baseRef += 21) {
                int span = Math.min(21, maxDel - baseRef + 1);
                long rChunk = refBank.word32Strand(refOffset, refLen, strand, endRef + baseRef);
                for (int j = 0; j < span; j++) {
                    if ((((int) rChunk) & anchorMask) == tgtAnchor12) {
                        return new int[] {0, baseRef + j};
                    }
                    rChunk >>>= 2;
                }
            }
        }

        // 2. Pure target insertion: dRef = 0, dTgt in [1 .. 12]
        int maxIns = Math.min(12, remTgt - MIN_SPLICE_ANCHOR);
        int maxCT = Math.min(6, remTgt - MIN_SPLICE_ANCHOR);
        int maxCR = Math.min(12, remRef - MIN_SPLICE_ANCHOR);
        if (maxIns < 1 && (maxCT < 1 || maxCR < 1)) {
            return null;
        }
        long refWord0 = refBank.word32Strand(refOffset, refLen, strand, endRef);
        if (maxIns >= 1) {
            int refAnchor12 = ((int) refWord0) & anchorMask;
            long tShifted = tgtWord0 >>> 2;
            for (int dTgt = 1; dTgt <= maxIns; dTgt++) {
                if ((((int) tShifted) & anchorMask) == refAnchor12) {
                    return new int[] {dTgt, 0};
                }
                tShifted >>>= 2;
            }
        }

        // 3. Compound micro-indel / substitution cluster: dTgt in [1..6], dRef in [1..12]
        if (maxCT >= 1 && maxCR >= 1) {
            long tWord = tgtWord0 >>> 2;
            for (int dt = 1; dt <= maxCT; dt++) {
                int t12 = ((int) tWord) & anchorMask;
                long rWord = refWord0 >>> 2;
                for (int dr = 1; dr <= maxCR; dr++) {
                    if ((((int) rWord) & anchorMask) == t12) {
                        return new int[] {dt, dr};
                    }
                    rWord >>>= 2;
                }
                tWord >>>= 2;
            }
        }
        return null;
    }

    private static double estimateHeaderBits(
            int kind,
            int donorId,
            int strand,
            int startRefPos,
            int startTgtPos,
            int[] mtf,
            int lastRef,
            int lastTgt) {
        if (kind == 2) {
            int dist = (strand == 0) ? (startTgtPos - startRefPos) : (startRefPos + 1);
            return 4.0 + approxCountBits(Math.max(0, dist));
        }
        int ds = (donorId << 1) | strand;
        double donorBits =
                (ds == mtf[0])
                        ? 1.2
                        : ((ds == mtf[1] || ds == mtf[2] || ds == mtf[3]) ? 4.0 : 7.5);
        double posBits;
        if (lastRef >= 0) {
            int d0 = Math.abs(startRefPos - lastRef);
            int d1 = Math.abs(startRefPos - (lastRef + (startTgtPos - lastTgt)));
            int minD = Math.min(d0, d1);
            posBits = minD <= 65536 ? (2.5 + approxCountBits(minD << 1)) : 23.0;
        } else {
            posBits = 2.0 + approxCountBits(startRefPos);
        }
        return 2.0 + donorBits + posBits;
    }

    private static double estimateEditBits(int len, int snpCount, MosaicConfig cfg) {
        if (snpCount == 0) {
            return 1.0; // single hasSnps=0 bit!
        }
        double avgGap = (double) len / (snpCount + 1);
        double gapBits = 1.8 + Math.max(1.5, Math.log(avgGap + 1.0) / Math.log(2.0));
        double subBits = cfg.allowKimuraV4 ? 1.32 : 2.00;
        return 1.5 + snpCount * (gapBits + subBits + 0.8);
    }

    private static double approxCountBits(int value) {
        if (value <= 1) return 2.5;
        int b = 31 - Integer.numberOfLeadingZeros(value);
        return 3.5 + b;
    }

    // =========================================================================
    // Symmetric Entropy Encoding & Decoding Helpers
    // =========================================================================

    private static void encodeSubblocks(
            BitRangeEnc enc,
            MosaicModels models,
            MosaicConfig cfg,
            Packed2BitSequence refBank,
            int refOffset,
            int refLen,
            int strand,
            int startRefPos,
            Packed2BitSequence target,
            int startTgtPos,
            Subblock[] subblocks) {
        int curRef = startRefPos;
        int curTgt = startTgtPos;
        int emaGap = 128;

        for (int b = 0; b < subblocks.length; b++) {
            Subblock sb = subblocks[b];
            int minL = (b == 0) ? MIN_MATCH_LEN : MIN_SPLICE_ANCHOR;
            models.writeCount(enc, models.matchLenSlots, b == 0 ? 0 : 1, 1, sb.len - minL);

            if (cfg.allowApprox) {
                int lenBin = sb.len < 64 ? 0 : (sb.len < 512 ? 1 : (sb.len < 4096 ? 2 : 3));
                int hasSnps = sb.snpCount > 0 ? 1 : 0;
                enc.encodeBit(models.hasSnps, ((b > 0 ? 1 : 0) << 2) | lenBin, hasSnps);

                if (hasSnps == 1) {
                    int localOff = 0;
                    int prevMod3Zero = 1;
                    for (int s = 0; s < sb.snpCount; s++) {
                        int snpOff = sb.snpOffsets[s];
                        int gap = snpOff - localOff;
                        int maxGap = (sb.len - 1) - localOff;

                        if (s > 0) {
                            // localOff < sb.len - 1 is guaranteed here because maxGap >= gap >= 0
                            // and snpOff <= sb.len - 2 (since final base of subblock is exact!)
                            int remBin = maxGap < 16 ? 0 : (maxGap < 128 ? 1 : 2);
                            int reg = gapRegime(emaGap);
                            enc.encodeBit(models.moreSnps, (reg << 2) | remBin, 1);
                        }

                        int reg = gapRegime(emaGap);
                        encodeSnpGap(
                                enc,
                                models,
                                cfg.allowWobbleMod3,
                                gap,
                                maxGap,
                                reg,
                                prevMod3Zero);
                        prevMod3Zero = ((gap + 1) % 3 == 0) ? 1 : 0;
                        emaGap = (emaGap * 3 + Math.min(2048, gap)) >>> 2;

                        int mRef = curRef + snpOff;
                        int mTgt = curTgt + snpOff;
                        int refB = refBank.baseAtStrand(refOffset, refLen, strand, mRef);
                        int tgtB = target.baseAt(mTgt);
                        if (cfg.allowKimuraV4) {
                            int prevRefB =
                                    (mRef > 0)
                                            ? refBank.baseAtStrand(refOffset, refLen, strand, mRef - 1)
                                            : 0;
                            int nextRefB =
                                    (mRef + 1 < refLen)
                                            ? refBank.baseAtStrand(refOffset, refLen, strand, mRef + 1)
                                            : 0;
                            int subCtx =
                                    (refB << 6)
                                            | (prevRefB << 4)
                                            | (nextRefB << 2)
                                            | ((reg <= 1 ? 1 : 0) << 1)
                                            | (cfg.allowWobbleMod3 ? prevMod3Zero : 0);
                            int delta = tgtB ^ refB; // guaranteed in {1, 2, 3}
                            if (delta == 2) {
                                enc.encodeBit(models.kimuraTi, subCtx, 1);
                            } else {
                                enc.encodeBit(models.kimuraTi, subCtx, 0);
                                enc.encodeBit(models.kimuraTv, subCtx, delta == 3 ? 1 : 0);
                            }
                        } else {
                            int prevTgtB = (mTgt > 0) ? target.baseAt(mTgt - 1) : 0;
                            models.writeFixedTree(
                                    enc, models.ablateBase4[prevTgtB], tgtB, 2);
                        }

                        localOff = snpOff + 1;
                    }
                    // Emit end-of-SNPs bit if another SNP could still legally fit inside [localOff .. sb.len - 1)
                    if (localOff < sb.len - 1) {
                        int maxGap = (sb.len - 1) - localOff;
                        int remBin = maxGap < 16 ? 0 : (maxGap < 128 ? 1 : 2);
                        int reg = gapRegime(emaGap);
                        enc.encodeBit(models.moreSnps, (reg << 2) | remBin, 0);
                    }
                }
            }

            int hasSnpsFlag = (cfg.allowApprox && sb.snpCount > 0) ? 1 : 0;
            curRef += sb.len;
            curTgt += sb.len;

            boolean canSplice =
                    cfg.allowIndelSplice
                            && (curTgt + MIN_SPLICE_ANCHOR <= target.length)
                            && (curRef + MIN_SPLICE_ANCHOR <= refLen);
            if (canSplice) {
                boolean hasSplice = (b + 1 < subblocks.length);
                int spCtx = ((b > 0 ? 1 : 0) << 2) | (hasSnpsFlag << 1) | (sb.len >= 256 ? 1 : 0);
                enc.encodeBit(models.hasSplice, spCtx, hasSplice ? 1 : 0);
                if (hasSplice) {
                    int dTgt = sb.dTgtAfter;
                    int dRef = sb.dRefAfter;
                    if (dTgt == 0) {
                        enc.encodeBit(models.spliceIsPureDel, 0, 1);
                        models.writeCount(enc, models.spliceDelSlots, 0, 7, dRef - 1);
                    } else {
                        enc.encodeBit(models.spliceIsPureDel, 0, 0);
                        if (dRef == 0) {
                            enc.encodeBit(models.spliceIsPureIns, 0, 1);
                            models.writeCount(enc, models.spliceInsSlots, 0, 6, dTgt - 1);
                        } else {
                            enc.encodeBit(models.spliceIsPureIns, 0, 0);
                            models.writeCount(enc, models.spliceInsSlots, 0, 6, dTgt - 1);
                            models.writeCount(enc, models.spliceDelSlots, 0, 7, dRef - 1);
                        }
                    }
                    curTgt += dTgt;
                    curRef += dRef;
                }
            }
        }
    }

    private static int[] decodeSubblocksFromPanel(
            BitRangeDec dec,
            MosaicModels models,
            MosaicConfig cfg,
            Packed2BitSequence refBank,
            int refOffset,
            int refLen,
            int strand,
            int startRefPos,
            byte[] out,
            int startTgtPos,
            int targetLen,
            byte[] noveltyBases,
            int novCursor,
            int noveltyLen,
            int donorId,
            MosaicObserver obs) {
        return decodeSubblocks(
                dec, models, cfg, refBank, refOffset, refLen, false, strand,
                startRefPos, out, startTgtPos, targetLen, noveltyBases, novCursor, noveltyLen,
                donorId, obs);
    }

    private static int[] decodeSubblocksFromSelf(
            BitRangeDec dec,
            MosaicModels models,
            MosaicConfig cfg,
            byte[] out,
            int selfRefLen,
            int strand,
            int startRefPos,
            int startTgtPos,
            int targetLen,
            byte[] noveltyBases,
            int novCursor,
            int noveltyLen,
            MosaicObserver obs) {
        return decodeSubblocks(
                dec, models, cfg, null, 0, selfRefLen, true, strand,
                startRefPos, out, startTgtPos, targetLen, noveltyBases, novCursor, noveltyLen,
                -1, obs);
    }

    private static int[] decodeSubblocks(
            BitRangeDec dec,
            MosaicModels models,
            MosaicConfig cfg,
            Packed2BitSequence refBank,
            int refOffset,
            int refLen,
            boolean isSelf,
            int strand,
            int startRefPos,
            byte[] out,
            int startTgtPos,
            int targetLen,
            byte[] noveltyBases,
            int novCursor,
            int noveltyLen,
            int donorId,
            MosaicObserver obs) {
        int curRef = startRefPos;
        int curTgt = startTgtPos;
        int emaGap = 128;
        int b = 0;

        while (true) {
            int minL = (b == 0) ? MIN_MATCH_LEN : MIN_SPLICE_ANCHOR;
            int subLen = models.readCount(dec, models.matchLenSlots, b == 0 ? 0 : 1, 1) + minL;
            if (curRef < 0 || curRef + subLen > refLen || curTgt + subLen > targetLen) {
                throw new IllegalStateException("Subblock exceeds donor or target bounds");
            }

            if (out != null) {
                if (!isSelf) {
                    copyRefToOut(refBank, refOffset, refLen, strand, curRef, out, curTgt, subLen);
                } else if (strand == 0) {
                    for (int k = 0; k < subLen; k++) {
                        out[curTgt + k] = out[curRef + k];
                    }
                } else {
                    int fwdTop = refLen - 1 - curRef;
                    for (int k = 0; k < subLen; k++) {
                        out[curTgt + k] = (byte) ((out[fwdTop - k] ^ 3) & 3);
                    }
                }
            }

            int spanIdx = -1;
            if (obs != null) {
                if (isSelf) obs.selfMatchBases += subLen;
                else obs.panelMatchBases += subLen;
                if (obs.collectSpans) {
                    spanIdx = obs.spans.size();
                    obs.spans.add(
                            new SpanEvent(isSelf ? 2 : 1, strand, false, donorId, curRef, curTgt, subLen));
                }
            }

            int hasSnpsFlag = 0;
            long subSnps = 0L;
            if (cfg.allowApprox) {
                int lenBin = subLen < 64 ? 0 : (subLen < 512 ? 1 : (subLen < 4096 ? 2 : 3));
                int hasSnps = dec.decodeBit(models.hasSnps, ((b > 0 ? 1 : 0) << 2) | lenBin);
                hasSnpsFlag = hasSnps;
                if (hasSnps == 1) {
                    int localOff = 0;
                    int prevMod3Zero = 1;
                    int s = 0;
                    while (true) {
                        int maxGap = (subLen - 1) - localOff;
                        if (s > 0) {
                            if (maxGap < 0) break;
                            int remBin = maxGap < 16 ? 0 : (maxGap < 128 ? 1 : 2);
                            int reg = gapRegime(emaGap);
                            if (dec.decodeBit(models.moreSnps, (reg << 2) | remBin) == 0) {
                                break;
                            }
                        }
                        int reg = gapRegime(emaGap);
                        int gap =
                                decodeSnpGap(
                                        dec,
                                        models,
                                        cfg.allowWobbleMod3,
                                        maxGap,
                                        reg,
                                        prevMod3Zero);
                        int snpOff = localOff + gap;
                        if (snpOff < 0 || snpOff >= subLen) {
                            throw new IllegalStateException("Decoded SNP offset out of subblock bounds");
                        }
                        prevMod3Zero = ((gap + 1) % 3 == 0) ? 1 : 0;
                        emaGap = (emaGap * 3 + Math.min(2048, gap)) >>> 2;

                        int mRef = curRef + snpOff;
                        int mTgt = curTgt + snpOff;
                        int refB =
                                (out != null)
                                        ? (out[mTgt] & 3)
                                        : refBank.baseAtStrand(refOffset, refLen, strand, mRef);
                        int tgtB;
                        if (cfg.allowKimuraV4) {
                            int prevRefB = 0;
                            int nextRefB = 0;
                            if (!isSelf) {
                                prevRefB =
                                        (mRef > 0)
                                                ? refBank.baseAtStrand(refOffset, refLen, strand, mRef - 1)
                                                : 0;
                                nextRefB =
                                        (mRef + 1 < refLen)
                                                ? refBank.baseAtStrand(refOffset, refLen, strand, mRef + 1)
                                                : 0;
                            } else {
                                prevRefB = (mRef > 0) ? selfBaseAt(out, refLen, strand, mRef - 1) : 0;
                                nextRefB =
                                        (mRef + 1 < refLen)
                                                ? selfBaseAt(out, refLen, strand, mRef + 1)
                                                : 0;
                            }
                            int subCtx =
                                    (refB << 6)
                                            | (prevRefB << 4)
                                            | (nextRefB << 2)
                                            | ((reg <= 1 ? 1 : 0) << 1)
                                            | (cfg.allowWobbleMod3 ? prevMod3Zero : 0);
                            int delta;
                            if (dec.decodeBit(models.kimuraTi, subCtx) == 1) {
                                delta = 2;
                            } else {
                                delta = dec.decodeBit(models.kimuraTv, subCtx) == 1 ? 3 : 1;
                            }
                            tgtB = refB ^ delta;
                        } else {
                            int prevTgtB = (mTgt > 0 && out != null) ? (out[mTgt - 1] & 3) : 0;
                            tgtB = models.readFixedTree(dec, models.ablateBase4[prevTgtB], 2);
                        }
                        if (out != null) out[mTgt] = (byte) tgtB;
                        if (obs != null) {
                            int delta = (refB ^ tgtB) & 3;
                            int isTi = (delta == 2) ? 1 : 0;
                            subSnps++;
                            obs.snps++;
                            if (isTi == 1) obs.transitions++;
                            else obs.transversions++;
                            if (prevMod3Zero != 0) obs.wobbleMod3Hits++;
                            if (obs.collectVariants && !isSelf) {
                                obs.variants.add(
                                        new VariantRecord(
                                                0,
                                                strand,
                                                isTi,
                                                prevMod3Zero,
                                                donorId,
                                                mRef,
                                                mTgt,
                                                1,
                                                1,
                                                new byte[] {(byte) refB},
                                                new byte[] {(byte) tgtB}));
                            }
                        }
                        localOff = snpOff + 1;
                        s++;
                        if (localOff >= subLen - 1) {
                            break;
                        }
                    }
                }
            }
            if (spanIdx >= 0 && subSnps > 0L) {
                SpanEvent prev = obs.spans.get(spanIdx);
                obs.spans.set(
                        spanIdx,
                        new SpanEvent(
                                prev.kind,
                                prev.strand,
                                true,
                                prev.donorId,
                                prev.refPos,
                                prev.tgtPos,
                                prev.tgtLen));
            }

            curRef += subLen;
            curTgt += subLen;

            boolean canSplice =
                    cfg.allowIndelSplice
                            && (curTgt + MIN_SPLICE_ANCHOR <= targetLen)
                            && (curRef + MIN_SPLICE_ANCHOR <= refLen);
            int spCtx = ((b > 0 ? 1 : 0) << 2) | (hasSnpsFlag << 1) | (subLen >= 256 ? 1 : 0);
            if (!canSplice || dec.decodeBit(models.hasSplice, spCtx) == 0) {
                break;
            }
            int dTgt;
            int dRef;
            if (dec.decodeBit(models.spliceIsPureDel, 0) == 1) {
                dTgt = 0;
                dRef = models.readCount(dec, models.spliceDelSlots, 0, 7) + 1;
            } else if (dec.decodeBit(models.spliceIsPureIns, 0) == 1) {
                dRef = 0;
                dTgt = models.readCount(dec, models.spliceInsSlots, 0, 6) + 1;
            } else {
                dTgt = models.readCount(dec, models.spliceInsSlots, 0, 6) + 1;
                dRef = models.readCount(dec, models.spliceDelSlots, 0, 7) + 1;
            }
            if (curTgt + dTgt > targetLen || novCursor + dTgt > noveltyLen || curRef + dRef > refLen) {
                throw new IllegalStateException("Splice insertion/deletion exceeds bounds");
            }
            int spliceRefPos = curRef;
            int spliceTgtPos = curTgt;
            int novStart = novCursor;
            if (dTgt > 0) {
                if (out != null) System.arraycopy(noveltyBases, novCursor, out, curTgt, dTgt);
                novCursor += dTgt;
                curTgt += dTgt;
            }
            curRef += dRef;
            if (obs != null) {
                obs.splices++;
                obs.insertedBases += dTgt;
                obs.deletedBases += dRef;
                if (obs.collectSpans && dTgt > 0) {
                    obs.spans.add(new SpanEvent(0, 0, true, -1, novStart, spliceTgtPos, dTgt));
                }
                if (obs.collectVariants && !isSelf) {
                    int kind = (dTgt > 0 && dRef == 0) ? 1 : ((dRef > 0 && dTgt == 0) ? 2 : 3);
                    byte[] rb = null;
                    if (dRef > 0 && dRef <= 65536) {
                        rb = new byte[dRef];
                        for (int k = 0; k < dRef; k++) {
                            rb[k] = (byte) refBank.baseAtStrand(refOffset, refLen, strand, spliceRefPos + k);
                        }
                    }
                    byte[] ab =
                            (dTgt > 0 && dTgt <= 65536)
                                    ? Arrays.copyOfRange(noveltyBases, novStart, novStart + dTgt)
                                    : null;
                    obs.variants.add(
                            new VariantRecord(
                                    kind, strand, 0, 0, donorId, spliceRefPos, spliceTgtPos, dRef, dTgt, rb, ab));
                }
            }
            b++;
        }
        return new int[] {curRef, curTgt, novCursor};
    }

    private static int selfBaseAt(byte[] out, int selfRefLen, int strand, int pos) {
        return strand == 0 ? (out[pos] & 3) : ((out[selfRefLen - 1 - pos] ^ 3) & 3);
    }

    private static void copyRefToOut(
            Packed2BitSequence refBank,
            int refOffset,
            int refLen,
            int strand,
            int startRefPos,
            byte[] out,
            int dstPos,
            int len) {
        int k = 0;
        if (strand == 0) {
            int rBase = refOffset + startRefPos;
            while (k + 32 <= len) {
                long w = refBank.word32Forward(rBase + k);
                for (int j = 0; j < 32; j++) {
                    out[dstPos + k + j] = (byte) (w & 3L);
                    w >>>= 2;
                }
                k += 32;
            }
            while (k < len) {
                out[dstPos + k] = (byte) refBank.baseAt(rBase + k);
                k++;
            }
        } else {
            int fwdTop = refOffset + refLen - 1 - startRefPos;
            while (k + 32 <= len) {
                long w = Packed2BitSequence.reverseComplementWord32(
                        refBank.word32Forward(fwdTop - k - 31));
                for (int j = 0; j < 32; j++) {
                    out[dstPos + k + j] = (byte) (w & 3L);
                    w >>>= 2;
                }
                k += 32;
            }
            while (k < len) {
                out[dstPos + k] = (byte) (refBank.baseAt(fwdTop - k) ^ 3);
                k++;
            }
        }
    }

    private static int gapRegime(int emaGap) {
        if (emaGap < 16) return 0;
        if (emaGap < 96) return 1;
        if (emaGap < 512) return 2;
        return 3;
    }

    private static void encodeSnpGap(
            BitRangeEnc enc,
            MosaicModels models,
            boolean wobbleMod3,
            int gap,
            int maxGap,
            int reg,
            int prevMod3Zero) {
        if (!wobbleMod3) {
            models.writeCount(enc, models.gapQuotSlots, reg, 5, gap);
            return;
        }
        int q = gap / 3;
        int remIdx = gap - 3 * q; // 0 -> r=1, 1 -> r=2, 2 -> r=0 (wobble!)
        models.writeCount(enc, models.gapQuotSlots, reg, 5, q);

        int maxRemIdx = Math.min(2, maxGap - 3 * q);
        if (maxRemIdx == 0) {
            return; // deterministic: remIdx must be 0
        }
        int qBin = q == 0 ? 0 : (q < 8 ? 1 : 2);
        int wobCtx = (reg << 3) | (qBin << 1) | prevMod3Zero;
        if (maxRemIdx == 1) {
            enc.encodeBit(models.wobbleRem01, wobCtx, remIdx);
            return;
        }
        if (remIdx == 2) {
            enc.encodeBit(models.wobbleIsThirdBase, wobCtx, 1);
        } else {
            enc.encodeBit(models.wobbleIsThirdBase, wobCtx, 0);
            enc.encodeBit(models.wobbleRem01, wobCtx, remIdx);
        }
    }

    private static int decodeSnpGap(
            BitRangeDec dec,
            MosaicModels models,
            boolean wobbleMod3,
            int maxGap,
            int reg,
            int prevMod3Zero) {
        if (!wobbleMod3) {
            return models.readCount(dec, models.gapQuotSlots, reg, 5);
        }
        int q = models.readCount(dec, models.gapQuotSlots, reg, 5);
        int maxRemIdx = Math.min(2, maxGap - 3 * q);
        if (maxRemIdx < 0) {
            throw new IllegalStateException("Decoded gap quotient exceeds subblock bound");
        }
        if (maxRemIdx == 0) {
            return 3 * q;
        }
        int qBin = q == 0 ? 0 : (q < 8 ? 1 : 2);
        int wobCtx = (reg << 3) | (qBin << 1) | prevMod3Zero;
        if (maxRemIdx == 1) {
            int remIdx = dec.decodeBit(models.wobbleRem01, wobCtx);
            return 3 * q + remIdx;
        }
        if (dec.decodeBit(models.wobbleIsThirdBase, wobCtx) == 1) {
            return 3 * q + 2;
        }
        int remIdx = dec.decodeBit(models.wobbleRem01, wobCtx);
        return 3 * q + remIdx;
    }

    private static void encodeMtfDonorStrand(
            BitRangeEnc enc, MosaicModels models, int[] mtf, int maxDonors, int ds, int afterLit) {
        int totalStates = maxDonors << 1;
        if (totalStates <= 1) {
            return;
        }
        if (ds == mtf[0]) {
            enc.encodeBit(models.mtfHit0, afterLit, 1);
            return;
        }
        enc.encodeBit(models.mtfHit0, afterLit, 0);
        if (totalStates == 2) {
            promoteMtf(mtf, ds);
            return;
        }
        int mtfSlot = -1;
        int maxSlot = Math.min(4, totalStates);
        for (int s = 1; s < maxSlot; s++) {
            if (mtf[s] == ds) {
                mtfSlot = s;
                break;
            }
        }
        if (mtfSlot >= 1) {
            if (totalStates > 4) {
                enc.encodeBit(models.mtfHit13, 0, 1);
            }
            if (maxSlot == 3) {
                enc.encodeBit(models.mtfChoice, 0, mtfSlot - 1);
            } else if (maxSlot >= 4) {
                if (mtfSlot == 1) {
                    enc.encodeBit(models.mtfChoice, 0, 0);
                } else {
                    enc.encodeBit(models.mtfChoice, 0, 1);
                    enc.encodeBit(models.mtfChoice, 1, mtfSlot - 2);
                }
            }
            promoteMtf(mtf, ds);
            return;
        }
        if (totalStates > 4) {
            enc.encodeBit(models.mtfHit13, 0, 0);
        }
        int donorId = ds >>> 1;
        int strand = ds & 1;
        models.writeBoundedInt(enc, models.donorTree, donorId, maxDonors);
        enc.encodeBit(models.panelStrand, 0, strand);
        promoteMtf(mtf, ds);
    }

    private static int decodeMtfDonorStrand(
            BitRangeDec dec, MosaicModels models, int[] mtf, int maxDonors, int afterLit) {
        int totalStates = maxDonors << 1;
        if (totalStates <= 1) {
            return 0;
        }
        if (dec.decodeBit(models.mtfHit0, afterLit) == 1) {
            return mtf[0];
        }
        if (totalStates == 2) {
            int ds = mtf[0] ^ 1;
            promoteMtf(mtf, ds);
            return ds;
        }
        int maxSlot = Math.min(4, totalStates);
        boolean inMtf13 = (totalStates <= 4) || (dec.decodeBit(models.mtfHit13, 0) == 1);
        if (inMtf13) {
            int slot;
            if (maxSlot == 3) {
                slot = 1 + dec.decodeBit(models.mtfChoice, 0);
            } else {
                if (dec.decodeBit(models.mtfChoice, 0) == 0) {
                    slot = 1;
                } else {
                    slot = 2 + dec.decodeBit(models.mtfChoice, 1);
                }
            }
            int ds = mtf[slot];
            promoteMtf(mtf, ds);
            return ds;
        }
        int donorId = models.readBoundedInt(dec, models.donorTree, maxDonors);
        int strand = dec.decodeBit(models.panelStrand, 0);
        int ds = (donorId << 1) | strand;
        promoteMtf(mtf, ds);
        return ds;
    }

    private static void encodePanelRefPos(
            BitRangeEnc enc,
            MosaicModels models,
            int refPos,
            int tgtPos,
            int donorLen,
            int lastRef,
            int lastTgt,
            int isMtf0) {
        if (lastRef >= 0) {
            int pred0 = Math.max(0, Math.min(donorLen - 1, lastRef));
            int pred1 = Math.max(0, Math.min(donorLen - 1, lastRef + (tgtPos - lastTgt)));
            int d0 = refPos - pred0;
            int d1 = refPos - pred1;
            int bestPred;
            int bestDelta;
            int useCol;
            if (Math.abs(d1) <= Math.abs(d0)) {
                bestPred = pred1;
                bestDelta = d1;
                useCol = 1;
            } else {
                bestPred = pred0;
                bestDelta = d0;
                useCol = 0;
            }
            if (Math.abs(bestDelta) <= 65536) {
                enc.encodeBit(models.usePosDelta, isMtf0, 1);
                if (pred0 != pred1) {
                    enc.encodeBit(models.useCollinearPred, isMtf0, useCol);
                }
                int zz = (bestDelta << 1) ^ (bestDelta >> 31);
                models.writeCount(enc, models.posDeltaSlots, useCol, 3, zz);
                bestPred &= 0;
                return;
            }
            enc.encodeBit(models.usePosDelta, isMtf0, 0);
        }
        models.writeCount(enc, models.rawPosSlots, 0, 4, refPos);
    }

    private static int decodePanelRefPos(
            BitRangeDec dec,
            MosaicModels models,
            int tgtPos,
            int donorLen,
            int lastRef,
            int lastTgt,
            int isMtf0) {
        if (lastRef >= 0 && dec.decodeBit(models.usePosDelta, isMtf0) == 1) {
            int pred0 = Math.max(0, Math.min(donorLen - 1, lastRef));
            int pred1 = Math.max(0, Math.min(donorLen - 1, lastRef + (tgtPos - lastTgt)));
            int pred = pred0;
            int useCol = 1;
            if (pred0 != pred1) {
                useCol = dec.decodeBit(models.useCollinearPred, isMtf0);
                pred = useCol == 1 ? pred1 : pred0;
            }
            int zz = models.readCount(dec, models.posDeltaSlots, useCol, 3);
            int delta = (zz >>> 1) ^ -(zz & 1);
            int refPos = pred + delta;
            if (refPos < 0 || refPos >= donorLen) {
                throw new IllegalStateException("Decoded delta position out of donor bounds");
            }
            return refPos;
        }
        int refPos = models.readCount(dec, models.rawPosSlots, 0, 4);
        if (refPos < 0 || refPos >= donorLen) {
            throw new IllegalStateException("Decoded raw position out of donor bounds");
        }
        return refPos;
    }

    private static int[] initMtf(int maxDonors) {
        int states = Math.max(4, maxDonors << 1);
        int[] mtf = new int[4];
        for (int i = 0; i < 4; i++) {
            mtf[i] = i % states;
        }
        return mtf;
    }

    private static void promoteMtf(int[] mtf, int val) {
        if (mtf[0] == val) return;
        int pos = 3;
        for (int i = 1; i < 3; i++) {
            if (mtf[i] == val) {
                pos = i;
                break;
            }
        }
        for (int i = pos; i > 0; i--) {
            mtf[i] = mtf[i - 1];
        }
        mtf[0] = val;
    }

    private static byte[] appendNovelBases(
            byte[] buf, int curLen, Packed2BitSequence seq, int start, int count) {
        if (curLen + count > buf.length) {
            int newCap = Math.max(curLen + count, buf.length + (buf.length >>> 1) + 1024);
            buf = Arrays.copyOf(buf, newCap);
        }
        for (int i = 0; i < count; i++) {
            buf[curLen + i] = (byte) seq.baseAt(start + i);
        }
        return buf;
    }

    private static int findNextSyncmerPos(TargetSyncmers ts, int startIdx, int minPos, int maxPos) {
        for (int i = startIdx; i < ts.count; i++) {
            int p = ts.positions[i];
            if (p > maxPos) return -1;
            if (p >= minPos) return p;
        }
        return -1;
    }

    private record TargetSyncmers(int[] positions, int[] canonKmers, byte[] strands, int count) {}

    private static TargetSyncmers extractTargetSyncmers(Packed2BitSequence seq, boolean primary) {
        int k = primary ? Packed2BitSequence.K_PRIMARY : Packed2BitSequence.K_SECONDARY;
        int s = primary ? Packed2BitSequence.S_PRIMARY : Packed2BitSequence.S_SECONDARY;
        int t = primary ? Packed2BitSequence.T_PRIMARY : Packed2BitSequence.T_SECONDARY;
        int cap = Math.max(32, seq.length / (k - s) + 32);
        int[][] pos = {new int[cap]};
        int[][] km = {new int[cap]};
        byte[][] st = {new byte[cap]};
        int[] cnt = {0};
        seq.forEachCanonicalOpenSyncmer(
                0,
                seq.length,
                k,
                s,
                t,
                (p, canon, strand) -> {
                    int idx = cnt[0];
                    if (idx == pos[0].length) {
                        int ncap = idx + (idx >>> 1) + 64;
                        pos[0] = Arrays.copyOf(pos[0], ncap);
                        km[0] = Arrays.copyOf(km[0], ncap);
                        st[0] = Arrays.copyOf(st[0], ncap);
                    }
                    pos[0][idx] = p;
                    km[0][idx] = (int) canon;
                    st[0][idx] = (byte) strand;
                    cnt[0] = idx + 1;
                });
        return new TargetSyncmers(pos[0], km[0], st[0], cnt[0]);
    }

    private static final class ParsedToken {
        final int kind; // 0 = LIT, 1 = PANEL, 2 = SELF, 3 = TANDEM
        final int litLen;
        final int donorId; // or period for TANDEM
        final int strand;
        final int startRefPos;
        final int totalTgtSpan;
        final int totalRefSpan;
        final Subblock[] subblocks;

        private ParsedToken(
                int kind,
                int litLen,
                int donorId,
                int strand,
                int startRefPos,
                int totalTgtSpan,
                int totalRefSpan,
                Subblock[] subblocks) {
            this.kind = kind;
            this.litLen = litLen;
            this.donorId = donorId;
            this.strand = strand;
            this.startRefPos = startRefPos;
            this.totalTgtSpan = totalTgtSpan;
            this.totalRefSpan = totalRefSpan;
            this.subblocks = subblocks;
        }

        static ParsedToken literal(int len) {
            return new ParsedToken(0, len, 0, 0, 0, len, 0, null);
        }

        static ParsedToken fromChain(int tgtPos, CandidateChain c) {
            tgtPos &= 0;
            return new ParsedToken(
                    c.kind,
                    0,
                    c.donorId,
                    c.strand,
                    c.startRefPos,
                    c.totalTgtSpan,
                    c.totalRefSpan,
                    c.subblocks);
        }
    }

    // =========================================================================
    // Adaptive Probability Tables & Binary Range Coder
    // =========================================================================

    private static final class MosaicModels {
        final int[] isLiteral = probs(4);
        final int[] isPanelMatch = probs(2);
        final int[] isTandem = probs(2);
        final int[] selfStrand = probs(2);
        final int[] panelStrand = probs(2);

        final int[] mtfHit0 = probs(2);
        final int[] mtfHit13 = probs(2);
        final int[] mtfChoice = probs(4);
        final int[] donorTree = probs(4096); // supports up to 4096 donors cleanly

        final int[] usePosDelta = probs(2);
        final int[] useCollinearPred = probs(2);

        final int[] hasSnps = probsInit(8, 600);
        final int[] moreSnps = probs(16);
        final int[] hasSplice = probsInit(8, 96);
        final int[] spliceIsPureDel = probs(2);
        final int[] spliceIsPureIns = probsInit(2, 3800);

        // Law I & II Coupled: Klein V_4 Kimura + CpG trinucleotide + mod-3 wobble contexts (256 contexts)
        final int[] kimuraTi = probsInit(256, 2700); // ~66% initial Kimura transition prior!
        final int[] kimuraTv = probs(256);
        final int[][] ablateBase4 = probs2d(4, 4);

        // Law II: Bijective Codon-Wobble phase models ((reg<<3) | (qBin<<1) | prevWobble < 32)
        final int[] wobbleIsThirdBase = probsInit(32, 1800);
        final int[] wobbleRem01 = probs(32);

        // Variable-length integer slot models (32 slots -> 6-bit tree of size 64, plus isolated low-bit refinements)
        final int[][] litLenSlots = probs2d(1, 64);
        final int[][] matchLenSlots = probs2d(2, 64);
        final int[][] selfDistSlots = probs2d(2, 64);
        final int[][] posDeltaSlots = probs2d(2, 64);
        final int[][] rawPosSlots = probs2d(1, 64);
        final int[][] gapQuotSlots = probs2d(4, 64);
        final int[][] spliceInsSlots = probs2d(1, 64);
        final int[][] spliceDelSlots = probs2d(1, 64);
        final int[][] lowBitProbs = probs2d(8, 128);

        final int[] tandemPeriodTree = probs(64);

        private static int[] probs(int n) {
            int[] a = new int[n];
            Arrays.fill(a, PROB_HALF);
            return a;
        }

        private static int[] probsInit(int n, int init) {
            int[] a = new int[n];
            Arrays.fill(a, init);
            return a;
        }

        private static int[][] probs2d(int rows, int cols) {
            int[][] m = new int[rows][cols];
            for (int r = 0; r < rows; r++) {
                Arrays.fill(m[r], PROB_HALF);
            }
            return m;
        }

        void writeFixedTree(BitRangeEnc enc, int[] tree, int sym, int bits) {
            int node = 1;
            for (int b = bits - 1; b >= 0; b--) {
                int bit = (sym >>> b) & 1;
                enc.encodeBit(tree, node, bit);
                node = (node << 1) | bit;
            }
        }

        int readFixedTree(BitRangeDec dec, int[] tree, int bits) {
            int node = 1;
            int val = 0;
            for (int b = 0; b < bits; b++) {
                int bit = dec.decodeBit(tree, node);
                val = (val << 1) | bit;
                node = (node << 1) | bit;
            }
            return val;
        }

        void writeBoundedInt(BitRangeEnc enc, int[] tree, int val, int maxExclusive) {
            if (maxExclusive <= 1) return;
            int bits = 32 - Integer.numberOfLeadingZeros(maxExclusive - 1);
            if (bits > 11) {
                writeCount(enc, rawPosSlots, 0, 4, val);
                return;
            }
            int node = 1;
            int prefix = 0;
            for (int b = bits - 1; b >= 0; b--) {
                int zeroUpper = (prefix << (b + 1)) | (1 << b);
                if (zeroUpper >= maxExclusive) {
                    // Upper branch >= maxExclusive is impossible; bit 0 is deterministic!
                    prefix = prefix << 1;
                    node = node << 1;
                    continue;
                }
                int bit = (val >>> b) & 1;
                enc.encodeBit(tree, node, bit);
                prefix = (prefix << 1) | bit;
                node = (node << 1) | bit;
            }
        }

        int readBoundedInt(BitRangeDec dec, int[] tree, int maxExclusive) {
            if (maxExclusive <= 1) return 0;
            int bits = 32 - Integer.numberOfLeadingZeros(maxExclusive - 1);
            if (bits > 11) {
                return readCount(dec, rawPosSlots, 0, 4);
            }
            int node = 1;
            int prefix = 0;
            for (int b = bits - 1; b >= 0; b--) {
                int zeroUpper = (prefix << (b + 1)) | (1 << b);
                if (zeroUpper >= maxExclusive) {
                    prefix = prefix << 1;
                    node = node << 1;
                    continue;
                }
                int bit = dec.decodeBit(tree, node);
                prefix = (prefix << 1) | bit;
                node = (node << 1) | bit;
            }
            return prefix;
        }

        void writeCount(BitRangeEnc enc, int[][] slotTable, int ctx, int lowCtx, int value) {
            int[] tree = slotTable[ctx];
            if (value < 4) {
                writeFixedTree(enc, tree, value, 6);
                return;
            }
            int msb = 31 - Integer.numberOfLeadingZeros(value);
            int extraBits = msb - 1;
            int secondBit = (value >>> extraBits) & 1;
            int slot = ((msb - 1) << 1) + 2 + secondBit; // 4..61 (fits in 6-bit tree!)
            writeFixedTree(enc, tree, slot, 6);
            int rem = value & ((1 << extraBits) - 1);
            int[] lowRow = lowBitProbs[lowCtx];
            for (int b = extraBits - 1; b >= 0; b--) {
                int bit = (rem >>> b) & 1;
                if (b < 4 && slot < 32) {
                    enc.encodeBit(lowRow, (slot << 2) | b, bit);
                } else {
                    enc.encodeDirectBit(bit);
                }
            }
        }

        int readCount(BitRangeDec dec, int[][] slotTable, int ctx, int lowCtx) {
            int slot = readFixedTree(dec, slotTable[ctx], 6);
            if (slot < 4) {
                return slot;
            }
            int extraBits = ((slot - 2) >>> 1);
            int secondBit = (slot - 2) & 1;
            int rem = 0;
            int[] lowRow = lowBitProbs[lowCtx];
            for (int b = extraBits - 1; b >= 0; b--) {
                int bit;
                if (b < 4 && slot < 32) {
                    bit = dec.decodeBit(lowRow, (slot << 2) | b);
                } else {
                    bit = dec.decodeDirectBit();
                }
                rem = (rem << 1) | bit;
            }
            return ((2 | secondBit) << extraBits) | rem;
        }
    }

    private static final class BitRangeEnc {
        private final ByteOut out;
        private long low = 0L;
        private int range = -1;
        private int cache = 0;
        private long cacheSize = 1L;

        BitRangeEnc(ByteOut out) {
            this.out = out;
        }

        void encodeBit(int[] probs, int idx, int bit) {
            int p = probs[idx];
            int bound = (range >>> PROB_BITS) * p;
            if (bit == 1) {
                range = bound;
                probs[idx] = p + ((PROB_SCALE - p + 20) >>> 5);
            } else {
                low += bound & 0xFFFFFFFFL;
                range -= bound;
                probs[idx] = p - ((p + 20) >>> 5);
            }
            while ((range & 0xFFFFFFFFL) < TOP) {
                range <<= 8;
                shiftLow();
            }
        }

        void encodeDirectBit(int bit) {
            range >>>= 1;
            if (bit == 1) {
                low += range & 0xFFFFFFFFL;
            }
            if ((range & 0xFFFFFFFFL) < TOP) {
                range <<= 8;
                shiftLow();
            }
        }

        private void shiftLow() {
            long lo = low;
            if ((int) (lo >>> 32) != 0 || (lo & 0xFFFFFFFFL) < 0xFF000000L) {
                int carry = (int) (lo >>> 32);
                int temp = cache;
                do {
                    out.writeByte(temp + carry);
                    temp = 0xFF;
                } while (--cacheSize != 0L);
                cache = (int) (lo >>> 24) & 0xFF;
            }
            cacheSize++;
            low = ((lo & 0xFFFFFFFFL) << 8) & 0xFFFFFFFFL;
        }

        void finish() {
            for (int i = 0; i < 5; i++) {
                shiftLow();
            }
        }
    }

    private static final class BitRangeDec {
        private final ByteIn in;
        private int range = -1;
        private int code = 0;

        BitRangeDec(ByteIn in) {
            this.in = in;
            for (int i = 0; i < 5; i++) {
                code = (code << 8) | nextByte();
            }
        }

        private int nextByte() {
            return in.hasRemaining() ? (in.readByte() & 0xFF) : 0;
        }

        int decodeBit(int[] probs, int idx) {
            int p = probs[idx];
            int bound = (range >>> PROB_BITS) * p;
            int bit;
            if ((code & 0xFFFFFFFFL) < (bound & 0xFFFFFFFFL)) {
                range = bound;
                probs[idx] = p + ((PROB_SCALE - p + 20) >>> 5);
                bit = 1;
            } else {
                code -= bound;
                range -= bound;
                probs[idx] = p - ((p + 20) >>> 5);
                bit = 0;
            }
            while ((range & 0xFFFFFFFFL) < TOP) {
                range <<= 8;
                code = (code << 8) | nextByte();
            }
            return bit;
        }

        int decodeDirectBit() {
            range >>>= 1;
            int bit;
            if ((code & 0xFFFFFFFFL) >= (range & 0xFFFFFFFFL)) {
                code -= range;
                bit = 1;
            } else {
                bit = 0;
            }
            if ((range & 0xFFFFFFFFL) < TOP) {
                range <<= 8;
                code = (code << 8) | nextByte();
            }
            return bit;
        }
    }
}
