package hypergrass;

import hypergrass.EimdlMosaicCodec.EncodedMosaic;
import hypergrass.EimdlMosaicCodec.MosaicConfig;
import hypergrass.EimdlMosaicCodec.MosaicStats;
import hypergrass.FastaStreamFactorizer.FactorizedChunk;
import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import java.util.zip.CRC32;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * Unified HyperGraSS + RC-MDL Production Engine ({@code HGS2} container format) implementing all
 * three compression regimes under the Evolutionary-Information MDL (EIMDL) &amp; Amdahl-Memory
 * Laws:
 *
 * <ul>
 *   <li><b>Mode S ({@code mode=0}):</b> Self-contained single genome / sequence (existing or
 *       synthetic/non-existing), combining bidirectional ({@code +/-}) approximate self-mosaic,
 *       tandem periodicity, and warm-primed parallel {@link FastNoveltyCoder}.</li>
 *   <li><b>Mode R ({@code mode=1}):</b> Reference-conditioned delta archive with 32-byte
 *       cryptographic SHA-256 Merkle manifest pinning over external reference panels.</li>
 *   <li><b>Mode C ({@code mode=2}):</b> Self-contained pangenome collection archive using
 *       syncmer-Jaccard phylogenetic Prim ordering, lock-free parallel multi-core rank mosaic
 *       encoding, bounded-memory clade-hub indexing, and a unified pan-genome novelty superstream.
 *       </li>
 * </ul>
 */
public final class HyperGraSSFusionEngine {

    public static final byte MAGIC_0 = 'H';
    public static final byte MAGIC_1 = 'G';
    public static final byte MAGIC_2 = 'S';
    public static final byte MAGIC_3 = '2';

    public static final int MODE_SINGLE = 0;
    public static final int MODE_REF = 1;
    public static final int MODE_COLLECTION = 2;

    private static final int FLAG_RAW_FLAT_ACGT = 1 << 2;
    private static final int FLAG_RAW_RNA_U = 1 << 3;

    /** Default target DNA size (4 MiB bases = 4 aligned 256 KiB pages) for headerless streams. */
    private static final int TARGET_COLLECTION_UNIT_BASES = Packed2BitSequence.PAGE_BASES * 4;

    private HyperGraSSFusionEngine() {}

    /** Configuration options for {@code HGS2} compression. */
    public static final class FusionOptions {
        public int mode = -1; // -1 = auto (MODE_REF if panel != null, MODE_SINGLE if singleFlag, else auto)
        public boolean singleGenome = false;
        public boolean boundedMemory = false;
        public int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
        public Path referencePanel = null; // directory of genomes or single reference file
        public Set<String> excludePanelNames = new HashSet<>();
        public List<String> explicitPanelNames = null;
        public int panelLimit = 0;
        public MosaicConfig mosaicConfig = MosaicConfig.DEFAULT;
    }

    /** End-to-end compression telemetry for CLI &amp; JSON benchmark harnesses. */
    public record CompressionReport(
            int mode,
            long rawBytes,
            long targetDnaBases,
            long archiveBytes,
            int panelGenomes,
            long panelBases,
            long panelMatchBases,
            long selfMatchBases,
            long tandemMatchBases,
            long noveltyBases,
            long panelTokens,
            long selfTokens,
            long tandemTokens,
            long literalRuns,
            long indelSplices,
            long substitutions,
            double compressSeconds) {}

    /** Byte-accounting summary of an {@code HGS2} archive for {@link ArchiveByteAccounting}. */
    public record Hgs2Accounting(
            int mode,
            long containerHeaderBytes,
            long panelPinningBytes,
            long metadataBytes,
            long mosaicBytes,
            long noveltyBytes,
            long totalBytes) {}

    public static boolean isHgs2Magic(byte[] prefix) {
        return prefix != null
                && prefix.length >= 4
                && prefix[0] == MAGIC_0
                && prefix[1] == MAGIC_1
                && prefix[2] == MAGIC_2
                && prefix[3] == MAGIC_3;
    }

    public static boolean isHgs2File(Path path) throws IOException {
        try (var in = Files.newInputStream(path)) {
            byte[] hdr = in.readNBytes(4);
            return isHgs2Magic(hdr);
        }
    }

    // =========================================================================
    // Public Compression Entrypoints
    // =========================================================================

    private record ParsedRawInput(
            long rawLength,
            int crcVal,
            int mode,
            boolean rawFlatAcgt,
            boolean rawRnaU,
            int totalRecords,
            List<InputUnit> units) {}

    public static CompressionReport compressFile(Path input, Path output, FusionOptions opts)
            throws IOException {
        long t0 = System.nanoTime();
        ParsedRawInput parsed = parseSource(input, null, opts);
        int initOutCap =
                (opts.boundedMemory || Runtime.getRuntime().maxMemory() <= 320L * 1024 * 1024)
                        ? 65536
                        : (int) Math.min(16_777_216L, Math.max(256L, parsed.rawLength / 16L));
        ByteOut archiveOut = new ByteOut(initOutCap);
        CompressionReport rep = compressParsed(parsed, archiveOut, opts, t0);
        try (OutputStream fos = Files.newOutputStream(output)) {
            fos.write(archiveOut.internalBuffer(), 0, archiveOut.size());
        }
        return new CompressionReport(
                rep.mode,
                rep.rawBytes,
                rep.targetDnaBases,
                archiveOut.size(),
                rep.panelGenomes,
                rep.panelBases,
                rep.panelMatchBases,
                rep.selfMatchBases,
                rep.tandemMatchBases,
                rep.noveltyBases,
                rep.panelTokens,
                rep.selfTokens,
                rep.tandemTokens,
                rep.literalRuns,
                rep.indelSplices,
                rep.substitutions,
                (System.nanoTime() - t0) / 1e9);
    }

    public static byte[] compress(byte[] raw, FusionOptions opts) throws IOException {
        long t0 = System.nanoTime();
        ParsedRawInput parsed = parseSource(null, raw, opts);
        ByteOut out = new ByteOut(Math.max(128, raw.length / 16));
        compressParsed(parsed, out, opts, t0);
        return out.toByteArray();
    }

    private static ParsedRawInput parseSource(Path inputPath, byte[] inMemRaw, FusionOptions opts)
            throws IOException {
        if (inMemRaw == null && inputPath != null) {
            ParsedRawInput streamed = tryStreamFlatFile(inputPath, opts);
            if (streamed != null) {
                return streamed;
            }
        }
        byte[] raw = (inMemRaw != null) ? inMemRaw : Files.readAllBytes(inputPath);
        CRC32 crc = new CRC32();
        crc.update(raw, 0, raw.length);
        int crcVal = (int) crc.getValue();

        int mode;
        if (opts.referencePanel != null) {
            mode = MODE_REF;
        } else if (opts.mode >= 0) {
            mode = opts.mode;
        } else if (opts.singleGenome || raw.length < 6_500_000) {
            mode = MODE_SINGLE;
        } else {
            mode = MODE_COLLECTION;
        }

        Packed2BitSequence flatSeq = Packed2BitSequence.tryFromFlatUppercaseAcgt(raw, raw.length);
        boolean rawFlatAcgt = (flatSeq != null && raw.length > 0);
        boolean rawRnaU = false;
        if (!rawFlatAcgt && raw.length > 0 && raw[0] != '>') {
            Packed2BitSequence uSeq = tryFromFlatUppercaseAcguOnly(raw, raw.length);
            if (uSeq != null) {
                flatSeq = uSeq;
                rawFlatAcgt = true;
                rawRnaU = true;
            }
        }

        List<InputUnit> units;
        int totalRecords;
        if (raw.length == 0) {
            units = List.of();
            totalRecords = 0;
            rawFlatAcgt = false;
        } else if (rawFlatAcgt) {
            totalRecords = 1;
            if (mode == MODE_COLLECTION && flatSeq.length > TARGET_COLLECTION_UNIT_BASES * 3 / 2) {
                units = sliceFlatSequenceIntoUnits(flatSeq, TARGET_COLLECTION_UNIT_BASES);
            } else {
                units = List.of(new InputUnit(flatSeq, null, 0, 1));
            }
        } else {
            if (mode == MODE_COLLECTION) {
                boolean boundedParse =
                        opts.boundedMemory
                                || Runtime.getRuntime().maxMemory() <= 320L * 1024 * 1024;
                units =
                        partitionFactorizedCollectionUnits(
                                raw,
                                TARGET_COLLECTION_UNIT_BASES,
                                Math.max(1, opts.threads),
                                boundedParse);
                int sumRec = 0;
                for (InputUnit u : units) {
                    sumRec += u.recordCount;
                }
                totalRecords = sumRec;
                if (units.size() <= 1) {
                    mode = MODE_SINGLE;
                }
            } else {
                FactorizedChunk chunk = FastaStreamFactorizer.factorize(raw, 0, raw.length);
                Packed2BitSequence seq =
                        Packed2BitSequence.fromTwoBitBytes(chunk.dnaBases, chunk.dnaLen);
                // Free 1-byte-per-base dnaBases immediately once 2-bit Packed2BitSequence is built!
                chunk.dnaBases = null;
                int recCnt = chunk.getRecordSpans().length;
                totalRecords = recCnt;
                units = List.of(new InputUnit(seq, chunk, 0, recCnt));
            }
        }
        return new ParsedRawInput(
                raw.length, crcVal, mode, rawFlatAcgt, rawRnaU, totalRecords, units);
    }

    private static CompressionReport compressParsed(
            ParsedRawInput parsed, ByteOut out, FusionOptions opts, long startNano)
            throws IOException {
        int threads = Math.max(1, opts.threads);
        boolean autoLowMem = Runtime.getRuntime().maxMemory() <= 320L * 1024 * 1024;
        boolean bounded = opts.boundedMemory || autoLowMem;

        final long rawLength = parsed.rawLength;
        final int crcVal = parsed.crcVal;
        int mode = parsed.mode;
        boolean rawFlatAcgt = parsed.rawFlatAcgt;
        boolean rawRnaU = parsed.rawRnaU;
        int totalRecords = parsed.totalRecords;
        List<InputUnit> units = parsed.units;
        parsed = null; // allow GC of ParsedRawInput wrapper once fields are extracted

        int headerFlags = (mode & 3);
        if (rawFlatAcgt) headerFlags |= FLAG_RAW_FLAT_ACGT;
        if (rawRnaU) headerFlags |= FLAG_RAW_RNA_U;

        // Write HGS2 header
        out.writeByte(MAGIC_0);
        out.writeByte(MAGIC_1);
        out.writeByte(MAGIC_2);
        out.writeByte(MAGIC_3);
        out.writeByte(headerFlags);
        out.writeLongLE(rawLength);
        out.writeIntLE(crcVal);
        out.writeVarInt(totalRecords);
        out.writeVarInt(units.size());

        // Mode R panel loading & cryptographic manifest pinning
        LoadedPanel loadedPanel = null;
        if (mode == MODE_REF) {
            loadedPanel =
                    loadReferencePanelForCompress(
                            opts.referencePanel,
                            opts.excludePanelNames,
                            opts.explicitPanelNames,
                            opts.panelLimit);
            out.writeVarInt(loadedPanel.names.length);
            byte[] compNames = compressPanelNames(loadedPanel.names);
            out.writeVarInt(compNames.length);
            out.writeBytes(compNames);
            out.writeBytes(loadedPanel.manifestSha256);
        }

        if (units.isEmpty()) {
            return new CompressionReport(
                    mode,
                    0L,
                    0L,
                    out.size(),
                    loadedPanel != null ? loadedPanel.names.length : 0,
                    loadedPanel != null ? loadedPanel.totalBases : 0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    (System.nanoTime() - startNano) / 1e9);
        }

        final int numUnits = units.size();

        // Determine rank ordering for Mode C (Prim's MST on canonical open-syncmer bottom sketches)
        int[] rankToUnit = new int[numUnits];
        int[] mstParentRank = new int[numUnits];
        Arrays.fill(mstParentRank, -1);
        int[][] rankDonorSimilarityOrder = new int[numUnits][];

        if (mode == MODE_COLLECTION && numUnits > 1) {
            buildPhylogeneticPrimOrder(
                    units, threads, rankToUnit, mstParentRank, rankDonorSimilarityOrder);
            for (int r = 0; r < numUnits; r++) {
                out.writeVarInt(rankToUnit[r]);
            }
            for (int r = 0; r < numUnits; r++) {
                out.writeVarInt(r); // maxAncestorRank = r (strictly causal < r)
            }
        } else {
            for (int u = 0; u < numUnits; u++) {
                rankToUnit[u] = u;
            }
        }

        // Write per-unit framing / 6-channel factorized metadata
        if (!rawFlatAcgt) {
            for (int u = 0; u < numUnits; u++) {
                InputUnit iu = units.get(u);
                out.writeVarInt(iu.firstRecordIdx);
                out.writeVarInt(iu.recordCount);
                ByteOut metaRaw = new ByteOut(256);
                FastaStreamFactorizer.writeAuxMetadata(metaRaw, iu.chunk);
                out.writeVarInt(metaRaw.size());
                out.writeBytes(metaRaw.internalBuffer(), 0, metaRaw.size());
                // Free bulky chunk arrays now that metadata is serialized & DNA is 2-bit packed!
                iu.chunk = null;
            }
        } else {
            for (int u = 0; u < numUnits; u++) {
                out.writeVarInt(units.get(u).seq.length);
            }
        }

        // Encode mosaic streams across units
        EncodedMosaic[] rankMosaics = new EncodedMosaic[numUnits];
        long panelBasesTotal = 0L;
        int panelGenomeCount = 0;
        long maxMemAtStart = Runtime.getRuntime().maxMemory();
        int safeStage1Threads =
                (maxMemAtStart <= 145L * 1024 * 1024)
                        ? 1
                        : (maxMemAtStart <= 205L * 1024 * 1024
                                ? Math.min(threads, 2)
                                : (bounded ? Math.min(threads, 4) : threads));

        if (mode == MODE_SINGLE) {
            InputUnit u0 = units.get(0);
            rankMosaics[0] =
                    EimdlMosaicCodec.encodeMosaic(
                            u0.seq, null, 0, -1, null, opts.mosaicConfig, safeStage1Threads);
        } else if (mode == MODE_REF) {
            panelGenomeCount = loadedPanel.seqs.length;
            panelBasesTotal = loadedPanel.totalBases;
            EimdlSyncmerIndex panelIdx =
                    loadedPanel.seqs.length > 0
                            ? EimdlSyncmerIndex.build(
                                    loadedPanel.seqs, null, !bounded, threads)
                            : null;

            // Compute fast similarity ranking of panel donors for each target unit
            final Packed2BitSequence[] panelSeqs = loadedPanel.seqs;
            long[][] panelSketches = new long[panelSeqs.length][];
            IntStream.range(0, panelSeqs.length)
                    .parallel()
                    .forEach(
                            d ->
                                    panelSketches[d] =
                                            panelSeqs[d].computeBottomSketch(
                                                    0, panelSeqs[d].length, 2048));

            for (int u = 0; u < numUnits; u++) {
                Packed2BitSequence tgtSeq = units.get(u).seq;
                int bestTwin = -1;
                int[] simOrder = null;
                if (panelIdx != null && loadedPanel.seqs.length > 0) {
                    long[] tSketch = tgtSeq.computeBottomSketch(0, tgtSeq.length, 2048);
                    double[] sims = new double[loadedPanel.seqs.length];
                    Integer[] orderBox = new Integer[loadedPanel.seqs.length];
                    for (int d = 0; d < loadedPanel.seqs.length; d++) {
                        sims[d] = Packed2BitSequence.sketchJaccard(tSketch, panelSketches[d]);
                        orderBox[d] = d;
                    }
                    Arrays.sort(
                            orderBox,
                            (a, b) -> {
                                int c = Double.compare(sims[b], sims[a]);
                                return c != 0 ? c : Integer.compare(a, b);
                            });
                    simOrder = new int[orderBox.length];
                    for (int d = 0; d < orderBox.length; d++) {
                        simOrder[d] = orderBox[d];
                    }
                    if (sims[simOrder[0]] > 0.05) {
                        bestTwin = simOrder[0];
                    }
                }
                rankMosaics[u] =
                        EimdlMosaicCodec.encodeMosaic(
                                tgtSeq,
                                panelIdx,
                                panelIdx != null ? panelIdx.numDonors : 0,
                                bestTwin,
                                simOrder,
                                opts.mosaicConfig,
                                safeStage1Threads);
                units.get(u).seq = null;
            }
            panelIdx = null;
            loadedPanel = null;
            units = null;
        } else {
            // MODE_COLLECTION: Rank-ordered parallel DAG mosaic encoding!
            Packed2BitSequence[] orderedDonors = new Packed2BitSequence[numUnits];
            for (int r = 0; r < numUnits; r++) {
                InputUnit iu = units.get(rankToUnit[r]);
                orderedDonors[r] = iu.seq;
                iu.seq = null;
            }
            units.clear();
            units = null;

            // Compute subsampling dynamically from (totalBankBases, maxMemory) so CSR + bank + workers
            // fit strictly inside tight heaps (-Xmx128m .. -Xmx1g) even at 300 Mbp .. 1 Gbp+!
            long maxMemBytes = Runtime.getRuntime().maxMemory();
            long totalBankBases = 0L;
            for (Packed2BitSequence seq : orderedDonors) {
                if (seq != null) {
                    totalBankBases += seq.length;
                }
            }
            long bankBytes = (totalBankBases + 3L) >>> 2;
            boolean ultraTight = maxMemBytes <= 205L * 1024 * 1024;
            boolean extremeTight = maxMemBytes <= 145L * 1024 * 1024;
            boolean bankHeavy = bankBytes * 5L > maxMemBytes * 2L;
            boolean memPressured = bounded || (bankBytes * 5L > maxMemBytes);

            int maxMosaicWorkers =
                    (extremeTight || bankHeavy)
                            ? 1
                            : (ultraTight
                                    ? Math.min(threads, 2)
                                    : (bounded ? Math.min(threads, 4) : threads));
            int workerThreads = Math.min(maxMosaicWorkers, numUnits);
            long workerReserve = (long) workerThreads * 6L * 1024 * 1024 + (4L << 20);
            long csrBudget =
                    Math.max(6L << 20, (long) (maxMemBytes * 0.72) - bankBytes - workerReserve);

            int subsampleShift = bounded ? (ultraTight ? 2 : 1) : 0;
            if (memPressured) {
                long estRawEntries =
                        Math.max(1024L, totalBankBases / Packed2BitSequence.W_PRIMARY);
                while (subsampleShift < 6) {
                    long ent = Math.max(1024L, estRawEntries >>> subsampleShift);
                    long estBucketBytes = Math.min(8L << 20, Math.max(64L << 10, ent * 2L));
                    if (ent * 6L + estBucketBytes <= csrBudget) {
                        break;
                    }
                    subsampleShift++;
                }
            }
            int buildThreads =
                    (extremeTight || bankHeavy)
                            ? 1
                            : (bounded ? Math.min(threads, 4) : threads);
            EimdlSyncmerIndex sharedDagIdx =
                    EimdlSyncmerIndex.buildWithSubsampling(
                            orderedDonors,
                            subsampleShift,
                            true,
                            !bounded && !memPressured && numUnits <= 16,
                            buildThreads);
            orderedDonors = null;
            boolean useIntraUnitChunking =
                    maxMosaicWorkers > 1
                            && (workerThreads <= 1
                                    || (numUnits <= 4
                                            && totalBankBases >= 16_000_000L
                                            && !memPressured
                                            && !bounded));
            if (workerThreads <= 1 || useIntraUnitChunking) {
                int chunkThreads = useIntraUnitChunking ? maxMosaicWorkers : 1;
                for (int r = 0; r < numUnits; r++) {
                    int dStart = sharedDagIdx.donorOffsets[r];
                    int dLen = sharedDagIdx.donorOffsets[r + 1] - dStart;
                    Packed2BitSequence targetSeq = sharedDagIdx.bank.sliceSeq(dStart, dLen);
                    if (r == 0) {
                        rankMosaics[0] =
                                EimdlMosaicCodec.encodeMosaic(
                                        targetSeq,
                                        null,
                                        0,
                                        -1,
                                        null,
                                        opts.mosaicConfig,
                                        chunkThreads);
                    } else {
                        rankMosaics[r] =
                                EimdlMosaicCodec.encodeMosaic(
                                        targetSeq,
                                        sharedDagIdx,
                                        r,
                                        mstParentRank[r],
                                        rankDonorSimilarityOrder[r],
                                        opts.mosaicConfig,
                                        chunkThreads);
                    }
                }
            } else {
                final EimdlSyncmerIndex dagIdx = sharedDagIdx;
                final int[] parentRanks = mstParentRank;
                final int[][] simOrders = rankDonorSimilarityOrder;
                final EncodedMosaic[] outMosaics = rankMosaics;
                ForkJoinPool pool = new ForkJoinPool(workerThreads);
                try {
                    pool.submit(
                                    () ->
                                            IntStream.range(0, numUnits)
                                                    .parallel()
                                                    .forEach(
                                                            r -> {
                                                                int dStart =
                                                                        dagIdx.donorOffsets[r];
                                                                int dLen =
                                                                        dagIdx.donorOffsets[r + 1]
                                                                                - dStart;
                                                                Packed2BitSequence targetSeq =
                                                                        dagIdx.bank.sliceSeq(
                                                                                dStart, dLen);
                                                                if (r == 0) {
                                                                    outMosaics[0] =
                                                                            EimdlMosaicCodec
                                                                                    .encodeMosaic(
                                                                                            targetSeq,
                                                                                            null,
                                                                                            0,
                                                                                            -1,
                                                                                            null,
                                                                                            opts.mosaicConfig);
                                                                } else {
                                                                    outMosaics[r] =
                                                                            EimdlMosaicCodec
                                                                                    .encodeMosaic(
                                                                                            targetSeq,
                                                                                            dagIdx,
                                                                                            r,
                                                                                            parentRanks[
                                                                                                    r],
                                                                                            simOrders[
                                                                                                    r],
                                                                                            opts.mosaicConfig);
                                                                }
                                                            }))
                            .join();
                } finally {
                    pool.shutdown();
                }
            }
            sharedDagIdx = null; // Free CSR index + 2-bit bank before Novelty Superstream coding!
            mstParentRank = null;
            rankDonorSimilarityOrder = null;
        }

        // Write mosaic streams & assemble Unified Pan-Genome Novelty Superstream
        long totalNoveltyLong = 0L;
        long totalTargetBases = 0L;
        long sumPanelMatch = 0L;
        long sumSelfMatch = 0L;
        long sumTandemMatch = 0L;
        long sumPanelTok = 0L;
        long sumSelfTok = 0L;
        long sumTandemTok = 0L;
        long sumLitRuns = 0L;
        long sumSplices = 0L;
        long sumSnps = 0L;

        byte[][] unitNovSlices = new byte[numUnits][];
        int[] unitNovLens = new int[numUnits];
        for (int r = 0; r < numUnits; r++) {
            EncodedMosaic em = rankMosaics[r];
            rankMosaics[r] = null; // release EncodedMosaic & mosaicBytes immediately after writing!
            out.writeVarInt(em.mosaicBytes().length);
            out.writeBytes(em.mosaicBytes());
            totalNoveltyLong += em.noveltyLen();
            unitNovSlices[r] = em.noveltyBases();
            unitNovLens[r] = em.noveltyLen();

            MosaicStats st = em.stats();
            totalTargetBases += st.targetBases();
            sumPanelMatch += st.panelMatchBases();
            sumSelfMatch += st.selfMatchBases();
            sumTandemMatch += st.tandemMatchBases();
            sumPanelTok += st.panelTokens();
            sumSelfTok += st.selfTokens();
            sumTandemTok += st.tandemTokens();
            sumLitRuns += st.literalRuns();
            sumSplices += st.indelSplices();
            sumSnps += st.substitutions();
        }
        rankMosaics = null;

        if (totalNoveltyLong > Integer.MAX_VALUE - 64) {
            throw new IllegalStateException("Novelty superstream exceeds 2.14 Gbp");
        }
        int totalNovelty = (int) totalNoveltyLong;
        byte[] noveltySuper = new byte[totalNovelty];
        int novDst = 0;
        for (int r = 0; r < numUnits; r++) {
            int nLen = unitNovLens[r];
            if (nLen > 0) {
                System.arraycopy(unitNovSlices[r], 0, noveltySuper, novDst, nLen);
                novDst += nLen;
            }
            unitNovSlices[r] = null; // release per-unit novelty buffer incrementally!
        }
        unitNovSlices = null;
        unitNovLens = null;

        byte[] noveltyComp = FastNoveltyCoder.encode(noveltySuper, totalNovelty, threads, bounded);
        noveltySuper = null;

        out.writeVarInt(totalNovelty);
        out.writeVarInt(noveltyComp.length);
        out.writeBytes(noveltyComp);

        double secs = (System.nanoTime() - startNano) / 1e9;
        return new CompressionReport(
                mode,
                rawLength,
                totalTargetBases,
                out.size(),
                panelGenomeCount,
                panelBasesTotal,
                sumPanelMatch,
                sumSelfMatch,
                sumTandemMatch,
                totalNovelty,
                sumPanelTok,
                sumSelfTok,
                sumTandemTok,
                sumLitRuns,
                sumSplices,
                sumSnps,
                secs);
    }

    // =========================================================================
    // Public Decompression, Streaming, Inspect, Query, & Byte Accounting
    // =========================================================================

    public static byte[] decompress(byte[] archive, Path referencePanel) throws IOException {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        decompressToStream(archive, referencePanel, baos);
        return baos.toByteArray();
    }

    public static void decompressFile(Path archivePath, Path outputPath, Path referencePanel)
            throws IOException {
        byte[] archive = Files.readAllBytes(archivePath);
        try (OutputStream out = new BufferedOutputStream(Files.newOutputStream(outputPath), 65536)) {
            decompressToStream(archive, referencePanel, out);
        }
    }

    public static void decompressToStream(
            byte[] archive, Path referencePanel, OutputStream outStream) throws IOException {
        ParsedHgs2 hdr = parseHeaderAndDirectory(archive);
        if (hdr.rawByteLen == 0L || hdr.numUnits == 0) {
            if (hdr.mode == MODE_REF && hdr.panelCount > 0) {
                verifyAndLoadPanelForDecode(referencePanel, hdr.panelNames, hdr.panelManifestSha256);
            }
            if (hdr.rawCrc32 != 0) {
                throw new IllegalStateException("CRC-32 mismatch on empty archive");
            }
            return;
        }

        int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
        boolean bounded = Runtime.getRuntime().maxMemory() <= 320L * 1024 * 1024;

        Packed2BitSequence refBank = null;
        int[] refOffsets = null;
        int refDonorCount = 0;
        if (hdr.mode == MODE_REF && hdr.panelCount > 0) {
            Packed2BitSequence[] panelSeqs =
                    verifyAndLoadPanelForDecode(
                            referencePanel, hdr.panelNames, hdr.panelManifestSha256);
            refDonorCount = panelSeqs.length;
            refOffsets = new int[refDonorCount + 1];
            refBank = Packed2BitSequence.concatenate(panelSeqs, refOffsets);
        }

        // Decode Unified Novelty DNA Superstream in parallel
        ByteIn novIn = new ByteIn(archive, hdr.noveltyOffset, hdr.noveltyPayloadLen);
        byte[] noveltySuper =
                FastNoveltyCoder.decode(novIn, hdr.totalNoveltyBases, threads, bounded);

        // Slice novelty & target lengths per rank
        int[] rankTargetLen = new int[hdr.numUnits];
        int[] rankNoveltyStart = new int[hdr.numUnits];
        int[] rankNoveltyLen = new int[hdr.numUnits];
        int nCur = 0;
        long totalTargetBasesLong = 0L;
        for (int r = 0; r < hdr.numUnits; r++) {
            byte[] mBytes =
                    Arrays.copyOfRange(
                            archive, hdr.mosaicOffsets[r], hdr.mosaicOffsets[r] + hdr.mosaicLens[r]);
            int[] tlNl = EimdlMosaicCodec.peekLengths(mBytes);
            rankTargetLen[r] = tlNl[0];
            rankNoveltyStart[r] = nCur;
            rankNoveltyLen[r] = tlNl[1];
            nCur += tlNl[1];
            totalTargetBasesLong += tlNl[0];
        }
        if (nCur != hdr.totalNoveltyBases) {
            throw new IllegalStateException("Novelty base sum mismatch in mosaic headers");
        }
        if (totalTargetBasesLong > Integer.MAX_VALUE - 64) {
            throw new IllegalStateException("Total DNA bases exceed 2.14 Gbp");
        }

        // Single 2-bit bank (0.25 B/base = 30 MB for 120 Mbp!) holds all decoded ranks with zero re-concatenation!
        int[] rankDnaOffsets = new int[hdr.numUnits + 1];
        for (int r = 0; r < hdr.numUnits; r++) {
            rankDnaOffsets[r + 1] = rankDnaOffsets[r] + rankTargetLen[r];
        }
        int totalDnaBases = rankDnaOffsets[hdr.numUnits];
        Packed2BitSequence decodedBank = Packed2BitSequence.allocatePaged(totalDnaBases);
        int[] unitToRank = new int[hdr.numUnits];

        if (hdr.mode == MODE_COLLECTION && hdr.numUnits > 1) {
            for (int r = 0; r < hdr.numUnits; r++) {
                unitToRank[hdr.rankToUnit[r]] = r;
                byte[] mBytes =
                        Arrays.copyOfRange(
                                archive,
                                hdr.mosaicOffsets[r],
                                hdr.mosaicOffsets[r] + hdr.mosaicLens[r]);
                byte[] novSlice =
                        Arrays.copyOfRange(
                                noveltySuper,
                                rankNoveltyStart[r],
                                rankNoveltyStart[r] + rankNoveltyLen[r]);
                byte[] twoBit =
                        (r == 0)
                                ? EimdlMosaicCodec.decodeMosaicToTwoBit(
                                        mBytes, novSlice, null, null, 0)
                                : EimdlMosaicCodec.decodeMosaicToTwoBit(
                                        mBytes,
                                        novSlice,
                                        decodedBank,
                                        rankDnaOffsets,
                                        hdr.maxAncestorRank[r]);
                decodedBank.writeTwoBitSlice(rankDnaOffsets[r], twoBit, twoBit.length);
            }
        } else {
            for (int u = 0; u < hdr.numUnits; u++) {
                unitToRank[u] = u;
                byte[] mBytes =
                        Arrays.copyOfRange(
                                archive,
                                hdr.mosaicOffsets[u],
                                hdr.mosaicOffsets[u] + hdr.mosaicLens[u]);
                byte[] novSlice =
                        Arrays.copyOfRange(
                                noveltySuper,
                                rankNoveltyStart[u],
                                rankNoveltyStart[u] + rankNoveltyLen[u]);
                byte[] twoBit =
                        EimdlMosaicCodec.decodeMosaicToTwoBit(
                                mBytes, novSlice, refBank, refOffsets, refDonorCount);
                decodedBank.writeTwoBitSlice(rankDnaOffsets[u], twoBit, twoBit.length);
            }
        }

        noveltySuper = null;
        refBank = null;

        // Reconstruct exact byte stream one unit at a time (~4.5 MB live buffer) & verify CRC-32
        CRC32 crc = new CRC32();
        long totalWritten = 0L;
        if (hdr.rawFlatAcgt) {
            byte[] map =
                    hdr.rawRnaU
                            ? new byte[] {'A', 'C', 'G', 'U'}
                            : new byte[] {'A', 'C', 'G', 'T'};
            byte[] ioBuf = new byte[65536];
            for (int u = 0; u < hdr.numUnits; u++) {
                int r = unitToRank[u];
                int start = rankDnaOffsets[r];
                int len = rankTargetLen[r];
                int idx = 0;
                while (idx < len) {
                    int take = Math.min(ioBuf.length, len - idx);
                    int j = 0;
                    while (j + 32 <= take) {
                        long w = decodedBank.word32Forward(start + idx + j);
                        for (int b = 0; b < 32; b++) {
                            ioBuf[j + b] = map[(int) (w & 3L)];
                            w >>>= 2;
                        }
                        j += 32;
                    }
                    while (j < take) {
                        ioBuf[j] = map[decodedBank.baseAt(start + idx + j)];
                        j++;
                    }
                    outStream.write(ioBuf, 0, take);
                    crc.update(ioBuf, 0, take);
                    totalWritten += take;
                    idx += take;
                }
            }
        } else {
            for (int u = 0; u < hdr.numUnits; u++) {
                int r = unitToRank[u];
                byte[] unitTwoBit =
                        decodedBank.extractTwoBitSlice(rankDnaOffsets[r], rankTargetLen[r]);
                ByteIn metaIn =
                        new ByteIn(archive, hdr.unitMetaOffsets[u], hdr.unitMetaLens[u]);
                FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(metaIn);
                totalWritten +=
                        FastaStreamFactorizer.reconstructTo(
                                chunk, unitTwoBit, outStream, crc);
            }
        }

        if (totalWritten != hdr.rawByteLen) {
            throw new IllegalStateException(
                    "Decompressed byte length mismatch: " + totalWritten + " vs " + hdr.rawByteLen);
        }
        if ((int) crc.getValue() != hdr.rawCrc32) {
            throw new IllegalStateException(
                    String.format(
                            "CRC-32 integrity failure: got 0x%08X, expected 0x%08X",
                            (int) crc.getValue(), hdr.rawCrc32));
        }
    }

    /** Random-access record / sub-interval query on an {@code HGS2} archive. */
    public static byte[] query(
            Path archivePath, int recordIdx, int qStart, int qEnd, Path referencePanel)
            throws IOException {
        byte[] archive = Files.readAllBytes(archivePath);
        ParsedHgs2 hdr = parseHeaderAndDirectory(archive);
        if (recordIdx < 0 || recordIdx >= hdr.totalRecords) {
            throw new IndexOutOfBoundsException(
                    "Record index " + recordIdx + " out of [0, " + hdr.totalRecords + ")");
        }

        int targetUnit = 0;
        int localRecIdx = recordIdx;
        if (!hdr.rawFlatAcgt) {
            boolean found = false;
            for (int u = 0; u < hdr.numUnits; u++) {
                if (recordIdx >= hdr.firstRecordIdx[u]
                        && recordIdx < hdr.firstRecordIdx[u] + hdr.recordCount[u]) {
                    targetUnit = u;
                    localRecIdx = recordIdx - hdr.firstRecordIdx[u];
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw new IllegalStateException("Record " + recordIdx + " not mapped to any unit");
            }
        }

        // If rawFlatAcgt and single record across units, reconstruct or slice exact overlapping units
        if (hdr.rawFlatAcgt) {
            byte[] full = decompress(archive, referencePanel);
            int s = (qStart >= 0) ? Math.max(0, Math.min(full.length, qStart)) : 0;
            int e = (qStart >= 0) ? ((qEnd < 0) ? full.length : Math.max(s, Math.min(full.length, qEnd))) : full.length;
            return Arrays.copyOfRange(full, s, e);
        }

        int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
        ByteIn novIn = new ByteIn(archive, hdr.noveltyOffset, hdr.noveltyPayloadLen);
        byte[] noveltySuper =
                FastNoveltyCoder.decode(novIn, hdr.totalNoveltyBases, threads, false);

        int[] rankTargetLen = new int[hdr.numUnits];
        int[] rankNoveltyStart = new int[hdr.numUnits];
        int[] rankNoveltyLen = new int[hdr.numUnits];
        int nCur = 0;
        int targetRank = targetUnit;
        for (int r = 0; r < hdr.numUnits; r++) {
            if (hdr.rankToUnit[r] == targetUnit) {
                targetRank = r;
            }
            byte[] mBytes =
                    Arrays.copyOfRange(
                            archive, hdr.mosaicOffsets[r], hdr.mosaicOffsets[r] + hdr.mosaicLens[r]);
            int[] tlNl = EimdlMosaicCodec.peekLengths(mBytes);
            rankTargetLen[r] = tlNl[0];
            rankNoveltyStart[r] = nCur;
            rankNoveltyLen[r] = tlNl[1];
            nCur += tlNl[1];
        }

        byte[] unitTwoBit;
        if (hdr.mode == MODE_REF && hdr.panelCount > 0) {
            Packed2BitSequence[] panelSeqs =
                    verifyAndLoadPanelForDecode(
                            referencePanel, hdr.panelNames, hdr.panelManifestSha256);
            int[] refOffsets = new int[panelSeqs.length + 1];
            Packed2BitSequence refBank = Packed2BitSequence.concatenate(panelSeqs, refOffsets);
            byte[] mBytes =
                    Arrays.copyOfRange(
                            archive,
                            hdr.mosaicOffsets[targetRank],
                            hdr.mosaicOffsets[targetRank] + hdr.mosaicLens[targetRank]);
            byte[] novSlice =
                    Arrays.copyOfRange(
                            noveltySuper,
                            rankNoveltyStart[targetRank],
                            rankNoveltyStart[targetRank] + rankNoveltyLen[targetRank]);
            unitTwoBit =
                    EimdlMosaicCodec.decodeMosaicToTwoBit(
                            mBytes, novSlice, refBank, refOffsets, panelSeqs.length);
        } else if (hdr.mode == MODE_COLLECTION && hdr.numUnits > 1) {
            int[] dOff = new int[targetRank + 1];
            for (int r = 0; r < targetRank; r++) {
                dOff[r + 1] = dOff[r] + rankTargetLen[r];
            }
            int ancBases = dOff[targetRank];
            Packed2BitSequence bank = Packed2BitSequence.allocatePaged(ancBases);
            byte[] latest = null;
            for (int r = 0; r <= targetRank; r++) {
                byte[] mBytes =
                        Arrays.copyOfRange(
                                archive,
                                hdr.mosaicOffsets[r],
                                hdr.mosaicOffsets[r] + hdr.mosaicLens[r]);
                byte[] novSlice =
                        Arrays.copyOfRange(
                                noveltySuper,
                                rankNoveltyStart[r],
                                rankNoveltyStart[r] + rankNoveltyLen[r]);
                if (r == 0) {
                    latest = EimdlMosaicCodec.decodeMosaicToTwoBit(mBytes, novSlice, null, null, 0);
                } else {
                    int anc = hdr.maxAncestorRank[r];
                    latest =
                            EimdlMosaicCodec.decodeMosaicToTwoBit(
                                    mBytes, novSlice, bank, dOff, anc);
                }
                if (r < targetRank) {
                    bank.writeTwoBitSlice(dOff[r], latest, latest.length);
                }
            }
            unitTwoBit = latest;
        } else {
            byte[] mBytes =
                    Arrays.copyOfRange(
                            archive,
                            hdr.mosaicOffsets[targetRank],
                            hdr.mosaicOffsets[targetRank] + hdr.mosaicLens[targetRank]);
            byte[] novSlice =
                    Arrays.copyOfRange(
                            noveltySuper,
                            rankNoveltyStart[targetRank],
                            rankNoveltyStart[targetRank] + rankNoveltyLen[targetRank]);
            unitTwoBit = EimdlMosaicCodec.decodeMosaicToTwoBit(mBytes, novSlice, null, null, 0);
        }

        ByteIn metaIn =
                new ByteIn(archive, hdr.unitMetaOffsets[targetUnit], hdr.unitMetaLens[targetUnit]);
        FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(metaIn);
        return FastaStreamFactorizer.extractRecordOrSlice(
                chunk, unitTwoBit, localRecIdx, qStart, qEnd);
    }

    public static HyperGraSS.ArchiveMetadata inspect(Path archivePath) throws IOException {
        byte[] archive = Files.readAllBytes(archivePath);
        ParsedHgs2 hdr = parseHeaderAndDirectory(archive);
        return new HyperGraSS.ArchiveMetadata(
                hdr.numUnits <= 1,
                hdr.mode != MODE_SINGLE,
                false,
                hdr.numUnits,
                hdr.rawByteLen,
                hdr.totalRecords,
                hdr.crc32UnsignedInt(),
                hdr.totalNoveltyBases,
                hdr.numUnits);
    }

    public static Hgs2Accounting accountArchive(byte[] archive) {
        ParsedHgs2 hdr = parseHeaderAndDirectory(archive);
        return new Hgs2Accounting(
                hdr.mode,
                hdr.headerBytes,
                hdr.panelPinningBytes,
                hdr.metadataBytes,
                hdr.mosaicSumBytes,
                hdr.noveltyFramingAndPayloadBytes,
                archive.length);
    }

    // =========================================================================
    // Phylogenetic Prim Ordering & Clade-Hub Selection
    // =========================================================================

    private static void buildPhylogeneticPrimOrder(
            List<InputUnit> units,
            int threads,
            int[] outRankToUnit,
            int[] outMstParentRank,
            int[][] outRankDonorOrder) {
        final int n = units.size();
        long[][] sketches = new long[n][];
        Runnable sketchTask =
                () ->
                        IntStream.range(0, n)
                                .parallel()
                                .forEach(
                                        i -> {
                                            Packed2BitSequence s = units.get(i).seq;
                                            sketches[i] = s.computeBottomSketch(0, s.length, 2048);
                                        });
        if (threads > 1) {
            ForkJoinPool pool = new ForkJoinPool(Math.min(threads, n));
            try {
                pool.submit(sketchTask).join();
            } finally {
                pool.shutdown();
            }
        } else {
            sketchTask.run();
        }

        double[][] sim = new double[n][n];
        for (int i = 0; i < n; i++) {
            sim[i][i] = 1.0;
            for (int j = i + 1; j < n; j++) {
                double jac = Packed2BitSequence.sketchJaccard(sketches[i], sketches[j]);
                sim[i][j] = jac;
                sim[j][i] = jac;
            }
        }

        // Choose root unit with maximum total similarity to all other units (medoid hub!)
        int root = 0;
        double bestSum = -1.0;
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) s += sim[i][j];
            }
            if (s > bestSum) {
                bestSum = s;
                root = i;
            }
        }

        boolean[] visited = new boolean[n];
        double[] maxSimToTree = new double[n];
        int[] bestTreeRank = new int[n];
        Arrays.fill(maxSimToTree, -1.0);
        Arrays.fill(bestTreeRank, 0);

        visited[root] = true;
        outRankToUnit[0] = root;
        outMstParentRank[0] = -1;

        for (int v = 0; v < n; v++) {
            if (!visited[v]) {
                maxSimToTree[v] = sim[root][v];
                bestTreeRank[v] = 0;
            }
        }

        for (int r = 1; r < n; r++) {
            int pick = -1;
            double bestSim = -1.0;
            for (int v = 0; v < n; v++) {
                if (!visited[v] && maxSimToTree[v] > bestSim) {
                    bestSim = maxSimToTree[v];
                    pick = v;
                }
            }
            visited[pick] = true;
            outRankToUnit[r] = pick;
            outMstParentRank[r] = bestTreeRank[pick];

            // Build similarity-sorted ancestor rank order [0 .. r-1] for unit pick
            Integer[] ancRanks = new Integer[r];
            for (int a = 0; a < r; a++) {
                ancRanks[a] = a;
            }
            final int uPick = pick;
            Arrays.sort(
                    ancRanks,
                    (a, b) -> {
                        int c =
                                Double.compare(
                                        sim[uPick][outRankToUnit[b]], sim[uPick][outRankToUnit[a]]);
                        return c != 0 ? c : Integer.compare(a, b);
                    });
            int[] ord = new int[r];
            for (int a = 0; a < r; a++) {
                ord[a] = ancRanks[a];
            }
            outRankDonorOrder[r] = ord;

            for (int v = 0; v < n; v++) {
                if (!visited[v] && sim[pick][v] > maxSimToTree[v]) {
                    maxSimToTree[v] = sim[pick][v];
                    bestTreeRank[v] = r;
                }
            }
        }
    }

    // =========================================================================
    // Input Partitioning & Reference Panel Pinning Helpers
    // =========================================================================

    private static final class InputUnit {
        Packed2BitSequence seq;
        FactorizedChunk chunk;
        final int firstRecordIdx;
        final int recordCount;

        InputUnit(Packed2BitSequence seq, FactorizedChunk chunk, int firstRecordIdx, int recordCount) {
            this.seq = seq;
            this.chunk = chunk;
            this.firstRecordIdx = firstRecordIdx;
            this.recordCount = recordCount;
        }
    }

    /**
     * Streams flat uppercase ACGT (or ACGU) `.seq` files directly into 2-bit {@link
     * Packed2BitSequence} slices using a 64 KiB I/O buffer, completely eliminating the `120 MB`
     * raw `byte[]` heap allocation so `-Xmx192m` / `-Xmx256m` bounded collections run with >50 MiB
     * headroom. Returns {@code null} immediately if the file contains FASTA headers, newlines, or
     * IUPAC ambiguity symbols.
     */
    private static ParsedRawInput tryStreamFlatFile(Path inputPath, FusionOptions opts)
            throws IOException {
        long fileSize = Files.size(inputPath);
        if (fileSize <= 0L || fileSize > Integer.MAX_VALUE - 64L) {
            return null;
        }
        int n = (int) fileSize;

        int mode;
        if (opts.referencePanel != null) {
            mode = MODE_REF;
        } else if (opts.mode >= 0) {
            mode = opts.mode;
        } else if (opts.singleGenome || n < 6_500_000) {
            mode = MODE_SINGLE;
        } else {
            mode = MODE_COLLECTION;
        }

        int numSlices =
                (mode == MODE_COLLECTION && n > TARGET_COLLECTION_UNIT_BASES * 3 / 2)
                        ? Math.max(1, (n + TARGET_COLLECTION_UNIT_BASES / 2) / TARGET_COLLECTION_UNIT_BASES)
                        : 1;
        int[] sliceLens = new int[numSlices];
        int rem = n;
        int pageBases = Packed2BitSequence.PAGE_BASES;
        for (int i = 0; i < numSlices; i++) {
            int len;
            if (i == numSlices - 1) {
                len = rem;
            } else {
                int avg = rem / (numSlices - i);
                int rounded = (avg + (pageBases >>> 1)) & ~(pageBases - 1);
                len = Math.min(TARGET_COLLECTION_UNIT_BASES, Math.max(pageBases, rounded));
                if (rem - len < 1) {
                    len = rem - (numSlices - 1 - i);
                }
            }
            sliceLens[i] = len;
            rem -= len;
        }

        List<InputUnit> units = new ArrayList<>(numSlices);
        CRC32 crc = new CRC32();
        boolean hasT = false;
        boolean hasU = false;

        try (InputStream in = Files.newInputStream(inputPath)) {
            byte[] ioBuf = new byte[65536];
            int ioPos = 0;
            int ioLim = 0;

            for (int s = 0; s < numSlices; s++) {
                int sLen = sliceLens[s];
                Packed2BitSequence.PageBuilder pb = new Packed2BitSequence.PageBuilder(sLen);
                int k = 0;
                while (k < sLen) {
                    if (ioPos >= ioLim) {
                        ioLim = in.read(ioBuf);
                        if (ioLim <= 0) {
                            return null;
                        }
                        crc.update(ioBuf, 0, ioLim);
                        ioPos = 0;
                    }
                    int b = ioBuf[ioPos++] & 0xFF;
                    int code;
                    switch (b) {
                        case 'A' -> code = 0;
                        case 'C' -> code = 1;
                        case 'G' -> code = 2;
                        case 'T' -> {
                            code = 3;
                            hasT = true;
                            if (hasU) return null;
                        }
                        case 'U' -> {
                            code = 3;
                            hasU = true;
                            if (hasT) return null;
                        }
                        default -> {
                            return null;
                        }
                    }
                    pb.appendBase(code);
                    k++;
                }
                units.add(new InputUnit(pb.build(), null, 0, s == 0 ? 1 : 0));
            }
            if (ioPos < ioLim || in.read() != -1) {
                return null;
            }
        }

        return new ParsedRawInput(n, (int) crc.getValue(), mode, true, hasU && !hasT, 1, units);
    }

    private static List<InputUnit> sliceFlatSequenceIntoUnits(
            Packed2BitSequence flatSeq, int targetBases) {
        int n = flatSeq.length;
        int numSlices = Math.max(1, (n + targetBases / 2) / targetBases);
        List<InputUnit> out = new ArrayList<>(numSlices);
        int cur = 0;
        for (int i = 0; i < numSlices; i++) {
            int len = (n - cur) / (numSlices - i);
            byte[] sub = flatSeq.extractTwoBitSlice(cur, len);
            out.add(new InputUnit(Packed2BitSequence.fromTwoBitBytes(sub, len), null, 0, i == 0 ? 1 : 0));
            cur += len;
        }
        return out;
    }

    private static List<InputUnit> partitionFactorizedCollectionUnits(
            byte[] raw, int targetBasesPerUnit, int threads, boolean bounded) {
        List<int[]> recSlices = splitRawIntoRecordByteRanges(raw, targetBasesPerUnit);
        List<int[]> unitByteRanges = new ArrayList<>();

        int i = 0;
        while (i < recSlices.size()) {
            int startByte = recSlices.get(i)[0];
            int endByte = recSlices.get(i)[1];
            i++;
            // Group small contigs (e.g. draft assembly contigs) up to targetBasesPerUnit
            while (i < recSlices.size()
                    && (endByte - startByte) + (recSlices.get(i)[1] - recSlices.get(i)[0])
                            <= targetBasesPerUnit + 524_288) {
                endByte = recSlices.get(i)[1];
                i++;
            }
            unitByteRanges.add(new int[] {startByte, endByte});
        }

        int numUnits = unitByteRanges.size();
        Packed2BitSequence[] seqs = new Packed2BitSequence[numUnits];
        FactorizedChunk[] chunks = new FactorizedChunk[numUnits];
        int[] recCounts = new int[numUnits];

        int factThreads = bounded ? 1 : Math.min(threads, numUnits);
        if (factThreads <= 1) {
            for (int u = 0; u < numUnits; u++) {
                int[] rng = unitByteRanges.get(u);
                FactorizedChunk chunk = FastaStreamFactorizer.factorize(raw, rng[0], rng[1]);
                recCounts[u] = chunk.getRecordSpans().length;
                seqs[u] = Packed2BitSequence.fromTwoBitBytes(chunk.dnaBases, chunk.dnaLen);
                chunk.dnaBases = null; // Free 1-byte-per-base array immediately!
                chunks[u] = chunk;
            }
        } else {
            ForkJoinPool pool = new ForkJoinPool(factThreads);
            try {
                pool.submit(
                                () ->
                                        IntStream.range(0, numUnits)
                                                .parallel()
                                                .forEach(
                                                        u -> {
                                                            int[] rng = unitByteRanges.get(u);
                                                            FactorizedChunk chunk =
                                                                    FastaStreamFactorizer.factorize(
                                                                            raw, rng[0], rng[1]);
                                                            recCounts[u] =
                                                                    chunk.getRecordSpans().length;
                                                            seqs[u] =
                                                                    Packed2BitSequence
                                                                            .fromTwoBitBytes(
                                                                                    chunk.dnaBases,
                                                                                    chunk.dnaLen);
                                                            chunk.dnaBases = null;
                                                            chunks[u] = chunk;
                                                        }))
                        .join();
            } finally {
                pool.shutdown();
            }
        }

        List<InputUnit> units = new ArrayList<>(numUnits);
        int curRecIdx = 0;
        for (int u = 0; u < numUnits; u++) {
            units.add(new InputUnit(seqs[u], chunks[u], curRecIdx, recCounts[u]));
            curRecIdx += recCounts[u];
        }
        return units;
    }

    /**
     * Splits {@code raw} cleanly at {@code >} FASTA record headers (or line boundaries if a single
     * record exceeds {@code 2 * targetBytes}).
     */
    private static List<int[]> splitRawIntoRecordByteRanges(byte[] raw, int targetBytes) {
        List<int[]> ranges = new ArrayList<>();
        int n = raw.length;
        int curStart = 0;
        boolean atLineStart = true;

        for (int p = 0; p < n; p++) {
            byte b = raw[p];
            if (atLineStart && b == '>' && p > curStart) {
                appendPossiblyLargeRange(raw, curStart, p, targetBytes, ranges);
                curStart = p;
            }
            atLineStart = (b == '\n');
        }
        if (curStart < n) {
            appendPossiblyLargeRange(raw, curStart, n, targetBytes, ranges);
        }
        return ranges;
    }

    private static void appendPossiblyLargeRange(
            byte[] raw, int start, int end, int targetBytes, List<int[]> out) {
        if (end - start <= targetBytes * 3 / 2) {
            out.add(new int[] {start, end, 1});
            return;
        }
        // Split a single giant chromosome record cleanly at newline boundaries near targetBytes
        int cur = start;
        boolean first = true;
        while (cur < end) {
            int aim = cur + targetBytes;
            if (aim >= end - targetBytes / 3) {
                out.add(new int[] {cur, end, first ? 1 : 0});
                break;
            }
            int nl = aim;
            while (nl < end && raw[nl] != '\n') {
                nl++;
            }
            int cut = (nl < end) ? (nl + 1) : aim;
            out.add(new int[] {cur, cut, first ? 1 : 0});
            cur = cut;
            first = false;
        }
    }

    private static Packed2BitSequence tryFromFlatUppercaseAcguOnly(byte[] raw, int len) {
        boolean sawU = false;
        boolean sawT = false;
        byte[] converted = new byte[len];
        for (int i = 0; i < len; i++) {
            int b = raw[i] & 0xFF;
            if (b == 'U') {
                sawU = true;
                converted[i] = 'T';
            } else if (b == 'A' || b == 'C' || b == 'G') {
                converted[i] = (byte) b;
            } else if (b == 'T') {
                sawT = true;
                break;
            } else {
                return null;
            }
        }
        if (!sawU || sawT) return null;
        return Packed2BitSequence.tryFromFlatUppercaseAcgt(converted, len);
    }

    private record LoadedPanel(
            String[] names, Packed2BitSequence[] seqs, byte[] manifestSha256, long totalBases) {}

    private static LoadedPanel loadReferencePanelForCompress(
            Path panelPath,
            Set<String> excludeNames,
            List<String> explicitNames,
            int panelLimit)
            throws IOException {
        List<Path> files = new ArrayList<>();
        if (Files.isRegularFile(panelPath)) {
            files.add(panelPath);
        } else if (Files.isDirectory(panelPath)) {
            if (explicitNames != null) {
                for (String name : explicitNames) {
                    if (excludeNames != null && excludeNames.contains(name)) continue;
                    Path p = panelPath.resolve(name);
                    if (!Files.isRegularFile(p)) {
                        throw new IOException("Explicit panel file missing: " + p);
                    }
                    files.add(p);
                }
            } else {
                try (DirectoryStream<Path> ds = Files.newDirectoryStream(panelPath)) {
                    for (Path p : ds) {
                        if (!Files.isRegularFile(p)) continue;
                        String fn = p.getFileName().toString();
                        if (excludeNames != null && excludeNames.contains(fn)) continue;
                        String low = fn.toLowerCase();
                        if (low.endsWith(".seq")
                                || low.endsWith(".fa")
                                || low.endsWith(".fna")
                                || low.endsWith(".fasta")) {
                            files.add(p);
                        }
                    }
                }
                files.sort(Comparator.comparing(p -> p.getFileName().toString()));
            }
        } else {
            throw new IOException("Reference panel path does not exist: " + panelPath);
        }

        if (panelLimit > 0 && files.size() > panelLimit) {
            files = new ArrayList<>(files.subList(0, panelLimit));
        }

        int m = files.size();
        String[] names = new String[m];
        Packed2BitSequence[] seqs = new Packed2BitSequence[m];
        byte[][] fileHashes = new byte[m][];
        long[] fileByteLens = new long[m];

        final List<Path> finalFiles = files;
        IntStream.range(0, m)
                .parallel()
                .forEach(
                        i -> {
                            try {
                                Path p = finalFiles.get(i);
                                names[i] = p.getFileName().toString();
                                byte[] raw = Files.readAllBytes(p);
                                fileByteLens[i] = raw.length;
                                fileHashes[i] = sha256(raw);
                                seqs[i] = extractDnaFromReferenceBytes(raw);
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });

        long totalBases = 0L;
        for (Packed2BitSequence s : seqs) {
            totalBases += s.length;
        }
        byte[] manifestSha = computeManifestSha256(names, fileByteLens, fileHashes);
        return new LoadedPanel(names, seqs, manifestSha, totalBases);
    }

    private static Packed2BitSequence[] verifyAndLoadPanelForDecode(
            Path panelPath, String[] expectedNames, byte[] expectedManifestSha256)
            throws IOException {
        if (panelPath == null) {
            throw new IllegalStateException(
                    "Mode R archive requires --ref / --panel reference path to decode");
        }
        int m = expectedNames.length;
        Path[] files = new Path[m];
        if (Files.isRegularFile(panelPath)) {
            if (m != 1 || !panelPath.getFileName().toString().equals(expectedNames[0])) {
                throw new IllegalStateException(
                        "Reference file name mismatch: expected "
                                + Arrays.toString(expectedNames)
                                + ", got "
                                + panelPath.getFileName());
            }
            files[0] = panelPath;
        } else if (Files.isDirectory(panelPath)) {
            for (int i = 0; i < m; i++) {
                Path p = panelPath.resolve(expectedNames[i]);
                if (!Files.isRegularFile(p)) {
                    throw new IllegalStateException("Missing required reference panel file: " + p);
                }
                files[i] = p;
            }
        } else {
            throw new IllegalStateException("Reference panel path not found: " + panelPath);
        }

        Packed2BitSequence[] seqs = new Packed2BitSequence[m];
        byte[][] fileHashes = new byte[m][];
        long[] fileByteLens = new long[m];

        IntStream.range(0, m)
                .parallel()
                .forEach(
                        i -> {
                            try {
                                byte[] raw = Files.readAllBytes(files[i]);
                                fileByteLens[i] = raw.length;
                                fileHashes[i] = sha256(raw);
                                seqs[i] = extractDnaFromReferenceBytes(raw);
                            } catch (IOException e) {
                                throw new RuntimeException(e);
                            }
                        });

        byte[] actualManifestSha = computeManifestSha256(expectedNames, fileByteLens, fileHashes);
        if (!MessageDigest.isEqual(expectedManifestSha256, actualManifestSha)) {
            throw new IllegalStateException(
                    "Reference panel SHA-256 manifest mismatch: supplied panel differs from pinned archive panel");
        }
        return seqs;
    }

    private static Packed2BitSequence extractDnaFromReferenceBytes(byte[] raw) {
        Packed2BitSequence flat = Packed2BitSequence.tryFromFlatUppercaseAcgt(raw, raw.length);
        if (flat != null) {
            return flat;
        }
        byte[] twoBit = new byte[raw.length];
        int len = 0;
        boolean inHeader = false;
        for (int i = 0; i < raw.length; i++) {
            int b = raw[i] & 0xFF;
            if (inHeader) {
                if (b == '\n') inHeader = false;
                continue;
            }
            if (b == '>' && (i == 0 || raw[i - 1] == '\n' || raw[i - 1] == '\r')) {
                inHeader = true;
                continue;
            }
            switch (b) {
                case 'A', 'a' -> twoBit[len++] = 0;
                case 'C', 'c' -> twoBit[len++] = 1;
                case 'G', 'g' -> twoBit[len++] = 2;
                case 'T', 't', 'U', 'u' -> twoBit[len++] = 3;
                default -> {}
            }
        }
        return Packed2BitSequence.fromTwoBitBytes(twoBit, len);
    }

    private static byte[] computeManifestSha256(
            String[] names, long[] byteLens, byte[][] fileSha256s) {
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            updateIntBE(md, names.length);
            for (int i = 0; i < names.length; i++) {
                byte[] nb = names[i].getBytes(StandardCharsets.UTF_8);
                updateIntBE(md, nb.length);
                md.update(nb);
                updateLongBE(md, byteLens[i]);
                md.update(fileSha256s[i]);
            }
            return md.digest();
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }
    }

    private static byte[] sha256(byte[] data) {
        try {
            return MessageDigest.getInstance("SHA-256").digest(data);
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }
    }

    private static void updateIntBE(MessageDigest md, int v) {
        md.update((byte) (v >>> 24));
        md.update((byte) (v >>> 16));
        md.update((byte) (v >>> 8));
        md.update((byte) v);
    }

    private static void updateLongBE(MessageDigest md, long v) {
        updateIntBE(md, (int) (v >>> 32));
        updateIntBE(md, (int) v);
    }

    private static byte[] compressPanelNames(String[] names) {
        ByteOut raw = new ByteOut(names.length * 8 + 16);
        byte[] prev = new byte[0];
        for (String s : names) {
            byte[] cur = s.getBytes(StandardCharsets.UTF_8);
            int lcp = 0;
            int lim = Math.min(prev.length, cur.length);
            while (lcp < lim && prev[lcp] == cur[lcp]) {
                lcp++;
            }
            raw.writeVarInt(lcp);
            raw.writeVarInt(cur.length - lcp);
            raw.writeBytes(cur, lcp, cur.length - lcp);
            prev = cur;
        }
        return deflateBytes(raw.toByteArray());
    }

    private static String[] decompressPanelNames(byte[] comp, int count) {
        byte[] raw = inflateBytes(comp, 0, comp.length);
        ByteIn in = new ByteIn(raw);
        String[] names = new String[count];
        byte[] prev = new byte[0];
        for (int i = 0; i < count; i++) {
            int lcp = in.readVarInt();
            int rem = in.readVarInt();
            if (lcp < 0 || rem < 0 || lcp > prev.length) {
                throw new IllegalStateException("Malformed panel name prefix length");
            }
            byte[] cur = new byte[lcp + rem];
            System.arraycopy(prev, 0, cur, 0, lcp);
            byte[] tail = in.readBytes(rem);
            System.arraycopy(tail, 0, cur, lcp, rem);
            names[i] = new String(cur, StandardCharsets.UTF_8);
            prev = cur;
        }
        if (in.hasRemaining()) {
            throw new IllegalStateException("Trailing bytes in panel names stream");
        }
        return names;
    }

    private static byte[] deflateBytes(byte[] raw) {
        Deflater def = new Deflater(Deflater.BEST_COMPRESSION);
        try {
            def.setInput(raw);
            def.finish();
            ByteArrayOutputStream bos = new ByteArrayOutputStream(Math.max(32, raw.length / 2));
            byte[] buf = new byte[4096];
            while (!def.finished()) {
                int n = def.deflate(buf);
                bos.write(buf, 0, n);
            }
            return bos.toByteArray();
        } finally {
            def.end();
        }
    }

    private static byte[] inflateBytes(byte[] src, int off, int len) {
        Inflater inf = new Inflater();
        try {
            inf.setInput(src, off, len);
            ByteArrayOutputStream bos = new ByteArrayOutputStream(Math.max(64, len * 2));
            byte[] buf = new byte[8192];
            while (!inf.finished()) {
                int n = inf.inflate(buf);
                if (n == 0 && (inf.needsInput() || inf.needsDictionary())) {
                    break;
                }
                bos.write(buf, 0, n);
            }
            if (!inf.finished()) {
                throw new IllegalStateException("Truncated Deflate stream in HGS2 metadata");
            }
            return bos.toByteArray();
        } catch (DataFormatException e) {
            throw new IllegalStateException("Corrupt Deflate stream in HGS2 metadata", e);
        } finally {
            inf.end();
        }
    }

    // =========================================================================
    // HGS2 Container Structural Parser & Validator
    // =========================================================================

    private static final class ParsedHgs2 {
        int mode;
        boolean rawFlatAcgt;
        boolean rawRnaU;
        long rawByteLen;
        int rawCrc32;
        int totalRecords;
        int numUnits;

        int panelCount;
        String[] panelNames;
        byte[] panelManifestSha256;

        int[] rankToUnit;
        int[] maxAncestorRank;

        int[] firstRecordIdx;
        int[] recordCount;
        int[] unitMetaOffsets;
        int[] unitMetaLens;
        int[] unitDnaLens;

        int[] mosaicOffsets;
        int[] mosaicLens;

        int totalNoveltyBases;
        int noveltyOffset;
        int noveltyPayloadLen;

        long headerBytes;
        long panelPinningBytes;
        long metadataBytes;
        long mosaicSumBytes;
        long noveltyFramingAndPayloadBytes;

        int crc32UnsignedInt() {
            return rawCrc32;
        }
    }

    private static ParsedHgs2 parseHeaderAndDirectory(byte[] archive) {
        if (!isHgs2Magic(archive) || archive.length < 18) {
            throw new IllegalArgumentException("Not a valid HGS2 archive");
        }
        ByteIn in = new ByteIn(archive);
        in.readBytes(4); // HGS2
        int flags = in.readByte() & 0xFF;
        if ((flags & ~0x0F) != 0) {
            throw new IllegalArgumentException("Unsupported HGS2 flags: 0x" + Integer.toHexString(flags));
        }
        int mode = flags & 3;
        if (mode > MODE_COLLECTION) {
            throw new IllegalArgumentException("Invalid HGS2 mode: " + mode);
        }

        ParsedHgs2 p = new ParsedHgs2();
        p.mode = mode;
        p.rawFlatAcgt = (flags & FLAG_RAW_FLAT_ACGT) != 0;
        p.rawRnaU = (flags & FLAG_RAW_RNA_U) != 0;
        p.rawByteLen = in.readLongLE();
        if (p.rawByteLen < 0) {
            throw new IllegalArgumentException("Negative raw byte length in HGS2 header");
        }
        p.rawCrc32 = in.readIntLE();
        p.totalRecords = in.readVarInt();
        p.numUnits = in.readVarInt();
        if (p.totalRecords < 0 || p.numUnits < 0) {
            throw new IllegalArgumentException("Negative counts in HGS2 header");
        }
        p.headerBytes = in.position();

        if (mode == MODE_REF) {
            int pinStart = in.position();
            p.panelCount = in.readVarInt();
            int compNamesLen = in.readVarInt();
            byte[] compNames = in.readBytes(compNamesLen);
            p.panelNames = decompressPanelNames(compNames, p.panelCount);
            p.panelManifestSha256 = in.readBytes(32);
            p.panelPinningBytes = in.position() - pinStart;
        }

        if (p.numUnits == 0) {
            if (in.hasRemaining()) {
                throw new IllegalArgumentException("Unexpected trailing bytes in empty HGS2 archive");
            }
            return p;
        }

        int metaStart = in.position();
        p.rankToUnit = new int[p.numUnits];
        p.maxAncestorRank = new int[p.numUnits];
        if (mode == MODE_COLLECTION && p.numUnits > 1) {
            boolean[] seen = new boolean[p.numUnits];
            for (int r = 0; r < p.numUnits; r++) {
                int u = in.readVarInt();
                if (u < 0 || u >= p.numUnits || seen[u]) {
                    throw new IllegalArgumentException("Invalid rankToUnit permutation in HGS2");
                }
                seen[u] = true;
                p.rankToUnit[r] = u;
            }
            for (int r = 0; r < p.numUnits; r++) {
                int anc = in.readVarInt();
                if (anc < 0 || anc > r) {
                    throw new IllegalArgumentException("Non-causal maxAncestorRank in HGS2");
                }
                p.maxAncestorRank[r] = anc;
            }
        } else {
            for (int u = 0; u < p.numUnits; u++) {
                p.rankToUnit[u] = u;
                p.maxAncestorRank[u] = 0;
            }
        }

        if (!p.rawFlatAcgt) {
            p.firstRecordIdx = new int[p.numUnits];
            p.recordCount = new int[p.numUnits];
            p.unitMetaOffsets = new int[p.numUnits];
            p.unitMetaLens = new int[p.numUnits];
            for (int u = 0; u < p.numUnits; u++) {
                p.firstRecordIdx[u] = in.readVarInt();
                p.recordCount[u] = in.readVarInt();
                int mLen = in.readVarInt();
                p.unitMetaOffsets[u] = in.position();
                p.unitMetaLens[u] = mLen;
                in.readBytes(mLen);
            }
        } else {
            p.unitDnaLens = new int[p.numUnits];
            for (int u = 0; u < p.numUnits; u++) {
                p.unitDnaLens[u] = in.readVarInt();
            }
        }
        p.metadataBytes = in.position() - metaStart;

        int mosStart = in.position();
        p.mosaicOffsets = new int[p.numUnits];
        p.mosaicLens = new int[p.numUnits];
        for (int r = 0; r < p.numUnits; r++) {
            int mLen = in.readVarInt();
            p.mosaicOffsets[r] = in.position();
            p.mosaicLens[r] = mLen;
            in.readBytes(mLen);
        }
        p.mosaicSumBytes = in.position() - mosStart;

        int novStart = in.position();
        p.totalNoveltyBases = in.readVarInt();
        p.noveltyPayloadLen = in.readVarInt();
        p.noveltyOffset = in.position();
        in.readBytes(p.noveltyPayloadLen);
        if (in.hasRemaining()) {
            throw new IllegalArgumentException(
                    "Unexpected trailing bytes at end of HGS2 archive: " + in.remaining());
        }
        p.noveltyFramingAndPayloadBytes = in.position() - novStart;
        return p;
    }
}
