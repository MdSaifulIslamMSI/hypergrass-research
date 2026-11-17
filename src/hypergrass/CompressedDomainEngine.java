package hypergrass;

import hypergrass.EimdlMosaicCodec.MosaicInspectResult;
import hypergrass.EimdlMosaicCodec.MosaicInspectStats;
import hypergrass.EimdlMosaicCodec.SpanEvent;
import hypergrass.EimdlMosaicCodec.VariantRecord;
import hypergrass.Hgs3OmicCodec.CdsSummary;
import hypergrass.Hgs3OmicCodec.ParsedHgs3Fastq;
import hypergrass.HyperGraSSFusionEngine.ParsedHgs2;
import hypergrass.VarIntCodec.ByteIn;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Compressed-Domain Genomic Compute Engine (Pure Java 25 LTS counterpart to Rust {@code compute.rs}).
 *
 * <p>Performs evolutionary inference ({@code analyze}), zero-FASTA VCF v4.2 variant extraction ({@code
 * variants}), exact k-mer search ({@code kmer}), and random-access region slicing ({@code query})
 * directly over {@code .hgs2} and {@code .hgs3} compressed archives.
 */
public final class CompressedDomainEngine {

    private CompressedDomainEngine() {}

    public record AnalyzeReport(
            String format,
            boolean isFastq,
            int mode,
            long rawBytes,
            long totalBases,
            long noveltyBases,
            long panelMatchBases,
            long selfMatchBases,
            long tandemMatchBases,
            long literalBases,
            long panelTokens,
            long snps,
            long transitions,
            long transversions,
            long wobbleMod3Hits,
            long splices,
            long insertedBases,
            long deletedBases,
            long recombinationSwitches,
            String[] donorNames,
            long[] donorMatchBases,
            int unitCount,
            CdsSummary cds,
            double seconds) {
        public double tiTv() {
            if (transversions == 0L) return transitions > 0L ? Double.POSITIVE_INFINITY : 0.0;
            return (double) transitions / transversions;
        }

        public double wobbleMod3Pct() {
            return snps == 0L ? 0.0 : (wobbleMod3Hits * 100.0) / snps;
        }

        public double panelConservationPct() {
            return totalBases == 0L ? 0.0 : (panelMatchBases * 100.0) / totalBases;
        }

        public String formatSummary() {
            StringBuilder sb = new StringBuilder();
            sb.append(
                    String.format(
                            Locale.ROOT,
                            "=== EIMDL Compressed-Domain Evolutionary Report ===%n"
                                    + "Container     : %s (mode=%d, units=%d)%n"
                                    + "Raw Bytes     : %d%n"
                                    + "Total Bases   : %d bp%n"
                                    + "Panel Copy    : %d bp (%.3f%%, %d tracts)%n"
                                    + "Self / Tandem : %d bp / %d bp%n"
                                    + "Literals      : %d bp (novelty=%d)%n"
                                    + "SNPs          : %d (Ti=%d, Tv=%d, Ti/Tv=%.4f)%n"
                                    + "Codon Wobble  : %d/%d (%.2f%% mod-3 spaced)%n"
                                    + "Indel Splices : %d (INS=%d bp, DEL=%d bp)%n"
                                    + "Switches      : %d inter-donor recombination switches%n"
                                    + "Elapsed       : %.4f s%n",
                            format,
                            mode,
                            unitCount,
                            rawBytes,
                            totalBases,
                            panelMatchBases,
                            panelConservationPct(),
                            panelTokens,
                            selfMatchBases,
                            tandemMatchBases,
                            literalBases,
                            noveltyBases,
                            snps,
                            transitions,
                            transversions,
                            tiTv(),
                            wobbleMod3Hits,
                            snps,
                            wobbleMod3Pct(),
                            splices,
                            insertedBases,
                            deletedBases,
                            recombinationSwitches,
                            seconds));
            if (cds != null) {
                sb.append(
                        String.format(
                                Locale.ROOT,
                                "Proteomics    : %d ORFs (>=30aa), longest=%d aa, GC=%.1f%%%n",
                                cds.orfCount(),
                                cds.longestOrfAa(),
                                cds.gcPermille() * 0.1));
            }
            return sb.toString();
        }

        public String toJson() {
            double tv = Double.isFinite(tiTv()) ? tiTv() : 0.0;
            return String.format(
                    Locale.ROOT,
                    "{\"format\":\"%s\",\"is_fastq\":%b,\"mode\":%d,\"raw_bytes\":%d,"
                            + "\"total_bases\":%d,\"panel_match_bases\":%d,\"self_match_bases\":%d,"
                            + "\"tandem_match_bases\":%d,\"literal_bases\":%d,\"panel_tokens\":%d,"
                            + "\"snps\":%d,\"transitions\":%d,\"transversions\":%d,\"ti_tv\":%.6f,"
                            + "\"wobble_mod3_hits\":%d,\"wobble_mod3_pct\":%.4f,\"splices\":%d,"
                            + "\"inserted_bases\":%d,\"deleted_bases\":%d,"
                            + "\"recombination_switches\":%d,\"elapsed_seconds\":%.6f}",
                    format,
                    isFastq,
                    mode,
                    rawBytes,
                    totalBases,
                    panelMatchBases,
                    selfMatchBases,
                    tandemMatchBases,
                    literalBases,
                    panelTokens,
                    snps,
                    transitions,
                    transversions,
                    tv,
                    wobbleMod3Hits,
                    wobbleMod3Pct(),
                    splices,
                    insertedBases,
                    deletedBases,
                    recombinationSwitches,
                    seconds);
        }
    }

    public record VcfResult(
            String vcfText,
            int totalVariants,
            int snvs,
            int insertions,
            int deletions,
            int complexIndels,
            int literalRunsSkipped,
            double seconds) {}

    public record KmerHit(int unitIndex, long tgtPos, int strand) {}

    public record KmerReport(
            int queryLen,
            List<KmerHit> hits,
            int exactZeroTargetSpans,
            int patchedOrLocalSpans,
            double seconds) {}

    private static final class LoadedContainer {
        String formatName;
        boolean isFastq;
        int mode;
        long rawBytes;
        String[] panelNames;
        Packed2BitSequence refBank;
        int[] refOffsets;
        int refCount;
        int[] rankToUnit;
        int[] maxAncestorRank;
        byte[][] rankMosaics;
        byte[][] rankNovelties;
        CdsSummary cds;
    }

    private static LoadedContainer loadContainer(Path archivePath, Path panelPath)
            throws IOException {
        byte[] bytes = Files.readAllBytes(archivePath);
        if (bytes.length < 5) {
            throw new IllegalArgumentException("Archive truncated");
        }
        int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
        boolean bounded = Runtime.getRuntime().maxMemory() <= 320L * 1024 * 1024;
        LoadedContainer c = new LoadedContainer();

        if (Hgs3OmicCodec.isHgs3Magic(bytes)
                && (bytes[4] & HyperGraSSFusionEngine.FLAG3_FASTQ) != 0) {
            ParsedHgs3Fastq p = Hgs3OmicCodec.parseFastqContainer(bytes, true);
            c.formatName = "HGS3-FASTQ";
            c.isFastq = true;
            c.mode = p.mode;
            c.rawBytes = p.rawBytes;
            c.panelNames = p.panelNames != null ? p.panelNames : new String[0];
            if (p.mode == HyperGraSSFusionEngine.MODE_REF && p.panelCount > 0) {
                Packed2BitSequence[] seqs =
                        HyperGraSSFusionEngine.verifyAndLoadPanelForDecode(
                                panelPath, p.panelNames, p.panelSha);
                c.refCount = seqs.length;
                c.refOffsets = new int[c.refCount + 1];
                c.refBank = Packed2BitSequence.concatenate(seqs, c.refOffsets);
            }
            if (p.rawBytes > 0 && p.totalRecords > 0) {
                c.rankToUnit = new int[] {0};
                c.maxAncestorRank = new int[] {0};
                c.rankMosaics = new byte[][] {p.mosaicBytes};
                c.rankNovelties =
                        new byte[][] {HgsAnsCoder.ansNoveltyDecode(p.ansNovComp, p.novLen, threads)};
            } else {
                c.rankToUnit = new int[0];
                c.maxAncestorRank = new int[0];
                c.rankMosaics = new byte[0][];
                c.rankNovelties = new byte[0][];
            }
            if (p.cdsComp != null && p.cdsComp.length > 0) {
                c.cds = Hgs3OmicCodec.deserializeCdsSummary(p.cdsComp);
            }
            return c;
        }

        ParsedHgs2 hdr = HyperGraSSFusionEngine.parseHeaderAndDirectory(bytes, true);
        c.formatName = hdr.isHgs3 ? "HGS3-FASTA" : "HGS2";
        c.isFastq = false;
        c.mode = hdr.mode;
        c.rawBytes = hdr.rawByteLen;
        c.panelNames = hdr.panelNames != null ? hdr.panelNames : new String[0];
        if (hdr.mode == HyperGraSSFusionEngine.MODE_REF && hdr.panelCount > 0) {
            Packed2BitSequence[] seqs =
                    HyperGraSSFusionEngine.verifyAndLoadPanelForDecode(
                            panelPath, hdr.panelNames, hdr.panelManifestSha256);
            c.refCount = seqs.length;
            c.refOffsets = new int[c.refCount + 1];
            c.refBank = Packed2BitSequence.concatenate(seqs, c.refOffsets);
        }
        if (hdr.cdsLen > 0) {
            byte[] cb = Arrays.copyOfRange(bytes, hdr.cdsOffset, hdr.cdsOffset + hdr.cdsLen);
            c.cds = Hgs3OmicCodec.deserializeCdsSummary(cb);
        }

        int n = hdr.numUnits;
        c.rankToUnit = hdr.rankToUnit != null ? hdr.rankToUnit : new int[0];
        c.maxAncestorRank = hdr.maxAncestorRank != null ? hdr.maxAncestorRank : new int[0];
        c.rankMosaics = new byte[n][];
        c.rankNovelties = new byte[n][];
        if (n == 0) {
            return c;
        }

        byte[] hgs3Nov = null;
        FastNoveltyCoder.NoveltyReader novReader = null;
        if (hdr.isHgs3) {
            byte[] novComp =
                    Arrays.copyOfRange(
                            bytes, hdr.noveltyOffset, hdr.noveltyOffset + hdr.noveltyPayloadLen);
            hgs3Nov = HgsAnsCoder.ansNoveltyDecode(novComp, hdr.totalNoveltyBases, threads);
        } else {
            ByteIn novIn = new ByteIn(bytes, hdr.noveltyOffset, hdr.noveltyPayloadLen);
            novReader =
                    FastNoveltyCoder.openNoveltyReader(novIn, hdr.totalNoveltyBases, threads, bounded);
        }

        int h3Off = 0;
        for (int r = 0; r < n; r++) {
            c.rankMosaics[r] =
                    Arrays.copyOfRange(
                            bytes, hdr.mosaicOffsets[r], hdr.mosaicOffsets[r] + hdr.mosaicLens[r]);
            int[] tlNl = EimdlMosaicCodec.peekLengths(c.rankMosaics[r]);
            int nlen = tlNl[1];
            if (hgs3Nov != null) {
                c.rankNovelties[r] = Arrays.copyOfRange(hgs3Nov, h3Off, h3Off + nlen);
                h3Off += nlen;
            } else {
                c.rankNovelties[r] = novReader.nextSlice(nlen);
            }
        }
        return c;
    }

    private record UnitInspect(int unitIndex, MosaicInspectResult res, byte[] novelty) {}

    private static List<UnitInspect> inspectAllRanks(
            LoadedContainer c, boolean collectVariants, boolean collectSpans) {
        int n = c.rankMosaics.length;
        List<UnitInspect> list = new ArrayList<>(n);
        if (c.mode == HyperGraSSFusionEngine.MODE_COLLECTION && n > 1) {
            int maxAnc = 0;
            for (int a : c.maxAncestorRank) if (a > maxAnc) maxAnc = a;
            List<Packed2BitSequence> ancSeqs = new ArrayList<>(maxAnc);
            int[] ancOffsets = new int[maxAnc + 1];
            Packed2BitSequence ancBank = null;
            for (int r = 0; r < n; r++) {
                int anc = c.maxAncestorRank[r];
                if (anc > 0) {
                    Packed2BitSequence[] sub = ancSeqs.subList(0, anc).toArray(new Packed2BitSequence[0]);
                    ancBank = Packed2BitSequence.concatenate(sub, ancOffsets);
                }
                MosaicInspectResult rep =
                        EimdlMosaicCodec.inspectMosaic(
                                c.rankMosaics[r],
                                c.rankNovelties[r],
                                ancBank,
                                ancOffsets,
                                anc,
                                collectVariants,
                                collectSpans);
                if (r < maxAnc) {
                    byte[] seq =
                            EimdlMosaicCodec.decodeMosaicToTwoBit(
                                    c.rankMosaics[r], c.rankNovelties[r], ancBank, ancOffsets, anc);
                    ancSeqs.add(Packed2BitSequence.fromTwoBitBytes(seq, seq.length));
                }
                list.add(new UnitInspect(c.rankToUnit[r], rep, c.rankNovelties[r]));
            }
        } else {
            for (int u = 0; u < n; u++) {
                MosaicInspectResult rep =
                        EimdlMosaicCodec.inspectMosaic(
                                c.rankMosaics[u],
                                c.rankNovelties[u],
                                c.refBank,
                                c.refOffsets,
                                c.refCount,
                                collectVariants,
                                collectSpans);
                list.add(new UnitInspect(u, rep, c.rankNovelties[u]));
            }
        }
        return list;
    }

    public static AnalyzeReport analyze(Path archive, Path panel) throws IOException {
        long t0 = System.nanoTime();
        LoadedContainer c = loadContainer(archive, panel);
        List<UnitInspect> units = inspectAllRanks(c, false, c.panelNames.length > 0);
        long totBases = 0L;
        long novBases = 0L;
        long pBases = 0L;
        long sBases = 0L;
        long tBases = 0L;
        long lBases = 0L;
        long pTok = 0L;
        long snps = 0L;
        long ti = 0L;
        long tv = 0L;
        long wob = 0L;
        long sp = 0L;
        long ins = 0L;
        long del = 0L;
        long sw = 0L;
        long[] donorBases = new long[c.panelNames.length];

        for (UnitInspect ui : units) {
            MosaicInspectStats st = ui.res.stats();
            totBases += st.targetLen();
            novBases += st.noveltyLen();
            pBases += st.panelMatchBases();
            sBases += st.selfMatchBases();
            tBases += st.tandemMatchBases();
            lBases += st.literalBases();
            pTok += st.panelTokens();
            snps += st.snps();
            ti += st.transitions();
            tv += st.transversions();
            wob += st.wobbleMod3Hits();
            sp += st.splices();
            ins += st.insertedBases();
            del += st.deletedBases();
            sw += st.recombinationSwitches();
            for (SpanEvent e : ui.res.spans()) {
                if (e.kind() == 1 && e.donorId() >= 0 && e.donorId() < donorBases.length) {
                    donorBases[e.donorId()] += e.tgtLen();
                }
            }
        }
        return new AnalyzeReport(
                c.formatName,
                c.isFastq,
                c.mode,
                c.rawBytes,
                totBases,
                novBases,
                pBases,
                sBases,
                tBases,
                lBases,
                pTok,
                snps,
                ti,
                tv,
                wob,
                sp,
                ins,
                del,
                sw,
                c.panelNames,
                donorBases,
                units.size(),
                c.cds,
                (System.nanoTime() - t0) / 1e9);
    }

    private static final char[] ACGT = {'A', 'C', 'G', 'T'};

    private static void appendBases(StringBuilder sb, byte[] bases) {
        for (byte b : bases) sb.append(ACGT[b & 3]);
    }

    public static VcfResult extractVariantsVcf(Path archive, Path panel, Path outVcf)
            throws IOException {
        long t0 = System.nanoTime();
        LoadedContainer c = loadContainer(archive, panel);
        List<UnitInspect> units = inspectAllRanks(c, true, false);
        StringBuilder vcf = new StringBuilder(4096);
        vcf.append("##fileformat=VCFv4.2\n");
        vcf.append("##source=EIMDL-HyperGraSS-CompressedDomain-Java25\n");
        vcf.append("##INFO=<ID=TGT_POS,Number=1,Type=Integer,Description=\"0-based coordinate in target archive\">\n");
        vcf.append("##INFO=<ID=UNIT,Number=1,Type=Integer,Description=\"Target archive unit index\">\n");
        vcf.append("##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Donor alignment strand (+ or -)\">\n");
        vcf.append("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant class (SNV, INS, DEL, COMPLEX)\">\n");
        vcf.append("##INFO=<ID=TI,Number=1,Type=Integer,Description=\"1 if transition SNP else 0\">\n");
        vcf.append("##INFO=<ID=MOD3,Number=1,Type=Integer,Description=\"1 if SNP step % 3 == 0 codon-wobble spaced\">\n");
        for (int i = 0; i < c.panelNames.length; i++) {
            int dlen = c.refOffsets != null ? (c.refOffsets[i + 1] - c.refOffsets[i]) : 0;
            vcf.append("##contig=<ID=").append(c.panelNames[i]).append(",length=").append(dlen).append(">\n");
        }
        vcf.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

        int total = 0;
        int snvs = 0;
        int ins = 0;
        int del = 0;
        int cplx = 0;
        int litSkip = 0;
        for (UnitInspect ui : units) {
            for (VariantRecord v : ui.res.variants()) {
                if (v.kind() == 4) {
                    litSkip++;
                    continue;
                }
                int did = v.donorId();
                if (did < 0 || did >= c.refCount) continue;
                int dOff = c.refOffsets[did];
                int dLen = c.refOffsets[did + 1] - dOff;
                String chrom = c.panelNames[did];

                int fwdStart =
                        (v.strand() == 0)
                                ? (int) v.refPos()
                                : (int) (dLen - v.refPos() - v.refLen());
                int fwdEnd = fwdStart + v.refLen();
                if (fwdStart < 0 || fwdEnd > dLen) continue;

                byte[] fwdRef = new byte[v.refLen()];
                for (int i = 0; i < v.refLen(); i++) {
                    fwdRef[i] = (byte) c.refBank.baseAt(dOff + fwdStart + i);
                }
                byte[] fwdAlt;
                if (v.altLen() == 0 || v.altBases() == null) {
                    fwdAlt = new byte[0];
                } else if (v.strand() == 0) {
                    fwdAlt = v.altBases();
                } else {
                    byte[] src = v.altBases();
                    fwdAlt = new byte[src.length];
                    for (int i = 0; i < src.length; i++) {
                        fwdAlt[i] = (byte) ((src[src.length - 1 - i] ^ 3) & 3);
                    }
                }

                long pos1;
                if (fwdRef.length == 0 || fwdAlt.length == 0) {
                    if (fwdStart > 0) {
                        byte anc = (byte) c.refBank.baseAt(dOff + fwdStart - 1);
                        byte[] nr = new byte[fwdRef.length + 1];
                        nr[0] = anc;
                        System.arraycopy(fwdRef, 0, nr, 1, fwdRef.length);
                        byte[] na = new byte[fwdAlt.length + 1];
                        na[0] = anc;
                        System.arraycopy(fwdAlt, 0, na, 1, fwdAlt.length);
                        fwdRef = nr;
                        fwdAlt = na;
                        pos1 = fwdStart;
                    } else if (fwdEnd < dLen) {
                        byte anc = (byte) c.refBank.baseAt(dOff + fwdEnd);
                        byte[] nr = Arrays.copyOf(fwdRef, fwdRef.length + 1);
                        nr[fwdRef.length] = anc;
                        byte[] na = Arrays.copyOf(fwdAlt, fwdAlt.length + 1);
                        na[fwdAlt.length] = anc;
                        fwdRef = nr;
                        fwdAlt = na;
                        pos1 = 1L;
                    } else {
                        pos1 = fwdStart + 1L;
                    }
                } else {
                    pos1 = fwdStart + 1L;
                }

                String kStr =
                        switch (v.kind()) {
                            case 0 -> {
                                snvs++;
                                yield "SNV";
                            }
                            case 1 -> {
                                ins++;
                                yield "INS";
                            }
                            case 2 -> {
                                del++;
                                yield "DEL";
                            }
                            default -> {
                                cplx++;
                                yield "COMPLEX";
                            }
                        };
                total++;
                vcf.append(chrom).append('\t').append(pos1).append("\t.\t");
                appendBases(vcf, fwdRef);
                vcf.append('\t');
                appendBases(vcf, fwdAlt);
                vcf.append("\t.\tPASS\tTGT_POS=")
                        .append(v.tgtPos())
                        .append(";UNIT=")
                        .append(ui.unitIndex)
                        .append(";STRAND=")
                        .append(v.strand() == 0 ? '+' : '-')
                        .append(";TYPE=")
                        .append(kStr)
                        .append(";TI=")
                        .append(v.isTransition())
                        .append(";MOD3=")
                        .append(v.prevMod3())
                        .append('\n');
            }
        }
        String text = vcf.toString();
        if (outVcf != null) {
            Files.writeString(outVcf, text, StandardCharsets.UTF_8);
        }
        return new VcfResult(
                text, total, snvs, ins, del, cplx, litSkip, (System.nanoTime() - t0) / 1e9);
    }

    private static byte[] parseQuery(String q) {
        byte[] out = new byte[q.length()];
        for (int i = 0; i < q.length(); i++) {
            char ch = q.charAt(i);
            out[i] =
                    switch (ch) {
                        case 'A', 'a' -> 0;
                        case 'C', 'c' -> 1;
                        case 'G', 'g' -> 2;
                        case 'T', 't', 'U', 'u' -> 3;
                        default -> throw new IllegalArgumentException("Non-ACGT base in query: " + ch);
                    };
        }
        if (out.length == 0) throw new IllegalArgumentException("Empty k-mer query");
        return out;
    }

    private static byte[] revComp(byte[] q) {
        byte[] r = new byte[q.length];
        for (int i = 0; i < q.length; i++) r[i] = (byte) ((q[q.length - 1 - i] ^ 3) & 3);
        return r;
    }

    private static int findSpanContaining(List<SpanEvent> spans, int curPos) {
        int lo = 0;
        int hi = spans.size();
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (spans.get(mid).tgtPos() <= curPos) lo = mid + 1;
            else hi = mid;
        }
        return Math.max(0, lo - 1);
    }

    private static int findSnpAlt(List<VariantRecord> snps, int pos) {
        int lo = 0;
        int hi = snps.size() - 1;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            long tp = snps.get(mid).tgtPos();
            if (tp == pos) return snps.get(mid).altBases()[0] & 3;
            if (tp < pos) lo = mid + 1;
            else hi = mid - 1;
        }
        return -1;
    }

    private static int donorBaseStrand(
            Packed2BitSequence bank, int[] offsets, int did, int strand, int mr) {
        int dOff = offsets[did];
        int dLen = offsets[did + 1] - dOff;
        return bank.baseAtStrand(dOff, dLen, strand, mr);
    }

    private static byte resolveUnitBaseFromSpans(
            int pos,
            List<SpanEvent> spans,
            List<VariantRecord> snps,
            byte[] novelty,
            Packed2BitSequence bank,
            int[] offsets) {
        int cur = pos;
        int flip = 0;
        for (int depth = 0; depth < 64; depth++) {
            SpanEvent sp = spans.get(findSpanContaining(spans, cur));
            int off = cur - (int) sp.tgtPos();
            switch (sp.kind()) {
                case 0 -> {
                    return (byte) ((novelty[(int) sp.refPos() + off] & 3) ^ flip);
                }
                case 1 -> {
                    if (sp.hasEdits()) {
                        int alt = findSnpAlt(snps, cur);
                        if (alt >= 0) return (byte) (alt ^ flip);
                    }
                    return (byte)
                            (donorBaseStrand(
                                            bank, offsets, sp.donorId(), sp.strand(), (int) sp.refPos() + off)
                                    ^ flip);
                }
                case 3 -> {
                    int period = Math.max(1, (int) sp.refPos());
                    cur = (int) sp.tgtPos() - period + (off % period);
                }
                default -> {
                    int mr = (int) sp.refPos() + off;
                    if (sp.strand() == 0) {
                        cur = mr;
                    } else {
                        cur = (int) sp.tgtPos() - 1 - mr;
                        flip ^= 3;
                    }
                }
            }
        }
        return 0;
    }

    private static byte resolveSpanPos(
            int pos,
            SpanEvent sp,
            List<SpanEvent> spans,
            List<VariantRecord> snps,
            byte[] novelty,
            Packed2BitSequence bank,
            int[] offsets) {
        int off = pos - (int) sp.tgtPos();
        return switch (sp.kind()) {
            case 0 -> (byte) (novelty[(int) sp.refPos() + off] & 3);
            case 1 -> {
                if (sp.hasEdits()) {
                    int alt = findSnpAlt(snps, pos);
                    if (alt >= 0) yield (byte) alt;
                }
                yield (byte)
                        donorBaseStrand(
                                bank, offsets, sp.donorId(), sp.strand(), (int) sp.refPos() + off);
            }
            default -> resolveUnitBaseFromSpans(pos, spans, snps, novelty, bank, offsets);
        };
    }

    private static boolean matchAt(byte[] arr, int off, byte[] pat) {
        for (int i = 0; i < pat.length; i++) {
            if (arr[off + i] != pat[i]) return false;
        }
        return true;
    }

    private static void scanWindow(
            byte[] win,
            int winLen,
            long basePos,
            int unitIdx,
            byte[] qFwd,
            byte[] qRc,
            List<KmerHit> hits) {
        int k = qFwd.length;
        if (winLen < k) return;
        for (int i = 0; i <= winLen - k; i++) {
            if (matchAt(win, i, qFwd)) hits.add(new KmerHit(unitIdx, basePos + i, 0));
            if (qRc != null && matchAt(win, i, qRc)) hits.add(new KmerHit(unitIdx, basePos + i, 1));
        }
    }

    public static KmerReport searchKmers(
            Path archive, Path panel, String queryAscii, boolean bothStrands) throws IOException {
        long t0 = System.nanoTime();
        byte[] qFwd = parseQuery(queryAscii);
        byte[] qRc = bothStrands ? revComp(qFwd) : null;
        int k = qFwd.length;

        LoadedContainer c = loadContainer(archive, panel);
        List<KmerHit> hits = new ArrayList<>();
        int zeroTargetSpans = 0;
        int patchedSpans = 0;

        if (c.mode == HyperGraSSFusionEngine.MODE_COLLECTION && c.rankMosaics.length > 1) {
            int n = c.rankMosaics.length;
            List<Packed2BitSequence> dec = new ArrayList<>(n);
            int[] offs = new int[n + 1];
            for (int r = 0; r < n; r++) {
                int anc = c.maxAncestorRank[r];
                Packed2BitSequence bank =
                        anc > 0
                                ? Packed2BitSequence.concatenate(
                                        dec.subList(0, anc).toArray(new Packed2BitSequence[0]), offs)
                                : null;
                byte[] dna =
                        EimdlMosaicCodec.decodeMosaicToTwoBit(
                                c.rankMosaics[r], c.rankNovelties[r], bank, offs, anc);
                scanWindow(dna, dna.length, 0L, c.rankToUnit[r], qFwd, qRc, hits);
                patchedSpans++;
                dec.add(Packed2BitSequence.fromTwoBitBytes(dna, dna.length));
            }
            dedupSortHits(hits);
            return new KmerReport(
                    k, hits, zeroTargetSpans, patchedSpans, (System.nanoTime() - t0) / 1e9);
        }

        List<UnitInspect> units = inspectAllRanks(c, true, true);
        for (UnitInspect ui : units) {
            boolean hasSelfEdits = false;
            for (SpanEvent sp : ui.res.spans()) {
                if (sp.kind() == 2 && sp.hasEdits()) {
                    hasSelfEdits = true;
                    break;
                }
            }
            if (hasSelfEdits) {
                byte[] dna =
                        EimdlMosaicCodec.decodeMosaicToTwoBit(
                                c.rankMosaics[ui.unitIndex],
                                ui.novelty,
                                c.refBank,
                                c.refOffsets,
                                c.refCount);
                scanWindow(dna, dna.length, 0L, ui.unitIndex, qFwd, qRc, hits);
                patchedSpans += ui.res.spans().size();
                continue;
            }

            List<VariantRecord> snps = new ArrayList<>();
            for (VariantRecord v : ui.res.variants()) {
                if (v.kind() == 0) snps.add(v);
            }
            snps.sort(Comparator.comparingLong(VariantRecord::tgtPos));

            int snpCursor = 0;
            byte[] spanBuf = new byte[0];
            byte[] tail = new byte[k > 1 ? k - 1 : 0];
            int tailLen = 0;
            byte[] jBuf = new byte[k * 2];

            for (SpanEvent sp : ui.res.spans()) {
                int spStart = (int) sp.tgtPos();
                int spLen = (int) sp.tgtLen();
                int spEnd = spStart + spLen;
                while (snpCursor < snps.size() && snps.get(snpCursor).tgtPos() < spStart) {
                    snpCursor++;
                }

                if (tailLen > 0 && k > 1) {
                    int takeR = Math.min(spLen, k - 1);
                    System.arraycopy(tail, 0, jBuf, 0, tailLen);
                    for (int off = 0; off < takeR; off++) {
                        jBuf[tailLen + off] =
                                resolveSpanPos(
                                        spStart + off,
                                        sp,
                                        ui.res.spans(),
                                        snps,
                                        ui.novelty,
                                        c.refBank,
                                        c.refOffsets);
                    }
                    int jLen = tailLen + takeR;
                    int jStart = spStart - tailLen;
                    if (jLen >= k) {
                        int lim = Math.min(tailLen, jLen - k + 1);
                        for (int i = 0; i < lim; i++) {
                            if (matchAt(jBuf, i, qFwd)) hits.add(new KmerHit(ui.unitIndex, jStart + i, 0));
                            if (qRc != null && matchAt(jBuf, i, qRc)) {
                                hits.add(new KmerHit(ui.unitIndex, jStart + i, 1));
                            }
                        }
                    }
                }

                if (spLen >= k) {
                    if (sp.kind() == 0) {
                        zeroTargetSpans++;
                        int rPos = (int) sp.refPos();
                        for (int off = 0; off <= spLen - k; off++) {
                            if (matchAt(ui.novelty, rPos + off, qFwd)) {
                                hits.add(new KmerHit(ui.unitIndex, spStart + off, 0));
                            }
                            if (qRc != null && matchAt(ui.novelty, rPos + off, qRc)) {
                                hits.add(new KmerHit(ui.unitIndex, spStart + off, 1));
                            }
                        }
                    } else if (sp.kind() == 1 && !sp.hasEdits()) {
                        zeroTargetSpans++;
                        int rs = (int) sp.refPos();
                        for (int off = 0; off <= spLen - k; off++) {
                            boolean mf = true;
                            boolean mr = qRc != null;
                            for (int j = 0; j < k; j++) {
                                int db =
                                        donorBaseStrand(
                                                c.refBank, c.refOffsets, sp.donorId(), sp.strand(), rs + off + j);
                                if (db != qFwd[j]) mf = false;
                                if (mr && db != qRc[j]) mr = false;
                                if (!mf && !mr) break;
                            }
                            if (mf) hits.add(new KmerHit(ui.unitIndex, spStart + off, 0));
                            if (mr) hits.add(new KmerHit(ui.unitIndex, spStart + off, 1));
                        }
                    } else {
                        patchedSpans++;
                        if (spanBuf.length < spLen) spanBuf = new byte[spLen];
                        if (sp.kind() == 1) {
                            int rs = (int) sp.refPos();
                            for (int off = 0; off < spLen; off++) {
                                spanBuf[off] =
                                        (byte)
                                                donorBaseStrand(
                                                        c.refBank,
                                                        c.refOffsets,
                                                        sp.donorId(),
                                                        sp.strand(),
                                                        rs + off);
                            }
                            int sc = snpCursor;
                            while (sc < snps.size() && snps.get(sc).tgtPos() < spEnd) {
                                VariantRecord sv = snps.get(sc++);
                                spanBuf[(int) sv.tgtPos() - spStart] = sv.altBases()[0];
                            }
                        } else {
                            for (int off = 0; off < spLen; off++) {
                                spanBuf[off] =
                                        resolveUnitBaseFromSpans(
                                                spStart + off,
                                                ui.res.spans(),
                                                snps,
                                                ui.novelty,
                                                c.refBank,
                                                c.refOffsets);
                            }
                        }
                        scanWindow(spanBuf, spLen, spStart, ui.unitIndex, qFwd, qRc, hits);
                    }
                } else {
                    zeroTargetSpans++;
                }

                if (k > 1) {
                    int keepPrev = Math.max(0, (k - 1) - spLen);
                    byte[] nextTail = new byte[k - 1];
                    int nt = 0;
                    if (keepPrev > 0 && tailLen > 0) {
                        int use = Math.min(keepPrev, tailLen);
                        System.arraycopy(tail, tailLen - use, nextTail, 0, use);
                        nt = use;
                    }
                    int takeThis = Math.min(spLen, (k - 1) - nt);
                    int startOff = spLen - takeThis;
                    for (int off = startOff; off < spLen; off++) {
                        nextTail[nt++] =
                                resolveSpanPos(
                                        spStart + off,
                                        sp,
                                        ui.res.spans(),
                                        snps,
                                        ui.novelty,
                                        c.refBank,
                                        c.refOffsets);
                    }
                    System.arraycopy(nextTail, 0, tail, 0, nt);
                    tailLen = nt;
                }
            }
        }
        dedupSortHits(hits);
        return new KmerReport(
                k, hits, zeroTargetSpans, patchedSpans, (System.nanoTime() - t0) / 1e9);
    }

    private static void dedupSortHits(List<KmerHit> hits) {
        hits.sort(
                Comparator.comparingInt(KmerHit::unitIndex)
                        .thenComparingLong(KmerHit::tgtPos)
                        .thenComparingInt(KmerHit::strand));
        for (int i = hits.size() - 1; i > 0; i--) {
            if (hits.get(i).equals(hits.get(i - 1))) hits.remove(i);
        }
    }

    public static String querySlice(
            Path archive, Path panel, int unitIdx, int start, int end) throws IOException {
        LoadedContainer c = loadContainer(archive, panel);
        int n = c.rankMosaics.length;
        if (unitIdx < 0 || unitIdx >= n) {
            throw new IndexOutOfBoundsException("Unit index " + unitIdx + " out of 0.." + n);
        }
        byte[] twoBit;
        if (c.mode == HyperGraSSFusionEngine.MODE_COLLECTION && n > 1) {
            List<Packed2BitSequence> dec = new ArrayList<>(n);
            int[] offs = new int[n + 1];
            twoBit = new byte[0];
            for (int r = 0; r < n; r++) {
                int anc = c.maxAncestorRank[r];
                Packed2BitSequence bank =
                        anc > 0
                                ? Packed2BitSequence.concatenate(
                                        dec.subList(0, anc).toArray(new Packed2BitSequence[0]), offs)
                                : null;
                byte[] seq =
                        EimdlMosaicCodec.decodeMosaicToTwoBit(
                                c.rankMosaics[r], c.rankNovelties[r], bank, offs, anc);
                if (c.rankToUnit[r] == unitIdx) {
                    twoBit = seq;
                    break;
                }
                dec.add(Packed2BitSequence.fromTwoBitBytes(seq, seq.length));
            }
        } else {
            twoBit =
                    EimdlMosaicCodec.decodeMosaicToTwoBit(
                            c.rankMosaics[unitIdx],
                            c.rankNovelties[unitIdx],
                            c.refBank,
                            c.refOffsets,
                            c.refCount);
        }
        int s = Math.max(0, Math.min(start, twoBit.length));
        int e = Math.max(s, Math.min(end, twoBit.length));
        StringBuilder sb = new StringBuilder(e - s);
        for (int i = s; i < e; i++) sb.append(ACGT[twoBit[i] & 3]);
        return sb.toString();
    }
}
