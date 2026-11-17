package hypergrass;

import hypergrass.CompressedDomainEngine.AnalyzeReport;
import hypergrass.CompressedDomainEngine.KmerHit;
import hypergrass.CompressedDomainEngine.KmerReport;
import hypergrass.CompressedDomainEngine.VcfResult;
import hypergrass.Hgs3OmicCodec.Hgs3Options;
import hypergrass.Hgs3OmicCodec.Hgs3Report;
import hypergrass.Hgs3OmicCodec.OmicFormat;
import hypergrass.HyperGraSSFusionEngine.FusionOptions;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * End-to-end verification suite for the three Next-Gen EIMDL Paradigm Pillars in Java 25 LTS:
 * (1) Compressed-domain analytics, VCFv4.2 variant extraction, exact k-mer search &amp; slicing;
 * (2) Unified {@code .hgs3} lossless multi-omic container (FASTQ + FASTA + CDS ORF digest);
 * (3) 4-stripe interleaved Byte-rANS entropy coder with Law IV noise bounds and Phred conditioning.
 */
public final class VerifyParadigm {

    private VerifyParadigm() {}

    private static byte[] lcgDna(long seed, int len) {
        long x = seed;
        byte[] out = new byte[len];
        for (int i = 0; i < len; i++) {
            x = x * 6364136223846793005L + 1442695040888963407L;
            out[i] = (byte) ((x >>> 32) & 3L);
        }
        return out;
    }

    private static String toAscii(byte[] twoBit) {
        char[] map = {'A', 'C', 'G', 'T'};
        StringBuilder sb = new StringBuilder(twoBit.length);
        for (byte b : twoBit) sb.append(map[b & 3]);
        return sb.toString();
    }

    private static void check(boolean cond, String msg) {
        if (!cond) throw new AssertionError(msg);
    }

    private static void verifyAns4() {
        byte[] motif = lcgDna(0x31415926L, 384);
        byte[] rep = new byte[384 * 250];
        for (int i = 0; i < 250; i++) System.arraycopy(motif, 0, rep, i * 384, 384);

        byte[] comp = HgsAnsCoder.ans4EncodeDnaSlice(rep, 0, rep.length);
        byte[] dec = HgsAnsCoder.ans4DecodeDnaSlice(comp, 0, comp.length, rep.length);
        check(Arrays.equals(dec, rep), "ANS4 DNA roundtrip mismatch");
        double repBpb = (comp.length * 8.0) / rep.length;
        check(repBpb < 1.2, "Expected ANS4 repeat bpb < 1.2, got " + repBpb);

        byte[] big = new byte[1_152_000];
        for (int i = 0; i < 3000; i++) System.arraycopy(motif, 0, big, i * 384, 384);
        byte[] b1 = HgsAnsCoder.ansNoveltyEncode(big, 1);
        byte[] b2 = HgsAnsCoder.ansNoveltyEncode(big, 2);
        check(Arrays.equals(b1, b2), "ANS novelty wire format must be thread-count deterministic");
        byte[] bd = HgsAnsCoder.ansNoveltyDecode(b2, big.length, 2);
        check(Arrays.equals(bd, big), "ANS multi-segment novelty decode mismatch");

        byte[] noise = lcgDna(0xFEEDFACECAFE0001L, 120_000);
        byte[] nc = HgsAnsCoder.ansNoveltyEncode(noise, 2);
        byte[] nd = HgsAnsCoder.ansNoveltyDecode(nc, noise.length, 2);
        check(Arrays.equals(nd, noise), "Law IV noise decode mismatch");
        double noiseBpb = (nc.length * 8.0) / noise.length;
        check(noiseBpb <= 2.02, "Expected Law IV ceiling <= 2.02 bpb, got " + noiseBpb);

        byte[] homo = HgsAnsCoder.computeHomoBins(Arrays.copyOf(rep, 40_000));
        byte[] mm = new byte[40_000];
        byte[] quals = new byte[40_000];
        for (int i = 0; i < quals.length; i++) {
            int cycle = i % 150;
            if (i % 37 == 0) {
                mm[i] = 1;
                quals[i] = '#';
            } else if (homo[i] >= 2) {
                quals[i] = '5';
            } else if (cycle > 135) {
                quals[i] = '?';
            } else if (i % 997 == 0) {
                quals[i] = '~'; // nibble escape (> 96)
            } else {
                quals[i] = 'I';
            }
        }
        byte[] qc = HgsAnsCoder.ans4EncodeQual(quals, mm, homo, 150);
        byte[] qd = HgsAnsCoder.ans4DecodeQual(qc, mm, homo, quals.length, 150);
        check(Arrays.equals(qd, quals), "ANS4 conditioned quality roundtrip mismatch");
        check(qc.length < quals.length / 4, "Expected high Phred compression");
    }

    private static void verifyHgs3FastqAndFasta(Path tmp) throws IOException {
        Path pdir = tmp.resolve("fq_panel");
        Files.createDirectories(pdir);
        byte[] refDna = lcgDna(0x42424242L, 6000);
        String refAscii = toAscii(refDna);
        String refFa = ">donor_ref\n" + refAscii + "\n";
        Files.writeString(pdir.resolve("ref0.fa"), refFa);

        StringBuilder fqSb = new StringBuilder();
        for (int r = 0; r < 80; r++) {
            String id = String.format("SRR999001.%d inst:1:1101:%d:%d length=150", r, 1000 + r, 2000 + r * 3);
            char[] read = refAscii.substring(r * 50, r * 50 + 150).toCharArray();
            char[] qual = new char[150];
            Arrays.fill(qual, 'I');
            if (r % 5 == 0) {
                read[30] = read[30] == 'T' ? 'A' : 'T';
                qual[30] = '!';
            }
            if (r % 9 == 0) {
                read[10] = 'N';
                read[11] = 'a';
            }
            String eol = (r % 2 == 0) ? "\n" : "\r\n";
            fqSb.append('@').append(id).append(eol);
            fqSb.append(read).append(eol);
            fqSb.append('+');
            if (r == 3) fqSb.append(id);
            else if (r == 7) fqSb.append("custom_comment");
            fqSb.append(eol);
            fqSb.append(qual);
            if (r + 1 < 80) fqSb.append(eol);
        }
        byte[] rawFq = fqSb.toString().getBytes(StandardCharsets.US_ASCII);
        Path fqIn = tmp.resolve("reads.fastq");
        Path arc3 = tmp.resolve("reads.hgs3");
        Path fqOut = tmp.resolve("decoded.fastq");
        Files.write(fqIn, rawFq);

        Hgs3Options opts = new Hgs3Options();
        opts.format = OmicFormat.AUTO;
        opts.emitCdsDigest = true;
        opts.fusion.referencePanel = pdir;
        Hgs3Report rep = Hgs3OmicCodec.compressOmics(fqIn, arc3, opts);
        check(rep.isFastq(), "Expected HGS3 FASTQ mode");
        check(rep.archiveBytes() < rawFq.length / 3L, "Expected >3x FASTQ compression");

        Hgs3OmicCodec.decompressOmics(arc3, fqOut, pdir, 2);
        check(Arrays.equals(Files.readAllBytes(fqOut), rawFq), "HGS3 FASTQ byte-exact roundtrip failed");

        AnalyzeReport ana = CompressedDomainEngine.analyze(arc3, pdir);
        check("HGS3-FASTQ".equals(ana.format()), "Expected HGS3-FASTQ analyze format");
        check(ana.panelMatchBases() > 8000, "Expected >8000 panel-matched read bases");

        // Panel tamper rejection
        byte[] bad = refFa.getBytes(StandardCharsets.US_ASCII);
        bad[15] = (byte) (bad[15] == 'A' ? 'C' : 'A');
        Files.write(pdir.resolve("ref0.fa"), bad);
        boolean rejected = false;
        try {
            Hgs3OmicCodec.decompressOmics(arc3, fqOut, pdir, 2);
        } catch (IllegalStateException ex) {
            rejected = true;
        }
        check(rejected, "Expected SHA-256 panel tamper rejection");

        // HGS3 FASTA + CDS ORF check
        String cdsFa = ">orf_test\nATG" + "GCT".repeat(40) + "TAA\n";
        Path faIn = tmp.resolve("cds.fa");
        Path faArc = tmp.resolve("cds.hgs3");
        Path faOut = tmp.resolve("cds_dec.fa");
        Files.writeString(faIn, cdsFa);
        Hgs3Options faOpts = new Hgs3Options();
        faOpts.format = OmicFormat.FASTA;
        faOpts.emitCdsDigest = true;
        faOpts.fusion.mode = HyperGraSSFusionEngine.MODE_SINGLE;
        faOpts.fusion.singleGenome = true;
        Hgs3Report faRep = Hgs3OmicCodec.compressOmics(faIn, faArc, faOpts);
        check(!faRep.isFastq() && faRep.cds() != null && faRep.cds().orfCount() >= 1, "Expected CDS ORF count >= 1");
        Hgs3OmicCodec.decompressOmics(faArc, faOut, null, 2);
        check(Files.readString(faOut).equals(cdsFa), "HGS3 FASTA roundtrip mismatch");
    }

    private static void verifyCompressedDomainGenomics(Path tmp) throws IOException {
        Path pdir = tmp.resolve("cd_panel");
        Files.createDirectories(pdir);
        byte[] d0 = lcgDna(0xAAAA1111L, 10_000);
        byte[] d1 = lcgDna(0xBBBB2222L, 10_000);
        Files.writeString(pdir.resolve("d0.fa"), ">donor0\n" + toAscii(d0) + "\n");
        Files.writeString(pdir.resolve("d1.fa"), ">donor1\n" + toAscii(d1) + "\n");

        byte[] segA = Arrays.copyOfRange(d0, 0, 3000);
        for (int i = 120; i < 2880; i += 45) {
            segA[i] ^= (byte) ((i % 90 == 0) ? 2 : 1); // alternating Ti/Tv at mod-3 step 45 bp
        }
        byte[] insLit = {0, 1, 2, 3, 2, 1, 0, 3};
        byte[] segC = Arrays.copyOfRange(d0, 3005, 5500);
        byte[] segD = new byte[3000];
        for (int i = 0; i < 3000; i++) segD[i] = (byte) ((d1[4999 - i] ^ 3) & 3);
        for (int i = 120; i < 2880; i += 48) {
            segD[i] ^= (byte) ((i % 96 == 0) ? 2 : 3);
        }

        byte[] tgt = new byte[segA.length + insLit.length + segC.length + segD.length];
        int cur = 0;
        for (byte[] part : new byte[][] {segA, insLit, segC, segD}) {
            System.arraycopy(part, 0, tgt, cur, part.length);
            cur += part.length;
        }
        String tgtAscii = toAscii(tgt);
        Path tgtFa = tmp.resolve("tgt.fa");
        Files.writeString(tgtFa, ">recomb\n" + tgtAscii + "\n");

        Path arc2 = tmp.resolve("tgt.hgs2");
        Path arc3 = tmp.resolve("tgt.hgs3");
        FusionOptions fopts = new FusionOptions();
        fopts.referencePanel = pdir;
        HyperGraSSFusionEngine.compressFile(tgtFa, arc2, fopts);
        Hgs3Options h3opts = new Hgs3Options();
        h3opts.format = OmicFormat.FASTA;
        h3opts.fusion = fopts;
        Hgs3OmicCodec.compressOmics(tgtFa, arc3, h3opts);

        for (Path arc : new Path[] {arc2, arc3}) {
            AnalyzeReport ana = CompressedDomainEngine.analyze(arc, pdir);
            check(ana.totalBases() == tgt.length, "Target length mismatch");
            check(ana.snps() >= 4, "Expected SNPs >= 4, got:\n" + ana.formatSummary());
            check(ana.recombinationSwitches() >= 1, "Expected recombination switches >= 1");

            Path vcfOut = tmp.resolve("calls.vcf");
            VcfResult vcf = CompressedDomainEngine.extractVariantsVcf(arc, pdir, vcfOut);
            check(vcf.vcfText().contains("##fileformat=VCFv4.2"), "Expected VCFv4.2 header");
            check(vcf.snvs() >= 4, "Expected VCF SNVs >= 4");
            check(vcf.vcfText().contains("STRAND=-"), "Expected RC donor strand variant call");

            String[] queries = {
                tgtAscii.substring(100, 115),
                tgtAscii.substring(498, 512),
                tgtAscii.substring(2995, 3015),
                tgtAscii.substring(6295, 6312),
                "ACGTACGT"
            };
            boolean[] bothFlags = {false, true, true, true, true};
            for (int qi = 0; qi < queries.length; qi++) {
                String q = queries[qi];
                boolean both = bothFlags[qi];
                String qRc = revCompStr(q);
                List<KmerHit> expected = new ArrayList<>();
                for (int p = 0; p <= tgtAscii.length() - q.length(); p++) {
                    if (tgtAscii.startsWith(q, p)) expected.add(new KmerHit(0, p, 0));
                    if (both && tgtAscii.startsWith(qRc, p)) expected.add(new KmerHit(0, p, 1));
                }
                expected.sort(Comparator.comparingLong(KmerHit::tgtPos).thenComparingInt(KmerHit::strand));
                for (int i = expected.size() - 1; i > 0; i--) {
                    if (expected.get(i).equals(expected.get(i - 1))) expected.remove(i);
                }
                KmerReport got = CompressedDomainEngine.searchKmers(arc, pdir, q, both);
                check(got.hits().equals(expected), "Java compressed k-mer mismatch on " + q);
            }

            String sliced = CompressedDomainEngine.querySlice(arc, pdir, 0, 2990, 3030);
            check(sliced.equals(tgtAscii.substring(2990, 3030)), "Query slice mismatch");
        }
    }

    private static String revCompStr(String s) {
        StringBuilder sb = new StringBuilder(s.length());
        for (int i = s.length() - 1; i >= 0; i--) {
            sb.append(
                    switch (s.charAt(i)) {
                        case 'A' -> 'T';
                        case 'C' -> 'G';
                        case 'G' -> 'C';
                        default -> 'A';
                    });
        }
        return sb.toString();
    }

    private static void deleteTree(Path dir) throws IOException {
        if (!Files.exists(dir)) return;
        try (var walk = Files.walk(dir)) {
            walk.sorted(Comparator.reverseOrder())
                    .forEach(
                            p -> {
                                try {
                                    Files.deleteIfExists(p);
                                } catch (IOException ignored) {
                                }
                            });
        }
    }

    public static void main(String[] args) throws Exception {
        Path tmp = Files.createTempDirectory("eimdl_paradigm_java_");
        try {
            verifyAns4();
            verifyHgs3FastqAndFasta(tmp);
            verifyCompressedDomainGenomics(tmp);
            System.out.println("VerifyParadigm PASS: 4x-rANS + HGS3 multi-omic + compressed-domain genomics");
        } finally {
            deleteTree(tmp);
        }
    }
}
