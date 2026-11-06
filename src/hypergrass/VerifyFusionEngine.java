package hypergrass;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * End-to-End Integration, Random-Access Query, Byte-Accounting, Pinning Security, &amp; Corruption
 * Rejection Suite for {@link HyperGraSSFusionEngine} ({@code HGS2}) and {@link HyperGraSS}.
 */
public final class VerifyFusionEngine {

    private static int checks = 0;

    private VerifyFusionEngine() {}

    public static void main(String[] args) throws Exception {
        System.out.println("=== HyperGraSS Fusion Engine (HGS2) Integration Verification ===");
        Path tmp = Files.createTempDirectory("verify-hgs2-");
        try {
            verifyModeSingleAndQuery(tmp);
            verifyModeRefAndCryptographicPinning(tmp);
            verifyModeCollectionPrimAndQuery(tmp);
            verifyCorruptArchiveRejection(tmp);
            verifyEimdlCliRoundTrip(tmp);
        } finally {
            deleteRecursively(tmp);
        }
        System.out.printf("=== ALL HGS2 FUSION ENGINE SUITES PASSED (%d checks) ===%n", checks);
    }

    private static void verifyModeSingleAndQuery(Path tmp) throws Exception {
        String fasta =
                ">chr1 synthetic chromosome 1\n"
                        + repeatWithMutations("ACGTACGTAGCTAGCTAGGCTAGCTAGCTAGC", 400, 11L)
                        + "\n>chr2lowercaseAndNs\r\n"
                        + "acgtacgtNNNNNNNNNNRYKMSWacgtacgt\r\n"
                        + ">chr3_tail\n"
                        + repeatWithMutations("TTAGGGTTAGGGTTAGGG", 300, 22L);
        byte[] raw = fasta.getBytes(StandardCharsets.UTF_8);
        Path in = tmp.resolve("single.fa");
        Path arch = tmp.resolve("single.hgs");
        Path out = tmp.resolve("single.out.fa");
        Files.write(in, raw);

        HyperGraSS.CompressOptions opts = HyperGraSS.CompressOptions.defaults();
        opts.fusionMode = true;
        opts.singleMode = true;
        HyperGraSS.compressFile(in, arch, opts);
        check(HyperGraSSFusionEngine.isHgs2File(arch), "Expected HGS2 magic");

        HyperGraSS.decompressFile(arch, out);
        check(Arrays.equals(raw, Files.readAllBytes(out)), "Mode S file round-trip mismatch");

        // Exact byte accounting sum == archive length
        byte[] archBytes = Files.readAllBytes(arch);
        HyperGraSSFusionEngine.Hgs2Accounting acc =
                HyperGraSSFusionEngine.accountArchive(archBytes);
        long sum =
                acc.containerHeaderBytes()
                        + acc.panelPinningBytes()
                        + acc.metadataBytes()
                        + acc.mosaicBytes()
                        + acc.noveltyBytes();
        check(sum == archBytes.length, "Mode S byte accounting sum != archive length");

        // Inspect & random-access query
        HyperGraSS.ArchiveMetadata meta = HyperGraSS.inspect(arch);
        check(meta.totalRawBytes == raw.length && meta.totalRecords == 3, "Inspect metadata mismatch");

        byte[] qChr2 = HyperGraSS.query(arch, 1, 8, 24);
        check(
                "NNNNNNNNNNRYKMSW".equals(new String(qChr2, StandardCharsets.UTF_8)),
                "Random-access query slice mismatch on chr2: " + new String(qChr2));
        System.out.println("  [PASS] Mode S round-trip, byte accounting, inspect, and sub-slice query");
    }

    private static void verifyModeRefAndCryptographicPinning(Path tmp) throws Exception {
        Path panelDir = tmp.resolve("panel");
        Files.createDirectories(panelDir);
        Random rng = new Random(33L);
        byte[] ref0 = randomDna(40_000, rng);
        byte[] ref1 = mutateDna(ref0, 80, rng);
        byte[] target = mutateDna(ref1, 25, rng);

        Files.write(panelDir.resolve("00.seq"), ref0);
        Files.write(panelDir.resolve("01.seq"), ref1);
        Files.write(panelDir.resolve("02.seq"), target);

        Path targetPath = panelDir.resolve("02.seq");
        Path arch = tmp.resolve("target_ref.hgs");
        Path out = tmp.resolve("target_ref.out");

        HyperGraSS.CompressOptions opts = HyperGraSS.CompressOptions.defaults();
        opts.fusionMode = true;
        opts.referencePanel = panelDir;
        opts.excludePanelNames.add("02.seq");
        HyperGraSS.compressFile(targetPath, arch, opts);

        // Must compress 40 kbp 25-SNP twin into a tiny delta (< 400 B!)
        long archSize = Files.size(arch);
        check(archSize < 400, "Mode R archive unexpectedly large: " + archSize);

        // Valid panel decodes identically
        HyperGraSS.decompressFile(arch, out, panelDir);
        check(Arrays.equals(target, Files.readAllBytes(out)), "Mode R round-trip mismatch");

        // Missing --ref must fail loudly
        boolean rejectedMissingRef = false;
        try {
            HyperGraSS.decompressFile(arch, out, null);
        } catch (Exception expected) {
            rejectedMissingRef = true;
        }
        check(rejectedMissingRef, "Mode R must reject missing reference panel");

        // Tampering with 1 byte in 01.seq must fail SHA-256 manifest verification!
        byte[] tampered = Arrays.copyOf(ref1, ref1.length);
        tampered[500] = (tampered[500] == 'A') ? (byte) 'G' : (byte) 'A';
        Files.write(panelDir.resolve("01.seq"), tampered);
        boolean rejectedTamper = false;
        try {
            HyperGraSS.decompressFile(arch, out, panelDir);
        } catch (Exception expected) {
            rejectedTamper = true;
        }
        check(rejectedTamper, "Mode R must reject 1-byte tampered reference panel");

        // Restore valid 01.seq and verify deletion of 00.seq is also rejected
        Files.write(panelDir.resolve("01.seq"), ref1);
        Files.delete(panelDir.resolve("00.seq"));
        boolean rejectedDeleted = false;
        try {
            HyperGraSS.decompressFile(arch, out, panelDir);
        } catch (Exception expected) {
            rejectedDeleted = true;
        }
        check(rejectedDeleted, "Mode R must reject missing panel file");
        Files.write(panelDir.resolve("00.seq"), ref0);

        System.out.printf(
                "  [PASS] Mode R delta (%d B for 40 kbp) + cryptographic SHA-256 manifest tamper rejection%n",
                archSize);
    }

    private static void verifyModeCollectionPrimAndQuery(Path tmp) throws Exception {
        Random rng = new Random(44L);
        // Build 6 strains with shuffled clade order so Prim MST reorders ranks != input order!
        byte[] cladeA = randomDna(35_000, rng);
        byte[] cladeB = randomDna(35_000, rng);
        String[][] records = {
            {"strain0_cladeA1", new String(mutateDna(cladeA, 40, rng), StandardCharsets.UTF_8)},
            {"strain1_cladeB1", new String(mutateDna(cladeB, 45, rng), StandardCharsets.UTF_8)},
            {"strain2_cladeA2", new String(mutateDna(cladeA, 15, rng), StandardCharsets.UTF_8)},
            {"strain3_cladeB2", new String(mutateDna(cladeB, 12, rng), StandardCharsets.UTF_8)},
            {"strain4_cladeA0", new String(cladeA, StandardCharsets.UTF_8)},
            {"strain5_cladeB0", new String(cladeB, StandardCharsets.UTF_8)}
        };

        StringBuilder sb = new StringBuilder();
        for (String[] rec : records) {
            sb.append('>').append(rec[0]).append('\n');
            String seq = rec[1];
            for (int i = 0; i < seq.length(); i += 80) {
                sb.append(seq, i, Math.min(seq.length(), i + 80)).append('\n');
            }
        }
        byte[] rawFasta = sb.toString().getBytes(StandardCharsets.UTF_8);
        Path colIn = tmp.resolve("pangenome6.fa");
        Path colArch = tmp.resolve("pangenome6.hgs");
        Path colOut = tmp.resolve("pangenome6.out.fa");
        Files.write(colIn, rawFasta);

        HyperGraSSFusionEngine.FusionOptions fOpts = new HyperGraSSFusionEngine.FusionOptions();
        fOpts.mode = HyperGraSSFusionEngine.MODE_COLLECTION;
        // Force per-record units on small synthetic test by compressing via multi-unit partition
        byte[] archBytes = HyperGraSSFusionEngine.compress(rawFasta, fOpts);
        Files.write(colArch, archBytes);

        HyperGraSS.decompressFile(colArch, colOut);
        check(
                Arrays.equals(rawFasta, Files.readAllBytes(colOut)),
                "Mode C multi-FASTA round-trip mismatch");

        // Query every strain record & sub-slice
        for (int r = 0; r < records.length; r++) {
            byte[] slice = HyperGraSS.query(colArch, r, 100, 250);
            String expected = records[r][1].substring(100, 250);
            check(
                    expected.equals(new String(slice, StandardCharsets.UTF_8)),
                    "Mode C query mismatch on record " + r);
        }
        System.out.printf(
                "  [PASS] Mode C pangenome round-trip (%d B -> %d B) + random-access query%n",
                rawFasta.length, archBytes.length);

        // Multi-unit flat .seq collection (> 6.75 Mbp so sliceFlatSequenceIntoUnits triggers)
        byte[] ancBlock = randomDna(300_000, new Random(909L));
        byte[] flatCollection = new byte[7_200_000];
        Random mr = new Random(910L);
        byte[] bases = {'A', 'C', 'G', 'T'};
        int writePos = 0;
        for (int copy = 0; copy < 24; copy++) {
            for (int i = 0; i < ancBlock.length && writePos < flatCollection.length; i++) {
                byte b = ancBlock[i];
                // Inject occasional adjacent-SNP + micro-deletion boundary patterns
                if (i == 12000 + copy * 13 && copy % 2 == 1) {
                    flatCollection[writePos++] = bases[(indexOf(b) + 1) & 3];
                    if (writePos < flatCollection.length) {
                        flatCollection[writePos++] = bases[(indexOf(ancBlock[i + 1]) + 2) & 3];
                    }
                    i += 5; // 4-bp micro-deletion right after 2 adjacent SNPs!
                    continue;
                }
                if (mr.nextInt(256) == 0) {
                    b = bases[(indexOf(b) + 2) & 3];
                }
                flatCollection[writePos++] = b;
            }
        }
        while (writePos < flatCollection.length) {
            flatCollection[writePos++] = 'A';
        }
        byte[] flatColArch = HyperGraSSFusionEngine.compress(flatCollection, fOpts);
        byte[] flatColDec = HyperGraSSFusionEngine.decompress(flatColArch, null);
        check(
                Arrays.equals(flatCollection, flatColDec),
                "Mode C multi-unit flat .seq round-trip mismatch");
        Path flatColPath = tmp.resolve("flat7m.hgs2");
        Files.write(flatColPath, flatColArch);
        byte[] qFlat = HyperGraSSFusionEngine.query(flatColPath, 0, 5_000_000, 5_000_250, null);
        check(
                Arrays.equals(Arrays.copyOfRange(flatCollection, 5_000_000, 5_000_250), qFlat),
                "Mode C multi-unit flat .seq cross-unit query mismatch");
        System.out.printf(
                "  [PASS] Mode C multi-unit flat .seq round-trip (7,200,000 B -> %d B)%n",
                flatColArch.length);
    }

    private static void verifyCorruptArchiveRejection(Path tmp) throws Exception {
        byte[] raw = randomDna(25_000, new Random(55L));
        HyperGraSSFusionEngine.FusionOptions opts = new HyperGraSSFusionEngine.FusionOptions();
        opts.mode = HyperGraSSFusionEngine.MODE_SINGLE;
        byte[] valid = HyperGraSSFusionEngine.compress(raw, opts);

        // 1. Truncated archives at 15 cut points must all fail
        for (int cut = 0; cut < valid.length; cut += Math.max(1, valid.length / 15)) {
            byte[] trunc = Arrays.copyOf(valid, cut);
            boolean failed = false;
            try {
                HyperGraSSFusionEngine.decompress(trunc, null);
            } catch (Throwable t) {
                failed = true;
            }
            check(failed, "Truncated archive at " + cut + " was not rejected");
        }

        // 2. Flipping CRC-32 byte in header must fail
        byte[] badCrc = Arrays.copyOf(valid, valid.length);
        badCrc[14] ^= 0x55;
        boolean crcFailed = false;
        try {
            HyperGraSSFusionEngine.decompress(badCrc, null);
        } catch (Throwable t) {
            crcFailed = true;
        }
        check(crcFailed, "Corrupted CRC-32 was not rejected");

        // 3. Invalid reserved flag bits must fail
        byte[] badFlags = Arrays.copyOf(valid, valid.length);
        badFlags[4] |= 0x40;
        boolean flagsFailed = false;
        try {
            HyperGraSSFusionEngine.decompress(badFlags, null);
        } catch (Throwable t) {
            flagsFailed = true;
        }
        check(flagsFailed, "Invalid HGS2 header flags were not rejected");

        System.out.println("  [PASS] Truncated, bad-CRC, and malformed HGS2 archives cleanly rejected");
    }

    private static void verifyEimdlCliRoundTrip(Path tmp) throws Exception {
        byte[] raw = randomDna(30_000, new Random(66L));
        Path in = tmp.resolve("cli_in.seq");
        Path arch = tmp.resolve("cli.hgs2");
        Path out = tmp.resolve("cli_out.seq");
        Files.write(in, raw);

        EimdlCli.main(
                new String[] {
                    "compress", "--target", in.toString(), "--out", arch.toString(), "--threads", "2"
                });
        EimdlCli.main(
                new String[] {"decompress", "--archive", arch.toString(), "--out", out.toString()});
        check(Arrays.equals(raw, Files.readAllBytes(out)), "EimdlCli round-trip mismatch");
        System.out.println("  [PASS] EimdlCli JSON telemetry adapter round-trip");
    }

    private static String repeatWithMutations(String motif, int reps, long seed) {
        Random rng = new Random(seed);
        StringBuilder sb = new StringBuilder(motif.length() * reps);
        char[] bases = {'A', 'C', 'G', 'T'};
        for (int r = 0; r < reps; r++) {
            for (int i = 0; i < motif.length(); i++) {
                char c = motif.charAt(i);
                if (rng.nextDouble() < 0.01) c = bases[rng.nextInt(4)];
                sb.append(c);
            }
            if ((r + 1) % 3 == 0) sb.append('\n');
        }
        return sb.toString();
    }

    private static byte[] randomDna(int len, Random rng) {
        byte[] out = new byte[len];
        byte[] map = {'A', 'C', 'G', 'T'};
        for (int i = 0; i < len; i++) out[i] = map[rng.nextInt(4)];
        return out;
    }

    private static byte[] mutateDna(byte[] src, int numSnps, Random rng) {
        byte[] out = Arrays.copyOf(src, src.length);
        byte[] map = {'A', 'C', 'G', 'T'};
        for (int i = 0; i < numSnps; i++) {
            int p = rng.nextInt(out.length);
            byte orig = out[p];
            byte next = map[(rng.nextInt(3) + 1 + indexOf(orig)) & 3];
            out[p] = next;
        }
        return out;
    }

    private static int indexOf(byte b) {
        return switch (b) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            default -> 3;
        };
    }

    private static void check(boolean cond, String msg) {
        if (!cond) throw new AssertionError(msg);
        checks++;
    }

    private static void deleteRecursively(Path root) throws Exception {
        if (root == null || !Files.exists(root)) return;
        try (var walk = Files.walk(root)) {
            walk.sorted((a, b) -> b.compareTo(a))
                    .forEach(
                            p -> {
                                try {
                                    Files.deleteIfExists(p);
                                } catch (Exception ignored) {
                                }
                            });
        }
    }
}
