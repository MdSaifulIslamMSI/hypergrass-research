package hypergrass;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.HexFormat;
import java.util.Random;

/**
 * Exhaustive SHA-256 byte-exact lossless verification suite covering all biological, formatting,
 * streaming (--lean), and O(1) random-access query edge cases.
 */
public final class VerifyAllEdgeCases {

    private static final Path GRASS_SRC = Path.of(
            "C:/Users/mdsai/.gemini/antigravity/brain/92f6b277-dea3-42a9-b123-90c65fe25a1c/scratch/GraSS/src");

    public static void main(String[] args) throws Exception {
        System.out.println("=== HyperGraSS Exhaustive SHA-256 Lossless Verification Suite ===");
        Path tmpDir = Files.createTempDirectory("hypergrass_verify_");
        int passed = 0;

        // 1. Trivial & boundary files
        passed += checkRoundTrip("EmptyFile (0 B)", new byte[0], tmpDir);
        passed += checkRoundTrip("SingleLF (1 B)", "\n".getBytes(StandardCharsets.UTF_8), tmpDir);
        passed += checkRoundTrip("SingleCRLF (2 B)", "\r\n".getBytes(StandardCharsets.UTF_8), tmpDir);
        passed += checkRoundTrip("Raw4BasesNoEOL", "ACGT".getBytes(StandardCharsets.UTF_8), tmpDir);
        passed += checkRoundTrip("Raw5BasesLF", "ACGTA\n".getBytes(StandardCharsets.UTF_8), tmpDir);
        passed += checkRoundTrip("Raw5BasesCRLF", "ACGTA\r\n".getBytes(StandardCharsets.UTF_8), tmpDir);

        // 2. Bundled GraSS files (In.fa, BuEb, GCF_001884535.1.fna CRLF + Linux LF)
        byte[] inFa = Files.readAllBytes(GRASS_SRC.resolve("In.fa"));
        passed += checkRoundTrip("GraSS In.fa (RNA + zy)", inFa, tmpDir);

        byte[] buebRaw = Files.readAllBytes(GRASS_SRC.resolve("BuEb"));
        String buebStr = new String(buebRaw, StandardCharsets.US_ASCII).replace("\r", "").replace("\n", "");
        byte[] bueb18940 = buebStr.getBytes(StandardCharsets.US_ASCII);
        passed += checkRoundTrip("Table 2 BuEb (18,940 B exact)", bueb18940, tmpDir);

        byte[] gcfCrlf = Files.readAllBytes(GRASS_SRC.resolve("GCF_001884535.1.fna"));
        passed += checkRoundTrip("Table 1 GCF_001884535.1 (CRLF disk)", gcfCrlf, tmpDir);

        byte[] gcfLf = new String(gcfCrlf, StandardCharsets.US_ASCII).replace("\r\n", "\n")
                .getBytes(StandardCharsets.US_ASCII);
        passed += checkRoundTrip("Table 1 GCF_001884535.1 (Linux LF 50,654 B)", gcfLf, tmpDir);

        // 3. Adversarial torture test (cases where GraSS fails/corrupts!)
        StringBuilder adv = new StringBuilder();
        adv.append(">empty_leading_record\n");
        adv.append(">chimeric_DNA_RNA_ragged_mixed_EOL\r\n");
        adv.append("acgtACGTuuuuUUUUttttTTTTRYSWKMBDHVryswkmbdhvZzYy-*.\n");
        adv.append("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\r\n");
        adv.append("\n"); // blank line inside sequence
        adv.append("A\n"); // ragged width 1
        adv.append("CCGGTTAA\r\n"); // ragged width 8
        adv.append(">hundred_thousand_Ns_scaffold\n");
        adv.append("ACGTACGT");
        adv.append("N".repeat(100_000));
        adv.append("tgcaTGCA\n");
        adv.append(">trailing_empty_record_no_newline");
        passed += checkRoundTrip("Adversarial Chimeric/Ragged/IUPAC/100k-N/EmptyRecs",
                adv.toString().getBytes(StandardCharsets.UTF_8), tmpDir);

        // 4. Multi-strain homologous + reverse-complement + SNP + indel + query test
        passed += verifyPangenomeAndRandomAccess(gcfLf, tmpDir);

        // 5. Giant multi-block single-chromosome slice query test
        passed += verifyMultiBlockChromosomeQuery(bueb18940, tmpDir);

        // 6. Multi-record input whose blocks start inside long records
        passed += verifyMultiRecordMidBlockQueries(bueb18940, tmpDir);

        // 7. Bounded file path: streamed archive assembly + parallel block window
        passed += verifyBoundedFilePaths(tmpDir);

        // 8. Opt-in single-genome profile: serial and segmented compact archives
        passed += verifySingleGenomeModes(tmpDir);

        System.out.println("\nALL " + passed + " VERIFICATION SUITES PASSED WITH 100% SHA-256 IDENTITY!");
    }

    private static int checkRoundTrip(String label, byte[] original, Path tmpDir) throws Exception {
        String expectedSha = sha256(original);

        // A. Standard mode
        byte[] compStd = HyperGraSS.compressBytes(original, HyperGraSS.CompressOptions.defaults());
        byte[] decStd = HyperGraSS.decompressBytes(compStd);
        if (!Arrays.equals(original, decStd) || !expectedSha.equals(sha256(decStd))) {
            throw new AssertionError("SHA-256 mismatch in Standard Mode on " + label);
        }

        // B. Streaming --lean file mode
        Path inFile = tmpDir.resolve("in.dat");
        Path hgsLean = tmpDir.resolve("out_lean.hgs");
        Path decFile = tmpDir.resolve("dec_lean.dat");
        Files.write(inFile, original);
        HyperGraSS.compressFile(inFile, hgsLean, HyperGraSS.CompressOptions.lean());
        HyperGraSS.decompressFile(hgsLean, decFile);
        byte[] decLean = Files.readAllBytes(decFile);
        if (!Arrays.equals(original, decLean) || !expectedSha.equals(sha256(decLean))) {
            throw new AssertionError("SHA-256 mismatch in --lean Streaming Mode on " + label);
        }

        double cr = compStd.length > 0 ? (double) original.length / compStd.length : 0.0;
        System.out.printf("  [PASS] %-48s | Raw: %,8d B -> HGS: %,7d B (CR: %6.2fx) | Lean: %,7d B | SHA256=%s%n",
                label, original.length, compStd.length, cr, Files.size(hgsLean), expectedSha.substring(0, 12));
        return 1;
    }

    private static int verifyPangenomeAndRandomAccess(byte[] gcfLf, Path tmpDir) throws Exception {
        // Extract pure sequence of GCF_001884535.1
        String txt = new String(gcfLf, StandardCharsets.US_ASCII);
        int nl = txt.indexOf('\n');
        String baseSeq = txt.substring(nl + 1).replace("\n", "");

        Random rng = new Random(42L);
        int numStrains = 60;
        String[] strainSeqs = new String[numStrains];
        StringBuilder multiFasta = new StringBuilder(numStrains * (baseSeq.length() + 120));

        for (int s = 0; s < numStrains; s++) {
            char[] arr = baseSeq.toCharArray();
            // Inject 25 SNPs (75% transitions, 25% transversions)
            for (int m = 0; m < 25; m++) {
                int pos = rng.nextInt(arr.length);
                arr[pos] = mutateBase(arr[pos], rng.nextDouble() < 0.75, rng);
            }
            String seq = new String(arr);
            // Circularly permute 1/3 of strains
            if (s % 3 == 1) {
                int cut = 5000 + rng.nextInt(30000);
                seq = seq.substring(cut) + seq.substring(0, cut);
            }
            // Reverse-complement 50% of strains (- strand!)
            if ((s & 1) == 1) {
                seq = revComp(seq);
            }
            // Add soft-masked lowercase run and N-gap in every 5th strain
            if (s % 5 == 0) {
                seq = seq.substring(0, 1000)
                        + seq.substring(1000, 1120).toLowerCase()
                        + "NNNNNNNNNN"
                        + seq.substring(1130);
            }
            strainSeqs[s] = seq;
            multiFasta.append(">GCF_001884535.1_strain_").append(String.format("%04d", s))
                    .append(" Gordonia phage GAL1 isolate\n");
            for (int i = 0; i < seq.length(); i += 70) {
                multiFasta.append(seq, i, Math.min(seq.length(), i + 70)).append('\n');
            }
        }

        byte[] raw = multiFasta.toString().getBytes(StandardCharsets.US_ASCII);
        Path inPath = tmpDir.resolve("pangenome60.fna");
        Path hgsPath = tmpDir.resolve("pangenome60.hgs");
        Path outPath = tmpDir.resolve("pangenome60.dec.fna");
        Files.write(inPath, raw);

        HyperGraSS.compressFile(inPath, hgsPath, HyperGraSS.CompressOptions.defaults());
        HyperGraSS.decompressFile(hgsPath, outPath);
        byte[] dec = Files.readAllBytes(outPath);
        if (!Arrays.equals(raw, dec)) {
            throw new AssertionError("SHA-256 mismatch on 60-strain mixed-strand pangenome!");
        }

        // Verify 50 random O(1) slice queries + full-record queries
        for (int q = 0; q < 50; q++) {
            int rec = rng.nextInt(numStrains);
            int start = rng.nextInt(strainSeqs[rec].length() - 200);
            int end = start + 1 + rng.nextInt(500);
            String expectedSub = strainSeqs[rec].substring(start, end);
            String gotSub = new String(HyperGraSS.query(hgsPath, rec, start, end), StandardCharsets.US_ASCII);
            if (!expectedSub.equals(gotSub)) {
                throw new AssertionError("Random-access slice mismatch on strain " + rec
                        + " [" + start + ".." + end + "): expected=" + expectedSub + ", got=" + gotSub);
            }
        }

        long compLen = Files.size(hgsPath);
        System.out.printf("  [PASS] %-48s | Raw: %,8d B -> HGS: %,7d B (CR: %6.2fx) | 50 O(1) Queries OK%n",
                "60-Strain (+/-) Permuted SNP Pangenome", raw.length, compLen, (double) raw.length / compLen);
        return 1;
    }

    private static int verifyMultiBlockChromosomeQuery(byte[] seedSeq, Path tmpDir) throws Exception {
        // Build a 1.5 MB single chromosome spanning multiple blocks
        String unit = new String(seedSeq, StandardCharsets.US_ASCII);
        StringBuilder chr = new StringBuilder(1_600_000);
        Random rng = new Random(99L);
        for (int c = 0; c < 80; c++) {
            char[] arr = unit.toCharArray();
            for (int m = 0; m < 40; m++) {
                int p = rng.nextInt(arr.length);
                arr[p] = mutateBase(arr[p], true, rng);
            }
            chr.append(arr);
        }
        String fullChr = chr.toString();
        StringBuilder fasta = new StringBuilder(fullChr.length() + 40_000);
        fasta.append(">chr1_synthetic_multiblock\n");
        for (int i = 0; i < fullChr.length(); i += 80) {
            fasta.append(fullChr, i, Math.min(fullChr.length(), i + 80)).append('\n');
        }
        byte[] raw = fasta.toString().getBytes(StandardCharsets.US_ASCII);
        Path inPath = tmpDir.resolve("chr1.fna");
        Path hgsPath = tmpDir.resolve("chr1.hgs");
        Files.write(inPath, raw);

        HyperGraSS.compressFile(inPath, hgsPath, HyperGraSS.CompressOptions.defaults());
        byte[] dec = HyperGraSS.decompressBytes(Files.readAllBytes(hgsPath));
        if (!Arrays.equals(raw, dec)) {
            throw new AssertionError("Multi-block chromosome round-trip mismatch!");
        }

        // Query across a block boundary (around 500,000 .. 505,000)
        int qStart = 498_500;
        int qEnd = 504_200;
        String got = new String(HyperGraSS.query(hgsPath, 0, qStart, qEnd), StandardCharsets.US_ASCII);
        String exp = fullChr.substring(qStart, qEnd);
        if (!exp.equals(got)) {
            throw new AssertionError("Cross-block chromosome query mismatch!");
        }
        System.out.printf("  [PASS] %-48s | Raw: %,8d B -> HGS: %,7d B | Cross-Block Query OK%n",
                "1.5 MB Multi-Block Single Chromosome", raw.length, Files.size(hgsPath));
        return 1;
    }

    private static int verifyMultiRecordMidBlockQueries(byte[] seedSeq, Path tmpDir) throws Exception {
        // Six near-identical records with a 16 KB block target. Each record is
        // longer than twice the target, so blocks start inside records and the
        // last record's header lands inside a mid-record block. The footer used
        // to count only blocks that start at a record boundary, which reported
        // 4 records for these 6 and rejected queries for the last two.
        int numRecords = 6;
        int recordLen = 40_000;
        int blockBytes = 16_384;
        String unit = new String(seedSeq, StandardCharsets.US_ASCII);
        Random rng = new Random(7L);
        String[] seqs = new String[numRecords];
        StringBuilder fasta = new StringBuilder(numRecords * recordLen + 16_000);
        for (int r = 0; r < numRecords; r++) {
            StringBuilder seq = new StringBuilder(recordLen + unit.length());
            while (seq.length() < recordLen) {
                seq.append(unit);
            }
            char[] arr = seq.substring(0, recordLen).toCharArray();
            for (int m = 0; m < 200; m++) {
                int p = rng.nextInt(arr.length);
                arr[p] = mutateBase(arr[p], rng.nextDouble() < 0.75, rng);
            }
            String record = new String(arr);
            seqs[r] = record;
            fasta.append(">synthetic_record_").append(r).append('\n');
            for (int i = 0; i < record.length(); i += 70) {
                fasta.append(record, i, Math.min(record.length(), i + 70)).append('\n');
            }
        }
        byte[] raw = fasta.toString().getBytes(StandardCharsets.US_ASCII);
        Path inPath = tmpDir.resolve("multi-record-mid-block.fna");
        Path hgsPath = tmpDir.resolve("multi-record-mid-block.hgs");
        Files.write(inPath, raw);

        HyperGraSS.CompressOptions opts = HyperGraSS.CompressOptions.defaults();
        opts.blockSizeBytes = blockBytes;
        HyperGraSS.compressFile(inPath, hgsPath, opts);
        if (!Arrays.equals(raw, HyperGraSS.decompressBytes(Files.readAllBytes(hgsPath)))) {
            throw new AssertionError("Multi-record mid-block round-trip mismatch!");
        }
        int counted = HyperGraSS.inspect(hgsPath).totalRecords;
        if (counted != numRecords) {
            throw new AssertionError("Footer reports " + counted + " records, expected " + numRecords);
        }

        for (int r = 0; r < numRecords; r++) {
            int start = (r == numRecords - 1) ? seqs[r].length() - 700 : 20_000;
            int end = start + 512;
            String expected = seqs[r].substring(start, end);
            String got = new String(HyperGraSS.query(hgsPath, r, start, end), StandardCharsets.US_ASCII);
            if (!expected.equals(got)) {
                throw new AssertionError("Mid-block record query mismatch on record " + r);
            }
        }
        System.out.printf("  [PASS] %-48s | Raw: %,8d B -> HGS: %,7d B | %d Records Queried%n",
                "Records with Mid-Record Block Starts", raw.length, Files.size(hgsPath), numRecords);
        return 1;
    }

    /** Bounded file path: streamed archive assembly and the parallel block window. */
    private static int verifyBoundedFilePaths(Path tmpDir) throws Exception {
        Random rng = new Random(20_261_103L);
        int numRecords = 6;
        int recordLen = 120_000;
        char[] template = new char[recordLen];
        for (int i = 0; i < recordLen; i++) {
            template[i] = "ACGT".charAt(rng.nextInt(4));
        }
        StringBuilder fasta = new StringBuilder();
        for (int r = 0; r < numRecords; r++) {
            char[] arr = template.clone();
            for (int m = 0; m < 400; m++) {
                int p = rng.nextInt(arr.length);
                arr[p] = mutateBase(arr[p], rng.nextDouble() < 0.75, rng);
            }
            fasta.append(">bounded_record_").append(r).append('\n');
            for (int i = 0; i < arr.length; i += 70) {
                fasta.append(arr, i, Math.min(70, arr.length - i)).append('\n');
            }
        }
        byte[] raw = fasta.toString().getBytes(StandardCharsets.US_ASCII);
        Path inPath = tmpDir.resolve("bounded.fna");
        Files.write(inPath, raw);

        HyperGraSS.CompressOptions seq = HyperGraSS.CompressOptions.defaults();
        seq.bounded = true;
        seq.threads = 1;
        seq.blockSizeBytes = 64 * 1024;
        seq.maxAnchorBases = 32 * 1024;
        Path seqPath = tmpDir.resolve("bounded-seq.hgs");
        HyperGraSS.compressFile(inPath, seqPath, seq);
        byte[] seqArchive = Files.readAllBytes(seqPath);
        if ((seqArchive[4] & 8) == 0) {
            throw new AssertionError("Single-threaded bounded archive must carry FLAG_SHARED_GAP");
        }
        if (!Arrays.equals(raw, HyperGraSS.decompressBytes(seqArchive))) {
            throw new AssertionError("Bounded streamed archive round-trip mismatch");
        }
        Path seqOut = tmpDir.resolve("bounded-seq.out");
        HyperGraSS.decompressFile(seqPath, seqOut);
        if (!Arrays.equals(raw, Files.readAllBytes(seqOut))) {
            throw new AssertionError("Bounded streamed file round-trip mismatch");
        }
        Path seqPath2 = tmpDir.resolve("bounded-seq2.hgs");
        HyperGraSS.compressFile(inPath, seqPath2, seq);
        if (!Arrays.equals(seqArchive, Files.readAllBytes(seqPath2))) {
            throw new AssertionError("Bounded streamed writer is not deterministic");
        }

        HyperGraSS.CompressOptions stride = HyperGraSS.CompressOptions.defaults();
        stride.bounded = true;
        stride.threads = 1;
        stride.blockSizeBytes = 64 * 1024;
        stride.maxAnchorBases = 32 * 1024;
        stride.seedStep = 8;
        Path stridePath = tmpDir.resolve("bounded-stride8.hgs");
        HyperGraSS.compressFile(inPath, stridePath, stride);
        if (!Arrays.equals(raw, HyperGraSS.decompressBytes(
                Files.readAllBytes(stridePath)))) {
            throw new AssertionError("Stride-8 bounded archive round-trip mismatch");
        }

        HyperGraSS.CompressOptions par = HyperGraSS.CompressOptions.defaults();
        par.bounded = true;
        par.threads = 3;
        par.blockSizeBytes = 64 * 1024;
        par.maxAnchorBases = 32 * 1024;
        Path parPath = tmpDir.resolve("bounded-par.hgs");
        HyperGraSS.compressFile(inPath, parPath, par);
        byte[] parArchive = Files.readAllBytes(parPath);
        if ((parArchive[4] & 8) != 0) {
            throw new AssertionError("Parallel bounded archive must clear FLAG_SHARED_GAP");
        }
        if (!Arrays.equals(raw, HyperGraSS.decompressBytes(parArchive))) {
            throw new AssertionError("Parallel bounded archive round-trip mismatch");
        }
        Path parOut = tmpDir.resolve("bounded-par.out");
        HyperGraSS.decompressFile(parPath, parOut);
        if (!Arrays.equals(raw, Files.readAllBytes(parOut))) {
            throw new AssertionError("Parallel bounded file round-trip mismatch");
        }

        HyperGraSS.CompressOptions cm = HyperGraSS.CompressOptions.defaults();
        cm.bounded = true;
        cm.threads = 1;
        cm.blockSizeBytes = 64 * 1024;
        cm.maxAnchorBases = 256 * 1024;
        cm.anchorCmPages = true;
        Path cmPath = tmpDir.resolve("bounded-cm-pages.hgs");
        HyperGraSS.compressFile(inPath, cmPath, cm);
        byte[] cmArchive = Files.readAllBytes(cmPath);
        if (!Arrays.equals(raw, HyperGraSS.decompressBytes(cmArchive))) {
            throw new AssertionError("Paged context-mixing anchor round-trip mismatch");
        }
        HyperGraSS.ArchiveMetadata cmMeta = HyperGraSS.inspect(cmPath);
        if (cmMeta.numAnchorPages < 2) {
            throw new AssertionError("Paged context-mixing anchor did not split into pages");
        }
        Path cmOut = tmpDir.resolve("bounded-cm-pages.out");
        HyperGraSS.decompressFile(cmPath, cmOut);
        if (!Arrays.equals(raw, Files.readAllBytes(cmOut))) {
            throw new AssertionError("Paged context-mixing anchor file round-trip mismatch");
        }

        int seqBlocks = HyperGraSS.inspect(seqPath).numBlocks;
        int parBlocks = HyperGraSS.inspect(parPath).numBlocks;
        if (seqBlocks != parBlocks || seqBlocks < 2) {
            throw new AssertionError(
                    "Bounded block count mismatch: " + seqBlocks + " vs " + parBlocks);
        }
        System.out.printf("  [PASS] %-48s | Raw: %,8d B -> HGS: %,7d B | %d blocks, 1 vs 3 threads%n",
                "Bounded streaming and parallel block paths", raw.length, seqArchive.length, seqBlocks);
        return 1;
    }

    /**
     * Exercises the opt-in single-genome container: the tuned flag, the optional
     * segmented DNA channel, random-access queries, byte accounting, and the fact
     * that the default compact format is untouched.
     */
    private static int verifySingleGenomeModes(Path tmpDir) throws Exception {
        Random rng = new Random(20_261_104L);
        int len = 600_000;
        char[] fwd = new char[len];
        for (int i = 0; i < len; i++) {
            fwd[i] = "ACGT".charAt(rng.nextInt(4));
        }
        // A duplicated reverse-complement block gives the repeat heads real work.
        char[] rc = revComp(new String(fwd, 0, 40_000)).toCharArray();
        StringBuilder fasta = new StringBuilder(len + 60_000);
        fasta.append(">chr_single description\n");
        for (int i = 0; i < len; i += 80) {
            fasta.append(fwd, i, Math.min(80, len - i)).append('\n');
        }
        fasta.append(">rc_block\n");
        for (int i = 0; i < rc.length; i += 80) {
            fasta.append(rc, i, Math.min(80, rc.length - i)).append('\n');
        }
        byte[] raw = fasta.toString().getBytes(StandardCharsets.US_ASCII);
        String expectedSlice = new String(fwd, 1000, 100);

        for (int segments : new int[] {1, 2, 5}) {
            HyperGraSS.CompressOptions opts = HyperGraSS.CompressOptions.defaults();
            opts.singleMode = true;
            opts.threads = segments;
            byte[] archive = HyperGraSS.compressBytes(raw, opts);
            if ((archive[4] & 16) == 0) {
                throw new AssertionError("Single-genome archive lacks the tuned flag");
            }
            boolean segmented = (archive[4] & 32) != 0;
            if (segmented != (segments > 1)) {
                throw new AssertionError("Segmented flag does not match thread count");
            }
            if (!Arrays.equals(raw, HyperGraSS.decompressBytes(archive))) {
                throw new AssertionError("Single-genome round-trip mismatch at "
                        + segments + " segment(s)");
            }
            ArchiveByteAccounting.Report account = ArchiveByteAccounting.report(archive);
            if (account.totalBytes() != archive.length) {
                throw new AssertionError("Single-genome byte accounting mismatch");
            }
            Path hgs = tmpDir.resolve("single-" + segments + ".hgs");
            Path out = tmpDir.resolve("single-" + segments + ".out");
            Files.write(hgs, archive);
            HyperGraSS.decompressFile(hgs, out);
            if (!Arrays.equals(raw, Files.readAllBytes(out))) {
                throw new AssertionError("Single-genome file round-trip mismatch at "
                        + segments + " segment(s)");
            }
            byte[] slice = HyperGraSS.query(hgs, 0, 1000, 1100);
            if (!expectedSlice.equals(new String(slice, StandardCharsets.US_ASCII))) {
                throw new AssertionError("Single-genome query mismatch at "
                        + segments + " segment(s)");
            }
            HyperGraSS.ArchiveMetadata meta = HyperGraSS.inspect(hgs);
            if (!meta.compactSingleBlock || meta.totalRawBytes != raw.length) {
                throw new AssertionError("Single-genome inspection mismatch");
            }
        }

        // The published compact path must remain byte-identical and decodable.
        HyperGraSS.CompressOptions published = HyperGraSS.CompressOptions.defaults();
        byte[] oldArchive = HyperGraSS.compressBytes(raw, published);
        if ((oldArchive[4] & 16) != 0) {
            throw new AssertionError("Default compact archive must not set the tuned flag");
        }
        if (!Arrays.equals(raw, HyperGraSS.decompressBytes(oldArchive))) {
            throw new AssertionError("Default compact archive no longer round-trips");
        }

        System.out.printf("  [PASS] %-48s | Raw: %,8d B -> tuned: %,7d B | 1/2/5 segments%n",
                "Single-genome tuned and segmented archives", raw.length, oldArchive.length);
        return 1;
    }

    static char mutateBase(char c, boolean transition, Random rng) {
        return switch (c) {
            case 'A' -> transition ? 'G' : (rng.nextBoolean() ? 'C' : 'T');
            case 'G' -> transition ? 'A' : (rng.nextBoolean() ? 'C' : 'T');
            case 'C' -> transition ? 'T' : (rng.nextBoolean() ? 'A' : 'G');
            case 'T' -> transition ? 'C' : (rng.nextBoolean() ? 'A' : 'G');
            default -> 'A';
        };
    }

    static String revComp(String s) {
        char[] out = new char[s.length()];
        for (int i = 0, j = s.length() - 1; i < s.length(); i++, j--) {
            out[i] = switch (s.charAt(j)) {
                case 'A' -> 'T';
                case 'T' -> 'A';
                case 'C' -> 'G';
                case 'G' -> 'C';
                case 'a' -> 't';
                case 't' -> 'a';
                case 'c' -> 'g';
                case 'g' -> 'c';
                default -> s.charAt(j);
            };
        }
        return new String(out);
    }

    static String sha256(byte[] data) throws Exception {
        MessageDigest md = MessageDigest.getInstance("SHA-256");
        return HexFormat.of().formatHex(md.digest(data));
    }
}
