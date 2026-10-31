package hypergrass;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

/** Saved-build size non-regression, exact reconstruction, queries and previous-decoder checks. */
public final class VerifyCostAwareEncoding {

    private VerifyCostAwareEncoding() {}

    private record TestCase(String name, int records, int length, int shared, int individual,
                            boolean groups, boolean rotate, boolean formatting, boolean lean,
                            boolean unrelated, int blockBytes) {}

    private static final TestCase[] CASES = {
        new TestCase("private", 64, 12_000, 0, 25, false, false, false, false, false, 524_288),
        new TestCase("shared", 64, 12_000, 48, 10, false, false, false, false, false, 524_288),
        new TestCase("subgroups", 80, 12_000, 40, 6, true, false, false, false, false, 524_288),
        new TestCase("dense", 48, 12_000, 120, 2, true, false, false, false, false, 524_288),
        new TestCase("small-collection", 8, 12_000, 50, 3, true, false, false, false, false, 524_288),
        new TestCase("rotations", 96, 12_000, 90, 10, true, true, false, false, false, 524_288),
        new TestCase("formatting", 40, 12_000, 40, 20, true, true, true, false, false, 524_288),
        new TestCase("no-variants", 24, 12_000, 0, 0, false, false, false, false, false, 524_288),
        new TestCase("unrelated", 8, 8_193, 0, 0, false, false, false, false, true, 524_288),
        new TestCase("lean", 16, 8_193, 40, 6, false, false, true, true, false, 65_536),
        new TestCase("large-collection", 192, 12_000, 120, 6, true, false, false, false, false, 524_288),
        new TestCase("many-blocks", 80, 8_193, 40, 6, true, false, false, false, false, 16_384),
    };

    public static void main(String[] args) throws Exception {
        if (args.length < 2 || args.length > 3
                || !("baseline".equals(args[0]) || "verify".equals(args[0]) || "decode".equals(args[0]))
                || (args.length == 3 && !"--bench".equals(args[2]))) {
            throw new IllegalArgumentException(
                    "Usage: VerifyCostAwareEncoding <baseline|verify|decode> <directory> [--bench]");
        }
        String mode = args[0];
        Path directory = Path.of(args[1]);
        boolean benchmark = args.length == 3;
        if ("baseline".equals(mode)) Files.createDirectories(directory);
        long saved = 0;
        int improved = 0;
        System.out.println("Synthetic regression corpus; not an independent real-data benchmark.");
        if (benchmark) System.out.println("Serial compression: 4 warmups, median of 7 measured runs.");
        for (int kind = 0; kind < CASES.length; kind++) {
            TestCase test = CASES[kind];
            String[] sequences = sequences(test, kind);
            byte[] raw = fasta(test, sequences);
            Path previous = directory.resolve(test.name + ".hgs");
            Path updated = directory.resolve(test.name + ".updated.hgs");
            if ("decode".equals(mode)) {
                verifyDecode(test, sequences, raw, updated);
                System.out.println("PASS previous decoder " + test.name);
                continue;
            }
            HyperGraSS.CompressOptions opts = test.lean
                    ? HyperGraSS.CompressOptions.lean() : HyperGraSS.CompressOptions.defaults();
            opts.blockSizeBytes = test.blockBytes;
            opts.parallel = false;
            if (!"baseline".equals(mode)) opts.optimize = true;
            long start = System.nanoTime();
            byte[] archive = HyperGraSS.compressBytes(raw, opts);
            double ms = (System.nanoTime() - start) / 1e6;
            if (benchmark) {
                long[] times = new long[7];
                for (int run = 0; run < 11; run++) {
                    start = System.nanoTime();
                    byte[] repeated = HyperGraSS.compressBytes(raw, opts);
                    long elapsed = System.nanoTime() - start;
                    require(Arrays.equals(archive, repeated), test.name + ": nondeterministic encoding");
                    if (run >= 4) times[run - 4] = elapsed;
                }
                Arrays.sort(times);
                ms = times[3] / 1e6;
            }
            ArchiveByteAccounting.Report report = ArchiveByteAccounting.report(archive);
            require(report.totalBytes() == archive.length, "Incomplete accounting");
            long oldSize = archive.length;
            if ("verify".equals(mode)) {
                byte[] old = Files.readAllBytes(previous);
                oldSize = old.length;
                require(archive.length <= old.length, test.name + ": archive grew");
                opts.optimize = false;
                require(Arrays.equals(old, HyperGraSS.compressBytes(raw, opts)),
                        test.name + ": default archive bytes changed");
                opts.optimize = true;
                verifyDecode(test, sequences, raw, previous);
                if (test.lean || "no-variants".equals(test.name) || test.unrelated) {
                    require(Arrays.equals(archive, old), test.name + ": unchanged path differs");
                }
                if (raw.length > 2 * 1024 * 1024 || "many-blocks".equals(test.name)) {
                    require(Arrays.equals(archive, old),
                            test.name + ": bounded cost selection changed an excluded archive");
                }
                saved += old.length - archive.length;
                if (archive.length < old.length) improved++;
                opts.parallel = true;
                require(Arrays.equals(archive, HyperGraSS.compressBytes(raw, opts)),
                        test.name + ": serial/parallel archives differ");
            }
            Path output = "baseline".equals(mode) ? previous : updated;
            Files.write(output, archive);
            verifyDecode(test, sequences, raw, output);
            System.out.printf(Locale.ROOT,
                    "PASS %s %s raw=%d old=%d archive=%d dictionary=%d %s=%.3f ms%n",
                    mode, test.name, raw.length, oldSize, archive.length,
                    report.categoryBytes().get("anchor.global_dictionary"), benchmark ? "median" : "cold", ms);
        }
        System.out.printf("ALL %d COST-AWARE ARCHIVE SUITES PASSED; improved=%d saved=%d bytes%n",
                CASES.length, improved, saved);
    }

    private static void verifyDecode(TestCase test, String[] sequences, byte[] raw, Path archive) throws Exception {
        require(Arrays.equals(raw, HyperGraSS.decompressBytes(Files.readAllBytes(archive))),
                test.name + ": byte-exact reconstruction failed");
        for (int record : new int[] {0, sequences.length / 2, sequences.length - 1}) {
            String seq = sequences[record];
            for (int start : new int[] {0, seq.length() / 2, seq.length() - 17}) {
                int end = Math.min(seq.length(), start + 257);
                require(Arrays.equals(seq.substring(start, end).getBytes(StandardCharsets.US_ASCII),
                        HyperGraSS.query(archive, record, start, end)), test.name + ": slice query failed");
            }
        }
    }

    private static String[] sequences(TestCase test, int kind) {
        Random random = new Random(0x71A9_2026L + kind);
        char[] initial = new char[test.length];
        for (int p = 0; p < initial.length; p++) initial[p] = "ACGT".charAt(random.nextInt(4));
        String[] sequences = new String[test.records];
        for (int rec = 0; rec < sequences.length; rec++) {
            char[] seq = initial.clone();
            if (test.unrelated) {
                for (int p = 0; p < seq.length; p++) seq[p] = "ACGT".charAt(random.nextInt(4));
            } else if (rec > 0) {
                for (int i = 0; i < test.shared; i++) mutate(seq, 43 + i * 71);
                if (test.groups) {
                    for (int p = 137 + (rec % 8) * 23; p < seq.length; p += 1091) mutate(seq, p);
                }
                for (int i = 0; i < test.individual; i++) mutate(seq, random.nextInt(seq.length));
            }
            String sequence = new String(seq);
            if (test.rotate) {
                int cut = random.nextInt(seq.length);
                sequence = sequence.substring(cut) + sequence.substring(0, cut);
                if ((rec & 1) != 0) sequence = reverseComplement(sequence);
            }
            if (test.formatting) {
                char[] formatted = sequence.toCharArray();
                formatted[101] = 'N';
                formatted[503] = 'R';
                sequence = new String(formatted);
                if ((rec & 1) != 0) sequence = sequence.toLowerCase(Locale.ROOT);
            }
            sequences[rec] = sequence;
        }
        return sequences;
    }

    private static byte[] fasta(TestCase test, String[] sequences) {
        StringBuilder text = new StringBuilder();
        for (int rec = 0; rec < sequences.length; rec++) {
            String eol = test.formatting && (rec & 1) != 0 ? "\r\n" : "\n";
            text.append('>').append(test.name).append('_').append(rec).append(eol);
            String sequence = sequences[rec];
            int width = test.formatting ? 65 + rec % 17 : 71;
            for (int p = 0; p < sequence.length(); p += width) {
                text.append(sequence, p, Math.min(sequence.length(), p + width));
                if (p + width < sequence.length() || rec + 1 < sequences.length) text.append(eol);
            }
        }
        return text.toString().getBytes(StandardCharsets.US_ASCII);
    }

    private static void mutate(char[] seq, int position) {
        seq[position] = "ACGT".charAt(("ACGT".indexOf(seq[position]) + 2) & 3);
    }

    private static String reverseComplement(String s) {
        char[] out = new char[s.length()];
        for (int p = 0; p < out.length; p++) out[p] = "TGCA".charAt("ACGT".indexOf(s.charAt(out.length - p - 1)));
        return new String(out);
    }

    private static void require(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
