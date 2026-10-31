package hypergrass;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.HexFormat;
import java.util.Locale;
import java.util.Random;

/** Cross-version archive-byte comparison using deterministic synthetic pangenomes. */
public final class VerifyDonorArchives {

    private VerifyDonorArchives() {}

    public static void main(String[] args) throws Exception {
        if (args.length < 2 || args.length > 4
                || !("baseline".equals(args[0]) || "verify".equals(args[0]))) {
            throw new IllegalArgumentException(
                    "Usage: VerifyDonorArchives <baseline|verify> <archive-directory> [--bench] [--optimize]");
        }
        boolean baseline = "baseline".equals(args[0]);
        boolean benchmark = false, optimize = false;
        for (int i = 2; i < args.length; i++) {
            if ("--bench".equals(args[i]) && !benchmark) benchmark = true;
            else if ("--optimize".equals(args[i]) && !optimize) optimize = true;
            else throw new IllegalArgumentException("Unknown or repeated option: " + args[i]);
        }
        Path dir = Path.of(args[1]);
        if (baseline) Files.createDirectories(dir);
        System.out.println("Synthetic archive compatibility check (not a real-genome benchmark).");
        if (benchmark) System.out.println("Serial compression: 4 warmups, median of 7 measured runs.");

        for (int kind = 0; kind < 4; kind++) {
            TestCase data = generateCase(kind);
            HyperGraSS.CompressOptions opts = HyperGraSS.CompressOptions.defaults();
            opts.parallel = false;
            if (optimize) opts.optimize = true;
            long start = System.nanoTime();
            byte[] archive = HyperGraSS.compressBytes(data.raw, opts);
            double ms = (System.nanoTime() - start) / 1e6;
            if (benchmark) {
                long[] times = new long[7];
                for (int run = 0; run < 11; run++) {
                    start = System.nanoTime();
                    byte[] repeated = HyperGraSS.compressBytes(data.raw, opts);
                    long elapsed = System.nanoTime() - start;
                    if (!Arrays.equals(repeated, archive)) {
                        throw new AssertionError("Nondeterministic archive on case " + kind);
                    }
                    if (run >= 4) times[run - 4] = elapsed;
                }
                Arrays.sort(times);
                ms = times[3] / 1e6;
            }
            if (!Arrays.equals(data.raw, HyperGraSS.decompressBytes(archive))) {
                throw new AssertionError("Round-trip mismatch on case " + kind);
            }
            Path path = dir.resolve("case" + kind + ".hgs");
            if (baseline) {
                Files.write(path, archive);
            } else if (!Arrays.equals(Files.readAllBytes(path), archive)) {
                throw new AssertionError("Archive bytes changed on case " + kind);
            }
            Random queries = new Random(9000 + kind);
            for (int q = 0; q < 15; q++) {
                int record = queries.nextInt(data.sequences.length);
                String seq = data.sequences[record];
                int from = queries.nextInt(seq.length() - 500);
                int to = from + 500;
                String got = new String(HyperGraSS.query(path, record, from, to), StandardCharsets.US_ASCII);
                if (!got.equals(seq.substring(from, to))) {
                    throw new AssertionError("Query mismatch on case " + kind + ", record " + record);
                }
            }
            String sha = HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(archive));
            System.out.printf(Locale.ROOT,
                    "PASS %s case%d raw=%d archive=%d %s=%.3f ms SHA256=%s%n",
                    args[0], kind, data.raw.length, archive.length, benchmark ? "median" : "cold", ms, sha);
        }
    }

    private static TestCase generateCase(int kind) {
        Random rng = new Random(4242 + kind);
        int length = 12_000, records = 128;
        char[] ancestor = new char[length];
        for (int p = 0; p < length; p++) ancestor[p] = "ACGT".charAt(rng.nextInt(4));
        StringBuilder text = new StringBuilder();
        String[] sequences = new String[records];
        for (int rec = 0; rec < records; rec++) {
            char[] seq = ancestor.clone();
            if (rec > 0) {
                int shared = kind == 0 ? 0 : kind == 1 ? 30 : 90;
                for (int p = 73; p < 73 + shared * 83; p += 83) seq[p] = mutate(seq[p]);
                if (kind == 2) {
                    for (int p = 131 + (rec % 8) * 19; p < length - 100; p += 1103) {
                        seq[p] = mutate(seq[p]);
                    }
                }
                int privateCount = kind == 3 ? 2 : 25;
                for (int m = 0; m < privateCount; m++) {
                    int p = rng.nextInt(length);
                    seq[p] = mutate(seq[p]);
                }
            }
            String sequence = new String(seq);
            if (rec % 3 == 1) {
                int cut = rng.nextInt(length);
                sequence = sequence.substring(cut) + sequence.substring(0, cut);
            }
            if ((rec & 1) != 0) sequence = reverseComplement(sequence);
            sequences[rec] = sequence;
            text.append(">strain_").append(rec).append('\n');
            for (int p = 0; p < length; p += 71) {
                text.append(sequence, p, Math.min(length, p + 71)).append('\n');
            }
        }
        return new TestCase(text.toString().getBytes(StandardCharsets.US_ASCII), sequences);
    }

    private static char mutate(char c) {
        return "ACGT".charAt(("ACGT".indexOf(c) + 2) & 3);
    }

    private static String reverseComplement(String s) {
        char[] out = new char[s.length()];
        for (int i = 0; i < out.length; i++) {
            out[i] = "TGCA".charAt("ACGT".indexOf(s.charAt(out.length - 1 - i)));
        }
        return new String(out);
    }

    private record TestCase(byte[] raw, String[] sequences) {}
}
