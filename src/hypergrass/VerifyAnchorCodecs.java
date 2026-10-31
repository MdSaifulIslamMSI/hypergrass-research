package hypergrass;

import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

/** Synthetic anchor-codec selection, size non-regression, and cross-version decoding checks. */
public final class VerifyAnchorCodecs {

    private VerifyAnchorCodecs() {}

    private record TestCase(String name, int length, int copies, boolean biased, boolean lean) {}

    private static final TestCase[] CASES = {
        new TestCase("uniform-high-reuse", 8_193, 16, false, false),
        new TestCase("biased-high-reuse", 8_193, 16, true, false),
        new TestCase("uniform-low-reuse", 8_193, 6, false, false),
        new TestCase("biased-low-reuse", 8_193, 6, true, false),
        new TestCase("short-unpadded", 509, 16, false, false),
        new TestCase("page-boundary", 65_536, 16, false, false),
        new TestCase("multipage-high-reuse", 65_537, 16, false, false),
        new TestCase("multipage-low-reuse", 65_537, 6, false, false),
        new TestCase("lean-no-anchor", 8_193, 6, false, true),
    };

    public static void main(String[] args) throws Exception {
        if (args.length != 0 && (args.length != 2
                || !("baseline".equals(args[0]) || "verify".equals(args[0]) || "decode".equals(args[0])))) {
            throw new IllegalArgumentException("Usage: VerifyAnchorCodecs [baseline|verify|decode <directory>]");
        }
        String mode = args.length == 0 ? "test" : args[0];
        boolean temporary = args.length == 0;
        Path directory = temporary ? Files.createTempDirectory("hypergrass_anchor_") : Path.of(args[1]);
        if (temporary || "baseline".equals(mode)) Files.createDirectories(directory);
        int packedChoices = 0, qcrmChoices = 0;
        try {
            for (TestCase test : CASES) {
                String sequence = sequence(test);
                byte[] raw = fasta(test, sequence);
                Path baseline = directory.resolve(test.name + ".hgs");
                Path updated = directory.resolve(test.name + ".updated.hgs");
                if ("decode".equals(mode)) {
                    verifyDecode(test, sequence, raw, updated);
                    System.out.println("PASS previous decoder " + test.name);
                    continue;
                }
                HyperGraSS.CompressOptions opts = test.lean
                        ? HyperGraSS.CompressOptions.lean() : HyperGraSS.CompressOptions.defaults();
                opts.blockSizeBytes = test.length + test.length / 71 + 256;
                opts.parallel = false;
                byte[] archive = HyperGraSS.compressBytes(raw, opts);
                if (!"baseline".equals(mode)) {
                    opts.parallel = true;
                    require(Arrays.equals(archive, HyperGraSS.compressBytes(raw, opts)),
                            test.name + ": parallel archive changed");
                }
                int choice = checkAnchor(test, archive, !"baseline".equals(mode));
                if (choice == 0) packedChoices++;
                if (choice == 1) qcrmChoices++;
                long oldSize = archive.length;
                if ("verify".equals(mode)) {
                    byte[] old = Files.readAllBytes(baseline);
                    oldSize = old.length;
                    verifyDecode(test, sequence, raw, baseline);
                    require(archive.length <= old.length, test.name + ": archive grew");
                    if (test.lean || test.length > CanonicalSyncmerEngine.ANCHOR_PAGE_SIZE) {
                        require(Arrays.equals(old, archive), test.name + ": unchanged path differs");
                    }
                }
                Path output = "baseline".equals(mode) ? baseline : updated;
                Files.write(output, archive);
                verifyDecode(test, sequence, raw, output);
                System.out.printf(Locale.ROOT, "PASS %s %s raw=%d old=%d archive=%d anchor=%s%n",
                        mode, test.name, raw.length, oldSize, archive.length,
                        choice < 0 ? "none" : choice == 0 ? "packed/Deflate" : "QCRM");
            }
            if (!"decode".equals(mode)) {
                require(packedChoices > 0 && qcrmChoices > 0, "Both existing anchor modes must be exercised");
            }
            System.out.println("ALL " + CASES.length + " ANCHOR CODEC SUITES PASSED");
        } finally {
            if (temporary) {
                for (TestCase test : CASES) Files.deleteIfExists(directory.resolve(test.name + ".updated.hgs"));
                Files.deleteIfExists(directory);
            }
        }
    }

    private static int checkAnchor(TestCase test, byte[] archive, boolean checkMinimum) {
        ByteIn in = new ByteIn(archive);
        in.readBytes(4);
        int flags = in.readByte();
        if (test.lean) {
            require((flags & 1) == 0, "Lean mode unexpectedly has an anchor");
            return -1;
        }
        require((flags & 1) != 0 && (flags & 4) == 0, test.name + ": expected multiblock anchor archive");
        int mode = in.readByte();
        int length = in.readVarInt();
        int pages = in.readVarInt();
        require(length > 0 && length <= test.length * test.copies,
                test.name + ": unexpected anchor length " + length);
        if (!test.biased && test.length >= 4_096) {
            require(length == test.length, test.name + ": expected complete first-record anchor");
        }
        boolean legacyPacked = (long) test.length * test.copies >= 12L * length;
        require(pages == (length + 65_535) / 65_536, test.name + ": unexpected page count");
        CanonicalSyncmerEngine.readGlobalDonors(in);
        int[] lengths = new int[pages], sizes = new int[pages];
        for (int p = 0; p < pages; p++) {
            lengths[p] = in.readVarInt();
            sizes[p] = in.readVarInt();
            require(lengths[p] == Math.min(65_536, length - p * 65_536),
                    test.name + ": unexpected page base length");
        }
        if (pages == 1 && checkMinimum) {
            byte[] payload = in.readBytes(sizes[0]);
            byte[] bases = mode == 0
                    ? VarIntCodec.unpack2Bit(VarIntCodec.readCompressedAux(new ByteIn(payload)), length)
                    : QuantumContextRepeatMixer.decode(new ByteIn(payload), length, false);
            byte[] packed = VarIntCodec.pack2Bit(bases, 0, length);
            ByteOut out = new ByteOut();
            VarIntCodec.writeCompressedAux(out, packed, packed.length);
            byte[] qcrm = QuantumContextRepeatMixer.encode(bases, length, false);
            int expected = out.size() == qcrm.length ? (legacyPacked ? 0 : 1)
                    : out.size() < qcrm.length ? 0 : 1;
            require(mode == expected, test.name + ": did not select the smaller existing codec");
            require(Arrays.equals(payload, mode == 0 ? out.toByteArray() : qcrm),
                    test.name + ": unexpected encoded anchor bytes");
        } else {
            require(mode == (legacyPacked ? 0 : 1), test.name + ": legacy selection changed");
        }
        return mode;
    }

    private static void verifyDecode(TestCase test, String sequence, byte[] raw, Path archive) throws Exception {
        require(Arrays.equals(raw, HyperGraSS.decompressBytes(Files.readAllBytes(archive))),
                test.name + ": byte-exact round-trip failed");
        for (int record : new int[] {0, test.copies / 2, test.copies - 1}) {
            String expected = (record & 1) == 0 ? sequence : sequence.toLowerCase(Locale.ROOT);
            for (int start : new int[] {0, test.length / 2, test.length - 13}) {
                int end = Math.min(test.length, start + 127);
                require(Arrays.equals(expected.substring(start, end).getBytes(StandardCharsets.US_ASCII),
                        HyperGraSS.query(archive, record, start, end)), test.name + ": slice mismatch");
            }
        }
    }

    private static String sequence(TestCase test) {
        Random random = new Random(8317);
        char[] bases = new char[test.length];
        for (int p = 0; p < bases.length; p++) {
            int value = test.biased && random.nextInt(10) < 7 ? 0 : random.nextInt(4);
            bases[p] = "ACGT".charAt(value);
        }
        return new String(bases);
    }

    private static byte[] fasta(TestCase test, String sequence) {
        StringBuilder text = new StringBuilder();
        for (int record = 0; record < test.copies; record++) {
            String eol = (record & 1) == 0 ? "\n" : "\r\n";
            String formatted = (record & 1) == 0 ? sequence : sequence.toLowerCase(Locale.ROOT);
            text.append('>').append(test.name).append('_').append(record).append(eol);
            int width = 71 + record % 3;
            for (int p = 0; p < formatted.length(); p += width) {
                text.append(formatted, p, Math.min(formatted.length(), p + width));
                if (p + width < formatted.length() || record + 1 < test.copies) text.append(eol);
            }
        }
        return text.toString().getBytes(StandardCharsets.US_ASCII);
    }

    private static void require(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
