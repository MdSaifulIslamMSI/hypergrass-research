package hypergrass;

import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.zip.CRC32;

/** Deterministic synthetic format/accounting checks; no external corpus or enabled assertions needed. */
public final class VerifyArchiveAccounting {

    private static int checks;
    private static int rejections;
    private static int suites;

    private VerifyArchiveAccounting() {}

    public static void main(String[] args) throws Exception {
        if (args.length != 0) throw new IllegalArgumentException("Usage: hypergrass.VerifyArchiveAccounting");
        Path directory = Files.createTempDirectory(Path.of("."), ".archive-accounting-verify-");
        try {
            exactEmpty();
            handBuiltModes(directory);
            realArchives(directory);
            malformedArchives();
            opaqueStreams();
            System.out.println("ALL " + suites + " ACCOUNTING SUITES PASSED (" + checks
                    + " checks, " + rejections + " malformed archives rejected)");
        } finally {
            try (var paths = Files.walk(directory)) {
                for (Path path : paths.sorted(Comparator.reverseOrder()).toList()) Files.delete(path);
            }
        }
    }

    private static void exactEmpty() {
        byte[] archive = compact(new byte[] {0, 0, 0, 0, 0, 0, 0});
        ArchiveByteAccounting.Report report = checkReport(archive);
        check(archive.length == 20, "Independent empty archive length");
        expect(report, "container.header", 13);
        expect(report, "metadata.framing", 4);
        expect(report, "metadata.line_widths", 1);
        expect(report, "metadata.line_endings", 1);
        expect(report, "dna.framing", 1);
        check(report.modeCounts().equals(Map.of(0, 1L, 1, 0L, 2, 0L, 3, 0L)), "Empty mode counts");
        for (var entry : report.categoryBytes().entrySet()) {
            if (!List.of("container.header", "metadata.framing", "metadata.line_widths",
                    "metadata.line_endings", "dna.framing").contains(entry.getKey())) {
                check(entry.getValue() == 0, "Absent empty category: " + entry.getKey());
            }
        }
        check(Arrays.equals(new byte[0], HyperGraSS.decompressBytes(archive)), "Hand-built empty decoder agreement");
        allTruncations(archive);
        pass("independent exact empty counts");
    }

    private static void handBuiltModes(Path directory) throws Exception {
        // These small streams are deliberately raw so the expected byte counts are independent
        // constants, not a second implementation of the accounting parser.
        for (int mode : new int[] {2, 3}) {
            for (boolean variants : new boolean[] {false, true}) {
                byte[] archive = matchArchive(mode, variants, 1, 0, 0);
                ArchiveByteAccounting.Report report = checkReport(archive);
                check(report.modeCounts().get(mode) == 1, "Hand-built match mode " + mode);
                expect(report, "container.header", 5);
                expect(report, "container.trailer", 12);
                expect(report, "anchor.header", 3);
                expect(report, "anchor.global_dictionary", 2);
                expect(report, "anchor.page_directory", 2);
                expect(report, "anchor.payload", 3);
                expect(report, "dna.framing", 1);
                expect(report, "primary.mode1_payload", 0);
                expect(report, "match.framing", 3);
                expect(report, "match.literal_lengths", 4);
                expect(report, "residual.framing", 3);
                expect(report, "residual.payload", 1);
                expect(report, "match.strands", 3);
                expect(report, "match.geometry", mode == 2 ? 4 : 5);
                expect(report, "variant.counts", variants ? 3 : 0);
                expect(report, "variant.gaps", variants ? 3 : 0);
                expect(report, "variant.differences", variants ? 3 : 0);
                expect(report, "donor.sweep_modes", variants && mode == 3 ? 3 : 0);
                expect(report, "donor.indices", variants && mode == 3 ? 1 : 0);
                expect(report, "donor.masks", variants && mode == 3 ? 1 : 0);
                byte[] expected = ascii(variants ? "TCCGT" : "TACGT");
                check(Arrays.equals(expected, HyperGraSS.decompressBytes(archive)), "Hand-built match decoder agreement");
                Path path = directory.resolve("fixture.hgs");
                Files.write(path, archive);
                check(report.equals(ArchiveByteAccounting.report(path)), "Path/byte[] API equivalence");
                allTruncations(archive);
                pass("Mode " + mode + " raw literals, variants=" + variants);
            }
        }
        for (int mode : new int[] {2, 3}) {
            for (int literals : new int[] {0, 64, 2048}) {
                int residualMode = literals >= 2048 ? 1 : 0;
                byte[] archive = matchArchive(mode, true, literals, residualMode, 0);
                ArchiveByteAccounting.Report report = checkReport(archive);
                if (literals == 0) {
                    expect(report, "residual.framing", 0);
                    expect(report, "residual.payload", 0);
                } else {
                    check(report.categoryBytes().get("residual.payload") > 0, "Nonempty residual bytes");
                }
                byte[] decoded = HyperGraSS.decompressBytes(archive);
                check(decoded.length == literals + 4, "Residual decoded length");
                check(decoded[literals] == 'C', "Variant follows residuals");
                allTruncations(archive);
                pass("Mode " + mode + " literal count=" + literals);
            }
        }
        for (int inheritance : new int[] {2, 3}) {
            byte[] archive = matchArchive(3, true, 0, 0, inheritance);
            ArchiveByteAccounting.Report report = checkReport(archive);
            expect(report, "donor.indices", 3);
            expect(report, "donor.masks", inheritance == 3 ? 3 : 1);
            expect(report, "variant.counts", 3);
            expect(report, "variant.gaps", 1);
            expect(report, "variant.differences", 1);
            check(report.categoryBytes().get("anchor.global_dictionary") > 2, "Global dictionary encoded bytes");
            check(Arrays.equals(ascii("CCGT"), HyperGraSS.decompressBytes(archive)), "Donor fixture decoder agreement");
            pass("Mode 3 donor inheritance=" + inheritance);
        }
    }

    private static void realArchives(Path directory) throws Exception {
        String[] texts = {"", "\n", "\r\n", "ACGT", ">empty\r\n>last", ">rna\nacguUUuNn\n",
                "raw-prefix\r\n>mixed\r\nacGTNnRYSWKMBDHV-*.\n\nUuTt\r\n>tail\nzz\r"};
        for (int i = 0; i < texts.length; i++) {
            byte[] raw = ascii(texts[i]);
            for (boolean lean : new boolean[] {false, true}) {
                HyperGraSS.CompressOptions options = lean
                        ? HyperGraSS.CompressOptions.lean() : HyperGraSS.CompressOptions.defaults();
                options.parallel = false;
                byte[] archive = HyperGraSS.compressBytes(raw, options);
                ArchiveByteAccounting.Report report = checkReport(archive);
                check(Arrays.equals(raw, HyperGraSS.decompressBytes(archive)), "Real edge round-trip");
                expect(report, "container.header", lean ? 5 : 13);
                check(report.blockCount() == 1, "Real edge block count");
                check(report.modeCounts().get(rawHasDna(raw) ? 1 : 0) == 1, "Real edge DNA mode");
                if (i == texts.length - 1) {
                    for (String channel : List.of("metadata.headers", "metadata.line_widths",
                            "metadata.line_endings", "metadata.case", "metadata.n_runs",
                            "metadata.exception_gaps", "metadata.exception_symbols")) {
                        check(report.categoryBytes().get(channel) > 0, "Mixed metadata channel: " + channel);
                    }
                }
                allTruncations(archive);
                pass("real edge " + i + (lean ? " lean" : " compact"));
            }
        }

        for (int kind = 0; kind < 3; kind++) {
            int length = kind == 2 ? 65_537 : 8_193;
            int copies = kind == 0 ? 16 : 6;
            byte[] raw = repeatedText(length, copies, kind == 1);
            HyperGraSS.CompressOptions options = HyperGraSS.CompressOptions.defaults();
            options.parallel = false;
            options.blockSizeBytes = length + length / 71 + 256;
            byte[] archive = HyperGraSS.compressBytes(raw, options);
            ArchiveByteAccounting.Report report = checkReport(archive);
            check((archive[4] & 5) == 1, "Real multiblock archive has anchor");
            check(report.blockCount() > 1, "Real multiblock block count");
            check(report.modeCounts().get(2) + report.modeCounts().get(3) > 0, "Real match mode used");
            check(Arrays.equals(raw, HyperGraSS.decompressBytes(archive)), "Real anchor archive round-trip");
            ByteIn anchor = new ByteIn(archive, 5, archive.length - 5);
            int anchorMode = anchor.readByte();
            anchor.readVarInt();
            int pages = anchor.readVarInt();
            CanonicalSyncmerEngine.readGlobalDonors(anchor);
            long directoryStart = anchor.position(), payloadBytes = 0;
            for (int page = 0; page < pages; page++) {
                anchor.readVarInt();
                payloadBytes += anchor.readVarInt();
            }
            expect(report, "anchor.page_directory", anchor.position() - directoryStart);
            expect(report, "anchor.payload", payloadBytes);
            check(kind != 2 || pages > 1, "Real multipage anchor");
            check(kind != 0 || anchorMode == 0, "Real packed anchor mode");
            check(kind != 1 || anchorMode == 1, "Real QCRM anchor mode");
            check(archive[footerOffset(archive)] == 2, "Real compressed footer directory");
            Path path = directory.resolve("real.hgs");
            Files.write(path, archive);
            check(report.equals(ArchiveByteAccounting.report(path)), "Real Path API");
            pass("real anchor kind=" + kind + " pages=" + pages + " codec=" + anchorMode);
        }

        Path input = directory.resolve("input.txt");
        Path archivePath = directory.resolve("streaming.hgs");
        Path decodedPath = directory.resolve("decoded.txt");
        for (byte[] raw : new byte[][] {new byte[0], repeatedText(8_193, 16, false)}) {
            Files.write(input, raw);
            HyperGraSS.compressStreamingLean(input, archivePath);
            byte[] archive = Files.readAllBytes(archivePath);
            ArchiveByteAccounting.Report report = checkReport(archive);
            check(report.equals(ArchiveByteAccounting.report(archivePath)), "Streaming Path API");
            expect(report, "anchor.payload", 0);
            expect(report, "container.header", 5);
            check(raw.length == 0 || report.blockCount() > 1, "Streaming multiblock count");
            HyperGraSS.decompressFile(archivePath, decodedPath);
            check(Arrays.equals(raw, Files.readAllBytes(decodedPath)), "Streaming round-trip");
            pass("real streaming lean raw=" + raw.length);
        }
    }

    private static void malformedArchives() {
        byte[] empty = compact(new byte[] {0, 0, 0, 0, 0, 0, 0});
        reject(changed(empty, 0, 'X'));
        reject(changed(empty, 3, 2));
        reject(changed(empty, 4, 12));
        reject(changed(empty, 4, 5));
        reject(changed(empty, 13, 32));
        reject(changed(empty, 14, 1));
        reject(changed(empty, 19, 4));
        reject(changed(empty, 19, 255));
        for (int mode : new int[] {3, 4, 255}) reject(changed(empty, 17, mode));
        reject(compact(new byte[] {0, 0, 0, 0, 1, 0})); // impossible raw auxiliary length
        reject(compact(new byte[] {0, 0, 0, 0, 2, 0, 1, 0}));
        reject(compact(new byte[] {0, 0, 0, 0, 2, 1, 0}));
        reject(compact(new byte[] {0, 0, 0, 0, 1, 127}));
        reject(compact(new byte[] {0, 0, 0, 0, 2, 1, 127}));
        for (byte[] badVarInt : new byte[][] {
                {(byte) 128}, {(byte) 128, (byte) 128, (byte) 128, (byte) 128, 8},
                {(byte) 128, (byte) 128, (byte) 128, (byte) 128, (byte) 128, 0},
                {(byte) 255, (byte) 255, (byte) 255, (byte) 255, 127}}) {
            ByteOut block = new ByteOut();
            block.writeByte(0);
            block.writeBytes(badVarInt);
            reject(compact(block.toByteArray()));
        }
        reject(Arrays.copyOf(empty, empty.length + 1));
        byte[] mode1 = compact(new byte[] {0, 1, 1, 0, 0, 0, 1, 127});
        reject(mode1);
        reject(compact(new byte[] {0, 1, 1, 0, 0, 0, 1, 0}));

        byte[] match = matchArchive(2, true, 1, 0, 0);
        reject(changed(match, 5, 2));
        reject(changed(match, 7, 2));
        reject(changed(match, 10, 3));
        reject(changed(match, 11, 127));
        reject(changed(match, 12, 3));
        reject(changed(match, 13, 2));
        reject(changed(match, match.length - 1, 0));
        for (long offset : new long[] {-1, 0, 4, match.length - 12L, Long.MAX_VALUE}) {
            byte[] bad = match.clone();
            putLongLE(bad, bad.length - 12, offset);
            reject(bad);
        }
        int footer = footerOffset(match);
        reject(changed(match, footer, 0));
        reject(changed(match, footer, 3));
        reject(changed(match, footer + 1, 127));
        // Locate the literal-mode byte in the independently hand-built raw fixture.
        ByteIn in = new ByteIn(match);
        in.setPosition(15);
        FastaStreamFactorizer.readAuxMetadata(in);
        int dnaOffset = in.position();
        reject(changed(match, dnaOffset, 4));
        in.readByte();
        in.readVarInt();
        in.readVarInt();
        in.readVarInt();
        VarIntCodec.readCompressedAux(in);
        int literalOffset = in.position();
        reject(changed(match, literalOffset, 2));
        reject(changed(match, literalOffset + 2, 0));
        reject(changed(match, literalOffset + 2, 127));
        reject(changed(match, literalOffset + 5, 2)); // wrong strand raw length
        byte[] noAnchor = match.clone();
        noAnchor[4] = 0;
        reject(noAnchor);
        pass("malformed magic, flags, modes, varints, lengths, pages and footer boundaries");
    }

    private static void opaqueStreams() {
        // Intentionally invalid Deflate and QCRM bytes must stay opaque to accounting.
        byte[] archive = compact(new byte[] {0, 1, 1, 1, 2, 100, 1, (byte) 255, 0, 1, 1, (byte) 255});
        ArchiveByteAccounting.Report report = checkReport(archive);
        expect(report, "metadata.line_widths", 4);
        expect(report, "primary.mode1_payload", 1);
        ByteOut footer = new ByteOut();
        footer.writeByte(2);
        footer.writeVarInt(Integer.MAX_VALUE);
        footer.writeVarInt(1);
        footer.writeByte(255);
        byte[] opaqueFooter = indexed(new byte[0], new byte[] {0, 0, 0, 0, 0, 0, 0}, 0,
                footer.toByteArray());
        expect(checkReport(opaqueFooter), "container.footer_directory", footer.size());
        // Preserve real framing while replacing every encoded QCRM anchor byte with junk.
        byte[] anchorPayload = {99, 98, 97};
        ByteOut anchor = new ByteOut();
        anchor.writeByte(1);
        anchor.writeVarInt(4);
        anchor.writeVarInt(1);
        anchor.writeVarInt(4);
        anchor.writeVarInt(0);
        anchor.writeVarInt(4);
        anchor.writeVarInt(anchorPayload.length);
        anchor.writeBytes(anchorPayload);
        byte[] opaqueAnchor = indexed(anchor.toByteArray(), new byte[] {0, 0, 0, 0, 0, 0, 0}, 0, null);
        expect(checkReport(opaqueAnchor), "anchor.payload", 3);
        pass("opaque compressed metadata, footer, primary and anchor payloads");
    }

    private static ArchiveByteAccounting.Report checkReport(byte[] archive) {
        byte[] before = archive.clone();
        ArchiveByteAccounting.Report report = ArchiveByteAccounting.report(archive);
        check(report.archiveBytes() == archive.length, "Archive length");
        check(report.totalBytes() == archive.length, "Exact encoded-byte sum");
        long sum = 0;
        for (long value : report.categoryBytes().values()) {
            check(value >= 0, "Nonnegative category");
            sum += value;
        }
        check(sum == archive.length, "Independent category sum");
        ArchiveByteAccounting.Report repeated = ArchiveByteAccounting.report(archive);
        check(report.equals(repeated), "Deterministic repeated report");
        check(report.toString().equals(repeated.toString()), "Deterministic iteration order");
        check(new ArrayList<>(report.modeCounts().keySet()).equals(List.of(0, 1, 2, 3)), "Ordered mode counts");
        check(Arrays.equals(before, archive), "Read-only archive bytes");
        try {
            report.categoryBytes().put("container.header", -1L);
            throw new AssertionError("Mutable category report");
        } catch (UnsupportedOperationException expected) {
            checks++;
        }
        try {
            report.modeCounts().put(0, -1L);
            throw new AssertionError("Mutable mode report");
        } catch (UnsupportedOperationException expected) {
            checks++;
        }
        return report;
    }

    private static void allTruncations(byte[] archive) {
        for (int length = 0; length < archive.length; length++) reject(Arrays.copyOf(archive, length));
    }

    private static void reject(byte[] archive) {
        try {
            ArchiveByteAccounting.report(archive);
        } catch (IllegalArgumentException expected) {
            check(expected.getMessage().contains("offset"), "Actionable malformed-archive offset");
            rejections++;
            return;
        }
        throw new AssertionError("Accepted malformed archive of " + archive.length + " bytes");
    }

    private static void expect(ArchiveByteAccounting.Report report, String category, long bytes) {
        check(report.categoryBytes().get(category) == bytes,
                category + ": expected " + bytes + ", got " + report.categoryBytes().get(category));
    }

    private static byte[] matchArchive(int mode, boolean variants, int literals, int residualMode, int inheritance) {
        ByteOut anchor = new ByteOut();
        anchor.writeByte(0);
        anchor.writeVarInt(4);
        anchor.writeVarInt(1);
        int[][] donors = inheritance == 0 ? new int[0][]
                : new int[][] {inheritance == 2 ? new int[] {1} : new int[] {1, 6}};
        CanonicalSyncmerEngine.writeGlobalDonors(anchor, 4, donors);
        anchor.writeVarInt(4);
        anchor.writeVarInt(3);
        anchor.writeBytes(new byte[] {1, 1, (byte) 228}); // packed ACGT

        byte[] literalBases = new byte[literals];
        Arrays.fill(literalBases, (byte) 3);
        byte[] raw = new byte[literals + 4];
        Arrays.fill(raw, 0, literals, (byte) 'T');
        System.arraycopy(ascii(variants ? "CCGT" : "ACGT"), 0, raw, literals, 4);
        var chunk = FastaStreamFactorizer.factorize(raw, 0, raw.length);
        ByteOut block = new ByteOut();
        FastaStreamFactorizer.writeAuxMetadata(block, chunk);
        block.writeByte(mode);
        block.writeVarInt(1);
        block.writeVarInt(literals);
        block.writeVarInt(variants ? 1 : 0);
        ByteOut lengths = new ByteOut();
        lengths.writeVarInt(literals);
        lengths.writeVarInt(0);
        VarIntCodec.writeCompressedAux(block, lengths);
        if (literals > 0) {
            block.writeByte(residualMode);
            if (residualMode == 0) {
                byte[] packed = VarIntCodec.pack2Bit(literalBases, 0, literals);
                VarIntCodec.writeCompressedAux(block, packed, packed.length);
            } else {
                byte[] encoded = QuantumContextRepeatMixer.encode(literalBases, literals, false);
                block.writeVarInt(encoded.length);
                block.writeBytes(encoded);
            }
        }
        block.writeBytes(new byte[] {1, 1, 0}); // strand bits
        block.writeBytes(mode == 2 ? new byte[] {1, 2, 0, 8} : new byte[] {1, 3, 4, 0, 0});
        if (variants) {
            if (mode == 3) {
                block.writeBytes(new byte[] {1, 1, (byte) (inheritance == 0 ? 1 : inheritance)});
                block.writeBytes(inheritance == 0 ? new byte[] {0} : new byte[] {1, 1, 0});
                block.writeBytes(inheritance == 3 ? new byte[] {1, 1, 1} : new byte[] {0});
            }
            block.writeBytes(new byte[] {1, 1, (byte) (inheritance == 0 ? 1 : 0)});
            block.writeBytes(inheritance == 0 ? new byte[] {1, 1, 0} : new byte[] {0});
            block.writeBytes(inheritance == 0 ? new byte[] {1, 1, 1} : new byte[] {0});
        }
        return indexed(anchor.toByteArray(), block.toByteArray(), raw.length, null, raw);
    }

    private static byte[] compact(byte[] block) {
        ByteOut out = new ByteOut();
        out.writeBytes(new byte[] {'H', 'G', 'S', 1, 4});
        out.writeIntLE(0);
        out.writeIntLE(0);
        out.writeBytes(block);
        return out.toByteArray();
    }

    private static byte[] indexed(byte[] anchor, byte[] block, int rawLength, byte[] footer) {
        return indexed(anchor, block, rawLength, footer, new byte[rawLength]);
    }

    private static byte[] indexed(byte[] anchor, byte[] block, int rawLength, byte[] footer, byte[] raw) {
        ByteOut out = new ByteOut();
        out.writeBytes(new byte[] {'H', 'G', 'S', 1, (byte) (anchor.length == 0 ? 0 : 1)});
        out.writeBytes(anchor);
        int blockOffset = out.size();
        out.writeBytes(block);
        int footerOffset = out.size();
        if (footer == null) {
            ByteOut body = new ByteOut();
            body.writeVarInt(1);
            body.writeVarLong(rawLength);
            body.writeVarInt(rawLength == 0 ? 0 : 1);
            CRC32 crc = new CRC32();
            crc.update(raw);
            body.writeIntLE((int) crc.getValue());
            body.writeVarLong(blockOffset);
            body.writeVarInt(block.length);
            body.writeVarInt(0);
            body.writeVarInt(1);
            body.writeVarLong(0);
            body.writeVarInt(rawLength);
            VarIntCodec.writeCompressedAux(out, body);
        } else {
            out.writeBytes(footer);
        }
        out.writeLongLE(footerOffset);
        out.writeIntLE(0x31464748);
        return out.toByteArray();
    }

    private static byte[] repeatedText(int length, int copies, boolean biased) {
        Random random = new Random(8317);
        char[] sequence = new char[length];
        for (int i = 0; i < length; i++) {
            int value = biased && random.nextInt(10) < 7 ? 0 : random.nextInt(4);
            sequence[i] = "ACGT".charAt(value);
        }
        String text = new String(sequence);
        StringBuilder raw = new StringBuilder();
        for (int record = 0; record < copies; record++) {
            String eol = (record & 1) == 0 ? "\n" : "\r\n";
            String formatted = (record & 1) == 0 ? text : text.toLowerCase(Locale.ROOT);
            raw.append(">record_").append(record).append(eol);
            for (int offset = 0; offset < length; offset += 71) {
                raw.append(formatted, offset, Math.min(length, offset + 71));
                if (offset + 71 < length || record + 1 < copies) raw.append(eol);
            }
        }
        return ascii(raw.toString());
    }

    private static boolean rawHasDna(byte[] raw) {
        return FastaStreamFactorizer.factorize(raw, 0, raw.length).dnaLen > 0;
    }

    private static int footerOffset(byte[] archive) {
        return (int) new ByteIn(archive, archive.length - 12, 12).readLongLE();
    }

    private static byte[] changed(byte[] archive, int offset, int value) {
        byte[] changed = archive.clone();
        changed[offset] = (byte) value;
        return changed;
    }

    private static void putLongLE(byte[] bytes, int offset, long value) {
        for (int i = 0; i < 8; i++) bytes[offset + i] = (byte) (value >>> (8 * i));
    }

    private static byte[] ascii(String text) {
        return text.getBytes(StandardCharsets.US_ASCII);
    }

    private static void check(boolean condition, String message) {
        checks++;
        if (!condition) throw new AssertionError(message);
    }

    private static void pass(String name) {
        suites++;
        System.out.println("PASS " + name);
    }
}
