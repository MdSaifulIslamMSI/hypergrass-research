package hypergrass;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Read-only encoded-byte accounting for version 1 .hgs archives.
 *
 * <p>Counts include each stream's mode/length prefix unless a separate framing category is
 * provided. In particular, primary and residual payload categories exclude their framing;
 * anchor payload includes the complete encoded pages. The footer's compressed directory is
 * indivisible and includes its auxiliary prefix, while its 12-byte trailer is separate.
 *
 * <p>No auxiliary stream is inflated and no primary payload is decoded. Blocks are scanned
 * sequentially to the trailer's footer offset, rather than using the compressed directory.
 * This checks framing and bounds, not compressed contents, directory entries, or CRC integrity.
 * The Path overload reads the archive into a byte array (and therefore requires it to fit one).
 */
public final class ArchiveByteAccounting {

    private static final int FOOTER_MAGIC = 0x31464748;
    private static final int ANCHOR_PAGE_BASES = 65_536;
    private static final List<String> CATEGORIES = List.of(
            "container.header", "container.footer_directory", "container.trailer",
            "anchor.header", "anchor.global_dictionary", "anchor.page_directory", "anchor.payload",
            "metadata.framing", "metadata.headers", "metadata.line_widths", "metadata.line_endings",
            "metadata.case", "metadata.n_runs", "metadata.exception_gaps", "metadata.exception_symbols",
            "dna.framing", "primary.mode1_payload", "residual.framing", "residual.payload",
            "match.framing", "match.literal_lengths", "match.strands", "match.geometry",
            "variant.counts", "variant.gaps", "variant.differences",
            "donor.sweep_modes", "donor.indices", "donor.masks");

    private ArchiveByteAccounting() {}

    /** Immutable, deterministically ordered totals; absent streams have zero-byte entries. */
    public record Report(long archiveBytes, Map<String, Long> categoryBytes, Map<Integer, Long> modeCounts) {
        public Report {
            categoryBytes = Collections.unmodifiableMap(new LinkedHashMap<>(categoryBytes));
            modeCounts = Collections.unmodifiableMap(new LinkedHashMap<>(modeCounts));
        }

        public long totalBytes() {
            long total = 0;
            for (long bytes : categoryBytes.values()) total += bytes;
            return total;
        }

        public long blockCount() {
            long total = 0;
            for (long count : modeCounts.values()) total += count;
            return total;
        }
    }

    public static Report report(Path archive) throws IOException {
        return report(Files.readAllBytes(archive));
    }

    /** @throws IllegalArgumentException if archive framing is invalid or truncated */
    public static Report report(byte[] archive) {
        return new Parser(Objects.requireNonNull(archive, "archive")).parse();
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 1) {
            throw new IllegalArgumentException(
                    "Usage: java -cp bin hypergrass.ArchiveByteAccounting <archive.hgs>");
        }
        Report report = report(Path.of(args[0]));
        System.out.println("category\tencoded_bytes");
        report.categoryBytes().forEach((category, bytes) -> System.out.println(category + "\t" + bytes));
        System.out.println("total\t" + report.totalBytes());
        System.out.println("archive_bytes\t" + report.archiveBytes());
        System.out.println("mode\tblocks");
        report.modeCounts().forEach((mode, count) -> System.out.println(mode + "\t" + count));
    }

    private static final class Parser {
        private final byte[] bytes;
        private final Map<String, Long> counts = new LinkedHashMap<>();
        private final Map<Integer, Long> modes = new LinkedHashMap<>();

        Parser(byte[] bytes) {
            this.bytes = bytes;
            for (String category : CATEGORIES) counts.put(category, 0L);
            for (int mode = 0; mode <= 6; mode++) modes.put(mode, 0L);
        }

        Report parse() {
            Cursor in = new Cursor(bytes, 0, bytes.length);
            if (in.u8() != 'H' || in.u8() != 'G' || in.u8() != 'S' || in.u8() != 1) {
                throw in.error("Invalid .hgs magic/version");
            }
            int flags = in.u8();
            if ((flags & ~63) != 0) throw in.error("Unknown container flags: " + flags);
            boolean hasAnchor = (flags & 1) != 0;
            boolean tuned = (flags & 16) != 0;
            boolean segmented = (flags & 32) != 0;
            if ((flags & 4) != 0) {
                if (hasAnchor) throw in.error("Compact archive cannot contain an anchor bank");
                if ((flags & 8) != 0) throw in.error("Compact archive cannot set the shared gap flag");
                if (segmented && !tuned) throw in.error("Segmented DNA channel requires the tuned flag");
                if (in.intLE() < 0) throw in.error("Negative compact raw length");
                in.intLE(); // CRC-32, not verified by accounting.
                add("container.header", in.pos);
                block(in, false, tuned);
                in.end();
            } else {
                if (tuned) throw in.error("Tuned DNA channel requires a compact archive");
                if (segmented) throw in.error("Segmented DNA channel requires a tuned compact archive");
                add("container.header", in.pos);
                if (bytes.length < 18) throw in.error("Truncated footer/trailer");
                int trailerOffset = bytes.length - 12;
                Cursor trailer = new Cursor(bytes, trailerOffset, bytes.length);
                long footerOffset = trailer.longLE();
                if (trailer.intLE() != FOOTER_MAGIC) throw trailer.error("Invalid .hgs footer magic");
                if (footerOffset < in.pos || footerOffset >= trailerOffset) {
                    throw trailer.error("Footer offset is outside archive body: " + footerOffset);
                }
                Cursor footer = new Cursor(bytes, (int) footerOffset, trailerOffset);
                if (aux(footer, "container.footer_directory") == 0) {
                    throw footer.error("Empty footer directory");
                }
                footer.end();
                add("container.trailer", 12);
                in = new Cursor(bytes, in.pos, (int) footerOffset);
                if (hasAnchor) anchor(in);
                while (in.pos < in.limit) block(in, hasAnchor, false);
                in.end();
            }
            Report report = new Report(bytes.length, counts, modes);
            if (report.totalBytes() != bytes.length) {
                throw new IllegalStateException("Accounting mismatch: " + report.totalBytes() + " != " + bytes.length);
            }
            return report;
        }

        private void anchor(Cursor in) {
            int start = in.pos;
            int mode = in.u8();
            if (mode > 1) throw in.error("Unknown anchor mode: " + mode);
            int length = in.positive("anchor length");
            int pages = in.positive("anchor page count");
            if (pages != (length - 1) / ANCHOR_PAGE_BASES + 1) {
                throw in.error("Anchor page count does not match anchor length");
            }
            add("anchor.header", in.pos - start);

            start = in.pos;
            in.varInt(); // primary seed length
            int donors = in.varInt();
            if (donors > 0) {
                skipAux(in);
                skipAux(in);
            }
            add("anchor.global_dictionary", in.pos - start);

            int tableStart = in.pos;
            int remainingBases = length;
            for (int page = 0; page < pages; page++) {
                int bases = in.positive("anchor page base length");
                if (bases != Math.min(ANCHOR_PAGE_BASES, remainingBases)) {
                    throw in.error("Invalid anchor page base length");
                }
                remainingBases -= bases;
                in.positive("anchor page byte length");
            }
            add("anchor.page_directory", in.pos - tableStart);
            Cursor table = new Cursor(bytes, tableStart, in.pos);
            for (int page = 0; page < pages; page++) {
                int bases = table.varInt();
                int encodedLength = table.varInt();
                in.require(encodedLength);
                if (mode == 0) {
                    Cursor payload = new Cursor(bytes, in.pos, in.pos + encodedLength);
                    int rawLength = skipAux(payload);
                    if (rawLength != (bases + 3) / 4) {
                        throw payload.error("Packed anchor length does not match page base length");
                    }
                    payload.end();
                }
                in.skip(encodedLength);
                add("anchor.payload", encodedLength);
            }
        }

        private void block(Cursor in, boolean hasAnchor, boolean tuned) {
            int start = in.pos;
            int flags = in.u8();
            if ((flags & ~31) != 0) throw in.error("Unknown metadata flags: " + flags);
            int dnaLength = in.varInt();
            int seqChars = in.varInt();
            in.varInt(); // sequence line count
            if (dnaLength > seqChars) throw in.error("DNA length exceeds sequence character count");
            add("metadata.framing", in.pos - start);
            if ((flags & 2) != 0) {
                start = in.pos;
                in.positive("header count");
                skipAux(in);
                add("metadata.headers", in.pos - start);
            }
            aux(in, "metadata.line_widths");
            aux(in, "metadata.line_endings");
            if ((flags & 4) != 0) {
                start = in.pos;
                in.positive("case interval count");
                in.positive("minimum case interval length");
                skipAux(in);
                skipAux(in);
                add("metadata.case", in.pos - start);
            }
            if ((flags & 8) != 0) {
                start = in.pos;
                in.positive("N-run count");
                skipAux(in);
                add("metadata.n_runs", in.pos - start);
            }
            if ((flags & 16) != 0) {
                start = in.pos;
                in.positive("exception count");
                skipAux(in);
                add("metadata.exception_gaps", in.pos - start);
                aux(in, "metadata.exception_symbols");
            }
            dna(in, dnaLength, hasAnchor, tuned);
        }

        private void dna(Cursor in, int dnaLength, boolean hasAnchor, boolean tuned) {
            int mode = in.u8();
            if (mode > 6) throw in.error("Unknown DNA mode: " + mode);
            if ((mode == 0) != (dnaLength == 0)) throw in.error("DNA mode does not match empty/nonempty length");
            add("dna.framing", 1);
            modes.put(mode, modes.get(mode) + 1);
            if (mode == 0) return;
            if (mode == 1) {
                if (tuned) {
                    int format = in.u8();
                    add("dna.framing", 1);
                    if (format != 1) throw in.error("Unknown tuned DNA format: " + format);
                }
                int start = in.pos;
                int length = in.positive("Mode 1 payload length");
                add("dna.framing", in.pos - start);
                in.skip(length);
                add("primary.mode1_payload", length);
                return;
            }
            if (!hasAnchor) throw in.error("Match mode requires an anchor bank");
            int start = in.pos;
            int matches = in.positive("match count");
            int literals = in.varInt();
            int variants = in.varInt();
            if (matches > dnaLength || literals > dnaLength || variants > dnaLength) {
                throw in.error("Match/literal/variant count exceeds DNA length");
            }
            add("match.framing", in.pos - start);
            aux(in, "match.literal_lengths");
            if (literals > 0) {
                start = in.pos;
                int literalMode = in.u8();
                int payloadLength;
                if (literalMode == 0) {
                    Aux header = auxHeader(in);
                    if (header.rawLength != (literals + 3L) / 4) {
                        throw in.error("Packed literal length does not match literal count");
                    }
                    payloadLength = header.payloadLength;
                } else if (literalMode == 1) {
                    payloadLength = in.positive("residual payload length");
                } else {
                    throw in.error("Unknown residual literal mode: " + literalMode);
                }
                add("residual.framing", in.pos - start);
                in.skip(payloadLength);
                add("residual.payload", payloadLength);
            }
            int strandLength = aux(in, "match.strands");
            if (strandLength != (matches + 7L) / 8) throw in.error("Strand length does not match match count");
            aux(in, "match.geometry");
            if (variants > 0) {
                if (mode >= 3) {
                    aux(in, "donor.sweep_modes");
                    aux(in, "donor.indices");
                    aux(in, "donor.masks");
                }
                aux(in, "variant.counts");
                aux(in, "variant.gaps");
                aux(in, "variant.differences");
            }
        }

        private int aux(Cursor in, String category) {
            int start = in.pos;
            int rawLength = skipAux(in);
            add(category, in.pos - start);
            return rawLength;
        }

        private void add(String category, long length) {
            counts.put(category, counts.get(category) + length);
        }
    }

    private record Aux(int rawLength, int payloadLength) {}

    private static Aux auxHeader(Cursor in) {
        int mode = in.u8();
        if (mode == 0) return new Aux(0, 0);
        if (mode < 1 || mode > 7) throw in.error("Unknown auxiliary mode: " + mode);
        int rawLength = in.positive("auxiliary raw length");
        if (mode == 4) {
            // The digit-plane body is self-delimiting and is consumed here, so the caller has
            // no trailing payload to skip.
            in.positive("plane value count");
            int planes = in.u8();
            for (int k = 0; k < planes; k++) {
                in.positive("plane length");
                skipAux(in);
            }
            return new Aux(rawLength, 0);
        }
        int payloadLength = mode == 1 ? rawLength : in.positive("auxiliary compressed length");
        return new Aux(rawLength, payloadLength);
    }

    private static int skipAux(Cursor in) {
        Aux header = auxHeader(in);
        in.skip(header.payloadLength);
        return header.rawLength;
    }

    private static final class Cursor {
        private final byte[] bytes;
        private int pos;
        private final int limit;

        Cursor(byte[] bytes, int pos, int limit) {
            this.bytes = bytes;
            this.pos = pos;
            this.limit = limit;
        }

        void require(int length) {
            if (length < 0 || length > limit - pos) throw error("Truncated stream: need " + length + " bytes");
        }

        int u8() {
            require(1);
            return bytes[pos++] & 255;
        }

        int varInt() {
            int value = 0;
            for (int shift = 0; shift <= 28; shift += 7) {
                int b = u8();
                if (shift == 28 && (b & 0xF8) != 0) throw error("Overflowing or negative length varint");
                value |= (b & 127) << shift;
                if ((b & 128) == 0) return value;
            }
            throw error("Malformed length varint");
        }

        int positive(String name) {
            int value = varInt();
            if (value == 0) throw error("Zero " + name);
            return value;
        }

        int intLE() {
            require(4);
            return u8() | (u8() << 8) | (u8() << 16) | (u8() << 24);
        }

        long longLE() {
            return (intLE() & 0xFFFFFFFFL) | ((long) intLE() << 32);
        }

        void skip(int length) {
            require(length);
            pos += length;
        }

        void end() {
            if (pos != limit) throw error("Unexpected trailing bytes: " + (limit - pos));
        }

        IllegalArgumentException error(String message) {
            return new IllegalArgumentException(message + " at archive offset " + pos);
        }
    }
}
