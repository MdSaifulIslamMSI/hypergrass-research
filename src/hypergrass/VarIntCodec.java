package hypergrass;

import java.util.Arrays;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * High-speed, zero-dependency binary serialization primitives, unsigned LEB128 / ZigZag
 * variable-length integers, 2-bit DNA bit-packing, and adaptive raw/Deflate auxiliary
 * stream compression for HyperGraSS.
 */
public final class VarIntCodec {

    private VarIntCodec() {}

    /**
     * Rolling window of raw auxiliary bytes shared across blocks by the bounded codec. The
     * encoder uses it as an immutable match dictionary while a block is being built; the raw
     * bytes are appended only after that block's representation has been chosen. The decoder
     * appends the plaintext it just recovered, so both sides see the same prefix history.
     */
    public static final class AuxHistory {
        private final int cap;
        private byte[] buf;
        private int len;

        public AuxHistory(int cap) {
            this.cap = Math.max(1, cap);
            this.buf = new byte[Math.min(this.cap, 1 << 16)];
        }

        public void append(byte[] raw, int rawLen) {
            if (rawLen <= 0) {
                return;
            }
            if (rawLen >= cap) {
                if (buf.length < cap) {
                    buf = Arrays.copyOf(buf, cap);
                }
                System.arraycopy(raw, rawLen - cap, buf, 0, cap);
                len = cap;
                return;
            }
            int keep = Math.min(len, cap - rawLen);
            if (buf.length < cap) {
                buf = Arrays.copyOf(buf, cap);
            }
            if (keep > 0 && keep < len) {
                System.arraycopy(buf, len - keep, buf, 0, keep);
            }
            System.arraycopy(raw, 0, buf, keep, rawLen);
            len = keep + rawLen;
        }

        public byte[] buffer() {
            return buf;
        }

        public int length() {
            return len;
        }
    }

    /** Growable byte output buffer with zero-boxing primitive writers. */
    public static final class ByteOut {
        private byte[] buf;
        private int pos;

        public ByteOut() {
            this(256);
        }

        public ByteOut(int initialCapacity) {
            this.buf = new byte[Math.max(16, initialCapacity)];
            this.pos = 0;
        }

        private void ensure(int extra) {
            if (pos + extra > buf.length) {
                int next = Math.max(buf.length << 1, pos + extra);
                buf = Arrays.copyOf(buf, next);
            }
        }

        public void writeByte(int b) {
            ensure(1);
            buf[pos++] = (byte) b;
        }

        public void writeBytes(byte[] src, int off, int len) {
            if (len <= 0) return;
            ensure(len);
            System.arraycopy(src, off, buf, pos, len);
            pos += len;
        }

        public void writeBytes(byte[] src) {
            writeBytes(src, 0, src.length);
        }

        public void writeIntLE(int v) {
            ensure(4);
            buf[pos++] = (byte) v;
            buf[pos++] = (byte) (v >>> 8);
            buf[pos++] = (byte) (v >>> 16);
            buf[pos++] = (byte) (v >>> 24);
        }

        public void writeLongLE(long v) {
            ensure(8);
            buf[pos++] = (byte) v;
            buf[pos++] = (byte) (v >>> 8);
            buf[pos++] = (byte) (v >>> 16);
            buf[pos++] = (byte) (v >>> 24);
            buf[pos++] = (byte) (v >>> 32);
            buf[pos++] = (byte) (v >>> 40);
            buf[pos++] = (byte) (v >>> 48);
            buf[pos++] = (byte) (v >>> 56);
        }

        public void writeVarInt(int value) {
            ensure(5);
            int v = value;
            while ((v & ~0x7F) != 0) {
                buf[pos++] = (byte) ((v & 0x7F) | 0x80);
                v >>>= 7;
            }
            buf[pos++] = (byte) v;
        }

        public void writeSignedVarInt(int value) {
            int zz = (value << 1) ^ (value >> 31);
            writeVarInt(zz);
        }

        public void writeVarLong(long value) {
            ensure(10);
            long v = value;
            while ((v & ~0x7FL) != 0L) {
                buf[pos++] = (byte) (((int) v & 0x7F) | 0x80);
                v >>>= 7;
            }
            buf[pos++] = (byte) v;
        }

        public int size() {
            return pos;
        }

        public void reset() {
            pos = 0;
        }

        public byte[] toByteArray() {
            return Arrays.copyOf(buf, pos);
        }

        public byte[] internalBuffer() {
            return buf;
        }
    }

    /** Zero-copy cursor over a byte array slice. */
    public static final class ByteIn {
        private final byte[] buf;
        private int pos;
        private final int limit;

        public ByteIn(byte[] buf) {
            this(buf, 0, buf.length);
        }

        public ByteIn(byte[] buf, int offset, int length) {
            this.buf = buf;
            this.pos = offset;
            this.limit = offset + length;
        }

        public boolean hasRemaining() {
            return pos < limit;
        }

        public int remaining() {
            return limit - pos;
        }

        public int position() {
            return pos;
        }

        public void setPosition(int newPos) {
            this.pos = newPos;
        }

        public int readByte() {
            if (pos >= limit) {
                throw new IllegalStateException("Unexpected EOF in ByteIn at " + pos);
            }
            return buf[pos++] & 0xFF;
        }

        public byte[] readBytes(int len) {
            if (len == 0) return new byte[0];
            if (pos + len > limit) {
                throw new IllegalStateException("Unexpected EOF reading " + len + " bytes at " + pos);
            }
            byte[] out = Arrays.copyOfRange(buf, pos, pos + len);
            pos += len;
            return out;
        }

        public void readBytes(byte[] dst, int dstOff, int len) {
            if (len == 0) return;
            if (pos + len > limit) {
                throw new IllegalStateException("Unexpected EOF reading " + len + " bytes at " + pos);
            }
            System.arraycopy(buf, pos, dst, dstOff, len);
            pos += len;
        }

        public int readIntLE() {
            if (pos + 4 > limit) {
                throw new IllegalStateException("Unexpected EOF reading int32 at " + pos);
            }
            int b0 = buf[pos++] & 0xFF;
            int b1 = buf[pos++] & 0xFF;
            int b2 = buf[pos++] & 0xFF;
            int b3 = buf[pos++] & 0xFF;
            return b0 | (b1 << 8) | (b2 << 16) | (b3 << 24);
        }

        public long readLongLE() {
            if (pos + 8 > limit) {
                throw new IllegalStateException("Unexpected EOF reading int64 at " + pos);
            }
            long b0 = buf[pos++] & 0xFFL;
            long b1 = buf[pos++] & 0xFFL;
            long b2 = buf[pos++] & 0xFFL;
            long b3 = buf[pos++] & 0xFFL;
            long b4 = buf[pos++] & 0xFFL;
            long b5 = buf[pos++] & 0xFFL;
            long b6 = buf[pos++] & 0xFFL;
            long b7 = buf[pos++] & 0xFFL;
            return b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
                    | (b4 << 32) | (b5 << 40) | (b6 << 48) | (b7 << 56);
        }

        public int readVarInt() {
            int result = 0;
            int shift = 0;
            while (true) {
                int b = readByte();
                result |= (b & 0x7F) << shift;
                if ((b & 0x80) == 0) {
                    return result;
                }
                shift += 7;
                if (shift > 28) {
                    throw new IllegalStateException("Malformed varint");
                }
            }
        }

        public int readSignedVarInt() {
            int zz = readVarInt();
            return (zz >>> 1) ^ -(zz & 1);
        }

        public long readVarLong() {
            long result = 0L;
            int shift = 0;
            while (true) {
                int b = readByte();
                result |= (long) (b & 0x7F) << shift;
                if ((b & 0x80) == 0) {
                    return result;
                }
                shift += 7;
                if (shift > 63) {
                    throw new IllegalStateException("Malformed varlong");
                }
            }
        }
    }

    /**
     * Writes an auxiliary channel choosing the strictly smallest of raw bytes, headerless Deflate,
     * an order-1 range coder, and a LEB128 digit-plane transform.
     * Mode 0: empty (1 byte total).
     * Mode 1: raw uncompressed.
     * Mode 2: raw Deflate (level 9, nowrap=true).
     * Mode 3: order-1 adaptive range coder (exponential-decay estimator).
     * Mode 4: LEB128 digit planes (self-delimiting continuation bits) coded separately.
     * Mode 5: LZ77 + range coder (256 KB window).
     * Mode 6: order-1 adaptive range coder (frequency-count estimator).
     * Mode 7: LZ77 + range coder with an explicit cross-block history (bounded anchor mode).
     */
    public static void writeCompressedAux(ByteOut out, byte[] raw, int rawLen) {
        writeCompressedAux(out, raw, rawLen, true, null, 0);
    }

    /**
     * Writes an auxiliary stream that may reference up to {@code histLen} bytes of preceding
     * context. The history is read-only here; the caller appends the raw bytes once the block
     * choice is final.
     */
    public static void writeCompressedAux(
            ByteOut out, byte[] raw, int rawLen, byte[] hist, int histLen) {
        writeCompressedAux(out, raw, rawLen, true, hist, histLen);
    }

    private static final boolean PLANES =
            !"false".equals(System.getProperty("hypergrass.planeAux"));
    private static final boolean LZ =
            !"false".equals(System.getProperty("hypergrass.lzAux"));
    private static final boolean AUX_DIAG =
            "true".equals(System.getProperty("hypergrass.auxDiag"));
    private static int diagDepth;

    private static String indent() {
        return "  ".repeat(diagDepth);
    }

    private static void writeCompressedAux(
            ByteOut out, byte[] raw, int rawLen, boolean allowPlanes) {
        writeCompressedAux(out, raw, rawLen, allowPlanes, null, 0);
    }

    private static void writeCompressedAux(
            ByteOut out, byte[] raw, int rawLen, boolean allowPlanes,
            byte[] hist, int histLen) {
        if (rawLen == 0) {
            out.writeByte(0);
            return;
        }
        if (rawLen <= 12) {
            out.writeByte(1);
            out.writeVarInt(rawLen);
            out.writeBytes(raw, 0, rawLen);
            return;
        }
        ByteOut simple = new ByteOut(rawLen / 2 + 64);
        writeSimpleAux(simple, raw, rawLen);
        byte[] bestBuf = Arrays.copyOf(simple.internalBuffer(), simple.size());
        if (allowPlanes && PLANES && rawLen >= 1024) {
            ByteOut planes = tryPlaneAux(raw, rawLen);
            if (AUX_DIAG) {
                System.err.println("AUXPLANE raw=" + rawLen + " simple=" + simple.size()
                        + " planes=" + (planes == null ? -1 : planes.size()));
            }
            if (planes != null && 1 + planes.size() < bestBuf.length) {
                bestBuf = new byte[1 + planes.size()];
                bestBuf[0] = 4;
                System.arraycopy(planes.internalBuffer(), 0, bestBuf, 1, planes.size());
            }
        }
        if (allowPlanes && LZ && rawLen >= 4096) {
            byte[] bestLz = LzRangeCoder.encode(raw, rawLen);
            int bestMode = 5;
            int bestTotal = 1 + varIntSize(rawLen) + varIntSize(bestLz.length) + bestLz.length;
            if (hist != null && histLen > 0) {
                byte[] dictLz = LzRangeCoder.encode(raw, rawLen, hist, histLen);
                int dictTotal = 1 + varIntSize(rawLen) + varIntSize(dictLz.length) + dictLz.length;
                if (AUX_DIAG) {
                    System.err.println("AUXDICT raw=" + rawLen + " hist=" + histLen
                            + " lz=" + bestTotal + " dict=" + dictTotal);
                }
                if (dictTotal < bestTotal) {
                    bestMode = 7;
                    bestLz = dictLz;
                    bestTotal = dictTotal;
                }
            } else if (AUX_DIAG) {
                System.err.println("AUXLZ raw=" + rawLen + " simple=" + simple.size()
                        + " lz=" + bestTotal);
            }
            if (bestTotal < bestBuf.length) {
                out.writeByte(bestMode);
                out.writeVarInt(rawLen);
                out.writeVarInt(bestLz.length);
                out.writeBytes(bestLz, 0, bestLz.length);
                return;
            }
        }
        out.writeBytes(bestBuf, 0, bestBuf.length);
        if (AUX_DIAG) {
            System.err.println(indent() + "AUXW raw=" + rawLen + " mode=" + (bestBuf[0] & 0xFF)
                    + " bytes=" + bestBuf.length);
        }
    }

    /**
     * Rewrites a canonical LEB128 stream into per-digit planes, where bit 7 of every symbol is a
     * "value continues in the next plane" flag. Splitting a value's low digit from its high digits
     * stops high-order noise from disturbing the very predictable low-digit context, which is worth
     * ~15% on variant-gap streams. Returns null when the input is not canonical LEB128.
     */
    private static ByteOut tryPlaneAux(byte[] raw, int rawLen) {
        int[] digitCounts = new int[9];
        int values = 0;
        int maxDigits = 0;
        long v = 0;
        int digits = 0;
        int shift = 0;
        for (int i = 0; i < rawLen; i++) {
            int b = raw[i] & 0xFF;
            if (digits >= 8) {
                return null;
            }
            v |= (long) (b & 0x7F) << shift;
            digits++;
            shift += 7;
            if ((b & 0x80) == 0) {
                if (varIntSizeLong(v) != digits) {
                    return null;
                }
                digitCounts[digits]++;
                values++;
                maxDigits = Math.max(maxDigits, digits);
                v = 0;
                digits = 0;
                shift = 0;
            }
        }
        if (digits != 0 || values < 2) {
            return null;
        }
        int[] planeLen = new int[maxDigits];
        for (int k = 0; k < maxDigits; k++) {
            int len = 0;
            for (int l = k + 1; l <= maxDigits; l++) {
                len += digitCounts[l];
            }
            planeLen[k] = len;
        }
        byte[][] planes = new byte[maxDigits][];
        for (int k = 0; k < maxDigits; k++) {
            planes[k] = new byte[planeLen[k]];
        }
        int[] cursor = new int[maxDigits];
        v = 0;
        digits = 0;
        shift = 0;
        for (int i = 0; i < rawLen; i++) {
            int b = raw[i] & 0xFF;
            v |= (long) (b & 0x7F) << shift;
            digits++;
            shift += 7;
            if ((b & 0x80) == 0) {
                for (int k = 0; k < digits; k++) {
                    int sym = (int) ((v >>> (7 * k)) & 0x7F);
                    if (k < digits - 1) {
                        sym |= 0x80;
                    }
                    planes[k][cursor[k]++] = (byte) sym;
                }
                v = 0;
                digits = 0;
                shift = 0;
            }
        }
        ByteOut body = new ByteOut(rawLen / 2 + 64);
        body.writeVarInt(rawLen);
        body.writeVarInt(values);
        body.writeByte(maxDigits);
        for (int k = 0; k < maxDigits; k++) {
            body.writeVarInt(planeLen[k]);
            diagDepth++;
            try {
                writeCompressedAux(body, planes[k], planeLen[k], false);
            } finally {
                diagDepth--;
            }
        }
        return body;
    }

    private static int varIntSizeLong(long value) {
        int n = 1;
        long v = value;
        while ((v & ~0x7FL) != 0) {
            n++;
            v >>>= 7;
        }
        return n;
    }

    private static void writeSimpleAux(ByteOut out, byte[] raw, int rawLen) {
        Deflater deflater = new Deflater(Deflater.BEST_COMPRESSION, true);
        int compLen;
        byte[] comp;
        try {
            deflater.setInput(raw, 0, rawLen);
            deflater.finish();
            comp = new byte[rawLen + 32];
            compLen = deflater.deflate(comp);
            if (!deflater.finished()) {
                compLen = rawLen;
            }
        } finally {
            deflater.end();
        }
        byte[] ranged = null;
        int rangedLen = 0;
        boolean countModel = RangeCoder.countModelDefault();
        if (rawLen >= 256 && !"false".equals(System.getProperty("hypergrass.rangeAux"))) {
            ranged = RangeCoder.encode(raw, rawLen, countModel);
            rangedLen = ranged.length;
        }
        if (ranged != null && rangedLen + varIntSize(rangedLen) + varIntSize(rawLen) + 1 < compLen
                && rangedLen < rawLen) {
            // Modes 3 and 6 differ only in the probability estimator, and the choice is recorded
            // here so an archive never depends on a JVM flag to decode.
            out.writeByte(countModel ? 6 : 3);
            out.writeVarInt(rawLen);
            out.writeVarInt(rangedLen);
            out.writeBytes(ranged, 0, rangedLen);
        } else if (compLen + varIntSize(compLen) < rawLen) {
            out.writeByte(2);
            out.writeVarInt(rawLen);
            out.writeVarInt(compLen);
            out.writeBytes(comp, 0, compLen);
        } else {
            out.writeByte(1);
            out.writeVarInt(rawLen);
            out.writeBytes(raw, 0, rawLen);
        }
    }

    public static void writeCompressedAux(ByteOut out, ByteOut payload) {
        writeCompressedAux(out, payload.internalBuffer(), payload.size());
    }

    public static void writeCompressedAux(
            ByteOut out, ByteOut payload, byte[] hist, int histLen) {
        writeCompressedAux(out, payload.internalBuffer(), payload.size(), true, hist, histLen);
    }

    public static byte[] readCompressedAux(ByteIn in) {
        int mode = in.readByte();
        if (mode == 0) {
            return new byte[0];
        }
        int rawLen = in.readVarInt();
        if (AUX_DIAG) {
            System.err.println(indent() + "AUXR mode=" + mode + " raw=" + rawLen);
        }
        diagDepth++;
        try {
            return readCompressedAuxBody(in, mode, rawLen);
        } finally {
            diagDepth--;
        }
    }

    /**
     * Reads a stream that may carry a mode-7 block and appends the recovered plaintext to
     * {@code hist}. Only the gap channel uses this entry point; nested plane streams go through
     * the plain method and never touch the shared history.
     */
    public static byte[] readCompressedAux(ByteIn in, AuxHistory hist) {
        int mode = in.readByte();
        if (mode == 0) {
            return new byte[0];
        }
        int rawLen = in.readVarInt();
        if (AUX_DIAG) {
            System.err.println(indent() + "AUXR mode=" + mode + " raw=" + rawLen);
        }
        diagDepth++;
        byte[] raw;
        try {
            if (mode == 7) {
                if (hist == null) {
                    throw new IllegalStateException("Mode 7 aux stream requires a history");
                }
                int compLen = in.readVarInt();
                raw = LzRangeCoder.decode(in.readBytes(compLen), rawLen,
                        hist.buffer(), hist.length());
            } else {
                raw = readCompressedAuxBody(in, mode, rawLen);
            }
        } finally {
            diagDepth--;
        }
        if (hist != null) {
            hist.append(raw, raw.length);
        }
        return raw;
    }

    private static byte[] readCompressedAuxBody(ByteIn in, int mode, int rawLen) {
        if (mode == 1) {
            return in.readBytes(rawLen);
        }
        if (mode == 3) {
            int compLen = in.readVarInt();
            return RangeCoder.decode(in.readBytes(compLen), rawLen, false);
        }
        if (mode == 6) {
            int compLen = in.readVarInt();
            return RangeCoder.decode(in.readBytes(compLen), rawLen, true);
        }
        if (mode == 5) {
            int compLen = in.readVarInt();
            return LzRangeCoder.decode(in.readBytes(compLen), rawLen);
        }
        if (mode == 4) {
            int values = in.readVarInt();
            int planes = in.readByte();
            byte[][] digits = new byte[planes][];
            for (int k = 0; k < planes; k++) {
                int len = in.readVarInt();
                digits[k] = readCompressedAux(in);
                if (digits[k].length != len) {
                    throw new IllegalStateException("Plane " + k + " length mismatch");
                }
            }
            if (AUX_DIAG) {
                StringBuilder sb = new StringBuilder();
                for (byte[] d : digits) {
                    sb.append(d.length).append(' ');
                }
                System.err.println(indent() + "AUXR4 values=" + values + " planes=" + planes
                        + " lens=" + sb);
            }
            ByteOut out = new ByteOut(rawLen);
            int[] cursor = new int[planes];
            for (int i = 0; i < values; i++) {
                long v = 0L;
                int shift = 0;
                for (int k = 0; k < planes; k++) {
                    if (cursor[k] >= digits[k].length) {
                        StringBuilder sb = new StringBuilder();
                        for (int z = 0; z < planes; z++) {
                            sb.append(digits[z].length).append('/').append(cursor[z]).append(' ');
                        }
                        throw new IllegalStateException("Plane " + k + " exhausted at value " + i
                                + "/" + values + " (len/cursor: " + sb + ")");
                    }
                    int sym = digits[k][cursor[k]++] & 0xFF;
                    v |= (long) (sym & 0x7F) << shift;
                    shift += 7;
                    if ((sym & 0x80) == 0) {
                        break;
                    }
                }
                while ((v & ~0x7FL) != 0) {
                    out.writeByte((int) (v & 0x7F) | 0x80);
                    v >>>= 7;
                }
                out.writeByte((int) v);
            }
            if (out.size() != rawLen) {
                throw new IllegalStateException("Plane aux decode mismatch: " + out.size() + " vs " + rawLen);
            }
            return Arrays.copyOf(out.internalBuffer(), out.size());
        }
        if (mode == 2) {
            int compLen = in.readVarInt();
            byte[] comp = in.readBytes(compLen);
            byte[] out = new byte[rawLen];
            Inflater inflater = new Inflater(true);
            try {
                inflater.setInput(comp, 0, compLen);
                int off = 0;
                while (off < rawLen && !inflater.finished()) {
                    int n = inflater.inflate(out, off, rawLen - off);
                    if (n == 0 && inflater.needsInput()) {
                        break;
                    }
                    off += n;
                }
                if (off != rawLen) {
                    throw new IllegalStateException("Deflate length mismatch: expected " + rawLen + ", got " + off);
                }
                return out;
            } catch (DataFormatException e) {
                throw new IllegalStateException("Corrupt Deflate auxiliary stream", e);
            } finally {
                inflater.end();
            }
        }
        throw new IllegalStateException("Unknown aux compression mode: " + mode);
    }

    public static int varIntSize(int value) {
        int v = value;
        int n = 1;
        while ((v & ~0x7F) != 0) {
            n++;
            v >>>= 7;
        }
        return n;
    }

    /** Packs 2-bit DNA bases (0=A, 1=C, 2=G, 3=T) at 4 bases per byte. */
    public static byte[] pack2Bit(byte[] bases, int off, int len) {
        byte[] packed = new byte[(len + 3) >>> 2];
        int fullBytes = len >>> 2;
        int p = off;
        for (int i = 0; i < fullBytes; i++) {
            packed[i] = (byte) ((bases[p] & 3)
                    | ((bases[p + 1] & 3) << 2)
                    | ((bases[p + 2] & 3) << 4)
                    | ((bases[p + 3] & 3) << 6));
            p += 4;
        }
        int rem = len & 3;
        if (rem > 0) {
            int last = 0;
            for (int r = 0; r < rem; r++) {
                last |= (bases[p + r] & 3) << (r << 1);
            }
            packed[fullBytes] = (byte) last;
        }
        return packed;
    }

    /** Unpacks 2-bit DNA bases (0=A, 1=C, 2=G, 3=T) from packed byte array. */
    public static byte[] unpack2Bit(byte[] packed, int len) {
        byte[] bases = new byte[len];
        int fullBytes = len >>> 2;
        int p = 0;
        for (int i = 0; i < fullBytes; i++) {
            int b = packed[i] & 0xFF;
            bases[p] = (byte) (b & 3);
            bases[p + 1] = (byte) ((b >>> 2) & 3);
            bases[p + 2] = (byte) ((b >>> 4) & 3);
            bases[p + 3] = (byte) ((b >>> 6) & 3);
            p += 4;
        }
        int rem = len & 3;
        if (rem > 0) {
            int b = packed[fullBytes] & 0xFF;
            for (int r = 0; r < rem; r++) {
                bases[p + r] = (byte) ((b >>> (r << 1)) & 3);
            }
        }
        return bases;
    }
}
