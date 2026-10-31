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
     * Writes an auxiliary channel choosing the strictly smaller of raw bytes vs headerless Deflate.
     * Mode 0: empty (1 byte total).
     * Mode 1: raw uncompressed.
     * Mode 2: raw Deflate (level 9, nowrap=true).
     */
    public static void writeCompressedAux(ByteOut out, byte[] raw, int rawLen) {
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
        Deflater deflater = new Deflater(Deflater.BEST_COMPRESSION, true);
        try {
            deflater.setInput(raw, 0, rawLen);
            deflater.finish();
            byte[] comp = new byte[rawLen + 32];
            int compLen = deflater.deflate(comp);
            if (deflater.finished() && compLen + varIntSize(compLen) < rawLen) {
                out.writeByte(2);
                out.writeVarInt(rawLen);
                out.writeVarInt(compLen);
                out.writeBytes(comp, 0, compLen);
            } else {
                out.writeByte(1);
                out.writeVarInt(rawLen);
                out.writeBytes(raw, 0, rawLen);
            }
        } finally {
            deflater.end();
        }
    }

    public static void writeCompressedAux(ByteOut out, ByteOut payload) {
        writeCompressedAux(out, payload.internalBuffer(), payload.size());
    }

    public static byte[] readCompressedAux(ByteIn in) {
        int mode = in.readByte();
        if (mode == 0) {
            return new byte[0];
        }
        int rawLen = in.readVarInt();
        if (mode == 1) {
            return in.readBytes(rawLen);
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
