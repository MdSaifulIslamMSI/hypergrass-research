package hypergrass;

import java.util.Arrays;

/**
 * Order-1 adaptive binary range coder used as an alternative to Deflate for auxiliary streams.
 * Each byte is coded as an 8-bit symbol whose per-node probabilities are conditioned on the
 * previous byte, which matches the strong local structure of LEB128 varint payloads.
 */
final class RangeCoder {

    private static final int TOP = 1 << 24;
    private static final int PROB_INIT = 1 << 10;
    /** Adaptation rate for the exponential-decay model; larger means slower forgetting. */
    private static final int MOVE_BITS =
            Math.max(1, Integer.getInteger("hypergrass.rangeMoveBits", 5));
    /** "count" uses a KT-style frequency estimator, "decay" the classic LZMA update. */
    private static final boolean COUNT_MODEL =
            "count".equals(System.getProperty("hypergrass.rangeModel", "decay"));

    private RangeCoder() {}

    static byte[] encode(byte[] raw, int len) {
        Encoder enc = new Encoder(len);
        int prev = 0;
        for (int i = 0; i < len; i++) {
            int sym = raw[i] & 0xFF;
            int base = prev << 8;
            int node = 1;
            for (int bit = 7; bit >= 0; bit--) {
                int b = (sym >>> bit) & 1;
                enc.encodeBit(base | node, b);
                node = (node << 1) | b;
            }
            prev = sym;
        }
        return enc.finish();
    }

    static byte[] decode(byte[] comp, int rawLen) {
        Decoder dec = new Decoder(comp, rawLen);
        byte[] out = new byte[rawLen];
        int prev = 0;
        for (int i = 0; i < rawLen; i++) {
            int base = prev << 8;
            int node = 1;
            for (int bit = 0; bit < 8; bit++) {
                int b = dec.decodeBit(base | node);
                node = (node << 1) | b;
            }
            int sym = node & 0xFF;
            out[i] = (byte) sym;
            prev = sym;
        }
        return out;
    }

    private static final class Model {
        final short[] probs = new short[1 << 16];
        final short[] c0 = new short[1 << 16];
        final short[] c1 = new short[1 << 16];

        Model() {
            Arrays.fill(probs, (short) PROB_INIT);
        }

        int prob(int idx) {
            if (!COUNT_MODEL) return probs[idx] & 0xFFFF;
            int zero = c0[idx] & 0xFFFF;
            int one = c1[idx] & 0xFFFF;
            int p = (int) (((long) (zero + 1) << 11) / (zero + one + 2));
            return Math.max(1, Math.min(2047, p));
        }

        void update(int idx, int bit) {
            if (!COUNT_MODEL) {
                int prob = probs[idx] & 0xFFFF;
                prob += bit == 0 ? (2048 - prob) >>> MOVE_BITS : -(prob >>> MOVE_BITS);
                probs[idx] = (short) prob;
                return;
            }
            int zero = c0[idx] & 0xFFFF;
            int one = c1[idx] & 0xFFFF;
            if (bit == 0) zero++; else one++;
            if (zero + one > 16384) {
                zero = (zero + 1) >>> 1;
                one = (one + 1) >>> 1;
            }
            c0[idx] = (short) zero;
            c1[idx] = (short) one;
        }
    }

    private static final class Encoder {
        private final Model model = new Model();
        private final VarIntCodec.ByteOut out;
        private final long[] low = {0L};
        private final int[] range = {-1};
        private int cache;
        private long cacheSize = 1L;

        Encoder(int rawLen) {
            out = new VarIntCodec.ByteOut(rawLen / 2 + 16);
        }

        void encodeBit(int idx, int bit) {
            int prob = model.prob(idx);
            int bound = (range[0] >>> 11) * prob;
            if (bit == 0) {
                range[0] = bound;
            } else {
                low[0] += bound & 0xFFFFFFFFL;
                range[0] -= bound;
            }
            model.update(idx, bit);
            while ((range[0] & 0xFFFFFFFFL) < TOP) {
                range[0] <<= 8;
                shiftLow();
            }
        }

        private void shiftLow() {
            long lo = low[0];
            if ((int) (lo >>> 32) != 0 || lo < 0xFF000000L) {
                int carry = (int) (lo >>> 32);
                int temp = cache;
                do {
                    out.writeByte(temp + carry);
                    temp = 0xFF;
                } while (--cacheSize != 0);
                cache = (int) (lo >>> 24) & 0xFF;
            }
            cacheSize++;
            // The reference algorithm shifts the *32-bit* low (a C unsigned cast), which keeps
            // low < 2^33 so the carry stays 0/1. Shifting the 64-bit value lets carry grow.
            low[0] = ((lo & 0xFFFFFFFFL) << 8) & 0xFFFFFFFFL;
        }

        byte[] finish() {
            for (int i = 0; i < 5; i++) {
                shiftLow();
            }
            return out.toByteArray();
        }
    }

    private static final class Decoder {
        private final Model model = new Model();
        private final byte[] src;
        private int pos;
        private int range = -1;
        private int code;

        Decoder(byte[] src, int rawLen) {
            this.src = src;
            for (int i = 0; i < 5; i++) {
                code = (code << 8) | nextByte();
            }
        }

        private int nextByte() {
            return pos < src.length ? src[pos++] & 0xFF : 0;
        }

        int decodeBit(int idx) {
            int prob = model.prob(idx);
            int bound = (range >>> 11) * prob;
            int bit;
            if ((code & 0xFFFFFFFFL) < (bound & 0xFFFFFFFFL)) {
                range = bound;
                bit = 0;
            } else {
                code -= bound;
                range -= bound;
                bit = 1;
            }
            model.update(idx, bit);
            while ((range & 0xFFFFFFFFL) < TOP) {
                range <<= 8;
                code = (code << 8) | nextByte();
            }
            return bit;
        }
    }
}
