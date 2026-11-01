package hypergrass;

import java.util.Arrays;

/**
 * LZ77 + adaptive binary range coder for auxiliary streams.
 *
 * <p>The variant-gap channel is LEB128 varints over genomic distances. An order-1 byte model
 * leaves ~250 KB of a 330 KB stream, while xz with a 256 KB window reaches 211 KB: the residual
 * redundancy is repeated varint sequences at arbitrary alignment, which a 32 KB Deflate window and
 * an order-1 model both miss.
 *
 * <p>Tokens are literal bytes and (length, distance) matches. Literals reuse the order-1 tree
 * model; match lengths use a 9-bit tree; distances use a slot (log2) tree conditioned on the
 * previous slot plus adaptive mantissa bits, so nearby repeated distances get cheap.
 */
final class LzRangeCoder {

    private static final int TOP = 1 << 24;
    private static final int PROB_INIT = 1 << 10;
    private static final int MOVE_BITS =
            Math.max(1, Integer.getInteger("hypergrass.lzMoveBits", 5));

    private static final int WINDOW = 1 << 18;
    private static final int MIN_MATCH = 4;
    private static final int MAX_MATCH = 512;
    private static final int MAX_CHAIN =
            Math.max(1, Integer.getInteger("hypergrass.lzChain", 64));
    private static final int HASH_BITS = 19;

    private static final int OFF_IS_MATCH = 0;
    private static final int OFF_LITERAL = 2;
    private static final int OFF_LEN = OFF_LITERAL + (256 << 8);
    private static final int OFF_SLOT = OFF_LEN + 1024;
    private static final int OFF_MANT = OFF_SLOT + (32 * 32);
    private static final int MODEL_SIZE = OFF_MANT + (32 * 32);

    private LzRangeCoder() {}

    private static final class Probs {
        final int[] p = new int[MODEL_SIZE];

        Probs() {
            Arrays.fill(p, PROB_INIT);
        }

        int prob(int idx) {
            return p[idx];
        }

        void update(int idx, int bit) {
            int v = p[idx];
            p[idx] = bit == 0 ? v + ((2048 - v) >>> MOVE_BITS) : v - (v >>> MOVE_BITS);
        }
    }

    static byte[] encode(byte[] raw, int len) {
        Encoder enc = new Encoder(len);
        int[] head = new int[1 << HASH_BITS];
        Arrays.fill(head, -1);
        int[] prev = new int[len];
        Arrays.fill(prev, -1);
        int i = 0;
        int prevByte = 0;
        int prevSlot = 0;
        boolean lastWasMatch = false;
        while (i < len) {
            int bestLen = 0;
            int bestDist = 0;
            if (i + MIN_MATCH <= len) {
                int h = hash(raw, i);
                int cand = head[h];
                int chain = 0;
                int limit = Math.min(MAX_MATCH, len - i);
                while (cand >= 0 && i - cand <= WINDOW && chain < MAX_CHAIN) {
                    int l = 0;
                    while (l < limit && raw[cand + l] == raw[i + l]) {
                        l++;
                    }
                    if (l > bestLen) {
                        bestLen = l;
                        bestDist = i - cand;
                        if (l >= limit) {
                            break;
                        }
                    }
                    cand = prev[cand];
                    chain++;
                }
            }
            if (bestLen >= MIN_MATCH) {
                enc.encodeBit(OFF_IS_MATCH + (lastWasMatch ? 1 : 0), 1);
                int lenVal = bestLen - MIN_MATCH;
                int node = 1;
                for (int bit = 8; bit >= 0; bit--) {
                    int b = (lenVal >>> bit) & 1;
                    enc.encodeBit(OFF_LEN + node, b);
                    node = (node << 1) | b;
                }
                int slot = 31 - Integer.numberOfLeadingZeros(bestDist);
                int sctx = prevSlot * 32;
                node = 1;
                for (int bit = 4; bit >= 0; bit--) {
                    int b = (slot >>> bit) & 1;
                    enc.encodeBit(OFF_SLOT + sctx + node, b);
                    node = (node << 1) | b;
                }
                int mant = bestDist - (1 << slot);
                for (int bit = slot - 1; bit >= 0; bit--) {
                    enc.encodeBit(OFF_MANT + slot * 32 + (slot - 1 - bit),
                            (mant >>> bit) & 1);
                }
                prevSlot = slot;
                lastWasMatch = true;
                for (int j = 0; j < bestLen; j++) {
                    if (i + j + MIN_MATCH <= len) {
                        int h = hash(raw, i + j);
                        prev[i + j] = head[h];
                        head[h] = i + j;
                    }
                }
                prevByte = raw[i + bestLen - 1] & 0xFF;
                i += bestLen;
            } else {
                int sym = raw[i] & 0xFF;
                enc.encodeBit(OFF_IS_MATCH + (lastWasMatch ? 1 : 0), 0);
                int base = OFF_LITERAL + (prevByte << 8);
                int node = 1;
                for (int bit = 7; bit >= 0; bit--) {
                    int b = (sym >>> bit) & 1;
                    enc.encodeBit(base + node, b);
                    node = (node << 1) | b;
                }
                if (i + MIN_MATCH <= len) {
                    int h = hash(raw, i);
                    prev[i] = head[h];
                    head[h] = i;
                }
                prevByte = sym;
                lastWasMatch = false;
                i++;
            }
        }
        return enc.finish();
    }

    static byte[] decode(byte[] comp, int rawLen) {
        Decoder dec = new Decoder(comp);
        byte[] out = new byte[rawLen];
        int i = 0;
        int prevByte = 0;
        int prevSlot = 0;
        boolean lastWasMatch = false;
        while (i < rawLen) {
            int isMatch = dec.decodeBit(OFF_IS_MATCH + (lastWasMatch ? 1 : 0));
            if (isMatch == 0) {
                int base = OFF_LITERAL + (prevByte << 8);
                int node = 1;
                for (int bit = 0; bit < 8; bit++) {
                    int b = dec.decodeBit(base + node);
                    node = (node << 1) | b;
                }
                int sym = node & 0xFF;
                out[i++] = (byte) sym;
                prevByte = sym;
                lastWasMatch = false;
            } else {
                int node = 1;
                for (int bit = 0; bit < 9; bit++) {
                    int b = dec.decodeBit(OFF_LEN + node);
                    node = (node << 1) | b;
                }
                int len = (node & 0x1FF) + MIN_MATCH;
                int sctx = prevSlot * 32;
                node = 1;
                for (int bit = 0; bit < 5; bit++) {
                    int b = dec.decodeBit(OFF_SLOT + sctx + node);
                    node = (node << 1) | b;
                }
                int slot = node & 0x1F;
                int mant = 0;
                for (int bit = slot - 1; bit >= 0; bit--) {
                    int b = dec.decodeBit(OFF_MANT + slot * 32 + (slot - 1 - bit));
                    mant |= b << bit;
                }
                int dist = (1 << slot) + mant;
                if (dist <= 0 || dist > i) {
                    throw new IllegalStateException("Invalid LZ distance " + dist + " at " + i);
                }
                if (i + len > rawLen) {
                    len = rawLen - i;
                }
                for (int j = 0; j < len; j++) {
                    out[i] = out[i - dist];
                    i++;
                }
                prevByte = out[i - 1] & 0xFF;
                prevSlot = slot;
                lastWasMatch = true;
            }
        }
        return out;
    }

    private static int hash(byte[] raw, int i) {
        int v = (raw[i] & 0xFF) | ((raw[i + 1] & 0xFF) << 8) | ((raw[i + 2] & 0xFF) << 16);
        v ^= v >>> 13;
        v *= 0x85EBCA6B;
        v ^= v >>> 16;
        return v & ((1 << HASH_BITS) - 1);
    }

    private static final class Encoder {
        private final Probs model = new Probs();
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
        private final Probs model = new Probs();
        private final byte[] src;
        private int pos;
        private int range = -1;
        private int code;

        Decoder(byte[] src) {
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
