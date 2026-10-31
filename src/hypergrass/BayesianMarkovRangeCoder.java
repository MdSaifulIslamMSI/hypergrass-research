package hypergrass;

import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;

/**
 * Allocation-Free Bayesian Multi-Order Markov + Adaptive Reading-Frame Codon-Phase +
 * Chargaff Reverse-Complement Symmetry Range Coder (BMM-RC).
 *
 * <p>Operates over pure 2-bit DNA/RNA sequences (0=A, 1=C, 2=G, 3=T/U) using 100% deterministic
 * FDLIBM/LUT arithmetic and a 64-bit Subbotin/LZMA carry-buffered range coder (14-bit cumulative
 * precision, TOT = 16384). Supports both High-Compression mode (orders 0..8 + dual codon-phase
 * experts, ~692 KB state) and --lean ultra-low-memory streaming mode (orders 0..6 + codon-phase,
 * ~46 KB state).
 */
public final class BayesianMarkovRangeCoder {

    private static final int PROB_BITS = 14;
    private static final int TOT = 1 << PROB_BITS; // 16384
    private static final long TOP_VALUE = 1L << 24;

    private static final int[] FULL_ORDERS = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    private static final int[] LEAN_ORDERS = {0, 1, 2, 3, 4, 5, 6};

    // Precomputed deterministic log table for quantized probabilities [1 .. 65536]
    private static final int LOG_SCALE = 65536;
    private static final double[] LOG_LUT = new double[LOG_SCALE + 1];

    static {
        LOG_LUT[0] = StrictMath.log(1.0 / LOG_SCALE);
        for (int i = 1; i <= LOG_SCALE; i++) {
            LOG_LUT[i] = StrictMath.log((double) i / LOG_SCALE);
        }
    }

    private BayesianMarkovRangeCoder() {}

    /** Encodes bases[0 .. len-1] into a compact BMM-RC byte array. */
    public static byte[] encode(byte[] bases, int len, boolean leanMode) {
        if (len == 0) {
            return new byte[0];
        }
        ByteOut out = new ByteOut((len >>> 2) + 64);
        AdaptiveGenomicModel model = new AdaptiveGenomicModel(leanMode, len);
        RangeEncoder rc = new RangeEncoder(out);
        int[] freqs = new int[4];

        for (int i = 0; i < len; i++) {
            model.predict(i, freqs);
            int sym = bases[i] & 3;
            int cumLow = 0;
            for (int s = 0; s < sym; s++) {
                cumLow += freqs[s];
            }
            rc.encodeSymbol(cumLow, freqs[sym]);
            model.update(i, sym);
        }
        rc.flush();
        return out.toByteArray();
    }

    /** Decodes exactly len bases from a BMM-RC compressed ByteIn slice. */
    public static byte[] decode(ByteIn in, int len, boolean leanMode) {
        if (len == 0) {
            return new byte[0];
        }
        byte[] out = new byte[len];
        AdaptiveGenomicModel model = new AdaptiveGenomicModel(leanMode, len);
        RangeDecoder rd = new RangeDecoder(in);
        int[] freqs = new int[4];

        for (int i = 0; i < len; i++) {
            model.predict(i, freqs);
            int target = rd.getThreshold();
            int cum = 0;
            int sym = 0;
            while (sym < 3 && cum + freqs[sym] <= target) {
                cum += freqs[sym];
                sym++;
            }
            rd.decodeSymbol(cum, freqs[sym]);
            out[i] = (byte) sym;
            model.update(i, sym);
        }
        return out;
    }

    /** Deterministic multi-expert genomic probability mixer. */
    private static final class AdaptiveGenomicModel {
        private final int[] orders;
        private final int numMarkov;
        private final int numExperts; // numMarkov + 2 codon-phase experts
        private final short[][] tables;
        private final int[] masks;
        private final double[] alphas;
        private final int[] maxSums;

        // Expert numMarkov: local reading-frame wobble expert (pos % 3, order-1 -> 12 contexts)
        private final short[] codon1Table = new short[3 * 4 * 4];
        // Expert numMarkov + 1: reading-frame hexamer/codon expert (pos % 3, order-3 -> 192 contexts)
        private final short[] codon3Table = new short[3 * 64 * 4];

        private final double[] logWeights;
        private final double[] weights;
        private final double[][] expertPreds;

        private long fwdHist = 0L;
        private long rcHist = 0L;

        AdaptiveGenomicModel(boolean leanMode, int seqLen) {
            int[] baseOrders = leanMode ? LEAN_ORDERS : FULL_ORDERS;
            // Include order-9 on large single/anchor sequences (>= 120,000 bp) in non-lean mode
            if (!leanMode && seqLen >= 120_000) {
                this.orders = new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
            } else {
                this.orders = baseOrders;
            }
            this.numMarkov = orders.length;
            this.numExperts = numMarkov + 2;
            this.tables = new short[numMarkov][];
            this.masks = new int[numMarkov];
            this.alphas = new double[numExperts];
            this.maxSums = new int[numMarkov];

            for (int m = 0; m < numMarkov; m++) {
                int k = orders[m];
                int numCtx = 1 << (k << 1);
                tables[m] = new short[numCtx << 2];
                masks[m] = numCtx - 1;
                if (k <= 2) {
                    alphas[m] = 2.0; // scaled by fwd increment = 2
                    maxSums[m] = 16384;
                } else if (k <= 4) {
                    alphas[m] = 1.0;
                    maxSums[m] = 8192;
                } else if (k <= 6) {
                    alphas[m] = 0.40;
                    maxSums[m] = 4096;
                } else {
                    alphas[m] = 0.20;
                    maxSums[m] = 2048;
                }
            }
            alphas[numMarkov] = 1.0;     // codon1
            alphas[numMarkov + 1] = 0.8; // codon3

            this.logWeights = new double[numExperts];
            this.weights = new double[numExperts];
            this.expertPreds = new double[numExperts][4];

            // Mild prior favoring low/mid orders on the very first bases
            for (int m = 0; m < numMarkov; m++) {
                int k = orders[m];
                logWeights[m] = (k <= 3) ? 0.0 : -1.5;
            }
            logWeights[numMarkov] = -0.5;
            logWeights[numMarkov + 1] = -0.8;
        }

        void predict(int pos, int[] outFreqs) {
            // Normalize expert weights via deterministic StrictMath.exp
            double maxLw = logWeights[0];
            for (int m = 1; m < numExperts; m++) {
                if (logWeights[m] > maxLw) maxLw = logWeights[m];
            }
            double wSum = 0.0;
            for (int m = 0; m < numExperts; m++) {
                double d = logWeights[m] - maxLw;
                double w = (d < -30.0) ? 1e-13 : StrictMath.exp(d);
                weights[m] = w;
                wSum += w;
            }
            double invWSum = 1.0 / wSum;

            double p0 = 0.0, p1 = 0.0, p2 = 0.0, p3 = 0.0;

            // Markov experts
            for (int m = 0; m < numMarkov; m++) {
                int k = orders[m];
                double[] ep = expertPreds[m];
                if (pos < k) {
                    ep[0] = 0.25;
                    ep[1] = 0.25;
                    ep[2] = 0.25;
                    ep[3] = 0.25;
                } else {
                    int ctx = (k == 0) ? 0 : (int) (fwdHist & masks[m]);
                    int baseIdx = ctx << 2;
                    short[] tab = tables[m];
                    int c0 = tab[baseIdx];
                    int c1 = tab[baseIdx + 1];
                    int c2 = tab[baseIdx + 2];
                    int c3 = tab[baseIdx + 3];
                    double alpha = alphas[m];
                    double denom = 1.0 / (c0 + c1 + c2 + c3 + 4.0 * alpha);
                    ep[0] = (c0 + alpha) * denom;
                    ep[1] = (c1 + alpha) * denom;
                    ep[2] = (c2 + alpha) * denom;
                    ep[3] = (c3 + alpha) * denom;
                }
                double wm = weights[m] * invWSum;
                p0 += wm * ep[0];
                p1 += wm * ep[1];
                p2 += wm * ep[2];
                p3 += wm * ep[3];
            }

            int phase = pos % 3;

            // Codon-phase order-1 expert (local reading-frame wobble tracker)
            {
                int m = numMarkov;
                double[] ep = expertPreds[m];
                if (pos < 1) {
                    ep[0] = ep[1] = ep[2] = ep[3] = 0.25;
                } else {
                    int ctx1 = (int) (fwdHist & 3L);
                    int baseIdx = ((phase << 2) | ctx1) << 2;
                    int c0 = codon1Table[baseIdx];
                    int c1 = codon1Table[baseIdx + 1];
                    int c2 = codon1Table[baseIdx + 2];
                    int c3 = codon1Table[baseIdx + 3];
                    double alpha = alphas[m];
                    double denom = 1.0 / (c0 + c1 + c2 + c3 + 4.0 * alpha);
                    ep[0] = (c0 + alpha) * denom;
                    ep[1] = (c1 + alpha) * denom;
                    ep[2] = (c2 + alpha) * denom;
                    ep[3] = (c3 + alpha) * denom;
                }
                double wm = weights[m] * invWSum;
                p0 += wm * ep[0];
                p1 += wm * ep[1];
                p2 += wm * ep[2];
                p3 += wm * ep[3];
            }

            // Codon-phase order-3 expert (in-frame codon/dicodon bias)
            {
                int m = numMarkov + 1;
                double[] ep = expertPreds[m];
                if (pos < 3) {
                    ep[0] = ep[1] = ep[2] = ep[3] = 0.25;
                } else {
                    int ctx3 = (int) (fwdHist & 63L);
                    int baseIdx = ((phase * 64) + ctx3) << 2;
                    int c0 = codon3Table[baseIdx];
                    int c1 = codon3Table[baseIdx + 1];
                    int c2 = codon3Table[baseIdx + 2];
                    int c3 = codon3Table[baseIdx + 3];
                    double alpha = alphas[m];
                    double denom = 1.0 / (c0 + c1 + c2 + c3 + 4.0 * alpha);
                    ep[0] = (c0 + alpha) * denom;
                    ep[1] = (c1 + alpha) * denom;
                    ep[2] = (c2 + alpha) * denom;
                    ep[3] = (c3 + alpha) * denom;
                }
                double wm = weights[m] * invWSum;
                p0 += wm * ep[0];
                p1 += wm * ep[1];
                p2 += wm * ep[2];
                p3 += wm * ep[3];
            }

            // Quantize [p0, p1, p2, p3] strictly into integer frequencies >= 1 summing to TOT (16384)
            int f0 = Math.max(1, (int) (p0 * TOT + 0.5));
            int f1 = Math.max(1, (int) (p1 * TOT + 0.5));
            int f2 = Math.max(1, (int) (p2 * TOT + 0.5));
            int f3 = Math.max(1, (int) (p3 * TOT + 0.5));
            int diff = TOT - (f0 + f1 + f2 + f3);
            if (diff != 0) {
                // Apply residual to largest frequency bucket so all frequencies remain >= 1
                if (f0 >= f1 && f0 >= f2 && f0 >= f3) f0 += diff;
                else if (f1 >= f2 && f1 >= f3) f1 += diff;
                else if (f2 >= f3) f2 += diff;
                else f3 += diff;
            }
            outFreqs[0] = f0;
            outFreqs[1] = f1;
            outFreqs[2] = f2;
            outFreqs[3] = f3;
        }

        void update(int pos, int sym) {
            // 1. Update expert weights via decayed log-loss (forgetting = 0.996, learning rate = 1.6)
            final double decay = 0.996;
            final double eta = 1.6;
            for (int m = 0; m < numExperts; m++) {
                int q = (int) (expertPreds[m][sym] * LOG_SCALE + 0.5);
                if (q < 1) q = 1;
                else if (q > LOG_SCALE) q = LOG_SCALE;
                logWeights[m] = decay * logWeights[m] + eta * LOG_LUT[q];
            }

            // Push reverse-complement of current base into top of rcHist BEFORE RC updates
            long prevFwdHist = fwdHist;
            rcHist = (rcHist >>> 2) | ((long) (sym ^ 3) << 60);

            // 2. Update Markov tables (forward weight = 2, plus Chargaff RC weight = 1 for k >= 4)
            for (int m = 0; m < numMarkov; m++) {
                int k = orders[m];
                if (pos >= k) {
                    short[] tab = tables[m];
                    int ctx = (k == 0) ? 0 : (int) (prevFwdHist & masks[m]);
                    incrementBucket(tab, ctx << 2, sym, 2, maxSums[m]);

                    if (k >= 4) {
                        int oldestBase = (int) ((prevFwdHist >>> ((k - 1) << 1)) & 3L);
                        int rcSym = oldestBase ^ 3;
                        int rcCtx = (int) ((rcHist >>> (62 - (k << 1))) & masks[m]);
                        incrementBucket(tab, rcCtx << 2, rcSym, 1, maxSums[m]);
                    }
                }
            }

            // 3. Update local reading-frame codon-phase experts (fast count halving to track ORF frame shifts)
            int phase = pos % 3;
            if (pos >= 1) {
                int ctx1 = (int) (prevFwdHist & 3L);
                incrementBucket(codon1Table, ((phase << 2) | ctx1) << 2, sym, 2, 128);
            }
            if (pos >= 3) {
                int ctx3 = (int) (prevFwdHist & 63L);
                incrementBucket(codon3Table, ((phase * 64) + ctx3) << 2, sym, 2, 256);
            }

            // 4. Advance forward history
            fwdHist = (prevFwdHist << 2) | sym;
        }

        private static void incrementBucket(short[] tab, int baseIdx, int sym, int inc, int maxSum) {
            int ns = tab[baseIdx + sym] + inc;
            tab[baseIdx + sym] = (short) ns;
            int sum = tab[baseIdx] + tab[baseIdx + 1] + tab[baseIdx + 2] + tab[baseIdx + 3];
            if (sum >= maxSum) {
                tab[baseIdx] = (short) (tab[baseIdx] >>> 1);
                tab[baseIdx + 1] = (short) (tab[baseIdx + 1] >>> 1);
                tab[baseIdx + 2] = (short) (tab[baseIdx + 2] >>> 1);
                tab[baseIdx + 3] = (short) (tab[baseIdx + 3] >>> 1);
            }
        }
    }

    /** 64-bit carry-buffered Range Encoder. */
    private static final class RangeEncoder {
        private final ByteOut out;
        private long low = 0L;
        private long range = 0xFFFFFFFFL;
        private int cache = -1;
        private long cacheSize = 0L;

        RangeEncoder(ByteOut out) {
            this.out = out;
        }

        void encodeSymbol(int cumLow, int symFreq) {
            range >>>= PROB_BITS;
            low += (cumLow & 0xFFFFFFFFL) * range;
            range *= symFreq;
            while (range < TOP_VALUE) {
                range <<= 8;
                shiftLow();
            }
        }

        private void shiftLow() {
            int carry = (int) (low >>> 32);
            if (carry != 0 || (low & 0xFF000000L) < 0xFF000000L) {
                if (cache >= 0) {
                    out.writeByte((cache + carry) & 0xFF);
                }
                while (cacheSize > 0) {
                    out.writeByte((0xFF + carry) & 0xFF);
                    cacheSize--;
                }
                cache = (int) ((low >>> 24) & 0xFF);
            } else {
                cacheSize++;
            }
            low = (low << 8) & 0xFFFFFFFFL;
        }

        void flush() {
            for (int i = 0; i < 5; i++) {
                shiftLow();
            }
        }
    }

    /** Matching 64-bit Range Decoder. */
    private static final class RangeDecoder {
        private final ByteIn in;
        private long range = 0xFFFFFFFFL;
        private long code = 0L;
        private long stepRange = 0L;

        RangeDecoder(ByteIn in) {
            this.in = in;
            for (int i = 0; i < 4; i++) {
                code = (code << 8) | (in.hasRemaining() ? (in.readByte() & 0xFFL) : 0L);
            }
        }

        int getThreshold() {
            stepRange = range >>> PROB_BITS;
            long t = code / stepRange;
            return (t >= TOT) ? (TOT - 1) : (int) t;
        }

        void decodeSymbol(int cumLow, int symFreq) {
            code -= (cumLow & 0xFFFFFFFFL) * stepRange;
            range = stepRange * symFreq;
            while (range < TOP_VALUE) {
                range <<= 8;
                code = ((code << 8) | (in.hasRemaining() ? (in.readByte() & 0xFFL) : 0L)) & 0xFFFFFFFFL;
            }
        }
    }
}
