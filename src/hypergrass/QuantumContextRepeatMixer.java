package hypergrass;

import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.util.Arrays;

/**
 * QuantumContextRepeatMixer (QCRM): 100% Deterministic Integer Fixed-Point Neural + Hierarchical
 * Markov + Spaced-Seed + Bidirectional Kimura-Repeat Context-Mixing Range Coder.
 *
 * <p>Architectural breakthroughs over first-generation BMM-RC, GeCo3, JARVIS3, DCom, and GraSS:
 * <ol>
 *   <li><b>Zero Floating-Point Hardware Drift:</b> Every probability, logarithm, exponential,
 *       neural logit, and SGD gradient is computed in pure 64-bit integer fixed-point arithmetic
 *       (zero {@code float}/{@code double} anywhere in the class), guaranteeing 100% bit-exact
 *       cross-CPU decompression across x86-64, ARM64, and RISC-V (solving GeCo3's FMA failure).</li>
 *   <li><b>Hierarchical Pitman-Yor / Dirichlet Base Measure:</b> High-order Markov, sparse-hash,
 *       and spaced-seed experts back off smoothly to the live Order-2 dinucleotide posterior
 *       instead of uniform 0.25, adapting instantaneously to extreme AT/GC skew (e.g. Plasmodium).</li>
 *   <li><b>High-Order 2-Way Cache-Associative Sparse Hash Contexts (k=12, k=15):</b> Captures long
 *       dispersed regulatory and protein motifs across megabase genomes in compact SWAR 64-bit slots.</li>
 *   <li><b>Wobble-Skipping Codon-Skeleton &amp; SNP-Bridging Spaced Seeds:</b> Conditions on
 *       {@code 110_110_110_11} (amino-acid skeleton across 4 codons ignoring synonymous 3rd-base
 *       wobble) and {@code 1111_0_1111} (1-SNP-tolerant 9-mer).</li>
 *   <li><b>Multi-Head Bidirectional Kimura-XOR Repeat Engine (Direct + &amp; Inverted -):</b>
 *       Tracks simultaneous homologous repeat cursors on both strands with 4-state Kimura
 *       transition/transversion XOR posteriors ({@code sym ^ templateBase}).</li>
 *   <li><b>Regime-Gated Online Multinomial Logistic Neural Mixer + Bayesian Meta-Pool:</b>
 *       Combines log-linear neural mixing (zero inactive-expert dilution!) with fast Bayesian
 *       exponential mixing across 8 genomic regime buckets.</li>
 * </ol>
 */
public final class QuantumContextRepeatMixer {

    private static final int PROB_BITS = 15;
    private static final int TOT = 1 << PROB_BITS; // 32768 range-coder scale
    private static final long TOP_VALUE = 1L << 24;

    private static final int P_BITS = 14;
    private static final int P_SCALE = 1 << P_BITS; // 16384 per-expert probability scale

    // Pure-integer precomputed natural log LUT: LN_LUT[p] = round(4096 * ln(p / 16384)) for p in [1..16384]
    private static final int LOG_Q = 4096; // Q12
    private static final int[] LN_LUT = new int[P_SCALE + 1];

    // Pure-integer precomputed exp LUT: EXP_LUT[x] = round(65536 * exp(-x / 4096)) for x in [0..65536]
    private static final int EXP_MAX_X = 65536; // covers down to exp(-16) ~= 1.1e-7
    private static final int[] EXP_NEG_LUT = new int[EXP_MAX_X + 1];

    static {
        initIntegerTables();
    }

    private static void initIntegerTables() {
        // 1. Exact integer bit-by-bit log2 computation (zero floating point!)
        // ln(x) = log2(x) * ln(2), with LN2_Q16 = 45426 (45426 / 65536 = 0.69314575...)
        final long LN2_Q16 = 45426L;
        for (int p = 1; p <= P_SCALE; p++) {
            int intLog2 = 0;
            long y = p;
            while (y < P_SCALE) {
                y <<= 1;
                intLog2--;
            }
            // Now y in [16384, 32768]; extract 16 fractional bits of log2(y / 16384)
            int fracBits = 0;
            for (int b = 15; b >= 0; b--) {
                y = (y * y) >>> P_BITS;
                if (y >= (P_SCALE << 1)) {
                    fracBits |= (1 << b);
                    y >>>= 1;
                }
            }
            long log2Q16 = ((long) intLog2 << 16) + fracBits; // <= 0
            long lnQ12 = (log2Q16 * LN2_Q16) >> 20; // Q16 * Q16 -> Q32 >> 20 -> Q12
            LN_LUT[p] = (int) lnQ12;
        }
        LN_LUT[0] = LN_LUT[1] - 2839;

        // 2. Pure-integer exp(-x / 4096) in Q16 [0..65536] via exp(-t) = 2^(-t / ln2)
        // 1 / ln(2) in Q16 = 94548 (94548 / 65536 = 1.442687988...)
        final long INV_LN2_Q16 = 94548L;
        for (int x = 0; x <= EXP_MAX_X; x++) {
            long uQ16 = ((long) x * INV_LN2_Q16) >>> 12; // x/4096 * (1/ln2) in Q16
            int k = (int) (uQ16 >>> 16); // integer part of negative exponent
            if (k >= 16) {
                EXP_NEG_LUT[x] = 1;
                continue;
            }
            // r in [0, 65535] represents fractional part in [0, 1); compute 2^(-r/65536) via 5-term minimax/Taylor
            // Let s = r * ln(2) in Q16: exp(-s) = 1 - s + s^2/2 - s^3/6 + s^4/24 - s^5/120
            long r = uQ16 & 0xFFFFL;
            long s = (r * LN2_Q16) >>> 16; // s in [0, 45426] (Q16)
            long s2 = (s * s) >>> 16;
            long s3 = (s2 * s) >>> 16;
            long s4 = (s3 * s) >>> 16;
            long s5 = (s4 * s) >>> 16;
            long expNegS = 65536L - s + (s2 >>> 1) - (s3 / 6L) + (s4 / 24L) - (s5 / 120L);
            int val = (int) (expNegS >>> k);
            EXP_NEG_LUT[x] = Math.max(1, val);
        }
    }

    private QuantumContextRepeatMixer() {}

    /** Encodes 2-bit DNA/RNA sequence bases[0 .. len-1] into a compact QCRM byte array. */
    public static byte[] encode(byte[] bases, int len, boolean leanMode) {
        if (len == 0) {
            return new byte[0];
        }
        ByteOut out = new ByteOut((len >>> 2) + 64);
        QcrmModel model = new QcrmModel(leanMode, len);
        RangeEncoder rc = new RangeEncoder(out);
        int[] freqs = new int[4];

        for (int i = 0; i < len; i++) {
            model.predict(bases, i, freqs);
            int sym = bases[i] & 3;
            int cumLow = 0;
            for (int s = 0; s < sym; s++) {
                cumLow += freqs[s];
            }
            rc.encodeSymbol(cumLow, freqs[sym]);
            model.update(bases, i, sym);
        }
        rc.flush();
        return out.toByteArray();
    }

    /** Decodes exactly len 2-bit bases from a QCRM compressed ByteIn stream. */
    public static byte[] decode(ByteIn in, int len, boolean leanMode) {
        if (len == 0) {
            return new byte[0];
        }
        byte[] out = new byte[len];
        QcrmModel model = new QcrmModel(leanMode, len);
        RangeDecoder rd = new RangeDecoder(in);
        int[] freqs = new int[4];

        for (int i = 0; i < len; i++) {
            model.predict(out, i, freqs);
            int target = rd.getThreshold();
            int cum = 0;
            int sym = 0;
            while (sym < 3 && cum + freqs[sym] <= target) {
                cum += freqs[sym];
                sym++;
            }
            rd.decodeSymbol(cum, freqs[sym]);
            out[i] = (byte) sym;
            model.update(out, i, sym);
        }
        return out;
    }

    /**
     * Compact 2-way cache-associative sparse high-order context table packing
     * [32-bit verification tag | c3(8b) | c2(8b) | c1(8b) | c0(8b)] into each 64-bit long.
     */
    private static final class SparseContextTable {
        private final long[] slots;
        private final int mask;

        SparseContextTable(int numSlotsPowerOf2) {
            this.slots = new long[numSlotsPowerOf2];
            this.mask = numSlotsPowerOf2 - 1;
        }

        int predict(long ctxKey, int[] priorP, int tau, int[] outP) {
            long mixed = mix64(ctxKey);
            int idx0 = (int) mixed & mask;
            int tag = (int) (mixed >>> 32);
            if (tag == 0) tag = 1; // 0 means empty slot

            long entry = slots[idx0];
            if ((int) (entry >>> 32) != tag) {
                entry = slots[idx0 ^ 1];
                if ((int) (entry >>> 32) != tag) {
                    System.arraycopy(priorP, 0, outP, 0, 4);
                    return 0;
                }
            }
            int c0 = (int) (entry & 0xFFL);
            int c1 = (int) ((entry >>> 8) & 0xFFL);
            int c2 = (int) ((entry >>> 16) & 0xFFL);
            int c3 = (int) ((entry >>> 24) & 0xFFL);
            int sum = c0 + c1 + c2 + c3;
            if (sum == 0) {
                System.arraycopy(priorP, 0, outP, 0, 4);
                return 0;
            }
            long denom = 4L * sum + tau; // tau in quarter-count units
            int p0 = Math.max(1, (int) (((long) c0 * 4L * P_SCALE + (long) tau * priorP[0] + (denom >> 1)) / denom));
            int p1 = Math.max(1, (int) (((long) c1 * 4L * P_SCALE + (long) tau * priorP[1] + (denom >> 1)) / denom));
            int p2 = Math.max(1, (int) (((long) c2 * 4L * P_SCALE + (long) tau * priorP[2] + (denom >> 1)) / denom));
            int p3 = Math.max(1, P_SCALE - (p0 + p1 + p2));
            normalizeP(p0, p1, p2, p3, outP);
            return sum;
        }

        void update(long ctxKey, int sym, int inc) {
            long mixed = mix64(ctxKey);
            int idx0 = (int) mixed & mask;
            int idx1 = idx0 ^ 1;
            int tag = (int) (mixed >>> 32);
            if (tag == 0) tag = 1;

            long e0 = slots[idx0];
            int t0 = (int) (e0 >>> 32);
            if (t0 == tag || t0 == 0) {
                slots[idx0] = bumpEntry(t0 == 0 ? (((long) tag) << 32) : e0, sym, inc);
                return;
            }
            long e1 = slots[idx1];
            int t1 = (int) (e1 >>> 32);
            if (t1 == tag || t1 == 0) {
                slots[idx1] = bumpEntry(t1 == 0 ? (((long) tag) << 32) : e1, sym, inc);
                return;
            }
            // LFU-lite replacement: only evict if existing bucket has weak count (<= 6)
            int s0 = sumCounts(e0);
            int s1 = sumCounts(e1);
            if (s0 <= s1 && s0 <= 6) {
                slots[idx0] = bumpEntry(((long) tag) << 32, sym, inc);
            } else if (s1 < s0 && s1 <= 6) {
                slots[idx1] = bumpEntry(((long) tag) << 32, sym, inc);
            }
        }

        private static int sumCounts(long e) {
            return (int) ((e & 0xFFL) + ((e >>> 8) & 0xFFL) + ((e >>> 16) & 0xFFL) + ((e >>> 24) & 0xFFL));
        }

        private static long bumpEntry(long entry, int sym, int inc) {
            int shift = sym << 3;
            int cur = (int) ((entry >>> shift) & 0xFFL);
            int next = cur + inc;
            if (next >= 248) {
                // SWAR parallel halving of all 4 8-bit counters in 1 instruction!
                long tagPart = entry & 0xFFFFFFFF00000000L;
                long halved = ((entry & 0xFFFFFFFFL) >>> 1) & 0x7F7F7F7FL;
                cur = (int) ((halved >>> shift) & 0xFFL);
                next = cur + Math.max(1, inc >> 1);
                entry = tagPart | halved;
            }
            long cleared = entry & ~(0xFFL << shift);
            return cleared | (((long) next & 0xFFL) << shift);
        }

        private static long mix64(long z) {
            z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
            z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
            return z ^ (z >>> 31);
        }
    }

    /** Bidirectional Repeat Position Cache (2-way recent-position ring per bucket for + strand, 1-way for -). */
    private static final class RepeatSeedIndex {
        private final int[] fwdPos0;
        private final int[] fwdPos1;
        private final short[] fwdTag0;
        private final short[] fwdTag1;

        private final int[] revCursor;
        private final short[] revTag;
        private final int mask;

        RepeatSeedIndex(int bucketsPowerOf2) {
            this.mask = bucketsPowerOf2 - 1;
            this.fwdPos0 = new int[bucketsPowerOf2];
            this.fwdPos1 = new int[bucketsPowerOf2];
            this.fwdTag0 = new short[bucketsPowerOf2];
            this.fwdTag1 = new short[bucketsPowerOf2];
            this.revCursor = new int[bucketsPowerOf2];
            this.revTag = new short[bucketsPowerOf2];
            Arrays.fill(this.fwdPos0, -1);
            Arrays.fill(this.fwdPos1, -1);
            Arrays.fill(this.revCursor, -1);
        }

        void lookupFwd(long kmer, int[] outTwoPositions) {
            long h = mixSeed(kmer);
            int slot = (int) h & mask;
            short tag = (short) ((h >>> 32) | 1);
            outTwoPositions[0] = (fwdTag0[slot] == tag) ? fwdPos0[slot] : -1;
            outTwoPositions[1] = (fwdTag1[slot] == tag) ? fwdPos1[slot] : -1;
        }

        int lookupRev(long kmer) {
            long h = mixSeed(kmer);
            int slot = (int) h & mask;
            short tag = (short) ((h >>> 32) | 1);
            return (revTag[slot] == tag) ? revCursor[slot] : -1;
        }

        void insertFwd(long kmer, int nextTargetPos) {
            long h = mixSeed(kmer);
            int slot = (int) h & mask;
            short tag = (short) ((h >>> 32) | 1);
            if (fwdTag0[slot] == tag) {
                // Keep previous occurrence in slot 1 if sufficiently separated (>= 4 bp)
                if (nextTargetPos - fwdPos0[slot] >= 4) {
                    fwdPos1[slot] = fwdPos0[slot];
                    fwdTag1[slot] = tag;
                }
                fwdPos0[slot] = nextTargetPos;
            } else {
                fwdPos1[slot] = fwdPos0[slot];
                fwdTag1[slot] = fwdTag0[slot];
                fwdPos0[slot] = nextTargetPos;
                fwdTag0[slot] = tag;
            }
        }

        void insertRev(long rcKmer, int predCursorBeforeSeed) {
            long h = mixSeed(rcKmer);
            int slot = (int) h & mask;
            short tag = (short) ((h >>> 32) | 1);
            revCursor[slot] = predCursorBeforeSeed;
            revTag[slot] = tag;
        }

        private static long mixSeed(long x) {
            x ^= (x >>> 23);
            x *= 0xFF51AFD7ED558CCDL;
            x ^= (x >>> 29);
            return x;
        }
    }

    /** Active Bidirectional Kimura-XOR Repeat Hypothesis Head. */
    private static final class RepeatHead {
        boolean active = false;
        int cursor = 0;
        int dir = 1; // +1 = direct repeat, -1 = inverted (reverse-complement) repeat
        int misses = 0;
        int maxMisses = 6;
        int streak = 0;
        int totalMatched = 0;
        // 4-state Kimura XOR counts: [0]=exact match, [1]=transv1, [2]=Kimura transition, [3]=transv2
        final int[] xorCnt = new int[4];

        void activate(int startCursor, int direction, int initMatchWeight, int maxConsecMiss) {
            this.active = true;
            this.cursor = startCursor;
            this.dir = direction;
            this.misses = 0;
            this.maxMisses = maxConsecMiss;
            this.streak = 0;
            this.totalMatched = 0;
            this.xorCnt[0] = initMatchWeight;
            this.xorCnt[1] = 2;
            this.xorCnt[2] = 4; // Biological Kimura transition prior (2x transversion!)
            this.xorCnt[3] = 2;
        }

        boolean predict(byte[] history, int pos, int[] fallbackP, int[] outP) {
            if (!active || cursor < 0 || cursor >= pos) {
                active = false;
                System.arraycopy(fallbackP, 0, outP, 0, 4);
                return false;
            }
            int tpl = (dir == 1) ? (history[cursor] & 3) : ((history[cursor] ^ 3) & 3);
            int cMatch = xorCnt[0] + Math.min(96, streak << 1);
            int c1 = xorCnt[1];
            int c2 = xorCnt[2];
            int c3 = xorCnt[3];
            int sum = cMatch + c1 + c2 + c3;

            int[] px = new int[4];
            px[0] = Math.max(1, (cMatch * P_SCALE + (sum >> 1)) / sum);
            px[1] = Math.max(1, (c1 * P_SCALE + (sum >> 1)) / sum);
            px[2] = Math.max(1, (c2 * P_SCALE + (sum >> 1)) / sum);
            px[3] = Math.max(1, (c3 * P_SCALE + (sum >> 1)) / sum);

            // Map XOR residual distribution back to nucleotide {A=0, C=1, G=2, T=3}: b = tpl ^ xor
            int p0 = px[0 ^ tpl];
            int p1 = px[1 ^ tpl];
            int p2 = px[2 ^ tpl];
            int p3 = px[3 ^ tpl];
            normalizeP(p0, p1, p2, p3, outP);
            return true;
        }

        void step(byte[] history, int pos, int actualSym) {
            if (!active) return;
            if (cursor < 0 || cursor >= pos) {
                active = false;
                return;
            }
            int tpl = (dir == 1) ? (history[cursor] & 3) : ((history[cursor] ^ 3) & 3);
            int xor = (actualSym ^ tpl) & 3;
            xorCnt[xor] += 8;
            int sum = xorCnt[0] + xorCnt[1] + xorCnt[2] + xorCnt[3];
            if (sum >= 256) {
                xorCnt[0] = Math.max(2, xorCnt[0] >>> 1);
                xorCnt[1] = Math.max(1, xorCnt[1] >>> 1);
                xorCnt[2] = Math.max(1, xorCnt[2] >>> 1);
                xorCnt[3] = Math.max(1, xorCnt[3] >>> 1);
            }
            if (xor == 0) {
                misses = 0;
                streak++;
                totalMatched++;
            } else {
                misses++;
                streak = 0;
                // Retire if consecutive misses exceed threshold (e.g. indel or end of repeat)
                if (misses > maxMisses) {
                    active = false;
                    return;
                }
            }
            cursor += dir;
            if (cursor < 0 || cursor > pos) {
                active = false;
            }
        }
    }
    static void normalizeP(int p0, int p1, int p2, int p3, int[] outP) {
        if (p0 < 1) p0 = 1;
        if (p1 < 1) p1 = 1;
        if (p2 < 1) p2 = 1;
        if (p3 < 1) p3 = 1;
        int diff = P_SCALE - (p0 + p1 + p2 + p3);
        if (diff != 0) {
            if (p0 >= p1 && p0 >= p2 && p0 >= p3) p0 += diff;
            else if (p1 >= p2 && p1 >= p3) p1 += diff;
            else if (p2 >= p3) p2 += diff;
            else p3 += diff;
        }
        outP[0] = p0;
        outP[1] = p1;
        outP[2] = p2;
        outP[3] = p3;
    }

    /** Full Deterministic Multi-Expert + Repeat + Neural/Bayesian Context-Mixing Model. */
    private static final class QcrmModel {
        private static final int K_FINE = 13;
        private static final long MASK_FINE = (1L << (K_FINE << 1)) - 1L;
        private static final int K_LONG = 18;
        private static final long MASK_LONG = (1L << (K_LONG << 1)) - 1L;
        private static final int K_LEAN = 12;
        private static final long MASK_LEAN = (1L << (K_LEAN << 1)) - 1L;

        private final boolean leanMode;
        private final int[] denseOrders;
        private final int numDense;
        private final short[][] denseTables;
        private final int[] denseMasks;
        private final int[] denseTaus;
        private final int[] denseMaxSums;

        // Codon reading-frame experts
        private final short[] codon1Table = new short[3 * 4 * 4];
        private final short[] codon3Table = new short[3 * 64 * 4];
        private final short[] codon5Table; // null in leanMode

        // Substitutional-tolerant spaced-seed experts (null on tiny/lean inputs)
        private final short[] wobbleSeedTable; // 110_110_110_11 (8 bases = 64K contexts)
        private final short[] gappedSeedTable; // 1111_0_1111 (8 bases = 64K contexts)

        // Sparse high-order hash context experts (k=12, k=15)
        private final SparseContextTable sparse12;
        private final SparseContextTable sparse15;
        private static final long MASK_12 = (1L << 24) - 1L;
        private static final long MASK_15 = (1L << 30) - 1L;

        // Bidirectional Repeat Engine
        private final RepeatSeedIndex repFine;
        private final RepeatSeedIndex repLong;
        private final RepeatHead[] heads;
        private final int[] tmpTwoPos = new int[2];

        // Expert Indexing & State
        private final int idxCodon1;
        private final int idxCodon3;
        private final int idxCodon5;
        private final int idxWobble;
        private final int idxGapped;
        private final int idxSparse12;
        private final int idxSparse15;
        private final int idxFirstHead;
        private final int numExperts;

        private final int[][] expertP;     // [numExperts][4] in P_SCALE (16384)
        private final int[][] expertLnP;   // [numExperts][4] in Q12 (4096)
        private final int[] uniformP = {4096, 4096, 4096, 4096};
        private final int[] pOrd0 = {4096, 4096, 4096, 4096};
        private final int[] pOrd2 = {4096, 4096, 4096, 4096};

        // Sub-Mixer 1: Fast Decayed Bayesian Log-Likelihood Pool (Q12 log-scores)
        private final long[] bayesLogScore;
        private final int[] bayesP = new int[4];

        // Sub-Mixer 2: 8-Regime Online Multinomial Logistic Neural Mixer (Q12 weights, 4096 = 1.0)
        private final int[][] neuralWeights; // [8][numExperts]
        private final int[] neuralP = new int[4];
        private int activeRegime = 0;

        // Meta-Mixer blending Sub-Mixer 1 (Bayesian) and Sub-Mixer 2 (Neural Logistic)
        private long metaScoreBayes = 0L;
        private long metaScoreNeural = -2048L; // slight initial advantage to Bayesian on first ~30 bases

        private long fwdHist = 0L;
        private long rcHist = 0L;

        QcrmModel(boolean leanMode, int seqLen) {
            this.leanMode = leanMode;
            if (leanMode) {
                this.denseOrders = new int[] {0, 1, 2, 3, 4, 5, 6};
            } else if (seqLen >= 300_000) {
                this.denseOrders = new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
            } else if (seqLen >= 32_768) {
                this.denseOrders = new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
            } else {
                this.denseOrders = new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8};
            }

            this.numDense = denseOrders.length;
            this.denseTables = new short[numDense][];
            this.denseMasks = new int[numDense];
            this.denseTaus = new int[numDense];
            this.denseMaxSums = new int[numDense];

            for (int m = 0; m < numDense; m++) {
                int k = denseOrders[m];
                int numCtx = 1 << (k << 1);
                denseTables[m] = new short[numCtx << 2];
                denseMasks[m] = numCtx - 1;
                if (k == 0) {
                    denseTaus[m] = 32;
                    denseMaxSums[m] = 16384;
                } else if (k <= 2) {
                    denseTaus[m] = 24;
                    denseMaxSums[m] = 16384;
                } else if (k <= 4) {
                    denseTaus[m] = 16;
                    denseMaxSums[m] = 8192;
                } else if (k <= 6) {
                    denseTaus[m] = 8;
                    denseMaxSums[m] = 4096;
                } else if (k <= 8) {
                    denseTaus[m] = 3;
                    denseMaxSums[m] = 2048;
                } else {
                    denseTaus[m] = 1;
                    denseMaxSums[m] = 2048;
                }
            }

            int nextIdx = numDense;
            this.idxCodon1 = nextIdx++;
            this.idxCodon3 = nextIdx++;

            if (!leanMode) {
                this.codon5Table = new short[3 * 1024 * 4];
                this.idxCodon5 = nextIdx++;
            } else {
                this.codon5Table = null;
                this.idxCodon5 = -1;
            }

            if (!leanMode && seqLen >= 32_768) {
                this.wobbleSeedTable = new short[65536 * 4];
                this.gappedSeedTable = new short[65536 * 4];
                this.idxWobble = nextIdx++;
                this.idxGapped = nextIdx++;
            } else {
                this.wobbleSeedTable = null;
                this.gappedSeedTable = null;
                this.idxWobble = -1;
                this.idxGapped = -1;
            }

            if (!leanMode && seqLen >= 90_000) {
                int slots = (seqLen >= 1_000_000) ? (1 << 20) : (1 << 18);
                this.sparse12 = new SparseContextTable(slots);
                this.sparse15 = new SparseContextTable(slots);
                this.idxSparse12 = nextIdx++;
                this.idxSparse15 = nextIdx++;
            } else {
                this.sparse12 = null;
                this.sparse15 = null;
                this.idxSparse12 = -1;
                this.idxSparse15 = -1;
            }

            // Configure Bidirectional Repeat Heads
            if (leanMode) {
                // Ultra-compact 2048-bucket (~28 KB total!) repeat index safe inside -Xmx8m!
                this.repFine = new RepeatSeedIndex(1 << 11);
                this.repLong = null;
                this.heads = new RepeatHead[2];
            } else {
                int fineBuckets = (seqLen >= 500_000) ? (1 << 19) : (1 << 16);
                this.repFine = new RepeatSeedIndex(fineBuckets);
                if (seqLen >= 25_000) {
                    int longBuckets = (seqLen >= 500_000) ? (1 << 19) : (1 << 16);
                    this.repLong = new RepeatSeedIndex(longBuckets);
                } else {
                    this.repLong = null;
                }
                this.heads = new RepeatHead[6];
            }
            for (int h = 0; h < heads.length; h++) {
                heads[h] = new RepeatHead();
            }
            this.idxFirstHead = nextIdx;
            nextIdx += heads.length;

            this.numExperts = nextIdx;
            this.expertP = new int[numExperts][4];
            this.expertLnP = new int[numExperts][4];
            this.bayesLogScore = new long[numExperts];
            this.neuralWeights = new int[8][numExperts];

            // Initialize priors for Bayesian & Neural weights
            for (int m = 0; m < numDense; m++) {
                int k = denseOrders[m];
                bayesLogScore[m] = (k <= 3) ? 0L : -6144L; // -1.5 * 4096
            }
            bayesLogScore[idxCodon1] = -2048L;
            bayesLogScore[idxCodon3] = -3276L;
            if (idxCodon5 >= 0) bayesLogScore[idxCodon5] = -4500L;
            if (idxWobble >= 0) bayesLogScore[idxWobble] = -4096L;
            if (idxGapped >= 0) bayesLogScore[idxGapped] = -4096L;
            if (idxSparse12 >= 0) bayesLogScore[idxSparse12] = -6144L;
            if (idxSparse15 >= 0) bayesLogScore[idxSparse15] = -7000L;
            for (int h = 0; h < heads.length; h++) {
                bayesLogScore[idxFirstHead + h] = -4096L;
            }

            // Initialize 8-regime neural logistic weights (sum across baseline experts ~= 4096 = 1.0)
            for (int r = 0; r < 8; r++) {
                boolean repRegime = (r & 4) != 0;
                boolean hiCtxRegime = (r & 1) != 0;
                for (int m = 0; m < numDense; m++) {
                    int k = denseOrders[m];
                    if (k >= 1 && k <= 4) {
                        neuralWeights[r][m] = hiCtxRegime ? 350 : 750;
                    } else if (k >= 5 && k <= 8) {
                        neuralWeights[r][m] = hiCtxRegime ? 650 : 250;
                    } else if (k >= 9) {
                        neuralWeights[r][m] = hiCtxRegime ? 800 : 150;
                    } else {
                        neuralWeights[r][m] = 150;
                    }
                }
                neuralWeights[r][idxCodon1] = 220;
                neuralWeights[r][idxCodon3] = 320;
                if (idxCodon5 >= 0) neuralWeights[r][idxCodon5] = 260;
                if (idxWobble >= 0) neuralWeights[r][idxWobble] = 300;
                if (idxGapped >= 0) neuralWeights[r][idxGapped] = 300;
                if (idxSparse12 >= 0) neuralWeights[r][idxSparse12] = hiCtxRegime ? 650 : 150;
                if (idxSparse15 >= 0) neuralWeights[r][idxSparse15] = hiCtxRegime ? 750 : 150;
                for (int h = 0; h < heads.length; h++) {
                    neuralWeights[r][idxFirstHead + h] = repRegime ? (h == 0 || h == 2 ? 1200 : 600) : 100;
                }
            }
        }

        private void spawnRepeatCandidates(byte[] history, int pos) {
            if (leanMode) {
                if (pos >= K_LEAN) {
                    long kmer = fwdHist & MASK_LEAN;
                    repFine.lookupFwd(kmer, tmpTwoPos);
                    int p0 = tmpTwoPos[0];
                    if (p0 > 0 && p0 < pos && !isCursorTracked(p0, 1)) {
                        if (!heads[0].active) {
                            heads[0].activate(p0, 1, 32, 5);
                        }
                    }
                    int rcCur = repFine.lookupRev(kmer);
                    if (rcCur >= 0 && rcCur < pos && !isCursorTracked(rcCur, -1)) {
                        if (!heads[1].active) {
                            heads[1].activate(rcCur, -1, 32, 5);
                        }
                    }
                }
                return;
            }

            // 1. Long homolog repeat seeds (k=18): ultra-specific (4^-18), high initial weight & mismatch tolerance
            if (repLong != null && pos >= K_LONG) {
                long kmer18 = fwdHist & MASK_LONG;
                repLong.lookupFwd(kmer18, tmpTwoPos);
                for (int c = 0; c < 2; c++) {
                    int candPos = tmpTwoPos[c];
                    if (candPos > 0 && candPos < pos && !isCursorTracked(candPos, 1)
                            && verifyExactMatchFwd(history, pos, candPos, K_LONG)) {
                        allocateHeadSlot(0, 1, candPos, 1, 56, 8);
                    }
                }
                int revCur18 = repLong.lookupRev(kmer18);
                if (revCur18 >= 0 && revCur18 + K_LONG < pos && !isCursorTracked(revCur18, -1)
                        && verifyExactMatchRev(history, pos, revCur18, K_LONG)) {
                    allocateHeadSlot(1, 0, revCur18, -1, 56, 8);
                }
            }

            // 2. Fine repeat seeds (k=13): catches shorter & diverged paralogs/hairpins
            if (pos >= K_FINE) {
                long kmer13 = fwdHist & MASK_FINE;
                repFine.lookupFwd(kmer13, tmpTwoPos);
                for (int c = 0; c < 2; c++) {
                    int candPos = tmpTwoPos[c];
                    if (candPos > 0 && candPos < pos && !isCursorTracked(candPos, 1)
                            && verifyExactMatchFwd(history, pos, candPos, K_FINE)) {
                        allocateHeadSlot(2, 3, candPos, 1, 28, 6);
                    }
                }
                int revCur13 = repFine.lookupRev(kmer13);
                if (revCur13 >= 0 && revCur13 + K_FINE < pos && !isCursorTracked(revCur13, -1)
                        && verifyExactMatchRev(history, pos, revCur13, K_FINE)) {
                    if (!heads[4].active || heads[4].streak == 0 && heads[4].misses >= 2) {
                        promoteVeteranIfStrong(heads[4]);
                        heads[4].activate(revCur13, -1, 28, 6);
                    }
                }
            }
        }

        private void allocateHeadSlot(int slotA, int slotB, int startCur, int dir, int initW, int maxMiss) {
            if (!heads[slotA].active) {
                heads[slotA].activate(startCur, dir, initW, maxMiss);
            } else if (!heads[slotB].active) {
                heads[slotB].activate(startCur, dir, initW, maxMiss);
            } else if (heads[slotB].misses >= 2 && heads[slotB].streak == 0) {
                promoteVeteranIfStrong(heads[slotB]);
                heads[slotB].activate(startCur, dir, initW, maxMiss);
            }
        }

        private void promoteVeteranIfStrong(RepeatHead evicted) {
            if (heads.length >= 6 && evicted.active && evicted.totalMatched >= 20) {
                RepeatHead vet = heads[5];
                if (!vet.active || evicted.totalMatched > vet.totalMatched) {
                    vet.active = true;
                    vet.cursor = evicted.cursor;
                    vet.dir = evicted.dir;
                    vet.misses = evicted.misses;
                    vet.maxMisses = 9;
                    vet.streak = evicted.streak;
                    vet.totalMatched = evicted.totalMatched;
                    System.arraycopy(evicted.xorCnt, 0, vet.xorCnt, 0, 4);
                }
            }
        }

        private boolean isCursorTracked(int candCursor, int dir) {
            for (RepeatHead h : heads) {
                if (h.active && h.dir == dir && h.cursor == candCursor) {
                    return true;
                }
            }
            return false;
        }

        private static boolean verifyExactMatchFwd(byte[] hist, int pos, int candNextPos, int k) {
            int a = pos - k;
            int b = candNextPos - k;
            if (b < 0) return false;
            for (int i = 0; i < k; i++) {
                if (hist[a + i] != hist[b + i]) return false;
            }
            return true;
        }

        private static boolean verifyExactMatchRev(byte[] hist, int pos, int candBeforeSeed, int k) {
            int a = pos - k;
            int bTop = candBeforeSeed + k;
            if (bTop >= pos) return false;
            for (int i = 0; i < k; i++) {
                if (hist[a + i] != (hist[bTop - i] ^ 3)) return false;
            }
            return true;
        }

        void predict(byte[] history, int pos, int[] outFreqs) {
            // 1. Spawn any newly triggered Direct (+) or Inverted (-) repeat heads from history[0..pos-1]
            spawnRepeatCandidates(history, pos);

            int highOrderCountSum = 0;

            // 2. Dense Hierarchical Markov Experts (k = 0 .. maxDense)
            for (int m = 0; m < numDense; m++) {
                int k = denseOrders[m];
                int[] ep = expertP[m];
                int[] prior = (k == 0) ? uniformP : (k <= 2 ? pOrd0 : pOrd2);
                if (pos < k) {
                    System.arraycopy(prior, 0, ep, 0, 4);
                } else {
                    int ctx = (k == 0) ? 0 : (int) (fwdHist & denseMasks[m]);
                    int sum = predictFromShortTable(denseTables[m], ctx << 2, prior, denseTaus[m], ep);
                    if (k >= 7) {
                        highOrderCountSum += sum;
                    }
                }
                if (k == 0) {
                    System.arraycopy(ep, 0, pOrd0, 0, 4);
                } else if (k == 2) {
                    System.arraycopy(ep, 0, pOrd2, 0, 4);
                }
            }

            // 3. Reading-Frame 3-Periodic Codon-Phase Experts
            int phase = pos % 3;
            if (pos < 1) {
                System.arraycopy(pOrd0, 0, expertP[idxCodon1], 0, 4);
            } else {
                int ctx1 = (int) (fwdHist & 3L);
                predictFromShortTable(codon1Table, ((phase << 2) | ctx1) << 2, pOrd0, 16, expertP[idxCodon1]);
            }

            if (pos < 3) {
                System.arraycopy(pOrd2, 0, expertP[idxCodon3], 0, 4);
            } else {
                int ctx3 = (int) (fwdHist & 63L);
                predictFromShortTable(codon3Table, ((phase * 64) + ctx3) << 2, pOrd2, 12, expertP[idxCodon3]);
            }

            if (idxCodon5 >= 0) {
                if (pos < 5) {
                    System.arraycopy(pOrd2, 0, expertP[idxCodon5], 0, 4);
                } else {
                    int ctx5 = (int) (fwdHist & 1023L);
                    predictFromShortTable(codon5Table, ((phase * 1024) + ctx5) << 2, pOrd2, 8, expertP[idxCodon5]);
                }
            }

            // 4. Wobble-Skipping (110_110_110_11) & SNP-Bridging (1111_0_1111) Spaced-Seed Experts
            if (idxWobble >= 0) {
                if (pos < 11) {
                    System.arraycopy(pOrd2, 0, expertP[idxWobble], 0, 4);
                } else {
                    int wCtx = extractWobbleSeed(fwdHist);
                    predictFromShortTable(wobbleSeedTable, wCtx << 2, pOrd2, 8, expertP[idxWobble]);
                }
                if (pos < 9) {
                    System.arraycopy(pOrd2, 0, expertP[idxGapped], 0, 4);
                } else {
                    int gCtx = extractGapped9Seed(fwdHist);
                    predictFromShortTable(gappedSeedTable, gCtx << 2, pOrd2, 8, expertP[idxGapped]);
                }
            }

            // 5. High-Order Cache-Associative Sparse Hash Contexts (k=12, k=15)
            if (idxSparse12 >= 0) {
                if (pos < 12) {
                    System.arraycopy(pOrd2, 0, expertP[idxSparse12], 0, 4);
                } else {
                    highOrderCountSum += sparse12.predict(fwdHist & MASK_12, pOrd2, 2, expertP[idxSparse12]);
                }
                if (pos < 15) {
                    System.arraycopy(pOrd2, 0, expertP[idxSparse15], 0, 4);
                } else {
                    highOrderCountSum += sparse15.predict(fwdHist & MASK_15, pOrd2, 1, expertP[idxSparse15]);
                }
            }

            // 6. Bidirectional Kimura-XOR Repeat Heads
            boolean anyStrongRepeat = false;
            for (int h = 0; h < heads.length; h++) {
                boolean active = heads[h].predict(history, pos, pOrd2, expertP[idxFirstHead + h]);
                if (active && (heads[h].streak >= 2 || heads[h].totalMatched >= 6)) {
                    anyStrongRepeat = true;
                }
            }

            // Precompute integer log-probabilities LN_LUT[p] (Q12) for all experts
            for (int m = 0; m < numExperts; m++) {
                int[] ep = expertP[m];
                int[] lp = expertLnP[m];
                lp[0] = LN_LUT[ep[0]];
                lp[1] = LN_LUT[ep[1]];
                lp[2] = LN_LUT[ep[2]];
                lp[3] = LN_LUT[ep[3]];
            }

            // 7. Sub-Mixer 1: Fast Decayed Bayesian Opinion Pool (pure integer EXP_NEG_LUT)
            long maxScore = bayesLogScore[0];
            for (int m = 1; m < numExperts; m++) {
                if (bayesLogScore[m] > maxScore) maxScore = bayesLogScore[m];
            }
            long bSum0 = 0L, bSum1 = 0L, bSum2 = 0L, bSum3 = 0L;
            long wTotal = 0L;
            for (int m = 0; m < numExperts; m++) {
                long diff = maxScore - bayesLogScore[m];
                int lutIdx = (diff >= EXP_MAX_X) ? EXP_MAX_X : (int) diff;
                long wm = EXP_NEG_LUT[lutIdx];
                wTotal += wm;
                int[] ep = expertP[m];
                bSum0 += wm * ep[0];
                bSum1 += wm * ep[1];
                bSum2 += wm * ep[2];
                bSum3 += wm * ep[3];
            }
            int bp0 = (int) ((bSum0 + (wTotal >> 1)) / wTotal);
            int bp1 = (int) ((bSum1 + (wTotal >> 1)) / wTotal);
            int bp2 = (int) ((bSum2 + (wTotal >> 1)) / wTotal);
            int bp3 = (int) ((bSum3 + (wTotal >> 1)) / wTotal);
            normalizeP(bp0, bp1, bp2, bp3, bayesP);

            // 8. Sub-Mixer 2: 8-Regime Online Multinomial Logistic Neural Mixer (zero inactive dilution!)
            activeRegime = (anyStrongRepeat ? 4 : 0)
                    | ((phase == 2) ? 2 : 0)
                    | ((highOrderCountSum >= 4) ? 1 : 0);
            int[] rw = neuralWeights[activeRegime];
            long z0 = 0L, z1 = 0L, z2 = 0L, z3 = 0L;
            for (int m = 0; m < numExperts; m++) {
                long wm = rw[m]; // Q12
                int[] lp = expertLnP[m]; // Q12
                z0 += wm * lp[0];
                z1 += wm * lp[1];
                z2 += wm * lp[2];
                z3 += wm * lp[3];
            }
            // Scale Q24 logits back to Q12
            z0 >>= 12;
            z1 >>= 12;
            z2 >>= 12;
            z3 >>= 12;
            long zMax = Math.max(Math.max(z0, z1), Math.max(z2, z3));
            long e0 = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, zMax - z0)];
            long e1 = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, zMax - z1)];
            long e2 = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, zMax - z2)];
            long e3 = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, zMax - z3)];
            long eSum = e0 + e1 + e2 + e3;
            int np0 = (int) ((e0 * P_SCALE + (eSum >> 1)) / eSum);
            int np1 = (int) ((e1 * P_SCALE + (eSum >> 1)) / eSum);
            int np2 = (int) ((e2 * P_SCALE + (eSum >> 1)) / eSum);
            int np3 = (int) ((e3 * P_SCALE + (eSum >> 1)) / eSum);
            normalizeP(np0, np1, np2, np3, neuralP);

            // 9. Exact Meta-Combiner between Bayesian Pool and Regime Neural Mixer -> TOT (32768)
            long mMax = Math.max(metaScoreBayes, metaScoreNeural);
            long wb = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, mMax - metaScoreBayes)];
            long wn = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, mMax - metaScoreNeural)];
            long wMeta = wb + wn;

            // Scale from P_SCALE (16384) to TOT (32768) = factor of 2
            int f0 = Math.max(1, (int) (((wb * bayesP[0] + wn * neuralP[0]) * 2L + (wMeta >> 1)) / wMeta));
            int f1 = Math.max(1, (int) (((wb * bayesP[1] + wn * neuralP[1]) * 2L + (wMeta >> 1)) / wMeta));
            int f2 = Math.max(1, (int) (((wb * bayesP[2] + wn * neuralP[2]) * 2L + (wMeta >> 1)) / wMeta));
            int f3 = Math.max(1, (int) (((wb * bayesP[3] + wn * neuralP[3]) * 2L + (wMeta >> 1)) / wMeta));
            int diff = TOT - (f0 + f1 + f2 + f3);
            if (diff != 0) {
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

        void update(byte[] history, int pos, int sym) {
            // 1. Update Meta-Mixer scores (decay = 255/256 ~= 0.9961, learning rate = 1.5)
            metaScoreBayes = ((metaScoreBayes * 255L) >> 8) + ((3L * LN_LUT[bayesP[sym]]) >> 1);
            metaScoreNeural = ((metaScoreNeural * 255L) >> 8) + ((3L * LN_LUT[neuralP[sym]]) >> 1);

            // 2. Update Bayesian Pool expert log-scores (decay = 255/256, eta = 1.625 = 13/8)
            for (int m = 0; m < numExperts; m++) {
                bayesLogScore[m] = ((bayesLogScore[m] * 255L) >> 8) + ((13L * expertLnP[m][sym]) >> 3);
            }

            // 3. Online Multinomial Logistic Neural SGD step on activeRegime
            // Gradient of -ln P_neural(sym) w.r.t. w_m is sum_b (P_neural(b) - I[b==sym]) * ln P_m(b)
            int r0 = (sym == 0 ? P_SCALE : 0) - neuralP[0]; // Q14
            int r1 = (sym == 1 ? P_SCALE : 0) - neuralP[1];
            int r2 = (sym == 2 ? P_SCALE : 0) - neuralP[2];
            int r3 = (sym == 3 ? P_SCALE : 0) - neuralP[3];
            int[] rw = neuralWeights[activeRegime];
            // Step size eta ~= 0.035: Q14 * Q12 = Q26; >> 19 gives Q7 ~= eta * Q12!
            for (int m = 0; m < numExperts; m++) {
                int[] lp = expertLnP[m];
                long grad = (long) r0 * lp[0] + (long) r1 * lp[1] + (long) r2 * lp[2] + (long) r3 * lp[3];
                int delta = (int) (grad >> 20);
                // Mild L2 shrinkage (1/2048) + non-negative weight clamp [0 .. 8192]
                int updated = rw[m] - (rw[m] >> 11) + delta;
                if (updated < 0) updated = 0;
                else if (updated > 8192) updated = 8192;
                rw[m] = updated;
            }

            // 4. Advance active Repeat Heads BEFORE inserting current k-mer into RepeatSeedIndex
            for (RepeatHead h : heads) {
                h.step(history, pos, sym);
            }

            // 5. Update Markov, Codon, Spaced-Seed, Sparse-Hash, and Repeat Seed Indices
            long prevFwdHist = fwdHist;
            rcHist = (rcHist >>> 2) | ((long) (sym ^ 3) << 60);

            // Insert k-mer ending at pos - 1 pointing to pos (Direct +), and RC k-mer ending at pos pointing to pos - k (Inverted -)
            if (leanMode) {
                if (pos >= K_LEAN) {
                    repFine.insertFwd(prevFwdHist & MASK_LEAN, pos);
                }
                if (pos >= K_LEAN) {
                    long rcKmer = (rcHist >>> (62 - (K_LEAN << 1))) & MASK_LEAN;
                    repFine.insertRev(rcKmer, pos - K_LEAN);
                }
            } else {
                if (pos >= K_FINE) {
                    repFine.insertFwd(prevFwdHist & MASK_FINE, pos);
                    long rc13 = (rcHist >>> (62 - (K_FINE << 1))) & MASK_FINE;
                    repFine.insertRev(rc13, pos - K_FINE);
                }
                if (repLong != null && pos >= K_LONG) {
                    repLong.insertFwd(prevFwdHist & MASK_LONG, pos);
                    long rc18 = (rcHist >>> (62 - (K_LONG << 1))) & MASK_LONG;
                    repLong.insertRev(rc18, pos - K_LONG);
                }
            }

            // Dense Markov updates (Forward +2, Chargaff RC +1 for k >= 4)
            for (int m = 0; m < numDense; m++) {
                int k = denseOrders[m];
                if (pos >= k) {
                    short[] tab = denseTables[m];
                    int ctx = (k == 0) ? 0 : (int) (prevFwdHist & denseMasks[m]);
                    incrementBucket(tab, ctx << 2, sym, 2, denseMaxSums[m]);

                    if (k >= 4) {
                        int oldestBase = (int) ((prevFwdHist >>> ((k - 1) << 1)) & 3L);
                        int rcSym = oldestBase ^ 3;
                        int rcCtx = (int) ((rcHist >>> (62 - (k << 1))) & denseMasks[m]);
                        incrementBucket(tab, rcCtx << 2, rcSym, 1, denseMaxSums[m]);
                    }
                }
            }

            // Codon-phase updates
            int phase = pos % 3;
            if (pos >= 1) {
                int ctx1 = (int) (prevFwdHist & 3L);
                incrementBucket(codon1Table, ((phase << 2) | ctx1) << 2, sym, 2, 160);
            }
            if (pos >= 3) {
                int ctx3 = (int) (prevFwdHist & 63L);
                incrementBucket(codon3Table, ((phase * 64) + ctx3) << 2, sym, 2, 320);
            }
            if (idxCodon5 >= 0 && pos >= 5) {
                int ctx5 = (int) (prevFwdHist & 1023L);
                incrementBucket(codon5Table, ((phase * 1024) + ctx5) << 2, sym, 2, 512);
            }

            // Spaced-seed updates
            if (idxWobble >= 0) {
                if (pos >= 11) {
                    int wCtx = extractWobbleSeed(prevFwdHist);
                    incrementBucket(wobbleSeedTable, wCtx << 2, sym, 2, 1024);
                }
                if (pos >= 9) {
                    int gCtx = extractGapped9Seed(prevFwdHist);
                    incrementBucket(gappedSeedTable, gCtx << 2, sym, 2, 1024);
                    // Symmetric palindrome mask 1111_0_1111 also supports Chargaff RC update!
                    int oldest9 = (int) ((prevFwdHist >>> 16) & 3L);
                    long rcWindow9 = rcHist >>> (62 - 18);
                    int rcGCtx = extractGapped9Seed(rcWindow9);
                    incrementBucket(gappedSeedTable, rcGCtx << 2, oldest9 ^ 3, 1, 1024);
                }
            }

            // Sparse high-order hash updates (Forward +2, Chargaff RC +1)
            if (idxSparse12 >= 0) {
                if (pos >= 12) {
                    sparse12.update(prevFwdHist & MASK_12, sym, 2);
                    int old12 = (int) ((prevFwdHist >>> 22) & 3L);
                    long rcCtx12 = (rcHist >>> (62 - 24)) & MASK_12;
                    sparse12.update(rcCtx12, old12 ^ 3, 1);
                }
                if (pos >= 15) {
                    sparse15.update(prevFwdHist & MASK_15, sym, 2);
                    int old15 = (int) ((prevFwdHist >>> 28) & 3L);
                    long rcCtx15 = (rcHist >>> (62 - 30)) & MASK_15;
                    sparse15.update(rcCtx15, old15 ^ 3, 1);
                }
            }

            fwdHist = (prevFwdHist << 2) | sym;
        }

        /** Extracts 8 bases (16 bits) from 11-bp history using codon-skeleton mask 110_110_110_11. */
        private static int extractWobbleSeed(long h) {
            int b01 = (int) (h & 0xFL);               // offsets 0, 1 (4 bits)
            int b34 = (int) ((h >>> 6) & 0xFL);       // offsets 3, 4 (4 bits)
            int b67 = (int) ((h >>> 12) & 0xFL);      // offsets 6, 7 (4 bits)
            int b910 = (int) ((h >>> 18) & 0xFL);     // offsets 9, 10 (4 bits)
            return b01 | (b34 << 4) | (b67 << 8) | (b910 << 12);
        }

        /** Extracts 8 bases (16 bits) from 9-bp history skipping middle SNP base (1111_0_1111). */
        private static int extractGapped9Seed(long h) {
            int lo4 = (int) (h & 0xFFL);              // offsets 0..3 (8 bits)
            int hi4 = (int) ((h >>> 10) & 0xFFL);     // offsets 5..8 (8 bits, skipping offset 4)
            return lo4 | (hi4 << 8);
        }

        private static int predictFromShortTable(short[] tab, int baseIdx, int[] priorP, int tau, int[] outP) {
            int c0 = tab[baseIdx];
            int c1 = tab[baseIdx + 1];
            int c2 = tab[baseIdx + 2];
            int c3 = tab[baseIdx + 3];
            int sum = c0 + c1 + c2 + c3;
            if (sum == 0) {
                System.arraycopy(priorP, 0, outP, 0, 4);
                return 0;
            }
            long denom = 4L * sum + tau; // tau in quarter-count units
            int p0 = Math.max(1, (int) (((long) c0 * 4L * P_SCALE + (long) tau * priorP[0] + (denom >> 1)) / denom));
            int p1 = Math.max(1, (int) (((long) c1 * 4L * P_SCALE + (long) tau * priorP[1] + (denom >> 1)) / denom));
            int p2 = Math.max(1, (int) (((long) c2 * 4L * P_SCALE + (long) tau * priorP[2] + (denom >> 1)) / denom));
            int p3 = Math.max(1, P_SCALE - (p0 + p1 + p2));
            normalizeP(p0, p1, p2, p3, outP);
            return sum;
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

    /** 64-bit carry-buffered Range Encoder (15-bit probability resolution, TOT = 32768). */
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

    /** Matching 64-bit carry-buffered Range Decoder. */
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
