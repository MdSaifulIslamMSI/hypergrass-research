package hypergrass;

import java.util.Arrays;

/**
 * High-throughput 2-bit packed DNA/RNA sequence store ({@code A=0, C=1, G=2, T/U=3}) with
 * branchless 64-bit SWAR (32-bases-per-register) forward and reverse-complement word extraction,
 * Klein four-group ({@code V_4}) XOR comparison, and strand-symmetric canonical open-syncmer
 * extraction.
 */
public final class Packed2BitSequence {

    /** Primary canonical open-syncmer parameters: k=16, s=10, window w=7, center offset t=3. */
    public static final int K_PRIMARY = 16;
    public static final int S_PRIMARY = 10;
    public static final int W_PRIMARY = K_PRIMARY - S_PRIMARY + 1; // 7
    public static final int T_PRIMARY = 3;

    /** Secondary diverged-island syncmer parameters: k=13, s=9, window w=5 (odd!), center t=2. */
    public static final int K_SECONDARY = 13;
    public static final int S_SECONDARY = 9;
    public static final int W_SECONDARY = K_SECONDARY - S_SECONDARY + 1; // 5
    public static final int T_SECONDARY = 2;

    private static final byte[] ASCII_TO_2BIT = new byte[256];
    private static final byte[] TWO_BIT_TO_ACGT = {'A', 'C', 'G', 'T'};
    private static final byte[] TWO_BIT_TO_ACGU = {'A', 'C', 'G', 'U'};

    static {
        Arrays.fill(ASCII_TO_2BIT, (byte) -1);
        ASCII_TO_2BIT['A'] = 0;
        ASCII_TO_2BIT['a'] = 0;
        ASCII_TO_2BIT['C'] = 1;
        ASCII_TO_2BIT['c'] = 1;
        ASCII_TO_2BIT['G'] = 2;
        ASCII_TO_2BIT['g'] = 2;
        ASCII_TO_2BIT['T'] = 3;
        ASCII_TO_2BIT['t'] = 3;
        ASCII_TO_2BIT['U'] = 3;
        ASCII_TO_2BIT['u'] = 3;
    }

    /**
     * Little-endian 2-bit lanes: base {@code i} resides in {@code words[i >>> 5]} at bit offset
     * {@code (i & 31) << 1}. Padded with at least two trailing zero words for branchless unaligned
     * 32-base loads.
     */
    public final long[] words;
    public final int length;

    public Packed2BitSequence(long[] words, int length) {
        this.words = words;
        this.length = length;
    }

    /** Packs an array of already-factorized 2-bit values ({@code 0..3} per byte). */
    public static Packed2BitSequence fromTwoBitBytes(byte[] bases, int len) {
        int numWords = ((len + 31) >>> 5) + 2;
        long[] words = new long[numWords];
        int fullWords = len >>> 5;
        int idx = 0;
        for (int w = 0; w < fullWords; w++) {
            long acc = 0L;
            for (int j = 0; j < 32; j++) {
                acc |= ((long) (bases[idx++] & 3)) << (j << 1);
            }
            words[w] = acc;
        }
        if (idx < len) {
            long acc = 0L;
            int shift = 0;
            while (idx < len) {
                acc |= ((long) (bases[idx++] & 3)) << shift;
                shift += 2;
            }
            words[fullWords] = acc;
        }
        return new Packed2BitSequence(words, len);
    }

    /**
     * Checks whether {@code raw[0..len)} consists purely of uppercase ASCII {@code A,C,G,T} with
     * no headers, whitespace, ambiguity codes, or lowercase characters, and packs it in one pass if
     * so; returns {@code null} otherwise so the caller can run full 6-channel factorization.
     */
    public static Packed2BitSequence tryFromFlatUppercaseAcgt(byte[] raw, int len) {
        int numWords = ((len + 31) >>> 5) + 2;
        long[] words = new long[numWords];
        int fullWords = len >>> 5;
        int idx = 0;
        for (int w = 0; w < fullWords; w++) {
            long acc = 0L;
            for (int j = 0; j < 32; j++) {
                int b = raw[idx++] & 0xFF;
                int v;
                if (b == 'A') v = 0;
                else if (b == 'C') v = 1;
                else if (b == 'G') v = 2;
                else if (b == 'T') v = 3;
                else return null;
                acc |= ((long) v) << (j << 1);
            }
            words[w] = acc;
        }
        if (idx < len) {
            long acc = 0L;
            int shift = 0;
            while (idx < len) {
                int b = raw[idx++] & 0xFF;
                int v;
                if (b == 'A') v = 0;
                else if (b == 'C') v = 1;
                else if (b == 'G') v = 2;
                else if (b == 'T') v = 3;
                else return null;
                acc |= ((long) v) << shift;
                shift += 2;
            }
            words[fullWords] = acc;
        }
        return new Packed2BitSequence(words, len);
    }

    /** Concatenates multiple {@link Packed2BitSequence} instances into one contiguous bank. */
    public static Packed2BitSequence concatenate(Packed2BitSequence[] seqs, int[] outOffsets) {
        return concatenateInternal(seqs, outOffsets, false);
    }

    /**
     * Concatenates {@code seqs} into one contiguous 2-bit bank while clearing {@code seqs[i] = null}
     * as each unit is copied so peak heap never duplicates the 2-bit bank under {@code -Xmx192m}.
     */
    public static Packed2BitSequence concatenateAndConsume(
            Packed2BitSequence[] seqs, int[] outOffsets) {
        return concatenateInternal(seqs, outOffsets, true);
    }

    private static Packed2BitSequence concatenateInternal(
            Packed2BitSequence[] seqs, int[] outOffsets, boolean consume) {
        long total = 0L;
        for (int i = 0; i < seqs.length; i++) {
            outOffsets[i] = (int) total;
            total += seqs[i].length;
            if (total > Integer.MAX_VALUE - 64) {
                throw new IllegalArgumentException("Concatenated sequence bank exceeds 2.14 Gbp");
            }
        }
        outOffsets[seqs.length] = (int) total;
        int totalLen = (int) total;
        long[] words = new long[((totalLen + 31) >>> 5) + 2];
        int dstBitPos = 0;
        for (int i = 0; i < seqs.length; i++) {
            Packed2BitSequence s = seqs[i];
            int rem = s.length;
            int srcPos = 0;
            while (rem > 0) {
                int take = Math.min(32, rem);
                long chunk = s.word32Forward(srcPos);
                if (take < 32) {
                    chunk &= (1L << (take << 1)) - 1L;
                }
                int wIdx = dstBitPos >>> 6;
                int bOff = dstBitPos & 63;
                words[wIdx] |= chunk << bOff;
                if (bOff + (take << 1) > 64) {
                    words[wIdx + 1] |= chunk >>> (64 - bOff);
                }
                dstBitPos += take << 1;
                srcPos += take;
                rem -= take;
            }
            if (consume) {
                seqs[i] = null;
            }
        }
        return new Packed2BitSequence(words, totalLen);
    }

    /** Extracts a 2-bit {@link Packed2BitSequence} sub-slice without allocating intermediate bytes. */
    public Packed2BitSequence sliceSeq(int startBase, int sliceLen) {
        long[] subWords = new long[((sliceLen + 31) >>> 5) + 2];
        int fullWords = sliceLen >>> 5;
        for (int w = 0; w < fullWords; w++) {
            subWords[w] = word32Forward(startBase + (w << 5));
        }
        int rem = sliceLen & 31;
        if (rem > 0) {
            subWords[fullWords] =
                    word32Forward(startBase + (fullWords << 5)) & ((1L << (rem << 1)) - 1L);
        }
        return new Packed2BitSequence(subWords, sliceLen);
    }

    /** Returns the single 2-bit base ({@code 0..3}) at {@code pos}. */
    public int baseAt(int pos) {
        return (int) ((words[pos >>> 5] >>> ((pos & 31) << 1)) & 3L);
    }

    /**
     * Returns base at coordinate {@code strandPos} within a sub-sequence {@code [seqOffset ..
     * seqOffset + seqLen)} on either forward ({@code strand == 0}) or reverse-complement ({@code
     * strand == 1}) orientation.
     */
    public int baseAtStrand(int seqOffset, int seqLen, int strand, int strandPos) {
        if (strand == 0) {
            return baseAt(seqOffset + strandPos);
        }
        int fwdPos = seqOffset + seqLen - 1 - strandPos;
        return baseAt(fwdPos) ^ 3;
    }

    /**
     * Branchless 64-bit load of up to 32 consecutive forward bases starting at {@code pos}. Caller
     * is responsible for clamping {@code min(32, length - pos)} valid lanes.
     */
    public long word32Forward(int pos) {
        int wIdx = pos >>> 5;
        int rem = pos & 31;
        if (rem == 0) {
            return words[wIdx];
        }
        int shift = rem << 1;
        return (words[wIdx] >>> shift) | (words[wIdx + 1] << (64 - shift));
    }

    /**
     * Reverses and Watson-Crick complements all 32 packed 2-bit bases in {@code u}: lane {@code j}
     * of output equals {@code (~lane_{31-j}(u)) & 3}. Executes in 5 register instructions with zero
     * memory loads.
     */
    public static long reverseComplementWord32(long u) {
        u = ~u;
        u = Long.reverseBytes(u);
        u = ((u & 0x0F0F0F0F0F0F0F0FL) << 4) | ((u >>> 4) & 0x0F0F0F0F0F0F0F0FL);
        u = ((u & 0x3333333333333333L) << 2) | ((u >>> 2) & 0x3333333333333333L);
        return u;
    }

    /**
     * Returns up to 32 consecutive bases on {@code strand} ({@code 0} = forward, {@code 1} =
     * reverse-complement) starting at {@code strandPos} within sub-sequence {@code [seqOffset ..
     * seqOffset + seqLen)}. Valid when {@code strandPos + 32 <= seqLen}; for tail positions where
     * {@code seqLen - strandPos < 32}, falls back cleanly to exact tail packing.
     */
    public long word32Strand(int seqOffset, int seqLen, int strand, int strandPos) {
        if (strand == 0) {
            return word32Forward(seqOffset + strandPos);
        }
        int rem = seqLen - strandPos;
        if (rem >= 32) {
            int fwdBot = seqOffset + seqLen - 32 - strandPos;
            return reverseComplementWord32(word32Forward(fwdBot));
        }
        long acc = 0L;
        int fwdTop = seqOffset + seqLen - 1 - strandPos;
        for (int j = 0; j < rem; j++) {
            long b = baseAt(fwdTop - j) ^ 3;
            acc |= b << (j << 1);
        }
        return acc;
    }

    /**
     * Counts consecutive exact matching bases between {@code ref} (within {@code [refOffset ..
     * refOffset + refLen)} on {@code strand} starting at {@code refStrandPos}) and {@code tgt}
     * (forward starting at {@code tgtPos}), up to {@code maxLen}, stepping 32 bases per CPU word.
     */
    public static int countExactMatch(
            Packed2BitSequence ref,
            int refOffset,
            int refLen,
            int strand,
            int refStrandPos,
            Packed2BitSequence tgt,
            int tgtPos,
            int maxLen) {
        int limit = Math.min(maxLen, Math.min(refLen - refStrandPos, tgt.length - tgtPos));
        int k = 0;
        if (strand == 0) {
            int rBase = refOffset + refStrandPos;
            while (k + 32 <= limit) {
                long diff = ref.word32Forward(rBase + k) ^ tgt.word32Forward(tgtPos + k);
                if (diff != 0L) {
                    return k + (Long.numberOfTrailingZeros(diff) >>> 1);
                }
                k += 32;
            }
            if (k < limit) {
                long diff = ref.word32Forward(rBase + k) ^ tgt.word32Forward(tgtPos + k);
                int extra = diff == 0L ? 32 : (Long.numberOfTrailingZeros(diff) >>> 1);
                return Math.min(limit, k + extra);
            }
            return limit;
        } else {
            int fwdTop = refOffset + refLen - 1 - refStrandPos;
            while (k + 32 <= limit) {
                long rWord = reverseComplementWord32(ref.word32Forward(fwdTop - k - 31));
                long diff = rWord ^ tgt.word32Forward(tgtPos + k);
                if (diff != 0L) {
                    return k + (Long.numberOfTrailingZeros(diff) >>> 1);
                }
                k += 32;
            }
            while (k < limit) {
                int rb = ref.baseAt(fwdTop - k) ^ 3;
                int tb = tgt.baseAt(tgtPos + k);
                if (rb != tb) {
                    break;
                }
                k++;
            }
            return k;
        }
    }

    /** Unpacks {@code [0..length)} into 2-bit bytes ({@code 0..3}). */
    public byte[] toTwoBitBytes() {
        return extractTwoBitSlice(0, length);
    }

    /** Writes {@code src[0..len)} (2-bit values {@code 0..3}) into this bank starting at {@code dstPos}. */
    public void writeTwoBitSlice(int dstPos, byte[] src, int len) {
        int i = 0;
        while (i < len && ((dstPos + i) & 31) != 0) {
            int p = dstPos + i;
            int wIdx = p >>> 5;
            int shift = (p & 31) << 1;
            words[wIdx] = (words[wIdx] & ~(3L << shift)) | (((long) (src[i] & 3)) << shift);
            i++;
        }
        int wIdx = (dstPos + i) >>> 5;
        while (i + 32 <= len) {
            long w = 0L;
            for (int j = 0; j < 32; j++) {
                w |= ((long) (src[i + j] & 3)) << (j << 1);
            }
            words[wIdx++] = w;
            i += 32;
        }
        while (i < len) {
            int p = dstPos + i;
            int remW = p >>> 5;
            int shift = (p & 31) << 1;
            words[remW] = (words[remW] & ~(3L << shift)) | (((long) (src[i] & 3)) << shift);
            i++;
        }
    }

    /** Extracts {@code [srcPos .. srcPos + len)} into a 2-bit byte array ({@code 0..3}). */
    public byte[] extractTwoBitSlice(int srcPos, int len) {
        byte[] out = new byte[len];
        int k = 0;
        while (k + 32 <= len) {
            long w = word32Forward(srcPos + k);
            for (int j = 0; j < 32; j++) {
                out[k + j] = (byte) (w & 3L);
                w >>>= 2;
            }
            k += 32;
        }
        while (k < len) {
            out[k] = (byte) baseAt(srcPos + k);
            k++;
        }
        return out;
    }

    /** Unpacks directly into uppercase ASCII {@code A,C,G,T} (or {@code A,C,G,U} if {@code rna}). */
    public byte[] toUppercaseAscii(boolean rna) {
        byte[] map = rna ? TWO_BIT_TO_ACGU : TWO_BIT_TO_ACGT;
        byte[] out = new byte[length];
        int fullWords = length >>> 5;
        int idx = 0;
        for (int w = 0; w < fullWords; w++) {
            long word = words[w];
            for (int j = 0; j < 32; j++) {
                out[idx++] = map[(int) (word & 3L)];
                word >>>= 2;
            }
        }
        if (idx < length) {
            long word = words[fullWords];
            while (idx < length) {
                out[idx++] = map[(int) (word & 3L)];
                word >>>= 2;
            }
        }
        return out;
    }

    /** Fast avalanche mixer for canonical k-mer / submer integer keys. */
    public static int mix32(long z) {
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        z = z ^ (z >>> 31);
        return (int) z;
    }

    public static long mix64(long z) {
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    /** Visitor callback for canonical open syncmers. */
    @FunctionalInterface
    public interface SyncmerConsumer {
        /**
         * @param pos 0-based start offset of the k-mer within the scanned interval
         * @param canonicalKmer unsigned canonical 2k-bit k-mer value {@code min(fwd, rc)}
         * @param canonicalStrand {@code 0} if {@code fwd <= rc}, {@code 1} if {@code rc < fwd}
         */
        void accept(int pos, long canonicalKmer, int canonicalStrand);
    }

    /**
     * Scans {@code [startOffset .. startOffset + segLen)} in one linear O(n) pass and emits all
     * strand-symmetric canonical open syncmers with parameters {@code (k, s, t)}, where {@code w =
     * k - s + 1} and {@code t = (w - 1) / 2} (strictly centered so reverse-complementation selects
     * the exact same genomic interval!).
     *
     * <p><b>Exact Strand-Symmetry Invariant:</b> When {@code w = k - s + 1} is odd and {@code t =
     * (w - 1) / 2}, reversing the order of the {@code w} canonical {@code s}-mer hashes maps index
     * {@code j} to {@code w - 1 - j}, fixing {@code t = w - 1 - t}. Ties between first and last
     * minimums are broken symmetrically by requiring strict minimum over {@code j != t}.
     */
    public void forEachCanonicalOpenSyncmer(
            int startOffset, int segLen, int k, int s, int t, SyncmerConsumer consumer) {
        if (segLen < k) {
            return;
        }
        final int w = k - s + 1;
        final long kMask = (k == 32) ? -1L : ((1L << (k << 1)) - 1L);
        final long sMask = (1L << (s << 1)) - 1L;
        final int rcKShift = (k - 1) << 1;
        final int rcSShift = (s - 1) << 1;

        int[] winHashes = new int[w]; // ring buffer of canonical s-mer hashes
        long fwdK = 0L;
        long rcK = 0L;
        long fwdS = 0L;
        long rcS = 0L;

        int end = startOffset + segLen;
        int relPos = 0;
        for (int p = startOffset; p < end; p++, relPos++) {
            long b = baseAt(p);
            long cb = b ^ 3L;

            fwdK = ((fwdK << 2) | b) & kMask;
            rcK = (rcK >>> 2) | (cb << rcKShift);

            fwdS = ((fwdS << 2) | b) & sMask;
            rcS = (rcS >>> 2) | (cb << rcSShift);

            if (relPos + 1 >= s) {
                int sIdx = relPos - s + 1;
                long canonS = Long.compareUnsigned(fwdS, rcS) <= 0 ? fwdS : rcS;
                winHashes[sIdx % w] = mix32(canonS);
            }

            if (relPos + 1 >= k) {
                int kStart = relPos - k + 1;
                int centerHash = winHashes[(kStart + t) % w];
                boolean isCenterUniqueMin = true;
                for (int j = 0; j < w; j++) {
                    if (j == t) continue;
                    if (Integer.compareUnsigned(winHashes[(kStart + j) % w], centerHash) <= 0) {
                        isCenterUniqueMin = false;
                        break;
                    }
                }
                if (isCenterUniqueMin) {
                    int cmp = Long.compareUnsigned(fwdK, rcK);
                    long canonK = cmp <= 0 ? fwdK : rcK;
                    int strand = cmp <= 0 ? 0 : 1;
                    consumer.accept(kStart, canonK, strand);
                }
            }
        }
    }

    /**
     * Computes a deterministic sorted bottom-{@code sketchSize} canonical syncmer hash sketch for
     * fast O(S) phylogenetic Jaccard estimation across genomes.
     */
    public long[] computeBottomSketch(int startOffset, int segLen, int sketchSize) {
        long[] heap = new long[sketchSize];
        int[] size = {0};
        // Max-heap on unsigned comparison via signed flip (x ^ Long.MIN_VALUE)
        forEachCanonicalOpenSyncmer(
                startOffset,
                segLen,
                K_PRIMARY,
                S_PRIMARY,
                T_PRIMARY,
                (pos, canonK, strand) -> {
                    long h = mix64(canonK) ^ Long.MIN_VALUE;
                    int n = size[0];
                    if (n < sketchSize) {
                        heap[n] = h;
                        siftUpMax(heap, n);
                        size[0] = n + 1;
                    } else if (h < heap[0]) {
                        heap[0] = h;
                        siftDownMax(heap, 0, n);
                    }
                });
        int n = size[0];
        long[] res = Arrays.copyOf(heap, n);
        Arrays.sort(res);
        // Deduplicate sorted sketch
        int u = 0;
        for (int i = 0; i < n; i++) {
            if (i == 0 || res[i] != res[i - 1]) {
                res[u++] = res[i];
            }
        }
        return u == res.length ? res : Arrays.copyOf(res, u);
    }

    /** Estimates canonical syncmer Jaccard similarity in [0, 1] between two sorted sketches. */
    public static double sketchJaccard(long[] a, long[] b) {
        if (a.length == 0 || b.length == 0) {
            return 0.0;
        }
        int i = 0;
        int j = 0;
        int inter = 0;
        int union = 0;
        int limit = Math.min(a.length, b.length);
        while (i < a.length && j < b.length && union < limit) {
            long va = a[i];
            long vb = b[j];
            if (va == vb) {
                inter++;
                i++;
                j++;
            } else if (va < vb) {
                i++;
            } else {
                j++;
            }
            union++;
        }
        return union == 0 ? 0.0 : (double) inter / union;
    }

    private static void siftUpMax(long[] heap, int idx) {
        long val = heap[idx];
        while (idx > 0) {
            int p = (idx - 1) >>> 1;
            if (heap[p] >= val) break;
            heap[idx] = heap[p];
            idx = p;
        }
        heap[idx] = val;
    }

    private static void siftDownMax(long[] heap, int idx, int n) {
        long val = heap[idx];
        int half = n >>> 1;
        while (idx < half) {
            int left = (idx << 1) + 1;
            int right = left + 1;
            int best = (right < n && heap[right] > heap[left]) ? right : left;
            if (val >= heap[best]) break;
            heap[idx] = heap[best];
            idx = best;
        }
        heap[idx] = val;
    }
}
