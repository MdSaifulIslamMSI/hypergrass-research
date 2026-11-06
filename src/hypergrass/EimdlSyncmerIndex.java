package hypergrass;

import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

/**
 * Memory-frugal Strand-Symmetric Canonical Open-Syncmer CSR Index (~1.1 B/base including 2-bit DNA)
 * plus a cache-resident 4-way Bidirectional Self-Repeat Syncmer Table.
 *
 * <p>Key guarantees:
 * <ol>
 *   <li><b>Exact 32-Bit Signature Pre-Filter:</b> Combines a $2^B$-bucket CSR directory with a
 *       16-bit signature {@code sig16} per entry so all 32 bits of every canonical 16-mer are
 *       verified sequentially in L1/L2 cache before touching sequence words.</li>
 *   <li><b>Zero Strand Duplication:</b> Stores one canonical entry per open syncmer with bit 31
 *       recording canonical orientation ({@code 0} or {@code 1}), halving RAM vs. dual-strand
 *       k-mer tables.</li>
 *   <li><b>Within-Donor High-Copy Cap + Distinct-Donor Preservation:</b> Caps pathological
 *       within-genome satellite floods while preserving distinct donor genomes so twin strains are
 *       never evicted (preserving Conditional Monotonicity, Law III).</li>
 * </ol>
 */
public final class EimdlSyncmerIndex {

    private static final int BLOCK_SHIFT = 14; // 16,384-base blocks for O(1) globalPos -> donorId
    private static final int MAX_WITHIN_DONOR_BUCKET_OCC = 16;

    public final Packed2BitSequence bank;
    public final int numDonors;
    public final int[] donorOffsets; // length numDonors + 1

    private final int[] posToDonorFloor;

    // Primary k=16 canonical open-syncmer CSR
    public final int primarySubsampleShift;
    public final int primarySubsampleMask;
    private final int bucketBits;
    private final int bucketMask;
    private final int[] bucketOffsets;
    private final int[] posAndStrand; // bit 31 = canonicalStrand, bits 0..30 = globalBasePos
    private final short[] sig16;

    // Optional secondary k=13 syncmer CSR for diverged homologous islands (null when disabled)
    private final int secBucketMask;
    private final int[] secBucketOffsets;
    private final int[] secPosAndStrand;
    private final short[] secSig16;

    private EimdlSyncmerIndex(
            Packed2BitSequence bank,
            int numDonors,
            int[] donorOffsets,
            int[] posToDonorFloor,
            int primarySubsampleShift,
            int bucketBits,
            int[] bucketOffsets,
            int[] posAndStrand,
            short[] sig16,
            int secBucketMask,
            int[] secBucketOffsets,
            int[] secPosAndStrand,
            short[] secSig16) {
        this.bank = bank;
        this.numDonors = numDonors;
        this.donorOffsets = donorOffsets;
        this.posToDonorFloor = posToDonorFloor;
        this.primarySubsampleShift = primarySubsampleShift;
        this.primarySubsampleMask = (1 << primarySubsampleShift) - 1;
        this.bucketBits = bucketBits;
        this.bucketMask = (1 << bucketBits) - 1;
        this.bucketOffsets = bucketOffsets;
        this.posAndStrand = posAndStrand;
        this.sig16 = sig16;
        this.secBucketMask = secBucketMask;
        this.secBucketOffsets = secBucketOffsets;
        this.secPosAndStrand = secPosAndStrand;
        this.secSig16 = secSig16;
    }

    /**
     * Builds a CSR canonical open-syncmer index over {@code donors[0..numDonors)}. If {@code
     * indexedDonorMask != null}, only donors with {@code indexedDonorMask[d] == true} are inserted
     * into the shared CSR tables (allowing selective donor indexing while keeping all donors
     * addressable in {@link #bank}).
     */
    public static EimdlSyncmerIndex build(
            Packed2BitSequence[] donors,
            boolean[] indexedDonorMask,
            boolean buildSecondary13,
            int threads) {
        return buildInternal(donors, indexedDonorMask, 0, false, buildSecondary13, threads);
    }

    /**
     * Builds a zero-materialization CSR canonical open-syncmer index across <b>all</b> donors using
     * deterministic hash subsampling ({@code 1 / 2^primarySubsampleShift}) and optional incremental
     * donor consumption ({@code consumeDonors == true}) so every donor remains searchable even
     * inside {@code -Xmx192m}.
     */
    public static EimdlSyncmerIndex buildWithSubsampling(
            Packed2BitSequence[] donors,
            int primarySubsampleShift,
            boolean consumeDonors,
            boolean buildSecondary13,
            int threads) {
        return buildInternal(
                donors, null, primarySubsampleShift, consumeDonors, buildSecondary13, threads);
    }

    private static EimdlSyncmerIndex buildInternal(
            Packed2BitSequence[] donors,
            boolean[] indexedDonorMask,
            int primarySubsampleShift,
            boolean consumeDonors,
            boolean buildSecondary13,
            int threads) {
        int n = donors.length;
        long indexedBases = 0L;
        for (int d = 0; d < n; d++) {
            if (indexedDonorMask == null || indexedDonorMask[d]) {
                indexedBases += donors[d].length;
            }
        }

        int[] donorOffsets = new int[n + 1];
        Packed2BitSequence bank =
                consumeDonors
                        ? Packed2BitSequence.concatenateAndConsume(donors, donorOffsets)
                        : Packed2BitSequence.concatenate(donors, donorOffsets);
        int totalBases = bank.length;

        int numBlocks = (totalBases >>> BLOCK_SHIFT) + 2;
        int[] posToDonorFloor = new int[numBlocks];
        int curD = 0;
        for (int b = 0; b < numBlocks; b++) {
            int bStart = b << BLOCK_SHIFT;
            while (curD + 1 < n && donorOffsets[curD + 1] <= bStart) {
                curD++;
            }
            posToDonorFloor[b] = curD;
        }

        long estPrimaryLong =
                Math.max(
                        1024L,
                        (indexedBases / Packed2BitSequence.W_PRIMARY) >>> primarySubsampleShift);
        int estPrimary = (int) Math.min(Integer.MAX_VALUE, estPrimaryLong);
        int maxPrimaryBits =
                primarySubsampleShift == 0 ? 20 : Math.max(16, 20 - primarySubsampleShift);
        int bucketBits = chooseBucketBits(estPrimary, 12, maxPrimaryBits);

        CsrTables primary =
                buildCsrPass(
                        bank,
                        donorOffsets,
                        indexedDonorMask,
                        Packed2BitSequence.K_PRIMARY,
                        Packed2BitSequence.S_PRIMARY,
                        Packed2BitSequence.T_PRIMARY,
                        bucketBits,
                        primarySubsampleShift,
                        threads);

        CsrTables secondary = null;
        int secMask = 0;
        if (buildSecondary13
                && primarySubsampleShift == 0
                && indexedBases > 0
                && indexedBases <= 80_000_000L) {
            int estSec = (int) Math.max(1024L, indexedBases / (Packed2BitSequence.W_SECONDARY * 2L));
            int secBits = chooseBucketBits(estSec, 12, 19);
            secMask = (1 << secBits) - 1;
            secondary =
                    buildCsrPass(
                            bank,
                            donorOffsets,
                            indexedDonorMask,
                            Packed2BitSequence.K_SECONDARY,
                            Packed2BitSequence.S_SECONDARY,
                            Packed2BitSequence.T_SECONDARY,
                            secBits,
                            1,
                            threads);
        }

        return new EimdlSyncmerIndex(
                bank,
                n,
                donorOffsets,
                posToDonorFloor,
                primarySubsampleShift,
                bucketBits,
                primary.bucketOffsets,
                primary.posAndStrand,
                primary.sig16,
                secMask,
                secondary != null ? secondary.bucketOffsets : null,
                secondary != null ? secondary.posAndStrand : null,
                secondary != null ? secondary.sig16 : null);
    }

    private record CsrTables(int[] bucketOffsets, int[] posAndStrand, short[] sig16) {}

    private record DonorSyncmers(int[] buckets, short[] sigs, int[] packedPosStrand, int count) {}

    private static CsrTables buildCsrPass(
            Packed2BitSequence bank,
            int[] donorOffsets,
            boolean[] indexedDonorMask,
            int k,
            int s,
            int t,
            int bucketBits,
            int subsampleShift,
            int threads) {
        boolean ampleHeap = Runtime.getRuntime().maxMemory() > 512L * 1024 * 1024;
        if (ampleHeap && threads > 1 && donorOffsets.length - 1 > 1) {
            return buildCsrPassParallel(
                    bank, donorOffsets, indexedDonorMask, k, s, t, bucketBits, subsampleShift, threads);
        }
        return buildCsrPassZeroAlloc(
                bank, donorOffsets, indexedDonorMask, k, s, t, bucketBits, subsampleShift);
    }

    /**
     * Multi-core parallel CSR builder for ample-heap runs: extracts syncmers per donor in parallel
     * in L1 cache and populates disjoint L3-cache-resident bucket stripes across cores, producing
     * CSR arrays bit-for-bit identical to {@link #buildCsrPassZeroAlloc}.
     */
    private static CsrTables buildCsrPassParallel(
            Packed2BitSequence bank,
            int[] donorOffsets,
            boolean[] indexedDonorMask,
            int k,
            int s,
            int t,
            int bucketBits,
            int subsampleShift,
            int threads) {
        final int n = donorOffsets.length - 1;
        final int numBuckets = 1 << bucketBits;
        final int mask = numBuckets - 1;
        final int subMask = (1 << subsampleShift) - 1;

        DonorSyncmers[] perDonor = new DonorSyncmers[n];
        ForkJoinPool pool = new ForkJoinPool(Math.min(threads, Math.max(n, 4)));
        try {
            pool.submit(
                            () ->
                                    IntStream.range(0, n)
                                            .parallel()
                                            .forEach(
                                                    d -> {
                                                        if (indexedDonorMask != null
                                                                && !indexedDonorMask[d]) {
                                                            perDonor[d] =
                                                                    new DonorSyncmers(
                                                                            new int[0],
                                                                            new short[0],
                                                                            new int[0],
                                                                            0);
                                                            return;
                                                        }
                                                        int start = donorOffsets[d];
                                                        int len = donorOffsets[d + 1] - start;
                                                        int cap =
                                                                Math.max(
                                                                        64,
                                                                        ((len / (k - s))
                                                                                        >>> subsampleShift)
                                                                                + 64);
                                                        int[][] bufs = {new int[cap], new int[cap]};
                                                        short[][] sbuf = {new short[cap]};
                                                        int[] cnt = {0};

                                                        bank.forEachCanonicalOpenSyncmer(
                                                                start,
                                                                len,
                                                                k,
                                                                s,
                                                                t,
                                                                (relPos, canonK, strand) -> {
                                                                    int h =
                                                                            Packed2BitSequence.mix32(
                                                                                    canonK);
                                                                    if ((h & subMask) != 0) {
                                                                        return;
                                                                    }
                                                                    int idx = cnt[0];
                                                                    if (idx == bufs[0].length) {
                                                                        int ncap =
                                                                                idx
                                                                                        + (idx >>> 1)
                                                                                        + 64;
                                                                        bufs[0] =
                                                                                Arrays.copyOf(
                                                                                        bufs[0],
                                                                                        ncap);
                                                                        bufs[1] =
                                                                                Arrays.copyOf(
                                                                                        bufs[1],
                                                                                        ncap);
                                                                        sbuf[0] =
                                                                                Arrays.copyOf(
                                                                                        sbuf[0],
                                                                                        ncap);
                                                                    }
                                                                    bufs[0][idx] =
                                                                            (h >>> subsampleShift)
                                                                                    & mask;
                                                                    sbuf[0][idx] =
                                                                            (short)
                                                                                    ((canonK >>> 16)
                                                                                            ^ (h
                                                                                                    >>> 16));
                                                                    bufs[1][idx] =
                                                                            (start + relPos)
                                                                                    | (strand << 31);
                                                                    cnt[0] = idx + 1;
                                                                });
                                                        perDonor[d] =
                                                                new DonorSyncmers(
                                                                        bufs[0],
                                                                        sbuf[0],
                                                                        bufs[1],
                                                                        cnt[0]);
                                                    }))
                    .join();

            final int stripes = Integer.highestOneBit(Math.min(threads, 4));
            final int stripeBuckets = numBuckets / stripes;
            final int[] bucketCounts = new int[numBuckets];

            pool.submit(
                            () ->
                                    IntStream.range(0, stripes)
                                            .parallel()
                                            .forEach(
                                                    st -> {
                                                        int bLo = st * stripeBuckets;
                                                        int bHi = bLo + stripeBuckets;
                                                        int[] lastDonor = new int[stripeBuckets];
                                                        byte[] occ = new byte[stripeBuckets];
                                                        Arrays.fill(lastDonor, -1);
                                                        for (int d = 0; d < n; d++) {
                                                            DonorSyncmers ds = perDonor[d];
                                                            int[] bArr = ds.buckets;
                                                            int c = ds.count;
                                                            for (int i = 0; i < c; i++) {
                                                                int b = bArr[i];
                                                                if (b < bLo || b >= bHi) continue;
                                                                int rel = b - bLo;
                                                                if (lastDonor[rel] != d) {
                                                                    lastDonor[rel] = d;
                                                                    occ[rel] = 1;
                                                                    bucketCounts[b]++;
                                                                } else if ((occ[rel] & 0xFF)
                                                                        < MAX_WITHIN_DONOR_BUCKET_OCC) {
                                                                    occ[rel]++;
                                                                    bucketCounts[b]++;
                                                                } else {
                                                                    bArr[i] = -1;
                                                                }
                                                            }
                                                        }
                                                    }))
                    .join();

            int[] bucketOffsets = new int[numBuckets + 1];
            int running = 0;
            for (int b = 0; b < numBuckets; b++) {
                bucketOffsets[b] = running;
                running += bucketCounts[b];
            }
            bucketOffsets[numBuckets] = running;

            int[] posAndStrand = new int[running];
            short[] sig16 = new short[running];

            pool.submit(
                            () ->
                                    IntStream.range(0, stripes)
                                            .parallel()
                                            .forEach(
                                                    st -> {
                                                        int bLo = st * stripeBuckets;
                                                        int bHi = bLo + stripeBuckets;
                                                        int[] cursor =
                                                                Arrays.copyOfRange(
                                                                        bucketOffsets, bLo, bHi);
                                                        for (int d = 0; d < n; d++) {
                                                            DonorSyncmers ds = perDonor[d];
                                                            int[] bArr = ds.buckets;
                                                            short[] sArr = ds.sigs;
                                                            int[] pArr = ds.packedPosStrand;
                                                            int c = ds.count;
                                                            for (int i = 0; i < c; i++) {
                                                                int b = bArr[i];
                                                                if (b < bLo || b >= bHi) continue;
                                                                int dst = cursor[b - bLo]++;
                                                                posAndStrand[dst] = pArr[i];
                                                                sig16[dst] = sArr[i];
                                                            }
                                                        }
                                                    }))
                    .join();

            return new CsrTables(bucketOffsets, posAndStrand, sig16);
        } finally {
            pool.shutdown();
        }
    }

    /**
     * Two-pass zero-materialization CSR builder: streams SWAR canonical open syncmers directly from
     * {@code bank} without allocating intermediate per-donor syncmer buffers (~340 MB saved on 120
     * Mbp collections while producing bit-identical CSR tables).
     */
    private static CsrTables buildCsrPassZeroAlloc(
            Packed2BitSequence bank,
            int[] donorOffsets,
            boolean[] indexedDonorMask,
            int k,
            int s,
            int t,
            int bucketBits,
            int subsampleShift) {
        final int n = donorOffsets.length - 1;
        final int numBuckets = 1 << bucketBits;
        final int mask = numBuckets - 1;
        final int subMask = (1 << subsampleShift) - 1;

        // Pass 1: Count entries per bucket while capping pathological within-donor satellite floods
        int[] bucketCounts = new int[numBuckets];
        int[] lastDonorSeen = new int[numBuckets];
        byte[] donorOccInBucket = new byte[numBuckets];
        Arrays.fill(lastDonorSeen, -1);

        for (int d = 0; d < n; d++) {
            if (indexedDonorMask != null && !indexedDonorMask[d]) {
                continue;
            }
            final int donorId = d;
            final int start = donorOffsets[d];
            final int len = donorOffsets[d + 1] - start;
            final int[] countsRef = bucketCounts;
            bank.forEachCanonicalOpenSyncmer(
                    start,
                    len,
                    k,
                    s,
                    t,
                    (relPos, canonK, strand) -> {
                        int h = Packed2BitSequence.mix32(canonK);
                        if ((h & subMask) != 0) {
                            return;
                        }
                        int b = (h >>> subsampleShift) & mask;
                        if (lastDonorSeen[b] != donorId) {
                            lastDonorSeen[b] = donorId;
                            donorOccInBucket[b] = 1;
                            countsRef[b]++;
                        } else if ((donorOccInBucket[b] & 0xFF) < MAX_WITHIN_DONOR_BUCKET_OCC) {
                            donorOccInBucket[b]++;
                            countsRef[b]++;
                        }
                    });
        }

        int[] bucketOffsets = new int[numBuckets + 1];
        int running = 0;
        for (int b = 0; b < numBuckets; b++) {
            bucketOffsets[b] = running;
            running += bucketCounts[b];
        }
        bucketOffsets[numBuckets] = running;
        bucketCounts = null; // allow GC before allocating CSR payload arrays

        int[] posAndStrand = new int[running];
        short[] sig16 = new short[running];
        int[] writeCursor = Arrays.copyOf(bucketOffsets, numBuckets);
        Arrays.fill(lastDonorSeen, -1);

        // Pass 2: Populate CSR payload directly in deterministic donor order
        for (int d = 0; d < n; d++) {
            if (indexedDonorMask != null && !indexedDonorMask[d]) {
                continue;
            }
            final int donorId = d;
            final int start = donorOffsets[d];
            final int len = donorOffsets[d + 1] - start;
            bank.forEachCanonicalOpenSyncmer(
                    start,
                    len,
                    k,
                    s,
                    t,
                    (relPos, canonK, strand) -> {
                        int h = Packed2BitSequence.mix32(canonK);
                        if ((h & subMask) != 0) {
                            return;
                        }
                        int b = (h >>> subsampleShift) & mask;
                        boolean accept = false;
                        if (lastDonorSeen[b] != donorId) {
                            lastDonorSeen[b] = donorId;
                            donorOccInBucket[b] = 1;
                            accept = true;
                        } else if ((donorOccInBucket[b] & 0xFF) < MAX_WITHIN_DONOR_BUCKET_OCC) {
                            donorOccInBucket[b]++;
                            accept = true;
                        }
                        if (accept) {
                            int dst = writeCursor[b]++;
                            posAndStrand[dst] = (start + relPos) | (strand << 31);
                            sig16[dst] = (short) ((canonK >>> 16) ^ (h >>> 16));
                        }
                    });
        }

        return new CsrTables(bucketOffsets, posAndStrand, sig16);
    }

    private static int chooseBucketBits(int estEntries, int minBits, int maxBits) {
        int bits = minBits;
        while (bits < maxBits && (1 << bits) < (estEntries >>> 1)) {
            bits++;
        }
        return bits;
    }

    /** Fast O(1) resolution of {@code globalBasePos} to {@code donorId} in {@code [0, numDonors)}. */
    public int donorOfGlobalPos(int globalPos) {
        int d = posToDonorFloor[globalPos >>> BLOCK_SHIFT];
        while (d + 1 < numDonors && globalPos >= donorOffsets[d + 1]) {
            d++;
        }
        return d;
    }

    public boolean acceptsPrimary(long canonK16) {
        return primarySubsampleMask == 0
                || (Packed2BitSequence.mix32(canonK16) & primarySubsampleMask) == 0;
    }

    public int bucketForPrimary(long canonK16) {
        return (Packed2BitSequence.mix32(canonK16) >>> primarySubsampleShift) & bucketMask;
    }

    public short sigForPrimary(long canonK16) {
        int h = Packed2BitSequence.mix32(canonK16);
        return (short) ((canonK16 >>> 16) ^ (h >>> 16));
    }

    public int bucketStart(int bucket) {
        return bucketOffsets[bucket];
    }

    public int bucketEnd(int bucket) {
        return bucketOffsets[bucket + 1];
    }

    public short entrySig(int idx) {
        return sig16[idx];
    }

    public int entryPackedPosAndStrand(int idx) {
        return posAndStrand[idx];
    }

    public boolean hasSecondary() {
        return secBucketOffsets != null;
    }

    public int secBucketFor(long canonK13) {
        int h = Packed2BitSequence.mix32(canonK13);
        if ((h & 1) != 0) {
            return -1;
        }
        return (h >>> 1) & secBucketMask;
    }

    public short secSigFor(long canonK13) {
        int h = Packed2BitSequence.mix32(canonK13);
        return (short) ((canonK13 >>> 16) ^ (h >>> 16));
    }

    public int secBucketStart(int bucket) {
        return secBucketOffsets[bucket];
    }

    public int secBucketEnd(int bucket) {
        return secBucketOffsets[bucket + 1];
    }

    public short secEntrySig(int idx) {
        return secSig16[idx];
    }

    public int secEntryPackedPosAndStrand(int idx) {
        return secPosAndStrand[idx];
    }

    /**
     * Ultra-compact L2/L3-cache-resident 4-way set-associative ring table for bidirectional
     * self-repeats ({@code + / -} strands) within the target sequence being parsed.
     *
     * <p>Footprint: {@code 2^15 * 4 = 131,072} slots (512 KB ints + 256 KB shorts = 768 KB per
     * worker thread!).
     */
    public static final class SelfSyncmerTable {
        public static final int WAYS = 4;
        private final int setMask;
        private final int[] packedPosStrand; // -1 = empty; bit 31 = canonicalStrand, bits 0..30 = pos
        private final short[] sigs;
        private final byte[] clockHand;

        public SelfSyncmerTable(int setBits) {
            int numSets = 1 << setBits;
            this.setMask = numSets - 1;
            this.packedPosStrand = new int[numSets * WAYS];
            Arrays.fill(this.packedPosStrand, -1);
            this.sigs = new short[numSets * WAYS];
            this.clockHand = new byte[numSets];
        }

        public void insert(long canonK, int pos, int canonicalStrand) {
            int h = Packed2BitSequence.mix32(canonK);
            int set = h & setMask;
            short sig = (short) ((canonK >>> 16) ^ (h >>> 16));
            int base = set << 2;
            // Update in-place if same signature was very recently inserted closely (< 64 bp)
            for (int w = 0; w < WAYS; w++) {
                int s = base + w;
                int prev = packedPosStrand[s];
                if (prev != -1 && sigs[s] == sig && pos - (prev & 0x7FFFFFFF) < 64) {
                    packedPosStrand[s] = pos | (canonicalStrand << 31);
                    return;
                }
            }
            int hand = clockHand[set] & 3;
            int slot = base + hand;
            packedPosStrand[slot] = pos | (canonicalStrand << 31);
            sigs[slot] = sig;
            clockHand[set] = (byte) ((hand + 1) & 3);
        }

        public int setBase(long canonK) {
            int h = Packed2BitSequence.mix32(canonK);
            return (h & setMask) << 2;
        }

        public short sigFor(long canonK) {
            int h = Packed2BitSequence.mix32(canonK);
            return (short) ((canonK >>> 16) ^ (h >>> 16));
        }

        public short slotSig(int slot) {
            return sigs[slot];
        }

        public int slotPackedPosStrand(int slot) {
            return packedPosStrand[slot];
        }
    }
}
