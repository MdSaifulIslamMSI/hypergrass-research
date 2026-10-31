package hypergrass;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Random;

/**
 * Deterministic, dependency-free differential verification of HaplotypeDonorIndex.
 * Run with --bench to also run a bounded synthetic intersection microbenchmark.
 */
public final class VerifyDonorIndex {

    private static final int MAX_DONORS = 512;
    private static final int[] EMPTY = new int[0];
    private static final long SEED = 0x5EED_D0A0_2026L;
    private static final Selection NO_DONOR = new Selection(-1, 0, 8);
    private static final int BENCH_REPETITIONS = 4;
    private static final int BENCH_WARMUPS = 4;
    private static final int BENCH_SAMPLES = 7;
    private static volatile long benchmarkSink;

    private long queriesChecked;
    private long donorCountsChecked;

    private VerifyDonorIndex() {
    }

    public static void main(String[] args) {
        if (args.length > 1 || (args.length == 1 && !"--bench".equals(args[0]))) {
            throw new IllegalArgumentException("Usage: java hypergrass.VerifyDonorIndex [--bench]");
        }

        VerifyDonorIndex verifier = new VerifyDonorIndex();
        long start = System.nanoTime();
        verifier.verifyPackedBoundaries();
        verifier.verifyExhaustiveSubsets();
        verifier.verifyIncrementalAndCap();
        verifier.verifyBitplaneCarriesAndLanes();
        verifier.verifyIndexBudgetTransition();
        verifier.verifyRandomizedPools();
        verifier.verifySyntheticPools();
        verifier.verifyScoring();
        System.out.printf(Locale.ROOT,
                "PASS: 8 deterministic suites; %,d queries, %,d donor counts and %,d selections checked "
                        + "(seed=0x%x, %.3f s).%n",
                verifier.queriesChecked, verifier.donorCountsChecked, verifier.queriesChecked,
                SEED, (System.nanoTime() - start) / 1e9);

        if (args.length == 1) {
            runBenchmarks();
        }
    }

    private void verifyPackedBoundaries() {
        int max = Integer.MAX_VALUE;
        int[][] cases = {
                EMPTY, {0}, {1}, {max}, {0, 1, 2, 3, 4, 5, 6, 7},
                {max - 3, max - 2, max - 1, max}, {0, 1, 3, 4, 5, max - 1, max},
                {4, 8, 12}, {5, 9, 13}, {6, 10, 14}, {7, 11, 15}
        };
        for (int capacity : new int[]{0, 1, 15, 16, 17, MAX_DONORS}) {
            TestPool pool = new TestPool(capacity);
            for (int[] query : cases) {
                checkQuery("empty pool, capacity=" + capacity, pool, query);
            }
            for (int i = 0; i < 24; i++) {
                pool.add(cases[i % cases.length]);
                for (int[] query : cases) {
                    checkQuery("packed boundaries, capacity=" + capacity + ", addition=" + i, pool, query);
                }
            }
        }
    }

    private void verifyExhaustiveSubsets() {
        int[] universe = {0, 1, 2, 3, Integer.MAX_VALUE - 1, Integer.MAX_VALUE};
        int[][] subsets = new int[1 << universe.length][];
        TestPool pool = new TestPool(MAX_DONORS);
        for (int mask = 0; mask < subsets.length; mask++) {
            int[] subset = new int[Integer.bitCount(mask)];
            int next = 0;
            for (int bit = 0; bit < universe.length; bit++) {
                if ((mask & (1 << bit)) != 0) {
                    subset[next++] = universe[bit];
                }
            }
            subsets[mask] = subset;
            pool.add(subset);
        }
        // The first query must catch up with an entire unqueried donor pool.
        for (int mask = 0; mask < subsets.length; mask++) {
            checkQuery("all subsets, mask=" + mask, pool, subsets[mask]);
        }
        for (int mask = 0; mask < subsets.length; mask++) {
            pool.add(subsets[mask]);
            checkQuery("subset incremental catchup, mask=" + mask, pool, subsets[mask]);
            checkQuery("subset complement, mask=" + mask, pool, subsets[mask ^ (subsets.length - 1)]);
            checkQuery("subset reset", pool, EMPTY);
        }
    }

    private void verifyIncrementalAndCap() {
        TestPool pool = new TestPool(MAX_DONORS);
        checkQuery("before first addition", pool, EMPTY);
        for (int i = 0; i < MAX_DONORS + 16; i++) {
            int[] donor = i % 13 == 0 ? EMPTY : new int[]{0, 1, (i + 8) << 2 | 1, Integer.MAX_VALUE};
            pool.add(donor);
            if (i >= 14) {
                String label = "incremental/cap, addition=" + i;
                checkQuery(label, pool, donor);
                checkQuery(label + " shared", pool, new int[]{0, 1, Integer.MAX_VALUE});
                checkQuery(label + " repeated", pool, donor);
                checkQuery(label + " reset", pool, (i & 1) == 0 ? EMPTY : new int[]{Integer.MAX_VALUE - 4});
            }
        }
        if (pool.index.size() != MAX_DONORS) {
            throw new AssertionError("The donor pool did not stop at 512 donors");
        }
    }

    private void verifyBitplaneCarriesAndLanes() {
        TestPool pool = new TestPool(MAX_DONORS);
        int[] longQuery = range(0, 1089);
        int[] boundaryDonors = {63, 64, 127, 511};
        for (int target : new int[]{63, 64, 65, 127, 128, 129, 511, MAX_DONORS}) {
            while (pool.donors.size() < target) {
                int donorId = pool.donors.size();
                int length = 1025 - donorId % 3;
                // Stagger memberships within and between words, with different carry chains per lane.
                int[] donor = Arrays.copyOf(range((donorId * 17) % 65, length), length + 1);
                donor[length] = ((4096 + donorId) << 2) | 1;
                pool.add(donor);
            }
            String label = "bitplane lanes/catchup, donors=" + target;
            checkQuery(label + " long", pool, longQuery);
            checkQuery(label + " short", pool, new int[]{63});
            checkQuery(label + " empty", pool, EMPTY);
            for (int donorId : boundaryDonors) {
                if (donorId < target) {
                    checkQuery(label + ", isolated donor=" + donorId, pool,
                            new int[]{((4096 + donorId) << 2) | 1});
                }
            }
        }
        // Donor zero contains every queried value, so both query length and exact count cross each power.
        for (int power = 2; power <= 1024; power <<= 1) {
            for (int length = power - 1; length <= power + 1; length++) {
                checkQuery("bitplane carry, query length=" + length, pool, range(0, length));
            }
        }
        checkQuery("bitplane high-carry reset short", pool, new int[]{0});
        checkQuery("bitplane high-carry reset empty", pool, EMPTY);
        checkQuery("bitplane high-carry repeat long", pool, longQuery);
    }

    private void verifyIndexBudgetTransition() {
        int[] sharedDonor = range(0, 4096);
        int[] extraDonor = {4096};
        int[] mixedQuery = {0, 2048, 4095, 4096};
        int[] disjointQuery = {4097, 4098, 4099, 4100};
        for (boolean incremental : new boolean[]{true, false}) {
            TestPool pool = new TestPool(33);
            String label = incremental ? "incremental index budget" : "all-at-once index budget";
            // Shared immutable input still counts as 4096 stored memberships for each donor.
            for (int d = 0; d < 31; d++) {
                pool.add(sharedDonor);
            }
            if (incremental) {
                for (int length = 1; length < 4; length++) {
                    checkQuery(label + " tiny before build, length=" + length, pool, range(0, length));
                }
                checkQuery(label + " build below budget", pool, sharedDonor);
            }
            pool.add(sharedDonor);
            if (incremental) {
                checkQuery(label + " catchup at 131072 memberships", pool, sharedDonor);
            }
            pool.add(extraDonor);
            // The new donor has a distinct variant, exposing stale counts after the fallback transition.
            int[][] queries = {mixedQuery, sharedDonor, extraDonor, EMPTY, disjointQuery, mixedQuery, sharedDonor};
            for (int q = 0; q < queries.length; q++) {
                checkQuery(label + " over budget, query=" + q, pool, queries[q]);
            }
        }
    }

    private void verifyRandomizedPools() {
        Random random = new Random(SEED);
        for (int domain : new int[]{32, 4096, Integer.MAX_VALUE}) {
            for (int capacity : new int[]{0, 1, 15, 16, 17, MAX_DONORS}) {
                TestPool pool = new TestPool(capacity);
                for (int i = 0; i < capacity + 8; i++) {
                    int[] donor = randomVariants(random, domain);
                    pool.add(donor);
                    String label = "random domain=" + domain + ", capacity=" + capacity + ", addition=" + i;
                    checkQuery(label + " self", pool, donor);
                    checkQuery(label + " independent", pool, randomVariants(random, domain));
                    if (!pool.donors.isEmpty()) {
                        checkQuery(label + " earlier donor", pool,
                                pool.donors.get(random.nextInt(pool.donors.size())));
                    }
                    if ((i & 15) == 0) {
                        checkQuery(label + " empty", pool, EMPTY);
                    }
                }
            }
        }
    }

    private void verifySyntheticPools() {
        for (PoolKind kind : PoolKind.values()) {
            Workload workload = createWorkload(kind);
            TestPool pool = new TestPool(MAX_DONORS);
            for (int target : new int[]{15, 16, 17, 64, 255, MAX_DONORS}) {
                while (pool.donors.size() < target) {
                    pool.add(workload.donors[pool.donors.size()]);
                }
                for (int[] query : workload.queries) {
                    checkQuery(kind.label + ", donors=" + target, pool, query);
                }
                checkQuery(kind.label + " self", pool, workload.donors[target - 1]);
                checkQuery(kind.label + " clear", pool, EMPTY);
                checkQuery(kind.label + " repeat", pool, workload.queries[0]);
            }
        }
    }

    private void verifyScoring() {
        checkSelection("no donors", new int[0][], new int[]{1}, NO_DONOR);
        checkSelection("no common variants", new int[][]{EMPTY, {2}}, new int[]{1}, NO_DONOR);
        checkSelection("score below threshold", new int[][]{{1, 2, 3}}, new int[]{1}, NO_DONOR);
        checkSelection("score equals threshold", new int[][]{{1, 2}}, new int[]{1}, NO_DONOR);
        checkSelection("score one above threshold", new int[][]{range(1, 19)}, new int[]{1, 2},
                new Selection(0, 2, 9));
        checkSelection("full donor omits mask", new int[][]{range(1, 20)}, range(1, 21),
                new Selection(0, 20, 352));
        checkSelection("equal scores, different variants", new int[][]{{1, 2}, {1, 3}}, new int[]{1, 2, 3},
                new Selection(0, 2, 28));
        checkSelection("equal scores, different common counts", new int[][]{{1}, range(1, 18)}, new int[]{1, 2},
                new Selection(0, 1, 10));
        checkSelection("equal scores, reverse order", new int[][]{range(1, 18), {1}}, new int[]{1, 2},
                new Selection(0, 2, 10));
        checkSelection("most common is not best score", new int[][]{range(1, 31), {1}}, new int[]{1, 2},
                new Selection(1, 1, 10));
        checkSelection("later strictly better donor", new int[][]{{1}, {1, 2}}, new int[]{1, 2},
                new Selection(1, 2, 28));

        int[][] acrossBoundary = new int[18][];
        Arrays.fill(acrossBoundary, EMPTY);
        acrossBoundary[15] = new int[]{1, 2};
        acrossBoundary[16] = new int[]{1, 3};
        checkSelection("first tie across index boundary", acrossBoundary, new int[]{1, 2, 3},
                new Selection(15, 2, 28));
    }

    private void checkSelection(String label, int[][] donors, int[] query, Selection expected) {
        for (int minimumSize : new int[]{0, 15, 16, 17, MAX_DONORS}) {
            TestPool pool = new TestPool(MAX_DONORS);
            for (int[] donor : donors) {
                pool.add(donor);
            }
            while (pool.donors.size() < minimumSize) {
                pool.add(EMPTY);
            }
            Selection actual = checkQuery(label + ", donors=" + pool.donors.size(), pool, query);
            if (!expected.equals(actual)) {
                throw new AssertionError(label + ": expected " + expected + ", got " + actual);
            }
        }
    }

    private Selection checkQuery(String label, TestPool pool, int[] query) {
        int[] original = query.clone();
        int size = pool.donors.size();
        int[] counts = pool.index.countIntersections(query);
        if (counts.length < size || pool.index.size() != size || !Arrays.equals(original, query)) {
            throw new AssertionError(label + ": invalid count buffer, changed size, or mutated query");
        }
        // Only the active prefix is specified; buffer identity need not survive an addition.
        if (pool.lastSize == size && counts != pool.lastCounts) {
            throw new AssertionError(label + ": counts buffer was not reused at a fixed pool size");
        }
        pool.lastCounts = counts;
        pool.lastSize = size;

        Selection expected = NO_DONOR;
        for (int d = 0; d < size; d++) {
            int[] donor = pool.donors.get(d);
            int common = intersection(query, donor);
            if (counts[d] != common || !Arrays.equals(donor, pool.index.get(d))) {
                throw new AssertionError(label + ": donor=" + d + ", expected count=" + common
                        + ", actual count=" + counts[d] + " (or donor contents changed)");
            }
            if (common != 0) {
                int score = common == donor.length ? 18 * common - 8 : 18 * common - donor.length - 8;
                if (score > expected.score) {
                    expected = new Selection(d, common, score);
                }
            }
        }
        Selection actual = indexedSelection(pool.index, counts);
        if (!expected.equals(actual)) {
            throw new AssertionError(label + ": expected selection " + expected + ", got " + actual);
        }
        queriesChecked++;
        donorCountsChecked += size;
        return actual;
    }

    private static Selection indexedSelection(HaplotypeDonorIndex index, int[] counts) {
        int bestDonor = -1;
        int bestCommon = 0;
        int bestScore = 8;
        for (int d = 0; d < index.size(); d++) {
            int common = counts[d];
            if (common == 0) {
                continue;
            }
            int length = index.get(d).length;
            int maskBits = common == length ? 0 : length;
            int score = common * 18 - maskBits - 8;
            if (score > bestScore) {
                bestDonor = d;
                bestCommon = common;
                bestScore = score;
            }
        }
        return new Selection(bestDonor, bestCommon, bestScore);
    }

    private static int intersection(int[] a, int[] b) {
        int i = 0;
        int j = 0;
        int common = 0;
        while (i < a.length && j < b.length) {
            if (a[i] < b[j]) {
                i++;
            } else if (a[i] > b[j]) {
                j++;
            } else {
                common++;
                i++;
                j++;
            }
        }
        return common;
    }

    private static int[] randomVariants(Random random, int domain) {
        int[] values = new int[random.nextInt(65)];
        for (int i = 0; i < values.length; i++) {
            values[i] = random.nextInt(domain);
        }
        Arrays.sort(values);
        int unique = 0;
        for (int value : values) {
            if (unique == 0 || values[unique - 1] != value) {
                values[unique++] = value;
            }
        }
        return Arrays.copyOf(values, unique);
    }

    private static int[] range(int start, int length) {
        int[] values = new int[length];
        for (int i = 0; i < length; i++) {
            values[i] = start + i;
        }
        return values;
    }

    private static Workload createWorkload(PoolKind kind) {
        int[][] donors = new int[MAX_DONORS][];
        int[][] queries = new int[16][];
        for (int d = 0; d < donors.length; d++) {
            donors[d] = syntheticVariants(kind, d, false);
        }
        for (int q = 0; q < queries.length; q++) {
            queries[q] = syntheticVariants(kind, (q * 37 + 11) % donors.length, true);
        }
        return new Workload(kind.label, donors, queries);
    }

    private static int[] syntheticVariants(PoolKind kind, int donorId, boolean query) {
        int shared = switch (kind) {
            case SPARSE_RELATED -> 24;
            case DENSE_SHARED -> 96;
            case DISJOINT -> 0;
        };
        int[] variants = new int[128];
        for (int i = 0; i < shared; i++) {
            int pos = kind == PoolKind.SPARSE_RELATED ? (donorId % 64) * 128 + i : i;
            variants[i] = (pos << 2) | 1;
        }
        for (int i = 0; i < variants.length - shared; i++) {
            int pos = 1_000_000 + donorId * 256 + i;
            if (query && kind != PoolKind.DISJOINT && i >= 8) {
                pos += 3_000_000;
            }
            int allele = query && kind == PoolKind.DISJOINT ? 2 : 1;
            variants[shared + i] = (pos << 2) | allele;
        }
        return variants;
    }

    private static void runBenchmarks() {
        System.out.println("\nSynthetic intersection microbenchmark (not end-to-end compression or real genomic data).");
        System.out.printf(Locale.ROOT, "%d warmup batches, %d measured batches; median times, alternating order.%n",
                BENCH_WARMUPS, BENCH_SAMPLES);
        System.out.println("Both paths consume every donor count with the same weighted checksum; oracle scans the full pool.");
        System.out.println("Warmed times exclude data generation, construction and first query; startup must be amortized.");
        System.out.println("No speed threshold is asserted; results depend on workload, JIT, GC and hardware.");

        for (PoolKind kind : PoolKind.values()) {
            Workload workload = createWorkload(kind);
            long start = System.nanoTime();
            HaplotypeDonorIndex index = new HaplotypeDonorIndex(MAX_DONORS);
            for (int[] donor : workload.donors) {
                index.add(donor);
            }
            long constructionNanos = System.nanoTime() - start;
            start = System.nanoTime();
            int[] firstCounts = index.countIntersections(workload.queries[0]);
            long firstQueryNanos = System.nanoTime() - start;
            benchmarkSink = checksum(firstCounts, index.size());

            long matchesPerPass = 0;
            long queryVariantsPerPass = 0;
            for (int[] query : workload.queries) {
                int[] counts = index.countIntersections(query);
                queryVariantsPerPass += query.length;
                for (int d = 0; d < workload.donors.length; d++) {
                    int common = intersection(query, workload.donors[d]);
                    if (counts[d] != common) {
                        throw new AssertionError("Synthetic benchmark count mismatch: " + workload.label + ", donor=" + d);
                    }
                    matchesPerPass += common;
                }
            }
            long expectedChecksum = oracleBatch(workload);
            for (int warmup = 0; warmup < BENCH_WARMUPS; warmup++) {
                checkBatch(timedBatch(workload, index, (warmup & 1) == 0), expectedChecksum);
                checkBatch(timedBatch(workload, index, (warmup & 1) != 0), expectedChecksum);
            }
            long[] oracleTimes = new long[BENCH_SAMPLES];
            long[] indexedTimes = new long[BENCH_SAMPLES];
            for (int sample = 0; sample < BENCH_SAMPLES; sample++) {
                boolean oracleFirst = (sample & 1) == 0;
                TimedBatch first = timedBatch(workload, index, oracleFirst);
                TimedBatch second = timedBatch(workload, index, !oracleFirst);
                checkBatch(first, expectedChecksum);
                checkBatch(second, expectedChecksum);
                oracleTimes[sample] = oracleFirst ? first.nanos : second.nanos;
                indexedTimes[sample] = oracleFirst ? second.nanos : first.nanos;
            }
            Arrays.sort(oracleTimes);
            Arrays.sort(indexedTimes);
            long queriesPerBatch = (long) workload.queries.length * BENCH_REPETITIONS;
            long donorVariants = 0;
            for (int[] donor : workload.donors) {
                donorVariants += donor.length;
            }
            System.out.printf(Locale.ROOT, "%n  Synthetic %-14s | %,d donors, %,d stored packed variants%n",
                    workload.label, index.size(), donorVariants);
            System.out.printf(Locale.ROOT,
                    "    Per algorithm/batch: %,d queries (%d vectors x %d repeats), %,d donor pairs/counts consumed,%n"
                            + "      %,d query packed values, %,d common variants summed across donor pairs.%n",
                    queriesPerBatch, workload.queries.length, BENCH_REPETITIONS, queriesPerBatch * index.size(),
                    queryVariantsPerPass * BENCH_REPETITIONS, matchesPerPass * BENCH_REPETITIONS);
            System.out.printf(Locale.ROOT,
                    "    One-shot construction/adds: %.3f ms; first query (includes lazy build): %.3f ms%n",
                    constructionNanos / 1e6, firstQueryNanos / 1e6);
            System.out.printf(Locale.ROOT,
                    "    Median/batch: oracle full scan=%.3f ms, countIntersections=%.3f ms; both checksums=%d%n",
                    oracleTimes[BENCH_SAMPLES / 2] / 1e6, indexedTimes[BENCH_SAMPLES / 2] / 1e6, expectedChecksum);
        }
    }

    private static TimedBatch timedBatch(Workload workload, HaplotypeDonorIndex index, boolean oracle) {
        long start = System.nanoTime();
        long result = oracle ? oracleBatch(workload) : indexedBatch(workload, index);
        long elapsed = System.nanoTime() - start;
        benchmarkSink = result;
        return new TimedBatch(elapsed, result);
    }

    private static void checkBatch(TimedBatch batch, long expectedChecksum) {
        if (batch.checksum != expectedChecksum) {
            throw new AssertionError("Microbenchmark checksum mismatch: expected " + expectedChecksum + ", got " + batch.checksum);
        }
    }

    private static long oracleBatch(Workload workload) {
        long result = 0;
        for (int repetition = 0; repetition < BENCH_REPETITIONS; repetition++) {
            for (int[] query : workload.queries) {
                for (int d = 0; d < workload.donors.length; d++) {
                    result += (d + 1L) * (intersection(query, workload.donors[d]) + 1L);
                }
            }
        }
        return result;
    }

    private static long indexedBatch(Workload workload, HaplotypeDonorIndex index) {
        long result = 0;
        for (int repetition = 0; repetition < BENCH_REPETITIONS; repetition++) {
            for (int[] query : workload.queries) {
                result += checksum(index.countIntersections(query), index.size());
            }
        }
        return result;
    }

    private static long checksum(int[] counts, int size) {
        long result = 0;
        for (int d = 0; d < size; d++) {
            result += (d + 1L) * (counts[d] + 1L);
        }
        return result;
    }

    private static final class TestPool {
        private final int capacity;
        private final HaplotypeDonorIndex index;
        private final List<int[]> donors = new ArrayList<>();
        private int[] lastCounts;
        private int lastSize = -1;

        private TestPool(int capacity) {
            this.capacity = capacity;
            index = new HaplotypeDonorIndex(capacity);
        }

        private void add(int[] variants) {
            int[] snapshot = variants.clone();
            index.add(variants);
            if (donors.size() < capacity) {
                donors.add(snapshot);
            }
            if (index.size() != donors.size() || !Arrays.equals(variants, snapshot)) {
                throw new AssertionError("Addition changed input or size/cap: capacity=" + capacity + ", expected size=" + donors.size());
            }
            for (int d = 0; d < donors.size(); d++) {
                if (!Arrays.equals(donors.get(d), index.get(d))) {
                    throw new AssertionError("Donor contents/order changed after addition: donor=" + d);
                }
            }
        }
    }

    private enum PoolKind {
        SPARSE_RELATED("sparse-related"), DENSE_SHARED("dense-shared"), DISJOINT("disjoint");

        private final String label;

        PoolKind(String label) {
            this.label = label;
        }
    }

    private record Selection(int donorId, int common, int score) {
    }

    private record Workload(String label, int[][] donors, int[][] queries) {
    }

    private record TimedBatch(long nanos, long checksum) {
    }
}
