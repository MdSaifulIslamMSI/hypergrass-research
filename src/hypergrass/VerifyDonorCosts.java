package hypergrass;

import hypergrass.VarIntCodec.ByteOut;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import java.util.TreeMap;

/** Independent pre-Deflate serialization oracle for the encoder-only donor-cost selector. */
public final class VerifyDonorCosts {

    private static final int[] EMPTY = new int[0];
    private static final int MAX_POSITION = Integer.MAX_VALUE >>> 2;
    private static final long SEED = 0xD0C0_2026_7139L;
    private long checks;
    private long selections;
    private long budgetChecks;
    private long serializedCandidates;
    private long randomizedCases;
    private long scratchSelections;
    private long donorSelections;
    private int suites;

    private VerifyDonorCosts() {}

    public static void main(String[] args) {
        if (args.length != 0) {
            throw new IllegalArgumentException("Usage: java -cp bin hypergrass.VerifyDonorCosts");
        }
        VerifyDonorCosts verifier = new VerifyDonorCosts();
        verifier.rawColumns();
        verifier.scratchAndInheritance();
        verifier.positionAndGapBoundaries();
        verifier.countBoundaries();
        verifier.maskBoundaries();
        verifier.donorIdsAndTies();
        verifier.budgets();
        verifier.exhaustiveSubsets();
        verifier.randomizedDifferential();
        System.out.printf(Locale.ROOT,
                "ALL %d DONOR COST SUITES PASSED (%d checks, %d selections, %d budget checks, "
                        + "%d serialized candidates, %d randomized cases; seed=0x%x)%n",
                verifier.suites, verifier.checks, verifier.selections, verifier.budgetChecks,
                verifier.serializedCandidates, verifier.randomizedCases, SEED);
    }

    private void rawColumns() {
        RawCandidate scratch = serialize(positions(0, 128), null, -1);
        equal(scratch.bytes(), 6, "scratch independent raw cost");
        bytes(scratch.modes, 1);
        bytes(scratch.ids);
        bytes(scratch.masks);
        bytes(scratch.counts, 2);
        bytes(scratch.gaps, 0, 127);
        bytes(scratch.differences, 1, 1);

        RawCandidate full = serialize(positions(0, 128), positions(0), 128);
        equal(full.bytes(), 7, "full inheritance with two-byte donor ID");
        bytes(full.modes, 2);
        bytes(full.ids, 128, 1);
        bytes(full.masks);
        bytes(full.counts, 1);
        bytes(full.gaps, 128, 1);
        bytes(full.differences, 1);

        int[] query = {pack(0, 1), pack(127, 2), pack(128, 3)};
        int[] donor = {pack(0, 2), pack(127, 2), pack(128, 1)};
        RawCandidate partial = serialize(query, donor, 0);
        equal(partial.bytes(), 8, "partial inheritance independent raw cost");
        equal(partial.common, 1, "alleles, not just coordinates, determine inheritance");
        bytes(partial.modes, 3);
        bytes(partial.ids, 0);
        bytes(partial.masks, 2);
        bytes(partial.counts, 2);
        bytes(partial.gaps, 0, 127);
        bytes(partial.differences, 1, 3);
        pass("independent raw columns and known byte totals");
    }

    private void scratchAndInheritance() {
        expected("empty query and pool", EMPTY, new int[0][], -1);
        expected("scratch without donors", positions(0, 127, 128), new int[0][], -1);
        expected("empty donors", positions(0, 1), new int[][] {EMPTY, EMPTY}, -1);
        expected("empty query with nonempty donors", EMPTY, new int[][] {positions(0, 1)}, -1);
        expected("disjoint donors", positions(0, 1), new int[][] {positions(2, 3)}, -1);
        expected("full identical donor", positions(0, 1, 2), new int[][] {positions(0, 1, 2)}, 0);
        expected("full subset donor", positions(0, 1, 2), new int[][] {positions(0, 1)}, 0);
        expected("partial retained donor", positions(0, 1, 2), new int[][] {positions(0, 1, 8)}, 0);
        expected("same coordinates, different alleles", positions(0, 1, 2),
                new int[][] {{pack(0, 2), pack(1, 2), pack(2, 2)}}, -1);
        expected("mixed allele mismatch tie", new int[] {pack(0, 1), pack(127, 2), pack(128, 3)},
                new int[][] {{pack(0, 2), pack(127, 2), pack(128, 1)}}, -1);
        expected("later strictly better donor", positions(0, 1, 2),
                new int[][] {positions(0), positions(0, 1)}, 1);
        pass("scratch, empty, disjoint, full, partial and allele mismatches");
    }

    private void positionAndGapBoundaries() {
        for (int boundary : new int[] {128, 16_384, 2_097_152, 268_435_456}) {
            int[] query = positions(boundary - 1, boundary);
            expected("added gap predictor resets at " + boundary, query,
                    new int[][] {positions(boundary - 1)}, -1);
            RawCandidate candidate = serialize(query, positions(boundary - 1), 0);
            ByteOut independentlyExpected = new ByteOut(16);
            independentlyExpected.writeVarInt(boundary);
            check(Arrays.equals(candidate.gaps.toByteArray(), independentlyExpected.toByteArray()),
                    "added gap starts at zero, not at last inherited coordinate");
        }
        int[] high = {0, 1, 126, 127, 128, 129, 16_382, 16_383, 16_384, 16_385,
                2_097_151, 2_097_152, 268_435_455, 268_435_456, MAX_POSITION - 1, MAX_POSITION};
        for (int position : high) {
            int[] query = positions(position);
            expected("high coordinate " + position, query, new int[][] {query}, 0);
        }
        int[] query = positions(MAX_POSITION - 1, MAX_POSITION);
        expected("largest nonnegative packed coordinate", query,
                new int[][] {positions(MAX_POSITION - 1)}, 0);
        RawCandidate highAdded = serialize(query, positions(MAX_POSITION - 1), 0);
        bytes(highAdded.gaps, 255, 255, 255, 255, 1);
        for (int gap : new int[] {0, 1, 126, 127, 128, 129, 16_383, 16_384, 2_097_151, 2_097_152}) {
            int[] spaced = positions(0, gap + 1, gap + 2);
            differential("gap boundary " + gap, spaced,
                    new int[][] {positions(0), positions(gap + 1), positions(0, gap + 2)});
        }
        pass("high coordinates and independently reset gap-varint boundaries");
    }

    private void countBoundaries() {
        for (int count : new int[] {0, 1, 2, 126, 127, 128, 129, 16_383, 16_384}) {
            int[] query = range(0, count);
            RawCandidate scratch = serialize(query, null, -1);
            ByteOut countBytes = new ByteOut(16);
            countBytes.writeVarInt(count);
            check(Arrays.equals(scratch.counts.toByteArray(), countBytes.toByteArray()), "scratch count varint");
            equal(scratch.bytes(), 1L + countBytes.size() + 2L * count, "contiguous scratch byte total");
            expected("scratch count=" + count, query, new int[0][], -1);
            expected("full count=" + count, query, new int[][] {query}, count == 0 ? -1 : 0);
            if (count > 0) {
                int[] withInheritedPrefix = range(0, count + 1);
                RawCandidate inherited = serialize(withInheritedPrefix, positions(0), 0);
                check(Arrays.equals(inherited.counts.toByteArray(), countBytes.toByteArray()), "added count varint");
                expected("added count=" + count, withInheritedPrefix, new int[][] {positions(0)}, 0);
            }
        }
        bytes(serialize(range(0, 127), null, -1).counts, 127);
        bytes(serialize(range(0, 128), null, -1).counts, 128, 1);
        bytes(serialize(range(0, 16_384), null, -1).counts, 128, 128, 1);
        pass("scratch and added-count varint boundaries");
    }

    private void maskBoundaries() {
        int[] lengths = {7, 8, 9, 55, 56, 57, 63, 64, 65, 71, 72, 73};
        int[] expectedBytes = {1, 1, 2, 7, 7, 8, 8, 8, 9, 9, 9, 10};
        for (int i = 0; i < lengths.length; i++) {
            int length = lengths[i];
            int[] donor = range(0, length);
            int[] query = positions(0, length - 1, length + 1);
            RawCandidate partial = serialize(query, donor, 0);
            equal(partial.masks.size(), expectedBytes[i], "rounded retained mask, donor length=" + length);
            byte[] expectedMask = new byte[expectedBytes[i]];
            expectedMask[0] |= 1;
            expectedMask[(length - 1) / 8] |= (byte) (1 << ((length - 1) % 8));
            check(Arrays.equals(partial.masks.toByteArray(), expectedMask), "first/last retained mask bits");
            differential("partial mask length=" + length, query, new int[][] {donor});
            RawCandidate full = serialize(donor, donor, 0);
            equal(full.masks.size(), 0, "full inheritance omits mask regardless of donor length");
            expected("full donor length=" + length, donor, new int[][] {donor}, 0);
        }
        for (int length : new int[] {7, 8, 9}) {
            int[] donor = new int[length];
            donor[0] = pack(0, 1);
            System.arraycopy(range(10, length - 1), 0, donor, 1, length - 1);
            long expectedCost = length == 9 ? 7 : 6;
            equal(serialize(positions(0, 1), donor, 0).bytes(), expectedCost, "one retained variant mask cost");
            expected("mask scratch tie/loss, length=" + length, positions(0, 1), new int[][] {donor}, -1);
        }
        pass("7/8/9 donor lengths and 7/8/9-byte mask rounding");
    }

    private void donorIdsAndTies() {
        for (int id : new int[] {0, 1, 63, 64, 126, 127, 128, 129, 255, 256, 511}) {
            int[][] donors = new int[id + 1][];
            Arrays.fill(donors, EMPTY);
            donors[id] = positions(0);
            expected("donor ID=" + id, positions(0), donors, id < 128 ? id : -1);
        }
        bytes(serialize(positions(0), positions(0), 127).ids, 127);
        bytes(serialize(positions(0), positions(0), 128).ids, 128, 1);
        int[] query = positions(0, 1, 2);
        expected("identical donor tie", query, new int[][] {query, query.clone()}, 0);
        expected("distinct donor tie", query, new int[][] {positions(0), positions(1)}, 0);
        expected("distinct donor tie reversed", query, new int[][] {positions(1), positions(0)}, 0);
        for (int repeat = 0; repeat < 16; repeat++) {
            expected("repeated partial scratch tie " + repeat, positions(0, 1),
                    new int[][] {positions(0, 2)}, -1);
        }
        int[][] donors = new int[129][];
        Arrays.fill(donors, EMPTY);
        donors[127] = positions(0, 1, 100);
        donors[128] = positions(0, 1);
        equal(serialize(query, donors[127], 127).bytes(), 6, "partial one-byte ID tie cost");
        equal(serialize(query, donors[128], 128).bytes(), 6, "full two-byte ID tie cost");
        expected("first donor wins across 127/128 ID boundary", query, donors, 127);
        pass("donor-ID varints, scratch ties and deterministic first-donor ties");
    }

    private void budgets() {
        equal(CostAwareDonorSelector.BUDGET_EXHAUSTED, -2, "budget exhaustion sentinel");
        for (long negative : new long[] {-1, Long.MIN_VALUE}) {
            try {
                new CostAwareDonorSelector(negative);
                throw new AssertionError("Accepted negative budget: " + negative);
            } catch (IllegalArgumentException expected) {
                checks++;
            }
        }
        int[] query = positions(0, 1, 2);
        HaplotypeDonorIndex emptyPool = index(new int[0][]);
        CostAwareDonorSelector zero = new CostAwareDonorSelector(0);
        budgetSelect("zero work at zero budget", zero, EMPTY, emptyPool, EMPTY, -1);
        budgetSelect("nonempty scratch at zero budget", zero, query, emptyPool, EMPTY, -2);
        budgetSelect("empty query after failed charge", zero, EMPTY, emptyPool, EMPTY, -1);

        int[][] unrelated = {range(20, 1024), {pack(0, 2), pack(1, 2), pack(2, 2)}};
        Oracle noCommon = oracle(query, unrelated);
        equal(noCommon.work, query.length, "zero-common donors are not charged");
        budgetSelect("only scratch scan is charged", new CostAwareDonorSelector(query.length), query,
                index(unrelated), noCommon.commonCounts, -1);

        int[][] donors = {query.clone(), positions(0), positions(0, 1, 10), EMPTY, range(100, 500)};
        Oracle expected = oracle(query, donors);
        equal(expected.work, 3 + 6 + 4 + 6, "scratch plus every positive-common donor scan");
        HaplotypeDonorIndex pool = index(donors);
        for (long budget : new long[] {0, 2, 3, 8, 9, 12, 13, 18, 19, 20, 1L << 32, Long.MAX_VALUE}) {
            budgetSelect("budget=" + budget, new CostAwareDonorSelector(budget), query, pool,
                    expected.commonCounts, budget < expected.work ? -2 : expected.donor);
        }
        CostAwareDonorSelector cumulative = new CostAwareDonorSelector(2 * expected.work);
        budgetSelect("first exact cumulative scan", cumulative, query, pool, expected.commonCounts, expected.donor);
        budgetSelect("second exact cumulative scan", cumulative, query, pool, expected.commonCounts, expected.donor);
        budgetSelect("empty cumulative scan is free", cumulative, EMPTY, pool, new int[donors.length], -1);
        budgetSelect("third cumulative scan exhausts", cumulative, query, pool, expected.commonCounts, -2);

        int[] second = positions(0, 1, 10, 127, 128);
        Oracle secondOracle = oracle(second, donors);
        long combined = expected.work + secondOracle.work;
        CostAwareDonorSelector exact = new CostAwareDonorSelector(combined);
        budgetSelect("mixed first scan", exact, query, pool, expected.commonCounts, expected.donor);
        budgetSelect("mixed second exact scan", exact, second, pool, secondOracle.commonCounts, secondOracle.donor);
        budgetSelect("mixed budget completely consumed", exact, positions(0), emptyPool, EMPTY, -2);
        CostAwareDonorSelector shortByOne = new CostAwareDonorSelector(combined - 1);
        budgetSelect("mixed first scan short budget", shortByOne, query, pool, expected.commonCounts, expected.donor);
        budgetSelect("mixed second scan short by one", shortByOne, second, pool, secondOracle.commonCounts, -2);
        pass("negative, zero, exact, exhausted and shared cumulative budgets");
    }

    private void exhaustiveSubsets() {
        int[] universe = {pack(0, 1), pack(1, 2), pack(126, 3), pack(127, 1),
                pack(128, 2), pack(16_383, 3), pack(16_384, 1), pack(MAX_POSITION, 3)};
        int[][] donors = new int[1 << universe.length][];
        for (int mask = 0; mask < donors.length; mask++) {
            int[] subset = new int[Integer.bitCount(mask)];
            int next = 0;
            for (int bit = 0; bit < universe.length; bit++) {
                if ((mask & (1 << bit)) != 0) subset[next++] = universe[bit];
            }
            donors[mask] = subset;
        }
        for (int mask = 0; mask < donors.length; mask++) {
            differential("exhaustive subset=" + mask, donors[mask], donors);
            differential("exhaustive changed alleles=" + mask, changedAlleles(donors[mask]), donors);
        }
        pass("all 256 boundary-coordinate subsets and allele-altered queries");
    }

    private void randomizedDifferential() {
        Random random = new Random(SEED);
        int[] poolSizes = {0, 1, 2, 7, 8, 9, 15, 16, 17, 63, 64, 65, 127, 128, 129, 255, 256, 512};
        long scratchBefore = scratchSelections, donorsBefore = donorSelections;
        for (int trial = 0; trial < 240; trial++) {
            int size = poolSizes[trial % poolSizes.length];
            int[] base = randomVariants(random, 8 + random.nextInt(153));
            int[][] donors = new int[size][];
            for (int d = 0; d < size; d++) {
                donors[d] = switch (d % 6) {
                    case 0 -> base.clone();
                    case 1 -> EMPTY;
                    case 2 -> changedAlleles(base);
                    case 3 -> randomVariants(random, random.nextInt(150));
                    default -> relatedVariants(random, base);
                };
            }
            int[][] queries = {EMPTY, base, changedAlleles(base), relatedVariants(random, base),
                    randomVariants(random, random.nextInt(160)), Arrays.copyOf(base, 1),
                    Arrays.copyOfRange(base, base.length / 2, base.length),
                    size == 0 ? base : donors[random.nextInt(size)],
                    randomVariants(random, 127), randomVariants(random, 128),
                    relatedVariants(random, size == 0 ? base : donors[random.nextInt(size)]),
                    randomVariants(random, 129)};
            Oracle[] results = new Oracle[queries.length];
            long cumulativeWork = 0;
            for (int q = 0; q < queries.length; q++) {
                results[q] = differential("random trial=" + trial + ", query=" + q, queries[q], donors);
                cumulativeWork += results[q].work;
                randomizedCases++;
            }
            HaplotypeDonorIndex pool = index(donors);
            CostAwareDonorSelector exact = new CostAwareDonorSelector(cumulativeWork);
            CostAwareDonorSelector shortByOne = new CostAwareDonorSelector(cumulativeWork - 1);
            for (int q = 0; q < queries.length; q++) {
                Oracle result = results[q];
                budgetSelect("random cumulative exact", exact, queries[q], pool, result.commonCounts, result.donor);
                budgetSelect("random cumulative short", shortByOne, queries[q], pool, result.commonCounts,
                        q + 1 == queries.length ? -2 : result.donor);
            }
            budgetSelect("random cumulative consumed", exact, positions(0), index(new int[0][]), EMPTY, -2);
        }
        check(scratchSelections > scratchBefore, "randomized oracle selects scratch");
        check(donorSelections > donorsBefore, "randomized oracle selects donors");
        pass("seeded randomized differential selection and cumulative budgets");
    }

    private void expected(String label, int[] query, int[][] donors, int expected) {
        Oracle result = differential(label, query, donors);
        equal(result.donor, expected, label + ": independent expected winner");
    }

    private Oracle differential(String label, int[] query, int[][] donors) {
        sorted(query);
        for (int[] donor : donors) sorted(donor);
        int[] originalQuery = query.clone();
        int[][] originalDonors = new int[donors.length][];
        for (int d = 0; d < donors.length; d++) originalDonors[d] = donors[d].clone();
        Oracle expected = oracle(query, donors);
        int[] originalCounts = expected.commonCounts.clone();
        HaplotypeDonorIndex pool = index(donors);
        CostAwareDonorSelector repeated = new CostAwareDonorSelector(Long.MAX_VALUE);
        select(label, repeated, query, pool, expected.commonCounts, expected.donor);
        select(label + " repeat", repeated, query, pool, expected.commonCounts, expected.donor);
        budgetSelect(label + " exact", new CostAwareDonorSelector(expected.work), query, pool,
                expected.commonCounts, expected.donor);
        if (expected.work > 0) {
            budgetSelect(label + " short by one", new CostAwareDonorSelector(expected.work - 1), query, pool,
                    expected.commonCounts, -2);
        }
        check(Arrays.equals(query, originalQuery), label + ": immutable query");
        check(Arrays.equals(expected.commonCounts, originalCounts), label + ": immutable common counts");
        for (int d = 0; d < donors.length; d++) {
            check(Arrays.equals(donors[d], originalDonors[d]), label + ": immutable donor=" + d);
        }
        if (expected.donor < 0) scratchSelections++;
        else donorSelections++;
        return expected;
    }

    private Oracle oracle(int[] query, int[][] donors) {
        long best = serialize(query, null, -1).bytes();
        long work = query.length;
        int winner = -1;
        int[] commonCounts = new int[donors.length];
        for (int d = 0; d < donors.length; d++) {
            RawCandidate candidate = serialize(query, donors[d], d);
            commonCounts[d] = candidate.common;
            if (candidate.common > 0) work += (long) query.length + donors[d].length;
            long cost = candidate.bytes();
            if (cost < best) {
                best = cost;
                winner = d;
            }
        }
        return new Oracle(winner, work, commonCounts);
    }

    /** Serialize candidate columns, not a cost formula: no varIntSize or compression shortcuts. */
    private RawCandidate serialize(int[] query, int[] donor, int donorId) {
        ByteOut modes = new ByteOut(16), ids = new ByteOut(16), masks = new ByteOut(16);
        ByteOut counts = new ByteOut(16), gaps = new ByteOut(16), differences = new ByteOut(16);
        int[] added = new int[query.length];
        int addedCount = 0;
        for (int variant : query) {
            if (donor == null || Arrays.binarySearch(donor, variant) < 0) added[addedCount++] = variant;
        }
        int common = query.length - addedCount;
        if (donor == null) {
            modes.writeByte(1);
        } else {
            boolean full = common == donor.length;
            modes.writeByte(full ? 2 : 3);
            ids.writeVarInt(donorId);
            if (!full) {
                // Build each mask byte from donor membership, independently of the added-list scan.
                for (int first = 0; first < donor.length; first += 8) {
                    int bits = 0;
                    for (int bit = 0; bit < 8 && first + bit < donor.length; bit++) {
                        if (Arrays.binarySearch(query, donor[first + bit]) >= 0) bits |= 1 << bit;
                    }
                    masks.writeByte(bits);
                }
            }
        }
        counts.writeVarInt(addedCount);
        int previousPosition = -1;
        for (int i = 0; i < addedCount; i++) {
            int position = added[i] >>> 2;
            gaps.writeVarInt(position - previousPosition - 1);
            differences.writeByte(added[i] & 3);
            previousPosition = position;
        }
        serializedCandidates++;
        return new RawCandidate(modes, ids, masks, counts, gaps, differences, common);
    }

    private record RawCandidate(ByteOut modes, ByteOut ids, ByteOut masks, ByteOut counts,
                                ByteOut gaps, ByteOut differences, int common) {
        long bytes() {
            return (long) modes.size() + ids.size() + masks.size() + counts.size()
                    + gaps.size() + differences.size();
        }
    }

    private record Oracle(int donor, long work, int[] commonCounts) {}

    private static HaplotypeDonorIndex index(int[][] donors) {
        HaplotypeDonorIndex index = new HaplotypeDonorIndex(donors.length);
        for (int[] donor : donors) index.add(donor);
        return index;
    }

    private void budgetSelect(String label, CostAwareDonorSelector selector, int[] query,
                              HaplotypeDonorIndex donors, int[] commonCounts, int expected) {
        budgetChecks++;
        select(label, selector, query, donors, commonCounts, expected);
    }

    private void select(String label, CostAwareDonorSelector selector, int[] query,
                        HaplotypeDonorIndex donors, int[] commonCounts, int expected) {
        int actual = selector.select(query, donors, commonCounts);
        selections++;
        equal(actual, expected, label);
    }

    private void sorted(int[] variants) {
        int previous = -1;
        for (int variant : variants) {
            int position = variant >>> 2;
            check(variant >= 0 && (variant & 3) != 0 && position > previous,
                    "fixture must have sorted, unique coordinates and nonzero alleles");
            previous = position;
        }
    }

    private static int[] positions(int... positions) {
        int[] variants = new int[positions.length];
        for (int i = 0; i < positions.length; i++) variants[i] = pack(positions[i], 1);
        return variants;
    }

    private static int[] range(int start, int length) {
        int[] variants = new int[length];
        for (int i = 0; i < length; i++) variants[i] = pack(start + i, 1);
        return variants;
    }

    private static int pack(int position, int allele) {
        if (position < 0 || position > MAX_POSITION || allele < 1 || allele > 3) {
            throw new IllegalArgumentException("Invalid test position/allele");
        }
        return (position << 2) | allele;
    }

    private static int[] changedAlleles(int[] variants) {
        int[] changed = variants.clone();
        for (int i = 0; i < changed.length; i++) changed[i] = pack(variants[i] >>> 2, variants[i] % 4 % 3 + 1);
        return changed;
    }

    private static int[] randomVariants(Random random, int length) {
        TreeMap<Integer, Integer> byPosition = new TreeMap<>();
        while (byPosition.size() < length) byPosition.put(randomPosition(random), 1 + random.nextInt(3));
        return packed(byPosition);
    }

    private static int[] relatedVariants(Random random, int[] base) {
        TreeMap<Integer, Integer> byPosition = new TreeMap<>();
        for (int variant : base) {
            if (random.nextInt(4) != 0) {
                int allele = (variant & 3);
                if (random.nextInt(5) == 0) allele = allele % 3 + 1;
                byPosition.put(variant >>> 2, allele);
            }
        }
        int extras = random.nextInt(24);
        for (int i = 0; i < extras; i++) byPosition.put(randomPosition(random), 1 + random.nextInt(3));
        return packed(byPosition);
    }

    private static int[] packed(TreeMap<Integer, Integer> byPosition) {
        int[] variants = new int[byPosition.size()];
        int next = 0;
        for (var entry : byPosition.entrySet()) variants[next++] = pack(entry.getKey(), entry.getValue());
        return variants;
    }

    private static int randomPosition(Random random) {
        int[] boundaries = {0, 1, 126, 127, 128, 129, 16_383, 16_384, 2_097_151,
                2_097_152, 268_435_455, 268_435_456, MAX_POSITION};
        return switch (random.nextInt(4)) {
            case 0 -> random.nextInt(512);
            case 1 -> boundaries[random.nextInt(boundaries.length)];
            case 2 -> MAX_POSITION - random.nextInt(512);
            default -> random.nextInt(MAX_POSITION + 1);
        };
    }

    private void bytes(ByteOut actual, int... expected) {
        byte[] bytes = new byte[expected.length];
        for (int i = 0; i < expected.length; i++) bytes[i] = (byte) expected[i];
        check(Arrays.equals(actual.toByteArray(), bytes), "independent expected raw column bytes");
    }

    private void equal(long actual, long expected, String label) {
        check(actual == expected, label + ": expected " + expected + ", got " + actual);
    }

    private void check(boolean condition, String label) {
        checks++;
        if (!condition) throw new AssertionError(label);
    }

    private void pass(String name) {
        suites++;
        System.out.println("PASS " + name);
    }
}
