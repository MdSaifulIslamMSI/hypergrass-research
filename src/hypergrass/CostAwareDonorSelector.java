package hypergrass;

/** Encoder-only raw-stream byte-cost model; final compressed bytes decide whether it is retained. */
final class CostAwareDonorSelector {

    static final int BUDGET_EXHAUSTED = -2;
    private long remainingWork;

    CostAwareDonorSelector(long workBudget) {
        if (workBudget < 0) throw new IllegalArgumentException("Negative donor-cost work budget");
        remainingWork = workBudget;
    }

    /** Returns -1 for scratch, a donor index, or BUDGET_EXHAUSTED to discard this encoding trial. */
    int select(int[] variants, HaplotypeDonorIndex donors, int[] commonCounts) {
        if (!charge(variants.length)) return BUDGET_EXHAUSTED;
        long bestCost = 1L + VarIntCodec.varIntSize(variants.length) + addedBytes(variants, null);
        int bestDonor = -1;
        for (int d = 0; d < donors.size(); d++) {
            int common = commonCounts[d];
            if (common == 0) continue;
            int[] donor = donors.get(d);
            if (!charge((long) variants.length + donor.length)) return BUDGET_EXHAUSTED;
            int added = variants.length - common;
            long maskBytes = common == donor.length ? 0 : (donor.length + 7L) / 8;
            long cost = 1L + VarIntCodec.varIntSize(d) + VarIntCodec.varIntSize(added)
                    + maskBytes + addedBytes(variants, donor);
            // Scratch wins ties; otherwise preserve the first donor attaining the lowest cost.
            if (cost < bestCost) {
                bestCost = cost;
                bestDonor = d;
            }
        }
        return bestDonor;
    }

    private boolean charge(long work) {
        if (work > remainingWork) return false;
        remainingWork -= work;
        return true;
    }

    /** Gap predictor resets per added list, just as in writeVariantList; differences cost one byte. */
    private static long addedBytes(int[] variants, int[] donor) {
        int j = 0, previousPosition = -1;
        long bytes = 0;
        for (int variant : variants) {
            if (donor != null) {
                while (j < donor.length && donor[j] < variant) j++;
                if (j < donor.length && donor[j] == variant) continue;
            }
            int position = variant >>> 2;
            bytes += VarIntCodec.varIntSize(position - previousPosition - 1) + 1L;
            previousPosition = position;
        }
        return bytes;
    }
}
