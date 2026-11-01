package hypergrass;

import java.util.HashMap;

/**
 * Encoder-only donor pool indexed by anchor span. Donor IDs are plain insertion order so the
 * decoder, which keeps a flat list, resolves them without any extra signaling. Only donors whose
 * anchor span overlaps the query span are scored, so a haplotype can inherit from the homologous
 * haplotype of any earlier genome instead of an arbitrary anchor window.
 */
final class RegionalDonorIndex {

    private static final int BIN_SHIFT = 13; // 8 kb anchor bins
    private static final int MAX_CANDIDATES =
            Math.max(1, Integer.getInteger("hypergrass.maxCandidates", 12));
    private static final boolean RECENCY =
            "true".equals(System.getProperty("hypergrass.donorRecency"));

    private int[][] donors = new int[256][];
    private int[] starts = new int[256];
    private int[] ends = new int[256];
    private int size;

    private final HashMap<Integer, int[]> bins = new HashMap<>();
    private final int[] candidateIds = new int[MAX_CANDIDATES];
    private final int[] candidateOverlap = new int[MAX_CANDIDATES];

    int size() {
        return size;
    }

    int[] get(int donorId) {
        return donors[donorId];
    }

    void add(int[] variants) {
        if (size == donors.length) {
            int grown = size * 2;
            donors = java.util.Arrays.copyOf(donors, grown);
            starts = java.util.Arrays.copyOf(starts, grown);
            ends = java.util.Arrays.copyOf(ends, grown);
        }
        int lo = variants[0] >>> 2;
        int hi = (variants[variants.length - 1] >>> 2) + 1;
        int id = size++;
        donors[id] = variants;
        starts[id] = lo;
        ends[id] = hi;
        int bin = lo >>> BIN_SHIFT;
        bins.put(bin, append(bins.get(bin), id));
    }

    private static int[] append(int[] ids, int id) {
        if (ids == null) {
            return new int[] {id};
        }
        int[] grown = java.util.Arrays.copyOf(ids, ids.length + 1);
        grown[ids.length] = id;
        return grown;
    }

    /**
     * Returns donor IDs whose anchor span overlaps {@code [lo, hi)}, ordered by decreasing overlap.
     * At most {@link #MAX_CANDIDATES} entries are returned.
     */
    int[] candidates(int lo, int hi) {
        int firstBin = lo >>> BIN_SHIFT;
        int lastBin = (hi - 1) >>> BIN_SHIFT;
        int found = 0;
        for (int bin = firstBin; bin <= lastBin; bin++) {
            int[] ids = bins.get(bin);
            if (ids == null) continue;
            if (RECENCY) {
                // Newest donors first: on equal span overlap the most recently processed haplotype
                // (adjacent genome) is the better recombination template than the oldest one.
                for (int k = ids.length - 1; k >= 0; k--) {
                    found = consider(ids[k], lo, hi, found);
                }
            } else {
                for (int id : ids) {
                    found = consider(id, lo, hi, found);
                }
            }
        }
        int[] out = new int[found];
        System.arraycopy(candidateIds, 0, out, 0, found);
        return out;
    }

    private int consider(int id, int lo, int hi, int found) {
        int overlap = Math.min(hi, ends[id]) - Math.max(lo, starts[id]);
        if (overlap <= 0) return found;
        return insert(id, overlap, found);
    }

    /** Keeps the candidate table ordered by decreasing overlap; insertion sort over <= 12 slots. */
    private int insert(int id, int overlap, int found) {
        if (found == MAX_CANDIDATES && overlap <= candidateOverlap[found - 1]) {
            return found;
        }
        int at = Math.min(found, MAX_CANDIDATES - 1);
        while (at > 0 && candidateOverlap[at - 1] < overlap) {
            candidateIds[at] = candidateIds[at - 1];
            candidateOverlap[at] = candidateOverlap[at - 1];
            at--;
        }
        candidateIds[at] = id;
        candidateOverlap[at] = overlap;
        return Math.min(found + 1, MAX_CANDIDATES);
    }
}
