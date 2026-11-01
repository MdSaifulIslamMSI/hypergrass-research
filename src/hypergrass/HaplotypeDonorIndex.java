package hypergrass;

import java.util.Arrays;

/**
 * Encoder-only inverted index from packed anchor variants to 64-donor membership masks. Bit-sliced
 * counters accumulate exact intersections for all lanes of a mask at once, including dense clades.
 * Variant arrays must remain sorted, unique and immutable. No index is serialized.
 */
final class HaplotypeDonorIndex {

    private static final int INDEX_THRESHOLD = 16;
    private static final int MAX_INDEX_VARIANTS = 131_072; // Bound per-pool index storage, not donor data.
    private static final int COUNTER_BITS = 31; // An intersection cannot exceed an int array length.

    private final int[][] donors;
    private final int[] counts;
    private int size;
    private long totalVariants;
    private int indexedDonors;

    // A zero head denotes an empty hash slot; posting links are one-based.
    private int[] keys;
    private int[] heads;
    private int distinctVariants;
    private int[] postingWords;
    private long[] membership;
    private int[] next;
    private int postings;
    private long[] counterPlanes;

    HaplotypeDonorIndex(int maxDonors) {
        donors = new int[maxDonors][];
        counts = new int[maxDonors];
    }

    int size() {
        return size;
    }

    int[] get(int donorId) {
        return donorId >= 0 && donorId < size ? donors[donorId] : EMPTY;
    }

    private static final int[] EMPTY = new int[0];

    void add(int[] variants) {
        if (size < donors.length) {
            donors[size++] = variants;
            totalVariants += variants.length;
        }
    }

    /** Returns a reused buffer of exact intersection sizes, indexed by donor insertion order. */
    int[] countIntersections(int[] variants) {
        Arrays.fill(counts, 0, size, 0);
        if (variants.length == 0) return counts;
        // Avoid an unamortized build for tiny queries and excessive index storage for large pools.
        if (size < INDEX_THRESHOLD || totalVariants > MAX_INDEX_VARIANTS
                || (keys == null && variants.length < 4)) {
            for (int d = 0; d < size; d++) {
                int[] donor = donors[d];
                int i = 0, j = 0, common = 0;
                while (i < variants.length && j < donor.length) {
                    if (variants[i] == donor[j]) {
                        common++;
                        i++;
                        j++;
                    } else if (variants[i] < donor[j]) {
                        i++;
                    } else {
                        j++;
                    }
                }
                counts[d] = common;
            }
            return counts;
        }

        if (keys == null) {
            keys = new int[64];
            heads = new int[64];
            postingWords = new int[128];
            membership = new long[128];
            next = new int[128];
            counterPlanes = new long[((donors.length + 63) >>> 6) * COUNTER_BITS];
        }
        while (indexedDonors < size) {
            indexDonor(indexedDonors++);
        }

        Arrays.fill(counterPlanes, 0L);
        int mask = keys.length - 1;
        for (int variant : variants) {
            int slot = hash(variant) & mask;
            while (heads[slot] != 0) {
                if (keys[slot] == variant) {
                    for (int p = heads[slot]; p != 0; p = next[p - 1]) {
                        long carry = membership[p - 1];
                        int plane = postingWords[p - 1] * COUNTER_BITS;
                        // Each bit is an independent donor lane; carry moves between planes.
                        while (carry != 0L) {
                            long old = counterPlanes[plane];
                            counterPlanes[plane++] = old ^ carry;
                            carry &= old;
                        }
                    }
                    break;
                }
                slot = (slot + 1) & mask;
            }
        }
        int numPlanes = 32 - Integer.numberOfLeadingZeros(variants.length);
        for (int word = 0; word < ((size + 63) >>> 6); word++) {
            for (int bit = 0; bit < numPlanes; bit++) {
                long lanes = counterPlanes[word * COUNTER_BITS + bit];
                while (lanes != 0L) {
                    int donorId = (word << 6) + Long.numberOfTrailingZeros(lanes);
                    counts[donorId] |= 1 << bit;
                    lanes &= lanes - 1;
                }
            }
        }
        return counts;
    }

    private void indexDonor(int donorId) {
        int word = donorId >>> 6;
        long lane = 1L << (donorId & 63);
        for (int variant : donors[donorId]) {
            int mask = keys.length - 1;
            int slot = hash(variant) & mask;
            while (heads[slot] != 0 && keys[slot] != variant) {
                slot = (slot + 1) & mask;
            }
            int head = heads[slot];
            // Append-only donor IDs put the newest 64-donor group at the list head.
            if (head != 0 && postingWords[head - 1] == word) {
                membership[head - 1] |= lane;
                continue;
            }
            if (postings == postingWords.length) {
                postingWords = Arrays.copyOf(postingWords, postings * 2);
                membership = Arrays.copyOf(membership, postings * 2);
                next = Arrays.copyOf(next, postings * 2);
            }
            if (head == 0) {
                if ((distinctVariants + 1) * 2 > keys.length) {
                    growHashTable();
                    mask = keys.length - 1;
                    slot = hash(variant) & mask;
                    while (heads[slot] != 0) {
                        slot = (slot + 1) & mask;
                    }
                }
                keys[slot] = variant;
                distinctVariants++;
            }
            postingWords[postings] = word;
            membership[postings] = lane;
            next[postings] = head;
            heads[slot] = ++postings;
        }
    }

    private void growHashTable() {
        int[] oldKeys = keys;
        int[] oldHeads = heads;
        keys = new int[oldKeys.length * 2];
        heads = new int[keys.length];
        int mask = keys.length - 1;
        for (int i = 0; i < oldKeys.length; i++) {
            if (oldHeads[i] == 0) continue;
            int slot = hash(oldKeys[i]) & mask;
            while (heads[slot] != 0) {
                slot = (slot + 1) & mask;
            }
            keys[slot] = oldKeys[i];
            heads[slot] = oldHeads[i];
        }
    }

    private static int hash(int variant) {
        variant ^= variant >>> 16;
        variant *= 0x7FEB352D;
        variant ^= variant >>> 15;
        variant *= 0x846CA68B;
        return variant ^ (variant >>> 16);
    }
}
