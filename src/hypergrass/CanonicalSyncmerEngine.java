package hypergrass;

import hypergrass.FastaStreamFactorizer.FactorizedChunk;
import hypergrass.FastaStreamFactorizer.RecordSpan;
import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Next-Generation Bidirectional Canonical K-mer + Anchor-Coordinate Haplotype Matrix &amp;
 * Virtual Clade-Ancestor Inference Engine (BC-LZG v2 + AC-HMC + VCAI).
 *
 * <p>Combines four orthogonal genomic compression breakthroughs:
 * <ul>
 *   <li><b>16-mer Dense Bidirectional Seed Index (step=4):</b> Guarantees hitting any &gt;= 19-bp
 *       exact homologous island on Watson (+) or Crick (-) strand, bridging across clustered SNPs
 *       down to ~94% ANI lineages.</li>
 *   <li><b>Circular-Genome Seam Predictor:</b> Interleaved circular boundary geometry prediction
 *       collapses 3 of 4 offset/length coordinates per circularly permuted genome to exact zero.</li>
 *   <li><b>Anchor-Coordinate Haplotype Matrix Coding (AC-HMC, Mode 3):</b> Projects substitutions
 *       from arbitrary circular offsets and (+/-) strands onto canonical sorted Anchor Loci
 *       {@code (anchorPos << 2) | xorDiff}, making clade haplotypes invariant to strand, circular
 *       rotation, and interspersed private mutations.</li>
 *   <li><b>Virtual Clade-Ancestor Inference (VCAI) + 1-Hop Global Clade Catalog:</b> Dynamically
 *       intersects sibling variant sets {@code V_a ∩ V_b} on both encoder and decoder to synthesize
 *       pure ancestral clade haplotypes (inheritable in Mode 2 with zero bitmask overhead) while
 *       preserving strict 1-hop O(1) block random access.</li>
 * </ul>
 */
public final class CanonicalSyncmerEngine {

    public static final int K = 16;
    private static final long KMER_MASK = (1L << (K << 1)) - 1L; // 32-bit exact 16-mer mask
    private static final int SEED_STEP = 4; // Any >= 19-bp exact island is guaranteed to hit a seed!
    private static final int MIN_MATCH_LEN = 22;
    private static final int NOVEL_SEGMENT_MIN = 64;
    public static final int ANCHOR_PAGE_SIZE = 65_536;
    public static final int DEFAULT_MAX_ANCHOR_BASES = 4 * 1024 * 1024; // 4 Mbp

    private static final int MAX_DONOR_POOL = 512;
    private static final int MAX_GLOBAL_DONORS = 96;
    private static final int[] EMPTY_INTS = new int[0];
    public static final int[][] EMPTY_DONORS = new int[0][];

    private CanonicalSyncmerEngine() {}

    /** Immutable indexed AnchorBank + 1-hop Global Clade-Haplotype Catalog. */
    public static final class AnchorBank {
        public final byte[] bases;
        public final int length;
        public final int primarySeedLen;
        public final int[][] globalDonors;
        final SeedHashTable index;

        public AnchorBank(byte[] bases, int length, SeedHashTable index) {
            this(bases, length, length, EMPTY_DONORS, index);
        }

        public AnchorBank(
                byte[] bases,
                int length,
                int primarySeedLen,
                int[][] globalDonors,
                SeedHashTable index) {
            this.bases = bases;
            this.length = length;
            this.primarySeedLen = primarySeedLen > 0 ? primarySeedLen : length;
            this.globalDonors = (globalDonors != null) ? globalDonors : EMPTY_DONORS;
            this.index = index;
        }

        /** Constructor used during decompression (no seed hash index needed). */
        public AnchorBank(byte[] bases, int length) {
            this(bases, length, length, EMPTY_DONORS, null);
        }

        public AnchorBank(byte[] bases, int length, int primarySeedLen, int[][] globalDonors) {
            this(bases, length, primarySeedLen, globalDonors, null);
        }
    }

    /** Flat open-addressed primitive hash table mapping 32-bit 16-mers to (strandBit << 30) | pos. */
    static final class SeedHashTable {
        private long[] keys;
        private int[] vals;
        private int mask;
        private int size;

        SeedHashTable(int initialCapPowerOf2) {
            int cap = 1;
            while (cap < initialCapPowerOf2) cap <<= 1;
            this.keys = new long[cap];
            this.vals = new int[cap];
            Arrays.fill(this.keys, -1L);
            this.mask = cap - 1;
            this.size = 0;
        }

        void insert(long kmer32, int strand, int pos) {
            if ((size + 1) * 2 >= keys.length) {
                rehash();
            }
            int slot = mix32(kmer32) & mask;
            int probes = 0;
            while (keys[slot] != -1L) {
                if (keys[slot] == kmer32) {
                    probes++;
                    if (probes >= 6) return;
                }
                slot = (slot + 1) & mask;
            }
            keys[slot] = kmer32;
            vals[slot] = (strand << 30) | pos;
            size++;
        }

        int firstSlot(long kmer32) {
            return mix32(kmer32) & mask;
        }

        long keyAt(int slot) {
            return keys[slot];
        }

        int valAt(int slot) {
            return vals[slot];
        }

        int nextSlot(int slot) {
            return (slot + 1) & mask;
        }

        private void rehash() {
            long[] oldK = keys;
            int[] oldV = vals;
            int newCap = oldK.length << 1;
            keys = new long[newCap];
            vals = new int[newCap];
            Arrays.fill(keys, -1L);
            mask = newCap - 1;
            size = 0;
            for (int i = 0; i < oldK.length; i++) {
                long k = oldK[i];
                if (k != -1L) {
                    int v = oldV[i];
                    insert(k, (v >>> 30) & 1, v & 0x3FFFFFFF);
                }
            }
        }

        private static int mix32(long x) {
            x ^= (x >>> 17);
            x *= 0xED5AD4BBAC4C2B51L;
            x ^= (x >>> 27);
            return (int) x;
        }
    }

    /**
     * Builds an incremental deduplicated AnchorBank across chunks AND harvests 1-hop Global
     * Virtual Clade Ancestor haplotypes when multi-record homologous clades are present.
     */
    public static AnchorBank buildAnchorBank(FactorizedChunk[] chunks, int maxAnchorBases) {
        int totalDna = 0;
        int totalRecords = 0;
        for (FactorizedChunk c : chunks) {
            totalDna += c.dnaLen;
            totalRecords += Math.max(1, c.numHeaders);
        }
        int capLimit = Math.min(maxAnchorBases, totalDna);
        if (capLimit < MIN_MATCH_LEN) {
            return new AnchorBank(new byte[0], 0, 0, EMPTY_DONORS, new SeedHashTable(64));
        }

        byte[] bank = new byte[Math.min(capLimit, 65536)];
        int bankLen = 0;
        int primarySeedLen = 0;
        SeedHashTable index = new SeedHashTable(32768);

        // Determine DNA length of the very first FASTA record so its complete sequence (including
        // short <256-bp 3' tail) forms a clean circular-ready primary archetype!
        int firstRecordTargetLen = 0;
        for (FactorizedChunk c : chunks) {
            if (c.dnaLen > 0) {
                RecordSpan[] spans = c.getRecordSpans();
                if (spans.length > 1) {
                    // Approximate DNA length of record 0 (exact when no N/exceptions, bounded otherwise)
                    firstRecordTargetLen = Math.min(c.dnaLen, spans[0].seqCharEnd - spans[0].seqCharStart);
                } else {
                    firstRecordTargetLen = c.dnaLen;
                }
                break;
            }
        }

        boolean completedFirstRecord = false;

        for (FactorizedChunk chunk : chunks) {
            if (bankLen >= maxAnchorBases) break;
            byte[] seq = chunk.dnaBases;
            int n = chunk.dnaLen;
            if (n == 0) continue;

            int pos = 0;
            int unalignedStart = 0;
            if (bankLen == 0) {
                int seedLen = Math.min(n, Math.min(maxAnchorBases, 4096));
                bank = ensureCapacity(bank, seedLen);
                System.arraycopy(seq, 0, bank, 0, seedLen);
                indexRange(bank, 0, seedLen, index);
                bankLen = seedLen;
                pos = seedLen;
                unalignedStart = seedLen;
                if (seedLen == n) {
                    primarySeedLen = bankLen;
                    completedFirstRecord = true;
                    continue;
                }
            }

            long rollingFwd = 0L;
            int validBases = 0;

            while (pos < n && bankLen < maxAnchorBases) {
                if (!completedFirstRecord && firstRecordTargetLen > 0 && pos >= firstRecordTargetLen) {
                    // Commit any remaining tail of the primary archetype record unconditionally
                    int rem = pos - unalignedStart;
                    if (rem > 0 && bankLen + rem <= maxAnchorBases) {
                        bank = ensureCapacity(bank, bankLen + rem);
                        System.arraycopy(seq, unalignedStart, bank, bankLen, rem);
                        indexRange(bank, Math.max(0, bankLen - K), bankLen + rem, index);
                        bankLen += rem;
                        unalignedStart = pos;
                    }
                    primarySeedLen = bankLen;
                    completedFirstRecord = true;
                }

                rollingFwd = ((rollingFwd << 2) | (seq[pos] & 3)) & KMER_MASK;
                validBases++;
                if (validBases >= K) {
                    int queryStart = pos - K + 1;
                    MatchResult m = findBestMatch(
                            seq, n, queryStart, unalignedStart, rollingFwd, bank, bankLen, primarySeedLen, index);
                    if (m != null && m.matchLen >= MIN_MATCH_LEN) {
                        int novelLen = m.targetStart - unalignedStart;
                        int minNovel = (!completedFirstRecord) ? 1 : NOVEL_SEGMENT_MIN;
                        if (novelLen >= minNovel && bankLen < maxAnchorBases) {
                            int add = Math.min(novelLen, maxAnchorBases - bankLen);
                            bank = ensureCapacity(bank, bankLen + add);
                            System.arraycopy(seq, unalignedStart, bank, bankLen, add);
                            int idxFrom = Math.max(0, bankLen - K);
                            indexRange(bank, idxFrom, bankLen + add, index);
                            bankLen += add;
                            MatchResult m2 = findBestMatch(
                                    seq, n, queryStart, m.targetStart, rollingFwd, bank, bankLen, primarySeedLen, index);
                            if (m2 != null && m2.matchLen > m.matchLen) {
                                m = m2;
                            }
                        }
                        pos = m.targetStart + m.matchLen;
                        unalignedStart = pos;
                        rollingFwd = 0L;
                        validBases = 0;
                        continue;
                    }
                }
                pos++;
                if (pos - unalignedStart >= 4096 && bankLen < maxAnchorBases) {
                    int safeCommit = (pos - unalignedStart) - K;
                    int add = Math.min(safeCommit, maxAnchorBases - bankLen);
                    if (add >= NOVEL_SEGMENT_MIN) {
                        bank = ensureCapacity(bank, bankLen + add);
                        System.arraycopy(seq, unalignedStart, bank, bankLen, add);
                        int idxFrom = Math.max(0, bankLen - K);
                        indexRange(bank, idxFrom, bankLen + add, index);
                        bankLen += add;
                        unalignedStart += add;
                    }
                }
            }

            int tailNovel = n - unalignedStart;
            int minTail = (!completedFirstRecord) ? 1 : NOVEL_SEGMENT_MIN;
            if (tailNovel >= minTail && bankLen < maxAnchorBases) {
                int add = Math.min(tailNovel, maxAnchorBases - bankLen);
                bank = ensureCapacity(bank, bankLen + add);
                System.arraycopy(seq, unalignedStart, bank, bankLen, add);
                int idxFrom = Math.max(0, bankLen - K);
                indexRange(bank, idxFrom, bankLen + add, index);
                bankLen += add;
            }
            if (!completedFirstRecord) {
                primarySeedLen = bankLen;
                completedFirstRecord = true;
            }
        }

        if (primarySeedLen <= 0) {
            primarySeedLen = bankLen;
        }
        byte[] trimmed = (bankLen == bank.length) ? bank : Arrays.copyOf(bank, bankLen);

        // Polish anchor bank to consensus majority allele across all chunks if multi-record cohort
        if (bankLen >= MIN_MATCH_LEN && totalRecords >= 2 && chunks.length > 1) {
            AnchorBank rawBank = new AnchorBank(trimmed, bankLen, primarySeedLen, EMPTY_DONORS, index);
            AnchorBank polished = consensusAnchor(rawBank, chunks);
            trimmed = polished.bases;
            index = polished.index;
        }

        // Harvest 1-hop Global Virtual Clade Ancestor haplotypes when multi-record cohorts exist
        int[][] globalDonors = EMPTY_DONORS;
        if (bankLen >= MIN_MATCH_LEN && totalRecords >= 6) {
            globalDonors = harvestGlobalCladeDonors(chunks, trimmed, bankLen, primarySeedLen, index);
        }

        return new AnchorBank(trimmed, bankLen, primarySeedLen, globalDonors, index);
    }

    /**
     * Returns a copy of {@code bank} whose bases are the per-position majority base of every chunk
     * aligned to it. This is an encoder-side transform only: the polished bases are what the archive
     * serializes, so decoders need no changes. Positions covered fewer than twice keep the original.
     */
    public static AnchorBank consensusAnchor(AnchorBank bank, FactorizedChunk[] chunks) {
        byte[] counts = new byte[bank.length << 2];
        for (FactorizedChunk c : chunks) {
            accumulateAnchorVotes(bank, c.dnaBases, c.dnaLen, counts);
        }
        return polishAnchor(bank, counts);
    }

    /**
     * Accumulates per-position base votes for one already-factorized chunk.
     */
    public static void accumulateAnchorVotes(
            AnchorBank bank, byte[] dna, int dnaLen, byte[] counts) {
        if (dnaLen < MIN_MATCH_LEN) return;
        ParsedAlignment aln = alignTargetToAnchor(
                dna, dnaLen, bank.bases, bank.length, bank.primarySeedLen, bank.index, false);
        if (aln == null || aln.numMatches == 0) return;
        int tPos = 0;
        for (int m = 0; m < aln.numMatches; m++) {
            tPos += aln.litLens[m];
            int aStart = aln.refOffsets[m];
            int mLen = aln.matchLens[m];
            if (aln.strands[m] == 0) {
                for (int d = 0; d < mLen; d++) {
                    bumpVote(counts, aStart + d, dna[tPos + d] & 3);
                }
            } else {
                for (int d = 0; d < mLen; d++) {
                    bumpVote(counts, aStart + mLen - 1 - d, (dna[tPos + d] ^ 3) & 3);
                }
            }
            tPos += mLen;
        }
    }

    private static void bumpVote(byte[] counts, int pos, int base) {
        if (pos < 0) return;
        int at = (pos << 2) | (base & 3);
        if (at >= counts.length) return;
        int cur = counts[at] & 0xFF;
        if (cur < 255) counts[at] = (byte) (cur + 1);
    }

    /** Applies accumulated votes and rebuilds the seed index. */
    public static AnchorBank polishAnchor(AnchorBank bank, byte[] counts) {
        int len = bank.length;
        byte[] bases = bank.bases.clone();
        int changed = 0;
        for (int p = 0; p < len; p++) {
            int at = p << 2;
            int cA = counts[at] & 0xFF;
            int cC = counts[at + 1] & 0xFF;
            int cG = counts[at + 2] & 0xFF;
            int cT = counts[at + 3] & 0xFF;
            int total = cA + cC + cG + cT;
            if (total < 2) continue;
            int best = 0, bestCount = cA;
            if (cC > bestCount) { best = 1; bestCount = cC; }
            if (cG > bestCount) { best = 2; bestCount = cG; }
            if (cT > bestCount) { best = 3; bestCount = cT; }
            int current = bases[p] & 3;
            if (bestCount > counts[at + current] && best != current) {
                bases[p] = (byte) best;
                changed++;
            }
        }
        if (changed == 0) {
            return bank;
        }
        int seedSlots = (int) Math.min(Integer.MAX_VALUE, (long) (len / SEED_STEP) * 3 + 1024);
        SeedHashTable index = new SeedHashTable(Math.max(1024, seedSlots));
        indexRange(bases, 0, len, index);
        return new AnchorBank(bases, len, bank.primarySeedLen, bank.globalDonors, index);
    }

    /**
     * Samples sweeps across chunks and applies pairwise Virtual Clade-Ancestor Inference (VCAI)
     * {@code I = V_a ∩ V_b} to extract shared ancestral clade SNP haplotypes for 1-hop global sharing.
     */
    private static int[][] harvestGlobalCladeDonors(
            FactorizedChunk[] chunks, byte[] anchor, int anchorLen, int primarySeedLen, SeedHashTable index) {
        List<int[]> sampleSweeps = new ArrayList<>(256);
        int chunkStep = Math.max(1, chunks.length / 16);
        for (int ci = 0; ci < chunks.length && sampleSweeps.size() < 240; ci += chunkStep) {
            FactorizedChunk c = chunks[ci];
            if (c.dnaLen < MIN_MATCH_LEN) continue;
            ParsedAlignment aln = alignTargetToAnchor(c.dnaBases, c.dnaLen, anchor, anchorLen, primarySeedLen, index, false);
            if (aln == null || aln.numMatches == 0 || aln.totalSnps == 0) continue;
            int[][] sweeps = buildSweepVariantSets(aln);
            for (int[] v : sweeps) {
                if (v.length >= 2 && v.length <= 8192) {
                    sampleSweeps.add(v);
                    if (sampleSweeps.size() >= 240) break;
                }
            }
        }

        if (sampleSweeps.size() < 2) {
            return EMPTY_DONORS;
        }

        List<int[]> global = new ArrayList<>();
        HaplotypeDonorIndex sampleIndex = new HaplotypeDonorIndex(sampleSweeps.size());
        sampleIndex.add(sampleSweeps.get(0));
        for (int i = 1; i < sampleSweeps.size() && global.size() < MAX_GLOBAL_DONORS; i++) {
            int[] vi = sampleSweeps.get(i);
            int[] commonCounts = sampleIndex.countIntersections(vi);
            int bestJ = -1;
            int bestCommon = 1; // require >= 2 shared SNPs
            for (int j = 0; j < i; j++) {
                int c = commonCounts[j];
                if (c > bestCommon) {
                    bestCommon = c;
                    bestJ = j;
                }
            }
            if (bestJ >= 0) {
                int[] inter = intersectSorted(vi, sampleSweeps.get(bestJ), bestCommon);
                if (!containsNearIdenticalDonor(global, inter)) {
                    global.add(inter);
                }
            }
            sampleIndex.add(vi);
        }
        return global.isEmpty() ? EMPTY_DONORS : global.toArray(new int[0][]);
    }

    private static boolean containsNearIdenticalDonor(List<int[]> list, int[] cand) {
        for (int[] d : list) {
            if (d.length == cand.length && Arrays.equals(d, cand)) {
                return true;
            }
            int common = countSortedIntersection(d, cand);
            int maxLen = Math.max(d.length, cand.length);
            if (maxLen > 0 && common * 10 >= maxLen * 9) {
                return true;
            }
        }
        return false;
    }

    /** Serializes global clade donors compactly into the .hgs AnchorBank header. */
    public static void writeGlobalDonors(ByteOut out, int primarySeedLen, int[][] globalDonors) {
        out.writeVarInt(primarySeedLen);
        int n = (globalDonors == null) ? 0 : globalDonors.length;
        out.writeVarInt(n);
        if (n == 0) return;

        ByteOut gapBuf = new ByteOut(n * 16);
        ByteOut diffBuf = new ByteOut(n * 8);
        for (int i = 0; i < n; i++) {
            int[] v = globalDonors[i];
            gapBuf.writeVarInt(v.length);
            int prevPos = -1;
            for (int x : v) {
                int pos = x >>> 2;
                int diff = x & 3;
                gapBuf.writeVarInt(pos - prevPos - 1);
                diffBuf.writeByte(diff);
                prevPos = pos;
            }
        }
        VarIntCodec.writeCompressedAux(out, gapBuf);
        VarIntCodec.writeCompressedAux(out, diffBuf);
    }

    /** Deserializes primarySeedLen + global clade donors from the .hgs AnchorBank header. */
    public static DecodedAnchorCatalog readGlobalDonors(ByteIn in) {
        int primarySeedLen = in.readVarInt();
        int n = in.readVarInt();
        if (n == 0) {
            return new DecodedAnchorCatalog(primarySeedLen, EMPTY_DONORS);
        }
        ByteIn gapIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        byte[] diffs = VarIntCodec.readCompressedAux(in);
        int dIdx = 0;
        int[][] donors = new int[n][];
        for (int i = 0; i < n; i++) {
            int len = gapIn.readVarInt();
            int[] v = new int[len];
            int prevPos = -1;
            for (int k = 0; k < len; k++) {
                int pos = prevPos + 1 + gapIn.readVarInt();
                int diff = diffs[dIdx++] & 3;
                v[k] = (pos << 2) | diff;
                prevPos = pos;
            }
            donors[i] = v;
        }
        return new DecodedAnchorCatalog(primarySeedLen, donors);
    }

    public static final class DecodedAnchorCatalog {
        public final int primarySeedLen;
        public final int[][] globalDonors;

        public DecodedAnchorCatalog(int primarySeedLen, int[][] globalDonors) {
            this.primarySeedLen = primarySeedLen;
            this.globalDonors = globalDonors;
        }
    }

    private static byte[] ensureCapacity(byte[] arr, int needed) {
        if (needed <= arr.length) return arr;
        int next = Math.max(arr.length << 1, needed);
        return Arrays.copyOf(arr, next);
    }

    private static void indexRange(byte[] bank, int startPos, int endPos, SeedHashTable index) {
        int firstSeed = ((startPos + SEED_STEP - 1) / SEED_STEP) * SEED_STEP;
        for (int p = firstSeed; p + K <= endPos; p += SEED_STEP) {
            long fwd = 0L;
            long rev = 0L;
            for (int k = 0; k < K; k++) {
                int b = bank[p + k] & 3;
                fwd = (fwd << 2) | b;
            }
            for (int k = K - 1; k >= 0; k--) {
                int b = (bank[p + k] ^ 3) & 3;
                rev = (rev << 2) | b;
            }
            index.insert(fwd, 0, p);
            index.insert(rev, 1, p);
        }
    }

    private static final class MatchResult {
        final int strand;
        final int anchorStart;
        final int targetStart;
        final int matchLen;

        MatchResult(int strand, int anchorStart, int targetStart, int matchLen) {
            this.strand = strand;
            this.anchorStart = anchorStart;
            this.targetStart = targetStart;
            this.matchLen = matchLen;
        }
    }

    private static MatchResult findBestMatch(
            byte[] target,
            int targetLen,
            int queryStart,
            int litFloor,
            long kmer32,
            byte[] anchor,
            int anchorLen,
            int primarySeedLen,
            SeedHashTable index) {
        int slot = index.firstSlot(kmer32);
        int bestStrand = 0;
        int bestAStart = 0;
        int bestTStart = 0;
        int bestLen = 0;
        int probes = 0;

        while (probes < 14) {
            long key = index.keyAt(slot);
            if (key == -1L) break;
            if (key == kmer32) {
                int val = index.valAt(slot);
                int strand = (val >>> 30) & 1;
                int seedPos = val & 0x3FFFFFFF;
                if (seedPos + K <= anchorLen && verifyExactSeed(target, queryStart, anchor, seedPos, strand)) {
                    MatchResult cand = extendWithSnpBridging(
                            target, targetLen, queryStart, litFloor, anchor, anchorLen, seedPos, strand);
                    if (cand.matchLen > bestLen) {
                        bestLen = cand.matchLen;
                        bestStrand = cand.strand;
                        bestAStart = cand.anchorStart;
                        bestTStart = cand.targetStart;
                        if (bestLen >= 8192 || bestTStart + bestLen == targetLen) {
                            break;
                        }
                    }
                }
            }
            slot = index.nextSlot(slot);
            probes++;
        }

        if (bestLen < MIN_MATCH_LEN) {
            return null;
        }
        return new MatchResult(bestStrand, bestAStart, bestTStart, bestLen);
    }

    private static boolean verifyExactSeed(byte[] target, int tPos, byte[] anchor, int aPos, int strand) {
        if (strand == 0) {
            for (int i = 0; i < K; i++) {
                if (target[tPos + i] != anchor[aPos + i]) return false;
            }
            return true;
        } else {
            for (int i = 0; i < K; i++) {
                if (target[tPos + i] != (anchor[aPos + K - 1 - i] ^ 3)) return false;
            }
            return true;
        }
    }

    private static MatchResult extendWithSnpBridging(
            byte[] target,
            int targetLen,
            int seedTStart,
            int litFloor,
            byte[] anchor,
            int anchorLen,
            int seedAPos,
            int strand) {
        int tStart = seedTStart;
        int tEnd = seedTStart + K;

        if (strand == 0) {
            int aStart = seedAPos;
            int aEnd = seedAPos + K;

            while (tStart > litFloor && aStart > 0) {
                if (target[tStart - 1] == anchor[aStart - 1]
                        || canBridgeLeftFwd(target, litFloor, tStart - 1, anchor, aStart - 1)) {
                    tStart--;
                    aStart--;
                } else {
                    break;
                }
            }
            while (tStart < seedTStart && target[tStart] != anchor[aStart]) {
                tStart++;
                aStart++;
            }

            while (tEnd < targetLen && aEnd < anchorLen) {
                if (target[tEnd] == anchor[aEnd]
                        || canBridgeRightFwd(target, targetLen, tEnd, anchor, anchorLen, aEnd)) {
                    tEnd++;
                    aEnd++;
                } else {
                    break;
                }
            }
            while (tEnd > seedTStart + K && target[tEnd - 1] != anchor[aEnd - 1]) {
                tEnd--;
                aEnd--;
            }
            return new MatchResult(0, aStart, tStart, tEnd - tStart);
        } else {
            int aStart = seedAPos;
            int aEnd = seedAPos + K;

            while (tStart > litFloor && aEnd < anchorLen) {
                if (target[tStart - 1] == (anchor[aEnd] ^ 3)
                        || canBridgeLeftRev(target, litFloor, tStart - 1, anchor, anchorLen, aEnd)) {
                    tStart--;
                    aEnd++;
                } else {
                    break;
                }
            }
            while (tStart < seedTStart && target[tStart] != (anchor[aEnd - 1] ^ 3)) {
                tStart++;
                aEnd--;
            }

            while (tEnd < targetLen && aStart > 0) {
                if (target[tEnd] == (anchor[aStart - 1] ^ 3)
                        || canBridgeRightRev(target, targetLen, tEnd, anchor, aStart - 1)) {
                    tEnd++;
                    aStart--;
                } else {
                    break;
                }
            }
            while (tEnd > seedTStart + K && target[tEnd - 1] != (anchor[aStart] ^ 3)) {
                tEnd--;
                aStart++;
            }
            return new MatchResult(1, aStart, tStart, tEnd - tStart);
        }
    }

    private static boolean canBridgeRightFwd(
            byte[] target, int targetLen, int tMismatch, byte[] anchor, int anchorLen, int aMismatch) {
        int avail = Math.min(16, Math.min(targetLen - tMismatch - 1, anchorLen - aMismatch - 1));
        if (avail < 3) return false;
        int matches = 0;
        for (int d = 1; d <= avail; d++) {
            if (target[tMismatch + d] == anchor[aMismatch + d]) matches++;
        }
        return (avail == 16) ? (matches >= 11) : (matches >= avail - (avail >= 8 ? 1 : 0));
    }

    private static boolean canBridgeLeftFwd(byte[] target, int litFloor, int tMismatch, byte[] anchor, int aMismatch) {
        int avail = Math.min(16, Math.min(tMismatch - litFloor, aMismatch));
        if (avail < 3) return false;
        int matches = 0;
        for (int d = 1; d <= avail; d++) {
            if (target[tMismatch - d] == anchor[aMismatch - d]) matches++;
        }
        return (avail == 16) ? (matches >= 11) : (matches >= avail - (avail >= 8 ? 1 : 0));
    }

    private static boolean canBridgeRightRev(byte[] target, int targetLen, int tMismatch, byte[] anchor, int aMismatch) {
        int avail = Math.min(16, Math.min(targetLen - tMismatch - 1, aMismatch));
        if (avail < 3) return false;
        int matches = 0;
        for (int d = 1; d <= avail; d++) {
            if (target[tMismatch + d] == (anchor[aMismatch - d] ^ 3)) matches++;
        }
        return (avail == 16) ? (matches >= 11) : (matches >= avail - (avail >= 8 ? 1 : 0));
    }

    private static boolean canBridgeLeftRev(
            byte[] target, int litFloor, int tMismatch, byte[] anchor, int anchorLen, int aMismatch) {
        int avail = Math.min(16, Math.min(tMismatch - litFloor, anchorLen - aMismatch - 1));
        if (avail < 3) return false;
        int matches = 0;
        for (int d = 1; d <= avail; d++) {
            if (target[tMismatch - d] == (anchor[aMismatch + d] ^ 3)) matches++;
        }
        return (avail == 16) ? (matches >= 11) : (matches >= avail - (avail >= 8 ? 1 : 0));
    }

    private static final class ParsedAlignment {
        int[] litLens;
        int tailLit;
        byte[] literals;
        int totalLiterals;
        byte[] strands;
        int[] refOffsets;
        int[] matchLens;
        int[] snpCounts;
        int numMatches;
        int[] snpGaps;
        byte[] snpDiffs;
        int totalSnps;
        int matchedBasesTotal;
    }

    private static ParsedAlignment alignTargetToAnchor(
            byte[] target,
            int targetLen,
            byte[] anchor,
            int anchorLen,
            int primarySeedLen,
            SeedHashTable index,
            boolean keepLiterals) {
        int cap = 32;
        int[] litLens = new int[cap];
        byte[] strands = new byte[cap];
        int[] refOffsets = new int[cap];
        int[] matchLens = new int[cap];
        int[] snpCounts = new int[cap];
        int numMatches = 0;

        int[] snpGaps = new int[128];
        byte[] snpDiffs = new byte[128];
        int totalSnps = 0;

        byte[] literals = keepLiterals ? new byte[Math.min(targetLen, 4096)] : EMPTY_BYTES;
        int totalLiterals = 0;
        int matchedBasesTotal = 0;

        int litStart = 0;
        int pos = 0;
        long rollingFwd = 0L;
        int validBases = 0;

        while (pos < targetLen) {
            rollingFwd = ((rollingFwd << 2) | (target[pos] & 3)) & KMER_MASK;
            validBases++;
            if (validBases >= K) {
                int queryStart = pos - K + 1;
                MatchResult m = findBestMatch(
                        target, targetLen, queryStart, litStart, rollingFwd, anchor, anchorLen, primarySeedLen, index);
                if (m != null && m.matchLen >= MIN_MATCH_LEN) {
                    int litLen = m.targetStart - litStart;
                    if (litLen > 0) {
                        if (keepLiterals) {
                            literals = ensureCapacity(literals, totalLiterals + litLen);
                            System.arraycopy(target, litStart, literals, totalLiterals, litLen);
                        }
                        totalLiterals += litLen;
                    }

                    if (numMatches == litLens.length) {
                        int nc = litLens.length << 1;
                        litLens = Arrays.copyOf(litLens, nc);
                        strands = Arrays.copyOf(strands, nc);
                        refOffsets = Arrays.copyOf(refOffsets, nc);
                        matchLens = Arrays.copyOf(matchLens, nc);
                        snpCounts = Arrays.copyOf(snpCounts, nc);
                    }

                    int matchSnps = 0;
                    int gapSinceLastSnp = 0;
                    int aStart = m.anchorStart;
                    int mLen = m.matchLen;
                    int tBase = m.targetStart;
                    if (m.strand == 0) {
                        for (int k = 0; k < mLen; k++) {
                            int xor = (target[tBase + k] ^ anchor[aStart + k]) & 3;
                            if (xor != 0) {
                                if (totalSnps == snpGaps.length) {
                                    snpGaps = Arrays.copyOf(snpGaps, snpGaps.length << 1);
                                    snpDiffs = Arrays.copyOf(snpDiffs, snpDiffs.length << 1);
                                }
                                snpGaps[totalSnps] = gapSinceLastSnp;
                                snpDiffs[totalSnps] = (byte) xor;
                                totalSnps++;
                                matchSnps++;
                                gapSinceLastSnp = 0;
                            } else {
                                gapSinceLastSnp++;
                            }
                        }
                    } else {
                        int tTop = tBase + mLen - 1;
                        for (int k = 0; k < mLen; k++) {
                            int refB = (anchor[aStart + k] ^ 3) & 3;
                            int xor = (target[tTop - k] ^ refB) & 3;
                            if (xor != 0) {
                                if (totalSnps == snpGaps.length) {
                                    snpGaps = Arrays.copyOf(snpGaps, snpGaps.length << 1);
                                    snpDiffs = Arrays.copyOf(snpDiffs, snpDiffs.length << 1);
                                }
                                snpGaps[totalSnps] = gapSinceLastSnp;
                                snpDiffs[totalSnps] = (byte) xor;
                                totalSnps++;
                                matchSnps++;
                                gapSinceLastSnp = 0;
                            } else {
                                gapSinceLastSnp++;
                            }
                        }
                    }

                    litLens[numMatches] = litLen;
                    strands[numMatches] = (byte) m.strand;
                    refOffsets[numMatches] = aStart;
                    matchLens[numMatches] = mLen;
                    snpCounts[numMatches] = matchSnps;
                    numMatches++;
                    matchedBasesTotal += mLen;

                    pos = m.targetStart + mLen;
                    litStart = pos;
                    rollingFwd = 0L;
                    validBases = 0;
                    continue;
                }
            }
            pos++;
        }

        int tailLit = targetLen - litStart;
        if (tailLit > 0) {
            if (keepLiterals) {
                literals = ensureCapacity(literals, totalLiterals + tailLit);
                System.arraycopy(target, litStart, literals, totalLiterals, tailLit);
            }
            totalLiterals += tailLit;
        }

        ParsedAlignment res = new ParsedAlignment();
        res.litLens = litLens;
        res.tailLit = tailLit;
        res.literals = literals;
        res.totalLiterals = totalLiterals;
        res.strands = strands;
        res.refOffsets = refOffsets;
        res.matchLens = matchLens;
        res.snpCounts = snpCounts;
        res.numMatches = numMatches;
        res.snpGaps = snpGaps;
        res.snpDiffs = snpDiffs;
        res.totalSnps = totalSnps;
        res.matchedBasesTotal = matchedBasesTotal;
        return res;
    }

    private static final byte[] EMPTY_BYTES = new byte[0];

    /**
     * Encodes target[0 .. targetLen) selecting the strictly smallest lossless representation among:
     * Mode 1 (QuantumContextRepeatMixer), Mode 2 (BC-LZG), and Mode 3 (AC-HMC + VCAI).
     */
    public static byte[] encodeBlockDna(
            byte[] target, int targetLen, AnchorBank bank, boolean leanMode) {
        return encodeBlockDna(target, targetLen, bank, leanMode, false);
    }

    static byte[] encodeBlockDna(
            byte[] target, int targetLen, AnchorBank bank, boolean leanMode, boolean costAware) {
        if (targetLen == 0) {
            return new byte[] {0}; // mode 0: empty
        }
        if (bank == null || bank.length < MIN_MATCH_LEN) {
            return encodeMode1Qcrm(target, targetLen, leanMode);
        }

        ParsedAlignment aln = alignTargetToAnchor(
                target, targetLen, bank.bases, bank.length, bank.primarySeedLen, bank.index, true);

        if (aln.numMatches == 0 || aln.matchedBasesTotal * 5 < targetLen) {
            return encodeMode1Qcrm(target, targetLen, leanMode);
        }

        // Encode shared prefix (literals + strands) once or compare Mode 2 vs Mode 3
        byte[] encodedLits = encodeResiduals(aln.literals, aln.totalLiterals, leanMode);
        byte[] mode2 = buildMode2Bytes(aln, encodedLits);
        byte[] best = mode2;

        if (aln.numMatches >= 1) {
            byte[] mode3 = buildMode3AcHmcBytes(aln, encodedLits, bank.primarySeedLen, bank.globalDonors, null);
            if (mode3.length < best.length) {
                best = mode3;
            }
            // Raw channel costs are only a proxy for Deflate. Bound the trial and retain it only
            // when its complete compressed representation beats the legacy Mode 2/3 winner.
            if (costAware && !leanMode && aln.totalSnps > 0 && aln.totalSnps <= 8_192 && aln.numMatches <= 256) {
                byte[] candidate = buildMode3AcHmcBytes(aln, encodedLits, bank.primarySeedLen,
                        bank.globalDonors, new CostAwareDonorSelector(1_048_576));
                if (candidate != null && candidate.length < best.length) best = candidate;
            }
        }

        // If < 55% of block aligned (e.g. single partially repetitive genome), also check pure QCRM!
        if (aln.matchedBasesTotal * 20 < targetLen * 11) {
            byte[] mode1 = encodeMode1Qcrm(target, targetLen, leanMode);
            if (mode1.length <= best.length) {
                return mode1;
            }
        }
        return best;
    }

    private static byte[] encodeMode1Qcrm(byte[] target, int targetLen, boolean leanMode) {
        ByteOut out = new ByteOut((targetLen >>> 2) + 16);
        out.writeByte(1); // Mode 1: pure QuantumContextRepeatMixer
        byte[] rc = QuantumContextRepeatMixer.encode(target, targetLen, leanMode);
        out.writeVarInt(rc.length);
        out.writeBytes(rc);
        return out.toByteArray();
    }

    private static byte[] encodeResiduals(byte[] literals, int totalLiterals, boolean leanMode) {
        ByteOut out = new ByteOut(Math.max(16, (totalLiterals >>> 2) + 16));
        if (totalLiterals >= 2048) {
            out.writeByte(1);
            byte[] litRc = QuantumContextRepeatMixer.encode(literals, totalLiterals, leanMode);
            out.writeVarInt(litRc.length);
            out.writeBytes(litRc);
        } else if (totalLiterals > 0) {
            out.writeByte(0);
            byte[] packed = VarIntCodec.pack2Bit(literals, 0, totalLiterals);
            VarIntCodec.writeCompressedAux(out, packed, packed.length);
        }
        return out.toByteArray();
    }

    private static byte[] buildMode2Bytes(ParsedAlignment aln, byte[] encodedLits) {
        int numMatches = aln.numMatches;
        ByteOut out = new ByteOut(256 + encodedLits.length);
        out.writeByte(2);
        out.writeVarInt(numMatches);
        out.writeVarInt(aln.totalLiterals);
        out.writeVarInt(aln.totalSnps);

        ByteOut litLenBuf = new ByteOut(numMatches * 2 + 4);
        for (int i = 0; i < numMatches; i++) {
            litLenBuf.writeVarInt(aln.litLens[i]);
        }
        litLenBuf.writeVarInt(aln.tailLit);
        VarIntCodec.writeCompressedAux(out, litLenBuf);

        out.writeBytes(encodedLits);

        byte[] strandBits = new byte[(numMatches + 7) >>> 3];
        for (int i = 0; i < numMatches; i++) {
            if (aln.strands[i] != 0) {
                strandBits[i >>> 3] |= (byte) (1 << (i & 7));
            }
        }
        VarIntCodec.writeCompressedAux(out, strandBits, strandBits.length);

        ByteOut geomBuf = new ByteOut(numMatches * 4);
        int prevRef = 0;
        for (int i = 0; i < numMatches; i++) {
            geomBuf.writeSignedVarInt(aln.refOffsets[i] - prevRef);
            prevRef = aln.refOffsets[i];
        }
        int prevLen = 0;
        for (int i = 0; i < numMatches; i++) {
            geomBuf.writeSignedVarInt(aln.matchLens[i] - prevLen);
            prevLen = aln.matchLens[i];
        }
        VarIntCodec.writeCompressedAux(out, geomBuf);

        if (aln.totalSnps > 0) {
            ByteOut snpCntBuf = new ByteOut(numMatches * 2);
            for (int i = 0; i < numMatches; i++) {
                snpCntBuf.writeVarInt(aln.snpCounts[i]);
            }
            VarIntCodec.writeCompressedAux(out, snpCntBuf);

            ByteOut snpGapBuf = new ByteOut(aln.totalSnps * 2);
            for (int i = 0; i < aln.totalSnps; i++) {
                snpGapBuf.writeVarInt(aln.snpGaps[i]);
            }
            VarIntCodec.writeCompressedAux(out, snpGapBuf);
            VarIntCodec.writeCompressedAux(out, aln.snpDiffs, aln.totalSnps);
        }
        return out.toByteArray();
    }

    // =========================================================================
    // MODE 3: ANCHOR-COORDINATE HAPLOTYPE MATRIX + VIRTUAL CLADE-ANCESTOR (AC-HMC)
    // =========================================================================

    /**
     * Deterministically partitions matches [0 .. numMatches) into disjoint-anchor Haplotype Sweeps
     * using ONLY (refOffsets, matchLens) so encoder and decoder agree with 0 framing bits!
     */
    private static int[] computeSweepStarts(int[] refOffsets, int[] matchLens, int numMatches) {
        int[] starts = new int[numMatches + 1];
        int numSweeps = 0;
        starts[numSweeps++] = 0;
        int curSweepFirst = 0;

        for (int m = 1; m < numMatches; m++) {
            int s = refOffsets[m];
            int e = s + matchLens[m];
            boolean overlaps = (m - curSweepFirst >= 16);
            if (!overlaps) {
                for (int prev = curSweepFirst; prev < m; prev++) {
                    int ps = refOffsets[prev];
                    int pe = ps + matchLens[prev];
                    if (s < pe && ps < e) {
                        overlaps = true;
                        break;
                    }
                }
            }
            if (overlaps) {
                starts[numSweeps++] = m;
                curSweepFirst = m;
            }
        }
        starts[numSweeps] = numMatches;
        return Arrays.copyOf(starts, numSweeps + 1);
    }

    /** Builds canonical coordinate-sorted variant arrays V_g = {(anchorPos << 2) | xorDiff} per sweep. */
    private static int[][] buildSweepVariantSets(ParsedAlignment aln) {
        int[] sweepStarts = computeSweepStarts(aln.refOffsets, aln.matchLens, aln.numMatches);
        int numSweeps = sweepStarts.length - 1;
        int[][] sweeps = new int[numSweeps][];
        int snpCursor = 0;

        for (int g = 0; g < numSweeps; g++) {
            int mStart = sweepStarts[g];
            int mEnd = sweepStarts[g + 1];
            int count = 0;
            for (int m = mStart; m < mEnd; m++) {
                count += aln.snpCounts[m];
            }
            if (count == 0) {
                sweeps[g] = EMPTY_INTS;
                continue;
            }
            int[] v = new int[count];
            int w = 0;
            for (int m = mStart; m < mEnd; m++) {
                int aStart = aln.refOffsets[m];
                int mSnps = aln.snpCounts[m];
                int relPos = 0;
                for (int s = 0; s < mSnps; s++) {
                    relPos += aln.snpGaps[snpCursor];
                    int anchorPos = aStart + relPos;
                    int diff = aln.snpDiffs[snpCursor] & 3;
                    v[w++] = (anchorPos << 2) | diff;
                    relPos++;
                    snpCursor++;
                }
            }
            if (mEnd - mStart > 1) {
                Arrays.sort(v); // strictly disjoint anchor intervals -> unique anchorPos order!
            }
            sweeps[g] = v;
        }
        return sweeps;
    }

    private static byte[] buildMode3AcHmcBytes(
            ParsedAlignment aln, byte[] encodedLits, int primarySeedLen, int[][] globalDonors,
            CostAwareDonorSelector costSelector) {
        int numMatches = aln.numMatches;
        ByteOut out = new ByteOut(256 + encodedLits.length);
        out.writeByte(3); // Mode 3: AC-HMC + VCAI + Circular Geometry
        out.writeVarInt(numMatches);
        out.writeVarInt(aln.totalLiterals);
        out.writeVarInt(aln.totalSnps);

        // 1. Literal lengths + tail literal
        ByteOut litLenBuf = new ByteOut(numMatches * 2 + 4);
        for (int i = 0; i < numMatches; i++) {
            litLenBuf.writeVarInt(aln.litLens[i]);
        }
        litLenBuf.writeVarInt(aln.tailLit);
        VarIntCodec.writeCompressedAux(out, litLenBuf);

        // 2. Residual literals
        out.writeBytes(encodedLits);

        // 3. Strand bits
        byte[] strandBits = new byte[(numMatches + 7) >>> 3];
        for (int i = 0; i < numMatches; i++) {
            if (aln.strands[i] != 0) {
                strandBits[i >>> 3] |= (byte) (1 << (i & 7));
            }
        }
        VarIntCodec.writeCompressedAux(out, strandBits, strandBits.length);

        // 4. Circular-aware interleaved geometry (collapses 3 of 4 circular coordinates to 0!)
        ByteOut geomBuf = new ByteOut(numMatches * 3 + 8);
        int coreEnd = Math.max(1, primarySeedLen);
        geomBuf.writeVarInt(coreEnd);
        for (int m = 0; m < numMatches; m++) {
            int predRef = predictCircularRefOffset(m, aln.refOffsets, aln.matchLens, coreEnd);
            geomBuf.writeSignedVarInt(aln.refOffsets[m] - predRef);
            int predLen = predictCircularMatchLen(m, aln.refOffsets, aln.matchLens, coreEnd);
            geomBuf.writeSignedVarInt(aln.matchLens[m] - predLen);
        }
        VarIntCodec.writeCompressedAux(out, geomBuf);

        // 5. Canonical Anchor-Coordinate Haplotype Sweeps with Dynamic VCAI
        if (aln.totalSnps > 0) {
            int[][] sweeps = buildSweepVariantSets(aln);
            int numSweeps = sweeps.length;
            byte[] sweepModes = new byte[numSweeps];
            ByteOut donorIdxBuf = new ByteOut(numSweeps * 2);
            ByteOut maskBuf = new ByteOut(numSweeps * 8);
            ByteOut addedCntBuf = new ByteOut(numSweeps * 2);
            ByteOut addedGapBuf = new ByteOut(aln.totalSnps);
            ByteOut addedDiffBuf = new ByteOut(aln.totalSnps);

            HaplotypeDonorIndex donorPool = new HaplotypeDonorIndex(MAX_DONOR_POOL);
            if (globalDonors != null) {
                for (int[] gd : globalDonors) {
                    if (donorPool.size() < MAX_DONOR_POOL) donorPool.add(gd);
                }
            }

            for (int g = 0; g < numSweeps; g++) {
                int[] v = sweeps[g];
                if (v.length == 0) {
                    sweepModes[g] = 0;
                    continue;
                }

                int[] commonCounts = donorPool.countIntersections(v);
                int bestDonor = -1;
                int bestCommon = 0;
                if (costSelector == null) {
                    int bestSavingsScore = 8; // require saving >= 1 byte over scratch mode 1
                    int poolSize = donorPool.size();
                    for (int d = 0; d < poolSize; d++) {
                        int[] cand = donorPool.get(d);
                        int common = commonCounts[d];
                        if (common == 0) continue;
                        int maskBits = (common == cand.length) ? 0 : cand.length;
                        int score = common * 18 - maskBits - 8; // ~18 bits saved per shared SNP
                        if (score > bestSavingsScore) {
                            bestSavingsScore = score;
                            bestDonor = d;
                            bestCommon = common;
                        }
                    }
                } else {
                    bestDonor = costSelector.select(v, donorPool, commonCounts);
                    if (bestDonor == CostAwareDonorSelector.BUDGET_EXHAUSTED) return null;
                    if (bestDonor >= 0) bestCommon = commonCounts[bestDonor];
                }

                if (bestDonor < 0) {
                    sweepModes[g] = 1; // Scratch sweep
                    addedCntBuf.writeVarInt(v.length);
                    writeVariantList(v, v.length, addedGapBuf, addedDiffBuf);
                    if (donorPool.size() < MAX_DONOR_POOL) {
                        donorPool.add(v);
                    }
                } else {
                    int[] dArr = donorPool.get(bestDonor);
                    if (bestCommon == dArr.length) {
                        // Mode 2: 100% of Donor SNPs inherited (ZERO bitmask!)
                        sweepModes[g] = 2;
                        donorIdxBuf.writeVarInt(bestDonor);
                        int addCnt = v.length - dArr.length;
                        addedCntBuf.writeVarInt(addCnt);
                        if (addCnt > 0) {
                            int[] added = subtractSortedSubset(v, dArr, addCnt);
                            writeVariantList(added, addCnt, addedGapBuf, addedDiffBuf);
                            if (donorPool.size() < MAX_DONOR_POOL) {
                                donorPool.add(v);
                            }
                        }
                    } else {
                        // Mode 3: Partial inheritance with retained bitmask + VCAI clade inference!
                        sweepModes[g] = 3;
                        donorIdxBuf.writeVarInt(bestDonor);
                        byte[] maskBytes = new byte[(dArr.length + 7) >>> 3];
                        int[] retained = new int[bestCommon];
                        int[] added = new int[v.length - bestCommon];
                        diffSortedLists(v, dArr, maskBytes, retained, added);

                        maskBuf.writeBytes(maskBytes);
                        addedCntBuf.writeVarInt(added.length);
                        if (added.length > 0) {
                            writeVariantList(added, added.length, addedGapBuf, addedDiffBuf);
                        }

                        if (donorPool.size() < MAX_DONOR_POOL) {
                            donorPool.add(v);
                        }
                        if (bestCommon >= 2 && bestCommon < dArr.length && bestCommon < v.length
                                && donorPool.size() < MAX_DONOR_POOL) {
                            donorPool.add(retained); // VCAI: synthesize Virtual Clade Ancestor!
                        }
                    }
                }
            }

            VarIntCodec.writeCompressedAux(out, sweepModes, numSweeps);
            VarIntCodec.writeCompressedAux(out, donorIdxBuf);
            VarIntCodec.writeCompressedAux(out, maskBuf);
            VarIntCodec.writeCompressedAux(out, addedCntBuf);
            VarIntCodec.writeCompressedAux(out, addedGapBuf);
            VarIntCodec.writeCompressedAux(out, addedDiffBuf);
        }

        return out.toByteArray();
    }

    private static int predictCircularRefOffset(int m, int[] refOffsets, int[] matchLens, int coreEnd) {
        if (m == 0) return 0;
        int prevRef = refOffsets[m - 1];
        int prevLen = matchLens[m - 1];
        if (prevRef > 0 && prevRef + prevLen >= coreEnd) {
            return 0; // Forward circular wraparound second half starts at 0
        }
        if (prevRef == 0 && prevLen < coreEnd) {
            return prevLen; // Reverse-complement circular wraparound second half starts at origin == prevLen
        }
        return 0; // Unpermuted genome starts at 0
    }

    private static int predictCircularMatchLen(int m, int[] refOffsets, int[] matchLens, int coreEnd) {
        int curRef = refOffsets[m];
        if (m > 0 && curRef == 0 && refOffsets[m - 1] > 0 && refOffsets[m - 1] + matchLens[m - 1] == coreEnd) {
            return refOffsets[m - 1]; // Forward circular wraparound second half has length == origin
        }
        if (curRef < coreEnd) {
            return coreEnd - curRef; // Extends cleanly to primary reference end
        }
        return 0;
    }

    private static void writeVariantList(int[] list, int len, ByteOut gapBuf, ByteOut diffBuf) {
        int prevPos = -1;
        for (int i = 0; i < len; i++) {
            int pos = list[i] >>> 2;
            int diff = list[i] & 3;
            gapBuf.writeVarInt(pos - prevPos - 1);
            diffBuf.writeByte(diff);
            prevPos = pos;
        }
    }

    private static int[] readVariantList(int len, ByteIn gapIn, byte[] diffs, int[] diffCursor) {
        int[] out = new int[len];
        int prevPos = -1;
        int dc = diffCursor[0];
        for (int i = 0; i < len; i++) {
            int pos = prevPos + 1 + gapIn.readVarInt();
            int diff = diffs[dc++] & 3;
            out[i] = (pos << 2) | diff;
            prevPos = pos;
        }
        diffCursor[0] = dc;
        return out;
    }

    private static int countSortedIntersection(int[] a, int[] b) {
        int i = 0, j = 0, c = 0;
        int na = a.length, nb = b.length;
        while (i < na && j < nb) {
            int va = a[i], vb = b[j];
            if (va == vb) {
                c++;
                i++;
                j++;
            } else if (va < vb) {
                i++;
            } else {
                j++;
            }
        }
        return c;
    }

    private static int[] intersectSorted(int[] a, int[] b, int commonLen) {
        int[] out = new int[commonLen];
        int i = 0, j = 0, w = 0;
        while (i < a.length && j < b.length) {
            int va = a[i], vb = b[j];
            if (va == vb) {
                out[w++] = va;
                i++;
                j++;
            } else if (va < vb) {
                i++;
            } else {
                j++;
            }
        }
        return out;
    }

    private static int[] subtractSortedSubset(int[] superset, int[] subset, int diffLen) {
        int[] out = new int[diffLen];
        int i = 0, j = 0, w = 0;
        while (i < superset.length) {
            if (j < subset.length && superset[i] == subset[j]) {
                i++;
                j++;
            } else {
                out[w++] = superset[i++];
            }
        }
        return out;
    }

    private static void diffSortedLists(
            int[] v, int[] donor, byte[] outMask, int[] outRetained, int[] outAdded) {
        int i = 0, j = 0, rw = 0, aw = 0;
        int nv = v.length, nd = donor.length;
        while (i < nv && j < nd) {
            int vv = v[i], dv = donor[j];
            if (vv == dv) {
                outMask[j >>> 3] |= (byte) (1 << (j & 7));
                outRetained[rw++] = dv;
                i++;
                j++;
            } else if (vv < dv) {
                outAdded[aw++] = vv;
                i++;
            } else {
                j++;
            }
        }
        while (i < nv) {
            outAdded[aw++] = v[i++];
        }
    }

    private static int[] mergeSortedDisjoint(int[] a, int[] b) {
        if (a.length == 0) return b;
        if (b.length == 0) return a;
        int[] out = new int[a.length + b.length];
        int i = 0, j = 0, w = 0;
        while (i < a.length && j < b.length) {
            if (a[i] < b[j]) {
                out[w++] = a[i++];
            } else {
                out[w++] = b[j++];
            }
        }
        while (i < a.length) out[w++] = a[i++];
        while (j < b.length) out[w++] = b[j++];
        return out;
    }

    private static int[][] decodeMode3SweepVariants(
            ByteIn in, int[] sweepStarts, int[][] globalDonors) {
        int numSweeps = sweepStarts.length - 1;
        byte[] sweepModes = VarIntCodec.readCompressedAux(in);
        ByteIn donorIdxIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        ByteIn maskIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        ByteIn addedCntIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        ByteIn addedGapIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        byte[] addedDiffs = VarIntCodec.readCompressedAux(in);
        int[] diffCursor = {0};

        List<int[]> donorPool = new ArrayList<>(MAX_DONOR_POOL);
        if (globalDonors != null) {
            for (int[] gd : globalDonors) {
                if (donorPool.size() < MAX_DONOR_POOL) donorPool.add(gd);
            }
        }

        int[][] sweeps = new int[numSweeps][];
        for (int g = 0; g < numSweeps; g++) {
            int sm = sweepModes[g] & 0xFF;
            if (sm == 0) {
                sweeps[g] = EMPTY_INTS;
            } else if (sm == 1) {
                int cnt = addedCntIn.readVarInt();
                int[] v = readVariantList(cnt, addedGapIn, addedDiffs, diffCursor);
                sweeps[g] = v;
                if (donorPool.size() < MAX_DONOR_POOL) {
                    donorPool.add(v);
                }
            } else if (sm == 2) {
                int dIdx = donorIdxIn.readVarInt();
                int[] dArr = donorPool.get(dIdx);
                int addCnt = addedCntIn.readVarInt();
                if (addCnt == 0) {
                    sweeps[g] = dArr;
                } else {
                    int[] added = readVariantList(addCnt, addedGapIn, addedDiffs, diffCursor);
                    int[] v = mergeSortedDisjoint(dArr, added);
                    sweeps[g] = v;
                    if (donorPool.size() < MAX_DONOR_POOL) {
                        donorPool.add(v);
                    }
                }
            } else if (sm == 3) {
                int dIdx = donorIdxIn.readVarInt();
                int[] dArr = donorPool.get(dIdx);
                int maskByteLen = (dArr.length + 7) >>> 3;
                byte[] maskBytes = maskIn.readBytes(maskByteLen);
                int retCnt = 0;
                for (int j = 0; j < dArr.length; j++) {
                    if (((maskBytes[j >>> 3] >>> (j & 7)) & 1) != 0) retCnt++;
                }
                int[] retained = new int[retCnt];
                int rw = 0;
                for (int j = 0; j < dArr.length; j++) {
                    if (((maskBytes[j >>> 3] >>> (j & 7)) & 1) != 0) {
                        retained[rw++] = dArr[j];
                    }
                }
                int addCnt = addedCntIn.readVarInt();
                int[] added = (addCnt > 0)
                        ? readVariantList(addCnt, addedGapIn, addedDiffs, diffCursor)
                        : EMPTY_INTS;
                int[] v = mergeSortedDisjoint(retained, added);
                sweeps[g] = v;

                if (donorPool.size() < MAX_DONOR_POOL) {
                    donorPool.add(v);
                }
                if (retCnt >= 2 && retCnt < dArr.length && retCnt < v.length
                        && donorPool.size() < MAX_DONOR_POOL) {
                    donorPool.add(retained);
                }
            } else {
                throw new IllegalStateException("Unknown AC-HMC sweep mode: " + sm);
            }
        }
        return sweeps;
    }

    /**
     * Inspects which AnchorBank page indices (each ANCHOR_PAGE_SIZE bases) are referenced by a
     * compressed Mode 2 or Mode 3 payload without materializing the output!
     */
    public static boolean[] findReferencedAnchorPages(ByteIn in, int numAnchorPages, boolean leanMode) {
        int savedPos = in.position();
        boolean[] needed = new boolean[numAnchorPages];
        int mode = in.readByte();
        if ((mode != 2 && mode != 3) || numAnchorPages == 0) {
            in.setPosition(savedPos);
            return needed;
        }
        int numMatches = in.readVarInt();
        int totalLiterals = in.readVarInt();
        in.readVarInt(); // totalSnps

        VarIntCodec.readCompressedAux(in); // litLenBuf
        skipResiduals(in, totalLiterals);
        VarIntCodec.readCompressedAux(in); // strandBits

        ByteIn geomIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        int[] refOffsets = new int[numMatches];
        int[] matchLens = new int[numMatches];
        if (mode == 2) {
            int prevRef = 0;
            for (int i = 0; i < numMatches; i++) {
                prevRef += geomIn.readSignedVarInt();
                refOffsets[i] = prevRef;
            }
            int prevLen = 0;
            for (int i = 0; i < numMatches; i++) {
                prevLen += geomIn.readSignedVarInt();
                matchLens[i] = prevLen;
            }
        } else {
            int coreEnd = geomIn.readVarInt();
            for (int m = 0; m < numMatches; m++) {
                int predRef = predictCircularRefOffset(m, refOffsets, matchLens, coreEnd);
                refOffsets[m] = predRef + geomIn.readSignedVarInt();
                int predLen = predictCircularMatchLen(m, refOffsets, matchLens, coreEnd);
                matchLens[m] = predLen + geomIn.readSignedVarInt();
            }
        }

        for (int i = 0; i < numMatches; i++) {
            int startPage = refOffsets[i] / ANCHOR_PAGE_SIZE;
            int endPage = Math.max(startPage, (refOffsets[i] + matchLens[i] - 1) / ANCHOR_PAGE_SIZE);
            for (int p = startPage; p <= endPage && p < numAnchorPages; p++) {
                if (p >= 0) needed[p] = true;
            }
        }
        in.setPosition(savedPos);
        return needed;
    }

    private static void skipResiduals(ByteIn in, int totalLiterals) {
        if (totalLiterals >= 2048) {
            int litMode = in.readByte();
            if (litMode == 1) {
                int rcLen = in.readVarInt();
                in.setPosition(in.position() + rcLen);
            } else {
                VarIntCodec.readCompressedAux(in);
            }
        } else if (totalLiterals > 0) {
            in.readByte();
            VarIntCodec.readCompressedAux(in);
        }
    }

    private static byte[] decodeResiduals(ByteIn in, int totalLiterals, boolean leanMode) {
        if (totalLiterals == 0) {
            return EMPTY_BYTES;
        }
        int litMode = in.readByte();
        if (litMode == 1) {
            int rcLen = in.readVarInt();
            ByteIn rcSlice = new ByteIn(in.readBytes(rcLen));
            return QuantumContextRepeatMixer.decode(rcSlice, totalLiterals, leanMode);
        } else {
            byte[] packed = VarIntCodec.readCompressedAux(in);
            return VarIntCodec.unpack2Bit(packed, totalLiterals);
        }
    }

    public static byte[] decodeBlockDna(ByteIn in, int targetLen, byte[] anchorBases, boolean leanMode) {
        return decodeBlockDna(in, targetLen, anchorBases, EMPTY_DONORS, leanMode);
    }

    /** Decodes a block's primary 2-bit DNA sequence from Mode 1 (QCRM), Mode 2 (BC-LZG), or Mode 3 (AC-HMC). */
    public static byte[] decodeBlockDna(
            ByteIn in, int targetLen, byte[] anchorBases, int[][] globalDonors, boolean leanMode) {
        int mode = in.readByte();
        if (mode == 0 || targetLen == 0) {
            return EMPTY_BYTES;
        }
        if (mode == 1) {
            int rcLen = in.readVarInt();
            ByteIn rcSlice = new ByteIn(in.readBytes(rcLen));
            return QuantumContextRepeatMixer.decode(rcSlice, targetLen, leanMode);
        }
        if (mode != 2 && mode != 3) {
            throw new IllegalStateException("Unknown DNA block mode: " + mode);
        }

        int numMatches = in.readVarInt();
        int totalLiterals = in.readVarInt();
        int totalSnps = in.readVarInt();

        ByteIn litLenIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        int[] litLens = new int[numMatches];
        for (int i = 0; i < numMatches; i++) {
            litLens[i] = litLenIn.readVarInt();
        }
        int tailLit = litLenIn.readVarInt();

        byte[] literals = decodeResiduals(in, totalLiterals, leanMode);
        byte[] strandBits = VarIntCodec.readCompressedAux(in);

        ByteIn geomIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        int[] refOffsets = new int[numMatches];
        int[] matchLens = new int[numMatches];
        if (mode == 2) {
            int prevRef = 0;
            for (int i = 0; i < numMatches; i++) {
                prevRef += geomIn.readSignedVarInt();
                refOffsets[i] = prevRef;
            }
            int prevLen = 0;
            for (int i = 0; i < numMatches; i++) {
                prevLen += geomIn.readSignedVarInt();
                matchLens[i] = prevLen;
            }
        } else {
            int coreEnd = geomIn.readVarInt();
            for (int m = 0; m < numMatches; m++) {
                int predRef = predictCircularRefOffset(m, refOffsets, matchLens, coreEnd);
                refOffsets[m] = predRef + geomIn.readSignedVarInt();
                int predLen = predictCircularMatchLen(m, refOffsets, matchLens, coreEnd);
                matchLens[m] = predLen + geomIn.readSignedVarInt();
            }
        }

        if (mode == 2) {
            int[] snpCounts = new int[numMatches];
            int[] snpGaps = null;
            byte[] snpDiffs = null;
            if (totalSnps > 0) {
                ByteIn snpCntIn = new ByteIn(VarIntCodec.readCompressedAux(in));
                for (int i = 0; i < numMatches; i++) {
                    snpCounts[i] = snpCntIn.readVarInt();
                }
                ByteIn snpGapIn = new ByteIn(VarIntCodec.readCompressedAux(in));
                snpGaps = new int[totalSnps];
                for (int i = 0; i < totalSnps; i++) {
                    snpGaps[i] = snpGapIn.readVarInt();
                }
                snpDiffs = VarIntCodec.readCompressedAux(in);
            }

            byte[] out = new byte[targetLen];
            int outPos = 0;
            int litPos = 0;
            int snpIdx = 0;
            for (int m = 0; m < numMatches; m++) {
                int lLen = litLens[m];
                if (lLen > 0) {
                    System.arraycopy(literals, litPos, out, outPos, lLen);
                    litPos += lLen;
                    outPos += lLen;
                }
                int strand = (strandBits[m >>> 3] >>> (m & 7)) & 1;
                int aStart = refOffsets[m];
                int mLen = matchLens[m];
                if (strand == 0) {
                    System.arraycopy(anchorBases, aStart, out, outPos, mLen);
                } else {
                    int aTop = aStart + mLen - 1;
                    for (int d = 0; d < mLen; d++) {
                        out[outPos + d] = (byte) (anchorBases[aTop - d] ^ 3);
                    }
                }
                int mSnps = snpCounts[m];
                if (mSnps > 0) {
                    int relPos = 0;
                    if (strand == 0) {
                        for (int s = 0; s < mSnps; s++) {
                            relPos += snpGaps[snpIdx];
                            out[outPos + relPos] ^= snpDiffs[snpIdx];
                            relPos++;
                            snpIdx++;
                        }
                    } else {
                        int outTop = outPos + mLen - 1;
                        for (int s = 0; s < mSnps; s++) {
                            relPos += snpGaps[snpIdx];
                            out[outTop - relPos] ^= snpDiffs[snpIdx];
                            relPos++;
                            snpIdx++;
                        }
                    }
                }
                outPos += mLen;
            }
            if (tailLit > 0) {
                System.arraycopy(literals, litPos, out, outPos, tailLit);
            }
            return out;
        }

        // Mode 3 reconstruction: copy template matches + apply canonical Sweep Haplotype Variants
        int[] sweepStarts = null;
        int[][] sweepVariants = null;
        if (totalSnps > 0) {
            sweepStarts = computeSweepStarts(refOffsets, matchLens, numMatches);
            sweepVariants = decodeMode3SweepVariants(in, sweepStarts, globalDonors);
        }

        byte[] out = new byte[targetLen];
        int[] matchTargetStarts = new int[numMatches];
        int outPos = 0;
        int litPos = 0;
        for (int m = 0; m < numMatches; m++) {
            int lLen = litLens[m];
            if (lLen > 0) {
                System.arraycopy(literals, litPos, out, outPos, lLen);
                litPos += lLen;
                outPos += lLen;
            }
            matchTargetStarts[m] = outPos;
            int strand = (strandBits[m >>> 3] >>> (m & 7)) & 1;
            int aStart = refOffsets[m];
            int mLen = matchLens[m];
            if (strand == 0) {
                System.arraycopy(anchorBases, aStart, out, outPos, mLen);
            } else {
                int aTop = aStart + mLen - 1;
                for (int d = 0; d < mLen; d++) {
                    out[outPos + d] = (byte) (anchorBases[aTop - d] ^ 3);
                }
            }
            outPos += mLen;
        }
        if (tailLit > 0) {
            System.arraycopy(literals, litPos, out, outPos, tailLit);
        }

        if (totalSnps > 0) {
            int numSweeps = sweepStarts.length - 1;
            for (int g = 0; g < numSweeps; g++) {
                int[] vList = sweepVariants[g];
                if (vList.length == 0) continue;
                int mStart = sweepStarts[g];
                int mEnd = sweepStarts[g + 1];
                for (int v : vList) {
                    int anchorPos = v >>> 2;
                    int xorDiff = v & 3;
                    for (int m = mStart; m < mEnd; m++) {
                        int aS = refOffsets[m];
                        int mL = matchLens[m];
                        if (anchorPos >= aS && anchorPos < aS + mL) {
                            int k = anchorPos - aS;
                            int strand = (strandBits[m >>> 3] >>> (m & 7)) & 1;
                            int tIdx = (strand == 0)
                                    ? (matchTargetStarts[m] + k)
                                    : (matchTargetStarts[m] + mL - 1 - k);
                            out[tIdx] ^= (byte) xorDiff;
                            break;
                        }
                    }
                }
            }
        }
        return out;
    }

    public static byte[] decodeBlockDnaRange(
            ByteIn in, int targetLen, int rangeStart, int rangeEnd, byte[] anchorBases, boolean leanMode) {
        return decodeBlockDnaRange(in, targetLen, rangeStart, rangeEnd, anchorBases, EMPTY_DONORS, leanMode);
    }

    /**
     * Sub-interval selective decoder for BOTH Mode 2 (BC-LZG) and Mode 3 (AC-HMC): materializes
     * ONLY target[rangeStart .. rangeEnd) without allocating or copying the rest of the block!
     */
    public static byte[] decodeBlockDnaRange(
            ByteIn in,
            int targetLen,
            int rangeStart,
            int rangeEnd,
            byte[] anchorBases,
            int[][] globalDonors,
            boolean leanMode) {
        if (rangeStart <= 0 && rangeEnd >= targetLen) {
            return decodeBlockDna(in, targetLen, anchorBases, globalDonors, leanMode);
        }
        if (rangeStart >= rangeEnd || targetLen == 0) {
            return EMPTY_BYTES;
        }
        int savedPos = in.position();
        int mode = in.readByte();
        if (mode != 2 && mode != 3) {
            in.setPosition(savedPos);
            byte[] full = decodeBlockDna(in, targetLen, anchorBases, globalDonors, leanMode);
            return Arrays.copyOfRange(full, rangeStart, rangeEnd);
        }

        int numMatches = in.readVarInt();
        int totalLiterals = in.readVarInt();
        int totalSnps = in.readVarInt();

        ByteIn litLenIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        int[] litLens = new int[numMatches];
        for (int i = 0; i < numMatches; i++) {
            litLens[i] = litLenIn.readVarInt();
        }
        int tailLit = litLenIn.readVarInt();

        byte[] literals = decodeResiduals(in, totalLiterals, leanMode);
        byte[] strandBits = VarIntCodec.readCompressedAux(in);

        ByteIn geomIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        int[] refOffsets = new int[numMatches];
        int[] matchLens = new int[numMatches];
        if (mode == 2) {
            int prevRef = 0;
            for (int i = 0; i < numMatches; i++) {
                prevRef += geomIn.readSignedVarInt();
                refOffsets[i] = prevRef;
            }
            int prevLen = 0;
            for (int i = 0; i < numMatches; i++) {
                prevLen += geomIn.readSignedVarInt();
                matchLens[i] = prevLen;
            }
        } else {
            int coreEnd = geomIn.readVarInt();
            for (int m = 0; m < numMatches; m++) {
                int predRef = predictCircularRefOffset(m, refOffsets, matchLens, coreEnd);
                refOffsets[m] = predRef + geomIn.readSignedVarInt();
                int predLen = predictCircularMatchLen(m, refOffsets, matchLens, coreEnd);
                matchLens[m] = predLen + geomIn.readSignedVarInt();
            }
        }

        int sliceLen = rangeEnd - rangeStart;
        byte[] out = new byte[sliceLen];

        if (mode == 2) {
            int[] snpCounts = new int[numMatches];
            int[] snpGaps = null;
            byte[] snpDiffs = null;
            if (totalSnps > 0) {
                ByteIn snpCntIn = new ByteIn(VarIntCodec.readCompressedAux(in));
                for (int i = 0; i < numMatches; i++) {
                    snpCounts[i] = snpCntIn.readVarInt();
                }
                ByteIn snpGapIn = new ByteIn(VarIntCodec.readCompressedAux(in));
                snpGaps = new int[totalSnps];
                for (int i = 0; i < totalSnps; i++) {
                    snpGaps[i] = snpGapIn.readVarInt();
                }
                snpDiffs = VarIntCodec.readCompressedAux(in);
            }

            int outPos = 0;
            int litPos = 0;
            int snpIdx = 0;
            for (int m = 0; m < numMatches; m++) {
                int lLen = litLens[m];
                if (lLen > 0) {
                    int segS = Math.max(outPos, rangeStart);
                    int segE = Math.min(outPos + lLen, rangeEnd);
                    if (segS < segE) {
                        System.arraycopy(literals, litPos + (segS - outPos), out, segS - rangeStart, segE - segS);
                    }
                    litPos += lLen;
                    outPos += lLen;
                }
                int mLen = matchLens[m];
                int mSnps = snpCounts[m];
                int segS = Math.max(outPos, rangeStart);
                int segE = Math.min(outPos + mLen, rangeEnd);

                if (segS < segE) {
                    int strand = (strandBits[m >>> 3] >>> (m & 7)) & 1;
                    int aStart = refOffsets[m];
                    int relS = segS - outPos;
                    int copyLen = segE - segS;
                    int dstOff = segS - rangeStart;
                    if (strand == 0) {
                        System.arraycopy(anchorBases, aStart + relS, out, dstOff, copyLen);
                    } else {
                        int aTop = aStart + mLen - 1;
                        for (int d = 0; d < copyLen; d++) {
                            out[dstOff + d] = (byte) (anchorBases[aTop - (relS + d)] ^ 3);
                        }
                    }
                    if (mSnps > 0) {
                        int relPos = 0;
                        if (strand == 0) {
                            for (int s = 0; s < mSnps; s++) {
                                relPos += snpGaps[snpIdx];
                                int absT = outPos + relPos;
                                if (absT >= segS && absT < segE) {
                                    out[absT - rangeStart] ^= snpDiffs[snpIdx];
                                }
                                relPos++;
                                snpIdx++;
                            }
                        } else {
                            int outTop = outPos + mLen - 1;
                            for (int s = 0; s < mSnps; s++) {
                                relPos += snpGaps[snpIdx];
                                int absT = outTop - relPos;
                                if (absT >= segS && absT < segE) {
                                    out[absT - rangeStart] ^= snpDiffs[snpIdx];
                                }
                                relPos++;
                                snpIdx++;
                            }
                        }
                    }
                } else {
                    snpIdx += mSnps;
                }
                outPos += mLen;
                if (outPos >= rangeEnd) {
                    return out;
                }
            }
            if (tailLit > 0 && outPos < rangeEnd) {
                int segS = Math.max(outPos, rangeStart);
                int segE = Math.min(outPos + tailLit, rangeEnd);
                if (segS < segE) {
                    System.arraycopy(literals, litPos + (segS - outPos), out, segS - rangeStart, segE - segS);
                }
            }
            return out;
        }

        // Mode 3 selective sub-interval reconstruction
        int[] sweepStarts = null;
        int[][] sweepVariants = null;
        if (totalSnps > 0) {
            sweepStarts = computeSweepStarts(refOffsets, matchLens, numMatches);
            sweepVariants = decodeMode3SweepVariants(in, sweepStarts, globalDonors);
        }

        int[] matchTargetStarts = new int[numMatches];
        int outPos = 0;
        int litPos = 0;
        for (int m = 0; m < numMatches; m++) {
            int lLen = litLens[m];
            if (lLen > 0) {
                int segS = Math.max(outPos, rangeStart);
                int segE = Math.min(outPos + lLen, rangeEnd);
                if (segS < segE) {
                    System.arraycopy(literals, litPos + (segS - outPos), out, segS - rangeStart, segE - segS);
                }
                litPos += lLen;
                outPos += lLen;
            }
            matchTargetStarts[m] = outPos;
            int mLen = matchLens[m];
            int segS = Math.max(outPos, rangeStart);
            int segE = Math.min(outPos + mLen, rangeEnd);
            if (segS < segE) {
                int strand = (strandBits[m >>> 3] >>> (m & 7)) & 1;
                int aStart = refOffsets[m];
                int relS = segS - outPos;
                int copyLen = segE - segS;
                int dstOff = segS - rangeStart;
                if (strand == 0) {
                    System.arraycopy(anchorBases, aStart + relS, out, dstOff, copyLen);
                } else {
                    int aTop = aStart + mLen - 1;
                    for (int d = 0; d < copyLen; d++) {
                        out[dstOff + d] = (byte) (anchorBases[aTop - (relS + d)] ^ 3);
                    }
                }
            }
            outPos += mLen;
        }
        if (tailLit > 0 && outPos < rangeEnd) {
            int segS = Math.max(outPos, rangeStart);
            int segE = Math.min(outPos + tailLit, rangeEnd);
            if (segS < segE) {
                System.arraycopy(literals, litPos + (segS - outPos), out, segS - rangeStart, segE - segS);
            }
        }

        if (totalSnps > 0) {
            int numSweeps = sweepStarts.length - 1;
            for (int g = 0; g < numSweeps; g++) {
                int mStart = sweepStarts[g];
                int mEnd = sweepStarts[g + 1];
                int sweepTargetStart = matchTargetStarts[mStart];
                int lastM = mEnd - 1;
                int sweepTargetEnd = matchTargetStarts[lastM] + matchLens[lastM];
                if (sweepTargetEnd <= rangeStart) continue;
                if (sweepTargetStart >= rangeEnd) break;

                int[] vList = sweepVariants[g];
                for (int v : vList) {
                    int anchorPos = v >>> 2;
                    int xorDiff = v & 3;
                    for (int m = mStart; m < mEnd; m++) {
                        int aS = refOffsets[m];
                        int mL = matchLens[m];
                        if (anchorPos >= aS && anchorPos < aS + mL) {
                            int k = anchorPos - aS;
                            int strand = (strandBits[m >>> 3] >>> (m & 7)) & 1;
                            int absT = (strand == 0)
                                    ? (matchTargetStarts[m] + k)
                                    : (matchTargetStarts[m] + mL - 1 - k);
                            if (absT >= rangeStart && absT < rangeEnd) {
                                out[absT - rangeStart] ^= (byte) xorDiff;
                            }
                            break;
                        }
                    }
                }
            }
        }
        return out;
    }
}
