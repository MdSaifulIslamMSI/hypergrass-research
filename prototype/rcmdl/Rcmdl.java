import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.stream.Stream;

/**
 * RC-MDL research prototype.
 *
 * <p>Reference-conditioned exact/approximate match parsing with an adaptive range-coded
 * token stream. Two modes are supported:
 *
 * <ul>
 *   <li><b>R</b> (reference-conditioned): the target is parsed into panel copies (forward or
 *       reverse-complement, with substitutions), self copies, and literals. The archive stores
 *       only the target; the panel is external and pinned by name, length, and SHA-256.</li>
 *   <li><b>S</b> (self-contained): no panel. The same parser uses only self copies and
 *       literals, which is the honest single-genome baseline for this architecture.</li>
 * </ul>
 *
 * <p>This is a flat-sequence prototype. It does not read FASTA and it does not write the
 * production .hgs container. It is meant to answer one question: how many bits per base does
 * an MDL-style mosaic parse recover on a real pangenome panel, and is that materially below
 * existing collection codecs once the panel is excluded from the archive?
 */
public final class Rcmdl {

    private static final int MAGIC = 0x52434D44; // "RCMD"
    private static final int VERSION = 2;
    private static final int MIN_MATCH = 16;
    private static final int MIN_SCORE = 32;
    private static final int MAX_CHAIN = 32;
    private static final int MAX_MATCH = 1 << 20;
    private static final int MAX_MISMATCHES = 64;
    private static final int MAX_DIST = 1 << 24;
    private static final int MAX_LITERAL_RUN = 1 << 16;
    private static final int MAX_PANELS = 64;
    private static final int HASH_BITS = 24;
    private static final int HASH_SIZE = 1 << HASH_BITS;
    private static final int TOP = 1 << 24;
    private static final int PROB_INIT = 1 << 10;

    private static final byte[] COMP = new byte[256];

    static {
        for (int i = 0; i < 256; i++) {
            COMP[i] = (byte) i;
        }
        COMP['A'] = 'T';
        COMP['C'] = 'G';
        COMP['G'] = 'C';
        COMP['T'] = 'A';
        COMP['a'] = 't';
        COMP['c'] = 'g';
        COMP['g'] = 'c';
        COMP['t'] = 'a';
        COMP['N'] = 'N';
        COMP['n'] = 'n';
    }

    private Rcmdl() {}

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            usage();
            return;
        }
        switch (args[0]) {
            case "compress" -> compress(args);
            case "decompress" -> decompress(args);
            default -> usage();
        }
    }

    private static void usage() {
        System.err.println("""
                usage:
                  Rcmdl compress   --target FILE --out FILE [--panel DIR] [--panel-list FILE]
                                   [--panel-limit N] [--k 16] [--step 8] [--max-chain 32] [--no-approx]
                  Rcmdl decompress --archive FILE [--panel DIR] --out FILE
                """);
    }

    // ---------------------------------------------------------------- records

    private record Genome(String name, byte[] seq, byte[] sha) {
        Genome(String name, byte[] seq) {
            this(name, seq, sha256(seq));
        }
    }

    private static final class Candidate {
        final boolean self;
        final int genome;
        final int pos;
        final int orient;
        final int len;
        final int mismatches;
        final int[] mmOff;
        final byte[] mmBase;
        final int score;

        Candidate(
                boolean self,
                int genome,
                int pos,
                int orient,
                int len,
                int mismatches,
                int[] mmOff,
                byte[] mmBase) {
            this.self = self;
            this.genome = genome;
            this.pos = pos;
            this.orient = orient;
            this.len = len;
            this.mismatches = mismatches;
            this.mmOff = mmOff;
            this.mmBase = mmBase;
            this.score = len - 4 * mismatches;
        }
    }

    private static final class Stats {
        long targetBases;
        long archiveBytes;
        long panelBases;
        long panelMatchBases;
        long selfMatchBases;
        long literalBases;
        long panelTokens;
        long selfTokens;
        long literalRuns;
        long mismatches;
        int panelGenomes;
        double compressSeconds;
        double decompressSeconds;
    }

    // ---------------------------------------------------------------- panel and index

    private static final class MatchIndex {
        final byte[] text;
        final int[] starts;
        final int[] ends;
        final int k;
        final int[] head;
        final int[] next;
        final int[] pos;
        final byte[] src;
        final byte[] orient;
        int count;

        MatchIndex(byte[] text, int[] starts, int[] ends, int k, int maxSeeds) {
            this.text = text;
            this.starts = starts;
            this.ends = ends;
            this.k = k;
            this.head = new int[HASH_SIZE];
            Arrays.fill(head, -1);
            this.next = new int[maxSeeds];
            this.pos = new int[maxSeeds];
            this.src = new byte[maxSeeds];
            this.orient = new byte[maxSeeds];
        }

        void insert(int genome, int localPos, int orientation, long code) {
            int h = hash(code);
            src[count] = (byte) genome;
            pos[count] = localPos;
            orient[count] = (byte) orientation;
            next[count] = head[h];
            head[h] = count++;
        }
    }

    private static MatchIndex buildPanelIndex(List<Genome> panels, int k, int step) {
        long total = 0;
        long seeds = 0;
        for (Genome g : panels) {
            total += g.seq.length;
            seeds += (g.seq.length / step + 2L) * 2L;
        }
        if (total > Integer.MAX_VALUE - 8L || seeds > Integer.MAX_VALUE - 8L) {
            throw new IllegalArgumentException(
                    "panel too large for this prototype: bases=" + total + " seeds=" + seeds);
        }
        byte[] text = new byte[(int) total];
        int[] starts = new int[panels.size()];
        int[] ends = new int[panels.size()];
        int off = 0;
        for (int gi = 0; gi < panels.size(); gi++) {
            Genome g = panels.get(gi);
            starts[gi] = off;
            System.arraycopy(g.seq, 0, text, off, g.seq.length);
            off += g.seq.length;
            ends[gi] = off;
        }
        MatchIndex index = new MatchIndex(text, starts, ends, k, (int) seeds);
        for (int gi = 0; gi < panels.size(); gi++) {
            byte[] s = panels.get(gi).seq;
            for (int p = 0; p + k <= s.length; p += step) {
                long code = kmerCode(s, p, k);
                if (code < 0) {
                    continue;
                }
                index.insert(gi, p, 0, code);
                index.insert(gi, p, 1, rcCode(code, k));
            }
        }
        return index;
    }

    private static MatchIndex buildSelfIndex(byte[] target, int k, int step) {
        int maxSeeds = target.length / step + 2;
        MatchIndex index =
                new MatchIndex(target, new int[] {0}, new int[] {target.length}, k, maxSeeds);
        for (int p = 0; p + k <= target.length; p += step) {
            long code = kmerCode(target, p, k);
            if (code >= 0) {
                index.insert(0, p, 0, code);
            }
        }
        return index;
    }

    private static long kmerCode(byte[] s, int p, int k) {
        long code = 0;
        for (int i = 0; i < k; i++) {
            int b = base(s[p + i]);
            if (b < 0) {
                return -1;
            }
            code = (code << 2) | b;
        }
        return code;
    }

    private static long rcCode(long code, int k) {
        long rc = 0;
        for (int i = 0; i < k; i++) {
            rc = (rc << 2) | (3 - (code & 3));
            code >>= 2;
        }
        return rc;
    }

    private static int base(byte b) {
        return switch (b) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            case 'T' -> 3;
            default -> -1;
        };
    }

    private static int hash(long code) {
        long x = code * 0x9E3779B97F4A7C15L;
        x ^= x >>> 29;
        return (int) (x & (HASH_SIZE - 1));
    }

    // ---------------------------------------------------------------- match finding

    private static Candidate findBest(
            byte[] target, int ti, MatchIndex index, boolean allowApprox, boolean self) {
        int k = index.k;
        if (ti + k > target.length) {
            return null;
        }
        long code = kmerCode(target, ti, k);
        if (code < 0) {
            return null;
        }
        int c = index.head[hash(code)];
        int chain = 0;
        Candidate best = null;
        int n = target.length;
        while (c >= 0 && chain < MAX_CHAIN) {
            int gi = index.src[c] & 0xFF;
            int p = index.pos[c];
            int o = index.orient[c];
            boolean eligible = self ? (p < ti && ti - p <= MAX_DIST) : true;
            if (eligible && seedMatches(target, ti, index, gi, p, o)) {
                int limit = Math.min(MAX_MATCH, n - ti);
                Candidate cand = extend(target, ti, index, gi, p, o, allowApprox, self, limit);
                if (cand != null && (best == null || cand.score > best.score)) {
                    best = cand;
                }
            }
            c = index.next[c];
            chain++;
        }
        return best;
    }

    private static boolean seedMatches(
            byte[] target, int ti, MatchIndex index, int gi, int p, int o) {
        int gs = index.starts[gi];
        int ge = index.ends[gi];
        int k = index.k;
        if (o == 0) {
            if (p < 0 || p + k > ge - gs) {
                return false;
            }
            for (int x = 0; x < k; x++) {
                if (target[ti + x] != index.text[gs + p + x]) {
                    return false;
                }
            }
        } else {
            if (p < 0 || p + k > ge - gs) {
                return false;
            }
            for (int x = 0; x < k; x++) {
                if (target[ti + x] != COMP[index.text[gs + p + k - 1 - x] & 0xFF]) {
                    return false;
                }
            }
        }
        return true;
    }

    private static Candidate extend(
            byte[] target,
            int ti,
            MatchIndex index,
            int gi,
            int p,
            int o,
            boolean allowMismatch,
            boolean self,
            int max) {
        int k = index.k;
        int gs = index.starts[gi];
        int ge = index.ends[gi];
        int x = 0;
        int m = 0;
        int run = 0;
        int lastGood = -1;
        int[] off = new int[MAX_MISMATCHES];
        byte[] base = new byte[MAX_MISMATCHES];
        while (x < max) {
            int q = o == 0 ? gs + p + x : gs + p + k - 1 - x;
            if (q < gs || q >= ge) {
                break;
            }
            byte tb = target[ti + x];
            byte pb = index.text[q];
            byte expected = o == 0 ? pb : COMP[pb & 0xFF];
            if (tb == expected) {
                lastGood = x;
                run = 0;
                x++;
                continue;
            }
            if (!allowMismatch || self || m >= MAX_MISMATCHES) {
                break;
            }
            off[m] = x;
            base[m] = tb;
            m++;
            run++;
            if (run >= 4) {
                break;
            }
            x++;
        }
        int len = run >= 4 ? lastGood + 1 : x;
        if (len < MIN_MATCH) {
            return null;
        }
        int mm = 0;
        while (mm < m && off[mm] < len) {
            mm++;
        }
        if (mm * 8 > len) {
            return null;
        }
        int start = o == 0 ? p : p + k - len;
        if (start < 0 || start + len > ge - gs) {
            return null;
        }
        return new Candidate(self, self ? -1 : gi, start, o, len, mm, off, base);
    }

    // ---------------------------------------------------------------- encoder

    private static byte[] encode(
            byte[] target,
            List<Genome> panels,
            MatchIndex panelIndex,
            MatchIndex selfIndex,
            boolean allowApprox,
            Stats stats) {
        Encoder enc = new Encoder();
        Models m = new Models();
        int n = target.length;
        int i = 0;
        int prevType = 2;
        int pLenPrev = 0;
        int pPosPrev = 0;
        int sLenPrev = 0;
        int sDistPrev = 0;
        int runPrev = 0;
        int[] lastPos = new int[panels.size()];
        ByteArrayOutputStream literals = new ByteArrayOutputStream(1 << 16);
        stats.targetBases = n;
        stats.panelGenomes = panels.size();
        for (Genome g : panels) {
            stats.panelBases += g.seq.length;
        }

        while (i < n) {
            Candidate panel = panelIndex == null ? null : findBest(target, i, panelIndex, allowApprox, false);
            Candidate best = panel;
            if (selfIndex != null) {
                Candidate self = findBest(target, i, selfIndex, false, true);
                if (self != null && (best == null || self.score > best.score)) {
                    best = self;
                }
            }
            if (best != null && best.score >= MIN_SCORE) {
                if (literals.size() > 0) {
                    prevType = writeType(enc, m, prevType, 2);
                    Slot rs = writeSlot(enc, m.runSlot, m.runMant, runPrev, literals.size() - 1);
                    runPrev = rs.slot();
                    for (byte b : literals.toByteArray()) {
                        writeLiteralByte(enc, m, b);
                    }
                    stats.literalBases += literals.size();
                    stats.literalRuns++;
                    literals.reset();
                }
                if (best.self) {
                    prevType = writeType(enc, m, prevType, 1);
                    Slot ds = writeSlot(enc, m.sDistSlot, m.sDistMant, sDistPrev, i - best.pos - 1);
                    sDistPrev = ds.slot();
                    Slot ls = writeSlot(enc, m.sLenSlot, m.sLenMant, sLenPrev, best.len - MIN_MATCH);
                    sLenPrev = ls.slot();
                    stats.selfMatchBases += best.len;
                    stats.selfTokens++;
                } else {
                    prevType = writeType(enc, m, prevType, 0);
                    writeTree(enc, m.source, 0, 6, best.genome);
                    enc.bit(m.orient, 0, best.orient);
                    Slot ls = writeSlot(enc, m.pLenSlot, m.pLenMant, pLenPrev, best.len - MIN_MATCH);
                    pLenPrev = ls.slot();
                    int delta = best.pos - lastPos[best.genome];
                    lastPos[best.genome] = best.pos;
                    Slot ps = writeSlot(enc, m.pPosSlot, m.pPosMant, pPosPrev, zigzag(delta));
                    pPosPrev = ps.slot();
                    enc.bit(m.mmFlag, 0, best.mismatches == 0 ? 0 : 1);
                    if (best.mismatches > 0) {
                        writeSlot(enc, m.mmCountSlot, m.mmCountMant, 0, best.mismatches - 1);
                        int prev = -1;
                        for (int j = 0; j < best.mismatches; j++) {
                            int gap = best.mmOff[j] - prev - 1;
                            writeSlot(enc, m.mmGapSlot, m.mmGapMant, 0, gap);
                            writeBase(enc, m, best.mmBase[j]);
                            prev = best.mmOff[j];
                        }
                        stats.mismatches += best.mismatches;
                    }
                    stats.panelMatchBases += best.len;
                    stats.panelTokens++;
                }
                i += best.len;
                updatePrev(m, target, i);
            } else {
                literals.write(target[i] & 0xFF);
                i++;
                if (literals.size() >= MAX_LITERAL_RUN) {
                    prevType = writeType(enc, m, prevType, 2);
                    Slot rs = writeSlot(enc, m.runSlot, m.runMant, runPrev, literals.size() - 1);
                    runPrev = rs.slot();
                    for (byte b : literals.toByteArray()) {
                        writeLiteralByte(enc, m, b);
                    }
                    stats.literalBases += literals.size();
                    stats.literalRuns++;
                    literals.reset();
                }
            }
        }
        if (literals.size() > 0) {
            writeType(enc, m, prevType, 2);
            Slot rs = writeSlot(enc, m.runSlot, m.runMant, runPrev, literals.size() - 1);
            for (byte b : literals.toByteArray()) {
                writeLiteralByte(enc, m, b);
            }
            stats.literalBases += literals.size();
            stats.literalRuns++;
        }
        return enc.finish();
    }

    // ---------------------------------------------------------------- decoder

    private static byte[] decode(
            byte[] tokens, int targetLen, List<Genome> panels, Stats stats) {
        Decoder dec = new Decoder(tokens);
        Models m = new Models();
        byte[] out = new byte[targetLen];
        int o = 0;
        int prevType = 2;
        int pLenPrev = 0;
        int pPosPrev = 0;
        int sLenPrev = 0;
        int sDistPrev = 0;
        int runPrev = 0;
        int[] lastPos = new int[panels.size()];

        while (o < targetLen) {
            int type = readType(dec, m, prevType);
            if (type == 2) {
                Slot rs = readSlot(dec, m.runSlot, m.runMant, runPrev);
                runPrev = rs.slot();
                int runLen = rs.value() + 1;
                if (o + runLen > targetLen) {
                    throw new IllegalStateException("literal run overruns target");
                }
                for (int i = 0; i < runLen; i++) {
                    out[o++] = readLiteralByte(dec, m);
                }
                stats.literalBases += runLen;
                stats.literalRuns++;
            } else if (type == 1) {
                Slot ds = readSlot(dec, m.sDistSlot, m.sDistMant, sDistPrev);
                sDistPrev = ds.slot();
                Slot ls = readSlot(dec, m.sLenSlot, m.sLenMant, sLenPrev);
                sLenPrev = ls.slot();
                int dist = ds.value() + 1;
                int len = ls.value() + MIN_MATCH;
                if (dist <= 0 || dist > o || o + len > targetLen) {
                    throw new IllegalStateException(
                            "invalid self match dist=" + dist + " len=" + len + " at " + o);
                }
                for (int j = 0; j < len; j++) {
                    out[o + j] = out[o + j - dist];
                }
                o += len;
                updatePrev(m, out, o);
                stats.selfMatchBases += len;
                stats.selfTokens++;
            } else if (type == 0) {
                int gi = readTree(dec, m.source, 0, 6);
                if (gi < 0 || gi >= panels.size()) {
                    throw new IllegalStateException("invalid panel source " + gi);
                }
                int orient = dec.bit(m.orient, 0);
                Slot ls = readSlot(dec, m.pLenSlot, m.pLenMant, pLenPrev);
                pLenPrev = ls.slot();
                Slot ps = readSlot(dec, m.pPosSlot, m.pPosMant, pPosPrev);
                pPosPrev = ps.slot();
                int len = ls.value() + MIN_MATCH;
                int delta = unzigzag(ps.value());
                int pos = lastPos[gi] + delta;
                lastPos[gi] = pos;
                int hasMismatch = dec.bit(m.mmFlag, 0);
                int mism = 0;
                int[] offs = new int[0];
                byte[] bases = new byte[0];
                if (hasMismatch == 1) {
                    Slot cs = readSlot(dec, m.mmCountSlot, m.mmCountMant, 0);
                    mism = cs.value() + 1;
                    offs = new int[mism];
                    bases = new byte[mism];
                    int prev = -1;
                    for (int j = 0; j < mism; j++) {
                        Slot gs = readSlot(dec, m.mmGapSlot, m.mmGapMant, 0);
                        int off = prev + 1 + gs.value();
                        if (off < 0 || off >= len) {
                            throw new IllegalStateException("invalid mismatch offset " + off);
                        }
                        offs[j] = off;
                        bases[j] = readBase(dec, m);
                        prev = off;
                    }
                }
                byte[] panel = panels.get(gi).seq;
                if (pos < 0 || pos + len > panel.length) {
                    throw new IllegalStateException(
                            "invalid panel match pos=" + pos + " len=" + len + " panel=" + panel.length);
                }
                if (o + len > targetLen) {
                    throw new IllegalStateException("panel match overruns target");
                }
                if (orient == 0) {
                    System.arraycopy(panel, pos, out, o, len);
                } else {
                    for (int j = 0; j < len; j++) {
                        out[o + j] = COMP[panel[pos + len - 1 - j] & 0xFF];
                    }
                }
                for (int j = 0; j < mism; j++) {
                    out[o + offs[j]] = bases[j];
                }
                o += len;
                updatePrev(m, out, o);
                stats.panelMatchBases += len;
                stats.panelTokens++;
                stats.mismatches += mism;
            } else {
                throw new IllegalStateException("unknown token type " + type);
            }
            prevType = type;
        }
        return out;
    }

    // ---------------------------------------------------------------- models

    private static final class Models {
        final Probs type = new Probs(12);
        final Probs source = new Probs(64);
        final Probs orient = new Probs(2);
        final Probs pLenSlot = new Probs(1024);
        final Probs pLenMant = new Probs(1024);
        final Probs pPosSlot = new Probs(1024);
        final Probs pPosMant = new Probs(1024);
        final Probs sLenSlot = new Probs(1024);
        final Probs sLenMant = new Probs(1024);
        final Probs sDistSlot = new Probs(1024);
        final Probs sDistMant = new Probs(1024);
        final Probs mmFlag = new Probs(2);
        final Probs mmCountSlot = new Probs(1024);
        final Probs mmCountMant = new Probs(1024);
        final Probs mmGapSlot = new Probs(1024);
        final Probs mmGapMant = new Probs(1024);
        final Probs mmAcgt = new Probs(1);
        final Probs mmBase = new Probs(4);
        final Probs mmByte = new Probs(256);
        final Probs litAcgt = new Probs(25);
        final Probs litBase = new Probs(100);
        final Probs litByte = new Probs(256);
        final Probs runSlot = new Probs(1024);
        final Probs runMant = new Probs(1024);
        int prev1 = 0;
        int prev2 = 0;
    }

    private static final class Probs {
        final int[] p;

        Probs(int size) {
            p = new int[size];
            Arrays.fill(p, PROB_INIT);
        }

        void update(int idx, int bit) {
            int v = p[idx];
            p[idx] = bit == 0 ? v + ((2048 - v) >>> 5) : v - (v >>> 5);
        }
    }

    private static final class Encoder {
        final ByteArrayOutputStream out = new ByteArrayOutputStream(1 << 20);
        long low = 0;
        int range = -1;
        int cache = 0;
        long cacheSize = 1;

        void bit(Probs model, int idx, int bit) {
            int prob = model.p[idx];
            int bound = (range >>> 11) * prob;
            if (bit == 0) {
                range = bound;
            } else {
                low += bound & 0xFFFFFFFFL;
                range -= bound;
            }
            model.update(idx, bit);
            while ((range & 0xFFFFFFFFL) < TOP) {
                range <<= 8;
                shiftLow();
            }
        }

        private void shiftLow() {
            long lo = low;
            if ((int) (lo >>> 32) != 0 || lo < 0xFF000000L) {
                int carry = (int) (lo >>> 32);
                int temp = cache;
                do {
                    out.write((temp + carry) & 0xFF);
                    temp = 0xFF;
                } while (--cacheSize != 0);
                cache = (int) (lo >>> 24) & 0xFF;
            }
            cacheSize++;
            low = (lo & 0xFFFFFFL) << 8;
        }

        byte[] finish() {
            for (int i = 0; i < 5; i++) {
                shiftLow();
            }
            return out.toByteArray();
        }
    }

    private static final class Decoder {
        final byte[] src;
        int pos;
        int range = -1;
        int code;

        Decoder(byte[] src) {
            this.src = src;
            for (int i = 0; i < 5; i++) {
                code = (code << 8) | next();
            }
        }

        private int next() {
            return pos < src.length ? src[pos++] & 0xFF : 0;
        }

        int bit(Probs model, int idx) {
            int prob = model.p[idx];
            int bound = (range >>> 11) * prob;
            int b;
            if ((code & 0xFFFFFFFFL) < (bound & 0xFFFFFFFFL)) {
                range = bound;
                b = 0;
            } else {
                code -= bound;
                range -= bound;
                b = 1;
            }
            model.update(idx, b);
            while ((range & 0xFFFFFFFFL) < TOP) {
                range <<= 8;
                code = (code << 8) | next();
            }
            return b;
        }
    }

    // ---------------------------------------------------------------- coding primitives

    private record Slot(int value, int slot) {}

    private static Slot writeSlot(
            Encoder enc, Probs slotP, Probs mantP, int prevSlot, int value) {
        int slot = value == 0 ? 0 : 32 - Integer.numberOfLeadingZeros(value);
        if (slot > 31) {
            throw new IllegalStateException("slot overflow for " + value);
        }
        writeTree(enc, slotP, prevSlot * 32, 5, slot);
        if (slot > 1) {
            int mant = value - (1 << (slot - 1));
            for (int bit = slot - 2; bit >= 0; bit--) {
                enc.bit(mantP, slot * 32 + (slot - 2 - bit), (mant >>> bit) & 1);
            }
        }
        return new Slot(value, slot);
    }

    private static Slot readSlot(Decoder dec, Probs slotP, Probs mantP, int prevSlot) {
        int slot = readTree(dec, slotP, prevSlot * 32, 5);
        if (slot == 0) {
            return new Slot(0, 0);
        }
        int value = 1 << (slot - 1);
        for (int bit = slot - 2; bit >= 0; bit--) {
            int b = dec.bit(mantP, slot * 32 + (slot - 2 - bit));
            value |= b << bit;
        }
        return new Slot(value, slot);
    }

    private static void writeTree(Encoder enc, Probs p, int base, int bits, int value) {
        int node = 1;
        for (int bit = bits - 1; bit >= 0; bit--) {
            int b = (value >>> bit) & 1;
            enc.bit(p, base + node, b);
            node = (node << 1) | b;
        }
    }

    private static int readTree(Decoder dec, Probs p, int base, int bits) {
        int node = 1;
        for (int i = 0; i < bits; i++) {
            node = (node << 1) | dec.bit(p, base + node);
        }
        return node - (1 << bits);
    }

    private static int writeType(Encoder enc, Models m, int prev, int type) {
        int base = prev * 4;
        int node = 1;
        for (int bit = 1; bit >= 0; bit--) {
            int b = (type >>> bit) & 1;
            enc.bit(m.type, base + node, b);
            node = (node << 1) | b;
        }
        return type;
    }

    private static int readType(Decoder dec, Models m, int prev) {
        int base = prev * 4;
        int node = 1;
        for (int i = 0; i < 2; i++) {
            node = (node << 1) | dec.bit(m.type, base + node);
        }
        return node - 4;
    }

    private static void writeBase(Encoder enc, Models m, byte b) {
        int code = base(b);
        if (code >= 0) {
            enc.bit(m.mmAcgt, 0, 1);
            writeTree(enc, m.mmBase, 0, 2, code);
        } else {
            enc.bit(m.mmAcgt, 0, 0);
            writeTree(enc, m.mmByte, 0, 8, b & 0xFF);
        }
    }

    private static byte readBase(Decoder dec, Models m) {
        if (dec.bit(m.mmAcgt, 0) == 1) {
            return switch (readTree(dec, m.mmBase, 0, 2)) {
                case 0 -> (byte) 'A';
                case 1 -> (byte) 'C';
                case 2 -> (byte) 'G';
                case 3 -> (byte) 'T';
                default -> throw new IllegalStateException("bad base code");
            };
        }
        return (byte) readTree(dec, m.mmByte, 0, 8);
    }

    private static void writeLiteralByte(Encoder enc, Models m, byte b) {
        int code = base(b);
        int ctx = m.prev1 * 5 + m.prev2;
        if (code >= 0) {
            enc.bit(m.litAcgt, ctx, 1);
            writeTree(enc, m.litBase, ctx * 4, 2, code);
        } else {
            enc.bit(m.litAcgt, ctx, 0);
            writeTree(enc, m.litByte, 0, 8, b & 0xFF);
        }
        m.prev1 = m.prev2;
        m.prev2 = code >= 0 ? code : 4;
    }

    private static byte readLiteralByte(Decoder dec, Models m) {
        int ctx = m.prev1 * 5 + m.prev2;
        byte b;
        if (dec.bit(m.litAcgt, ctx) == 1) {
            int code = readTree(dec, m.litBase, ctx * 4, 2);
            b = switch (code) {
                case 0 -> (byte) 'A';
                case 1 -> (byte) 'C';
                case 2 -> (byte) 'G';
                case 3 -> (byte) 'T';
                default -> throw new IllegalStateException("bad literal base");
            };
            m.prev1 = m.prev2;
            m.prev2 = code;
            return b;
        }
        b = (byte) readTree(dec, m.litByte, 0, 8);
        m.prev1 = m.prev2;
        m.prev2 = 4;
        return b;
    }

    private static void updatePrev(Models m, byte[] data, int end) {
        int last = 0;
        int prev = 0;
        if (end >= 1) {
            int b = base(data[end - 1]);
            last = b >= 0 ? b : 4;
        }
        if (end >= 2) {
            int b = base(data[end - 2]);
            prev = b >= 0 ? b : 4;
        }
        m.prev1 = prev;
        m.prev2 = last;
    }

    private static int zigzag(int v) {
        return (v << 1) ^ (v >> 31);
    }

    private static int unzigzag(int u) {
        return (u >>> 1) ^ -(u & 1);
    }

    // ---------------------------------------------------------------- archive and CLI

    private static void compress(String[] args) throws Exception {
        long start = System.nanoTime();
        Path targetPath = null;
        Path panelPath = null;
        Path panelListPath = null;
        Path outPath = null;
        int k = 16;
        int step = 8;
        int panelLimit = 0;
        String exclude = null;
        boolean allowApprox = true;
        for (int i = 1; i < args.length; i++) {
            switch (args[i]) {
                case "--target" -> targetPath = Path.of(args[++i]);
                case "--panel" -> panelPath = Path.of(args[++i]);
                case "--panel-list" -> panelListPath = Path.of(args[++i]);
                case "--out" -> outPath = Path.of(args[++i]);
                case "--k" -> k = Integer.parseInt(args[++i]);
                case "--step" -> step = Integer.parseInt(args[++i]);
                case "--panel-limit" -> panelLimit = Integer.parseInt(args[++i]);
                case "--exclude" -> exclude = args[++i];
                case "--no-approx" -> allowApprox = false;
                default -> throw new IllegalArgumentException("unknown option " + args[i]);
            }
        }
        if (targetPath == null || outPath == null) {
            usage();
            throw new IllegalArgumentException("--target and --out are required");
        }
        if (k < 8 || k > 31) {
            throw new IllegalArgumentException("k must be 8..31");
        }
        if (step < 1 || step > 255) {
            throw new IllegalArgumentException("step must be 1..255");
        }
        if (panelLimit < 0 || panelLimit > MAX_PANELS) {
            throw new IllegalArgumentException("panel-limit must be 0.." + MAX_PANELS);
        }
        byte[] target = Files.readAllBytes(targetPath);
        List<Genome> panels = new ArrayList<>();
        if (panelListPath != null) {
            if (panelPath == null) {
                throw new IllegalArgumentException("--panel-list requires --panel");
            }
            List<String> names = new ArrayList<>();
            for (String line : Files.readAllLines(panelListPath)) {
                String name = line.trim();
                if (!name.isEmpty() && !name.startsWith("#")) {
                    names.add(name);
                }
            }
            if (panelLimit > 0 && panelLimit < names.size()) {
                names = names.subList(0, panelLimit);
            }
            for (String name : names) {
                Path p = panelPath.resolve(name);
                panels.add(new Genome(p.getFileName().toString(), Files.readAllBytes(p)));
            }
        } else if (panelPath != null) {
            int used = 0;
            for (Path p : panelFiles(panelPath, 0)) {
                if (exclude != null && p.getFileName().toString().equals(exclude)) {
                    continue;
                }
                if (panelLimit > 0 && used >= panelLimit) {
                    break;
                }
                panels.add(new Genome(p.getFileName().toString(), Files.readAllBytes(p)));
                used++;
            }
        }
        MatchIndex panelIndex = panels.isEmpty() ? null : buildPanelIndex(panels, k, step);
        MatchIndex selfIndex = buildSelfIndex(target, k, step);
        Stats stats = new Stats();
        byte[] tokens = encode(target, panels, panelIndex, selfIndex, allowApprox, stats);
        byte[] archive = buildArchive(target, panels, tokens, k, step, allowApprox ? 1 : 0);
        Files.write(outPath, archive);
        stats.archiveBytes = archive.length;
        stats.compressSeconds = (System.nanoTime() - start) / 1e9;
        printStats(stats);
    }

    private static void decompress(String[] args) throws Exception {
        long start = System.nanoTime();
        Path archivePath = null;
        Path panelPath = null;
        Path outPath = null;
        for (int i = 1; i < args.length; i++) {
            switch (args[i]) {
                case "--archive" -> archivePath = Path.of(args[++i]);
                case "--panel" -> panelPath = Path.of(args[++i]);
                case "--out" -> outPath = Path.of(args[++i]);
                default -> throw new IllegalArgumentException("unknown option " + args[i]);
            }
        }
        if (archivePath == null || outPath == null) {
            usage();
            throw new IllegalArgumentException("--archive and --out are required");
        }
        byte[] all = Files.readAllBytes(archivePath);
        java.io.DataInputStream in =
                new java.io.DataInputStream(new java.io.ByteArrayInputStream(all));
        if (in.readInt() != MAGIC) {
            throw new IllegalStateException("bad magic");
        }
        int version = in.readUnsignedByte();
        if (version != VERSION) {
            throw new IllegalStateException("unsupported version " + version);
        }
        int mode = in.readUnsignedByte();
        int targetLen = in.readInt();
        int panelCount = in.readUnsignedShort();
        List<Genome> panels = new ArrayList<>();
        for (int i = 0; i < panelCount; i++) {
            int nameLen = in.readUnsignedShort();
            byte[] nameBytes = new byte[nameLen];
            in.readFully(nameBytes);
            String name = new String(nameBytes, StandardCharsets.UTF_8);
            int len = in.readInt();
            byte[] hash = new byte[32];
            in.readFully(hash);
            if (panelPath == null) {
                throw new IllegalStateException("archive needs a panel but --panel was not given");
            }
            byte[] seq = Files.readAllBytes(panelPath.resolve(name));
            if (seq.length != len) {
                throw new IllegalStateException("panel length mismatch for " + name);
            }
            if (!Arrays.equals(sha256(seq), hash)) {
                throw new IllegalStateException("panel digest mismatch for " + name);
            }
            panels.add(new Genome(name, seq, hash));
        }
        byte[] targetHash = new byte[32];
        in.readFully(targetHash);
        int k = in.readUnsignedByte();
        int step = in.readUnsignedByte();
        int flags = in.readUnsignedByte();
        int tokenLen = in.readInt();
        byte[] tokens = new byte[tokenLen];
        in.readFully(tokens);
        Stats stats = new Stats();
        byte[] out = decode(tokens, targetLen, panels, stats);
        if (!Arrays.equals(sha256(out), targetHash)) {
            throw new IllegalStateException("target digest mismatch");
        }
        Files.write(outPath, out);
        stats.targetBases = targetLen;
        stats.archiveBytes = all.length;
        stats.panelGenomes = panels.size();
        for (Genome g : panels) {
            stats.panelBases += g.seq.length;
        }
        stats.decompressSeconds = (System.nanoTime() - start) / 1e9;
        System.err.println("decoded mode=" + mode + " k=" + k + " step=" + step + " flags=" + flags);
        printStats(stats);
    }

    private static List<Path> panelFiles(Path dir, int limit) throws IOException {
        try (Stream<Path> stream = Files.list(dir)) {
            List<Path> files = stream
                    .filter(p -> p.getFileName().toString().endsWith(".seq"))
                    .sorted()
                    .toList();
            if (limit > 0 && limit < files.size()) {
                return files.subList(0, limit);
            }
            return files;
        }
    }

    private static byte[] buildArchive(
            byte[] target, List<Genome> panels, byte[] tokens, int k, int step, int flags)
            throws IOException {
        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        DataOutputStream out = new DataOutputStream(buf);
        out.writeInt(MAGIC);
        out.writeByte(VERSION);
        out.writeByte(panels.isEmpty() ? 0 : 1);
        out.writeInt(target.length);
        out.writeShort(panels.size());
        for (Genome g : panels) {
            byte[] name = g.name.getBytes(StandardCharsets.UTF_8);
            out.writeShort(name.length);
            out.write(name);
            out.writeInt(g.seq.length);
            out.write(g.sha);
        }
        out.write(sha256(target));
        out.writeByte(k);
        out.writeByte(step);
        out.writeByte(flags);
        out.writeInt(tokens.length);
        out.write(tokens);
        out.flush();
        return buf.toByteArray();
    }

    private static void printStats(Stats s) {
        double bpb = s.targetBases == 0 ? 0.0 : s.archiveBytes * 8.0 / s.targetBases;
        System.out.printf(Locale.ROOT,
                "{\"target_bytes\":%d,\"archive_bytes\":%d,\"bpb\":%.6f,"
                        + "\"panel_genomes\":%d,\"panel_bases\":%d,\"panel_match_bases\":%d,"
                        + "\"self_match_bases\":%d,\"literal_bases\":%d,\"panel_tokens\":%d,"
                        + "\"self_tokens\":%d,\"literal_runs\":%d,\"mismatches\":%d,"
                        + "\"compress_seconds\":%.3f,\"decompress_seconds\":%.3f}%n",
                s.targetBases, s.archiveBytes, bpb, s.panelGenomes, s.panelBases,
                s.panelMatchBases, s.selfMatchBases, s.literalBases, s.panelTokens,
                s.selfTokens, s.literalRuns, s.mismatches,
                s.compressSeconds, s.decompressSeconds);
    }

    private static byte[] sha256(byte[] data) {
        try {
            return MessageDigest.getInstance("SHA-256").digest(data);
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }
    }
}
