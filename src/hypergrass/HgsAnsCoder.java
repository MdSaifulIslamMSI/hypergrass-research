package hypergrass;

import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

/**
 * Pure Java 25 LTS bit-exact counterpart to C23 {@code native/c_kernels/hgs_ans.c} and C++23
 * {@code hgs_cpp_ans_novelty_encode/decode}: 4-stripe interleaved Byte-rANS ({@code PROB_BITS=12},
 * {@code L=1<<23}) for DNA novelty streams (Order-0/2/4/6 + Law IV raw 2-bit ceiling) and
 * pangenome-mismatch + homopolymer + cycle-conditioned Phred quality scores.
 */
public final class HgsAnsCoder {

    private static final int ANS_PROB_BITS = 12;
    private static final int ANS_SCALE = 1 << ANS_PROB_BITS; // 4096
    private static final long ANS_L = 1L << 23; // 8,388,608
    private static final long ANS_XMAX_MUL = (ANS_L >>> ANS_PROB_BITS) << 8; // 524,288

    private static final int ANS_MODE_RAW_2BIT = 0;
    private static final int ANS_MODE_RANS4_DNA = 1;
    private static final int ANS_MODE_RAW_QUAL = 2;
    private static final int ANS_MODE_RANS4_QUAL = 3;

    public static final int FLAG_ANS4_NOVELTY = 0x10;
    public static final int ANS4_SEG_BASES = 1_000_000;

    private HgsAnsCoder() {}

    private static void normalizeCdf(int[] weights, int n, int[] outCum, int[] outFreq) {
        long tot = 0;
        int best = 0;
        int bestW = 0;
        for (int i = 0; i < n; i++) {
            int w = weights[i] + 1;
            tot += w;
            if (w > bestW) {
                bestW = w;
                best = i;
            }
        }
        int budget = ANS_SCALE - n;
        int sum = 0;
        for (int i = 0; i < n; i++) {
            int f = 1 + (int) (((long) (weights[i] + 1) * budget) / tot);
            outFreq[i] = f;
            sum += f;
        }
        outFreq[best] += (ANS_SCALE - sum);
        int c = 0;
        for (int i = 0; i < n; i++) {
            outCum[i] = c;
            c += outFreq[i];
        }
    }

    private static final class ByteRev {
        byte[] data = new byte[1024];
        int len = 0;

        void push(byte b) {
            if (len == data.length) {
                data = Arrays.copyOf(data, data.length << 1);
            }
            data[len++] = b;
        }
    }

    private static long encodeOne(long x, ByteRev rev, int cum, int freq) {
        long xmax = ANS_XMAX_MUL * freq;
        while (x >= xmax) {
            rev.push((byte) (x & 0xFFL));
            x >>>= 8;
        }
        return ((x / freq) << ANS_PROB_BITS) + (x % freq) + cum;
    }

    private static final class AnsDecCursor {
        long state;
        int pos;
        final int end;
        final byte[] buf;

        AnsDecCursor(long initState, int startPos, int endPos, byte[] buf) {
            this.state = initState;
            this.pos = startPos;
            this.end = endPos;
            this.buf = buf;
        }

        int step(int[] cums, int[] freqs, int nSyms) {
            long x = state;
            int slot = (int) (x & (ANS_SCALE - 1));
            int sym = 0;
            while (sym + 1 < nSyms && slot >= cums[sym + 1]) {
                sym++;
            }
            int f = freqs[sym];
            int c = cums[sym];
            x = f * (x >>> ANS_PROB_BITS) + (slot - c);
            while (x < ANS_L) {
                if (pos >= end) {
                    throw new IllegalStateException("Truncated rANS compressed stripe");
                }
                x = (x << 8) | (buf[pos++] & 0xFFL);
            }
            state = x;
            return sym;
        }
    }

    private static void writeU32LE(ByteOut out, long v) {
        out.writeByte((int) (v & 0xFFL));
        out.writeByte((int) ((v >>> 8) & 0xFFL));
        out.writeByte((int) ((v >>> 16) & 0xFFL));
        out.writeByte((int) ((v >>> 24) & 0xFFL));
    }

    private static long readU32LE(byte[] buf, int off) {
        return (buf[off] & 0xFFL)
                | ((buf[off + 1] & 0xFFL) << 8)
                | ((buf[off + 2] & 0xFFL) << 16)
                | ((buf[off + 3] & 0xFFL) << 24);
    }

    // ---------------- Adaptive 4-Stripe DNA Predictor (Order 0/2/4/6) ----------------

    private static final class DnaLaneModel {
        final int[][] o2 = new int[16][4];
        final int[][] o4 = new int[256][4];
        final int[][] o6 = new int[4096][4];
        int ctx2 = 0;
        int ctx4 = 0;
        int ctx6 = 0;

        void cdf(int[] cum, int[] freq, int[] w) {
            int[] p2 = o2[ctx2];
            int[] p4 = o4[ctx4];
            int[] p6 = o6[ctx6];
            for (int s = 0; s < 4; s++) {
                w[s] = p2[s] + (p4[s] << 2) + (p6[s] << 4);
            }
            normalizeCdf(w, 4, cum, freq);
        }

        void update(int b) {
            int[] p2 = o2[ctx2];
            int[] p4 = o4[ctx4];
            int[] p6 = o6[ctx6];
            if (++p2[b] > 2048) {
                for (int s = 0; s < 4; s++) p2[s] >>>= 1;
            }
            if (++p4[b] > 2048) {
                for (int s = 0; s < 4; s++) p4[s] >>>= 1;
            }
            if (++p6[b] > 2048) {
                for (int s = 0; s < 4; s++) p6[s] >>>= 1;
            }
            ctx2 = ((ctx2 << 2) | b) & 15;
            ctx4 = ((ctx4 << 2) | b) & 255;
            ctx6 = ((ctx6 << 2) | b) & 4095;
        }
    }

    private static byte[] packRaw2Bit(byte[] bases, int off, int len) {
        int nbytes = (len + 3) >>> 2;
        byte[] out = new byte[1 + nbytes];
        out[0] = (byte) ANS_MODE_RAW_2BIT;
        int i = 0;
        int op = 1;
        while (i + 4 <= len) {
            out[op++] =
                    (byte)
                            ((bases[off + i] & 3)
                                    | ((bases[off + i + 1] & 3) << 2)
                                    | ((bases[off + i + 2] & 3) << 4)
                                    | ((bases[off + i + 3] & 3) << 6));
            i += 4;
        }
        if (i < len) {
            int v = 0;
            int sh = 0;
            while (i < len) {
                v |= (bases[off + i++] & 3) << sh;
                sh += 2;
            }
            out[op] = (byte) v;
        }
        return out;
    }

    public static byte[] ans4EncodeDnaSlice(byte[] twoBit, int startOff, int len) {
        if (len == 0) {
            return new byte[] {(byte) ANS_MODE_RAW_2BIT};
        }
        if (len < 64) {
            return packRaw2Bit(twoBit, startOff, len);
        }
        int baseQ = len >>> 2;
        int rem = len & 3;
        int[] laneLen = new int[4];
        int[] laneOff = new int[4];
        int cur = 0;
        DnaLaneModel[] models = new DnaLaneModel[4];
        int[][] cumRecs = new int[4][];
        int[][] freqRecs = new int[4][];
        ByteRev[] rev = new ByteRev[4];
        for (int l = 0; l < 4; l++) {
            laneLen[l] = baseQ + (l < rem ? 1 : 0);
            laneOff[l] = startOff + cur;
            cur += laneLen[l];
            models[l] = new DnaLaneModel();
            cumRecs[l] = new int[laneLen[l]];
            freqRecs[l] = new int[laneLen[l]];
            rev[l] = new ByteRev();
        }

        int[] cum = new int[4];
        int[] freq = new int[4];
        int[] w = new int[4];
        for (int l = 0; l < 4; l++) {
            int sOff = laneOff[l];
            int n = laneLen[l];
            DnaLaneModel m = models[l];
            for (int i = 0; i < n; i++) {
                int b = twoBit[sOff + i] & 3;
                m.cdf(cum, freq, w);
                cumRecs[l][i] = cum[b];
                freqRecs[l][i] = freq[b];
                m.update(b);
            }
        }

        long[] state = {ANS_L, ANS_L, ANS_L, ANS_L};
        for (int l = 0; l < 4; l++) {
            int n = laneLen[l];
            for (int i = n; i > 0; i--) {
                state[l] = encodeOne(state[l], rev[l], cumRecs[l][i - 1], freqRecs[l][i - 1]);
            }
        }

        int rawSize = 1 + ((len + 3) >>> 2);
        int ransSize = 33;
        for (int l = 0; l < 4; l++) {
            ransSize += rev[l].len;
        }
        if (ransSize >= rawSize) {
            return packRaw2Bit(twoBit, startOff, len);
        }

        ByteOut out = new ByteOut(ransSize);
        out.writeByte(ANS_MODE_RANS4_DNA);
        for (int l = 0; l < 4; l++) writeU32LE(out, state[l]);
        for (int l = 0; l < 4; l++) writeU32LE(out, rev[l].len);
        for (int l = 0; l < 4; l++) {
            byte[] rd = rev[l].data;
            for (int i = rev[l].len; i > 0; i--) {
                out.writeByte(rd[i - 1]);
            }
        }
        return out.toByteArray();
    }

    public static byte[] ans4DecodeDnaSlice(byte[] payload, int off, int payloadLen, int expectedLen) {
        if (payloadLen < 1) {
            throw new IllegalStateException("Truncated ANS4 DNA block");
        }
        int mode = payload[off] & 0xFF;
        if (expectedLen == 0) {
            if (mode != ANS_MODE_RAW_2BIT) {
                throw new IllegalStateException("Unexpected mode on empty DNA block");
            }
            return new byte[0];
        }
        byte[] out = new byte[expectedLen];
        if (mode == ANS_MODE_RAW_2BIT) {
            int need = 1 + ((expectedLen + 3) >>> 2);
            if (payloadLen < need) {
                throw new IllegalStateException("Truncated RAW_2BIT ANS block");
            }
            for (int i = 0; i < expectedLen; i++) {
                out[i] = (byte) ((payload[off + 1 + (i >>> 2)] >>> ((i & 3) << 1)) & 3);
            }
            return out;
        }
        if (mode != ANS_MODE_RANS4_DNA || payloadLen < 33) {
            throw new IllegalStateException("Invalid ANS4 DNA payload header");
        }

        int baseQ = expectedLen >>> 2;
        int rem = expectedLen & 3;
        int[] laneLen = new int[4];
        int[] laneOff = new int[4];
        int cur = 0;
        DnaLaneModel[] models = new DnaLaneModel[4];
        AnsDecCursor[] decs = new AnsDecCursor[4];
        int pos = off + 33;
        int endBound = off + payloadLen;
        for (int l = 0; l < 4; l++) {
            laneLen[l] = baseQ + (l < rem ? 1 : 0);
            laneOff[l] = cur;
            cur += laneLen[l];
            models[l] = new DnaLaneModel();
            long st = readU32LE(payload, off + 1 + l * 4);
            int bl = (int) readU32LE(payload, off + 17 + l * 4);
            if (bl < 0 || pos + bl > endBound) {
                throw new IllegalStateException("Corrupt ANS4 stripe lengths");
            }
            decs[l] = new AnsDecCursor(st, pos, pos + bl, payload);
            pos += bl;
        }

        int maxQ = laneLen[0];
        int[][] cum = new int[4][4];
        int[][] freq = new int[4][4];
        int[] w = new int[4];
        for (int k = 0; k < maxQ; k++) {
            for (int l = 0; l < 4; l++) {
                if (k < laneLen[l]) models[l].cdf(cum[l], freq[l], w);
            }
            for (int l = 0; l < 4; l++) {
                if (k < laneLen[l]) {
                    int sym = decs[l].step(cum[l], freq[l], 4);
                    out[laneOff[l] + k] = (byte) sym;
                    models[l].update(sym);
                }
            }
        }
        for (int l = 0; l < 4; l++) {
            if (decs[l].state != ANS_L || decs[l].pos != decs[l].end) {
                throw new IllegalStateException("rANS DNA stripe final state mismatch on lane " + l);
            }
        }
        return out;
    }

    public static byte[] ansNoveltyEncode(byte[] twoBitBases, int threads) {
        int len = twoBitBases.length;
        if (len == 0) {
            return new byte[0];
        }
        int totalSegs = (len + ANS4_SEG_BASES - 1) / ANS4_SEG_BASES;
        int[] segStarts = new int[totalSegs];
        int[] segLens = new int[totalSegs];
        int cur = 0;
        for (int i = 0; i < totalSegs; i++) {
            segStarts[i] = cur;
            int chunk = Math.min(ANS4_SEG_BASES, len - cur);
            segLens[i] = chunk;
            cur += chunk;
        }
        byte[][] chunks = new byte[totalSegs][];
        int workers = Math.min(Math.max(1, threads), totalSegs);
        if (workers <= 1) {
            for (int i = 0; i < totalSegs; i++) {
                chunks[i] = ans4EncodeDnaSlice(twoBitBases, segStarts[i], segLens[i]);
            }
        } else {
            ForkJoinPool pool = new ForkJoinPool(workers);
            try {
                pool.submit(
                                () ->
                                        IntStream.range(0, totalSegs)
                                                .parallel()
                                                .forEach(
                                                        i ->
                                                                chunks[i] =
                                                                        ans4EncodeDnaSlice(
                                                                                twoBitBases,
                                                                                segStarts[i],
                                                                                segLens[i])))
                        .join();
            } finally {
                pool.shutdown();
            }
        }

        int dataBytes = 0;
        for (byte[] c : chunks) dataBytes += c.length;
        ByteOut out = new ByteOut(32 + totalSegs * 8 + dataBytes);
        out.writeByte(FLAG_ANS4_NOVELTY);
        out.writeVarInt(totalSegs);
        for (int i = 0; i < totalSegs; i++) {
            out.writeVarInt(segLens[i]);
            out.writeVarInt(chunks[i].length);
        }
        for (int i = 0; i < totalSegs; i++) {
            out.writeBytes(chunks[i]);
        }
        return out.toByteArray();
    }

    public static byte[] ansNoveltyDecode(byte[] payload, int expectedLen, int threads) {
        if (expectedLen == 0) {
            if (payload.length != 0) {
                throw new IllegalStateException("Non-empty ANS payload for empty novelty stream");
            }
            return new byte[0];
        }
        ByteIn in = new ByteIn(payload);
        int flags = in.readByte() & 0xFF;
        if (flags != FLAG_ANS4_NOVELTY) {
            throw new IllegalStateException("Expected FLAG_ANS4_NOVELTY (0x10), got 0x" + Integer.toHexString(flags));
        }
        int totalSegs = in.readVarInt();
        if (totalSegs <= 0) {
            throw new IllegalStateException("Invalid ANS4 segment count");
        }
        int[] segStarts = new int[totalSegs];
        int[] segLens = new int[totalSegs];
        int[] compLens = new int[totalSegs];
        int sumBases = 0;
        for (int i = 0; i < totalSegs; i++) {
            int sl = in.readVarInt();
            int cl = in.readVarInt();
            segStarts[i] = sumBases;
            segLens[i] = sl;
            compLens[i] = cl;
            sumBases += sl;
        }
        if (sumBases != expectedLen) {
            throw new IllegalStateException("ANS4 decoded base length mismatch");
        }
        int[] compOffs = new int[totalSegs];
        int curOff = in.position();
        for (int i = 0; i < totalSegs; i++) {
            compOffs[i] = curOff;
            curOff += compLens[i];
        }
        if (curOff != payload.length) {
            throw new IllegalStateException("ANS4 payload trailing or truncated bytes");
        }

        byte[] dst = new byte[expectedLen];
        int workers = Math.min(Math.max(1, threads), totalSegs);
        if (workers <= 1) {
            for (int i = 0; i < totalSegs; i++) {
                byte[] seg = ans4DecodeDnaSlice(payload, compOffs[i], compLens[i], segLens[i]);
                System.arraycopy(seg, 0, dst, segStarts[i], segLens[i]);
            }
        } else {
            ForkJoinPool pool = new ForkJoinPool(workers);
            try {
                pool.submit(
                                () ->
                                        IntStream.range(0, totalSegs)
                                                .parallel()
                                                .forEach(
                                                        i -> {
                                                            byte[] seg =
                                                                    ans4DecodeDnaSlice(
                                                                            payload,
                                                                            compOffs[i],
                                                                            compLens[i],
                                                                            segLens[i]);
                                                            System.arraycopy(
                                                                    seg,
                                                                    0,
                                                                    dst,
                                                                    segStarts[i],
                                                                    segLens[i]);
                                                        }))
                        .join();
            } finally {
                pool.shutdown();
            }
        }
        return dst;
    }

    // ---------------- Pangenome-Conditioned Phred Quality 4-Stripe Coder ----------------

    public static byte[] computeHomoBins(byte[] twoBitBases) {
        int len = twoBitBases.length;
        byte[] bins = new byte[len];
        if (len == 0) {
            return bins;
        }
        int run = 1;
        int prev = twoBitBases[0] & 3;
        bins[0] = 0;
        for (int i = 1; i < len; i++) {
            int b = twoBitBases[i] & 3;
            if (b == prev) {
                run++;
            } else {
                run = 1;
                prev = b;
            }
            bins[i] = (byte) ((run == 1) ? 0 : (run <= 3 ? 1 : (run <= 6 ? 2 : 3)));
        }
        return bins;
    }

    private static final class QualModel {
        final int[][] same = new int[128][2];
        final int[][] coarse = new int[128][9];
        final int[][] fine = new int[128][8];
    }

    private static int qualCtx(int prevQ, int mismatch, int homoBin, int pos, int readLenHint) {
        int pqBin = (prevQ >= 33 && prevQ <= 96) ? ((prevQ - 33) >>> 4) : 3;
        int cycBin = (readLenHint > 0) ? (((pos % readLenHint) * 2) / readLenHint) : 0;
        return ((mismatch & 1) << 5) | ((homoBin & 3) << 3) | (pqBin << 1) | (cycBin & 1);
    }

    private static int qualFineCtx(int cBin, int mismatch, int prevQ) {
        int pq = (prevQ >= 33 && prevQ <= 96) ? ((prevQ - 33) >>> 5) : 0;
        return ((cBin & 7) << 2) | ((mismatch & 1) << 1) | (pq & 1);
    }

    private static final class SymRecVec {
        int[] cums = new int[4096];
        int[] freqs = new int[4096];
        int len = 0;

        void push(int cum, int freq) {
            if (len == cums.length) {
                cums = Arrays.copyOf(cums, cums.length << 1);
                freqs = Arrays.copyOf(freqs, freqs.length << 1);
            }
            cums[len] = cum;
            freqs[len] = freq;
            len++;
        }
    }

    public static byte[] ans4EncodeQual(
            byte[] qualsRaw, byte[] mismatchMask, byte[] homoBins, int readLenHint) {
        int len = qualsRaw.length;
        if (len == 0) {
            return new byte[] {(byte) ANS_MODE_RAW_QUAL};
        }
        int[] u16Cum = new int[16];
        int[] u16Freq = new int[16];
        normalizeCdf(new int[16], 16, u16Cum, u16Freq);

        int baseQ = len >>> 2;
        int rem = len & 3;
        int[] laneLen = new int[4];
        int[] laneOff = new int[4];
        int cur = 0;
        QualModel[] models = new QualModel[4];
        SymRecVec[] recs = new SymRecVec[4];
        ByteRev[] rev = new ByteRev[4];
        for (int l = 0; l < 4; l++) {
            laneLen[l] = baseQ + (l < rem ? 1 : 0);
            laneOff[l] = cur;
            cur += laneLen[l];
            models[l] = new QualModel();
            recs[l] = new SymRecVec();
            rev[l] = new ByteRev();
        }

        int[] scum = new int[2];
        int[] sfreq = new int[2];
        int[] ccum = new int[9];
        int[] cfreq = new int[9];
        int[] fcum = new int[8];
        int[] ffreq = new int[8];

        for (int l = 0; l < 4; l++) {
            int prevQ = 66;
            QualModel m = models[l];
            for (int k = 0; k < laneLen[l]; k++) {
                int gpos = laneOff[l] + k;
                int q = qualsRaw[gpos] & 0xFF;
                int mm = mismatchMask != null ? (mismatchMask[gpos] & 1) : 0;
                int hb = homoBins != null ? (homoBins[gpos] & 3) : 0;
                int ctx = qualCtx(prevQ, mm, hb, gpos, readLenHint);

                normalizeCdf(m.same[ctx], 2, scum, sfreq);
                int isSame = (q == prevQ) ? 1 : 0;
                recs[l].push(scum[isSame], sfreq[isSame]);
                if (++m.same[ctx][isSame] > 4096) {
                    m.same[ctx][0] >>>= 1;
                    m.same[ctx][1] >>>= 1;
                }

                if (isSame == 0) {
                    normalizeCdf(m.coarse[ctx], 9, ccum, cfreq);
                    if (q >= 33 && q <= 96) {
                        int rel = q - 33;
                        int cbin = rel >>> 3;
                        int fbin = rel & 7;
                        recs[l].push(ccum[cbin], cfreq[cbin]);
                        if (++m.coarse[ctx][cbin] > 4096) {
                            for (int c = 0; c < 9; c++) m.coarse[ctx][c] >>>= 1;
                        }
                        int fctx = qualFineCtx(cbin, mm, prevQ);
                        normalizeCdf(m.fine[fctx], 8, fcum, ffreq);
                        recs[l].push(fcum[fbin], ffreq[fbin]);
                        if (++m.fine[fctx][fbin] > 4096) {
                            for (int f = 0; f < 8; f++) m.fine[fctx][f] >>>= 1;
                        }
                    } else {
                        recs[l].push(ccum[8], cfreq[8]);
                        if (++m.coarse[ctx][8] > 4096) {
                            for (int c = 0; c < 9; c++) m.coarse[ctx][c] >>>= 1;
                        }
                        int hi = (q >>> 4) & 0xF;
                        int lo = q & 0xF;
                        recs[l].push(u16Cum[hi], u16Freq[hi]);
                        recs[l].push(u16Cum[lo], u16Freq[lo]);
                    }
                }
                prevQ = q;
            }
        }

        long[] state = {ANS_L, ANS_L, ANS_L, ANS_L};
        for (int l = 0; l < 4; l++) {
            SymRecVec rv = recs[l];
            for (int i = rv.len; i > 0; i--) {
                state[l] = encodeOne(state[l], rev[l], rv.cums[i - 1], rv.freqs[i - 1]);
            }
        }

        int ransSize = 33;
        for (int l = 0; l < 4; l++) ransSize += rev[l].len;
        if (ransSize >= len + 1) {
            byte[] rawOut = new byte[1 + len];
            rawOut[0] = (byte) ANS_MODE_RAW_QUAL;
            System.arraycopy(qualsRaw, 0, rawOut, 1, len);
            return rawOut;
        }

        ByteOut out = new ByteOut(ransSize);
        out.writeByte(ANS_MODE_RANS4_QUAL);
        for (int l = 0; l < 4; l++) writeU32LE(out, state[l]);
        for (int l = 0; l < 4; l++) writeU32LE(out, rev[l].len);
        for (int l = 0; l < 4; l++) {
            byte[] rd = rev[l].data;
            for (int i = rev[l].len; i > 0; i--) {
                out.writeByte(rd[i - 1]);
            }
        }
        return out.toByteArray();
    }

    public static byte[] ans4DecodeQual(
            byte[] payload, byte[] mismatchMask, byte[] homoBins, int expectedLen, int readLenHint) {
        if (payload.length < 1) {
            throw new IllegalStateException("Truncated ANS4 quality block");
        }
        int mode = payload[0] & 0xFF;
        if (expectedLen == 0) {
            if (mode != ANS_MODE_RAW_QUAL) {
                throw new IllegalStateException("Invalid mode on empty quality stream");
            }
            return new byte[0];
        }
        byte[] out = new byte[expectedLen];
        if (mode == ANS_MODE_RAW_QUAL) {
            if (payload.length < 1 + expectedLen) {
                throw new IllegalStateException("Truncated RAW_QUAL block");
            }
            System.arraycopy(payload, 1, out, 0, expectedLen);
            return out;
        }
        if (mode != ANS_MODE_RANS4_QUAL || payload.length < 33) {
            throw new IllegalStateException("Invalid ANS4 quality header");
        }

        int[] u16Cum = new int[16];
        int[] u16Freq = new int[16];
        normalizeCdf(new int[16], 16, u16Cum, u16Freq);

        int baseQ = expectedLen >>> 2;
        int rem = expectedLen & 3;
        int[] laneLen = new int[4];
        int[] laneOff = new int[4];
        int cur = 0;
        QualModel[] models = new QualModel[4];
        AnsDecCursor[] decs = new AnsDecCursor[4];
        int pos = 33;
        for (int l = 0; l < 4; l++) {
            laneLen[l] = baseQ + (l < rem ? 1 : 0);
            laneOff[l] = cur;
            cur += laneLen[l];
            models[l] = new QualModel();
            long st = readU32LE(payload, 1 + l * 4);
            int bl = (int) readU32LE(payload, 17 + l * 4);
            if (bl < 0 || pos + bl > payload.length) {
                throw new IllegalStateException("Corrupt ANS4 quality stripe lengths");
            }
            decs[l] = new AnsDecCursor(st, pos, pos + bl, payload);
            pos += bl;
        }

        int[] scum = new int[2];
        int[] sfreq = new int[2];
        int[] ccum = new int[9];
        int[] cfreq = new int[9];
        int[] fcum = new int[8];
        int[] ffreq = new int[8];

        for (int l = 0; l < 4; l++) {
            int prevQ = 66;
            QualModel m = models[l];
            for (int k = 0; k < laneLen[l]; k++) {
                int gpos = laneOff[l] + k;
                int mm = mismatchMask != null ? (mismatchMask[gpos] & 1) : 0;
                int hb = homoBins != null ? (homoBins[gpos] & 3) : 0;
                int ctx = qualCtx(prevQ, mm, hb, gpos, readLenHint);

                normalizeCdf(m.same[ctx], 2, scum, sfreq);
                int isSame = decs[l].step(scum, sfreq, 2);
                if (++m.same[ctx][isSame] > 4096) {
                    m.same[ctx][0] >>>= 1;
                    m.same[ctx][1] >>>= 1;
                }

                int q = prevQ;
                if (isSame == 0) {
                    normalizeCdf(m.coarse[ctx], 9, ccum, cfreq);
                    int cbin = decs[l].step(ccum, cfreq, 9);
                    if (++m.coarse[ctx][cbin] > 4096) {
                        for (int c = 0; c < 9; c++) m.coarse[ctx][c] >>>= 1;
                    }
                    if (cbin < 8) {
                        int fctx = qualFineCtx(cbin, mm, prevQ);
                        normalizeCdf(m.fine[fctx], 8, fcum, ffreq);
                        int fbin = decs[l].step(fcum, ffreq, 8);
                        if (++m.fine[fctx][fbin] > 4096) {
                            for (int f = 0; f < 8; f++) m.fine[fctx][f] >>>= 1;
                        }
                        q = 33 + (cbin << 3) + fbin;
                    } else {
                        int hi = decs[l].step(u16Cum, u16Freq, 16);
                        int lo = decs[l].step(u16Cum, u16Freq, 16);
                        q = (hi << 4) | lo;
                    }
                }
                out[gpos] = (byte) q;
                prevQ = q;
            }
        }
        for (int l = 0; l < 4; l++) {
            if (decs[l].state != ANS_L || decs[l].pos != decs[l].end) {
                throw new IllegalStateException("rANS quality stripe final state mismatch on lane " + l);
            }
        }
        return out;
    }
}
