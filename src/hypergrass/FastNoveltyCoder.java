package hypergrass;

import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

/**
 * Blazing-fast, 100% deterministic integer SWAR Chargaff-Symmetric Multi-Order + Codon-Frame +
 * Wobble/Gapped Spaced-Seed + Adaptive Pan-Superstream Sparse-12 + 6-Head Multi-Scale Bidirectional
 * Kimura-XOR Repeat Context Mixer with Secondary Symbol Estimation (SSE), Warm-Primed Parallel
 * Segmentation, and a strict Raw 2-Bit Incompressible Safety Bypass (Law IV &amp; Law V).
 *
 * <p>Automatically adapts capacity to information scale (Single-Genome Anti-Overfitting Law):
 *
 * <ul>
 *   <li><b>Single-Genome / Bounded Profile ({@code len <= 6,000,000} or {@code boundedMem}):</b>
 *       100% L2/L3-cache-resident dense SWAR {@code Order-0..7} (87 KiB) + {@code Codon-1/3/5} +
 *       {@code Wobble-8} + {@code Gapped-8} + 6 Kimura-XOR Repeat Heads (~1.1 MiB/worker), avoiding
 *       high-order {@code Order-8/9/Sparse-12} single-chromosome overfitting while running 2.5x
 *       faster.
 *   <li><b>Pan-Genome Superstream Profile ({@code len > 6,000,000} and {@code !boundedMem}):</b>
 *       Activates {@code Order-8}, {@code Order-9}, and {@code Sparse-12} across multi-megabase
 *       accessory gene superstreams.
 * </ul>
 */
public final class FastNoveltyCoder {

    private static final int PROB_BITS = 15;
    private static final int TOT = 1 << PROB_BITS; // 32768
    private static final long TOP_VALUE = 1L << 24;

    private static final int P_BITS = 14;
    private static final int P_SCALE = 1 << P_BITS; // 16384

    private static final int PAN_SUPER_THRESHOLD = 6_000_000;
    private static final int TURBO_THRESHOLD = 8_000_000;
    private static final int TURBO_SEG_BASES = 4_000_000;
    private static final int FLAG_TURBO_NOVELTY = 4;

    private static final byte SEG_RAW_2BIT = 0;
    private static final byte SEG_CRM = 1;

    // Pure-integer natural log LUT: LN_LUT[p] = round(4096 * ln(p / 16384)) for p in [1..16384]
    private static final int[] LN_LUT = new int[P_SCALE + 1];
    // Pure-integer exp LUT: EXP_NEG_LUT[x] = round(65536 * exp(-x / 4096)) for x in [0..65536]
    private static final int EXP_MAX_X = 65536;
    private static final int[] EXP_NEG_LUT = new int[EXP_MAX_X + 1];
    // Precomputed exact reciprocals (1 << 30) / d for d in [1..65599]
    private static final int[] RECIP = new int[65600];

    static {
        initTables();
    }

    private static void initTables() {
        final long LN2_Q16 = 45426L;
        for (int p = 1; p <= P_SCALE; p++) {
            int intLog2 = 0;
            long y = p;
            while (y < P_SCALE) {
                y <<= 1;
                intLog2--;
            }
            int fracBits = 0;
            for (int b = 15; b >= 0; b--) {
                y = (y * y) >>> P_BITS;
                if (y >= (P_SCALE << 1)) {
                    fracBits |= (1 << b);
                    y >>>= 1;
                }
            }
            long log2Q16 = ((long) intLog2 << 16) + fracBits;
            LN_LUT[p] = (int) ((log2Q16 * LN2_Q16) >> 20);
        }
        LN_LUT[0] = LN_LUT[1] - 2839;

        final long INV_LN2_Q16 = 94548L;
        for (int x = 0; x <= EXP_MAX_X; x++) {
            long uQ16 = ((long) x * INV_LN2_Q16) >>> 12;
            int k = (int) (uQ16 >>> 16);
            if (k >= 16) {
                EXP_NEG_LUT[x] = 1;
                continue;
            }
            long r = uQ16 & 0xFFFFL;
            long s = (r * LN2_Q16) >>> 16;
            long s2 = (s * s) >>> 16;
            long s3 = (s2 * s) >>> 16;
            long s4 = (s3 * s) >>> 16;
            long s5 = (s4 * s) >>> 16;
            long expNegS = 65536L - s + (s2 >>> 1) - (s3 / 6L) + (s4 / 24L) - (s5 / 120L);
            EXP_NEG_LUT[x] = Math.max(1, (int) (expNegS >>> k));
        }

        for (int d = 1; d < RECIP.length; d++) {
            RECIP[d] = (int) ((1L << 30) / d);
        }
    }

    private FastNoveltyCoder() {}

    static int divRoundSmall(long num, int denom) {
        long n = num + (denom >> 1);
        if (denom < RECIP.length) {
            long q = (n * (RECIP[denom] & 0xFFFFFFFFL)) >>> 30;
            while (n - q * denom >= denom) {
                q++;
            }
            return (int) q;
        }
        return (int) (n / denom);
    }

    /**
     * Encodes {@code bases[0..len)} (values {@code 0..3}) with optional warm-primed parallelism and
     * per-segment incompressible raw-2-bit safety cap.
     */
    public static byte[] encode(byte[] bases, int len, int threads, boolean boundedMem) {
        if (len == 0) {
            ByteOut empty = new ByteOut(4);
            empty.writeVarInt(0);
            return empty.toByteArray();
        }

        int turboCutoff = boundedMem ? TURBO_THRESHOLD : 64_000_000;
        if (len >= turboCutoff) {
            return encodeTurbo(bases, len, threads, boundedMem);
        }

        boolean panSuper = len > PAN_SUPER_THRESHOLD;
        int modelFlags = (boundedMem ? 1 : 0) | (panSuper ? 2 : 0);

        int warmBases = panSuper ? 131_072 : (len >= 1_500_000 ? 393_216 : 65_536);
        int minTailBases = panSuper ? 524_288 : (len >= 1_500_000 ? 2_000_000 : 131_072);

        int tail = len - warmBases;
        int maxWorkers = panSuper ? Math.min(threads, 8) : threads;
        int workerCount =
                (maxWorkers <= 1 || tail < minTailBases)
                        ? 0
                        : Math.min(maxWorkers, tail / minTailBases);

        if (workerCount <= 1) {
            FastCrmModel model = new FastCrmModel(null, boundedMem, panSuper, 0, len, 1);
            byte[] seg = encodeSegmentWithSafety(bases, 0, len, model);
            ByteOut out = new ByteOut(seg.length + 16);
            out.writeVarInt(1);
            out.writeByte(modelFlags);
            out.writeVarInt(len);
            out.writeVarInt(seg.length);
            out.writeBytes(seg);
            return out.toByteArray();
        }

        // Step 1: Encode warm-up prefix [0 .. warmBases) once and snapshot model + prefix repeat rings
        int totalSegs = 1 + workerCount;
        FastCrmModel warmModel =
                new FastCrmModel(null, boundedMem, panSuper, warmBases, len, totalSegs);
        byte[] warmSeg = encodeSegmentWithSafety(bases, 0, warmBases, warmModel);
        WarmSnapshot snapshot = warmModel.snapshot();
        warmModel = null; // Free warm-up model before spawning parallel segment workers!

        // Step 2: Partition remaining [warmBases .. len) across workerCount parallel chunks
        int[] segLens = new int[totalSegs];
        int[] segStarts = new int[totalSegs];
        segLens[0] = warmBases;
        segStarts[0] = 0;

        int cursor = warmBases;
        for (int w = 0; w < workerCount; w++) {
            int sLen = tail / workerCount + (w < (tail % workerCount) ? 1 : 0);
            segStarts[w + 1] = cursor;
            segLens[w + 1] = sLen;
            cursor += sLen;
        }

        byte[][] encoded = new byte[totalSegs][];
        encoded[0] = warmSeg;

        final int finalWarmBases = warmBases;
        int poolThreads = boundedMem ? Math.min(workerCount, 6) : workerCount;
        ForkJoinPool pool = new ForkJoinPool(poolThreads);
        try {
            pool.submit(
                            () ->
                                    IntStream.range(0, workerCount)
                                            .parallel()
                                            .forEach(
                                                    w -> {
                                                        int idx = w + 1;
                                                        FastCrmModel m =
                                                                new FastCrmModel(
                                                                        snapshot,
                                                                        boundedMem,
                                                                        panSuper,
                                                                        finalWarmBases,
                                                                        len,
                                                                        totalSegs);
                                                        encoded[idx] =
                                                                encodeSegmentWithSafety(
                                                                        bases,
                                                                        segStarts[idx],
                                                                        segLens[idx],
                                                                        m);
                                                    }))
                    .join();
        } finally {
            pool.shutdown();
        }

        int totalComp = 0;
        for (byte[] e : encoded) {
            totalComp += e.length;
        }
        ByteOut out = new ByteOut(totalComp + totalSegs * 8 + 16);
        out.writeVarInt(totalSegs);
        out.writeByte(modelFlags);
        for (int i = 0; i < totalSegs; i++) {
            out.writeVarInt(segLens[i]);
        }
        for (int i = 0; i < totalSegs; i++) {
            out.writeVarInt(encoded[i].length);
        }
        for (int i = 0; i < totalSegs; i++) {
            out.writeBytes(encoded[i]);
        }
        return out.toByteArray();
    }

    /** Decodes {@code expectedLen} 2-bit bases produced by {@link #encode}. */
    public static byte[] decode(ByteIn in, int expectedLen, int threads, boolean boundedMem) {
        int totalSegs = in.readVarInt();
        if (expectedLen == 0) {
            if (totalSegs != 0) {
                throw new IllegalStateException("Expected 0 novelty segments for empty stream");
            }
            return new byte[0];
        }
        if (totalSegs <= 0 || totalSegs > expectedLen / 1024 + 2) {
            throw new IllegalStateException("Malformed novelty segment count: " + totalSegs);
        }
        int modelFlags = in.readByte() & 0xFF;
        if ((modelFlags & ~7) != 0) {
            throw new IllegalStateException("Unsupported novelty model flags: " + modelFlags);
        }
        boolean encodedBounded = (modelFlags & 1) != 0;
        boolean encodedPanSuper = (modelFlags & 2) != 0;
        boolean encodedTurbo = (modelFlags & FLAG_TURBO_NOVELTY) != 0;

        int[] segLens = new int[totalSegs];
        int sumLen = 0;
        for (int i = 0; i < totalSegs; i++) {
            segLens[i] = in.readVarInt();
            if (segLens[i] < 0 || sumLen + segLens[i] > expectedLen) {
                throw new IllegalStateException("Invalid novelty segment length");
            }
            sumLen += segLens[i];
        }
        if (sumLen != expectedLen) {
            throw new IllegalStateException(
                    "Novelty segment length sum mismatch: " + sumLen + " vs " + expectedLen);
        }
        int[] compLens = new int[totalSegs];
        for (int i = 0; i < totalSegs; i++) {
            compLens[i] = in.readVarInt();
        }
        byte[][] payloads = new byte[totalSegs][];
        for (int i = 0; i < totalSegs; i++) {
            payloads[i] = in.readBytes(compLens[i]);
        }

        byte[] out = new byte[expectedLen];
        if (encodedTurbo) {
            decodeTurbo(payloads, out, segLens, threads, boundedMem);
            return out;
        }
        if (totalSegs == 1) {
            FastCrmModel model =
                    new FastCrmModel(null, encodedBounded, encodedPanSuper, 0, expectedLen, 1);
            decodeSegmentWithSafety(payloads[0], out, 0, segLens[0], model, false);
            return out;
        }

        int warmBases = segLens[0];
        FastCrmModel warmModel =
                new FastCrmModel(
                        null, encodedBounded, encodedPanSuper, warmBases, expectedLen, totalSegs);
        decodeSegmentWithSafety(payloads[0], out, 0, warmBases, warmModel, true);
        WarmSnapshot snapshot = warmModel.snapshot();
        warmModel = null; // Free warm-up model before spawning parallel segment workers!

        int[] segStarts = new int[totalSegs];
        int cur = 0;
        for (int i = 0; i < totalSegs; i++) {
            segStarts[i] = cur;
            cur += segLens[i];
        }

        int maxSafeWorkers =
                boundedMem
                        ? Math.min(Math.max(1, threads), 6)
                        : (encodedPanSuper ? Math.min(Math.max(1, threads), 8) : Math.max(1, threads));
        int workers = Math.min(maxSafeWorkers, totalSegs - 1);
        if (workers <= 1) {
            for (int i = 1; i < totalSegs; i++) {
                FastCrmModel m =
                        new FastCrmModel(
                                snapshot,
                                encodedBounded,
                                encodedPanSuper,
                                warmBases,
                                expectedLen,
                                totalSegs);
                decodeSegmentWithSafety(payloads[i], out, segStarts[i], segLens[i], m, false);
            }
        } else {
            ForkJoinPool pool = new ForkJoinPool(workers);
            try {
                pool.submit(
                                () ->
                                        IntStream.range(1, totalSegs)
                                                .parallel()
                                                .forEach(
                                                        i -> {
                                                            FastCrmModel m =
                                                                    new FastCrmModel(
                                                                            snapshot,
                                                                            encodedBounded,
                                                                            encodedPanSuper,
                                                                            warmBases,
                                                                            expectedLen,
                                                                            totalSegs);
                                                            decodeSegmentWithSafety(
                                                                    payloads[i],
                                                                    out,
                                                                    segStarts[i],
                                                                    segLens[i],
                                                                    m,
                                                                    false);
                                                        }))
                        .join();
            } finally {
                pool.shutdown();
            }
        }
        return out;
    }

    private static byte[] encodeSegmentWithSafety(
            byte[] bases, int start, int len, FastCrmModel model) {
        int rawBytesLen = (len + 3) >>> 2;
        ByteOut crmOut = new ByteOut(rawBytesLen + 32);
        crmOut.writeByte(SEG_CRM);
        RangeEnc rc = new RangeEnc(crmOut);
        int[] freqs = new int[4];

        int end = start + len;
        for (int i = start; i < end; i++) {
            int relPos = i - start;
            model.predict(bases, start, relPos, freqs);
            int sym = bases[i] & 3;
            int cumLow = 0;
            for (int s = 0; s < sym; s++) {
                cumLow += freqs[s];
            }
            rc.encodeSymbol(cumLow, freqs[sym]);
            model.update(bases, start, relPos, sym);
        }
        rc.flush();

        if (crmOut.size() <= 1 + rawBytesLen) {
            return crmOut.toByteArray();
        }
        return packRaw2BitSegment(bases, start, len);
    }

    private static byte[] packRaw2BitSegment(byte[] bases, int start, int len) {
        int rawBytesLen = (len + 3) >>> 2;
        int end = start + len;
        byte[] raw = new byte[1 + rawBytesLen];
        raw[0] = SEG_RAW_2BIT;
        int full = len >>> 2;
        int idx = start;
        for (int b = 0; b < full; b++) {
            int v =
                    (bases[idx] & 3)
                            | ((bases[idx + 1] & 3) << 2)
                            | ((bases[idx + 2] & 3) << 4)
                            | ((bases[idx + 3] & 3) << 6);
            raw[1 + b] = (byte) v;
            idx += 4;
        }
        if (idx < end) {
            int v = 0;
            int shift = 0;
            while (idx < end) {
                v |= (bases[idx++] & 3) << shift;
                shift += 2;
            }
            raw[1 + full] = (byte) v;
        }
        return raw;
    }

    private static byte[] encodeTurbo(byte[] bases, int len, int threads, boolean boundedMem) {
        int totalSegs = (len + TURBO_SEG_BASES - 1) / TURBO_SEG_BASES;
        int[] segStarts = new int[totalSegs];
        int[] segLens = new int[totalSegs];
        for (int i = 0; i < totalSegs; i++) {
            segStarts[i] = i * TURBO_SEG_BASES;
            segLens[i] = Math.min(TURBO_SEG_BASES, len - segStarts[i]);
        }
        byte[][] encoded = new byte[totalSegs][];
        int workers = Math.max(1, Math.min(threads, Math.min(totalSegs, boundedMem ? 4 : 16)));
        if (workers <= 1) {
            for (int i = 0; i < totalSegs; i++) {
                encoded[i] = encodeTurboSegment(bases, segStarts[i], segLens[i]);
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
                                                                encoded[i] =
                                                                        encodeTurboSegment(
                                                                                bases,
                                                                                segStarts[i],
                                                                                segLens[i])))
                        .join();
            } finally {
                pool.shutdown();
            }
        }
        int totalComp = 0;
        for (byte[] e : encoded) {
            totalComp += e.length;
        }
        int modelFlags = (boundedMem ? 1 : 0) | FLAG_TURBO_NOVELTY;
        ByteOut out = new ByteOut(totalComp + totalSegs * 8 + 16);
        out.writeVarInt(totalSegs);
        out.writeByte(modelFlags);
        for (int i = 0; i < totalSegs; i++) {
            out.writeVarInt(segLens[i]);
        }
        for (int i = 0; i < totalSegs; i++) {
            out.writeVarInt(encoded[i].length);
        }
        for (int i = 0; i < totalSegs; i++) {
            out.writeBytes(encoded[i]);
        }
        return out.toByteArray();
    }

    /**
     * Bounded-memory encoder over chunked per-unit novelty slices, progressively releasing
     * {@code slices[i] = null} so large chromosomes and collections avoid allocating a second
     * 1-byte/base superstream under tight {@code -Xmx256m} heaps.
     */
    public static byte[] encodeFromSlices(
            byte[][] slices, int[] sliceLens, int totalLen, int threads, boolean boundedMem) {
        if (totalLen <= 0) {
            ByteOut empty = new ByteOut(4);
            empty.writeVarInt(0);
            return empty.toByteArray();
        }
        int turboCutoff = boundedMem ? TURBO_THRESHOLD : 64_000_000;
        if (totalLen < turboCutoff) {
            byte[] flat = new byte[totalLen];
            int dst = 0;
            for (int i = 0; i < slices.length; i++) {
                int sl = sliceLens[i];
                if (sl > 0 && slices[i] != null) {
                    System.arraycopy(slices[i], 0, flat, dst, sl);
                    dst += sl;
                }
                slices[i] = null;
            }
            return encode(flat, totalLen, threads, boundedMem);
        }

        int totalSegs = (totalLen + TURBO_SEG_BASES - 1) / TURBO_SEG_BASES;
        int[] segLens = new int[totalSegs];
        for (int i = 0; i < totalSegs; i++) {
            segLens[i] = Math.min(TURBO_SEG_BASES, totalLen - i * TURBO_SEG_BASES);
        }

        byte[][] encoded = new byte[totalSegs][];
        int sliceIdx = 0;
        int sliceOff = 0;
        int batchWorkers = boundedMem ? 1 : Math.max(1, Math.min(threads, 4));

        if (batchWorkers <= 1) {
            byte[] scratch = new byte[Math.min(TURBO_SEG_BASES, totalLen)];
            for (int s = 0; s < totalSegs; s++) {
                int need = segLens[s];
                int filled = 0;
                while (filled < need && sliceIdx < slices.length) {
                    int sl = sliceLens[sliceIdx];
                    int avail = sl - sliceOff;
                    if (avail <= 0) {
                        slices[sliceIdx] = null;
                        sliceIdx++;
                        sliceOff = 0;
                        continue;
                    }
                    int take = Math.min(avail, need - filled);
                    System.arraycopy(slices[sliceIdx], sliceOff, scratch, filled, take);
                    filled += take;
                    sliceOff += take;
                    if (sliceOff >= sl) {
                        slices[sliceIdx] = null;
                        sliceIdx++;
                        sliceOff = 0;
                    }
                }
                encoded[s] = encodeTurboSegment(scratch, 0, need);
            }
        } else {
            byte[][] batchBufs = new byte[batchWorkers][];
            int[] batchLens = new int[batchWorkers];
            for (int w = 0; w < batchWorkers; w++) {
                batchBufs[w] = new byte[Math.min(TURBO_SEG_BASES, totalLen)];
            }
            ForkJoinPool pool = new ForkJoinPool(batchWorkers);
            try {
                for (int baseSeg = 0; baseSeg < totalSegs; baseSeg += batchWorkers) {
                    final int count = Math.min(batchWorkers, totalSegs - baseSeg);
                    for (int w = 0; w < count; w++) {
                        int need = segLens[baseSeg + w];
                        batchLens[w] = need;
                        int filled = 0;
                        while (filled < need && sliceIdx < slices.length) {
                            int sl = sliceLens[sliceIdx];
                            int avail = sl - sliceOff;
                            if (avail <= 0) {
                                slices[sliceIdx] = null;
                                sliceIdx++;
                                sliceOff = 0;
                                continue;
                            }
                            int take = Math.min(avail, need - filled);
                            System.arraycopy(slices[sliceIdx], sliceOff, batchBufs[w], filled, take);
                            filled += take;
                            sliceOff += take;
                            if (sliceOff >= sl) {
                                slices[sliceIdx] = null;
                                sliceIdx++;
                                sliceOff = 0;
                            }
                        }
                    }
                    final int startSeg = baseSeg;
                    pool.submit(
                                    () ->
                                            IntStream.range(0, count)
                                                    .parallel()
                                                    .forEach(
                                                            w ->
                                                                    encoded[startSeg + w] =
                                                                            encodeTurboSegment(
                                                                                    batchBufs[w],
                                                                                    0,
                                                                                    batchLens[w])))
                            .join();
                }
            } finally {
                pool.shutdown();
            }
        }

        long totalComp = 0L;
        for (byte[] e : encoded) {
            totalComp += e.length;
        }
        int modelFlags = (boundedMem ? 1 : 0) | FLAG_TURBO_NOVELTY;
        ByteOut out = new ByteOut((int) Math.min(Integer.MAX_VALUE - 64L, totalComp + totalSegs * 8L + 16L));
        out.writeVarInt(totalSegs);
        out.writeByte(modelFlags);
        for (int i = 0; i < totalSegs; i++) {
            out.writeVarInt(segLens[i]);
        }
        for (int i = 0; i < totalSegs; i++) {
            out.writeVarInt(encoded[i].length);
        }
        for (int i = 0; i < totalSegs; i++) {
            out.writeBytes(encoded[i]);
            encoded[i] = null;
        }
        return out.toByteArray();
    }

    /** Streaming / lazy reader over novelty payloads keeping <= 4 MB live during Turbo bounded decode. */
    public static final class NoveltyReader {
        private final boolean lazyTurbo;
        private final int[] segLens;
        private final byte[][] payloads;
        private int segIdx = -1;
        private byte[] curSegBuf;
        private int curSegLen = 0;
        private int curSegOff = 0;

        private final byte[] fullBuf;
        private int fullOff = 0;

        private NoveltyReader(int[] segLens, byte[][] payloads) {
            this.lazyTurbo = true;
            this.segLens = segLens;
            this.payloads = payloads;
            this.fullBuf = null;
        }

        private NoveltyReader(byte[] fullBuf) {
            this.lazyTurbo = false;
            this.segLens = null;
            this.payloads = null;
            this.fullBuf = fullBuf;
        }

        public byte[] nextSlice(int len) {
            if (len <= 0) {
                return new byte[0];
            }
            if (!lazyTurbo) {
                byte[] out = Arrays.copyOfRange(fullBuf, fullOff, fullOff + len);
                fullOff += len;
                return out;
            }
            byte[] out = new byte[len];
            int filled = 0;
            while (filled < len) {
                if (segIdx < 0 || curSegOff >= curSegLen) {
                    segIdx++;
                    if (segIdx >= payloads.length) {
                        throw new IllegalStateException("Novelty stream exhausted early");
                    }
                    int sl = segLens[segIdx];
                    if (curSegBuf == null || curSegBuf.length < sl) {
                        curSegBuf = new byte[sl];
                    }
                    decodeTurboSegment(payloads[segIdx], curSegBuf, 0, sl);
                    payloads[segIdx] = null; // release compressed payload incrementally!
                    curSegLen = sl;
                    curSegOff = 0;
                }
                int avail = curSegLen - curSegOff;
                int take = Math.min(avail, len - filled);
                System.arraycopy(curSegBuf, curSegOff, out, filled, take);
                filled += take;
                curSegOff += take;
            }
            return out;
        }
    }

    public static NoveltyReader openNoveltyReader(
            ByteIn in, int expectedLen, int threads, boolean boundedMem) {
        if (expectedLen == 0) {
            int totalSegs = in.readVarInt();
            if (totalSegs != 0) {
                throw new IllegalStateException("Expected 0 novelty segments for empty stream");
            }
            return new NoveltyReader(new byte[0]);
        }
        int mark = in.position();
        int totalSegs = in.readVarInt();
        if (totalSegs <= 0 || totalSegs > expectedLen / 1024 + 2) {
            throw new IllegalStateException("Malformed novelty segment count: " + totalSegs);
        }
        int modelFlags = in.readByte() & 0xFF;
        if ((modelFlags & ~7) != 0) {
            throw new IllegalStateException("Unsupported novelty model flags: " + modelFlags);
        }
        boolean encodedTurbo = (modelFlags & FLAG_TURBO_NOVELTY) != 0;
        if (encodedTurbo && (boundedMem || expectedLen >= 32_000_000)) {
            int[] segLens = new int[totalSegs];
            int sumLen = 0;
            for (int i = 0; i < totalSegs; i++) {
                segLens[i] = in.readVarInt();
                if (segLens[i] < 0 || sumLen + segLens[i] > expectedLen) {
                    throw new IllegalStateException("Invalid novelty segment length");
                }
                sumLen += segLens[i];
            }
            if (sumLen != expectedLen) {
                throw new IllegalStateException(
                        "Novelty segment length sum mismatch: " + sumLen + " vs " + expectedLen);
            }
            int[] compLens = new int[totalSegs];
            for (int i = 0; i < totalSegs; i++) {
                compLens[i] = in.readVarInt();
            }
            byte[][] payloads = new byte[totalSegs][];
            for (int i = 0; i < totalSegs; i++) {
                payloads[i] = in.readBytes(compLens[i]);
            }
            return new NoveltyReader(segLens, payloads);
        }
        in.setPosition(mark);
        return new NoveltyReader(decode(in, expectedLen, threads, boundedMem));
    }

    private static void decodeTurbo(
            byte[][] payloads, byte[] out, int[] segLens, int threads, boolean boundedMem) {
        int totalSegs = payloads.length;
        int[] segStarts = new int[totalSegs];
        int cur = 0;
        for (int i = 0; i < totalSegs; i++) {
            segStarts[i] = cur;
            cur += segLens[i];
        }
        int workers = Math.max(1, Math.min(threads, Math.min(totalSegs, boundedMem ? 4 : 16)));
        if (workers <= 1) {
            for (int i = 0; i < totalSegs; i++) {
                decodeTurboSegment(payloads[i], out, segStarts[i], segLens[i]);
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
                                                                decodeTurboSegment(
                                                                        payloads[i],
                                                                        out,
                                                                        segStarts[i],
                                                                        segLens[i])))
                        .join();
            } finally {
                pool.shutdown();
            }
        }
    }

    /**
     * L2-cache-resident Turbo SWAR Order-4/6/8 + Fuzzy Repeat context model (~401 KiB/worker),
     * fusing prediction, 4-symbol range coding, and count adaptation in one tight loop.
     */
    private static byte[] encodeTurboSegment(byte[] bases, int start, int len) {
        int rawBytesLen = (len + 3) >>> 2;
        ByteOut crmOut = new ByteOut(rawBytesLen + 64);
        crmOut.writeByte(SEG_CRM);
        RangeEnc rc = new RangeEnc(crmOut);

        int[] tab4 = new int[256];
        int[] tab6 = new int[4096];
        int[] tab8 = new int[65536];
        int[] repTab = new int[32768];
        int ctx4 = 0, ctx6 = 0, ctx8 = 0, ctx12 = 0;
        int repPos = -1, repStreak = 0, repMisses = 0;

        for (int r = 0; r < len; r++) {
            int c4 = tab4[ctx4];
            int c6 = tab6[ctx6];
            int c8 = tab8[ctx8];
            int w0 = 2 + (c4 & 0x7F) + ((c6 & 0x7F) << 1) + ((c8 & 0x7F) << 2);
            int w1 = 2 + ((c4 >>> 8) & 0x7F) + (((c6 >>> 8) & 0x7F) << 1) + (((c8 >>> 8) & 0x7F) << 2);
            int w2 = 2 + ((c4 >>> 16) & 0x7F) + (((c6 >>> 16) & 0x7F) << 1) + (((c8 >>> 16) & 0x7F) << 2);
            int w3 = 2 + ((c4 >>> 24) & 0x7F) + (((c6 >>> 24) & 0x7F) << 1) + (((c8 >>> 24) & 0x7F) << 2);
            int repPred = -1;
            if (repPos >= 0) {
                repPred = bases[start + repPos] & 3;
                int boost = Math.min(2048, repStreak << 6);
                switch (repPred) {
                    case 0 -> w0 += boost;
                    case 1 -> w1 += boost;
                    case 2 -> w2 += boost;
                    default -> w3 += boost;
                }
            }

            int tot = w0 + w1 + w2 + w3;
            long inv = RECIP[tot];
            int f0 = (int) ((w0 * inv) >>> 15);
            int f1 = (int) ((w1 * inv) >>> 15);
            int f2 = (int) ((w2 * inv) >>> 15);
            int f3 = TOT - f0 - f1 - f2;

            int sym = bases[start + r] & 3;
            int cumLow;
            int freq;
            switch (sym) {
                case 0 -> { cumLow = 0; freq = f0; }
                case 1 -> { cumLow = f0; freq = f1; }
                case 2 -> { cumLow = f0 + f1; freq = f2; }
                default -> { cumLow = f0 + f1 + f2; freq = f3; }
            }
            rc.encodeSymbol(cumLow, freq);

            int shift = sym << 3;
            int nv4 = c4 + (1 << shift);
            tab4[ctx4] = ((nv4 & 0x80808080) != 0) ? ((nv4 & 0xFEFEFEFE) >>> 1) : nv4;
            int nv6 = c6 + (1 << shift);
            tab6[ctx6] = ((nv6 & 0x80808080) != 0) ? ((nv6 & 0xFEFEFEFE) >>> 1) : nv6;
            int nv8 = c8 + (1 << shift);
            tab8[ctx8] = ((nv8 & 0x80808080) != 0) ? ((nv8 & 0xFEFEFEFE) >>> 1) : nv8;

            if (repPos >= 0) {
                if (sym == repPred) {
                    repStreak++;
                    repMisses = 0;
                    repPos++;
                } else {
                    repMisses++;
                    if (repMisses > 2) {
                        repPos = -1;
                        repStreak = 0;
                    } else {
                        repPos++;
                    }
                }
            }
            if (r >= 12) {
                int h12 = ((ctx12 * 0x9E3779B1) >>> 17) & 0x7FFF;
                if (repPos < 0) {
                    int cand = repTab[h12] - 1;
                    if (cand >= 12
                            && (r - cand) > 16
                            && (bases[start + cand] & 3) == sym
                            && bases[start + cand - 1] == bases[start + r - 1]
                            && bases[start + cand - 2] == bases[start + r - 2]
                            && bases[start + cand - 3] == bases[start + r - 3]) {
                        repPos = cand + 1;
                        repStreak = 4;
                        repMisses = 0;
                    }
                }
                repTab[h12] = r + 1;
            }

            ctx4 = ((ctx4 << 2) | sym) & 0xFF;
            ctx6 = ((ctx6 << 2) | sym) & 0xFFF;
            ctx8 = ((ctx8 << 2) | sym) & 0xFFFF;
            ctx12 = ((ctx12 << 2) | sym) & 0xFFFFFF;
        }
        rc.flush();

        if (crmOut.size() <= 1 + rawBytesLen) {
            return crmOut.toByteArray();
        }
        return packRaw2BitSegment(bases, start, len);
    }

    private static void decodeTurboSegment(byte[] payload, byte[] dst, int dstStart, int len) {
        if (payload.length == 0) {
            throw new IllegalStateException("Empty turbo novelty segment payload");
        }
        byte mode = payload[0];
        if (mode == SEG_RAW_2BIT) {
            int expectedBytes = 1 + ((len + 3) >>> 2);
            if (payload.length != expectedBytes) {
                throw new IllegalStateException("Malformed raw 2-bit novelty payload length");
            }
            int idx = dstStart;
            int end = dstStart + len;
            for (int b = 1; b < payload.length && idx < end; b++) {
                int v = payload[b] & 0xFF;
                for (int j = 0; j < 4 && idx < end; j++) {
                    dst[idx++] = (byte) (v & 3);
                    v >>>= 2;
                }
            }
            return;
        }
        if (mode != SEG_CRM) {
            throw new IllegalStateException("Unknown turbo novelty segment mode: " + mode);
        }
        ByteIn in = new ByteIn(payload, 1, payload.length - 1);
        RangeDec rd = new RangeDec(in);

        int[] tab4 = new int[256];
        int[] tab6 = new int[4096];
        int[] tab8 = new int[65536];
        int[] repTab = new int[32768];
        int ctx4 = 0, ctx6 = 0, ctx8 = 0, ctx12 = 0;
        int repPos = -1, repStreak = 0, repMisses = 0;

        for (int r = 0; r < len; r++) {
            int c4 = tab4[ctx4];
            int c6 = tab6[ctx6];
            int c8 = tab8[ctx8];
            int w0 = 2 + (c4 & 0x7F) + ((c6 & 0x7F) << 1) + ((c8 & 0x7F) << 2);
            int w1 = 2 + ((c4 >>> 8) & 0x7F) + (((c6 >>> 8) & 0x7F) << 1) + (((c8 >>> 8) & 0x7F) << 2);
            int w2 = 2 + ((c4 >>> 16) & 0x7F) + (((c6 >>> 16) & 0x7F) << 1) + (((c8 >>> 16) & 0x7F) << 2);
            int w3 = 2 + ((c4 >>> 24) & 0x7F) + (((c6 >>> 24) & 0x7F) << 1) + (((c8 >>> 24) & 0x7F) << 2);
            int repPred = -1;
            if (repPos >= 0) {
                repPred = dst[dstStart + repPos] & 3;
                int boost = Math.min(2048, repStreak << 6);
                switch (repPred) {
                    case 0 -> w0 += boost;
                    case 1 -> w1 += boost;
                    case 2 -> w2 += boost;
                    default -> w3 += boost;
                }
            }

            int tot = w0 + w1 + w2 + w3;
            long inv = RECIP[tot];
            int f0 = (int) ((w0 * inv) >>> 15);
            int f1 = (int) ((w1 * inv) >>> 15);
            int f2 = (int) ((w2 * inv) >>> 15);
            int f3 = TOT - f0 - f1 - f2;

            int target = rd.getThreshold();
            int cum;
            int freq;
            int sym;
            if (target < f0) {
                sym = 0;
                cum = 0;
                freq = f0;
            } else if (target < f0 + f1) {
                sym = 1;
                cum = f0;
                freq = f1;
            } else if (target < f0 + f1 + f2) {
                sym = 2;
                cum = f0 + f1;
                freq = f2;
            } else {
                sym = 3;
                cum = f0 + f1 + f2;
                freq = f3;
            }
            rd.decodeSymbol(cum, freq);
            dst[dstStart + r] = (byte) sym;

            int shift = sym << 3;
            int nv4 = c4 + (1 << shift);
            tab4[ctx4] = ((nv4 & 0x80808080) != 0) ? ((nv4 & 0xFEFEFEFE) >>> 1) : nv4;
            int nv6 = c6 + (1 << shift);
            tab6[ctx6] = ((nv6 & 0x80808080) != 0) ? ((nv6 & 0xFEFEFEFE) >>> 1) : nv6;
            int nv8 = c8 + (1 << shift);
            tab8[ctx8] = ((nv8 & 0x80808080) != 0) ? ((nv8 & 0xFEFEFEFE) >>> 1) : nv8;

            if (repPos >= 0) {
                if (sym == repPred) {
                    repStreak++;
                    repMisses = 0;
                    repPos++;
                } else {
                    repMisses++;
                    if (repMisses > 2) {
                        repPos = -1;
                        repStreak = 0;
                    } else {
                        repPos++;
                    }
                }
            }
            if (r >= 12) {
                int h12 = ((ctx12 * 0x9E3779B1) >>> 17) & 0x7FFF;
                if (repPos < 0) {
                    int cand = repTab[h12] - 1;
                    if (cand >= 12
                            && (r - cand) > 16
                            && (dst[dstStart + cand] & 3) == sym
                            && dst[dstStart + cand - 1] == dst[dstStart + r - 1]
                            && dst[dstStart + cand - 2] == dst[dstStart + r - 2]
                            && dst[dstStart + cand - 3] == dst[dstStart + r - 3]) {
                        repPos = cand + 1;
                        repStreak = 4;
                        repMisses = 0;
                    }
                }
                repTab[h12] = r + 1;
            }

            ctx4 = ((ctx4 << 2) | sym) & 0xFF;
            ctx6 = ((ctx6 << 2) | sym) & 0xFFF;
            ctx8 = ((ctx8 << 2) | sym) & 0xFFFF;
            ctx12 = ((ctx12 << 2) | sym) & 0xFFFFFF;
        }
    }

    private static void decodeSegmentWithSafety(
            byte[] payload,
            byte[] dst,
            int dstStart,
            int len,
            FastCrmModel model,
            boolean buildWarmState) {
        if (payload.length == 0) {
            throw new IllegalStateException("Empty novelty segment payload");
        }
        byte mode = payload[0];
        if (mode == SEG_RAW_2BIT) {
            int expectedBytes = 1 + ((len + 3) >>> 2);
            if (payload.length != expectedBytes) {
                throw new IllegalStateException("Malformed raw 2-bit novelty payload length");
            }
            int idx = dstStart;
            int end = dstStart + len;
            for (int b = 1; b < payload.length && idx < end; b++) {
                int v = payload[b] & 0xFF;
                for (int j = 0; j < 4 && idx < end; j++) {
                    dst[idx++] = (byte) (v & 3);
                    v >>>= 2;
                }
            }
            if (buildWarmState) {
                int[] dummy = new int[4];
                for (int r = 0; r < len; r++) {
                    model.predict(dst, 0, r, dummy);
                    model.update(dst, 0, r, dst[r] & 3);
                }
            }
            return;
        }
        if (mode != SEG_CRM) {
            throw new IllegalStateException("Unknown novelty segment mode: " + mode);
        }
        ByteIn in = new ByteIn(payload, 1, payload.length - 1);
        RangeDec rd = new RangeDec(in);
        int[] freqs = new int[4];
        for (int r = 0; r < len; r++) {
            model.predict(dst, dstStart, r, freqs);
            int target = rd.getThreshold();
            int cum = 0;
            int sym = 0;
            while (sym < 3 && cum + freqs[sym] <= target) {
                cum += freqs[sym];
                sym++;
            }
            rd.decodeSymbol(cum, freqs[sym]);
            dst[dstStart + r] = (byte) sym;
            model.update(dst, dstStart, r, sym);
        }
    }

    private record WarmSnapshot(
            short[] ord0,
            short[] ord1,
            short[] ord2,
            short[] ord3,
            short[] ord4,
            short[] ord5,
            short[] ord6,
            short[] ord7,
            int[] ord8,
            int[] ord9,
            short[] orfPlus2,
            short[] orfMinus2,
            short[] orfPlus3,
            short[] orfMinus3,
            short[] orfPlus5,
            short[] orfMinus5,
            short[] wobble8,
            short[] gapped8,
            long[] sparse12,
            long[] repFineFwd,
            long[] repFineRev,
            long[] repLongFwd,
            long[] repLongRev,
            short[] sseCounts,
            int[][] weights,
            long[] bayesScores,
            long metaBayes,
            long metaNeural) {}

    /**
     * L2/L3-cache-resident Chargaff-Symmetric Context + 6-State Causal Open-Reading-Frame HMM +
     * Spaced-Seed + 6-Head Bidirectional Kimura-XOR Repeat Model (19 experts in single/bounded
     * profile; 22 experts in pan-superstream profile).
     */
    private static final class FastCrmModel {
        private final boolean panSuper;
        private final boolean useOrd8;
        private final int warmPrefixLimit;
        private final int numExperts;

        // Dense 16-bit Markov tables (ord0..ord7: 174.7 KiB total, L2-resident!) + SWAR ord8/ord9
        private final short[] ord0;
        private final short[] ord1;
        private final short[] ord2;
        private final short[] ord3;
        private final short[] ord4;
        private final short[] ord5;
        private final short[] ord6;
        private final short[] ord7;
        private final int[] ord8; // non-null only when panSuper
        private final int[] ord9; // non-null only when panSuper

        // 6-State Causal Open-Reading-Frame HMM tables: cache-line grouped [((ctx * 3) + c) << 2]
        // Tracks {+0, +1, +2, -0, -1, -2} reading frames dynamically across gene boundaries!
        private final short[] orfPlus2;
        private final short[] orfMinus2;
        private final short[] orfPlus3;
        private final short[] orfMinus3;
        private final short[] orfPlus5;
        private final short[] orfMinus5;
        private final long[] frameLog = new long[6];
        private final int[] frameW = new int[6];
        private final int[][] pPlus2 = new int[3][4];
        private final int[][] pMinus2 = new int[3][4];
        private final int[][] pPlus3 = new int[3][4];
        private final int[][] pMinus3 = new int[3][4];
        private final int[][] pPlus5 = new int[3][4];
        private final int[][] pMinus5 = new int[3][4];
        private int bestFrame = 0;
        private int frameWsum = 6;

        // Spaced-seed tables
        private final short[] wobble8;
        private final short[] gapped8;

        // Pan-superstream Sparse-12 table (non-null only when panSuper)
        private final long[] sparse12;
        private final int sparseMask;

        // O(1) 32-bit-tagged Repeat Seed Rings (zero false positives, zero byte-verification loops!)
        // Entry packing in long: [tag32 (bits 63..32) | pos32 (bits 31..0)] (tag == 0 means empty)
        private static final int K_FINE = 13;
        private static final long MASK_FINE = (1L << (K_FINE << 1)) - 1L;
        private static final int K_LONG = 18;
        private static final long MASK_LONG = (1L << (K_LONG << 1)) - 1L;

        private final long[] repFineFwd; // 2-way associative (idx, idx^1)
        private final long[] repFineRev; // 1-way
        private final int repFineMask;

        private final long[] repLongFwd; // 2-way associative (idx, idx^1)
        private final long[] repLongRev; // 1-way
        private final int repLongMask;

        private final MicroRepeatHead[] heads = new MicroRepeatHead[6];
        private final int idxFirstHead;

        private final int[][] weights;
        private final long[] bayesScores;

        private final short[] sseCounts;
        private int lastSseState = -1;

        private final int[][] expP;
        private final int[][] expLnP;
        private final int[] uniformP = {4096, 4096, 4096, 4096};
        private final int[] p0Buf = {4096, 4096, 4096, 4096};
        private final int[] p2Buf = {4096, 4096, 4096, 4096};
        private final int[] neuralP = new int[4];
        private final int[] bayesP = new int[4];
        private long metaBayes;
        private long metaNeural;

        private long fwdHist = 0L;
        private long rcHist = 0L;
        private int activeRegime = 0;

        private int cachedSparseIdx;
        private int cachedSparseTag;

        FastCrmModel(
                WarmSnapshot snap,
                boolean boundedMem,
                boolean panSuper,
                int warmPrefixLimit,
                int totalStreamLen,
                int totalSegs) {
            this.panSuper = panSuper;
            this.useOrd8 = panSuper;
            this.warmPrefixLimit = warmPrefixLimit;
            // Base experts:
            // 0..7: ord0..ord7 (8)
            // 8,9,10: codon1, codon3, codon5 (3)
            // 11,12: wobble8, gapped8 (2)
            // [if useOrd8: 13=ord8 (+1)]
            // [if panSuper: 14=ord9, 15=sparse12 (+2)]
            // last 6: heads[0..5] (6)
            this.idxFirstHead = panSuper ? 16 : (useOrd8 ? 14 : 13);
            this.numExperts = idxFirstHead + 6;

            for (int h = 0; h < 6; h++) {
                this.heads[h] = new MicroRepeatHead();
            }

            int spSize = panSuper ? (1 << 18) : 0;
            this.sparseMask = spSize - 1;

            int ringBuckets =
                    (panSuper && totalSegs >= 4)
                            ? (1 << 16)
                            : (boundedMem
                                    ? (1 << 14)
                                    : (totalStreamLen >= 500_000 ? (1 << 19) : (1 << 16)));
            this.repFineMask = ringBuckets - 1;
            this.repLongMask = ringBuckets - 1;

            this.expP = new int[numExperts][4];
            this.expLnP = new int[numExperts][4];

            if (snap == null) {
                this.ord0 = new short[4];
                this.ord1 = new short[4 * 4];
                this.ord2 = new short[16 * 4];
                this.ord3 = new short[64 * 4];
                this.ord4 = new short[256 * 4];
                this.ord5 = new short[1024 * 4];
                this.ord6 = new short[4096 * 4];
                this.ord7 = new short[16384 * 4];
                this.ord8 = useOrd8 ? new int[65536] : null;
                this.ord9 = panSuper ? new int[262144] : null;

                this.orfPlus2 = new short[16 * 3 * 4];
                this.orfMinus2 = new short[16 * 3 * 4];
                this.orfPlus3 = new short[64 * 3 * 4];
                this.orfMinus3 = new short[64 * 3 * 4];
                this.orfPlus5 = new short[1024 * 3 * 4];
                this.orfMinus5 = new short[1024 * 3 * 4];
                seedUniversalStopPriors(this.orfPlus2, this.orfMinus2, this.orfPlus3, this.orfMinus3);
                this.wobble8 = new short[65536 * 4];
                this.gapped8 = new short[65536 * 4];
                this.sparse12 = panSuper ? new long[spSize] : null;

                this.repFineFwd = new long[ringBuckets << 1];
                this.repFineRev = new long[ringBuckets];
                this.repLongFwd = new long[ringBuckets << 1];
                this.repLongRev = new long[ringBuckets];
                this.sseCounts = new short[1536 * 4];

                this.bayesScores = new long[numExperts];
                for (int k = 0; k <= 7; k++) {
                    this.bayesScores[k] = (k <= 3) ? 0L : -6144L;
                }
                this.bayesScores[8] = -2048L;
                this.bayesScores[9] = -3276L;
                this.bayesScores[10] = -4500L;
                this.bayesScores[11] = -4096L;
                this.bayesScores[12] = -4096L;
                if (useOrd8) {
                    this.bayesScores[13] = -6144L;
                }
                if (panSuper) {
                    this.bayesScores[14] = -6144L;
                    this.bayesScores[15] = -6144L;
                }
                for (int h = 0; h < 6; h++) {
                    this.bayesScores[idxFirstHead + h] = -4096L;
                }

                this.weights = new int[8][numExperts];
                for (int r = 0; r < 8; r++) {
                    boolean repReg = (r & 4) != 0;
                    boolean hiReg = (r & 1) != 0;
                    this.weights[r][0] = 150;
                    for (int k = 1; k <= 4; k++) {
                        this.weights[r][k] = hiReg ? 350 : 750;
                    }
                    for (int k = 5; k <= 7; k++) {
                        this.weights[r][k] = hiReg ? 650 : 250;
                    }
                    this.weights[r][8] = 220;
                    this.weights[r][9] = 320;
                    this.weights[r][10] = 360;
                    this.weights[r][11] = 300;
                    this.weights[r][12] = 300;
                    if (useOrd8) {
                        this.weights[r][13] = hiReg ? 650 : 220;
                    }
                    if (panSuper) {
                        this.weights[r][14] = hiReg ? 780 : 150;
                        this.weights[r][15] = hiReg ? 650 : 150;
                    }
                    for (int h = 0; h < 6; h++) {
                        this.weights[r][idxFirstHead + h] =
                                repReg ? (h == 0 || h == 2 ? 1200 : 600) : 100;
                    }
                }
                this.metaBayes = 0L;
                this.metaNeural = -2048L;
            } else {
                this.ord0 = snap.ord0.clone();
                this.ord1 = snap.ord1.clone();
                this.ord2 = snap.ord2.clone();
                this.ord3 = snap.ord3.clone();
                this.ord4 = snap.ord4.clone();
                this.ord5 = snap.ord5.clone();
                this.ord6 = snap.ord6.clone();
                this.ord7 = snap.ord7.clone();
                this.ord8 = snap.ord8 == null ? null : snap.ord8.clone();
                this.ord9 = snap.ord9 == null ? null : snap.ord9.clone();
                this.orfPlus2 = snap.orfPlus2.clone();
                this.orfMinus2 = snap.orfMinus2.clone();
                this.orfPlus3 = snap.orfPlus3.clone();
                this.orfMinus3 = snap.orfMinus3.clone();
                this.orfPlus5 = snap.orfPlus5.clone();
                this.orfMinus5 = snap.orfMinus5.clone();
                this.wobble8 = snap.wobble8.clone();
                this.gapped8 = snap.gapped8.clone();
                this.sparse12 = snap.sparse12 == null ? null : snap.sparse12.clone();
                this.repFineFwd = snap.repFineFwd.clone();
                this.repFineRev = snap.repFineRev.clone();
                this.repLongFwd = snap.repLongFwd.clone();
                this.repLongRev = snap.repLongRev.clone();
                this.sseCounts = snap.sseCounts.clone();
                this.weights = new int[8][];
                for (int r = 0; r < 8; r++) {
                    this.weights[r] = snap.weights[r].clone();
                }
                this.bayesScores = snap.bayesScores.clone();
                this.metaBayes = snap.metaBayes;
                this.metaNeural = snap.metaNeural;
            }
        }

        private static void seedUniversalStopPriors(
                short[] plus2, short[] minus2, short[] plus3, short[] minus3) {
            // Universal stop-codon phase symmetry breaker at c = 2:
            // + strand stops: TAA (ctx2=12, sym=0), TAG (ctx2=12, sym=2), TGA (ctx2=14, sym=0)
            // - strand RC stops on + strand: TTA (ctx2=15, sym=0), CTA (ctx2=7, sym=0), TCA (ctx2=13, sym=0)
            for (int ctx2 = 0; ctx2 < 16; ctx2++) {
                int bp2 = ((ctx2 * 3) + 2) << 2;
                if (ctx2 == 12) {
                    plus2[bp2 + 1] += 4;
                    plus2[bp2 + 3] += 4;
                } else if (ctx2 == 14) {
                    plus2[bp2 + 1] += 3;
                    plus2[bp2 + 2] += 3;
                    plus2[bp2 + 3] += 3;
                }
                int bm2 = ((ctx2 * 3) + 2) << 2;
                if (ctx2 == 15 || ctx2 == 7 || ctx2 == 13) {
                    minus2[bm2 + 1] += 3;
                    minus2[bm2 + 2] += 3;
                    minus2[bm2 + 3] += 3;
                }
                for (int b0 = 0; b0 < 4; b0++) {
                    int ctx3 = (b0 << 4) | ctx2;
                    int bp3 = ((ctx3 * 3) + 2) << 2;
                    if (ctx2 == 12) {
                        plus3[bp3 + 1] += 4;
                        plus3[bp3 + 3] += 4;
                    } else if (ctx2 == 14) {
                        plus3[bp3 + 1] += 3;
                        plus3[bp3 + 2] += 3;
                        plus3[bp3 + 3] += 3;
                    }
                    int bm3 = ((ctx3 * 3) + 2) << 2;
                    if (ctx2 == 15 || ctx2 == 7 || ctx2 == 13) {
                        minus3[bm3 + 1] += 3;
                        minus3[bm3 + 2] += 3;
                        minus3[bm3 + 3] += 3;
                    }
                }
            }
        }

        WarmSnapshot snapshot() {
            int[][] wCopy = new int[8][];
            for (int r = 0; r < 8; r++) {
                wCopy[r] = weights[r].clone();
            }
            return new WarmSnapshot(
                    ord0.clone(),
                    ord1.clone(),
                    ord2.clone(),
                    ord3.clone(),
                    ord4.clone(),
                    ord5.clone(),
                    ord6.clone(),
                    ord7.clone(),
                    ord8 == null ? null : ord8.clone(),
                    ord9 == null ? null : ord9.clone(),
                    orfPlus2.clone(),
                    orfMinus2.clone(),
                    orfPlus3.clone(),
                    orfMinus3.clone(),
                    orfPlus5.clone(),
                    orfMinus5.clone(),
                    wobble8.clone(),
                    gapped8.clone(),
                    sparse12 == null ? null : sparse12.clone(),
                    repFineFwd.clone(),
                    repFineRev.clone(),
                    repLongFwd.clone(),
                    repLongRev.clone(),
                    sseCounts.clone(),
                    wCopy,
                    bayesScores.clone(),
                    metaBayes,
                    metaNeural);
        }

        void predict(byte[] buf, int segStart, int relPos, int[] outFreqs) {
            int globalPos = segStart + relPos;
            spawnRepeatCandidates(segStart, relPos, globalPos);

            int highOrderCountSum = 0;

            // 0..7: L1/L2-cache-resident 16-bit Dense Markov Hierarchy with Hierarchical Pitman-Yor Backoff
            shortPredict(ord0, 0, uniformP, 32, expP[0]);

            if (relPos >= 1) {
                shortPredict(ord1, ((int) (fwdHist & 3L)) << 2, expP[0], 24, expP[1]);
            } else {
                System.arraycopy(expP[0], 0, expP[1], 0, 4);
            }

            if (relPos >= 2) {
                shortPredict(ord2, ((int) (fwdHist & 15L)) << 2, expP[1], 24, expP[2]);
            } else {
                System.arraycopy(expP[1], 0, expP[2], 0, 4);
            }

            if (relPos >= 3) {
                shortPredict(ord3, ((int) (fwdHist & 63L)) << 2, expP[2], 16, expP[3]);
            } else {
                System.arraycopy(expP[2], 0, expP[3], 0, 4);
            }

            if (relPos >= 4) {
                shortPredict(ord4, ((int) (fwdHist & 255L)) << 2, expP[3], 16, expP[4]);
            } else {
                System.arraycopy(expP[3], 0, expP[4], 0, 4);
            }

            if (relPos >= 5) {
                shortPredict(ord5, ((int) (fwdHist & 1023L)) << 2, expP[4], 10, expP[5]);
            } else {
                System.arraycopy(expP[4], 0, expP[5], 0, 4);
            }

            if (relPos >= 6) {
                shortPredict(ord6, ((int) (fwdHist & 4095L)) << 2, expP[5], 8, expP[6]);
            } else {
                System.arraycopy(expP[5], 0, expP[6], 0, 4);
            }

            if (relPos >= 7) {
                highOrderCountSum +=
                        shortPredict(ord7, ((int) (fwdHist & 16383L)) << 2, expP[6], 4, expP[7]);
            } else {
                System.arraycopy(expP[6], 0, expP[7], 0, 4);
            }

            // 8, 9, 10: 6-State Causal Open-Reading-Frame Hidden Markov Model (6F-ORF-HMM)
            long maxFL = frameLog[0];
            for (int s = 1; s < 6; s++) {
                if (frameLog[s] > maxFL) {
                    maxFL = frameLog[s];
                }
            }
            int fwSum = 0;
            int bestS = 0;
            for (int s = 0; s < 6; s++) {
                if (maxFL - frameLog[s] > 26624L) {
                    frameLog[s] = maxFL - 26624L; // ~6.5 nats gene-boundary switching floor
                }
                int w = Math.max(1, EXP_NEG_LUT[(int) (maxFL - frameLog[s])] >>> 3);
                frameW[s] = w;
                fwSum += w;
                if (w > frameW[bestS]) {
                    bestS = s;
                }
            }
            this.bestFrame = bestS;
            this.frameWsum = fwSum;

            int rmod3 = relPos % 3;
            int phase = rmod3;
            if (frameW[bestS] * 2 >= fwSum) {
                int cBest = rmod3 - (bestS < 3 ? bestS : bestS - 3);
                if (cBest < 0) cBest += 3;
                phase = (bestS < 3) ? cBest : (2 - cBest);
            }

            int ctx2Orf = (int) (fwdHist & 15L);
            int ctx3Orf = (int) (fwdHist & 63L);
            int ctx5Orf = (int) (fwdHist & 1023L);
            for (int c = 0; c < 3; c++) {
                if (relPos >= 2) {
                    shortPredict(orfPlus2, ((ctx2Orf * 3) + c) << 2, expP[1], 16, pPlus2[c]);
                    shortPredict(orfMinus2, ((ctx2Orf * 3) + c) << 2, expP[1], 16, pMinus2[c]);
                } else {
                    System.arraycopy(expP[1], 0, pPlus2[c], 0, 4);
                    System.arraycopy(expP[1], 0, pMinus2[c], 0, 4);
                }
                if (relPos >= 3) {
                    shortPredict(orfPlus3, ((ctx3Orf * 3) + c) << 2, pPlus2[c], 12, pPlus3[c]);
                    shortPredict(orfMinus3, ((ctx3Orf * 3) + c) << 2, pMinus2[c], 12, pMinus3[c]);
                } else {
                    System.arraycopy(pPlus2[c], 0, pPlus3[c], 0, 4);
                    System.arraycopy(pMinus2[c], 0, pMinus3[c], 0, 4);
                }
                if (relPos >= 5) {
                    shortPredict(orfPlus5, ((ctx5Orf * 3) + c) << 2, pPlus3[c], 8, pPlus5[c]);
                    shortPredict(orfMinus5, ((ctx5Orf * 3) + c) << 2, pMinus3[c], 8, pMinus5[c]);
                } else {
                    System.arraycopy(pPlus3[c], 0, pPlus5[c], 0, 4);
                    System.arraycopy(pMinus3[c], 0, pMinus5[c], 0, 4);
                }
            }

            long a20 = 0L, a21 = 0L, a22 = 0L, a23 = 0L;
            long a30 = 0L, a31 = 0L, a32 = 0L, a33 = 0L;
            long a50 = 0L, a51 = 0L, a52 = 0L, a53 = 0L;
            for (int s = 0; s < 6; s++) {
                long ws = frameW[s];
                int c = rmod3 - (s < 3 ? s : s - 3);
                if (c < 0) c += 3;
                int[] q2 = (s < 3) ? pPlus2[c] : pMinus2[c];
                int[] q3 = (s < 3) ? pPlus3[c] : pMinus3[c];
                int[] q5 = (s < 3) ? pPlus5[c] : pMinus5[c];
                a20 += ws * q2[0];
                a21 += ws * q2[1];
                a22 += ws * q2[2];
                a23 += ws * q2[3];
                a30 += ws * q3[0];
                a31 += ws * q3[1];
                a32 += ws * q3[2];
                a33 += ws * q3[3];
                a50 += ws * q5[0];
                a51 += ws * q5[1];
                a52 += ws * q5[2];
                a53 += ws * q5[3];
            }
            QuantumContextRepeatMixer.normalizeP(
                    divRoundSmall(a20, fwSum),
                    divRoundSmall(a21, fwSum),
                    divRoundSmall(a22, fwSum),
                    divRoundSmall(a23, fwSum),
                    expP[8]);
            QuantumContextRepeatMixer.normalizeP(
                    divRoundSmall(a30, fwSum),
                    divRoundSmall(a31, fwSum),
                    divRoundSmall(a32, fwSum),
                    divRoundSmall(a33, fwSum),
                    expP[9]);
            QuantumContextRepeatMixer.normalizeP(
                    divRoundSmall(a50, fwSum),
                    divRoundSmall(a51, fwSum),
                    divRoundSmall(a52, fwSum),
                    divRoundSmall(a53, fwSum),
                    expP[10]);

            // 11, 12: Wobble (110_110_110_11) & Palindromic Gapped-9 (1111_0_1111) Spaced Seeds
            if (relPos >= 11) {
                shortPredict(wobble8, extractWobbleKey(fwdHist) << 2, expP[4], 8, expP[11]);
            } else {
                System.arraycopy(expP[4], 0, expP[11], 0, 4);
            }
            if (relPos >= 9) {
                shortPredict(gapped8, extractGapped9Key(fwdHist) << 2, expP[4], 8, expP[12]);
            } else {
                System.arraycopy(expP[4], 0, expP[12], 0, 4);
            }

            // Optional Megabase / Pan-Superstream High-Order Experts (13, 14, 15)
            if (useOrd8) {
                if (relPos >= 8) {
                    int p8 = ord8[(int) (fwdHist & 65535L)];
                    highOrderCountSum += swarSum(p8);
                    swarPredict(p8, expP[7], 3, expP[13]);
                } else {
                    System.arraycopy(expP[7], 0, expP[13], 0, 4);
                }
            }
            if (panSuper) {
                if (relPos >= 9) {
                    int p9 = ord9[(int) (fwdHist & 262143L)];
                    highOrderCountSum += swarSum(p9);
                    swarPredict(p9, expP[13], 2, expP[14]);
                } else {
                    System.arraycopy(expP[13], 0, expP[14], 0, 4);
                }
                if (relPos >= 12) {
                    highOrderCountSum += predictSparse12(fwdHist & 0xFFFFFFL, expP[14], expP[15]);
                } else {
                    System.arraycopy(expP[14], 0, expP[15], 0, 4);
                }
            }

            // 6 Multi-Scale Bidirectional Kimura-XOR Repeat Heads
            boolean anyStrongRepeat = false;
            for (int h = 0; h < 6; h++) {
                MicroRepeatHead mh = heads[h];
                if (mh.predict(buf, segStart, globalPos, warmPrefixLimit, expP[idxFirstHead + h])) {
                    if (mh.streak >= 2 || mh.totalMatched >= 6) {
                        anyStrongRepeat = true;
                    }
                } else {
                    System.arraycopy(expP[4], 0, expP[idxFirstHead + h], 0, 4);
                }
            }

            activeRegime =
                    (anyStrongRepeat ? 4 : 0)
                            | ((phase == 2) ? 2 : 0)
                            | ((highOrderCountSum >= 4) ? 1 : 0);

            // Fused Log-Probability + Bayesian Opinion Pool + 8-Regime Multinomial Logistic Mix
            long maxBayes = bayesScores[0];
            for (int e = 1; e < numExperts; e++) {
                if (bayesScores[e] > maxBayes) {
                    maxBayes = bayesScores[e];
                }
            }

            long bAcc0 = 0L, bAcc1 = 0L, bAcc2 = 0L, bAcc3 = 0L;
            long bWsum = 0L;
            int[] wReg = weights[activeRegime];
            long logit0 = 0L, logit1 = 0L, logit2 = 0L, logit3 = 0L;

            for (int e = 0; e < numExperts; e++) {
                int[] pe = expP[e];
                int[] lne = expLnP[e];
                int lp0 = LN_LUT[pe[0]];
                int lp1 = LN_LUT[pe[1]];
                int lp2 = LN_LUT[pe[2]];
                int lp3 = LN_LUT[pe[3]];
                lne[0] = lp0;
                lne[1] = lp1;
                lne[2] = lp2;
                lne[3] = lp3;

                long diff = maxBayes - bayesScores[e];
                long bw = EXP_NEG_LUT[diff >= EXP_MAX_X ? EXP_MAX_X : (int) diff];
                bAcc0 += bw * pe[0];
                bAcc1 += bw * pe[1];
                bAcc2 += bw * pe[2];
                bAcc3 += bw * pe[3];
                bWsum += bw;

                long nw = wReg[e];
                logit0 += nw * lp0;
                logit1 += nw * lp1;
                logit2 += nw * lp2;
                logit3 += nw * lp3;
            }

            long bHalf = bWsum >> 1;
            long bInv = (1L << 36) / bWsum;
            int bp0 = divFast(bAcc0, bWsum, bHalf, bInv, 36);
            int bp1 = divFast(bAcc1, bWsum, bHalf, bInv, 36);
            int bp2 = divFast(bAcc2, bWsum, bHalf, bInv, 36);
            int bp3 = divFast(bAcc3, bWsum, bHalf, bInv, 36);
            QuantumContextRepeatMixer.normalizeP(bp0, bp1, bp2, bp3, bayesP);

            long l0 = logit0 >> 12, l1 = logit1 >> 12, l2 = logit2 >> 12, l3 = logit3 >> 12;
            long maxL = Math.max(Math.max(l0, l1), Math.max(l2, l3));
            int e0 = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, maxL - l0)];
            int e1 = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, maxL - l1)];
            int e2 = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, maxL - l2)];
            int e3 = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, maxL - l3)];
            int eSum = e0 + e1 + e2 + e3;
            long eHalf = eSum >> 1;
            long eInv = (1L << 32) / eSum;
            QuantumContextRepeatMixer.normalizeP(
                    divFast((long) e0 * P_SCALE, eSum, eHalf, eInv, 32),
                    divFast((long) e1 * P_SCALE, eSum, eHalf, eInv, 32),
                    divFast((long) e2 * P_SCALE, eSum, eHalf, eInv, 32),
                    divFast((long) e3 * P_SCALE, eSum, eHalf, eInv, 32),
                    neuralP);

            long maxMeta = Math.max(metaBayes, metaNeural);
            long mb = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, maxMeta - metaBayes)];
            long mn = EXP_NEG_LUT[(int) Math.min(EXP_MAX_X, maxMeta - metaNeural)];
            long mSum = mb + mn;
            long mHalf = mSum >> 1;
            long mInv = (1L << 34) / mSum;

            int f0 = Math.max(1, divFast((mb * bayesP[0] + mn * neuralP[0]) << 1, mSum, mHalf, mInv, 34));
            int f1 = Math.max(1, divFast((mb * bayesP[1] + mn * neuralP[1]) << 1, mSum, mHalf, mInv, 34));
            int f2 = Math.max(1, divFast((mb * bayesP[2] + mn * neuralP[2]) << 1, mSum, mHalf, mInv, 34));
            int f3 = Math.max(1, divFast((mb * bayesP[3] + mn * neuralP[3]) << 1, mSum, mHalf, mInv, 34));
            int fix = TOT - (f0 + f1 + f2 + f3);
            if (fix != 0) {
                if (f0 >= f1 && f0 >= f2 && f0 >= f3) f0 += fix;
                else if (f1 >= f2 && f1 >= f3) f1 += fix;
                else if (f2 >= f3) f2 += fix;
                else f3 += fix;
            }

            // Secondary Symbol Estimation (SSE) Calibration (12 KB L1 table)
            if (relPos >= 2) {
                int sBest = 0;
                int fMax = f0;
                if (f1 > fMax) {
                    sBest = 1;
                    fMax = f1;
                }
                if (f2 > fMax) {
                    sBest = 2;
                    fMax = f2;
                }
                if (f3 > fMax) {
                    sBest = 3;
                    fMax = f3;
                }
                int c2 = (int) (fwdHist & 15L);
                int bin = Math.min(7, Math.max(0, (fMax - 8192) / 3072));
                int sseState = (((phase * 16 + c2) << 2) | sBest) * 8 + bin;
                lastSseState = sseState;

                int sseBase = sseState << 2;
                int sc0 = sseCounts[sseBase];
                int sc1 = sseCounts[sseBase + 1];
                int sc2 = sseCounts[sseBase + 2];
                int sc3 = sseCounts[sseBase + 3];
                int scSum = sc0 + sc1 + sc2 + sc3;
                if (scSum > 0) {
                    int w = Math.min(8, scSum >> 1);
                    int sf0 = divRoundSmall((long) sc0 * TOT, scSum);
                    int sf1 = divRoundSmall((long) sc1 * TOT, scSum);
                    int sf2 = divRoundSmall((long) sc2 * TOT, scSum);
                    int sf3 = divRoundSmall((long) sc3 * TOT, scSum);

                    f0 = Math.max(1, ((32 - w) * f0 + w * sf0 + 16) >> 5);
                    f1 = Math.max(1, ((32 - w) * f1 + w * sf1 + 16) >> 5);
                    f2 = Math.max(1, ((32 - w) * f2 + w * sf2 + 16) >> 5);
                    f3 = Math.max(1, ((32 - w) * f3 + w * sf3 + 16) >> 5);

                    fix = TOT - (f0 + f1 + f2 + f3);
                    if (fix != 0) {
                        if (f0 >= f1 && f0 >= f2 && f0 >= f3) f0 += fix;
                        else if (f1 >= f2 && f1 >= f3) f1 += fix;
                        else if (f2 >= f3) f2 += fix;
                        else f3 += fix;
                    }
                }
            } else {
                lastSseState = -1;
            }

            outFreqs[0] = f0;
            outFreqs[1] = f1;
            outFreqs[2] = f2;
            outFreqs[3] = f3;
        }

        void update(byte[] buf, int segStart, int relPos, int sym) {
            int globalPos = segStart + relPos;

            metaBayes = ((metaBayes * 255L) >> 8) + ((3L * LN_LUT[bayesP[sym]]) >> 1);
            metaNeural = ((metaNeural * 255L) >> 8) + ((3L * LN_LUT[neuralP[sym]]) >> 1);

            for (int e = 0; e < numExperts; e++) {
                bayesScores[e] = ((bayesScores[e] * 255L) >> 8) + ((13L * expLnP[e][sym]) >> 3);
            }

            int r0 = (sym == 0 ? P_SCALE : 0) - neuralP[0];
            int r1 = (sym == 1 ? P_SCALE : 0) - neuralP[1];
            int r2 = (sym == 2 ? P_SCALE : 0) - neuralP[2];
            int r3 = (sym == 3 ? P_SCALE : 0) - neuralP[3];
            int[] rw = weights[activeRegime];
            for (int e = 0; e < numExperts; e++) {
                int[] lne = expLnP[e];
                long grad =
                        (long) r0 * lne[0]
                                + (long) r1 * lne[1]
                                + (long) r2 * lne[2]
                                + (long) r3 * lne[3];
                int updated = rw[e] - (rw[e] >> 11) + (int) (grad >> 20);
                if (updated < 0) updated = 0;
                else if (updated > 8192) updated = 8192;
                rw[e] = updated;
            }

            if (lastSseState >= 0) {
                int sseBase = lastSseState << 2;
                sseCounts[sseBase + sym] = (short) (sseCounts[sseBase + sym] + 2);
                int sSum =
                        sseCounts[sseBase]
                                + sseCounts[sseBase + 1]
                                + sseCounts[sseBase + 2]
                                + sseCounts[sseBase + 3];
                if (sSum >= 128) {
                    sseCounts[sseBase] = (short) (sseCounts[sseBase] >>> 1);
                    sseCounts[sseBase + 1] = (short) (sseCounts[sseBase + 1] >>> 1);
                    sseCounts[sseBase + 2] = (short) (sseCounts[sseBase + 2] >>> 1);
                    sseCounts[sseBase + 3] = (short) (sseCounts[sseBase + 3] >>> 1);
                }
            }

            // Step active repeat heads BEFORE inserting current k-mer into rings
            for (int h = 0; h < 6; h++) {
                heads[h].step(buf, segStart, globalPos, warmPrefixLimit, sym);
            }

            long prevFwd = fwdHist;
            rcHist = (rcHist >>> 2) | (((long) (sym ^ 3) & 3L) << 60);

            // Insert k=13 & k=18 forward/reverse repeat seeds using globalPos
            if (relPos >= K_FINE) {
                insertRingFwd(repFineFwd, repFineMask, prevFwd & MASK_FINE, globalPos);
                long rc13 = (rcHist >>> (62 - (K_FINE << 1))) & MASK_FINE;
                insertRingRev(repFineRev, repFineMask, rc13, globalPos - K_FINE);
            }
            if (relPos >= K_LONG) {
                insertRingFwd(repLongFwd, repLongMask, prevFwd & MASK_LONG, globalPos);
                long rc18 = (rcHist >>> (62 - (K_LONG << 1))) & MASK_LONG;
                insertRingRev(repLongRev, repLongMask, rc18, globalPos - K_LONG);
            }

            // Chargaff-Symmetric (+2 forward, +1 reverse-complement) Context Updates
            shortBump(ord0, 0, sym, 2, 16384);
            if (relPos >= 1) {
                shortBump(ord1, ((int) (prevFwd & 3L)) << 2, sym, 2, 16384);
            }
            if (relPos >= 2) {
                shortBump(ord2, ((int) (prevFwd & 15L)) << 2, sym, 2, 16384);
            }
            if (relPos >= 3) {
                shortBump(ord3, ((int) (prevFwd & 63L)) << 2, sym, 2, 8192);
            }
            if (relPos >= 4) {
                shortBump(ord4, ((int) (prevFwd & 255L)) << 2, sym, 2, 8192);
                int rcSym4 = ((int) (prevFwd >>> 6) & 3) ^ 3;
                int rcCtx4 = (int) ((rcHist >>> (62 - 8)) & 255L);
                shortBump(ord4, rcCtx4 << 2, rcSym4, 1, 8192);
            }
            if (relPos >= 5) {
                shortBump(ord5, ((int) (prevFwd & 1023L)) << 2, sym, 2, 4096);
                int rcSym5 = ((int) (prevFwd >>> 8) & 3) ^ 3;
                int rcCtx5 = (int) ((rcHist >>> (62 - 10)) & 1023L);
                shortBump(ord5, rcCtx5 << 2, rcSym5, 1, 4096);
            }
            if (relPos >= 6) {
                shortBump(ord6, ((int) (prevFwd & 4095L)) << 2, sym, 2, 4096);
                int rcSym6 = ((int) (prevFwd >>> 10) & 3) ^ 3;
                int rcCtx6 = (int) ((rcHist >>> (62 - 12)) & 4095L);
                shortBump(ord6, rcCtx6 << 2, rcSym6, 1, 4096);
            }
            if (relPos >= 7) {
                shortBump(ord7, ((int) (prevFwd & 16383L)) << 2, sym, 2, 2048);
                int rcSym7 = ((int) (prevFwd >>> 12) & 3) ^ 3;
                int rcCtx7 = (int) ((rcHist >>> (62 - 14)) & 16383L);
                shortBump(ord7, rcCtx7 << 2, rcSym7, 1, 2048);
            }

            int rmod3 = relPos % 3;
            for (int s = 0; s < 6; s++) {
                int c = rmod3 - (s < 3 ? s : s - 3);
                if (c < 0) c += 3;
                int p3s = (s < 3) ? pPlus3[c][sym] : pMinus3[c][sym];
                int p5s = (s < 3) ? pPlus5[c][sym] : pMinus5[c][sym];
                int pMix = (3 * p3s + 5 * p5s) >> 3;
                frameLog[s] = ((frameLog[s] * 255L) >> 8) + LN_LUT[pMix];
            }
            if (frameW[bestFrame] * 9 >= frameWsum * 4) {
                int c = rmod3 - (bestFrame < 3 ? bestFrame : bestFrame - 3);
                if (c < 0) c += 3;
                boolean plus = bestFrame < 3;
                if (relPos >= 2) {
                    shortBump(
                            plus ? orfPlus2 : orfMinus2,
                            ((((int) (prevFwd & 15L)) * 3) + c) << 2,
                            sym,
                            2,
                            4096);
                }
                if (relPos >= 3) {
                    shortBump(
                            plus ? orfPlus3 : orfMinus3,
                            ((((int) (prevFwd & 63L)) * 3) + c) << 2,
                            sym,
                            2,
                            4096);
                }
                if (relPos >= 5) {
                    shortBump(
                            plus ? orfPlus5 : orfMinus5,
                            ((((int) (prevFwd & 1023L)) * 3) + c) << 2,
                            sym,
                            2,
                            4096);
                }
            }

            if (relPos >= 11) {
                shortBump(wobble8, extractWobbleKey(prevFwd) << 2, sym, 2, 1024);
            }
            if (relPos >= 9) {
                shortBump(gapped8, extractGapped9Key(prevFwd) << 2, sym, 2, 1024);
                int rcSym9 = ((int) (prevFwd >>> 16) & 3) ^ 3;
                int rcG = extractGapped9Key(rcHist >>> (62 - 18));
                shortBump(gapped8, rcG << 2, rcSym9, 1, 1024);
            }

            if (useOrd8 && relPos >= 8) {
                int i8 = (int) (prevFwd & 65535L);
                ord8[i8] = swarBump(ord8[i8], sym, 2);
                int rcSym8 = ((int) (prevFwd >>> 14) & 3) ^ 3;
                int rcCtx8 = (int) ((rcHist >>> (62 - 16)) & 65535L);
                ord8[rcCtx8] = swarBump(ord8[rcCtx8], rcSym8, 1);
            }
            if (panSuper) {
                if (relPos >= 9) {
                    int i9 = (int) (prevFwd & 262143L);
                    ord9[i9] = swarBump(ord9[i9], sym, 2);
                    int rcSym9 = ((int) (prevFwd >>> 16) & 3) ^ 3;
                    int rcCtx9 = (int) ((rcHist >>> (62 - 18)) & 262143L);
                    ord9[rcCtx9] = swarBump(ord9[rcCtx9], rcSym9, 1);
                }
                if (relPos >= 12) {
                    updateSparse12Cached(sym, 2);
                    int rcSym12 = ((int) (prevFwd >>> 22) & 3) ^ 3;
                    long rcCtx12 = (rcHist >>> (62 - 24)) & 0xFFFFFFL;
                    updateSparse12Key(rcCtx12, rcSym12, 1);
                }
            }

            fwdHist = (prevFwd << 2) | (sym & 3L);
        }

        private void spawnRepeatCandidates(int segStart, int relPos, int globalPos) {
            if (relPos >= K_LONG) {
                long kmer18 = fwdHist & MASK_LONG;
                long h = mixSeed(kmer18);
                int slot = (int) h & repLongMask;
                int fSlot = slot << 1;
                int tag = (int) (h >>> 32);
                if (tag == 0) tag = 1;
                long e0 = repLongFwd[fSlot];
                if ((int) (e0 >>> 32) == tag) {
                    int c0 = (int) e0;
                    if (isValidRepeatCursor(c0, 1, segStart, globalPos) && !isCursorTracked(c0, 1)) {
                        allocateHeadSlot(0, 1, c0, 1, 56, 8);
                    }
                }
                long e1 = repLongFwd[fSlot + 1];
                if ((int) (e1 >>> 32) == tag) {
                    int c1 = (int) e1;
                    if (isValidRepeatCursor(c1, 1, segStart, globalPos) && !isCursorTracked(c1, 1)) {
                        allocateHeadSlot(0, 1, c1, 1, 56, 8);
                    }
                }
                long er = repLongRev[slot];
                if ((int) (er >>> 32) == tag) {
                    int cr = (int) er;
                    if (cr + K_LONG < globalPos
                            && isValidRepeatCursor(cr, -1, segStart, globalPos)
                            && !isCursorTracked(cr, -1)) {
                        allocateHeadSlot(1, 0, cr, -1, 56, 8);
                    }
                }
            }

            if (relPos >= K_FINE) {
                long kmer13 = fwdHist & MASK_FINE;
                long h = mixSeed(kmer13);
                int slot = (int) h & repFineMask;
                int fSlot = slot << 1;
                int tag = (int) (h >>> 32);
                if (tag == 0) tag = 1;
                long e0 = repFineFwd[fSlot];
                if ((int) (e0 >>> 32) == tag) {
                    int c0 = (int) e0;
                    if (isValidRepeatCursor(c0, 1, segStart, globalPos) && !isCursorTracked(c0, 1)) {
                        allocateHeadSlot(2, 3, c0, 1, 28, 6);
                    }
                }
                long e1 = repFineFwd[fSlot + 1];
                if ((int) (e1 >>> 32) == tag) {
                    int c1 = (int) e1;
                    if (isValidRepeatCursor(c1, 1, segStart, globalPos) && !isCursorTracked(c1, 1)) {
                        allocateHeadSlot(2, 3, c1, 1, 28, 6);
                    }
                }
                long er = repFineRev[slot];
                if ((int) (er >>> 32) == tag) {
                    int cr = (int) er;
                    if (cr + K_FINE < globalPos
                            && isValidRepeatCursor(cr, -1, segStart, globalPos)
                            && !isCursorTracked(cr, -1)) {
                        if (!heads[4].active || (heads[4].streak == 0 && heads[4].misses >= 2)) {
                            promoteVeteranIfStrong(heads[4]);
                            heads[4].activate(cr, -1, 28, 6);
                        }
                    }
                }
            }
        }

        private boolean isValidRepeatCursor(int cur, int dir, int segStart, int globalPos) {
            if (cur < 0 || cur >= globalPos) return false;
            if (segStart == 0) return true;
            // In parallel worker (segStart >= warmPrefixLimit), safe regions are [0 .. warmPrefixLimit) and [segStart .. globalPos)
            return cur < warmPrefixLimit || cur >= segStart;
        }

        private void allocateHeadSlot(
                int slotA, int slotB, int startCur, int dir, int initW, int maxMiss) {
            if (!heads[slotA].active) {
                heads[slotA].activate(startCur, dir, initW, maxMiss);
            } else if (!heads[slotB].active) {
                heads[slotB].activate(startCur, dir, initW, maxMiss);
            } else if (heads[slotB].misses >= 2 && heads[slotB].streak == 0) {
                promoteVeteranIfStrong(heads[slotB]);
                heads[slotB].activate(startCur, dir, initW, maxMiss);
            }
        }

        private void promoteVeteranIfStrong(MicroRepeatHead evicted) {
            if (evicted.active && evicted.totalMatched >= 20) {
                MicroRepeatHead vet = heads[5];
                if (!vet.active || evicted.totalMatched > vet.totalMatched) {
                    vet.copyFrom(evicted);
                    vet.maxMisses = 9;
                }
            }
        }

        private boolean isCursorTracked(int candCursor, int dir) {
            for (int h = 0; h < 6; h++) {
                MicroRepeatHead mh = heads[h];
                if (mh.active && mh.dir == dir && mh.cursor == candCursor) {
                    return true;
                }
            }
            return false;
        }

        private static void insertRingFwd(long[] ring, int mask, long kmer, int nextGlobalPos) {
            long h = mixSeed(kmer);
            int slot0 = ((int) h & mask) << 1;
            int slot1 = slot0 + 1;
            int tag = (int) (h >>> 32);
            if (tag == 0) tag = 1;
            long packed = (((long) tag) << 32) | (nextGlobalPos & 0xFFFFFFFFL);
            long e0 = ring[slot0];
            if ((int) (e0 >>> 32) == tag) {
                if (nextGlobalPos - (int) e0 >= 4) {
                    ring[slot1] = e0;
                }
                ring[slot0] = packed;
            } else {
                ring[slot1] = e0;
                ring[slot0] = packed;
            }
        }

        private static void insertRingRev(long[] ring, int mask, long rcKmer, int predBeforeSeed) {
            long h = mixSeed(rcKmer);
            int slot = (int) h & mask;
            int tag = (int) (h >>> 32);
            if (tag == 0) tag = 1;
            ring[slot] = (((long) tag) << 32) | (predBeforeSeed & 0xFFFFFFFFL);
        }

        private static long mixSeed(long x) {
            x ^= (x >>> 23);
            x *= 0xFF51AFD7ED558CCDL;
            x ^= (x >>> 29);
            return x;
        }

        private static int extractWobbleKey(long h) {
            int b01 = (int) (h & 0xFL);
            int b34 = (int) ((h >>> 6) & 0xFL);
            int b67 = (int) ((h >>> 12) & 0xFL);
            int b9A = (int) ((h >>> 18) & 0xFL);
            return b01 | (b34 << 4) | (b67 << 8) | (b9A << 12);
        }

        private static int extractGapped9Key(long h) {
            int lo4 = (int) (h & 0xFFL);
            int hi4 = (int) ((h >>> 10) & 0xFFL);
            return lo4 | (hi4 << 8);
        }

        private int predictSparse12(long ctx24, int[] priorP, int[] outP) {
            long mixed = Packed2BitSequence.mix64(ctx24);
            int idx0 = ((int) mixed & sparseMask) & ~1;
            int tag = (int) (mixed >>> 32);
            if (tag == 0) tag = 1;
            cachedSparseIdx = idx0;
            cachedSparseTag = tag;

            long e = sparse12[idx0];
            if ((int) (e >>> 32) != tag) {
                e = sparse12[idx0 + 1];
                if ((int) (e >>> 32) != tag) {
                    System.arraycopy(priorP, 0, outP, 0, 4);
                    return 0;
                }
            }
            int low = (int) e;
            swarPredict(low, priorP, 2, outP);
            return swarSum(low);
        }

        private void updateSparse12Cached(int sym, int inc) {
            updateSparse12Slot(cachedSparseIdx, cachedSparseTag, sym, inc);
        }

        private void updateSparse12Key(long ctx24, int sym, int inc) {
            long mixed = Packed2BitSequence.mix64(ctx24);
            int idx0 = ((int) mixed & sparseMask) & ~1;
            int tag = (int) (mixed >>> 32);
            if (tag == 0) tag = 1;
            updateSparse12Slot(idx0, tag, sym, inc);
        }

        private void updateSparse12Slot(int idx0, int tag, int sym, int inc) {
            int idx1 = idx0 + 1;
            long e0 = sparse12[idx0];
            int t0 = (int) (e0 >>> 32);
            if (t0 == tag || t0 == 0) {
                int low = swarBump(t0 == 0 ? 0 : (int) e0, sym, inc);
                sparse12[idx0] = (((long) tag) << 32) | (low & 0xFFFFFFFFL);
                return;
            }
            long e1 = sparse12[idx1];
            int t1 = (int) (e1 >>> 32);
            if (t1 == tag || t1 == 0) {
                int low = swarBump(t1 == 0 ? 0 : (int) e1, sym, inc);
                sparse12[idx1] = (((long) tag) << 32) | (low & 0xFFFFFFFFL);
                return;
            }
            int s0 = swarSum((int) e0);
            int s1 = swarSum((int) e1);
            if (s0 <= s1 && s0 <= 6) {
                sparse12[idx0] = (((long) tag) << 32) | (swarBump(0, sym, inc) & 0xFFFFFFFFL);
            } else if (s1 < s0 && s1 <= 6) {
                sparse12[idx1] = (((long) tag) << 32) | (swarBump(0, sym, inc) & 0xFFFFFFFFL);
            }
        }

        private static int divFast(long num, long denom, long half, long inv, int shift) {
            long n = num + half;
            long q = (n * inv) >>> shift;
            if (n - q * denom >= denom) {
                q++;
            }
            return (int) q;
        }

        private static int divRecip(long num, int denom, int half, long r) {
            long n = num + half;
            long q = (n * r) >>> 30;
            if (n - q * denom >= denom) {
                q++;
            }
            return (int) q;
        }

        private static int shortPredict(
                short[] tab, int baseIdx, int[] priorP, int tau, int[] outP) {
            int c0 = tab[baseIdx];
            int c1 = tab[baseIdx + 1];
            int c2 = tab[baseIdx + 2];
            int c3 = tab[baseIdx + 3];
            int sum = c0 + c1 + c2 + c3;
            if (sum == 0) {
                System.arraycopy(priorP, 0, outP, 0, 4);
                return 0;
            }
            int denom = (sum << 2) + tau;
            int half = denom >> 1;
            long r = RECIP[denom] & 0xFFFFFFFFL;
            int p0 =
                    Math.max(
                            1,
                            divRecip(
                                    ((long) c0 << 2) * P_SCALE + (long) tau * priorP[0],
                                    denom,
                                    half,
                                    r));
            int p1 =
                    Math.max(
                            1,
                            divRecip(
                                    ((long) c1 << 2) * P_SCALE + (long) tau * priorP[1],
                                    denom,
                                    half,
                                    r));
            int p2 =
                    Math.max(
                            1,
                            divRecip(
                                    ((long) c2 << 2) * P_SCALE + (long) tau * priorP[2],
                                    denom,
                                    half,
                                    r));
            int p3 = Math.max(1, P_SCALE - (p0 + p1 + p2));
            QuantumContextRepeatMixer.normalizeP(p0, p1, p2, p3, outP);
            return sum;
        }

        private static void shortBump(short[] tab, int baseIdx, int sym, int inc, int maxSum) {
            int ns = tab[baseIdx + sym] + inc;
            tab[baseIdx + sym] = (short) ns;
            int sum = tab[baseIdx] + tab[baseIdx + 1] + tab[baseIdx + 2] + tab[baseIdx + 3];
            if (sum >= maxSum) {
                tab[baseIdx] = (short) (tab[baseIdx] >>> 1);
                tab[baseIdx + 1] = (short) (tab[baseIdx + 1] >>> 1);
                tab[baseIdx + 2] = (short) (tab[baseIdx + 2] >>> 1);
                tab[baseIdx + 3] = (short) (tab[baseIdx + 3] >>> 1);
            }
        }

        private static void swarPredict(int packed, int[] priorP, int tau, int[] outP) {
            if (packed == 0) {
                System.arraycopy(priorP, 0, outP, 0, 4);
                return;
            }
            int c0 = packed & 0xFF;
            int c1 = (packed >>> 8) & 0xFF;
            int c2 = (packed >>> 16) & 0xFF;
            int c3 = (packed >>> 24) & 0xFF;
            int sum = c0 + c1 + c2 + c3;
            int denom = (sum << 2) + tau;
            int half = denom >> 1;
            long r = RECIP[denom] & 0xFFFFFFFFL;
            int p0 =
                    Math.max(
                            1,
                            divRecip(
                                    ((long) c0 << 2) * P_SCALE + (long) tau * priorP[0],
                                    denom,
                                    half,
                                    r));
            int p1 =
                    Math.max(
                            1,
                            divRecip(
                                    ((long) c1 << 2) * P_SCALE + (long) tau * priorP[1],
                                    denom,
                                    half,
                                    r));
            int p2 =
                    Math.max(
                            1,
                            divRecip(
                                    ((long) c2 << 2) * P_SCALE + (long) tau * priorP[2],
                                    denom,
                                    half,
                                    r));
            int p3 = Math.max(1, P_SCALE - (p0 + p1 + p2));
            QuantumContextRepeatMixer.normalizeP(p0, p1, p2, p3, outP);
        }

        private static int swarSum(int packed) {
            return (packed & 0xFF)
                    + ((packed >>> 8) & 0xFF)
                    + ((packed >>> 16) & 0xFF)
                    + ((packed >>> 24) & 0xFF);
        }

        private static int swarBump(int packed, int sym, int inc) {
            int shift = sym << 3;
            int c = (packed >>> shift) & 0xFF;
            if (c + inc >= 250) {
                packed = (packed >>> 1) & 0x7F7F7F7F;
            }
            return packed + (inc << shift);
        }
    }

    private static final class MicroRepeatHead {
        boolean active;
        int cursor;
        int dir;
        int misses;
        int maxMisses;
        int streak;
        int totalMatched;
        final int[] px = new int[4];
        final int[] xorCnt = {28, 2, 4, 2};

        void activate(int cursor, int dir, int initMatch, int maxMisses) {
            this.active = true;
            this.cursor = cursor;
            this.dir = dir;
            this.misses = 0;
            this.maxMisses = maxMisses;
            this.streak = 0;
            this.totalMatched = 0;
            this.xorCnt[0] = initMatch;
            this.xorCnt[1] = 2;
            this.xorCnt[2] = 4;
            this.xorCnt[3] = 2;
        }

        void copyFrom(MicroRepeatHead src) {
            this.active = src.active;
            this.cursor = src.cursor;
            this.dir = src.dir;
            this.misses = src.misses;
            this.maxMisses = src.maxMisses;
            this.streak = src.streak;
            this.totalMatched = src.totalMatched;
            System.arraycopy(src.xorCnt, 0, this.xorCnt, 0, 4);
        }

        private boolean cursorAccessible(int segStart, int globalPos, int warmLimit) {
            if (cursor < 0 || cursor >= globalPos) return false;
            if (segStart == 0) return true;
            return cursor < warmLimit || cursor >= segStart;
        }

        boolean predict(byte[] buf, int segStart, int globalPos, int warmLimit, int[] outP) {
            if (!active || !cursorAccessible(segStart, globalPos, warmLimit)) {
                active = false;
                return false;
            }
            int tpl = buf[cursor] & 3;
            if (dir < 0) tpl ^= 3;
            int cMatch =
                    xorCnt[0]
                            + Math.min(
                                    96,
                                    (streak << 1)
                                            + (misses == 0 ? Math.min(48, totalMatched) : 0));
            int c1 = xorCnt[1];
            int c2 = xorCnt[2];
            int c3 = xorCnt[3];
            int sum = cMatch + c1 + c2 + c3;
            px[0] = Math.max(1, divRoundSmall((long) cMatch * P_SCALE, sum));
            px[1] = Math.max(1, divRoundSmall((long) c1 * P_SCALE, sum));
            px[2] = Math.max(1, divRoundSmall((long) c2 * P_SCALE, sum));
            px[3] = Math.max(1, P_SCALE - (px[0] + px[1] + px[2]));
            QuantumContextRepeatMixer.normalizeP(
                    px[tpl], px[1 ^ tpl], px[2 ^ tpl], px[3 ^ tpl], outP);
            return true;
        }

        void step(byte[] buf, int segStart, int globalPos, int warmLimit, int actualSym) {
            if (!active) return;
            if (!cursorAccessible(segStart, globalPos, warmLimit)) {
                active = false;
                return;
            }
            int tpl = buf[cursor] & 3;
            if (dir < 0) tpl ^= 3;
            int delta = (actualSym ^ tpl) & 3;
            xorCnt[delta] += 8;
            if (xorCnt[0] + xorCnt[1] + xorCnt[2] + xorCnt[3] >= 256) {
                xorCnt[0] = Math.max(2, xorCnt[0] >>> 1);
                xorCnt[1] = Math.max(1, xorCnt[1] >>> 1);
                xorCnt[2] = Math.max(1, xorCnt[2] >>> 1);
                xorCnt[3] = Math.max(1, xorCnt[3] >>> 1);
            }
            if (delta == 0) {
                misses = 0;
                streak++;
                totalMatched++;
            } else {
                misses++;
                streak = 0;
                if (misses > maxMisses) {
                    active = false;
                    return;
                }
            }
            cursor += dir;
            if (!cursorAccessible(segStart, globalPos + 1, warmLimit)) {
                active = false;
            }
        }
    }

    /** 64-bit carry-safe Range Encoder matching QCRM's proven renormalization protocol. */
    private static final class RangeEnc {
        private final ByteOut out;
        private long low = 0L;
        private long range = 0xFFFFFFFFL;
        private int cache = -1;
        private long cacheSize = 0L;

        RangeEnc(ByteOut out) {
            this.out = out;
        }

        void encodeSymbol(int cumLow, int freq) {
            long r = range >>> PROB_BITS;
            low += cumLow * r;
            range = freq * r;
            while (range < TOP_VALUE) {
                shiftLow();
                range = (range << 8) & 0xFFFFFFFFL;
            }
        }

        private void shiftLow() {
            long carry = low >>> 32;
            if (carry != 0L || (low & 0xFFFFFFFFL) < 0xFF000000L) {
                if (cache >= 0) {
                    out.writeByte((int) (cache + carry));
                }
                int fill = (int) (0xFFL + carry);
                while (cacheSize > 0) {
                    out.writeByte(fill);
                    cacheSize--;
                }
                cache = (int) ((low >>> 24) & 0xFFL);
            } else {
                cacheSize++;
            }
            low = (low << 8) & 0xFFFFFFFFL;
        }

        void flush() {
            for (int i = 0; i < 5; i++) {
                shiftLow();
            }
        }
    }

    /** Matching Range Decoder. */
    private static final class RangeDec {
        private final ByteIn in;
        private long range = 0xFFFFFFFFL;
        private long code = 0L;

        RangeDec(ByteIn in) {
            this.in = in;
            for (int i = 0; i < 4; i++) {
                code = (code << 8) | readByteSafe();
            }
        }

        private long readByteSafe() {
            return in.hasRemaining() ? (in.readByte() & 0xFFL) : 0L;
        }

        int getThreshold() {
            long r = range >>> PROB_BITS;
            int t = (int) (code / r);
            return t >= TOT ? (TOT - 1) : t;
        }

        void decodeSymbol(int cumLow, int freq) {
            long r = range >>> PROB_BITS;
            code -= cumLow * r;
            range = freq * r;
            while (range < TOP_VALUE) {
                code = ((code << 8) | readByteSafe()) & 0xFFFFFFFFL;
                range = (range << 8) & 0xFFFFFFFFL;
            }
        }
    }
}
