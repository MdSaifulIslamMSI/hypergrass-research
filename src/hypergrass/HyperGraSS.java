package hypergrass;

import hypergrass.CanonicalSyncmerEngine.AnchorBank;
import hypergrass.FastaStreamFactorizer.BlockSlice;
import hypergrass.FastaStreamFactorizer.FactorizedChunk;
import hypergrass.FastaStreamFactorizer.RecordSpan;
import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;
import java.util.zip.CRC32;

/**
 * HyperGraSS: Next-Generation Lossless Genomic & Pangenomic Sequence Compression Engine.
 *
 * <p>Combines:
 * <ol>
 *   <li>6-channel universal lossless FASTA/RNA/DNA stream factorization (zero parity avalanche)</li>
 *   <li>BC-LZG bidirectional canonical 20-mer + homologous SNP-bridging referential engine</li>
 *   <li>BMM-RC allocation-free Bayesian multi-order Markov + codon-phase + RC range coder</li>
 *   <li>1-hop page-indexed random-access .hgs container with parallel block execution and
 *       strictly bounded (< 1 MB) --lean streaming mode</li>
 * </ol>
 */
public final class HyperGraSS {

    public static final String VERSION = "1.0.0";

    private static final byte[] MAGIC = {'H', 'G', 'S', 0x01};
    private static final int FOOTER_MAGIC = 0x31464748; // "HGF1" little-endian
    private static final int DEFAULT_BLOCK_BYTES = 524_288; // 512 KB
    private static final int LEAN_BLOCK_BYTES = 65_536;     // 64 KB bounded streaming chunk

    private static final int FLAG_HAS_ANCHOR = 1;
    private static final int FLAG_LEAN_MODE = 2;
    private static final int FLAG_COMPACT_SINGLE = 4;
    /** Blocks carry a mode-7 gap channel whose LZ dictionary spans block boundaries. */
    private static final int FLAG_SHARED_GAP = 8;
    private static final int GAP_HISTORY_BYTES = 256 * 1024;

    private HyperGraSS() {}

    public static final class CompressOptions {
        public boolean leanMode = false;
        public boolean parallel = true;
        public boolean optimize = false;
        public int blockSizeBytes = DEFAULT_BLOCK_BYTES;
        public int maxAnchorBases = CanonicalSyncmerEngine.DEFAULT_MAX_ANCHOR_BASES;
        public boolean anchorSinglePage = false;
        public boolean bounded = false;

        public static CompressOptions defaults() {
            return new CompressOptions();
        }

        public static CompressOptions lean() {
            CompressOptions o = new CompressOptions();
            o.leanMode = true;
            o.parallel = false;
            o.blockSizeBytes = LEAN_BLOCK_BYTES;
            o.maxAnchorBases = 0;
            return o;
        }
    }

    public static final class ArchiveMetadata {
        public final boolean compactSingleBlock;
        public final boolean hasAnchorBank;
        public final boolean leanMode;
        public final int numBlocks;
        public final long totalRawBytes;
        public final int totalRecords;
        public final int crc32;
        public final int anchorTotalBases;
        public final int numAnchorPages;

        public ArchiveMetadata(
                boolean compactSingleBlock,
                boolean hasAnchorBank,
                boolean leanMode,
                int numBlocks,
                long totalRawBytes,
                int totalRecords,
                int crc32,
                int anchorTotalBases,
                int numAnchorPages) {
            this.compactSingleBlock = compactSingleBlock;
            this.hasAnchorBank = hasAnchorBank;
            this.leanMode = leanMode;
            this.numBlocks = numBlocks;
            this.totalRawBytes = totalRawBytes;
            this.totalRecords = totalRecords;
            this.crc32 = crc32;
            this.anchorTotalBases = anchorTotalBases;
            this.numAnchorPages = numAnchorPages;
        }
    }

    private static final class BlockDirEntry {
        long fileOffset;
        int compressedLen;
        int firstRecordIdx;
        int recordCount;
        long firstRecordBaseOffset;
        int blockSeqChars;
    }

    /** Compresses an input file to a .hgs archive on disk. */
    public static void compressFile(Path input, Path output, CompressOptions opts) throws IOException {
        if (opts.leanMode) {
            compressStreamingLean(input, output);
            return;
        }
        if (opts.bounded) {
            compressFileBounded(input, output, opts);
            return;
        }
        byte[] raw = Files.readAllBytes(input);
        byte[] compressed = compressBytes(raw, opts);
        Files.write(output, compressed);
    }

    /**
     * Two-pass bounded-memory compressor for the anchor-bank format. Pass one slices and
     * factorizes the file one block at a time to build the anchor bank; pass two re-reads and
     * encodes each block. The input and the full chunk set are never resident together.
     */
    public static void compressFileBounded(Path input, Path output, CompressOptions opts) throws IOException {
        int blockBytes = (!opts.leanMode
                && opts.blockSizeBytes == DEFAULT_BLOCK_BYTES
                && Files.size(input) >= 16 * 1024 * 1024)
                ? (2 * 1024 * 1024)
                : opts.blockSizeBytes;
        FastaStreamFactorizer.FilePartition part = FastaStreamFactorizer.partitionFile(input, blockBytes);
        List<BlockSlice> slices = part.slices;
        int numBlocks = slices.size();

        try (FileChannel ch = FileChannel.open(input, StandardOpenOption.READ)) {
            CanonicalSyncmerEngine.AnchorBankBuilder builder =
                    new CanonicalSyncmerEngine.AnchorBankBuilder(opts.maxAnchorBases, part.size);
            long totalDna = 0L;
            int totalRecords = 0;
            for (BlockSlice s : slices) {
                FactorizedChunk chunk = factorizeSlice(ch, s);
                totalDna += chunk.dnaLen;
                totalRecords = Math.max(totalRecords, s.firstRecordIndex + s.recordCount);
                builder.addChunk(chunk);
            }
            if (part.size == 0) {
                totalRecords = 0;
            }
            boolean useAnchor = opts.maxAnchorBases > 0
                    && totalDna >= CanonicalSyncmerEngine.K * 2L;
            AnchorBank bank = useAnchor ? builder.finish() : null;
            // The builder shares the bank's seed index; dropping it here lets the index die before
            // polishing allocates a replacement, which is what keeps the bounded peak under 256 MiB.
            builder = null;
            // Match the in-memory path: polish the anchor to the cross-genome majority base.
            // Votes are accumulated by streaming the slices again, so no chunk is retained.
            // The seed index is built after the last vote buffer dies: building it inside the
            // pass would hold the vote array and the ~56 MB index alive at the same time.
            // The default packs four counters per position (18.9 MB on this bank, 4 bits each);
            // -Dhypergrass.boundedVotes8=true selects the exact 8-bit form (37.7 MB) for the
            // saturation ablation. Saturation at 15 votes costs 2,218 B on the pangenome.
            if (bank != null && bank.length > 0 && numBlocks > 1
                    && CanonicalSyncmerEngine.consensusEnabled()) {
                boolean exactVotes = "true".equals(
                        System.getProperty("hypergrass.boundedVotes8"));
                for (int pass = 0; pass < CanonicalSyncmerEngine.CONSENSUS_PASSES; pass++) {
                    byte[] votes = new byte[bank.length << (exactVotes ? 2 : 1)];
                    for (BlockSlice s : slices) {
                        FactorizedChunk chunk = factorizeSlice(ch, s);
                        if (exactVotes) {
                            CanonicalSyncmerEngine.accumulateAnchorVotes(
                                    bank, chunk.dnaBases, chunk.dnaLen, votes);
                        } else {
                            CanonicalSyncmerEngine.accumulateAnchorVotesPacked(
                                    bank, chunk.dnaBases, chunk.dnaLen, votes);
                        }
                    }
                    bank = CanonicalSyncmerEngine.polishAnchor(bank, votes, !exactVotes, false);
                }
                bank = CanonicalSyncmerEngine.withSeedIndex(bank);
            }

            List<BlockDirEntry> dir = new ArrayList<>(numBlocks);
            byte[][] encodedBlocks = new byte[numBlocks][];
            VarIntCodec.AuxHistory gapHistory =
                    new VarIntCodec.AuxHistory(GAP_HISTORY_BYTES);
            for (int i = 0; i < numBlocks; i++) {
                BlockSlice s = slices.get(i);
                FactorizedChunk chunk = factorizeSlice(ch, s);
                encodedBlocks[i] = encodeChunkPayload(chunk, bank, false, false, gapHistory);
                BlockDirEntry e = new BlockDirEntry();
                e.firstRecordIdx = s.firstRecordIndex;
                e.recordCount = s.recordCount;
                e.firstRecordBaseOffset = s.firstRecordBaseOffset;
                e.blockSeqChars = chunk.totalSeqChars;
                dir.add(e);
            }
            byte[] archive = serializeMultiBlockArchive(
                    part.size, part.crc32, totalRecords, totalDna, dir, encodedBlocks,
                    bank, bank != null && bank.length > 0, true, opts);
            Files.write(output, archive);
        }
    }

    private static FactorizedChunk factorizeSlice(FileChannel ch, BlockSlice slice) throws IOException {
        int len = slice.endByte - slice.startByte;
        if (len == 0) {
            return FastaStreamFactorizer.factorize(new byte[0], 0, 0);
        }
        ByteBuffer buf = ByteBuffer.allocate(len);
        readFully(ch, slice.startByte, buf);
        return FastaStreamFactorizer.factorize(buf.array(), 0, len);
    }

    /** Decompresses a .hgs archive on disk to output path, verifying CRC-32 integrity. */
    public static void decompressFile(Path archive, Path output) throws IOException {
        try (FileChannel ch = FileChannel.open(archive, StandardOpenOption.READ)) {
            ByteBuffer hdr = ByteBuffer.allocate(5);
            readFully(ch, 0L, hdr);
            hdr.flip();
            verifyMagic(hdr);
            int modeFlags = hdr.get() & 0xFF;
            if ((modeFlags & FLAG_LEAN_MODE) != 0 && (modeFlags & FLAG_COMPACT_SINGLE) == 0) {
                decompressStreaming(ch, output, modeFlags);
                return;
            }
            if ((modeFlags & FLAG_COMPACT_SINGLE) == 0) {
                decompressMultiBlockStreaming(ch, output, modeFlags);
                return;
            }
        }
        byte[] comp = Files.readAllBytes(archive);
        byte[] raw = decompressBytes(comp);
        Files.write(output, raw);
    }

    /**
     * Streaming decompressor for multi-block archives: decodes the anchor once, then decodes and
     * writes one block at a time. Peak heap is one compressed block plus its DNA and raw bytes.
     */
    private static void decompressMultiBlockStreaming(
            FileChannel ch, Path output, int modeFlags) throws IOException {
        boolean hasAnchor = (modeFlags & FLAG_HAS_ANCHOR) != 0;
        boolean leanMode = (modeFlags & FLAG_LEAN_MODE) != 0;
        boolean sharedGap = (modeFlags & FLAG_SHARED_GAP) != 0;
        ParsedFooter footer = readFooterFromChannel(ch);
        if (footer.entries.isEmpty()) {
            throw new IllegalStateException("Archive has no blocks");
        }
        long firstBlockOffset = footer.entries.get(0).fileOffset;

        byte[] anchorBases = null;
        int[][] globalDonors = CanonicalSyncmerEngine.EMPTY_DONORS;
        if (hasAnchor) {
            int headerLen = (int) (firstBlockOffset - 5L);
            ByteBuffer headerBuf = ByteBuffer.allocate(headerLen);
            readFully(ch, 5L, headerBuf);
            ByteIn in = new ByteIn(headerBuf.array());
            boolean fast2Bit = (in.readByte() == 0);
            int anchorLen = in.readVarInt();
            int numPages = in.readVarInt();
            globalDonors = CanonicalSyncmerEngine.readGlobalDonors(in).globalDonors;
            int[] pageBaseLens = new int[numPages];
            int[] pageByteLens = new int[numPages];
            for (int p = 0; p < numPages; p++) {
                pageBaseLens[p] = in.readVarInt();
                pageByteLens[p] = in.readVarInt();
            }
            anchorBases = new byte[anchorLen];
            int anchorOffset = 0;
            for (int p = 0; p < numPages; p++) {
                byte[] pagePayload = in.readBytes(pageByteLens[p]);
                byte[] dec = decodeAnchorPage(pagePayload, pageBaseLens[p], fast2Bit, leanMode);
                System.arraycopy(dec, 0, anchorBases, anchorOffset, pageBaseLens[p]);
                anchorOffset += pageBaseLens[p];
            }
        }

        CRC32 crc = new CRC32();
        long produced = 0L;
        VarIntCodec.AuxHistory gapHistory =
                sharedGap ? new VarIntCodec.AuxHistory(GAP_HISTORY_BYTES) : null;
        try (OutputStream out = new BufferedOutputStream(Files.newOutputStream(output), 1 << 20)) {
            for (BlockDirEntry e : footer.entries) {
                ByteBuffer blockBuf = ByteBuffer.allocate(e.compressedLen);
                readFully(ch, e.fileOffset, blockBuf);
                ByteIn bIn = new ByteIn(blockBuf.array());
                FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(bIn);
                byte[] dna = CanonicalSyncmerEngine.decodeBlockDna(
                        bIn, chunk.dnaLen, anchorBases, globalDonors, leanMode, gapHistory);
                byte[] recBytes = FastaStreamFactorizer.reconstructBytes(chunk, dna);
                crc.update(recBytes, 0, recBytes.length);
                out.write(recBytes);
                produced += recBytes.length;
            }
        }
        if (produced != footer.totalRawBytes) {
            throw new IllegalStateException(
                    "Raw byte count mismatch: expected " + footer.totalRawBytes + ", got " + produced);
        }
        if ((int) crc.getValue() != footer.crc32) {
            throw new IllegalStateException("CRC-32 integrity check failed!");
        }
    }

    /**
     * Strictly bounded-memory (< 1 MB live buffers) streaming compressor (--lean mode).
     * Never loads the full input or output file into heap memory; safe under -Xmx8m on 100+ MB genomes!
     */
    public static void compressStreamingLean(Path input, Path output) throws IOException {
        List<BlockDirEntry> dir = new ArrayList<>();
        CRC32 crc = new CRC32();
        long totalRawBytes = 0L;
        int globalRecordIdx = -1;
        boolean sawAnyBytes = false;
        boolean startsWithHeader = false;
        long currentRecordBaseAcc = 0L;

        try (BufferedInputStream in = new BufferedInputStream(Files.newInputStream(input), 32_768);
             BufferedOutputStream out = new BufferedOutputStream(Files.newOutputStream(output), 32_768)) {

            long writtenBytes = 0L;
            out.write(MAGIC);
            out.write(FLAG_LEAN_MODE); // multi-block, no anchor bank, lean mode
            writtenBytes += 5L;

            ByteArrayOutputStream chunkBuf = new ByteArrayOutputStream(LEAN_BLOCK_BYTES + 4096);
            ByteArrayOutputStream lineBuf = new ByteArrayOutputStream(256);

            int blockFirstRecord = 0;
            int blockRecordCount = 0;
            long blockFirstRecOffset = 0L;

            while (true) {
                lineBuf.reset();
                int b;
                boolean eof = false;
                while (true) {
                    b = in.read();
                    if (b == -1) {
                        eof = true;
                        break;
                    }
                    lineBuf.write(b);
                    if (b == '\n' || lineBuf.size() >= LEAN_BLOCK_BYTES) {
                        break;
                    }
                }

                byte[] line = lineBuf.toByteArray();
                if (line.length > 0) {
                    if (!sawAnyBytes) {
                        sawAnyBytes = true;
                        startsWithHeader = (line[0] == '>');
                        if (!startsWithHeader) {
                            globalRecordIdx = 0;
                            blockFirstRecord = 0;
                            blockRecordCount = 1;
                        }
                    }
                    crc.update(line, 0, line.length);
                    totalRawBytes += line.length;

                    boolean isHeader = (line[0] == '>' && (line.length > 1 || line[0] == '>'));
                    if (isHeader && chunkBuf.size() >= LEAN_BLOCK_BYTES && blockRecordCount > 0) {
                        byte[] rawBlock = chunkBuf.toByteArray();
                        chunkBuf.reset();
                        byte[] encBlock = encodeSingleChunkBytes(rawBlock, null, true);
                        out.write(encBlock);

                        BlockDirEntry e = new BlockDirEntry();
                        e.fileOffset = writtenBytes;
                        e.compressedLen = encBlock.length;
                        e.firstRecordIdx = blockFirstRecord;
                        e.recordCount = blockRecordCount;
                        e.firstRecordBaseOffset = blockFirstRecOffset;
                        e.blockSeqChars = countSeqCharsQuick(rawBlock);
                        dir.add(e);
                        writtenBytes += encBlock.length;

                        blockFirstRecord = globalRecordIdx + 1;
                        blockRecordCount = 0;
                        blockFirstRecOffset = 0L;
                    } else if (!isHeader && chunkBuf.size() >= LEAN_BLOCK_BYTES) {
                        byte[] rawBlock = chunkBuf.toByteArray();
                        chunkBuf.reset();
                        byte[] encBlock = encodeSingleChunkBytes(rawBlock, null, true);
                        out.write(encBlock);

                        BlockDirEntry e = new BlockDirEntry();
                        e.fileOffset = writtenBytes;
                        e.compressedLen = encBlock.length;
                        e.firstRecordIdx = Math.max(0, blockFirstRecord);
                        e.recordCount = Math.max(1, blockRecordCount);
                        e.firstRecordBaseOffset = blockFirstRecOffset;
                        int seqChars = countSeqCharsQuick(rawBlock);
                        e.blockSeqChars = seqChars;
                        dir.add(e);
                        writtenBytes += encBlock.length;

                        blockFirstRecord = Math.max(0, globalRecordIdx);
                        blockRecordCount = 1;
                        blockFirstRecOffset = currentRecordBaseAcc;
                    }

                    if (isHeader) {
                        globalRecordIdx++;
                        if (blockRecordCount == 0) {
                            blockFirstRecord = globalRecordIdx;
                            blockFirstRecOffset = 0L;
                        }
                        blockRecordCount++;
                        currentRecordBaseAcc = 0L;
                    } else {
                        int contentChars = line.length;
                        if (contentChars > 0 && line[contentChars - 1] == '\n') {
                            contentChars--;
                            if (contentChars > 0 && line[contentChars - 1] == '\r') {
                                contentChars--;
                            }
                        }
                        currentRecordBaseAcc += contentChars;
                    }
                    chunkBuf.write(line, 0, line.length);
                }

                if (eof) {
                    if (chunkBuf.size() > 0 || !sawAnyBytes) {
                        byte[] rawBlock = chunkBuf.toByteArray();
                        chunkBuf.reset();
                        byte[] encBlock = encodeSingleChunkBytes(rawBlock, null, true);
                        out.write(encBlock);

                        BlockDirEntry e = new BlockDirEntry();
                        e.fileOffset = writtenBytes;
                        e.compressedLen = encBlock.length;
                        e.firstRecordIdx = Math.max(0, blockFirstRecord);
                        e.recordCount = Math.max(1, blockRecordCount);
                        e.firstRecordBaseOffset = blockFirstRecOffset;
                        e.blockSeqChars = countSeqCharsQuick(rawBlock);
                        dir.add(e);
                        writtenBytes += encBlock.length;
                    }
                    break;
                }
            }

            int totalRecords = !sawAnyBytes ? 0 : (startsWithHeader ? (globalRecordIdx + 1) : Math.max(1, globalRecordIdx + 1));
            long footerOffset = writtenBytes;
            ByteOut footer = buildFooterBytes(dir, totalRawBytes, totalRecords, (int) crc.getValue(), footerOffset);
            out.write(footer.internalBuffer(), 0, footer.size());
        }
    }

    private static int countSeqCharsQuick(byte[] raw) {
        int p = 0;
        int end = raw.length;
        int chars = 0;
        while (p < end) {
            int lineStart = p;
            int nl = p;
            while (nl < end && raw[nl] != '\n') nl++;
            int contentEnd = (nl < end && nl > lineStart && raw[nl - 1] == '\r') ? (nl - 1) : nl;
            if (!(contentEnd > lineStart && raw[lineStart] == '>')) {
                chars += (contentEnd - lineStart);
            }
            p = (nl < end) ? (nl + 1) : end;
        }
        return chars;
    }

    /** Streaming bounded-memory block-by-block decompressor. */
    private static void decompressStreaming(FileChannel ch, Path output, int modeFlags) throws IOException {
        boolean leanMode = (modeFlags & FLAG_LEAN_MODE) != 0;
        ParsedFooter footer = readFooterFromChannel(ch);
        CRC32 crc = new CRC32();
        long produced = 0L;

        try (OutputStream out = new BufferedOutputStream(Files.newOutputStream(output), 32_768)) {
            for (BlockDirEntry b : footer.entries) {
                ByteBuffer buf = ByteBuffer.allocate(b.compressedLen);
                readFully(ch, b.fileOffset, buf);
                ByteIn blockIn = new ByteIn(buf.array());
                FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(blockIn);
                byte[] dna = CanonicalSyncmerEngine.decodeBlockDna(blockIn, chunk.dnaLen, null, leanMode);
                byte[] recBytes = FastaStreamFactorizer.reconstructBytes(chunk, dna);
                crc.update(recBytes, 0, recBytes.length);
                out.write(recBytes);
                produced += recBytes.length;
            }
        }
        if (produced != footer.totalRawBytes) {
            throw new IllegalStateException("Raw byte count mismatch: expected " + footer.totalRawBytes + ", got " + produced);
        }
        if ((int) crc.getValue() != footer.crc32) {
            throw new IllegalStateException("CRC-32 integrity check failed!");
        }
    }

    /** Compresses raw input bytes in-memory using full 6-channel + BC-LZG + BMM-RC + Parallelism. */
    public static byte[] compressBytes(byte[] raw, CompressOptions opts) {
        CRC32 crc = new CRC32();
        crc.update(raw, 0, raw.length);
        int crcVal = (int) crc.getValue();

        int effectiveBlockBytes = (!opts.leanMode
                && opts.blockSizeBytes == DEFAULT_BLOCK_BYTES
                && raw.length >= 16 * 1024 * 1024)
                ? (2 * 1024 * 1024)
                : opts.blockSizeBytes;
        List<BlockSlice> slices = FastaStreamFactorizer.partitionIntoBlocks(raw, effectiveBlockBytes);
        int numBlocks = slices.size();

        // Factorize all blocks (in parallel when numBlocks > 1 and opts.parallel)
        FactorizedChunk[] chunks = new FactorizedChunk[numBlocks];
        if (opts.parallel && numBlocks > 1) {
            IntStream.range(0, numBlocks).parallel().forEach(i -> {
                BlockSlice s = slices.get(i);
                chunks[i] = FastaStreamFactorizer.factorize(raw, s.startByte, s.endByte);
            });
        } else {
            for (int i = 0; i < numBlocks; i++) {
                BlockSlice s = slices.get(i);
                chunks[i] = FastaStreamFactorizer.factorize(raw, s.startByte, s.endByte);
            }
        }

        long totalDna = 0L;
        int totalRecords = 0;
        for (int i = 0; i < numBlocks; i++) {
            totalDna += chunks[i].dnaLen;
            BlockSlice s = slices.get(i);
            // A block can start inside a record when a long record is split at a
            // line boundary. Its recordCount still includes that in-progress
            // record, so the maximum over every block is the true record count.
            // Counting only blocks that start at a record boundary undercounts
            // records that begin inside a mid-record block.
            totalRecords = Math.max(totalRecords, s.firstRecordIndex + s.recordCount);
        }
        if (raw.length == 0) {
            totalRecords = 0;
        }

        // Single-genome / low-redundancy assembly path (<= 32 MB): if an AnchorBank does NOT find
        // substantial cross-record/self deduplication (bank covers > 25% of totalDna), encode the
        // entire genome in one unified QuantumContextRepeatMixer (QCRM) superblock so high-order
        // Markov, Spaced-Seed, Sparse-Hash, and Bidirectional Repeat models never reset at 512 KB boundaries!
        if (!opts.leanMode && raw.length <= 32 * 1024 * 1024) {
            if (numBlocks == 1 && chunks[0].getRecordSpans().length <= 1 && chunks[0].dnaLen < 65_536) {
                FactorizedChunk fullChunk = chunks[0];
                byte[] singlePayload = encodeChunkPayload(fullChunk, null, false);
                return wrapCompactSingleArchive(raw.length, crcVal, singlePayload);
            }
            // Probe AnchorBank compression density quickly
            boolean checkMulti = opts.maxAnchorBases > 0 && totalDna >= CanonicalSyncmerEngine.K * 2;
            AnchorBank probeBank = checkMulti
                    ? CanonicalSyncmerEngine.buildAnchorBank(chunks, opts.maxAnchorBases)
                    : null;
            boolean highlyRedundantCollection = (probeBank != null
                    && probeBank.length > 0
                    && (long) probeBank.length * 4L <= totalDna);

            if (!highlyRedundantCollection) {
                FactorizedChunk unified = (numBlocks == 1)
                        ? chunks[0]
                        : FastaStreamFactorizer.factorize(raw, 0, raw.length);
                byte[] unifiedPayload = encodeChunkPayload(unified, null, false);
                return wrapCompactSingleArchive(raw.length, crcVal, unifiedPayload);
            }

            // Highly redundant collection / pangenome / tandem-repeat chromosome: use prebuilt probeBank!
            byte[] multiArchive = buildMultiBlockArchiveWithBank(
                    raw.length, crcVal, slices, chunks, totalDna, totalRecords, opts, probeBank);
            if (numBlocks == 1) {
                byte[] singlePayload = encodeChunkPayload(chunks[0], null, false);
                if (13 + singlePayload.length <= multiArchive.length) {
                    return wrapCompactSingleArchive(raw.length, crcVal, singlePayload);
                }
            }
            return multiArchive;
        }

        return buildMultiBlockArchive(raw.length, crcVal, slices, chunks, totalDna, totalRecords, opts);
    }

    private static byte[] wrapCompactSingleArchive(int rawLen, int crcVal, byte[] singlePayload) {
        ByteOut out = new ByteOut(13 + singlePayload.length);
        out.writeBytes(MAGIC);
        out.writeByte(FLAG_COMPACT_SINGLE);
        out.writeIntLE(rawLen);
        out.writeIntLE(crcVal);
        out.writeBytes(singlePayload);
        return out.toByteArray();
    }

    private static byte[] buildMultiBlockArchive(
            long rawByteLen,
            int crcVal,
            List<BlockSlice> slices,
            FactorizedChunk[] chunks,
            long totalDna,
            int totalRecords,
            CompressOptions opts) {
        boolean useAnchor = !opts.leanMode && opts.maxAnchorBases > 0 && totalDna >= CanonicalSyncmerEngine.K * 2;
        AnchorBank bank = useAnchor
                ? CanonicalSyncmerEngine.buildAnchorBank(chunks, opts.maxAnchorBases)
                : null;
        return buildMultiBlockArchiveWithBank(
                rawByteLen, crcVal, slices, chunks, totalDna, totalRecords, opts, bank);
    }

    private static byte[] buildMultiBlockArchiveWithBank(
            long rawByteLen,
            int crcVal,
            List<BlockSlice> slices,
            FactorizedChunk[] chunks,
            long totalDna,
            int totalRecords,
            CompressOptions opts,
            AnchorBank bank) {
        byte[] best = encodeMultiBlockArchiveWithBank(
                rawByteLen, crcVal, slices, chunks, totalDna, totalRecords, opts, bank, false);
        // Global reuse must pay for the dictionary itself, not just shorten individual blocks.
        // Bound the extra full encoding to small archives and preserve the old choice on ties.
        if (opts.optimize && !opts.leanMode && rawByteLen <= 2 * 1024 * 1024 && slices.size() <= 16
                && bank != null && bank.length <= CanonicalSyncmerEngine.ANCHOR_PAGE_SIZE
                && bank.globalDonors.length > 0) {
            byte[] withDictionary = encodeMultiBlockArchiveWithBank(
                    rawByteLen, crcVal, slices, chunks, totalDna, totalRecords, opts, bank, true);
            if (withDictionary.length < best.length) best = withDictionary;
            AnchorBank localOnly = new AnchorBank(bank.bases, bank.length, bank.primarySeedLen,
                    CanonicalSyncmerEngine.EMPTY_DONORS, bank.index);
            byte[] candidate = encodeMultiBlockArchiveWithBank(
                    rawByteLen, crcVal, slices, chunks, totalDna, totalRecords, opts, localOnly, true);
            if (candidate.length < best.length) return candidate;
        }
        return best;
    }

    private static byte[] encodeMultiBlockArchiveWithBank(
            long rawByteLen,
            int crcVal,
            List<BlockSlice> slices,
            FactorizedChunk[] chunks,
            long totalDna,
            int totalRecords,
            CompressOptions opts,
            AnchorBank bank,
            boolean costAware) {
        int numBlocks = slices.size();
        boolean useAnchor = (bank != null && bank.length > 0);
        if (!useAnchor) {
            bank = null;
        }

        // Encode all blocks in parallel against the immutable AnchorBank + Global Clade Catalog!
        final AnchorBank finalBank = bank;
        byte[][] encodedBlocks = new byte[numBlocks][];
        if (opts.parallel && numBlocks > 1) {
            IntStream.range(0, numBlocks).parallel().forEach(i -> {
                encodedBlocks[i] = encodeChunkPayload(chunks[i], finalBank, opts.leanMode, costAware);
            });
        } else {
            for (int i = 0; i < numBlocks; i++) {
                encodedBlocks[i] = encodeChunkPayload(chunks[i], finalBank, opts.leanMode, costAware);
            }
        }

        List<BlockDirEntry> dir = new ArrayList<>(numBlocks);
        for (int i = 0; i < numBlocks; i++) {
            BlockSlice s = slices.get(i);
            BlockDirEntry e = new BlockDirEntry();
            e.firstRecordIdx = s.firstRecordIndex;
            e.recordCount = s.recordCount;
            e.firstRecordBaseOffset = s.firstRecordBaseOffset;
            e.blockSeqChars = chunks[i].totalSeqChars;
            dir.add(e);
        }
        return serializeMultiBlockArchive(
                rawByteLen, crcVal, totalRecords, totalDna, dir, encodedBlocks, bank, useAnchor,
                false, opts);
    }

    private static byte[] serializeMultiBlockArchive(
            long rawByteLen,
            int crcVal,
            int totalRecords,
            long totalDna,
            List<BlockDirEntry> dir,
            byte[][] encodedBlocks,
            AnchorBank bank,
            boolean useAnchor,
            boolean sharedGap,
            CompressOptions opts) {
        int numBlocks = encodedBlocks.length;
        int modeFlags = 0;
        if (useAnchor) modeFlags |= FLAG_HAS_ANCHOR;
        if (opts.leanMode) modeFlags |= FLAG_LEAN_MODE;
        if (sharedGap) modeFlags |= FLAG_SHARED_GAP;

        ByteOut out = new ByteOut(4096);
        out.writeBytes(MAGIC);
        out.writeByte(modeFlags);

        // Serialize page-partitioned AnchorBank + 1-hop Global Clade Catalog
        if (useAnchor) {
            int anchorLen = bank.length;
            int pageSize = opts.anchorSinglePage
                    ? Math.max(1, anchorLen)
                    : CanonicalSyncmerEngine.ANCHOR_PAGE_SIZE;
            int numPages = (anchorLen + pageSize - 1) / pageSize;
            final boolean preferFast2Bit = (totalDna >= 12L * anchorLen);

            byte[][] encPages = new byte[numPages][];
            int[] pageBaseLens = new int[numPages];
            if (opts.parallel && numPages > 1) {
                IntStream.range(0, numPages).parallel().forEach(p -> {
                    int pStart = p * pageSize;
                    int pLen = Math.min(pageSize, anchorLen - pStart);
                    pageBaseLens[p] = pLen;
                    encPages[p] = encodeAnchorPage(bank.bases, pStart, pLen, preferFast2Bit, opts.leanMode);
                });
            } else {
                for (int p = 0; p < numPages; p++) {
                    int pStart = p * pageSize;
                    int pLen = Math.min(pageSize, anchorLen - pStart);
                    pageBaseLens[p] = pLen;
                    encPages[p] = encodeAnchorPage(bank.bases, pStart, pLen, preferFast2Bit, opts.leanMode);
                }
            }

            boolean fast2BitAnchor = preferFast2Bit;
            // Bound the extra context-mixing work to one page; reuse alone does not predict entropy.
            // A smaller payload also has a no-larger length varint. Keep the legacy choice on ties.
            if (!opts.leanMode && numPages == 1) {
                byte[] alternate = encodeAnchorPage(bank.bases, 0, anchorLen, !preferFast2Bit, false);
                if (alternate.length < encPages[0].length) {
                    encPages[0] = alternate;
                    fast2BitAnchor = !preferFast2Bit;
                }
            }
            out.writeByte(fast2BitAnchor ? 0 : 1);
            out.writeVarInt(anchorLen);
            out.writeVarInt(numPages);
            CanonicalSyncmerEngine.writeGlobalDonors(out, bank.primarySeedLen, bank.globalDonors);

            for (int p = 0; p < numPages; p++) {
                out.writeVarInt(pageBaseLens[p]);
                out.writeVarInt(encPages[p].length);
            }
            for (int p = 0; p < numPages; p++) {
                out.writeBytes(encPages[p]);
            }
        }

        for (int i = 0; i < numBlocks; i++) {
            BlockDirEntry e = dir.get(i);
            e.fileOffset = out.size();
            e.compressedLen = encodedBlocks[i].length;
            out.writeBytes(encodedBlocks[i]);
        }

        long footerOffset = out.size();
        ByteOut footer = buildFooterBytes(dir, rawByteLen, totalRecords, crcVal, footerOffset);
        out.writeBytes(footer.internalBuffer(), 0, footer.size());
        return out.toByteArray();
    }

    private static byte[] encodeAnchorPage(
            byte[] anchorBases, int offset, int len, boolean fast2Bit, boolean leanMode) {
        if (fast2Bit) {
            byte[] packed = VarIntCodec.pack2Bit(anchorBases, offset, len);
            ByteOut tmp = new ByteOut(packed.length + 16);
            VarIntCodec.writeCompressedAux(tmp, packed, packed.length);
            return tmp.toByteArray();
        } else {
            byte[] slice = (offset == 0 && len == anchorBases.length)
                    ? anchorBases
                    : Arrays.copyOfRange(anchorBases, offset, offset + len);
            return QuantumContextRepeatMixer.encode(slice, len, leanMode);
        }
    }

    private static byte[] decodeAnchorPage(byte[] pagePayload, int pageBaseLen, boolean fast2Bit, boolean leanMode) {
        ByteIn in = new ByteIn(pagePayload);
        if (fast2Bit) {
            byte[] packed = VarIntCodec.readCompressedAux(in);
            return VarIntCodec.unpack2Bit(packed, pageBaseLen);
        } else {
            return QuantumContextRepeatMixer.decode(in, pageBaseLen, leanMode);
        }
    }

    private static byte[] encodeSingleChunkBytes(byte[] rawChunk, AnchorBank bank, boolean leanMode) {
        FactorizedChunk chunk = FastaStreamFactorizer.factorize(rawChunk, 0, rawChunk.length);
        return encodeChunkPayload(chunk, bank, leanMode);
    }

    private static byte[] encodeChunkPayload(FactorizedChunk chunk, AnchorBank bank, boolean leanMode) {
        return encodeChunkPayload(chunk, bank, leanMode, false, null);
    }

    private static byte[] encodeChunkPayload(
            FactorizedChunk chunk, AnchorBank bank, boolean leanMode, boolean costAware) {
        return encodeChunkPayload(chunk, bank, leanMode, costAware, null);
    }

    private static byte[] encodeChunkPayload(
            FactorizedChunk chunk, AnchorBank bank, boolean leanMode, boolean costAware,
            VarIntCodec.AuxHistory gapHistory) {
        ByteOut out = new ByteOut(256);
        FastaStreamFactorizer.writeAuxMetadata(out, chunk);
        byte[] dnaEnc = CanonicalSyncmerEngine.encodeBlockDna(
                chunk.dnaBases, chunk.dnaLen, bank, leanMode, costAware, gapHistory);
        out.writeBytes(dnaEnc);
        return out.toByteArray();
    }

    /** Decompresses a full .hgs archive in memory and verifies exact CRC-32 checksum. */
    public static byte[] decompressBytes(byte[] archive) {
        ByteIn in = new ByteIn(archive);
        if (archive.length < 13) {
            throw new IllegalArgumentException("Truncated .hgs archive");
        }
        for (int i = 0; i < 4; i++) {
            if ((byte) in.readByte() != MAGIC[i]) {
                throw new IllegalArgumentException("Invalid .hgs magic header");
            }
        }
        int modeFlags = in.readByte();
        boolean hasAnchor = (modeFlags & FLAG_HAS_ANCHOR) != 0;
        boolean leanMode = (modeFlags & FLAG_LEAN_MODE) != 0;
        boolean compactSingle = (modeFlags & FLAG_COMPACT_SINGLE) != 0;
        boolean sharedGap = (modeFlags & FLAG_SHARED_GAP) != 0;

        if (compactSingle) {
            int expectedRawLen = in.readIntLE();
            int expectedCrc = in.readIntLE();
            FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(in);
            byte[] dna = CanonicalSyncmerEngine.decodeBlockDna(in, chunk.dnaLen, null, leanMode);
            byte[] raw = FastaStreamFactorizer.reconstructBytes(chunk, dna);
            verifyIntegrity(raw, expectedRawLen, expectedCrc);
            return raw;
        }

        // Read Footer Directory from trailer
        ParsedFooter footer = parseFooterFromBytes(archive);

        // Decode AnchorBank + 1-hop Global Clade Catalog (if present) in parallel across pages
        byte[] anchorBases = null;
        int[][] globalDonors = CanonicalSyncmerEngine.EMPTY_DONORS;
        if (hasAnchor) {
            boolean fast2Bit = (in.readByte() == 0);
            int anchorLen = in.readVarInt();
            int numPages = in.readVarInt();
            globalDonors = CanonicalSyncmerEngine.readGlobalDonors(in).globalDonors;
            int[] pageBaseLens = new int[numPages];
            int[] pageByteLens = new int[numPages];
            for (int p = 0; p < numPages; p++) {
                pageBaseLens[p] = in.readVarInt();
                pageByteLens[p] = in.readVarInt();
            }
            byte[][] pageBytes = new byte[numPages][];
            for (int p = 0; p < numPages; p++) {
                pageBytes[p] = in.readBytes(pageByteLens[p]);
            }
            anchorBases = new byte[anchorLen];
            final byte[] targetAnchor = anchorBases;
            final int[] pageOffsets = new int[numPages + 1];
            for (int p = 0; p < numPages; p++) {
                pageOffsets[p + 1] = pageOffsets[p] + pageBaseLens[p];
            }
            if (numPages > 1) {
                IntStream.range(0, numPages).parallel().forEach(p -> {
                    byte[] dec = decodeAnchorPage(pageBytes[p], pageBaseLens[p], fast2Bit, leanMode);
                    System.arraycopy(dec, 0, targetAnchor, pageOffsets[p], pageBaseLens[p]);
                });
            } else if (numPages == 1) {
                byte[] dec = decodeAnchorPage(pageBytes[0], pageBaseLens[0], fast2Bit, leanMode);
                System.arraycopy(dec, 0, targetAnchor, 0, pageBaseLens[0]);
            }
        }

        // Decode all blocks in parallel
        int numBlocks = footer.entries.size();
        byte[][] decodedBlocks = new byte[numBlocks][];
        final byte[] finalAnchor = anchorBases;
        final int[][] finalGlobalDonors = globalDonors;
        if (sharedGap) {
            VarIntCodec.AuxHistory gapHistory = new VarIntCodec.AuxHistory(GAP_HISTORY_BYTES);
            for (int b = 0; b < numBlocks; b++) {
                BlockDirEntry e = footer.entries.get(b);
                ByteIn bIn = new ByteIn(archive, (int) e.fileOffset, e.compressedLen);
                FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(bIn);
                byte[] dna = CanonicalSyncmerEngine.decodeBlockDna(
                        bIn, chunk.dnaLen, finalAnchor, finalGlobalDonors, leanMode, gapHistory);
                decodedBlocks[b] = FastaStreamFactorizer.reconstructBytes(chunk, dna);
            }
        } else if (numBlocks > 1) {
            IntStream.range(0, numBlocks).parallel().forEach(b -> {
                BlockDirEntry e = footer.entries.get(b);
                ByteIn bIn = new ByteIn(archive, (int) e.fileOffset, e.compressedLen);
                FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(bIn);
                byte[] dna = CanonicalSyncmerEngine.decodeBlockDna(
                        bIn, chunk.dnaLen, finalAnchor, finalGlobalDonors, leanMode);
                decodedBlocks[b] = FastaStreamFactorizer.reconstructBytes(chunk, dna);
            });
        } else if (numBlocks == 1) {
            BlockDirEntry e = footer.entries.get(0);
            ByteIn bIn = new ByteIn(archive, (int) e.fileOffset, e.compressedLen);
            FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(bIn);
            byte[] dna = CanonicalSyncmerEngine.decodeBlockDna(
                    bIn, chunk.dnaLen, finalAnchor, finalGlobalDonors, leanMode);
            decodedBlocks[0] = FastaStreamFactorizer.reconstructBytes(chunk, dna);
        }

        if (footer.totalRawBytes > Integer.MAX_VALUE) {
            throw new IllegalStateException("Archive > 2 GB; use decompressFile streaming API");
        }
        byte[] full = new byte[(int) footer.totalRawBytes];
        int pos = 0;
        for (int b = 0; b < numBlocks; b++) {
            byte[] part = decodedBlocks[b];
            System.arraycopy(part, 0, full, pos, part.length);
            pos += part.length;
        }
        verifyIntegrity(full, footer.totalRawBytes, footer.crc32);
        return full;
    }

    /**
     * O(1) 1-Hop Random-Access Record & Subsequence Query directly from disk via FileChannel!
     *
     * @param archive path to .hgs archive
     * @param recordIndex 0-based FASTA record index
     * @param qStart 0-based start base within record, or -1 to extract the entire formatted FASTA record
     * @param qEnd 0-based exclusive end base within record, or -1 for end of record
     */
    public static byte[] query(Path archive, int recordIndex, int qStart, int qEnd) throws IOException {
        try (FileChannel ch = FileChannel.open(archive, StandardOpenOption.READ)) {
            long fileSize = ch.size();
            ByteBuffer hdr = ByteBuffer.allocate((int) Math.min(fileSize, 4096));
            readFully(ch, 0L, hdr);
            hdr.flip();
            verifyMagic(hdr);
            int modeFlags = hdr.get() & 0xFF;
            boolean hasAnchor = (modeFlags & FLAG_HAS_ANCHOR) != 0;
            boolean leanMode = (modeFlags & FLAG_LEAN_MODE) != 0;
            boolean compactSingle = (modeFlags & FLAG_COMPACT_SINGLE) != 0;

            if (compactSingle) {
                byte[] all = Files.readAllBytes(archive);
                ByteIn in = new ByteIn(all, 13, all.length - 13);
                FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(in);
                byte[] dna = CanonicalSyncmerEngine.decodeBlockDna(in, chunk.dnaLen, null, leanMode);
                return FastaStreamFactorizer.extractRecordOrSlice(chunk, dna, recordIndex, qStart, qEnd);
            }

            ParsedFooter footer = readFooterFromChannel(ch);
            if (recordIndex < 0 || recordIndex >= footer.totalRecords) {
                throw new IndexOutOfBoundsException(
                        "Record index " + recordIndex + " out of [0, " + footer.totalRecords + ")");
            }

            // Parse AnchorBank page table + Global Clade Catalog header (only metadata—do NOT read unreferenced pages!)
            boolean fast2BitAnchor = false;
            int anchorLen = 0;
            int numAnchorPages = 0;
            int[][] globalDonors = CanonicalSyncmerEngine.EMPTY_DONORS;
            int[] pageBaseLens = null;
            int[] pageByteLens = null;
            long[] pageFileOffsets = null;

            if (hasAnchor) {
                int anchorHdrCap = (int) Math.min(fileSize - 5L, 262_144L);
                ByteBuffer abuf = ByteBuffer.allocate(anchorHdrCap);
                readFully(ch, 5L, abuf);
                ByteIn aIn = new ByteIn(abuf.array());
                fast2BitAnchor = (aIn.readByte() == 0);
                anchorLen = aIn.readVarInt();
                numAnchorPages = aIn.readVarInt();
                globalDonors = CanonicalSyncmerEngine.readGlobalDonors(aIn).globalDonors;
                pageBaseLens = new int[numAnchorPages];
                pageByteLens = new int[numAnchorPages];
                pageFileOffsets = new long[numAnchorPages];
                for (int p = 0; p < numAnchorPages; p++) {
                    pageBaseLens[p] = aIn.readVarInt();
                    pageByteLens[p] = aIn.readVarInt();
                }
                long curPageOff = 5L + aIn.position();
                for (int p = 0; p < numAnchorPages; p++) {
                    pageFileOffsets[p] = curPageOff;
                    curPageOff += pageByteLens[p];
                }
            }

            // Binary-search first block containing recordIndex
            List<BlockDirEntry> entries = footer.entries;
            int firstBlock = findFirstBlockForRecord(entries, recordIndex);

            // Check whether recordIndex is contained inside a single block (99.99% of multi-FASTA queries)
            BlockDirEntry b0 = entries.get(firstBlock);
            boolean spansMultipleBlocks = (firstBlock + 1 < entries.size()
                    && entries.get(firstBlock + 1).firstRecordIdx == recordIndex
                    && entries.get(firstBlock + 1).firstRecordBaseOffset > 0L);

            if (!spansMultipleBlocks) {
                int localIdx = recordIndex - b0.firstRecordIdx;
                byte[] blockBytes = readBlockBytes(ch, b0);
                ByteIn bIn = new ByteIn(blockBytes);
                FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(bIn);

                byte[] anchorBases = null;
                if (hasAnchor) {
                    boolean[] needed = CanonicalSyncmerEngine.findReferencedAnchorPages(
                            bIn, numAnchorPages, leanMode);
                    anchorBases = loadSelectedAnchorPages(
                            ch, anchorLen, numAnchorPages, needed,
                            pageFileOffsets, pageByteLens, pageBaseLens, fast2BitAnchor, leanMode);
                }
                return FastaStreamFactorizer.extractRecordOrSliceSelective(
                        chunk, bIn, anchorBases, globalDonors, leanMode, localIdx, qStart, qEnd);
            }

            // Multi-block giant chromosome record: touch ONLY the blocks overlapping [qStart, qEnd)!
            ByteOut multiOut = new ByteOut(qStart >= 0 && qEnd >= qStart ? (qEnd - qStart) : 65536);
            byte[] anchorBases = null;
            for (int bi = firstBlock; bi < entries.size(); bi++) {
                BlockDirEntry be = entries.get(bi);
                if (bi > firstBlock && (be.firstRecordIdx != recordIndex || be.firstRecordBaseOffset == 0L)) {
                    break;
                }
                int localIdx = (bi == firstBlock) ? (recordIndex - be.firstRecordIdx) : 0;
                long recStartInGenome = (bi == firstBlock && localIdx > 0) ? 0L : be.firstRecordBaseOffset;
                long approxEndInGenome = recStartInGenome + be.blockSeqChars;

                if (qStart >= 0) {
                    if (approxEndInGenome <= qStart) continue;
                    if (qEnd >= 0 && recStartInGenome >= qEnd) break;
                }

                byte[] blockBytes = readBlockBytes(ch, be);
                ByteIn bIn = new ByteIn(blockBytes);
                FactorizedChunk chunk = FastaStreamFactorizer.readAuxMetadata(bIn);
                if (hasAnchor) {
                    boolean[] needed = CanonicalSyncmerEngine.findReferencedAnchorPages(
                            bIn, numAnchorPages, leanMode);
                    if (anchorBases == null) anchorBases = new byte[anchorLen];
                    populateSelectedAnchorPages(
                            ch, anchorBases, numAnchorPages, needed,
                            pageFileOffsets, pageByteLens, pageBaseLens, fast2BitAnchor, leanMode);
                }
                byte[] dna = CanonicalSyncmerEngine.decodeBlockDna(
                        bIn, chunk.dnaLen, anchorBases, globalDonors, leanMode);
                if (qStart >= 0) {
                    int localS = (int) Math.max(0L, qStart - recStartInGenome);
                    int localE = (qEnd < 0)
                            ? -1
                            : (int) Math.max(0L, qEnd - recStartInGenome);
                    byte[] sub = FastaStreamFactorizer.extractRecordOrSlice(chunk, dna, localIdx, localS, localE);
                    multiOut.writeBytes(sub);
                } else {
                    RecordSpan span = chunk.getRecordSpans()[localIdx];
                    if (bi == firstBlock && span.headerBytes != null) {
                        multiOut.writeByte('>');
                        multiOut.writeBytes(span.headerBytes);
                        multiOut.writeByte('\n');
                    }
                    byte[] seqChars = FastaStreamFactorizer.reconstructSequenceChars(chunk, dna);
                    int cur = span.seqCharStart;
                    for (int line = span.firstSeqLine; line < span.endSeqLine; line++) {
                        int len = chunk.seqLineLengths[line];
                        if (len > 0) {
                            multiOut.writeBytes(seqChars, cur, len);
                            cur += len;
                        }
                        multiOut.writeByte('\n');
                    }
                }
            }
            return multiOut.toByteArray();
        }
    }

    /** Inspects archive structural metadata in O(1) without decompressing blocks. */
    public static ArchiveMetadata inspect(Path archive) throws IOException {
        try (FileChannel ch = FileChannel.open(archive, StandardOpenOption.READ)) {
            ByteBuffer hdr = ByteBuffer.allocate(13);
            readFully(ch, 0L, hdr);
            hdr.flip();
            verifyMagic(hdr);
            int modeFlags = hdr.get() & 0xFF;
            boolean hasAnchor = (modeFlags & FLAG_HAS_ANCHOR) != 0;
            boolean leanMode = (modeFlags & FLAG_LEAN_MODE) != 0;
            boolean compactSingle = (modeFlags & FLAG_COMPACT_SINGLE) != 0;

            if (compactSingle) {
                ByteIn in = new ByteIn(hdr.array(), 5, 8);
                int rawLen = in.readIntLE();
                int crc = in.readIntLE();
                return new ArchiveMetadata(true, false, leanMode, 1, rawLen, rawLen > 0 ? 1 : 0, crc, 0, 0);
            }

            ParsedFooter footer = readFooterFromChannel(ch);
            int anchorLen = 0;
            int numPages = 0;
            if (hasAnchor) {
                ByteBuffer abuf = ByteBuffer.allocate((int) Math.min(64L, ch.size() - 5L));
                readFully(ch, 5L, abuf);
                ByteIn aIn = new ByteIn(abuf.array());
                aIn.readByte();
                anchorLen = aIn.readVarInt();
                numPages = aIn.readVarInt();
            }
            return new ArchiveMetadata(
                    false, hasAnchor, leanMode, footer.entries.size(),
                    footer.totalRawBytes, footer.totalRecords, footer.crc32, anchorLen, numPages);
        }
    }

    private static int findFirstBlockForRecord(List<BlockDirEntry> entries, int recordIndex) {
        int lo = 0;
        int hi = entries.size() - 1;
        int ans = 0;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            BlockDirEntry e = entries.get(mid);
            int lastRecInBlock = e.firstRecordIdx + e.recordCount - 1;
            if (lastRecInBlock < recordIndex) {
                lo = mid + 1;
            } else if (e.firstRecordIdx > recordIndex) {
                hi = mid - 1;
            } else {
                ans = mid;
                // Walk left if giant chromosome spans earlier blocks with same recordIndex
                hi = mid - 1;
            }
        }
        return ans;
    }

    private static byte[] readBlockBytes(FileChannel ch, BlockDirEntry e) throws IOException {
        ByteBuffer buf = ByteBuffer.allocate(e.compressedLen);
        readFully(ch, e.fileOffset, buf);
        return buf.array();
    }

    private static byte[] loadSelectedAnchorPages(
            FileChannel ch,
            int anchorLen,
            int numPages,
            boolean[] needed,
            long[] offsets,
            int[] byteLens,
            int[] baseLens,
            boolean fast2Bit,
            boolean leanMode) throws IOException {
        boolean any = false;
        for (boolean b : needed) {
            if (b) {
                any = true;
                break;
            }
        }
        if (!any) return null;
        byte[] anchorBases = new byte[anchorLen];
        populateSelectedAnchorPages(
                ch, anchorBases, numPages, needed, offsets, byteLens, baseLens, fast2Bit, leanMode);
        return anchorBases;
    }

    private static void populateSelectedAnchorPages(
            FileChannel ch,
            byte[] anchorBases,
            int numPages,
            boolean[] needed,
            long[] offsets,
            int[] byteLens,
            int[] baseLens,
            boolean fast2Bit,
            boolean leanMode) throws IOException {
        int pageOffset = 0;
        for (int p = 0; p < numPages; p++) {
            if (!needed[p]) {
                pageOffset += baseLens[p];
                continue;
            }
            ByteBuffer pBuf = ByteBuffer.allocate(byteLens[p]);
            readFully(ch, offsets[p], pBuf);
            byte[] dec = decodeAnchorPage(pBuf.array(), baseLens[p], fast2Bit, leanMode);
            System.arraycopy(dec, 0, anchorBases, pageOffset, baseLens[p]);
            pageOffset += baseLens[p];
        }
    }

    private static ByteOut buildFooterBytes(
            List<BlockDirEntry> dir, long totalRawBytes, int totalRecords, int crcVal, long footerOffset) {
        ByteOut body = new ByteOut(dir.size() * 12 + 32);
        body.writeVarInt(dir.size());
        body.writeVarLong(totalRawBytes);
        body.writeVarInt(totalRecords);
        body.writeIntLE(crcVal);

        long prevEnd = 5L;
        int prevFirstRec = 0;
        for (int i = 0; i < dir.size(); i++) {
            BlockDirEntry e = dir.get(i);
            long offDelta = (i == 0) ? e.fileOffset : (e.fileOffset - prevEnd);
            body.writeVarLong(offDelta);
            body.writeVarInt(e.compressedLen);
            body.writeVarInt(e.firstRecordIdx - prevFirstRec);
            body.writeVarInt(e.recordCount);
            body.writeVarLong(e.firstRecordBaseOffset);
            body.writeVarInt(e.blockSeqChars);
            prevEnd = e.fileOffset + e.compressedLen;
            prevFirstRec = e.firstRecordIdx;
        }

        ByteOut out = new ByteOut(body.size() + 20);
        VarIntCodec.writeCompressedAux(out, body);
        out.writeLongLE(footerOffset);
        out.writeIntLE(FOOTER_MAGIC);
        return out;
    }

    private static final class ParsedFooter {
        final long totalRawBytes;
        final int totalRecords;
        final int crc32;
        final List<BlockDirEntry> entries;

        ParsedFooter(long totalRawBytes, int totalRecords, int crc32, List<BlockDirEntry> entries) {
            this.totalRawBytes = totalRawBytes;
            this.totalRecords = totalRecords;
            this.crc32 = crc32;
            this.entries = entries;
        }
    }

    private static ParsedFooter parseFooterFromBytes(byte[] archive) {
        ByteIn trailer = new ByteIn(archive, archive.length - 12, 12);
        long footerOff = trailer.readLongLE();
        int fMagic = trailer.readIntLE();
        if (fMagic != FOOTER_MAGIC) {
            throw new IllegalStateException("Invalid .hgs footer magic");
        }
        ByteIn fIn = new ByteIn(archive, (int) footerOff, archive.length - 12 - (int) footerOff);
        return decodeFooterBody(new ByteIn(VarIntCodec.readCompressedAux(fIn)));
    }

    private static ParsedFooter readFooterFromChannel(FileChannel ch) throws IOException {
        long size = ch.size();
        ByteBuffer tr = ByteBuffer.allocate(12);
        readFully(ch, size - 12L, tr);
        ByteIn trIn = new ByteIn(tr.array());
        long footerOff = trIn.readLongLE();
        int fMagic = trIn.readIntLE();
        if (fMagic != FOOTER_MAGIC) {
            throw new IllegalStateException("Invalid .hgs footer magic");
        }
        int footerLen = (int) (size - 12L - footerOff);
        ByteBuffer fb = ByteBuffer.allocate(footerLen);
        readFully(ch, footerOff, fb);
        ByteIn fIn = new ByteIn(fb.array());
        return decodeFooterBody(new ByteIn(VarIntCodec.readCompressedAux(fIn)));
    }

    private static ParsedFooter decodeFooterBody(ByteIn body) {
        int numBlocks = body.readVarInt();
        long totalRawBytes = body.readVarLong();
        int totalRecords = body.readVarInt();
        int crc = body.readIntLE();
        List<BlockDirEntry> list = new ArrayList<>(numBlocks);

        long prevEnd = 0L;
        int prevFirstRec = 0;
        for (int i = 0; i < numBlocks; i++) {
            BlockDirEntry e = new BlockDirEntry();
            long offDelta = body.readVarLong();
            e.fileOffset = (i == 0) ? offDelta : (prevEnd + offDelta);
            e.compressedLen = body.readVarInt();
            e.firstRecordIdx = prevFirstRec + body.readVarInt();
            e.recordCount = body.readVarInt();
            e.firstRecordBaseOffset = body.readVarLong();
            e.blockSeqChars = body.readVarInt();
            prevEnd = e.fileOffset + e.compressedLen;
            prevFirstRec = e.firstRecordIdx;
            list.add(e);
        }
        return new ParsedFooter(totalRawBytes, totalRecords, crc, list);
    }

    private static void verifyMagic(ByteBuffer buf) {
        for (int i = 0; i < 4; i++) {
            if (buf.get() != MAGIC[i]) {
                throw new IllegalArgumentException("Invalid .hgs magic header");
            }
        }
    }

    private static void verifyIntegrity(byte[] raw, long expectedLen, int expectedCrc) {
        if (raw.length != expectedLen) {
            throw new IllegalStateException(
                    "Decoded length mismatch: expected " + expectedLen + ", got " + raw.length);
        }
        CRC32 crc = new CRC32();
        crc.update(raw, 0, raw.length);
        if ((int) crc.getValue() != expectedCrc) {
            throw new IllegalStateException("CRC-32 checksum mismatch on decompressed output!");
        }
    }

    private static void readFully(FileChannel ch, long pos, ByteBuffer dst) throws IOException {
        long cur = pos;
        while (dst.hasRemaining()) {
            int n = ch.read(dst, cur);
            if (n < 0) break;
            cur += n;
        }
    }

    /** Command-line interface for HyperGraSS. */
    public static void main(String[] args) throws Exception {
        if (args.length >= 1
                && ("--version".equals(args[0]) || "version".equals(args[0]) || "-V".equals(args[0]))) {
            System.out.println("HyperGraSS " + VERSION);
            return;
        }
        if (args.length < 2) {
            printUsage();
            return;
        }
        String cmd = args[0].toLowerCase();
        switch (cmd) {
            case "compress", "c" -> {
                CompressOptions opts = CompressOptions.defaults();
                int idx = 1;
                while (idx < args.length && args[idx].startsWith("--")) {
                    if ("--lean".equals(args[idx])) {
                        opts = CompressOptions.lean();
                    } else if ("--single-thread".equals(args[idx])) {
                        opts.parallel = false;
                    } else if ("--optimize".equals(args[idx])) {
                        opts.optimize = true;
                    } else if ("--pangenome".equals(args[idx])) {
                        // Scale-tuned profile: 16 MiB anchor bank, 8 MiB blocks.
                        opts.maxAnchorBases = 16 * 1024 * 1024;
                        opts.blockSizeBytes = 8 * 1024 * 1024;
                    } else if ("--anchor-single-page".equals(args[idx])) {
                        opts.anchorSinglePage = true;
                    } else if ("--bounded".equals(args[idx])) {
                        opts.bounded = true;
                    } else if ("--anchor-mib".equals(args[idx])) {
                        if (idx + 1 >= args.length) {
                            printUsage();
                            return;
                        }
                        int mib = Integer.parseInt(args[++idx]);
                        if (mib <= 0) {
                            System.err.println("--anchor-mib must be positive");
                            return;
                        }
                        opts.maxAnchorBases = mib * 1024 * 1024;
                    } else if ("--block-mib".equals(args[idx])) {
                        if (idx + 1 >= args.length) {
                            printUsage();
                            return;
                        }
                        int mib = Integer.parseInt(args[++idx]);
                        if (mib <= 0) {
                            System.err.println("--block-mib must be positive");
                            return;
                        }
                        opts.blockSizeBytes = mib * 1024 * 1024;
                    }
                    idx++;
                }
                if (idx + 1 >= args.length) {
                    printUsage();
                    return;
                }
                Path in = Path.of(args[idx]);
                Path out = Path.of(args[idx + 1]);
                long t0 = System.nanoTime();
                compressFile(in, out, opts);
                long t1 = System.nanoTime();
                long inBytes = Files.size(in);
                long outBytes = Files.size(out);
                double cr = outBytes > 0 ? (double) inBytes / outBytes : 0.0;
                System.out.printf(
                        "Compressed %s (%,d B) -> %s (%,d B) | CR = %.2fx | %.2f ms%n",
                        in.getFileName(), inBytes, out.getFileName(), outBytes, cr, (t1 - t0) / 1e6);
            }
            case "decompress", "d" -> {
                Path in = Path.of(args[1]);
                Path out = Path.of(args[2]);
                long t0 = System.nanoTime();
                decompressFile(in, out);
                long t1 = System.nanoTime();
                System.out.printf(
                        "Decompressed %s -> %s (%,d B) | CRC-32 OK | %.2f ms%n",
                        in.getFileName(), out.getFileName(), Files.size(out), (t1 - t0) / 1e6);
            }
            case "query", "q" -> {
                Path archive = Path.of(args[1]);
                int recIdx = Integer.parseInt(args[2]);
                int qStart = (args.length >= 4) ? Integer.parseInt(args[3]) : -1;
                int qEnd = (args.length >= 5) ? Integer.parseInt(args[4]) : -1;
                long t0 = System.nanoTime();
                byte[] slice = query(archive, recIdx, qStart, qEnd);
                long t1 = System.nanoTime();
                System.err.printf("[Query record=%d [%d..%d) -> %,d B in %.3f ms]%n",
                        recIdx, qStart, qEnd, slice.length, (t1 - t0) / 1e6);
                System.out.write(slice);
                if (slice.length > 0 && slice[slice.length - 1] != '\n') {
                    System.out.println();
                }
            }
            case "profile", "p" -> ArchiveByteAccounting.main(new String[] {args[1]});
            case "inspect", "i" -> {
                Path archive = Path.of(args[1]);
                ArchiveMetadata meta = inspect(archive);
                long compBytes = Files.size(archive);
                System.out.printf("Archive          : %s (%,d B)%n", archive.getFileName(), compBytes);
                System.out.printf("Raw Size         : %,d B (CR = %.2fx)%n",
                        meta.totalRawBytes, compBytes > 0 ? (double) meta.totalRawBytes / compBytes : 0.0);
                System.out.printf("Records / Blocks : %,d records across %,d block(s)%n",
                        meta.totalRecords, meta.numBlocks);
                System.out.printf("Mode             : compactSingle=%s, anchorBank=%s (%,d bp, %d pages), lean=%s%n",
                        meta.compactSingleBlock, meta.hasAnchorBank, meta.anchorTotalBases,
                        meta.numAnchorPages, meta.leanMode);
                System.out.printf("CRC-32           : 0x%08X%n", meta.crc32);
            }
            default -> printUsage();
        }
    }

    private static void printUsage() {
        System.out.println("HyperGraSS CLI (Java 21):");
        System.out.println("  compress [--lean] [--single-thread] [--optimize] [--bounded] <input.fa> <output.hgs>");
        System.out.println("    --optimize: bounded extra size search for small archives; costs more CPU");
        System.out.println("    --bounded: two-pass bounded-memory compression for large anchored archives");
        System.out.println("    --pangenome: scale-tuned 16 MiB anchor bank and 8 MiB blocks");
        System.out.println("    --anchor-mib N / --block-mib N: explicit tuning for large inputs");
        System.out.println("  decompress <input.hgs> <output.fa>");
        System.out.println("  query <input.hgs> <recordIdx> [startBase endBase]");
        System.out.println("  inspect <input.hgs>");
        System.out.println("  profile <input.hgs>  (exact encoded-byte accounting; reads archive into memory)");
    }
}
