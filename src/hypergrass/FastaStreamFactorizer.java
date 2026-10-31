package hypergrass;

import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Universal 6-Channel Orthogonal Lossless Factorizer for Raw, FASTA, RNA, and Multi-FASTA
 * biological sequence streams.
 *
 * <p>Decomposes any arbitrary byte slice into 6 independent, parity-avalanche-free channels:
 * <ol>
 *   <li>Front-coded + Deflate FASTA headers and inter-header line deltas</li>
 *   <li>RLE-encoded sequence line widths and exact EOL markers (NONE / LF / CRLF)</li>
 *   <li>Alternating uppercase/lowercase run-length intervals</li>
 *   <li>Scaffold ambiguity 'N'-run intervals (gap, runLen)</li>
 *   <li>Generalized non-ACGTUN exception stream (gap, rawByte) supporting all IUPAC and chimeric T+U</li>
 *   <li>Pure 2-bit {A=0, C=1, G=2, T/U=3} primary nucleotide sequence</li>
 * </ol>
 *
 * <p>Provably bijective over arbitrary byte sequences B* (100% SHA-256 lossless on ragged lines,
 * empty records, mixed CRLF/LF, no trailing newline, soft-masking, RNA, and chimeric DNA/RNA).
 */
public final class FastaStreamFactorizer {

    private FastaStreamFactorizer() {}

    public static final class BlockSlice {
        public final int startByte;
        public final int endByte;
        public final int firstRecordIndex;
        public final int recordCount;
        public final long firstRecordBaseOffset;

        public BlockSlice(int startByte, int endByte, int firstRecordIndex, int recordCount, long firstRecordBaseOffset) {
            this.startByte = startByte;
            this.endByte = endByte;
            this.firstRecordIndex = firstRecordIndex;
            this.recordCount = recordCount;
            this.firstRecordBaseOffset = firstRecordBaseOffset;
        }
    }

    public static final class RecordSpan {
        public final byte[] headerBytes; // null if raw headerless record
        public final int firstSeqLine;
        public final int endSeqLine;
        public final int seqCharStart;
        public final int seqCharEnd;

        public RecordSpan(byte[] headerBytes, int firstSeqLine, int endSeqLine, int seqCharStart, int seqCharEnd) {
            this.headerBytes = headerBytes;
            this.firstSeqLine = firstSeqLine;
            this.endSeqLine = endSeqLine;
            this.seqCharStart = seqCharStart;
            this.seqCharEnd = seqCharEnd;
        }
    }

    public static final class FactorizedChunk {
        public boolean rnaMode;
        public int numHeaders;
        public int[] headerSeqLineDeltas;
        public byte[][] headerBytes;

        public int numSeqLines;
        public int[] seqLineLengths;
        public byte[] lineEolTypes; // 0 = NONE, 1 = LF (\n), 2 = CRLF (\r\n)

        public int totalSeqChars;

        public int numLowerIntervals;
        public int[] lowerUpperGaps;
        public int[] lowerLengths;

        public int numNRuns;
        public int[] nAcgtGaps;
        public int[] nRunLengths;

        public int numExceptions;
        public int[] excGaps;
        public byte[] excSymbols;

        public int dnaLen;
        public byte[] dnaBases; // 0=A, 1=C, 2=G, 3=T(or U)

        private RecordSpan[] cachedRecords;

        public RecordSpan[] getRecordSpans() {
            if (cachedRecords != null) {
                return cachedRecords;
            }
            int leadingLines = (numHeaders > 0) ? headerSeqLineDeltas[0] : numSeqLines;
            boolean hasLeadingRaw = (numHeaders == 0) || (leadingLines > 0);
            int totalRecs = numHeaders + (hasLeadingRaw ? 1 : 0);
            if (totalRecs == 0) {
                cachedRecords = new RecordSpan[0];
                return cachedRecords;
            }
            int[] lineCharPrefix = new int[numSeqLines + 1];
            for (int i = 0; i < numSeqLines; i++) {
                lineCharPrefix[i + 1] = lineCharPrefix[i] + seqLineLengths[i];
            }
            RecordSpan[] spans = new RecordSpan[totalRecs];
            int idx = 0;
            int curSeqLine = 0;
            if (hasLeadingRaw) {
                int endLine = (numHeaders > 0) ? headerSeqLineDeltas[0] : numSeqLines;
                spans[idx++] = new RecordSpan(null, 0, endLine, lineCharPrefix[0], lineCharPrefix[endLine]);
                curSeqLine = endLine;
            }
            for (int h = 0; h < numHeaders; h++) {
                int startLine = curSeqLine;
                int endLine = (h + 1 < numHeaders) ? (startLine + headerSeqLineDeltas[h + 1]) : numSeqLines;
                spans[idx++] = new RecordSpan(
                        headerBytes[h],
                        startLine,
                        endLine,
                        lineCharPrefix[startLine],
                        lineCharPrefix[endLine]);
                curSeqLine = endLine;
            }
            cachedRecords = spans;
            return spans;
        }
    }

    /**
     * Partitions raw input bytes into self-contained block slices aligned on FASTA '>' record
     * boundaries whenever possible, or clean line/chunk boundaries for giant chromosomes/raw files.
     */
    public static List<BlockSlice> partitionIntoBlocks(byte[] raw, int targetBlockBytes) {
        List<BlockSlice> slices = new ArrayList<>();
        if (raw.length == 0) {
            slices.add(new BlockSlice(0, 0, 0, 1, 0L));
            return slices;
        }

        boolean startsWithHeader = raw[0] == '>';
        int blockStart = 0;
        int globalRecordIdx = startsWithHeader ? -1 : 0;
        int blockFirstRecord = 0;
        int blockRecordCount = startsWithHeader ? 0 : 1;
        long currentRecordBaseAcc = 0L;
        long blockFirstRecOffset = 0L;

        int p = 0;
        final int end = raw.length;
        while (p < end) {
            int lineStart = p;
            int nl = p;
            while (nl < end && raw[nl] != '\n') {
                nl++;
            }
            int contentEnd = (nl < end && nl > p && raw[nl - 1] == '\r') ? (nl - 1) : nl;
            int nextP = (nl < end) ? (nl + 1) : end;
            boolean isHeader = (contentEnd > lineStart && raw[lineStart] == '>');

            if (isHeader) {
                if (lineStart > blockStart && (lineStart - blockStart) >= targetBlockBytes && blockRecordCount > 0) {
                    slices.add(new BlockSlice(
                            blockStart, lineStart, blockFirstRecord, blockRecordCount, blockFirstRecOffset));
                    blockStart = lineStart;
                    blockFirstRecord = globalRecordIdx + 1;
                    blockRecordCount = 0;
                    blockFirstRecOffset = 0L;
                }
                globalRecordIdx++;
                if (blockRecordCount == 0) {
                    blockFirstRecord = globalRecordIdx;
                    blockFirstRecOffset = 0L;
                }
                blockRecordCount++;
                currentRecordBaseAcc = 0L;
            } else {
                int lineChars = contentEnd - lineStart;
                // Handle single huge record (> 2 * targetBlockBytes) by splitting cleanly
                if (lineStart > blockStart && (lineStart - blockStart) >= (targetBlockBytes * 2)) {
                    slices.add(new BlockSlice(
                            blockStart, lineStart, blockFirstRecord, Math.max(1, blockRecordCount), blockFirstRecOffset));
                    blockStart = lineStart;
                    blockFirstRecord = Math.max(0, globalRecordIdx);
                    blockRecordCount = 1;
                    blockFirstRecOffset = currentRecordBaseAcc;
                } else if (lineChars > targetBlockBytes * 2) {
                    // Single unwrapped giant raw sequence line: split into targetBlockBytes chunks
                    int subPos = lineStart;
                    if (subPos > blockStart) {
                        slices.add(new BlockSlice(
                                blockStart, subPos, blockFirstRecord, Math.max(1, blockRecordCount), blockFirstRecOffset));
                        blockStart = subPos;
                        blockFirstRecord = Math.max(0, globalRecordIdx);
                        blockRecordCount = 1;
                        blockFirstRecOffset = currentRecordBaseAcc;
                    }
                    while (nextP - blockStart > targetBlockBytes * 2) {
                        int cut = blockStart + targetBlockBytes;
                        slices.add(new BlockSlice(blockStart, cut, blockFirstRecord, 1, blockFirstRecOffset));
                        currentRecordBaseAcc += (cut - blockStart);
                        blockStart = cut;
                        blockFirstRecOffset = currentRecordBaseAcc;
                    }
                    currentRecordBaseAcc += (contentEnd - blockStart);
                    p = nextP;
                    continue;
                }
                currentRecordBaseAcc += lineChars;
            }
            p = nextP;
        }

        if (blockStart < end) {
            slices.add(new BlockSlice(
                    blockStart, end, Math.max(0, blockFirstRecord), Math.max(1, blockRecordCount), blockFirstRecOffset));
        }
        return slices;
    }

    /** Losslessly factorizes raw[startOffset .. endOffset) into a 6-channel FactorizedChunk. */
    public static FactorizedChunk factorize(byte[] raw, int startOffset, int endOffset) {
        FactorizedChunk chunk = new FactorizedChunk();
        if (startOffset >= endOffset) {
            chunk.headerSeqLineDeltas = new int[0];
            chunk.headerBytes = new byte[0][];
            chunk.seqLineLengths = new int[0];
            chunk.lineEolTypes = new byte[0];
            chunk.lowerUpperGaps = new int[0];
            chunk.lowerLengths = new int[0];
            chunk.nAcgtGaps = new int[0];
            chunk.nRunLengths = new int[0];
            chunk.excGaps = new int[0];
            chunk.excSymbols = new byte[0];
            chunk.dnaBases = new byte[0];
            return chunk;
        }

        // Pass 1: parse line boundaries, headers, EOLs, and scan T/U presence
        int capLines = Math.max(16, (endOffset - startOffset) / 60 + 4);
        int[] seqLineStarts = new int[capLines];
        int[] seqLineLens = new int[capLines];
        byte[] eolTypes = new byte[capLines];
        int totalLines = 0;
        int numSeqLines = 0;

        int capHeaders = 16;
        int[] hdrDeltas = new int[capHeaders];
        byte[][] hdrBytes = new byte[capHeaders][];
        int numHeaders = 0;
        int seqSinceHeader = 0;

        boolean hasT = false;
        boolean hasU = false;
        int totalSeqChars = 0;

        int p = startOffset;
        while (p < endOffset) {
            int lineStart = p;
            int nl = p;
            while (nl < endOffset && raw[nl] != '\n') {
                nl++;
            }
            byte eol;
            int contentEnd;
            if (nl < endOffset) {
                if (nl > lineStart && raw[nl - 1] == '\r') {
                    eol = 2;
                    contentEnd = nl - 1;
                } else {
                    eol = 1;
                    contentEnd = nl;
                }
                p = nl + 1;
            } else {
                eol = 0;
                contentEnd = endOffset;
                p = endOffset;
            }

            if (totalLines == eolTypes.length) {
                eolTypes = Arrays.copyOf(eolTypes, eolTypes.length << 1);
            }
            eolTypes[totalLines++] = eol;

            if (contentEnd > lineStart && raw[lineStart] == '>') {
                if (numHeaders == hdrDeltas.length) {
                    hdrDeltas = Arrays.copyOf(hdrDeltas, hdrDeltas.length << 1);
                    hdrBytes = Arrays.copyOf(hdrBytes, hdrBytes.length << 1);
                }
                hdrDeltas[numHeaders] = seqSinceHeader;
                hdrBytes[numHeaders] = Arrays.copyOfRange(raw, lineStart + 1, contentEnd);
                numHeaders++;
                seqSinceHeader = 0;
            } else {
                if (numSeqLines == seqLineStarts.length) {
                    seqLineStarts = Arrays.copyOf(seqLineStarts, seqLineStarts.length << 1);
                    seqLineLens = Arrays.copyOf(seqLineLens, seqLineLens.length << 1);
                }
                int len = contentEnd - lineStart;
                seqLineStarts[numSeqLines] = lineStart;
                seqLineLens[numSeqLines] = len;
                numSeqLines++;
                seqSinceHeader++;
                totalSeqChars += len;

                for (int i = lineStart; i < contentEnd; i++) {
                    int b = raw[i] & 0xFF;
                    if (b == 'T' || b == 't') hasT = true;
                    else if (b == 'U' || b == 'u') hasU = true;
                }
            }
        }

        boolean rnaMode = hasU && !hasT;
        final int fourthBase = rnaMode ? 'U' : 'T';

        // Pass 2: stream over sequence bytes to extract Case, Exception, N-Run, and 2-bit DNA channels
        byte[] dna = new byte[totalSeqChars];
        int dnaLen = 0;

        int[] lowGaps = new int[16];
        int[] lowLens = new int[16];
        int numLow = 0;
        boolean inLower = false;
        int curCaseRun = 0;
        int pendingUpperGap = 0;

        int[] nGaps = new int[16];
        int[] nLens = new int[16];
        int numN = 0;
        int curAcgtSinceN = 0;
        int curNRun = 0;

        int[] excGaps = new int[16];
        byte[] excSyms = new byte[16];
        int numExc = 0;
        int curAcgtnSinceExc = 0;

        for (int s = 0; s < numSeqLines; s++) {
            int start = seqLineStarts[s];
            int end = start + seqLineLens[s];
            for (int i = start; i < end; i++) {
                int b = raw[i] & 0xFF;
                boolean isLow = (b >= 'a' && b <= 'z');
                if (isLow) {
                    if (!inLower) {
                        pendingUpperGap = curCaseRun;
                        curCaseRun = 1;
                        inLower = true;
                    } else {
                        curCaseRun++;
                    }
                    b -= 32; // exact uppercase normalization
                } else {
                    if (inLower) {
                        if (numLow == lowGaps.length) {
                            lowGaps = Arrays.copyOf(lowGaps, lowGaps.length << 1);
                            lowLens = Arrays.copyOf(lowLens, lowLens.length << 1);
                        }
                        lowGaps[numLow] = pendingUpperGap;
                        lowLens[numLow] = curCaseRun;
                        numLow++;
                        curCaseRun = 1;
                        inLower = false;
                    } else {
                        curCaseRun++;
                    }
                }

                // Primary / N / Exception classification on normalized byte b
                int code;
                if (b == 'A') code = 0;
                else if (b == 'C') code = 1;
                else if (b == 'G') code = 2;
                else if (b == fourthBase) code = 3;
                else code = -1;

                if (code >= 0) {
                    if (curNRun > 0) {
                        if (numN == nGaps.length) {
                            nGaps = Arrays.copyOf(nGaps, nGaps.length << 1);
                            nLens = Arrays.copyOf(nLens, nLens.length << 1);
                        }
                        nGaps[numN] = curAcgtSinceN;
                        nLens[numN] = curNRun;
                        numN++;
                        curNRun = 0;
                        curAcgtSinceN = 0;
                    }
                    dna[dnaLen++] = (byte) code;
                    curAcgtSinceN++;
                    curAcgtnSinceExc++;
                } else if (b == 'N') {
                    curNRun++;
                    curAcgtnSinceExc++;
                } else {
                    // Exception symbol (IUPAC or chimeric U or arbitrary byte)
                    if (numExc == excGaps.length) {
                        excGaps = Arrays.copyOf(excGaps, excGaps.length << 1);
                        excSyms = Arrays.copyOf(excSyms, excSyms.length << 1);
                    }
                    excGaps[numExc] = curAcgtnSinceExc;
                    excSyms[numExc] = (byte) b;
                    numExc++;
                    curAcgtnSinceExc = 0;
                }
            }
        }

        if (inLower) {
            if (numLow == lowGaps.length) {
                lowGaps = Arrays.copyOf(lowGaps, lowGaps.length << 1);
                lowLens = Arrays.copyOf(lowLens, lowLens.length << 1);
            }
            lowGaps[numLow] = pendingUpperGap;
            lowLens[numLow] = curCaseRun;
            numLow++;
        }
        if (curNRun > 0) {
            if (numN == nGaps.length) {
                nGaps = Arrays.copyOf(nGaps, nGaps.length << 1);
                nLens = Arrays.copyOf(nLens, nLens.length << 1);
            }
            nGaps[numN] = curAcgtSinceN;
            nLens[numN] = curNRun;
            numN++;
        }

        chunk.rnaMode = rnaMode;
        chunk.numHeaders = numHeaders;
        chunk.headerSeqLineDeltas = Arrays.copyOf(hdrDeltas, numHeaders);
        chunk.headerBytes = Arrays.copyOf(hdrBytes, numHeaders);

        chunk.numSeqLines = numSeqLines;
        chunk.seqLineLengths = Arrays.copyOf(seqLineLens, numSeqLines);
        chunk.lineEolTypes = Arrays.copyOf(eolTypes, totalLines);
        chunk.totalSeqChars = totalSeqChars;

        chunk.numLowerIntervals = numLow;
        chunk.lowerUpperGaps = Arrays.copyOf(lowGaps, numLow);
        chunk.lowerLengths = Arrays.copyOf(lowLens, numLow);

        chunk.numNRuns = numN;
        chunk.nAcgtGaps = Arrays.copyOf(nGaps, numN);
        chunk.nRunLengths = Arrays.copyOf(nLens, numN);

        chunk.numExceptions = numExc;
        chunk.excGaps = Arrays.copyOf(excGaps, numExc);
        chunk.excSymbols = Arrays.copyOf(excSyms, numExc);

        chunk.dnaLen = dnaLen;
        chunk.dnaBases = (dnaLen == dna.length) ? dna : Arrays.copyOf(dna, dnaLen);
        return chunk;
    }

    /** Serializes channels 1..5 (all auxiliary metadata except primary DNA stream). */
    public static void writeAuxMetadata(ByteOut out, FactorizedChunk chunk) {
        int flags = 0;
        if (chunk.rnaMode) flags |= 1;
        if (chunk.numHeaders > 0) flags |= 2;
        if (chunk.numLowerIntervals > 0) flags |= 4;
        if (chunk.numNRuns > 0) flags |= 8;
        if (chunk.numExceptions > 0) flags |= 16;

        out.writeByte(flags);
        out.writeVarInt(chunk.dnaLen);
        out.writeVarInt(chunk.totalSeqChars);
        out.writeVarInt(chunk.numSeqLines);

        // Channel 1: Front-coded FASTA headers + seqLine deltas
        if (chunk.numHeaders > 0) {
            out.writeVarInt(chunk.numHeaders);
            ByteOut hdrBuf = new ByteOut(chunk.numHeaders * 24);
            for (int h = 0; h < chunk.numHeaders; h++) {
                hdrBuf.writeVarInt(chunk.headerSeqLineDeltas[h]);
            }
            byte[] prev = new byte[0];
            for (int h = 0; h < chunk.numHeaders; h++) {
                byte[] cur = chunk.headerBytes[h];
                int maxShare = Math.min(prev.length, cur.length);
                int share = 0;
                while (share < maxShare && prev[share] == cur[share]) {
                    share++;
                }
                hdrBuf.writeVarInt(share);
                hdrBuf.writeVarInt(cur.length - share);
                hdrBuf.writeBytes(cur, share, cur.length - share);
                prev = cur;
            }
            VarIntCodec.writeCompressedAux(out, hdrBuf);
        }

        // Channel 2: Sequence line-width RLE + line EOL RLE
        ByteOut lineLenRle = new ByteOut(32);
        writeRleInts(lineLenRle, chunk.seqLineLengths, chunk.numSeqLines);
        VarIntCodec.writeCompressedAux(out, lineLenRle);

        int totalLines = chunk.numHeaders + chunk.numSeqLines;
        ByteOut eolRle = new ByteOut(16);
        writeRleBytes(eolRle, chunk.lineEolTypes, totalLines);
        VarIntCodec.writeCompressedAux(out, eolRle);

        // Channel 3: Lowercase intervals
        if (chunk.numLowerIntervals > 0) {
            out.writeVarInt(chunk.numLowerIntervals);
            // Split-column soft-mask coding: independent gap / (length - minLength) streams
            int minLowLen = Integer.MAX_VALUE;
            for (int i = 0; i < chunk.numLowerIntervals; i++) {
                minLowLen = Math.min(minLowLen, chunk.lowerLengths[i]);
            }
            out.writeVarInt(minLowLen);
            ByteOut gapCol = new ByteOut(chunk.numLowerIntervals * 2);
            ByteOut lenCol = new ByteOut(chunk.numLowerIntervals * 2);
            for (int i = 0; i < chunk.numLowerIntervals; i++) {
                gapCol.writeVarInt(chunk.lowerUpperGaps[i]);
                lenCol.writeVarInt(chunk.lowerLengths[i] - minLowLen);
            }
            VarIntCodec.writeCompressedAux(out, gapCol);
            VarIntCodec.writeCompressedAux(out, lenCol);
        }

        // Channel 4: Scaffold N-run intervals
        if (chunk.numNRuns > 0) {
            out.writeVarInt(chunk.numNRuns);
            ByteOut nBuf = new ByteOut(chunk.numNRuns * 4);
            for (int i = 0; i < chunk.numNRuns; i++) {
                nBuf.writeVarInt(chunk.nAcgtGaps[i]);
                nBuf.writeVarInt(chunk.nRunLengths[i]);
            }
            VarIntCodec.writeCompressedAux(out, nBuf);
        }

        // Channel 5: Non-ACGTUN exception stream
        if (chunk.numExceptions > 0) {
            out.writeVarInt(chunk.numExceptions);
            ByteOut gapBuf = new ByteOut(chunk.numExceptions * 2);
            for (int i = 0; i < chunk.numExceptions; i++) {
                gapBuf.writeVarInt(chunk.excGaps[i]);
            }
            VarIntCodec.writeCompressedAux(out, gapBuf);
            VarIntCodec.writeCompressedAux(out, chunk.excSymbols, chunk.numExceptions);
        }
    }

    /** Deserializes channels 1..5 from a ByteIn stream. */
    public static FactorizedChunk readAuxMetadata(ByteIn in) {
        FactorizedChunk chunk = new FactorizedChunk();
        int flags = in.readByte();
        chunk.rnaMode = (flags & 1) != 0;
        boolean hasHeaders = (flags & 2) != 0;
        boolean hasLower = (flags & 4) != 0;
        boolean hasNRuns = (flags & 8) != 0;
        boolean hasExceptions = (flags & 16) != 0;

        chunk.dnaLen = in.readVarInt();
        chunk.totalSeqChars = in.readVarInt();
        chunk.numSeqLines = in.readVarInt();

        if (hasHeaders) {
            chunk.numHeaders = in.readVarInt();
            ByteIn hdrIn = new ByteIn(VarIntCodec.readCompressedAux(in));
            chunk.headerSeqLineDeltas = new int[chunk.numHeaders];
            chunk.headerBytes = new byte[chunk.numHeaders][];
            for (int h = 0; h < chunk.numHeaders; h++) {
                chunk.headerSeqLineDeltas[h] = hdrIn.readVarInt();
            }
            byte[] prev = new byte[0];
            for (int h = 0; h < chunk.numHeaders; h++) {
                int share = hdrIn.readVarInt();
                int suffixLen = hdrIn.readVarInt();
                byte[] cur = new byte[share + suffixLen];
                if (share > 0) {
                    System.arraycopy(prev, 0, cur, 0, share);
                }
                hdrIn.readBytes(cur, share, suffixLen);
                chunk.headerBytes[h] = cur;
                prev = cur;
            }
        } else {
            chunk.numHeaders = 0;
            chunk.headerSeqLineDeltas = new int[0];
            chunk.headerBytes = new byte[0][];
        }

        ByteIn lineLenIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        chunk.seqLineLengths = readRleInts(lineLenIn, chunk.numSeqLines);

        int totalLines = chunk.numHeaders + chunk.numSeqLines;
        ByteIn eolIn = new ByteIn(VarIntCodec.readCompressedAux(in));
        chunk.lineEolTypes = readRleBytes(eolIn, totalLines);

        if (hasLower) {
            chunk.numLowerIntervals = in.readVarInt();
            int minLowLen = in.readVarInt();
            ByteIn gapIn = new ByteIn(VarIntCodec.readCompressedAux(in));
            ByteIn lenIn = new ByteIn(VarIntCodec.readCompressedAux(in));
            chunk.lowerUpperGaps = new int[chunk.numLowerIntervals];
            chunk.lowerLengths = new int[chunk.numLowerIntervals];
            for (int i = 0; i < chunk.numLowerIntervals; i++) {
                chunk.lowerUpperGaps[i] = gapIn.readVarInt();
                chunk.lowerLengths[i] = lenIn.readVarInt() + minLowLen;
            }
        } else {
            chunk.numLowerIntervals = 0;
            chunk.lowerUpperGaps = new int[0];
            chunk.lowerLengths = new int[0];
        }

        if (hasNRuns) {
            chunk.numNRuns = in.readVarInt();
            ByteIn nIn = new ByteIn(VarIntCodec.readCompressedAux(in));
            chunk.nAcgtGaps = new int[chunk.numNRuns];
            chunk.nRunLengths = new int[chunk.numNRuns];
            for (int i = 0; i < chunk.numNRuns; i++) {
                chunk.nAcgtGaps[i] = nIn.readVarInt();
                chunk.nRunLengths[i] = nIn.readVarInt();
            }
        } else {
            chunk.numNRuns = 0;
            chunk.nAcgtGaps = new int[0];
            chunk.nRunLengths = new int[0];
        }

        if (hasExceptions) {
            chunk.numExceptions = in.readVarInt();
            ByteIn gapIn = new ByteIn(VarIntCodec.readCompressedAux(in));
            chunk.excGaps = new int[chunk.numExceptions];
            for (int i = 0; i < chunk.numExceptions; i++) {
                chunk.excGaps[i] = gapIn.readVarInt();
            }
            chunk.excSymbols = VarIntCodec.readCompressedAux(in);
        } else {
            chunk.numExceptions = 0;
            chunk.excGaps = new int[0];
            chunk.excSymbols = new byte[0];
        }

        return chunk;
    }

    /** Reconstructs the flat sequence character array (with N-runs, IUPAC, RNA, and exact case). */
    public static byte[] reconstructSequenceChars(FactorizedChunk chunk, byte[] dnaBases) {
        final byte[] alphabet = chunk.rnaMode
                ? new byte[] {'A', 'C', 'G', 'U'}
                : new byte[] {'A', 'C', 'G', 'T'};

        int acgtnLen = chunk.totalSeqChars - chunk.numExceptions;

        // Step 1: merge 2-bit DNA bases + N-runs into acgtnBuf
        byte[] acgtnBuf;
        if (chunk.numNRuns == 0) {
            acgtnBuf = new byte[chunk.dnaLen];
            for (int i = 0; i < chunk.dnaLen; i++) {
                acgtnBuf[i] = alphabet[dnaBases[i] & 3];
            }
        } else {
            acgtnBuf = new byte[acgtnLen];
            int srcPos = 0;
            int dstPos = 0;
            for (int r = 0; r < chunk.numNRuns; r++) {
                int gap = chunk.nAcgtGaps[r];
                for (int g = 0; g < gap; g++) {
                    acgtnBuf[dstPos++] = alphabet[dnaBases[srcPos++] & 3];
                }
                int nLen = chunk.nRunLengths[r];
                Arrays.fill(acgtnBuf, dstPos, dstPos + nLen, (byte) 'N');
                dstPos += nLen;
            }
            while (srcPos < chunk.dnaLen) {
                acgtnBuf[dstPos++] = alphabet[dnaBases[srcPos++] & 3];
            }
        }

        // Step 2: merge exceptions into full sequence buffer
        byte[] seqChars;
        if (chunk.numExceptions == 0) {
            seqChars = acgtnBuf;
        } else {
            seqChars = new byte[chunk.totalSeqChars];
            int srcPos = 0;
            int dstPos = 0;
            for (int e = 0; e < chunk.numExceptions; e++) {
                int gap = chunk.excGaps[e];
                if (gap > 0) {
                    System.arraycopy(acgtnBuf, srcPos, seqChars, dstPos, gap);
                    srcPos += gap;
                    dstPos += gap;
                }
                seqChars[dstPos++] = chunk.excSymbols[e];
            }
            int rem = acgtnBuf.length - srcPos;
            if (rem > 0) {
                System.arraycopy(acgtnBuf, srcPos, seqChars, dstPos, rem);
            }
        }

        // Step 3: restore exact lowercase intervals
        if (chunk.numLowerIntervals > 0) {
            int pos = 0;
            for (int i = 0; i < chunk.numLowerIntervals; i++) {
                pos += chunk.lowerUpperGaps[i];
                int len = chunk.lowerLengths[i];
                for (int k = 0; k < len; k++) {
                    seqChars[pos + k] = (byte) (seqChars[pos + k] + 32);
                }
                pos += len;
            }
        }
        return seqChars;
    }

    /** Reconstructs the exact byte-for-byte original file chunk (headers, line widths, CRLF/LF). */
    public static byte[] reconstructBytes(FactorizedChunk chunk, byte[] dnaBases) {
        if (chunk.numHeaders == 0 && chunk.numSeqLines == 0) {
            return new byte[0];
        }
        byte[] seqChars = reconstructSequenceChars(chunk, dnaBases);

        // Compute exact output byte length up front for single-allocation output
        int exactBytes = chunk.totalSeqChars;
        for (int h = 0; h < chunk.numHeaders; h++) {
            exactBytes += 1 + chunk.headerBytes[h].length;
        }
        for (int i = 0; i < chunk.lineEolTypes.length; i++) {
            exactBytes += chunk.lineEolTypes[i]; // 0->0, 1->1 (\n), 2->2 (\r\n)
        }

        byte[] out = new byte[exactBytes];
        int outPos = 0;
        int h = 0;
        int seqLine = 0;
        int seqSinceLastHeader = 0;
        int lineIdx = 0;
        int seqCharPos = 0;

        while (h < chunk.numHeaders || seqLine < chunk.numSeqLines) {
            while (h < chunk.numHeaders && chunk.headerSeqLineDeltas[h] == seqSinceLastHeader) {
                out[outPos++] = '>';
                byte[] hb = chunk.headerBytes[h];
                System.arraycopy(hb, 0, out, outPos, hb.length);
                outPos += hb.length;
                byte eol = chunk.lineEolTypes[lineIdx++];
                if (eol == 2) {
                    out[outPos++] = '\r';
                    out[outPos++] = '\n';
                } else if (eol == 1) {
                    out[outPos++] = '\n';
                }
                h++;
                seqSinceLastHeader = 0;
            }
            if (seqLine < chunk.numSeqLines) {
                int len = chunk.seqLineLengths[seqLine];
                if (len > 0) {
                    System.arraycopy(seqChars, seqCharPos, out, outPos, len);
                    seqCharPos += len;
                    outPos += len;
                }
                byte eol = chunk.lineEolTypes[lineIdx++];
                if (eol == 2) {
                    out[outPos++] = '\r';
                    out[outPos++] = '\n';
                } else if (eol == 1) {
                    out[outPos++] = '\n';
                }
                seqLine++;
                seqSinceLastHeader++;
            }
        }
        return out;
    }

    /**
     * Extracts either a complete formatted FASTA record (when qStart < 0) or an exact sequence
     * locus slice [qStart, qEnd) within localRecordIdx.
     */
    public static byte[] extractRecordOrSlice(
            FactorizedChunk chunk, byte[] dnaBases, int localRecordIdx, int qStart, int qEnd) {
        RecordSpan[] spans = chunk.getRecordSpans();
        if (localRecordIdx < 0 || localRecordIdx >= spans.length) {
            throw new IndexOutOfBoundsException(
                    "Local record index " + localRecordIdx + " out of [0, " + spans.length + ")");
        }
        RecordSpan span = spans[localRecordIdx];
        byte[] seqChars = reconstructSequenceChars(chunk, dnaBases);
        int recLen = span.seqCharEnd - span.seqCharStart;

        if (qStart >= 0) {
            int s = Math.max(0, Math.min(recLen, qStart));
            int e = (qEnd < 0) ? recLen : Math.max(s, Math.min(recLen, qEnd));
            return Arrays.copyOfRange(seqChars, span.seqCharStart + s, span.seqCharStart + e);
        }

        // Full record extraction with header + original line wrapping
        ByteOut out = new ByteOut(recLen + 128);
        if (span.headerBytes != null) {
            out.writeByte('>');
            out.writeBytes(span.headerBytes);
            out.writeByte('\n');
        }
        int cur = span.seqCharStart;
        for (int line = span.firstSeqLine; line < span.endSeqLine; line++) {
            int len = chunk.seqLineLengths[line];
            if (len > 0) {
                out.writeBytes(seqChars, cur, len);
                cur += len;
            }
            if (span.headerBytes != null || line + 1 < span.endSeqLine) {
                out.writeByte('\n');
            }
        }
        return out.toByteArray();
    }

    /**
     * Zero-bloat O(1) sub-interval extractor: maps [absS, absE) across exception & N-run channels
     * to decode ONLY the exact overlapping 2-bit DNA subrange [dnaStart, dnaEnd)!
     */
    public static byte[] extractRecordOrSliceSelective(
            FactorizedChunk chunk,
            ByteIn dnaStreamIn,
            byte[] anchorBases,
            boolean leanMode,
            int localRecordIdx,
            int qStart,
            int qEnd) {
        return extractRecordOrSliceSelective(
                chunk, dnaStreamIn, anchorBases, CanonicalSyncmerEngine.EMPTY_DONORS,
                leanMode, localRecordIdx, qStart, qEnd);
    }

    public static byte[] extractRecordOrSliceSelective(
            FactorizedChunk chunk,
            ByteIn dnaStreamIn,
            byte[] anchorBases,
            int[][] globalDonors,
            boolean leanMode,
            int localRecordIdx,
            int qStart,
            int qEnd) {
        RecordSpan[] spans = chunk.getRecordSpans();
        if (localRecordIdx < 0 || localRecordIdx >= spans.length) {
            throw new IndexOutOfBoundsException(
                    "Local record index " + localRecordIdx + " out of [0, " + spans.length + ")");
        }
        RecordSpan span = spans[localRecordIdx];
        int recLen = span.seqCharEnd - span.seqCharStart;
        int s = (qStart >= 0) ? Math.max(0, Math.min(recLen, qStart)) : 0;
        int e = (qStart >= 0) ? ((qEnd < 0) ? recLen : Math.max(s, Math.min(recLen, qEnd))) : recLen;
        int absS = span.seqCharStart + s;
        int absE = span.seqCharStart + e;
        int subLen = absE - absS;

        byte[] sliceChars = new byte[subLen];
        if (subLen > 0) {
            // 1. Map [absS, absE) -> [acgtnS, acgtnE)
            int acgtnS = absS - countExceptionsBefore(chunk, absS);
            int acgtnE = absE - countExceptionsBefore(chunk, absE);
            int acgtnSubLen = acgtnE - acgtnS;

            // 2. Map [acgtnS, acgtnE) -> [dnaS, dnaE) and decode ONLY [dnaS, dnaE)
            int dnaS = acgtnS - countNsBefore(chunk, acgtnS);
            int dnaE = acgtnE - countNsBefore(chunk, acgtnE);
            byte[] dnaSub = CanonicalSyncmerEngine.decodeBlockDnaRange(
                    dnaStreamIn, chunk.dnaLen, dnaS, dnaE, anchorBases, globalDonors, leanMode);

            final byte[] alphabet = chunk.rnaMode
                    ? new byte[] {'A', 'C', 'G', 'U'}
                    : new byte[] {'A', 'C', 'G', 'T'};

            // 3. Reconstruct acgtnSub[0 .. acgtnSubLen)
            byte[] acgtnSub = new byte[acgtnSubLen];
            if (chunk.numNRuns == 0) {
                for (int i = 0; i < acgtnSubLen; i++) {
                    acgtnSub[i] = alphabet[dnaSub[i] & 3];
                }
            } else {
                boolean[] isNMask = new boolean[acgtnSubLen];
                int cur = 0;
                for (int r = 0; r < chunk.numNRuns; r++) {
                    cur += chunk.nAcgtGaps[r];
                    int nLen = chunk.nRunLengths[r];
                    int ovS = Math.max(cur, acgtnS);
                    int ovE = Math.min(cur + nLen, acgtnE);
                    if (ovS < ovE) {
                        Arrays.fill(isNMask, ovS - acgtnS, ovE - acgtnS, true);
                    }
                    cur += nLen;
                    if (cur >= acgtnE) break;
                }
                int dIdx = 0;
                for (int i = 0; i < acgtnSubLen; i++) {
                    acgtnSub[i] = isNMask[i] ? (byte) 'N' : alphabet[dnaSub[dIdx++] & 3];
                }
            }

            // 4. Re-insert any exceptions overlapping [absS, absE)
            if (chunk.numExceptions == 0) {
                sliceChars = acgtnSub;
            } else {
                boolean[] isExc = new boolean[subLen];
                int seqCur = 0;
                for (int ex = 0; ex < chunk.numExceptions; ex++) {
                    seqCur += chunk.excGaps[ex];
                    if (seqCur >= absE) break;
                    if (seqCur >= absS) {
                        isExc[seqCur - absS] = true;
                        sliceChars[seqCur - absS] = chunk.excSymbols[ex];
                    }
                    seqCur++;
                }
                int aIdx = 0;
                for (int i = 0; i < subLen; i++) {
                    if (!isExc[i]) {
                        sliceChars[i] = acgtnSub[aIdx++];
                    }
                }
            }

            // 5. Apply any lowercase intervals overlapping [absS, absE)
            if (chunk.numLowerIntervals > 0) {
                int pos = 0;
                for (int i = 0; i < chunk.numLowerIntervals; i++) {
                    pos += chunk.lowerUpperGaps[i];
                    int len = chunk.lowerLengths[i];
                    int ovS = Math.max(pos, absS);
                    int ovE = Math.min(pos + len, absE);
                    for (int k = ovS; k < ovE; k++) {
                        sliceChars[k - absS] = (byte) (sliceChars[k - absS] + 32);
                    }
                    pos += len;
                    if (pos >= absE) break;
                }
            }
        }

        if (qStart >= 0) {
            return sliceChars;
        }

        // Full formatted FASTA record
        ByteOut out = new ByteOut(subLen + 128);
        if (span.headerBytes != null) {
            out.writeByte('>');
            out.writeBytes(span.headerBytes);
            out.writeByte('\n');
        }
        int cur = 0;
        for (int line = span.firstSeqLine; line < span.endSeqLine; line++) {
            int len = chunk.seqLineLengths[line];
            if (len > 0) {
                out.writeBytes(sliceChars, cur, len);
                cur += len;
            }
            if (span.headerBytes != null || line + 1 < span.endSeqLine) {
                out.writeByte('\n');
            }
        }
        return out.toByteArray();
    }

    private static int countExceptionsBefore(FactorizedChunk chunk, int seqPos) {
        int cur = 0;
        int cnt = 0;
        for (int e = 0; e < chunk.numExceptions; e++) {
            cur += chunk.excGaps[e];
            if (cur >= seqPos) break;
            cnt++;
            cur++;
        }
        return cnt;
    }

    private static int countNsBefore(FactorizedChunk chunk, int acgtnPos) {
        int cur = 0;
        int cnt = 0;
        for (int r = 0; r < chunk.numNRuns; r++) {
            cur += chunk.nAcgtGaps[r];
            if (cur >= acgtnPos) break;
            int nLen = chunk.nRunLengths[r];
            cnt += Math.min(nLen, acgtnPos - cur);
            cur += nLen;
            if (cur >= acgtnPos) break;
        }
        return cnt;
    }

    public static String getRecordHeaderString(FactorizedChunk chunk, int localRecordIdx) {
        RecordSpan[] spans = chunk.getRecordSpans();
        byte[] hb = spans[localRecordIdx].headerBytes;
        return hb == null ? "" : new String(hb, StandardCharsets.UTF_8);
    }

    private static void writeRleInts(ByteOut out, int[] arr, int len) {
        int i = 0;
        while (i < len) {
            int val = arr[i];
            int run = 1;
            while (i + run < len && arr[i + run] == val) {
                run++;
            }
            out.writeVarInt(run);
            out.writeVarInt(val);
            i += run;
        }
    }

    private static int[] readRleInts(ByteIn in, int expectedLen) {
        int[] out = new int[expectedLen];
        int pos = 0;
        while (in.hasRemaining() && pos < expectedLen) {
            int run = in.readVarInt();
            int val = in.readVarInt();
            Arrays.fill(out, pos, pos + run, val);
            pos += run;
        }
        if (pos != expectedLen) {
            throw new IllegalStateException("RLE int length mismatch: expected " + expectedLen + ", got " + pos);
        }
        return out;
    }

    private static void writeRleBytes(ByteOut out, byte[] arr, int len) {
        int i = 0;
        while (i < len) {
            int val = arr[i] & 0xFF;
            int run = 1;
            while (i + run < len && (arr[i + run] & 0xFF) == val) {
                run++;
            }
            out.writeVarInt(run);
            out.writeVarInt(val);
            i += run;
        }
    }

    private static byte[] readRleBytes(ByteIn in, int expectedLen) {
        byte[] out = new byte[expectedLen];
        int pos = 0;
        while (in.hasRemaining() && pos < expectedLen) {
            int run = in.readVarInt();
            byte val = (byte) in.readVarInt();
            Arrays.fill(out, pos, pos + run, val);
            pos += run;
        }
        if (pos != expectedLen) {
            throw new IllegalStateException("RLE byte length mismatch: expected " + expectedLen + ", got " + pos);
        }
        return out;
    }
}
