package hypergrass;

import hypergrass.EimdlMosaicCodec.EncodedMosaic;
import hypergrass.EimdlMosaicCodec.MosaicConfig;
import hypergrass.EimdlMosaicCodec.MosaicInspectResult;
import hypergrass.EimdlMosaicCodec.VariantRecord;
import hypergrass.HyperGraSSFusionEngine.CompressionReport;
import hypergrass.HyperGraSSFusionEngine.FusionOptions;
import hypergrass.HyperGraSSFusionEngine.LoadedPanel;
import hypergrass.VarIntCodec.ByteIn;
import hypergrass.VarIntCodec.ByteOut;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.zip.CRC32;

/**
 * Unified Multi-Omic {@code .hgs3} (v3.0) Container Codec in Pure Java 25 LTS.
 *
 * <p>100% wire-compatible with native Rust/C23/C++23 {@code hgs3.rs} + {@code hgs_ans.c}:
 *
 * <ul>
 *   <li><b>Lossless FASTQ Mode:</b> Strips common instrument ID prefixes &amp; sorts LCP read IDs,
 *       packs CRLF/LF &amp; {@code +} comments, isolates ambiguity/soft-masked bases cleanly, models
 *       DNA via reference-conditioned Mosaic + parallel 4-stripe Byte-rANS, and couples Phred
 *       quality scores directly to Mosaic mismatch &amp; homopolymer bins at zero extra bits.
 *   <li><b>FASTA Mode:</b> Delegates to {@link HyperGraSSFusionEngine} upgraded with parallel
 *       4-stripe Byte-rANS DNA novelty ({@code HGS3}) + optional 6-frame ORF proteome digest.
 * </ul>
 */
public final class Hgs3OmicCodec {

    public enum OmicFormat {
        AUTO,
        FASTA,
        FASTQ
    }

    public static final class Hgs3Options {
        public OmicFormat format = OmicFormat.AUTO;
        public FusionOptions fusion = new FusionOptions();
        public boolean emitCdsDigest = false;
    }

    public record CdsSummary(
            long orfCount, int longestOrfAa, long totalOrfAa, int gcPermille, byte[] proteomeSha256) {}

    public record Hgs3Report(
            boolean isFastq,
            int mode,
            long rawBytes,
            long archiveBytes,
            int totalRecords,
            long dnaBases,
            long qualBytes,
            int metaCompBytes,
            int mosaicCompBytes,
            int noveltyCompBytes,
            int qualCompBytes,
            CdsSummary cds,
            double seconds) {
        public double bpb() {
            return rawBytes == 0 ? 0.0 : (archiveBytes * 8.0) / rawBytes;
        }

        public double ratio() {
            return archiveBytes == 0 ? 0.0 : ((double) rawBytes) / archiveBytes;
        }
    }

    static final class ParsedHgs3Fastq {
        int mode;
        long rawBytes;
        int expectedCrc;
        int totalRecords;
        int panelCount;
        String[] panelNames;
        byte[] panelSha;
        int totalDnaBases;
        int totalQualBytes;
        int readLenHint;
        byte[] metaComp;
        byte[] mosaicBytes;
        int novLen;
        byte[] ansNovComp;
        byte[] qualComp;
        byte[] cdsComp;
    }

    private Hgs3OmicCodec() {}

    public static boolean isHgs3Magic(byte[] b) {
        return b != null
                && b.length >= 4
                && b[0] == 'H'
                && b[1] == 'G'
                && b[2] == 'S'
                && b[3] == '3';
    }

    private static final byte[] CODON_TABLE =
            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
                    .getBytes(StandardCharsets.US_ASCII);

    public static CdsSummary scanCdsSummary(byte[] twoBit) {
        if (twoBit == null || twoBit.length == 0) {
            return null;
        }
        int n = twoBit.length;
        long gc = 0L;
        for (byte b : twoBit) {
            int v = b & 3;
            if (v == 1 || v == 2) gc++;
        }
        int gcPermille = (int) Math.min(1000L, (gc * 1000L) / n);

        MessageDigest md;
        try {
            md = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }

        long orfCount = 0L;
        int longest = 0;
        long totalAa = 0L;
        ByteOut pep = new ByteOut(256);

        for (int rev = 0; rev <= 1; rev++) {
            for (int frame = 0; frame < 3; frame++) {
                int pos = frame;
                boolean inOrf = false;
                pep.reset();
                while (pos + 3 <= n) {
                    int b0;
                    int b1;
                    int b2;
                    if (rev == 0) {
                        b0 = twoBit[pos] & 3;
                        b1 = twoBit[pos + 1] & 3;
                        b2 = twoBit[pos + 2] & 3;
                    } else {
                        int i = n - 1 - pos;
                        b0 = (twoBit[i] ^ 3) & 3;
                        b1 = (twoBit[i - 1] ^ 3) & 3;
                        b2 = (twoBit[i - 2] ^ 3) & 3;
                    }
                    byte aa = CODON_TABLE[(b0 << 4) | (b1 << 2) | b2];
                    if (!inOrf) {
                        if (aa == 'M') {
                            inOrf = true;
                            pep.reset();
                            pep.writeByte('M');
                        }
                    } else if (aa == '*') {
                        int plen = pep.size();
                        if (plen >= 30) {
                            orfCount++;
                            totalAa += plen;
                            if (plen > longest) longest = plen;
                            md.update(pep.internalBuffer(), 0, plen);
                            md.update((byte) '*');
                        }
                        inOrf = false;
                    } else {
                        pep.writeByte(aa);
                    }
                    pos += 3;
                }
            }
        }
        return new CdsSummary(orfCount, longest, totalAa, gcPermille, md.digest());
    }

    static byte[] serializeCdsSummary(CdsSummary s) {
        ByteOut w = new ByteOut(64);
        w.writeVarLong(s.orfCount);
        w.writeVarInt(s.longestOrfAa);
        w.writeVarLong(s.totalOrfAa);
        w.writeVarInt(s.gcPermille);
        w.writeBytes(s.proteomeSha256);
        return HyperGraSSFusionEngine.deflateBytes(w.toByteArray());
    }

    static CdsSummary deserializeCdsSummary(byte[] comp) {
        byte[] raw = HyperGraSSFusionEngine.inflateBytes(comp, 0, comp.length);
        ByteIn r = new ByteIn(raw);
        long orf = r.readVarLong();
        int longest = r.readVarInt();
        long tot = r.readVarLong();
        int gc = r.readVarInt();
        byte[] sha = r.readBytes(32);
        return new CdsSummary(orf, longest, tot, gc, sha);
    }

    static byte[] buildMismatchMask(int len, List<VariantRecord> variants) {
        byte[] mask = new byte[len];
        for (VariantRecord v : variants) {
            int start = (int) v.tgtPos();
            if (start < 0 || start >= len) continue;
            int span = switch (v.kind()) {
                case 0 -> 1;
                case 1, 3, 4 -> Math.max(1, v.altLen());
                default -> 0;
            };
            int end = Math.min(len, start + span);
            for (int i = start; i < end; i++) {
                mask[i] = 1;
            }
        }
        return mask;
    }

    public static OmicFormat detectInputFormat(Path path, OmicFormat explicit) throws IOException {
        if (explicit != OmicFormat.AUTO) return explicit;
        String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        if (name.endsWith(".fq") || name.endsWith(".fastq")) {
            return OmicFormat.FASTQ;
        }
        try (InputStream in = Files.newInputStream(path)) {
            byte[] probe = in.readNBytes(64);
            for (byte b : probe) {
                if (Character.isWhitespace(b & 0xFF)) continue;
                return (b == '@') ? OmicFormat.FASTQ : OmicFormat.FASTA;
            }
        }
        return OmicFormat.FASTA;
    }

    static ParsedHgs3Fastq parseFastqContainer(byte[] bytes, boolean decodePanel) {
        if (!isHgs3Magic(bytes) || bytes.length < 18) {
            throw new IllegalArgumentException("Not a valid HGS3 archive");
        }
        ByteIn r = new ByteIn(bytes);
        r.readBytes(4); // HGS3
        int flags = r.readByte() & 0xFF;
        int mode = flags & 3;
        if ((flags & HyperGraSSFusionEngine.FLAG3_FASTQ) == 0) {
            throw new IllegalArgumentException("HGS3 archive does not carry FASTQ payload");
        }
        boolean hasCds = (flags & HyperGraSSFusionEngine.FLAG3_HAS_CDS) != 0;
        ParsedHgs3Fastq p = new ParsedHgs3Fastq();
        p.mode = mode;
        p.rawBytes = r.readLongLE();
        p.expectedCrc = r.readIntLE();
        p.totalRecords = r.readVarInt();
        int numUnits = r.readVarInt();

        if (mode == HyperGraSSFusionEngine.MODE_REF) {
            p.panelCount = r.readVarInt();
            int cnLen = r.readVarInt();
            byte[] cn = r.readBytes(cnLen);
            if (decodePanel) {
                p.panelNames = HyperGraSSFusionEngine.decompressPanelNames(cn, p.panelCount);
            }
            p.panelSha = r.readBytes(32);
        }

        if (p.rawBytes == 0 || numUnits == 0 || p.totalRecords == 0) {
            p.readLenHint = 1;
            p.cdsComp = new byte[0];
            return p;
        }

        p.totalDnaBases = r.readVarInt();
        p.totalQualBytes = r.readVarInt();
        p.readLenHint = Math.max(1, r.readVarInt());

        p.metaComp = r.readBytes(r.readVarInt());
        p.mosaicBytes = r.readBytes(r.readVarInt());
        p.novLen = r.readVarInt();
        p.ansNovComp = r.readBytes(r.readVarInt());
        p.qualComp = r.readBytes(r.readVarInt());
        if (hasCds && r.hasRemaining()) {
            p.cdsComp = r.readBytes(r.readVarInt());
        } else {
            p.cdsComp = new byte[0];
        }
        return p;
    }

    private static int[] scanLine(byte[] raw, int pos) {
        int start = pos;
        int n = raw.length;
        while (pos < n) {
            if (raw[pos] == '\n') {
                return new int[] {start, pos, 0, pos + 1};
            }
            if (raw[pos] == '\r' && pos + 1 < n && raw[pos + 1] == '\n') {
                return new int[] {start, pos, 1, pos + 2};
            }
            pos++;
        }
        return new int[] {start, n, 2, n};
    }

    private static final byte[] EOL_LF = {'\n'};
    private static final byte[] EOL_CRLF = {'\r', '\n'};
    private static final byte[] EOL_NONE = {};

    private static byte[] eolBytes(int code) {
        return switch (code) {
            case 0 -> EOL_LF;
            case 1 -> EOL_CRLF;
            default -> EOL_NONE;
        };
    }

    private static int lcpLen(byte[] a, int aS, int aE, byte[] b, int bS, int bE) {
        int lim = Math.min(aE - aS, bE - bS);
        int i = 0;
        while (i < lim && a[aS + i] == b[bS + i]) i++;
        return i;
    }

    public static Hgs3Report compressOmics(Path input, Path output, Hgs3Options opts)
            throws IOException {
        OmicFormat fmt = detectInputFormat(input, opts.format);
        if (fmt == OmicFormat.FASTQ) {
            return compressFastqToHgs3(input, output, opts);
        }

        long t0 = System.nanoTime();
        CdsSummary cds = null;
        byte[] cdsComp = null;
        if (opts.emitCdsDigest) {
            byte[] raw = Files.readAllBytes(input);
            byte[] two = extractFastaTwoBitSimple(raw);
            cds = scanCdsSummary(two);
            if (cds != null) {
                cdsComp = serializeCdsSummary(cds);
            }
        }
        CompressionReport rep =
                HyperGraSSFusionEngine.compressFileToHgs3Fasta(input, output, opts.fusion, cdsComp);
        return new Hgs3Report(
                false,
                rep.mode(),
                rep.rawBytes(),
                rep.archiveBytes(),
                0,
                rep.targetDnaBases(),
                0L,
                0,
                0,
                0,
                0,
                cds,
                (System.nanoTime() - t0) / 1e9);
    }

    private static byte[] extractFastaTwoBitSimple(byte[] raw) {
        ByteOut two = new ByteOut(raw.length);
        boolean inHdr = false;
        boolean atSol = true;
        for (byte rawB : raw) {
            int b = rawB & 0xFF;
            if (inHdr) {
                if (b == '\n') {
                    inHdr = false;
                    atSol = true;
                }
                continue;
            }
            if (b == '>' && atSol) {
                inHdr = true;
                atSol = false;
                continue;
            }
            atSol = (b == '\n' || b == '\r');
            switch (b) {
                case 'A', 'a' -> two.writeByte(0);
                case 'C', 'c' -> two.writeByte(1);
                case 'G', 'g' -> two.writeByte(2);
                case 'T', 't', 'U', 'u' -> two.writeByte(3);
                default -> {}
            }
        }
        return two.toByteArray();
    }

    private static Hgs3Report compressFastqToHgs3(Path input, Path output, Hgs3Options opts)
            throws IOException {
        long t0 = System.nanoTime();
        byte[] raw = Files.readAllBytes(input);
        long rawBytes = raw.length;
        CRC32 crc = new CRC32();
        crc.update(raw, 0, raw.length);
        int rawCrc = (int) crc.getValue();
        int threads = Math.max(1, opts.fusion.threads);
        int mode =
                (opts.fusion.referencePanel != null)
                        ? HyperGraSSFusionEngine.MODE_REF
                        : HyperGraSSFusionEngine.MODE_SINGLE;

        if (raw.length == 0) {
            ByteOut out = new ByteOut(64);
            out.writeBytes("HGS3".getBytes(StandardCharsets.US_ASCII));
            out.writeByte((mode & 3) | HyperGraSSFusionEngine.FLAG3_FASTQ);
            out.writeLongLE(0L);
            out.writeIntLE(0);
            out.writeVarInt(0);
            out.writeVarInt(0);
            if (mode == HyperGraSSFusionEngine.MODE_REF) {
                LoadedPanel lp =
                        HyperGraSSFusionEngine.loadReferencePanelForCompress(
                                opts.fusion.referencePanel,
                                opts.fusion.excludePanelNames,
                                opts.fusion.explicitPanelNames,
                                opts.fusion.panelLimit);
                out.writeVarInt(lp.names().length);
                byte[] cn = HyperGraSSFusionEngine.compressPanelNames(lp.names());
                out.writeVarInt(cn.length);
                out.writeBytes(cn);
                out.writeBytes(lp.manifestSha256());
            }
            byte[] emptyBytes = out.toByteArray();
            Files.write(output, emptyBytes);
            return new Hgs3Report(
                    true, mode, 0L, emptyBytes.length, 0, 0L, 0L, 0, 0, 0, 0, null,
                    (System.nanoTime() - t0) / 1e9);
        }

        int nRecCap = Math.max(16, raw.length / 120);
        int[] idS = new int[nRecCap];
        int[] idE = new int[nRecCap];
        int[] plusS = new int[nRecCap];
        int[] plusE = new int[nRecCap];
        int[] seqS = new int[nRecCap];
        int[] seqE = new int[nRecCap];
        int[] qualS = new int[nRecCap];
        int[] qualE = new int[nRecCap];
        int[] eolFlat = new int[nRecCap * 4];
        int nRecs = 0;
        int pos = 0;
        int totalDna = 0;
        int totalQual = 0;

        while (pos < raw.length) {
            if (raw[pos] != '@') {
                throw new IllegalArgumentException(
                        "Invalid FASTQ start byte at offset " + pos);
            }
            pos++;
            if (nRecs == idS.length) {
                int nc = idS.length << 1;
                idS = Arrays.copyOf(idS, nc);
                idE = Arrays.copyOf(idE, nc);
                plusS = Arrays.copyOf(plusS, nc);
                plusE = Arrays.copyOf(plusE, nc);
                seqS = Arrays.copyOf(seqS, nc);
                seqE = Arrays.copyOf(seqE, nc);
                qualS = Arrays.copyOf(qualS, nc);
                qualE = Arrays.copyOf(qualE, nc);
                eolFlat = Arrays.copyOf(eolFlat, nc * 4);
            }
            int[] l1 = scanLine(raw, pos);
            pos = l1[3];
            int[] l2 = scanLine(raw, pos);
            pos = l2[3];
            if (pos >= raw.length || raw[pos] != '+') {
                throw new IllegalArgumentException("Expected '+' FASTQ separator line at " + pos);
            }
            pos++;
            int[] l3 = scanLine(raw, pos);
            pos = l3[3];
            int[] l4 = scanLine(raw, pos);
            pos = l4[3];

            idS[nRecs] = l1[0];
            idE[nRecs] = l1[1];
            seqS[nRecs] = l2[0];
            seqE[nRecs] = l2[1];
            plusS[nRecs] = l3[0];
            plusE[nRecs] = l3[1];
            qualS[nRecs] = l4[0];
            qualE[nRecs] = l4[1];
            eolFlat[nRecs * 4] = l1[2];
            eolFlat[nRecs * 4 + 1] = l2[2];
            eolFlat[nRecs * 4 + 2] = l3[2];
            eolFlat[nRecs * 4 + 3] = l4[2];
            totalDna += l2[1] - l2[0];
            totalQual += l4[1] - l4[0];
            nRecs++;
        }

        int firstEol = eolFlat[0];
        boolean uniformInt = firstEol < 2;
        for (int i = 0; i < nRecs && uniformInt; i++) {
            int o = i * 4;
            if (eolFlat[o] != firstEol
                    || eolFlat[o + 1] != firstEol
                    || eolFlat[o + 2] != firstEol
                    || (i + 1 < nRecs && eolFlat[o + 3] != firstEol)) {
                uniformInt = false;
            }
        }
        int eolMode = uniformInt ? firstEol : 2;
        int finalEol = eolFlat[(nRecs - 1) * 4 + 3];

        boolean allPlusEmpty = true;
        boolean allPlusEqId = true;
        for (int i = 0; i < nRecs; i++) {
            if (plusS[i] != plusE[i]) allPlusEmpty = false;
            if (!Arrays.equals(raw, plusS[i], plusE[i], raw, idS[i], idE[i])) {
                allPlusEqId = false;
            }
        }
        int plusMode = allPlusEmpty ? 0 : (allPlusEqId ? 1 : 2);

        int firstSeqLen = seqE[0] - seqS[0];
        boolean fixedSeqLen = true;
        boolean fixedQualLen = true;
        int[] sortedLens = new int[nRecs];
        for (int i = 0; i < nRecs; i++) {
            int sl = seqE[i] - seqS[i];
            int ql = qualE[i] - qualS[i];
            sortedLens[i] = sl;
            if (sl != firstSeqLen) fixedSeqLen = false;
            if (ql != sl) fixedQualLen = false;
        }
        Arrays.sort(sortedLens);
        int readLenHint = Math.max(1, sortedLens[nRecs / 2]);

        int cpLen = idE[0] - idS[0];
        for (int i = 1; i < nRecs && cpLen > 0; i++) {
            cpLen = lcpLen(raw, idS[0], idS[0] + cpLen, raw, idS[i], idE[i]);
        }

        byte[] dna2Bit = new byte[totalDna];
        byte[] qualsAll = new byte[totalQual];
        int[] excPos = new int[Math.min(totalDna, 1024)];
        byte[] excByte = new byte[excPos.length];
        int excCount = 0;
        int dCur = 0;
        int qCur = 0;
        byte[] acgtUp = {'A', 'C', 'G', 'T'};

        for (int r = 0; r < nRecs; r++) {
            for (int p = seqS[r]; p < seqE[r]; p++) {
                byte ch = raw[p];
                byte tb =
                        switch (ch) {
                            case 'A', 'a' -> 0;
                            case 'C', 'c' -> 1;
                            case 'G', 'g' -> 2;
                            case 'T', 't', 'U', 'u' -> 3;
                            default -> 0;
                        };
                dna2Bit[dCur] = tb;
                if (ch != acgtUp[tb]) {
                    if (excCount == excPos.length) {
                        int nc = excPos.length << 1;
                        excPos = Arrays.copyOf(excPos, nc);
                        excByte = Arrays.copyOf(excByte, nc);
                    }
                    excPos[excCount] = dCur;
                    excByte[excCount] = ch;
                    excCount++;
                }
                dCur++;
            }
            int ql = qualE[r] - qualS[r];
            System.arraycopy(raw, qualS[r], qualsAll, qCur, ql);
            qCur += ql;
        }

        ByteOut mw = new ByteOut(nRecs * 12 + 256);
        int uflags =
                (eolMode & 3)
                        | ((plusMode & 3) << 2)
                        | (fixedSeqLen ? 16 : 0)
                        | (fixedQualLen ? 32 : 0);
        mw.writeByte(uflags);
        mw.writeByte(finalEol);
        mw.writeVarInt(cpLen);
        mw.writeBytes(raw, idS[0], cpLen);
        if (eolMode == 2) {
            for (int i = 0; i < nRecs; i++) {
                for (int j = 0; j < 4; j++) mw.writeByte(eolFlat[i * 4 + j]);
            }
        }
        if (fixedSeqLen) {
            mw.writeVarInt(firstSeqLen);
        } else {
            for (int i = 0; i < nRecs; i++) mw.writeVarInt(seqE[i] - seqS[i]);
        }
        if (!fixedQualLen) {
            for (int i = 0; i < nRecs; i++) mw.writeVarInt(qualE[i] - qualS[i]);
        }
        int prevSufS = 0;
        int prevSufE = 0;
        for (int i = 0; i < nRecs; i++) {
            int sS = idS[i] + cpLen;
            int sE = idE[i];
            int l = lcpLen(raw, prevSufS, prevSufE, raw, sS, sE);
            mw.writeVarInt(l);
            mw.writeVarInt((sE - sS) - l);
            mw.writeBytes(raw, sS + l, (sE - sS) - l);
            prevSufS = sS;
            prevSufE = sE;
        }
        if (plusMode == 2) {
            for (int i = 0; i < nRecs; i++) {
                if (plusS[i] == plusE[i]) {
                    mw.writeByte(0);
                } else if (Arrays.equals(raw, plusS[i], plusE[i], raw, idS[i], idE[i])) {
                    mw.writeByte(1);
                } else {
                    mw.writeByte(2);
                    mw.writeVarInt(plusE[i] - plusS[i]);
                    mw.writeBytes(raw, plusS[i], plusE[i] - plusS[i]);
                }
            }
        }
        if (excCount == 0) {
            mw.writeByte(0);
        } else if (excCount * 3 <= totalDna) {
            mw.writeByte(1);
            mw.writeVarInt(excCount);
            int last = 0;
            for (int i = 0; i < excCount; i++) {
                mw.writeVarInt(excPos[i] - last);
                mw.writeByte(excByte[i]);
                last = excPos[i];
            }
        } else {
            mw.writeByte(2);
            byte[] dense = new byte[totalDna];
            for (int i = 0; i < excCount; i++) {
                int ep = excPos[i];
                dense[ep] = (byte) (excByte[i] ^ acgtUp[dna2Bit[ep]]);
            }
            mw.writeBytes(dense);
        }
        byte[] metaComp = HyperGraSSFusionEngine.deflateBytes(mw.toByteArray());

        LoadedPanel lp = null;
        Packed2BitSequence refBank = null;
        int[] refOffsets = null;
        int refCount = 0;
        if (mode == HyperGraSSFusionEngine.MODE_REF) {
            lp =
                    HyperGraSSFusionEngine.loadReferencePanelForCompress(
                            opts.fusion.referencePanel,
                            opts.fusion.excludePanelNames,
                            opts.fusion.explicitPanelNames,
                            opts.fusion.panelLimit);
            refCount = lp.seqs().length;
            refOffsets = new int[refCount + 1];
            refBank = Packed2BitSequence.concatenate(lp.seqs(), refOffsets);
        }

        Packed2BitSequence targetSeq = Packed2BitSequence.fromTwoBitBytes(dna2Bit, totalDna);
        EncodedMosaic em;
        if (lp != null && lp.seqs().length > 0) {
            int subShift = lp.totalBases() > 300_000_000L ? 2 : 0;
            boolean buildSec = lp.totalBases() < 250_000_000L;
            EimdlSyncmerIndex idx =
                    EimdlSyncmerIndex.buildWithSubsampling(
                            lp.seqs(), subShift, false, buildSec, threads);
            em =
                    EimdlMosaicCodec.encodeMosaic(
                            targetSeq, idx, refCount, -1, null, opts.fusion.mosaicConfig, threads);
        } else {
            em =
                    EimdlMosaicCodec.encodeMosaic(
                            targetSeq, null, 0, -1, null, opts.fusion.mosaicConfig, threads);
        }

        int novLen = em.noveltyLen();
        byte[] novExact = Arrays.copyOf(em.noveltyBases(), novLen);
        byte[] ansNovComp = HgsAnsCoder.ansNoveltyEncode(novExact, threads);

        MosaicInspectResult inspectRep =
                EimdlMosaicCodec.inspectMosaic(
                        em.mosaicBytes(), novExact, refBank, refOffsets, refCount, true, false);
        byte[] mm = buildMismatchMask(totalDna, inspectRep.variants());
        byte[] hb = HgsAnsCoder.computeHomoBins(dna2Bit);
        if (mm.length != totalQual) mm = Arrays.copyOf(mm, totalQual);
        if (hb.length != totalQual) hb = Arrays.copyOf(hb, totalQual);
        byte[] qualComp = HgsAnsCoder.ans4EncodeQual(qualsAll, mm, hb, readLenHint);

        CdsSummary cds = opts.emitCdsDigest ? scanCdsSummary(dna2Bit) : null;

        ByteOut out = new ByteOut(128 + metaComp.length + em.mosaicBytes().length);
        int hflags = (mode & 3) | HyperGraSSFusionEngine.FLAG3_FASTQ;
        if (cds != null) hflags |= HyperGraSSFusionEngine.FLAG3_HAS_CDS;
        out.writeBytes("HGS3".getBytes(StandardCharsets.US_ASCII));
        out.writeByte(hflags);
        out.writeLongLE(rawBytes);
        out.writeIntLE(rawCrc);
        out.writeVarInt(nRecs);
        out.writeVarInt(1);

        if (lp != null) {
            out.writeVarInt(lp.names().length);
            byte[] cn = HyperGraSSFusionEngine.compressPanelNames(lp.names());
            out.writeVarInt(cn.length);
            out.writeBytes(cn);
            out.writeBytes(lp.manifestSha256());
        }

        out.writeVarInt(totalDna);
        out.writeVarInt(totalQual);
        out.writeVarInt(readLenHint);

        out.writeVarInt(metaComp.length);
        out.writeBytes(metaComp);
        out.writeVarInt(em.mosaicBytes().length);
        out.writeBytes(em.mosaicBytes());
        out.writeVarInt(novLen);
        out.writeVarInt(ansNovComp.length);
        out.writeBytes(ansNovComp);
        out.writeVarInt(qualComp.length);
        out.writeBytes(qualComp);
        if (cds != null) {
            byte[] cb = serializeCdsSummary(cds);
            out.writeVarInt(cb.length);
            out.writeBytes(cb);
        }

        byte[] finalBytes = out.toByteArray();
        Files.write(output, finalBytes);
        return new Hgs3Report(
                true,
                mode,
                rawBytes,
                finalBytes.length,
                nRecs,
                totalDna,
                totalQual,
                metaComp.length,
                em.mosaicBytes().length,
                ansNovComp.length,
                qualComp.length,
                cds,
                (System.nanoTime() - t0) / 1e9);
    }

    public static double decompressOmics(Path archive, Path output, Path panel, int threads)
            throws IOException {
        long t0 = System.nanoTime();
        byte[] bytes = Files.readAllBytes(archive);
        if (bytes.length >= 4 && HyperGraSSFusionEngine.isHgs2Magic(bytes)) {
            HyperGraSSFusionEngine.decompressFile(archive, output, panel);
            return (System.nanoTime() - t0) / 1e9;
        }
        if (!isHgs3Magic(bytes) || bytes.length < 5) {
            throw new IllegalArgumentException("Archive is neither HGS2 nor HGS3");
        }
        if ((bytes[4] & HyperGraSSFusionEngine.FLAG3_FASTQ) != 0) {
            decompressFastqHgs3(bytes, output, panel, threads);
            return (System.nanoTime() - t0) / 1e9;
        }
        HyperGraSSFusionEngine.decompressFile(archive, output, panel);
        return (System.nanoTime() - t0) / 1e9;
    }

    private static void decompressFastqHgs3(byte[] bytes, Path output, Path panel, int threads)
            throws IOException {
        ParsedHgs3Fastq p = parseFastqContainer(bytes, true);
        Packed2BitSequence refBank = null;
        int[] refOffsets = null;
        int refCount = 0;
        if (p.mode == HyperGraSSFusionEngine.MODE_REF && p.panelCount > 0) {
            Packed2BitSequence[] seqs =
                    HyperGraSSFusionEngine.verifyAndLoadPanelForDecode(
                            panel, p.panelNames, p.panelSha);
            refCount = seqs.length;
            refOffsets = new int[refCount + 1];
            refBank = Packed2BitSequence.concatenate(seqs, refOffsets);
        }
        if (p.rawBytes == 0 || p.totalRecords == 0) {
            if (p.expectedCrc != 0) {
                throw new IllegalStateException("CRC mismatch on empty HGS3 archive");
            }
            Files.write(output, new byte[0]);
            return;
        }

        byte[] novDna = HgsAnsCoder.ansNoveltyDecode(p.ansNovComp, p.novLen, Math.max(1, threads));
        byte[] dna2Bit =
                EimdlMosaicCodec.decodeMosaicToTwoBit(
                        p.mosaicBytes, novDna, refBank, refOffsets, refCount);
        if (dna2Bit.length != p.totalDnaBases) {
            throw new IllegalStateException("HGS3 FASTQ DNA length mismatch");
        }

        MosaicInspectResult inspectRep =
                EimdlMosaicCodec.inspectMosaic(
                        p.mosaicBytes, novDna, refBank, refOffsets, refCount, true, false);
        byte[] mm = buildMismatchMask(dna2Bit.length, inspectRep.variants());
        byte[] hb = HgsAnsCoder.computeHomoBins(dna2Bit);
        if (mm.length != p.totalQualBytes) mm = Arrays.copyOf(mm, p.totalQualBytes);
        if (hb.length != p.totalQualBytes) hb = Arrays.copyOf(hb, p.totalQualBytes);
        byte[] quals = HgsAnsCoder.ans4DecodeQual(p.qualComp, mm, hb, p.totalQualBytes, p.readLenHint);

        byte[] metaRaw = HyperGraSSFusionEngine.inflateBytes(p.metaComp, 0, p.metaComp.length);
        ByteIn mr = new ByteIn(metaRaw);
        int uflags = mr.readByte() & 0xFF;
        int eolMode = uflags & 3;
        int plusMode = (uflags >>> 2) & 3;
        boolean fixedSeqLen = (uflags & 16) != 0;
        boolean fixedQualLen = (uflags & 32) != 0;
        int finalEol = mr.readByte() & 0xFF;
        int cpLen = mr.readVarInt();
        byte[] commonPrefix = mr.readBytes(cpLen);

        int nRecs = p.totalRecords;
        int[] perRecEols = eolMode == 2 ? new int[nRecs * 4] : null;
        if (eolMode == 2) {
            for (int i = 0; i < nRecs * 4; i++) perRecEols[i] = mr.readByte() & 0xFF;
        }

        int[] seqLens = new int[nRecs];
        if (fixedSeqLen) {
            Arrays.fill(seqLens, mr.readVarInt());
        } else {
            for (int i = 0; i < nRecs; i++) seqLens[i] = mr.readVarInt();
        }
        int[] qualLens = new int[nRecs];
        if (fixedQualLen) {
            System.arraycopy(seqLens, 0, qualLens, 0, nRecs);
        } else {
            for (int i = 0; i < nRecs; i++) qualLens[i] = mr.readVarInt();
        }

        byte[][] ids = new byte[nRecs][];
        byte[] prevSuf = new byte[0];
        for (int i = 0; i < nRecs; i++) {
            int l = mr.readVarInt();
            int rem = mr.readVarInt();
            byte[] suf = new byte[l + rem];
            System.arraycopy(prevSuf, 0, suf, 0, l);
            byte[] tail = mr.readBytes(rem);
            System.arraycopy(tail, 0, suf, l, rem);
            byte[] full = new byte[cpLen + suf.length];
            System.arraycopy(commonPrefix, 0, full, 0, cpLen);
            System.arraycopy(suf, 0, full, cpLen, suf.length);
            prevSuf = suf;
            ids[i] = full;
        }

        byte[][] plusLines = plusMode == 2 ? new byte[nRecs][] : null;
        if (plusMode == 2) {
            for (int i = 0; i < nRecs; i++) {
                int tag = mr.readByte() & 0xFF;
                plusLines[i] =
                        switch (tag) {
                            case 0 -> new byte[0];
                            case 1 -> ids[i];
                            case 2 -> mr.readBytes(mr.readVarInt());
                            default -> throw new IllegalStateException("Invalid FASTQ plus tag " + tag);
                        };
            }
        }

        byte[] acgtUp = {'A', 'C', 'G', 'T'};
        byte[] seqAscii = new byte[dna2Bit.length];
        for (int i = 0; i < dna2Bit.length; i++) seqAscii[i] = acgtUp[dna2Bit[i] & 3];
        int maskMode = mr.readByte() & 0xFF;
        switch (maskMode) {
            case 0 -> {}
            case 1 -> {
                int cnt = mr.readVarInt();
                int cur = 0;
                for (int i = 0; i < cnt; i++) {
                    cur += mr.readVarInt();
                    seqAscii[cur] = (byte) mr.readByte();
                }
            }
            case 2 -> {
                byte[] dense = mr.readBytes(dna2Bit.length);
                for (int i = 0; i < dna2Bit.length; i++) {
                    if (dense[i] != 0) seqAscii[i] ^= dense[i];
                }
            }
            default -> throw new IllegalStateException("Invalid ambiguity mask mode " + maskMode);
        }

        ByteArrayOutputStream bos = new ByteArrayOutputStream((int) Math.min(p.rawBytes, 16_777_216));
        CRC32 crc = new CRC32();
        int dOff = 0;
        int qOff = 0;
        for (int i = 0; i < nRecs; i++) {
            int e0 = eolMode == 2 ? perRecEols[i * 4] : eolMode;
            int e1 = eolMode == 2 ? perRecEols[i * 4 + 1] : eolMode;
            int e2 = eolMode == 2 ? perRecEols[i * 4 + 2] : eolMode;
            int e3 = eolMode == 2 ? perRecEols[i * 4 + 3] : (i + 1 == nRecs ? finalEol : eolMode);

            bos.write('@');
            crc.update('@');
            bos.write(ids[i], 0, ids[i].length);
            crc.update(ids[i], 0, ids[i].length);
            byte[] eb0 = eolBytes(e0);
            bos.write(eb0, 0, eb0.length);
            crc.update(eb0, 0, eb0.length);

            int sl = seqLens[i];
            bos.write(seqAscii, dOff, sl);
            crc.update(seqAscii, dOff, sl);
            dOff += sl;
            byte[] eb1 = eolBytes(e1);
            bos.write(eb1, 0, eb1.length);
            crc.update(eb1, 0, eb1.length);

            bos.write('+');
            crc.update('+');
            if (plusMode == 1) {
                bos.write(ids[i], 0, ids[i].length);
                crc.update(ids[i], 0, ids[i].length);
            } else if (plusMode == 2) {
                bos.write(plusLines[i], 0, plusLines[i].length);
                crc.update(plusLines[i], 0, plusLines[i].length);
            }
            byte[] eb2 = eolBytes(e2);
            bos.write(eb2, 0, eb2.length);
            crc.update(eb2, 0, eb2.length);

            int ql = qualLens[i];
            bos.write(quals, qOff, ql);
            crc.update(quals, qOff, ql);
            qOff += ql;
            byte[] eb3 = eolBytes(e3);
            bos.write(eb3, 0, eb3.length);
            crc.update(eb3, 0, eb3.length);
        }

        byte[] recOut = bos.toByteArray();
        if (recOut.length != p.rawBytes) {
            throw new IllegalStateException("HGS3 FASTQ length mismatch");
        }
        if ((int) crc.getValue() != p.expectedCrc) {
            throw new IllegalStateException("HGS3 FASTQ CRC-32 mismatch");
        }
        Files.write(output, recOut);
    }
}
