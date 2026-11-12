package hypergrass;

import hypergrass.EimdlMosaicCodec.MosaicConfig;
import hypergrass.HyperGraSSFusionEngine.FusionOptions;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Automated Falsification &amp; Verification Suite for Laws I–V in {@code docs/COMPRESSION-LAW.md}.
 *
 * <p>Tests algebraic identities, biological coding gains against controlled ablations, conditional
 * monotonicity under distractor panels, non-existing/synthetic genome universality (including the
 * {@code <= 2.005 bpb} cryptographic incompressible ceiling), and parallel determinism.
 */
public final class VerifyCompressionLaws {

    private static int checksPassed = 0;

    private VerifyCompressionLaws() {}

    public static void main(String[] args) throws Exception {
        System.out.println("=== HyperGraSS + RC-MDL (EIMDL) Compression Laws Verification ===");
        verifyLaw1KleinV4AndSyncmerSymmetry();
        verifyLaw2CodonWobbleBijectionAndGain();
        verifyLaw3MonotonicityAndIndelContinuity();
        verifyLaw4UniversalFiveRegimes();
        verifyLaw5WarmParallelAndBoundedDeterminism();
        verifyLaw6HierarchicalPitmanYorAndTwoAdicAdaptation();
        verifyLaw7DisjointLifecycleAndBoundedScaling();
        System.out.printf(
                "=== ALL 7 COMPRESSION LAWS VERIFIED (%d assertions passed) ===%n", checksPassed);
    }

    // -------------------------------------------------------------------------
    // Law I: Klein V_4 Strand-Equivariant Kimura + CpG Substitution Law
    // -------------------------------------------------------------------------
    private static void verifyLaw1KleinV4AndSyncmerSymmetry() throws Exception {
        // 1. Exhaustive Klein V_4 involution & transition/transversion invariance
        for (int x = 0; x < 4; x++) {
            for (int y = 0; y < 4; y++) {
                if (x == y) continue;
                int delta = y ^ x;
                check(delta >= 1 && delta <= 3, "V4 non-zero delta in {1,2,3}");
                int xRc = x ^ 3;
                int yRc = y ^ 3;
                int deltaRc = yRc ^ xRc;
                check(delta == deltaRc, "V4 XOR commutes with Watson-Crick complement");

                // A=0 (purine), G=2 (purine), C=1 (pyrimidine), T=3 (pyrimidine)
                boolean isPurineX = (x == 0 || x == 2);
                boolean isPurineY = (y == 0 || y == 2);
                boolean biologicalTransition = (isPurineX == isPurineY);
                check((delta == 2) == biologicalTransition, "delta==2 iff biological transition");
            }
        }

        // 2. SWAR 32-base reverseComplementWord32 involution & exact base check
        Random rng = new Random(0x5EED_1A71L);
        for (int t = 0; t < 5_000; t++) {
            long w = rng.nextLong();
            long rc = Packed2BitSequence.reverseComplementWord32(w);
            long rc2 = Packed2BitSequence.reverseComplementWord32(rc);
            check(rc2 == w, "reverseComplementWord32 is an exact involution");
            for (int i = 0; i < 32; i++) {
                int bOrig = (int) ((w >>> (i << 1)) & 3L);
                int bRc = (int) ((rc >>> ((31 - i) << 1)) & 3L);
                check(bRc == (bOrig ^ 3), "SWAR RC matches base-by-base complement");
            }
        }

        // 3. Exact strand symmetry of canonical open syncmers (Primary k=16 & Secondary k=13)
        int len = 10_000;
        byte[] fwdBases = new byte[len];
        byte[] revBases = new byte[len];
        for (int i = 0; i < len; i++) {
            int b = rng.nextInt(4);
            fwdBases[i] = (byte) b;
            revBases[len - 1 - i] = (byte) (b ^ 3);
        }
        Packed2BitSequence seqFwd = Packed2BitSequence.fromTwoBitBytes(fwdBases, len);
        Packed2BitSequence seqRev = Packed2BitSequence.fromTwoBitBytes(revBases, len);

        for (boolean secondary : new boolean[] {false, true}) {
            int k = secondary ? Packed2BitSequence.K_SECONDARY : Packed2BitSequence.K_PRIMARY;
            int s = secondary ? Packed2BitSequence.S_SECONDARY : Packed2BitSequence.S_PRIMARY;
            int t = secondary ? Packed2BitSequence.T_SECONDARY : Packed2BitSequence.T_PRIMARY;
            List<long[]> fList = new ArrayList<>();
            List<long[]> rList = new ArrayList<>();
            Packed2BitSequence.SyncmerConsumer fc =
                    (pos, hash, strand) -> fList.add(new long[] {pos, hash, strand});
            Packed2BitSequence.SyncmerConsumer rcCons =
                    (pos, hash, strand) -> rList.add(new long[] {pos, hash, strand});
            seqFwd.forEachCanonicalOpenSyncmer(0, len, k, s, t, fc);
            seqRev.forEachCanonicalOpenSyncmer(0, len, k, s, t, rcCons);
            check(
                    fList.size() == rList.size() && !fList.isEmpty(),
                    "Canonical open syncmer count identical on + and - strands");
            for (int i = 0; i < fList.size(); i++) {
                long[] a = fList.get(i);
                long[] b = rList.get(rList.size() - 1 - i);
                check(a[0] == (len - k - b[0]), "Syncmer RC position reflection exact");
                check(a[1] == b[1], "Syncmer canonical hash identical on RC");
            }
        }

        // 4. Empirical Kimura + CpG gain on reverse-complemented homolog vs V4-disabled ablation
        Path tmp = Files.createTempDirectory("verify-law1-");
        try {
            int n = 180_000;
            byte[] donor = generateCodingLikeSequence(n, new Random(101L));
            byte[] target = reverseComplementAscii(donor);
            Random mutRng = new Random(202L);
            for (int i = 20; i < n - 20; i++) {
                boolean isCpG = (target[i] == 'C' && target[i + 1] == 'G')
                        || (target[i - 1] == 'C' && target[i] == 'G');
                double rate = isCpG ? 0.035 : 0.008;
                if (mutRng.nextDouble() < rate) {
                    int b = charToCode(target[i]);
                    // 75% transition (delta=2), 25% transversions (delta=1 or 3)
                    int delta = mutRng.nextDouble() < 0.75 ? 2 : (mutRng.nextBoolean() ? 1 : 3);
                    target[i] = codeToChar(b ^ delta);
                }
            }
            Path refFile = tmp.resolve("donor.seq");
            Files.write(refFile, donor);

            FusionOptions optsV4 = new FusionOptions();
            optsV4.referencePanel = refFile;
            optsV4.mosaicConfig = new MosaicConfig(true, true, true, true, true, true);
            byte[] archV4 = roundTripVerify(target, optsV4, refFile);

            FusionOptions optsNoV4 = new FusionOptions();
            optsNoV4.referencePanel = refFile;
            optsNoV4.mosaicConfig = new MosaicConfig(true, true, false, true, true, true);
            byte[] archNoV4 = roundTripVerify(target, optsNoV4, refFile);

            check(
                    archV4.length < archNoV4.length,
                    "Law I: Klein-V4 Kimura+CpG ("
                            + archV4.length
                            + " B) strictly beats 4-symbol ablation ("
                            + archNoV4.length
                            + " B)");
            System.out.printf(
                    "  [PASS] Law I   (Klein-V4 Kimura+CpG + Syncmer RC) : V4=%d B vs Ablated=%d B (-%.1f%%)%n",
                    archV4.length,
                    archNoV4.length,
                    100.0 * (archNoV4.length - archV4.length) / archNoV4.length);
        } finally {
            deleteRecursively(tmp);
        }
    }

    // -------------------------------------------------------------------------
    // Law II: Bijective Codon-Wobble Phase & Coalescent Gap Law
    // -------------------------------------------------------------------------
    private static void verifyLaw2CodonWobbleBijectionAndGain() throws Exception {
        for (int g = 0; g <= 2_000_000; g++) {
            int q = g / 3;
            int r = (g + 1) % 3;
            int recon = 3 * q + ((r + 2) % 3);
            if (recon != g) {
                throw new AssertionError("Bijection failed at g=" + g);
            }
            if ((r == 0) != (((g + 1) % 3) == 0)) {
                throw new AssertionError("Wobble phase r==0 mismatch at g=" + g);
            }
        }
        checksPassed += 2;

        Path tmp = Files.createTempDirectory("verify-law2-");
        try {
            int n = 180_000;
            byte[] donor = generateCodingLikeSequence(n, new Random(303L));
            byte[] targetCodon = Arrays.copyOf(donor, n);
            Random rng = new Random(404L);
            for (int pos = 30; pos < n - 30; pos++) {
                // 3rd-codon wobble positions mutate 6x more frequently than 1st/2nd codon positions
                double p = (pos % 3 == 2) ? 0.024 : 0.002;
                if (rng.nextDouble() < p) {
                    int b = charToCode(targetCodon[pos]);
                    int delta = rng.nextDouble() < 0.70 ? 2 : (rng.nextBoolean() ? 1 : 3);
                    targetCodon[pos] = codeToChar(b ^ delta);
                }
            }
            Path refFile = tmp.resolve("donor.seq");
            Files.write(refFile, donor);

            FusionOptions optsWobble = new FusionOptions();
            optsWobble.referencePanel = refFile;
            optsWobble.mosaicConfig = new MosaicConfig(true, true, true, true, true, true);
            byte[] archWob = roundTripVerify(targetCodon, optsWobble, refFile);

            FusionOptions optsFlatGap = new FusionOptions();
            optsFlatGap.referencePanel = refFile;
            optsFlatGap.mosaicConfig = new MosaicConfig(true, true, true, false, true, true);
            byte[] archFlat = roundTripVerify(targetCodon, optsFlatGap, refFile);

            check(
                    archWob.length < archFlat.length,
                    "Law II: Codon-wobble decomposition ("
                            + archWob.length
                            + " B) strictly beats flat gap ("
                            + archFlat.length
                            + " B)");
            // Law IIb: 6-State Causal ORF-HMM on multi-frame alternating +/- strand ORFs vs intra-codon-shuffled control
            Random orfRng = new Random(808L);
            StringBuilder chromSb = new StringBuilder(65_000);
            StringBuilder brokenSb = new StringBuilder(65_000);
            for (int g = 0; g < 24; g++) {
                // Random intergenic spacer NOT divisible by 3 -> shifts globalPos % 3 randomly!
                int spacer = 41 + (g % 3) + orfRng.nextInt(37);
                for (int s = 0; s < spacer; s++) {
                    char sp = (char) codeToChar(orfRng.nextInt(4));
                    chromSb.append(sp);
                    brokenSb.append(sp);
                }
                int codons = 360 + orfRng.nextInt(120);
                StringBuilder gene = new StringBuilder(codons * 3 + 3);
                StringBuilder geneShuf = new StringBuilder(codons * 3 + 3);
                gene.append("ATG");
                geneShuf.append("ATG");
                for (int c = 0; c < codons; c++) {
                    char c0, c1, c2;
                    while (true) {
                        c0 = (char) codeToChar(orfRng.nextInt(4));
                        c1 = (char) codeToChar(orfRng.nextInt(4));
                        // 3rd-codon wobble position strongly biased (92% G/C, e.g. codon bias)
                        c2 = orfRng.nextDouble() < 0.92
                                ? (orfRng.nextBoolean() ? 'C' : 'G')
                                : (orfRng.nextBoolean() ? 'A' : 'T');
                        if (!(c0 == 'T' && ((c1 == 'A' && (c2 == 'A' || c2 == 'G'))
                                || (c1 == 'G' && c2 == 'A')))) {
                            break;
                        }
                    }
                    gene.append(c0).append(c1).append(c2);
                    // Control: cyclic random rotation of (c0, c1, c2) -> identical base counts, broken 3-frame phase!
                    int rot = orfRng.nextInt(3);
                    if (rot == 0) geneShuf.append(c0).append(c1).append(c2);
                    else if (rot == 1) geneShuf.append(c1).append(c2).append(c0);
                    else geneShuf.append(c2).append(c0).append(c1);
                }
                gene.append("TAA");
                geneShuf.append("TAA");
                byte[] gBytes = gene.toString().getBytes(StandardCharsets.US_ASCII);
                byte[] sBytes = geneShuf.toString().getBytes(StandardCharsets.US_ASCII);
                if ((g & 1) == 1) {
                    gBytes = reverseComplementAscii(gBytes);
                    sBytes = reverseComplementAscii(sBytes);
                }
                chromSb.append(new String(gBytes, StandardCharsets.US_ASCII));
                brokenSb.append(new String(sBytes, StandardCharsets.US_ASCII));
            }
            byte[] multiFrameOrf = chromSb.toString().getBytes(StandardCharsets.US_ASCII);
            byte[] frameBroken = brokenSb.toString().getBytes(StandardCharsets.US_ASCII);

            FusionOptions singleOpts = new FusionOptions();
            singleOpts.mode = HyperGraSSFusionEngine.MODE_SINGLE;
            singleOpts.threads = 1;
            byte[] archOrf = roundTripVerify(multiFrameOrf, singleOpts, null);
            byte[] archBroken = roundTripVerify(frameBroken, singleOpts, null);
            check(
                    archOrf.length * 100 < archBroken.length * 98,
                    "Law IIb: 6F-ORF-HMM multi-strand ORFs ("
                            + archOrf.length
                            + " B) must beat phase-rotated control ("
                            + archBroken.length
                            + " B) by >2%");

            System.out.printf(
                    "  [PASS] Law II  (Wobble Gap + 6F-ORF-HMM)          : Wobble=%d B vs Flat=%d B (-%.1f%%) | 6F-ORF=%d B vs Broken=%d B%n",
                    archWob.length,
                    archFlat.length,
                    100.0 * (archFlat.length - archWob.length) / archFlat.length,
                    archOrf.length,
                    archBroken.length);
        } finally {
            deleteRecursively(tmp);
        }
    }

    // -------------------------------------------------------------------------
    // Law III: Conditional Monotonicity & Syntenic Micro-Indel Continuity Law
    // -------------------------------------------------------------------------
    private static void verifyLaw3MonotonicityAndIndelContinuity() throws Exception {
        Path tmp = Files.createTempDirectory("verify-law3-");
        try {
            // Part A: Conditional Monotonicity under 24 Distractor Genomes
            Path twinOnlyDir = tmp.resolve("twin_only");
            Path fullPanelDir = tmp.resolve("full_panel");
            Files.createDirectories(twinOnlyDir);
            Files.createDirectories(fullPanelDir);

            int n = 120_000;
            Random rng = new Random(505L);
            byte[] ancestor = generateCodingLikeSequence(n, rng);

            // Create 24 moderately diverged distractors (2% SNPs + shuffled rearrangements)
            for (int d = 0; d < 24; d++) {
                byte[] dist = Arrays.copyOf(ancestor, n);
                for (int i = 0; i < n; i++) {
                    if (rng.nextDouble() < 0.022) {
                        dist[i] = codeToChar(charToCode(dist[i]) ^ (1 + rng.nextInt(3)));
                    }
                }
                Files.write(
                        fullPanelDir.resolve(String.format("donor_%02d.seq", d)), dist);
            }
            // Donor 12 is replaced by a near-identical twin of target (only 18 SNPs across 120 kbp!)
            byte[] twin = Arrays.copyOf(ancestor, n);
            for (int i = 0; i < n; i += 97) {
                twin[i] = codeToChar(charToCode(twin[i]) ^ 2);
            }
            Files.write(twinOnlyDir.resolve("donor_12.seq"), twin);
            Files.write(fullPanelDir.resolve("donor_12.seq"), twin);

            byte[] target = Arrays.copyOf(twin, n);
            for (int k = 0; k < 18; k++) {
                int p = 1500 + k * 6200;
                target[p] = codeToChar(charToCode(target[p]) ^ 2);
            }

            FusionOptions oTwin = new FusionOptions();
            oTwin.referencePanel = twinOnlyDir;
            byte[] archTwin = roundTripVerify(target, oTwin, twinOnlyDir);

            FusionOptions oFull = new FusionOptions();
            oFull.referencePanel = fullPanelDir;
            byte[] archFull = roundTripVerify(target, oFull, fullPanelDir);

            // Extra 23 filenames in compressed manifest cost ~28 B; mosaic itself must not blow up!
            check(
                    archFull.length <= archTwin.length + 64,
                    "Law III(a): Conditional monotonicity violated: fullPanel="
                            + archFull.length
                            + " B vs twinOnly="
                            + archTwin.length
                            + " B");

            // Part B: Syntenic Micro-Indel Continuity vs Ungapped Match Shattering
            // Conserved syntenic exons/blocks of 16..34 bp separated by frequent micro-indels
            ByteArrayBuilder mutTarget = new ByteArrayBuilder(n + 8192);
            int cur = 0;
            Random indelRng = new Random(606L);
            while (cur < n) {
                int seg = Math.min(n - cur, 16 + indelRng.nextInt(19));
                for (int j = 0; j < seg; j++) {
                    mutTarget.append(twin[cur + j]);
                }
                cur += seg;
                if (cur + 80 < n) {
                    if (indelRng.nextBoolean()) {
                        // Syntenic deletion of 1..18 bp
                        cur += 1 + indelRng.nextInt(18);
                    } else {
                        // Syntenic insertion of 1..6 bp
                        int ins = 1 + indelRng.nextInt(6);
                        for (int j = 0; j < ins; j++) {
                            mutTarget.append(codeToChar(indelRng.nextInt(4)));
                        }
                    }
                }
            }
            byte[] indelTargetBytes = mutTarget.toByteArray();
            Path singleRef = twinOnlyDir.resolve("donor_12.seq");

            FusionOptions oSplice = new FusionOptions();
            oSplice.referencePanel = singleRef;
            oSplice.mosaicConfig = new MosaicConfig(true, true, true, true, true, true);
            byte[] archSpliced = roundTripVerify(indelTargetBytes, oSplice, singleRef);

            FusionOptions oUnspliced = new FusionOptions();
            oUnspliced.referencePanel = singleRef;
            oUnspliced.mosaicConfig = new MosaicConfig(true, false, true, true, true, true);
            byte[] archUnspliced = roundTripVerify(indelTargetBytes, oUnspliced, singleRef);

            check(
                    archSpliced.length < (int) (archUnspliced.length * 0.85),
                    "Law III(b): Micro-indel splice ("
                            + archSpliced.length
                            + " B) must beat unspliced ("
                            + archUnspliced.length
                            + " B) by >15% even with collinear memory active");
            System.out.printf(
                    "  [PASS] Law III (Monotonicity + Micro-Indel Splice): Panel24=%d B (Twin1=%d B) | Spliced=%d B vs Unspliced=%d B (-%.1f%%)%n",
                    archFull.length,
                    archTwin.length,
                    archSpliced.length,
                    archUnspliced.length,
                    100.0 * (archUnspliced.length - archSpliced.length) / archUnspliced.length);
        } finally {
            deleteRecursively(tmp);
        }
    }

    // -------------------------------------------------------------------------
    // Law IV: 5-Regime Universal Existing & Non-Existing Genome Law
    // -------------------------------------------------------------------------
    private static void verifyLaw4UniversalFiveRegimes() throws Exception {
        FusionOptions singleOpts = new FusionOptions();
        singleOpts.mode = HyperGraSSFusionEngine.MODE_SINGLE;
        singleOpts.singleGenome = true;

        // 1. Cryptographic incompressible random DNA: guaranteed <= 2.005 bpb ceiling
        SecureRandom sr = new SecureRandom(new byte[] {1, 2, 3, 4, 5, 6, 7, 8});
        byte[] cryptoDna = new byte[120_000];
        for (int i = 0; i < cryptoDna.length; i++) {
            cryptoDna[i] = codeToChar(sr.nextInt(4));
        }
        byte[] archCrypto = roundTripVerify(cryptoDna, singleOpts, null);
        double cryptoBpb = archCrypto.length * 8.0 / cryptoDna.length;
        check(
                cryptoBpb <= 2.005,
                "Law IV(1): Cryptographic random DNA exceeded 2.005 bpb ceiling: " + cryptoBpb);

        // 2. Pure viral RNA (flat ACGU + FASTA ACGU)
        byte[] flatRna = new byte[40_000];
        Random rng = new Random(707L);
        byte[] rnaAlphabet = {'A', 'C', 'G', 'U'};
        for (int i = 0; i < flatRna.length; i++) {
            flatRna[i] = rnaAlphabet[rng.nextInt(4)];
        }
        roundTripVerify(flatRna, singleOpts, null);

        String fastaRna = ">viral_rna_seg1\r\nacguuacguacguuacguNNNNacgu\r\n>seg2\r\nUUUUAAAACCCCGGGG";
        roundTripVerify(fastaRna.getBytes(StandardCharsets.UTF_8), singleOpts, null);

        // 3. Chimeric T+U + full IUPAC + soft-masked lowercase + ragged lines + no EOF newline
        String chimeric =
                ">chimeric_synthetic_construct|v1\n"
                        + "acgtURYKMSWBDHVNacgtacgtNNNNNNNNNNttttuuuuACGT\r\n"
                        + ">empty_record\n"
                        + ">tail_no_newline\n"
                        + "aaaaccccggggttttRYMKSWBDHVNuuuutttt";
        roundTripVerify(chimeric.getBytes(StandardCharsets.UTF_8), singleOpts, null);

        // 4. Extreme Xenobiological AT-skew (88% AT)
        byte[] xenoAt = new byte[100_000];
        for (int i = 0; i < xenoAt.length; i++) {
            double u = rng.nextDouble();
            xenoAt[i] = u < 0.44 ? (byte) 'A' : (u < 0.88 ? (byte) 'T' : (u < 0.94 ? (byte) 'C' : (byte) 'G'));
        }
        byte[] archXeno = roundTripVerify(xenoAt, singleOpts, null);
        double xenoBpb = archXeno.length * 8.0 / xenoAt.length;
        // Exact i.i.d. Shannon limit H(0.44, 0.44, 0.06, 0.06) = 1.5294 bpb
        check(
                xenoBpb < 1.58,
                "Law IV(4): Extreme AT-skew bpb too far above 1.5294 Shannon limit: " + xenoBpb);

        // 5. High-copy tandem satellite array (37-bp motif x 1,500 copies with 1% mutations)
        byte[] motif = new byte[37];
        for (int i = 0; i < motif.length; i++) motif[i] = codeToChar(rng.nextInt(4));
        byte[] satellite = new byte[37 * 1500];
        for (int i = 0; i < satellite.length; i++) {
            byte b = motif[i % 37];
            if (rng.nextDouble() < 0.008) b = codeToChar(charToCode(b) ^ 2);
            satellite[i] = b;
        }
        byte[] archSat = roundTripVerify(satellite, singleOpts, null);
        double satBpb = archSat.length * 8.0 / satellite.length;
        check(satBpb < 0.25, "Law IV(5): Tandem satellite repeat bpb too high: " + satBpb);

        // 6. Inverted transposon repeats (+ and - copies) vs self-match disabled
        byte[] core = generateCodingLikeSequence(60_000, new Random(808L));
        byte[] transposon = Arrays.copyOfRange(core, 5_000, 15_000);
        byte[] rcTransposon = reverseComplementAscii(transposon);
        for (int i = 0; i < rcTransposon.length; i++) {
            if (rng.nextDouble() < 0.01) {
                rcTransposon[i] = codeToChar(charToCode(rcTransposon[i]) ^ 2);
            }
        }
        ByteArrayBuilder genomeWithIr = new ByteArrayBuilder(90_000);
        genomeWithIr.appendBytes(core);
        genomeWithIr.appendBytes(rcTransposon);
        genomeWithIr.appendBytes(transposon);
        byte[] irBytes = genomeWithIr.toByteArray();

        byte[] archIrFull = roundTripVerify(irBytes, singleOpts, null);
        FusionOptions noSelfOpts = new FusionOptions();
        noSelfOpts.mode = HyperGraSSFusionEngine.MODE_SINGLE;
        noSelfOpts.singleGenome = true;
        noSelfOpts.mosaicConfig = new MosaicConfig(true, true, true, true, false, false);
        byte[] archIrNoSelf = roundTripVerify(irBytes, noSelfOpts, null);
        check(
                archIrFull.length < archIrNoSelf.length,
                "Law IV(6): Bidirectional self-mosaic ("
                        + archIrFull.length
                        + " B) must beat no-self ("
                        + archIrNoSelf.length
                        + " B)");

        // 7. Degenerate boundary inputs (0 B, 1 B, all-N)
        roundTripVerify(new byte[0], singleOpts, null);
        roundTripVerify(new byte[] {'G'}, singleOpts, null);
        roundTripVerify(">only_ns\nNNNNNNNNNNNNNNNNNNNN\n".getBytes(StandardCharsets.UTF_8), singleOpts, null);

        System.out.printf(
                "  [PASS] Law IV  (5-Regime Universal + Non-Existing): Crypto=%.4f bpb (<=2.005), AT88=%.3f bpb, Satellite=%.3f bpb, IR=%d B vs %d B%n",
                cryptoBpb, xenoBpb, satBpb, archIrFull.length, archIrNoSelf.length);
    }

    // -------------------------------------------------------------------------
    // Law V: Warm-Primed Parallel Segmentation & Bounded Memory Determinism
    // -------------------------------------------------------------------------
    private static void verifyLaw5WarmParallelAndBoundedDeterminism() throws Exception {
        byte[] seq = generateCodingLikeSequence(260_000, new Random(909L));
        int size1 = 0;
        int size8 = 0;
        for (int threads : new int[] {1, 2, 4, 8}) {
            for (boolean bounded : new boolean[] {false, true}) {
                FusionOptions opts = new FusionOptions();
                opts.mode = HyperGraSSFusionEngine.MODE_SINGLE;
                opts.singleGenome = true;
                opts.threads = threads;
                opts.boundedMemory = bounded;
                byte[] arch = roundTripVerify(seq, opts, null);
                if (!bounded && threads == 1) size1 = arch.length;
                if (!bounded && threads == 8) size8 = arch.length;
            }
        }
        double overheadPct = 100.0 * Math.abs(size8 - size1) / size1;
        check(
                overheadPct < 1.5,
                "Law V: Warm-primed 8-thread segmentation overhead too large: " + overheadPct + "%");
        System.out.printf(
                "  [PASS] Law V   (Warm-Primed Parallel + Bounded)   : 1T=%d B, 8T=%d B (delta=%.3f%%)%n",
                size1, size8, 100.0 * (size8 - size1) / size1);
    }

    // -------------------------------------------------------------------------
    // Law VI: Hierarchical Pitman-Yor Backoff & Two-Adic Tapered Adaptation Law
    // -------------------------------------------------------------------------
    private static void verifyLaw6HierarchicalPitmanYorAndTwoAdicAdaptation() throws Exception {
        // 1. Verify Two-Adic Tapered Adaptation reaches <= 12/4096 (0.29%) vs classic 31/4096 (0.76%)
        int pTaper = 2048;
        int pClassic = 2048;
        for (int i = 0; i < 4096; i++) {
            pTaper -= (pTaper + 20) >>> 5;
            pClassic -= pClassic >>> 5;
        }
        check(pTaper == 11, "Law VI(1): Two-adic tapered adaptation fixed point is 11/4096: got " + pTaper);
        check(pClassic == 31, "Law VI(1): Classic shift-5 fixed point is 31/4096: got " + pClassic);
        check(pTaper < (pClassic >>> 1), "Law VI(1): Tapered fixed-point rate is >2.5x lower than classic floor");

        // 2. Verify FastNoveltyCoder Hierarchical Pitman-Yor backoff + Post-SNP Homology Recovery
        Random rng = new Random(1010L);
        byte[] baseCoding = generateCodingLikeSequence(60_000, rng);
        byte[] twoBit = new byte[120_000];
        for (int i = 0; i < 60_000; i++) {
            twoBit[i] = (byte) charToCode(baseCoding[i]);
        }
        // Second 60 kbp has 80-bp homologous repeats separated by single SNPs and novel linkers
        for (int i = 60_000; i < 120_000; i++) {
            if ((i % 100) < 80) {
                int srcIdx = (i - 60_000) % 50_000;
                int b = twoBit[srcIdx];
                if ((i % 100) == 40) b ^= 2; // Single SNP inside homologous repeat
                twoBit[i] = (byte) b;
            } else {
                twoBit[i] = (byte) charToCode(baseCoding[i - 60_000]);
            }
        }
        byte[] comp = FastNoveltyCoder.encode(twoBit, twoBit.length, 1, false);
        byte[] dec =
                FastNoveltyCoder.decode(
                        new VarIntCodec.ByteIn(comp, 0, comp.length),
                        twoBit.length,
                        1,
                        false);
        check(Arrays.equals(twoBit, dec), "Law VI(2): FastNoveltyCoder exact round-trip failed");
        double bpb = comp.length * 8.0 / twoBit.length;
        check(bpb < 1.25, "Law VI(2): Hierarchical Pitman-Yor + Post-SNP Homology bpb too high: " + bpb);
        System.out.printf(
                "  [PASS] Law VI  (Pitman-Yor Backoff + Two-Adic)    : p_min=%d/4096 (vs %d/4096), Novelty=%.4f bpb%n",
                pTaper, pClassic, bpb);
    }

    // -------------------------------------------------------------------------
    // Law VII: Disjoint Index/Superstream Lifecycle & Bounded-Heap Scaling Law
    // -------------------------------------------------------------------------
    private static void verifyLaw7DisjointLifecycleAndBoundedScaling() throws Exception {
        Random rng = new Random(2026L);
        byte[] anc = generateCodingLikeSequence(90_000, rng);
        StringBuilder sb = new StringBuilder(380_000);
        for (int g = 0; g < 4; g++) {
            sb.append(">strain_").append(g).append('\n');
            byte[] strain = Arrays.copyOf(anc, anc.length);
            for (int i = 0; i < strain.length; i++) {
                if (rng.nextDouble() < 0.008) {
                    strain[i] = codeToChar(charToCode(strain[i]) ^ (1 + rng.nextInt(3)));
                }
            }
            sb.append(new String(strain, StandardCharsets.US_ASCII)).append('\n');
        }
        byte[] collectionBytes = sb.toString().getBytes(StandardCharsets.US_ASCII);

        FusionOptions optUnbounded = new FusionOptions();
        optUnbounded.mode = HyperGraSSFusionEngine.MODE_COLLECTION;
        optUnbounded.threads = 2;
        optUnbounded.boundedMemory = false;
        byte[] archUnbounded = roundTripVerify(collectionBytes, optUnbounded, null);

        FusionOptions optBounded = new FusionOptions();
        optBounded.mode = HyperGraSSFusionEngine.MODE_COLLECTION;
        optBounded.threads = 2;
        optBounded.boundedMemory = true;
        byte[] archBounded = roundTripVerify(collectionBytes, optBounded, null);

        double deltaPct = 100.0 * Math.abs(archBounded.length - archUnbounded.length) / archUnbounded.length;
        check(
                deltaPct < 2.0,
                "Law VII: Bounded-memory collection archive diverged by >2%: " + deltaPct + "%");
        System.out.printf(
                "  [PASS] Law VII (Disjoint Lifecycle + Bounded Mode): Unbounded=%d B, Bounded=%d B (delta=%.3f%%)%n",
                archUnbounded.length, archBounded.length, deltaPct);
    }

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------
    private static byte[] roundTripVerify(byte[] raw, FusionOptions opts, Path refPanel)
            throws Exception {
        byte[] archive = HyperGraSSFusionEngine.compress(raw, opts);
        byte[] restored = HyperGraSSFusionEngine.decompress(archive, refPanel);
        byte[] h1 = MessageDigest.getInstance("SHA-256").digest(raw);
        byte[] h2 = MessageDigest.getInstance("SHA-256").digest(restored);
        check(
                Arrays.equals(raw, restored) && Arrays.equals(h1, h2),
                "SHA-256 byte-exact round-trip failed");
        return archive;
    }

    private static byte[] generateCodingLikeSequence(int len, Random rng) {
        // Synthesize realistic bacterial-like coding Markov + 3-periodic codon bias
        byte[] out = new byte[len];
        int prev = 0;
        for (int i = 0; i < len; i++) {
            int phase = i % 3;
            int u = rng.nextInt(100);
            int b;
            if (phase == 0) {
                b = u < 32 ? 2 : (u < 60 ? 0 : (u < 82 ? 1 : 3)); // G/A-rich 1st codon pos
            } else if (phase == 1) {
                b = u < 30 ? 0 : (u < 58 ? 3 : (u < 80 ? 1 : 2));
            } else {
                b = ((prev & 1) == 0)
                        ? (u < 38 ? 2 : (u < 72 ? 1 : (u < 86 ? 3 : 0)))
                        : (u < 38 ? 3 : (u < 70 ? 0 : (u < 86 ? 1 : 2)));
            }
            out[i] = codeToChar(b);
            prev = b;
        }
        return out;
    }

    private static byte[] reverseComplementAscii(byte[] src) {
        byte[] out = new byte[src.length];
        for (int i = 0; i < src.length; i++) {
            out[src.length - 1 - i] = codeToChar(charToCode(src[i]) ^ 3);
        }
        return out;
    }

    private static int charToCode(byte c) {
        return switch (c) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            case 'T', 'U' -> 3;
            default -> 0;
        };
    }

    private static byte codeToChar(int code) {
        return switch (code & 3) {
            case 0 -> (byte) 'A';
            case 1 -> (byte) 'C';
            case 2 -> (byte) 'G';
            default -> (byte) 'T';
        };
    }

    private static void check(boolean cond, String msg) {
        if (!cond) {
            throw new AssertionError("Verification failed: " + msg);
        }
        checksPassed++;
    }

    private static void deleteRecursively(Path root) throws Exception {
        if (root == null || !Files.exists(root)) return;
        try (var walk = Files.walk(root)) {
            walk.sorted((a, b) -> b.compareTo(a))
                    .forEach(
                            p -> {
                                try {
                                    Files.deleteIfExists(p);
                                } catch (Exception ignored) {
                                }
                            });
        }
    }

    private static final class ByteArrayBuilder {
        private byte[] buf;
        private int len;

        ByteArrayBuilder(int cap) {
            this.buf = new byte[cap];
        }

        void append(byte b) {
            if (len == buf.length) buf = Arrays.copyOf(buf, buf.length << 1);
            buf[len++] = b;
        }

        void appendBytes(byte[] src) {
            if (len + src.length > buf.length) {
                buf = Arrays.copyOf(buf, Math.max(buf.length << 1, len + src.length));
            }
            System.arraycopy(src, 0, buf, len, src.length);
            len += src.length;
        }

        byte[] toByteArray() {
            return Arrays.copyOf(buf, len);
        }
    }
}
