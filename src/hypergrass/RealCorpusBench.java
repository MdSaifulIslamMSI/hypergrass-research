package hypergrass;

import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.HexFormat;

/** Real-genome benchmark through the full public pipeline with SHA-256 round-trip verification. */
public final class RealCorpusBench {
    public static void main(String[] args) throws Exception {
        Path dir = Path.of(args[0]);
        Path tmp = Files.createTempDirectory("hgs_real_");
        // name, best published SOTA bytes (specialized compressors), who
        Object[][] cases = {
            {"BuEb", 4690, "Jarvis"}, {"AgPh", 10708, "GeCo2"}, {"YeMi", 16798, "GeCo2"},
            {"AeCa", 378374, "GeCo3"}, {"HePy", 374362, "Jarvis"}, {"HaHi", 898031, "GeCo3"},
            {"EsCo", 1095028, "GeCo3"}, {"PlFa", 1908680, "GeCo3"}, {"ScPo", 2511563, "GeCo3"},
            {"GCF_001884535.1.fna", 11983, "DCom"}, {"GCA_000398605.1.fna", 120179, "DCom"},
            {"GCA_000211355.2.fna", 384952, "GraSS"}, {"GCA_000988165.1.fna", 1344131, "LFQC"},
            {"GCA_000165345.1.fna", 2181328, "LFQC"},
            {"Real_Ebolavirus_GenBank_200.fasta", -1, "n/a"},
        };
        int wins = 0, total = 0;
        for (Object[] c : cases) {
            String name = (String) c[0];
            int sota = (Integer) c[1];
            Path in = dir.resolve(name);
            Path hgs = tmp.resolve(name + ".hgs");
            Path dec = tmp.resolve(name + ".dec");
            long t0 = System.nanoTime();
            HyperGraSS.compressFile(in, hgs, HyperGraSS.CompressOptions.defaults());
            long t1 = System.nanoTime();
            HyperGraSS.decompressFile(hgs, dec);
            long t2 = System.nanoTime();
            String a = sha(in), b = sha(dec);
            boolean ok = a.equals(b);
            long raw = Files.size(in), comp = Files.size(hgs);
            String verdict = sota < 0 ? "-" : (comp < sota ? "WIN " + (sota - comp) : "gap +" + (comp - sota));
            if (sota > 0) { total++; if (comp < sota) wins++; }
            System.out.printf("%-34s raw=%,11d hgs=%,10d CR=%7.3fx sota=%,10d(%s) %s | SHA %s | comp %.1fs dec %.1fs%n",
                    name, raw, comp, (double) raw / comp, sota, c[2], verdict, ok ? "OK" : "MISMATCH!!",
                    (t1 - t0) / 1e9, (t2 - t1) / 1e9);
            System.out.flush();
            if (!ok) throw new AssertionError("SHA mismatch " + name);
            Files.deleteIfExists(dec);
        }
        System.out.println("WINS vs best published specialized compressor: " + wins + " / " + total);
    }

    static String sha(Path p) throws Exception {
        MessageDigest md = MessageDigest.getInstance("SHA-256");
        try (var in = Files.newInputStream(p)) {
            byte[] buf = new byte[1 << 20];
            int n;
            while ((n = in.read(buf)) > 0) md.update(buf, 0, n);
        }
        return HexFormat.of().formatHex(md.digest());
    }
}
