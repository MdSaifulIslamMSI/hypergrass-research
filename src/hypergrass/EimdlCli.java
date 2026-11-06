package hypergrass;

import hypergrass.EimdlMosaicCodec.MosaicConfig;
import hypergrass.HyperGraSSFusionEngine.CompressionReport;
import hypergrass.HyperGraSSFusionEngine.FusionOptions;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Drop-in JSON-telemetry CLI adapter for {@link HyperGraSSFusionEngine}, compatible with all
 * benchmark harnesses ({@code run.py}, {@code collection.py}, {@code cluster_loo.py}) as well as
 * ablation experiments across Laws I–V.
 */
public final class EimdlCli {

    private EimdlCli() {}

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            printUsage();
            System.exit(1);
        }
        String sub = args[0];
        switch (sub) {
            case "compress" -> runCompress(args);
            case "decompress" -> runDecompress(args);
            default -> {
                printUsage();
                System.exit(1);
            }
        }
    }

    private static void runCompress(String[] args) throws IOException {
        Path target = null;
        Path out = null;
        FusionOptions opts = new FusionOptions();

        boolean allowApprox = true;
        boolean allowSplice = true;
        boolean allowV4 = true;
        boolean allowWobble = true;
        boolean allowSelf = true;
        boolean allowTandem = true;

        int i = 1;
        while (i < args.length) {
            String flag = args[i];
            switch (flag) {
                case "--target", "--input" -> target = Path.of(args[++i]);
                case "--out", "--output" -> out = Path.of(args[++i]);
                case "--panel", "--ref" -> opts.referencePanel = Path.of(args[++i]);
                case "--exclude" -> opts.excludePanelNames.add(args[++i]);
                case "--panel-list" -> opts.explicitPanelNames = readPanelList(Path.of(args[++i]));
                case "--panel-limit" -> opts.panelLimit = Integer.parseInt(args[++i]);
                case "--threads" -> opts.threads = Math.max(1, Integer.parseInt(args[++i]));
                case "--bounded" -> opts.boundedMemory = true;
                case "--collection" -> opts.mode = HyperGraSSFusionEngine.MODE_COLLECTION;
                case "--single" -> {
                    opts.mode = HyperGraSSFusionEngine.MODE_SINGLE;
                    opts.singleGenome = true;
                }
                case "--k", "--step" -> {
                    // Accepted for seamless CLI compatibility with Rcmdl Python harnesses
                    i++;
                }
                case "--no-approx" -> allowApprox = false;
                case "--no-splice" -> allowSplice = false;
                case "--no-v4" -> allowV4 = false;
                case "--no-wobble" -> allowWobble = false;
                case "--no-self" -> allowSelf = false;
                case "--no-tandem" -> allowTandem = false;
                default -> throw new IllegalArgumentException("Unknown compress flag: " + flag);
            }
            i++;
        }
        if (target == null || out == null) {
            throw new IllegalArgumentException("compress requires --target <path> and --out <path>");
        }
        opts.mosaicConfig =
                new MosaicConfig(
                        allowApprox, allowSplice, allowV4, allowWobble, allowSelf, allowTandem);

        CompressionReport rep = HyperGraSSFusionEngine.compressFile(target, out, opts);
        long selfAndTandemBases = rep.selfMatchBases() + rep.tandemMatchBases();
        long selfAndTandemTokens = rep.selfTokens() + rep.tandemTokens();

        String json =
                String.format(
                        Locale.ROOT,
                        "{\"mode\":%d,\"raw_bytes\":%d,\"target_bases\":%d,\"archive_bytes\":%d,"
                                + "\"panel_genomes\":%d,\"panel_bases\":%d,\"panel_match_bases\":%d,"
                                + "\"self_match_bases\":%d,\"tandem_match_bases\":%d,\"literal_bases\":%d,"
                                + "\"panel_tokens\":%d,\"self_tokens\":%d,\"tandem_tokens\":%d,"
                                + "\"literal_runs\":%d,\"indel_splices\":%d,\"mismatches\":%d,"
                                + "\"compress_seconds\":%.6f}",
                        rep.mode(),
                        rep.rawBytes(),
                        rep.targetDnaBases(),
                        rep.archiveBytes(),
                        rep.panelGenomes(),
                        rep.panelBases(),
                        rep.panelMatchBases(),
                        selfAndTandemBases,
                        rep.tandemMatchBases(),
                        rep.noveltyBases(),
                        rep.panelTokens(),
                        selfAndTandemTokens,
                        rep.tandemTokens(),
                        rep.literalRuns(),
                        rep.indelSplices(),
                        rep.substitutions(),
                        rep.compressSeconds());
        System.out.println(json);
    }

    private static void runDecompress(String[] args) throws IOException {
        Path archive = null;
        Path panel = null;
        Path out = null;

        int i = 1;
        while (i < args.length) {
            String flag = args[i];
            switch (flag) {
                case "--archive", "--input" -> archive = Path.of(args[++i]);
                case "--panel", "--ref" -> panel = Path.of(args[++i]);
                case "--out", "--output" -> out = Path.of(args[++i]);
                default -> throw new IllegalArgumentException("Unknown decompress flag: " + flag);
            }
            i++;
        }
        if (archive == null || out == null) {
            throw new IllegalArgumentException(
                    "decompress requires --archive <path> and --out <path>");
        }
        long t0 = System.nanoTime();
        HyperGraSSFusionEngine.decompressFile(archive, out, panel);
        double secs = (System.nanoTime() - t0) / 1e9;
        long archBytes = Files.size(archive);
        long outBytes = Files.size(out);
        System.out.printf(
                Locale.ROOT,
                "{\"archive_bytes\":%d,\"output_bytes\":%d,\"decompress_seconds\":%.6f}%n",
                archBytes,
                outBytes,
                secs);
    }

    private static List<String> readPanelList(Path path) throws IOException {
        List<String> out = new ArrayList<>();
        for (String line : Files.readAllLines(path, StandardCharsets.UTF_8)) {
            String t = line.trim();
            if (!t.isEmpty()) out.add(t);
        }
        return out;
    }

    private static void printUsage() {
        System.err.println(
                "Usage:\n"
                        + "  java -cp bin hypergrass.EimdlCli compress --target <file> --out <archive> "
                        + "[--panel <dir|file>] [--exclude <name>] [--panel-list <file>] "
                        + "[--panel-limit <N>] [--collection] [--single] [--bounded] [--threads <N>]\n"
                        + "  java -cp bin hypergrass.EimdlCli decompress --archive <archive> "
                        + "[--panel <dir|file>] --out <file>");
    }
}
