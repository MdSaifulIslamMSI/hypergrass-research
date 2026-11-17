package hypergrass;

import hypergrass.CompressedDomainEngine.AnalyzeReport;
import hypergrass.CompressedDomainEngine.KmerHit;
import hypergrass.CompressedDomainEngine.KmerReport;
import hypergrass.CompressedDomainEngine.VcfResult;
import hypergrass.EimdlMosaicCodec.MosaicConfig;
import hypergrass.Hgs3OmicCodec.Hgs3Options;
import hypergrass.Hgs3OmicCodec.Hgs3Report;
import hypergrass.Hgs3OmicCodec.OmicFormat;
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
 * Unified JSON/Text CLI adapter for {@link HyperGraSSFusionEngine}, {@link Hgs3OmicCodec}, and
 * {@link CompressedDomainEngine}.
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
            case "decompress", "decompress-omics" -> runDecompress(args);
            case "compress-omics" -> runCompressOmics(args);
            case "analyze" -> runAnalyze(args);
            case "variants" -> runVariants(args);
            case "kmer" -> runKmer(args);
            case "query" -> runQuery(args);
            default -> {
                printUsage();
                System.exit(1);
            }
        }
    }

    private static final class CompressArgs {
        Path target;
        Path out;
        final FusionOptions opts = new FusionOptions();
        boolean allowApprox = true;
        boolean allowSplice = true;
        boolean allowV4 = true;
        boolean allowWobble = true;
        boolean allowSelf = true;
        boolean allowTandem = true;
        OmicFormat omicFormat = OmicFormat.AUTO;
        boolean cds = false;

        void finalizeMosaic() {
            opts.mosaicConfig =
                    new MosaicConfig(
                            allowApprox, allowSplice, allowV4, allowWobble, allowSelf, allowTandem);
        }
    }

    private static CompressArgs parseCompressArgs(String[] args, boolean omicsMode) throws IOException {
        CompressArgs ca = new CompressArgs();
        int i = 1;
        while (i < args.length) {
            String flag = args[i];
            switch (flag) {
                case "--target", "--input" -> ca.target = Path.of(args[++i]);
                case "--out", "--output" -> ca.out = Path.of(args[++i]);
                case "--panel", "--ref" -> ca.opts.referencePanel = Path.of(args[++i]);
                case "--exclude" -> ca.opts.excludePanelNames.add(args[++i]);
                case "--panel-list" -> ca.opts.explicitPanelNames = readPanelList(Path.of(args[++i]));
                case "--panel-limit" -> ca.opts.panelLimit = Integer.parseInt(args[++i]);
                case "--threads" -> ca.opts.threads = Math.max(1, Integer.parseInt(args[++i]));
                case "--bounded" -> ca.opts.boundedMemory = true;
                case "--collection" -> ca.opts.mode = HyperGraSSFusionEngine.MODE_COLLECTION;
                case "--single" -> {
                    ca.opts.mode = HyperGraSSFusionEngine.MODE_SINGLE;
                    ca.opts.singleGenome = true;
                }
                case "--k", "--step" -> i++;
                case "--no-approx" -> ca.allowApprox = false;
                case "--no-splice" -> ca.allowSplice = false;
                case "--no-v4" -> ca.allowV4 = false;
                case "--no-wobble" -> ca.allowWobble = false;
                case "--no-self" -> ca.allowSelf = false;
                case "--no-tandem" -> ca.allowTandem = false;
                case "--format" -> {
                    if (!omicsMode) throw new IllegalArgumentException("Unknown compress flag: " + flag);
                    String fv = args[++i].toLowerCase(Locale.ROOT);
                    ca.omicFormat =
                            switch (fv) {
                                case "auto" -> OmicFormat.AUTO;
                                case "fasta", "fa", "fna" -> OmicFormat.FASTA;
                                case "fastq", "fq" -> OmicFormat.FASTQ;
                                default -> throw new IllegalArgumentException("Unknown --format: " + fv);
                            };
                }
                case "--cds" -> {
                    if (!omicsMode) throw new IllegalArgumentException("Unknown compress flag: " + flag);
                    ca.cds = true;
                }
                default -> throw new IllegalArgumentException("Unknown compress flag: " + flag);
            }
            i++;
        }
        if (ca.target == null || ca.out == null) {
            throw new IllegalArgumentException("Requires --target <path> and --out <path>");
        }
        ca.finalizeMosaic();
        return ca;
    }

    private static void runCompress(String[] args) throws IOException {
        CompressArgs ca = parseCompressArgs(args, false);
        CompressionReport rep = HyperGraSSFusionEngine.compressFile(ca.target, ca.out, ca.opts);
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

    private static void runCompressOmics(String[] args) throws IOException {
        CompressArgs ca = parseCompressArgs(args, true);
        Hgs3Options h3 = new Hgs3Options();
        h3.format = ca.omicFormat;
        h3.fusion = ca.opts;
        h3.emitCdsDigest = ca.cds;
        Hgs3Report rep = Hgs3OmicCodec.compressOmics(ca.target, ca.out, h3);
        System.out.printf(
                Locale.ROOT,
                "{\"container\":\"HGS3\",\"is_fastq\":%b,\"mode\":%d,\"engine\":\"java-25-ans4\","
                        + "\"compress_seconds\":%.6f,\"raw_bytes\":%d,\"archive_bytes\":%d,"
                        + "\"bpb\":%.6f,\"ratio\":%.4f,\"records\":%d,\"dna_bases\":%d,"
                        + "\"qual_bytes\":%d,\"meta_comp_bytes\":%d,\"mosaic_comp_bytes\":%d,"
                        + "\"novelty_comp_bytes\":%d,\"qual_comp_bytes\":%d,\"cds_orfs\":%d}%n",
                rep.isFastq(),
                rep.mode(),
                rep.seconds(),
                rep.rawBytes(),
                rep.archiveBytes(),
                rep.bpb(),
                rep.ratio(),
                rep.totalRecords(),
                rep.dnaBases(),
                rep.qualBytes(),
                rep.metaCompBytes(),
                rep.mosaicCompBytes(),
                rep.noveltyCompBytes(),
                rep.qualCompBytes(),
                rep.cds() != null ? rep.cds().orfCount() : 0L);
    }

    private static void runDecompress(String[] args) throws IOException {
        Path archive = null;
        Path panel = null;
        Path out = null;
        int threads = Math.max(1, Runtime.getRuntime().availableProcessors());

        int i = 1;
        while (i < args.length) {
            String flag = args[i];
            switch (flag) {
                case "--archive", "--input" -> archive = Path.of(args[++i]);
                case "--panel", "--ref" -> panel = Path.of(args[++i]);
                case "--out", "--output" -> out = Path.of(args[++i]);
                case "--threads" -> threads = Math.max(1, Integer.parseInt(args[++i]));
                case "--bounded" -> {}
                default -> throw new IllegalArgumentException("Unknown decompress flag: " + flag);
            }
            i++;
        }
        if (archive == null || out == null) {
            throw new IllegalArgumentException(
                    "decompress requires --archive <path> and --out <path>");
        }
        double secs = Hgs3OmicCodec.decompressOmics(archive, out, panel, threads);
        long archBytes = Files.size(archive);
        long outBytes = Files.size(out);
        System.out.printf(
                Locale.ROOT,
                "{\"archive_bytes\":%d,\"output_bytes\":%d,\"decompress_seconds\":%.6f}%n",
                archBytes,
                outBytes,
                secs);
    }

    private static void runAnalyze(String[] args) throws IOException {
        Path archive = null;
        Path panel = null;
        boolean json = false;
        int i = 1;
        while (i < args.length) {
            switch (args[i]) {
                case "--archive", "--input" -> archive = Path.of(args[++i]);
                case "--panel", "--ref" -> panel = Path.of(args[++i]);
                case "--json" -> json = true;
                case "--bounded" -> {}
                case "--threads" -> i++;
                default -> throw new IllegalArgumentException("Unknown analyze flag: " + args[i]);
            }
            i++;
        }
        if (archive == null) throw new IllegalArgumentException("analyze requires --archive <path>");
        AnalyzeReport rep = CompressedDomainEngine.analyze(archive, panel);
        System.out.print(json ? rep.toJson() + "\n" : rep.formatSummary());
    }

    private static void runVariants(String[] args) throws IOException {
        Path archive = null;
        Path panel = null;
        Path out = null;
        boolean json = false;
        int i = 1;
        while (i < args.length) {
            switch (args[i]) {
                case "--archive", "--input" -> archive = Path.of(args[++i]);
                case "--panel", "--ref" -> panel = Path.of(args[++i]);
                case "--out", "--output" -> out = Path.of(args[++i]);
                case "--json" -> json = true;
                case "--bounded" -> {}
                case "--threads" -> i++;
                default -> throw new IllegalArgumentException("Unknown variants flag: " + args[i]);
            }
            i++;
        }
        if (archive == null) throw new IllegalArgumentException("variants requires --archive <path>");
        VcfResult res = CompressedDomainEngine.extractVariantsVcf(archive, panel, out);
        if (json) {
            System.out.printf(
                    Locale.ROOT,
                    "{\"total_variants\":%d,\"snvs\":%d,\"insertions\":%d,\"deletions\":%d,"
                            + "\"complex_indels\":%d,\"literal_runs_skipped\":%d,\"elapsed_seconds\":%.6f}%n",
                    res.totalVariants(), res.snvs(), res.insertions(), res.deletions(),
                    res.complexIndels(), res.literalRunsSkipped(), res.seconds());
        } else if (out != null) {
            System.out.printf(
                    Locale.ROOT,
                    "VCF emitted: %d variants (SNV=%d, INS=%d, DEL=%d, COMPLEX=%d) in %.4f s%n",
                    res.totalVariants(), res.snvs(), res.insertions(), res.deletions(),
                    res.complexIndels(), res.seconds());
        } else {
            System.out.print(res.vcfText());
        }
    }

    private static void runKmer(String[] args) throws IOException {
        Path archive = null;
        Path panel = null;
        String query = null;
        boolean both = false;
        boolean json = false;
        int i = 1;
        while (i < args.length) {
            switch (args[i]) {
                case "--archive", "--input" -> archive = Path.of(args[++i]);
                case "--panel", "--ref" -> panel = Path.of(args[++i]);
                case "--query", "-q" -> query = args[++i];
                case "--both-strands" -> both = true;
                case "--json" -> json = true;
                case "--bounded" -> {}
                case "--threads" -> i++;
                default -> throw new IllegalArgumentException("Unknown kmer flag: " + args[i]);
            }
            i++;
        }
        if (archive == null || query == null) {
            throw new IllegalArgumentException("kmer requires --archive <path> and --query <ACGT>");
        }
        KmerReport rep = CompressedDomainEngine.searchKmers(archive, panel, query, both);
        if (json) {
            System.out.printf(
                    Locale.ROOT,
                    "{\"query_len\":%d,\"hit_count\":%d,\"exact_zero_target_spans\":%d,"
                            + "\"patched_or_local_spans\":%d,\"elapsed_seconds\":%.6f}%n",
                    rep.queryLen(), rep.hits().size(), rep.exactZeroTargetSpans(),
                    rep.patchedOrLocalSpans(), rep.seconds());
        } else {
            System.out.printf(
                    Locale.ROOT,
                    "Found %d exact hits for %d-mer (zero-target-spans=%d, patched-spans=%d, %.4f s):%n",
                    rep.hits().size(), rep.queryLen(), rep.exactZeroTargetSpans(),
                    rep.patchedOrLocalSpans(), rep.seconds());
            int shown = 0;
            for (KmerHit h : rep.hits()) {
                if (shown++ >= 50) break;
                System.out.printf(
                        "  unit=%d\tpos=%d\tstrand=%c%n",
                        h.unitIndex(), h.tgtPos(), h.strand() == 0 ? '+' : '-');
            }
        }
    }

    private static void runQuery(String[] args) throws IOException {
        Path archive = null;
        Path panel = null;
        int unit = 0;
        int start = -1;
        int end = -1;
        int i = 1;
        while (i < args.length) {
            switch (args[i]) {
                case "--archive", "--input" -> archive = Path.of(args[++i]);
                case "--panel", "--ref" -> panel = Path.of(args[++i]);
                case "--unit" -> unit = Integer.parseInt(args[++i]);
                case "--start" -> start = Integer.parseInt(args[++i]);
                case "--end" -> end = Integer.parseInt(args[++i]);
                case "--bounded" -> {}
                case "--threads" -> i++;
                default -> throw new IllegalArgumentException("Unknown query flag: " + args[i]);
            }
            i++;
        }
        if (archive == null || start < 0 || end < 0) {
            throw new IllegalArgumentException("query requires --archive <path> --start <pos> --end <pos>");
        }
        System.out.println(CompressedDomainEngine.querySlice(archive, panel, unit, start, end));
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
                        + "  java -cp bin hypergrass.EimdlCli compress --target <file> --out <archive.hgs2> "
                        + "[--panel <dir|file>] [--collection] [--single] [--bounded] [--threads <N>]\n"
                        + "  java -cp bin hypergrass.EimdlCli decompress --archive <archive.hgs2|hgs3> "
                        + "[--panel <dir|file>] --out <file>\n"
                        + "  java -cp bin hypergrass.EimdlCli compress-omics --target <fasta|fastq> "
                        + "--out <archive.hgs3> [--format auto|fasta|fastq] [--panel <dir|file>] [--cds]\n"
                        + "  java -cp bin hypergrass.EimdlCli decompress-omics --archive <archive.hgs3|hgs2> "
                        + "[--panel <dir|file>] --out <file>\n"
                        + "  java -cp bin hypergrass.EimdlCli analyze --archive <archive> [--panel <dir|file>] [--json]\n"
                        + "  java -cp bin hypergrass.EimdlCli variants --archive <archive> --panel <dir|file> [--out <vcf>] [--json]\n"
                        + "  java -cp bin hypergrass.EimdlCli kmer --archive <archive> [--panel <dir|file>] --query <ACGT> [--both-strands] [--json]\n"
                        + "  java -cp bin hypergrass.EimdlCli query --archive <archive> [--panel <dir|file>] [--unit <N>] --start <pos> --end <pos>");
    }
}
