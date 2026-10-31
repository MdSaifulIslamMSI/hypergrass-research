package ablation;

import hypergrass.HyperGraSS;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;

/**
 * One ablation point: compress the input with a chosen anchor-bank cap and
 * block size, restore it, and report size, time, and exactness as one JSON line.
 *
 * This is a benchmark harness, not a production entry point. Each invocation is
 * a fresh JVM so the caller can record peak RSS for that configuration alone.
 */
public final class AnchorBankAblation {
    private AnchorBankAblation() {}

    public static void main(String[] args) throws Exception {
        if (args.length != 5) {
            System.err.println(
                    "usage: AnchorBankAblation INPUT BANK_MB BLOCK_MB ARCHIVE RESTORED");
            System.exit(2);
        }
        Path input = Path.of(args[0]);
        int bankMb = Integer.parseInt(args[1]);
        int blockMb = Integer.parseInt(args[2]);
        Path archive = Path.of(args[3]);
        Path restored = Path.of(args[4]);

        HyperGraSS.CompressOptions opts = HyperGraSS.CompressOptions.defaults();
        opts.maxAnchorBases = bankMb * 1024 * 1024;
        opts.blockSizeBytes = blockMb * 1024 * 1024;

        long compressStart = System.nanoTime();
        HyperGraSS.compressFile(input, archive, opts);
        double compressSeconds = (System.nanoTime() - compressStart) / 1_000_000_000.0;

        long decompressStart = System.nanoTime();
        HyperGraSS.decompressFile(archive, restored);
        double decompressSeconds = (System.nanoTime() - decompressStart) / 1_000_000_000.0;

        boolean exact = Files.mismatch(input, restored) == -1;
        HyperGraSS.ArchiveMetadata meta = HyperGraSS.inspect(archive);
        System.out.printf(
                Locale.ROOT,
                "{\"bank_mb\":%d,\"block_mb\":%d,\"archive_bytes\":%d,"
                        + "\"compress_seconds\":%.3f,\"decompress_seconds\":%.3f,"
                        + "\"exact\":%s,\"records\":%d,\"blocks\":%d,"
                        + "\"anchor_bases\":%d,\"anchor_pages\":%d}%n",
                bankMb, blockMb, Files.size(archive),
                compressSeconds, decompressSeconds, exact,
                meta.totalRecords, meta.numBlocks,
                meta.anchorTotalBases, meta.numAnchorPages);
    }
}
