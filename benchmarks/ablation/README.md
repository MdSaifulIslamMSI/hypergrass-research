# Anchor-bank and block-size ablation

This benchmark tests the hypothesis recorded in [`docs/BENCHMARKS.md`](../../docs/BENCHMARKS.md):
the default 4 MiB anchor bank cannot represent a full ~5 Mbp pangenome record,
and 2 MiB blocks cannot exploit the same long-range matches as xz's 64 MiB
dictionary. The question is whether a larger bank, larger blocks, or both close
the 8.9% size gap to `xz -9` on the 121.8 MB pangenome.

## Protocol

- Input: the same 24-record, 121,825,884-byte FASTA used by
  `pangenome-scale-2`, verified against its manifest SHA-256.
- Configurations: anchor-bank caps of 4, 8, and 16 MiB crossed with block sizes
  of 2, 4, and 8 MiB.
- Every configuration runs in a fresh JVM inside the locked Linux image with
  2 CPUs and a 2 GiB container limit.
- Each run compresses, decompresses, and checks the restored file byte-for-byte
  with `Files.mismatch`.
- Archive size, internal compression/decompression time, process wall/CPU time,
  and peak RSS are recorded.
- The `xz -9` baseline is the already-published 8,394,900-byte archive from
  `pangenome-scale-2`; xz is not rerun, because the input and tool version are
  unchanged.

This is a controlled ablation, not a parameter sweep. It has one input and one
host; every configuration was measured once, and the two best configurations
were repeated with byte-identical archives. A configuration that wins here
still needs a second pangenome before it can become a default.

## Run

```sh
# Compile the harness (bin/ is not committed)
javac -cp bin -d bin benchmarks/ablation/AnchorBankAblation.java

python -B benchmarks/ablation/run.py --run-id anchor-ablation-1
```

The report is written to `benchmarks/ablation/runs/<run-id>/report.json`.
Run directories are gitignored. The two published `report.json` files are
force-added as evidence; raw work directories and archives are not.

## Reading the result

- `archive_bytes` — final `.hgs` size.
- `delta_vs_xz_percent` — positive means larger than xz; negative means smaller.
- `beats_xz` — true only when the archive is strictly smaller.
- `exact` — byte-exact round trip.
- `peak_rss_bytes` — peak child RSS for the whole configuration process.
- `wall_seconds` — container-measured process wall time, including JVM startup.

A larger bank or block that reduces size but exceeds the 2 GiB container limit
is an infeasible configuration, not a win.
