# End-to-end query workload

This benchmark asks a narrower question than “which compressor is smallest?”:

> For a stored 121.8 MB pangenome, how long does it take to serve many small
> record slices, and what does one full restore cost?

It uses the archives and input already verified by `pangenome-scale-2` and
`pangenome-scale-3`. No archive is rebuilt, so the size and round-trip evidence
stays the same.

## Workload

- Input: 24-record, 121,825,884-byte *E. coli* pangenome.
- Queries: 13 deterministic 2,000-base slices, one every second record plus the
  last record. Offsets are 37% into each record, which spreads reads across the
  whole archive.
- Each query runs in a fresh process inside the locked Linux image.
- HyperGraSS answers with its page-indexed `query` command.
- gzip and xz decompress the stream up to the requested record and stop once
  the slice is complete. A late-record query therefore pays for almost the
  whole archive.
- Every returned slice is compared byte-for-byte with the FASTA source.
- One full decompression per profile is taken from the already-published
  `pangenome-scale-2` or `pangenome-scale-3` report.

## Metrics

- archive bytes,
- per-query wall time: median, p95, maximum,
- query peak RSS,
- one full-restore wall time and peak RSS,
- exact query count,
- external dependencies.

The analysis computes query speedups, the number of queries needed for one
restore plus Q queries to favour HyperGraSS on wall time, and the Pareto
frontier over archive size, query latency, restore time, and query memory.

The break-even calculation ignores storage price and concurrency. It is a
query-intensity threshold, not a total cost of ownership model.

## Run

```sh
python -B benchmarks/endtoend/run.py --run-id endtoend-1
```

The report is written to `benchmarks/endtoend/runs/<run-id>/report.json`. Run
directories are gitignored; the published `endtoend-1` report is force-added as
evidence, while the raw work directories are not.

## Limits

This is one pangenome, one host, one repetition, and a single-query process
model. Batch queries inside one JVM would reduce HyperGraSS startup overhead and
are not measured here. gzip/xz could be combined with a separate index or
cached decompression; that is a different system and is not what this benchmark
measures. Timings on this host varied by up to 2.5x between identical scale
runs, so treat latency as indicative.
