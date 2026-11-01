# Large query workload

## Question

The 13-query workload found the tuned `--pangenome` profile tied with `xz -9` on
per-query latency and ahead only on storage. Thirteen queries is a small sample,
and it queried one offset in every other record. Does the result change when the
workload is larger and spreads queries across every record?

This run serves 96 deterministic 2,000-base slices: all 24 records, four offsets
per record. It compares the tuned profile with sequential gzip and xz
extraction, using the same archives as the scale test.

## Protocol

- Input: the same 24-record, 121,825,884-byte *E. coli* pangenome used by the
  scale test, SHA-256
  `f03f0dd10c541b99ddf63d18053e3a9c6009910a81700ba18f15e24dcb7d2285`.
- Queries: all 24 records, at 20%, 40%, 60%, and 80% of each record length,
  2,000 bases each (96 slices total). The slice length is reduced only if a
  record is shorter than 2,000 bases.
- Profiles: `hypergrass-pangenome`, `xz-9`, and `gzip-9`. The default and
  `--lean` profiles are not repeated; the 13-query workload already shows the
  default is dominated and `--lean` trades a much larger archive for latency.
- Each query runs in the locked Linux container with 2 CPUs and the profile's
  memory limit, in a fresh process. HyperGraSS answers from its block/page
  index; gzip and xz decompress sequentially until the slice is complete.
- Archives are the exact artifacts from the scale reports, checked by size and
  SHA-256 before querying. Full-restore wall time and peak RSS come from those
  reports.
- Output bytes are compared with the source FASTA. Process wall time, child
  CPU time, and peak RSS are measured inside the container. JVM or Python
  startup is included in every query.

<!-- E2E-TABLES:START -->
<!-- E2E-TABLES:END -->

## Interpretation

_(Filled from the measured report.)_

## Limits

- One pangenome, one host, one run, 96 queries, and four fixed offsets per
  record. No repeated-query cache, concurrency, or mixed read/write workload is
  measured.
- Query latency varies strongly by record and offset; the median is not
  evidence of constant-time access.
- JVM or Python startup is included in every query, so absolute latencies are
  higher than the in-process work.
- Full-restore numbers are reused from the scale reports rather than remeasured
  in this run; archive hashes were checked before reuse.
- The result is a workload study on one input, not a claim that indexed access
  is always faster or slower than sequential decompression.

## Reproduce

```sh
python -B benchmarks/endtoend/run.py --run-id endtoend-96 \
  --queries-per-record 4 --profiles hypergrass-pangenome,xz-9,gzip-9
python -B benchmarks/endtoend/render.py \
  --report benchmarks/endtoend/runs/endtoend-96/report.json \
  --doc docs/END-TO-END-LARGE.md
```

Reports are written under `benchmarks/endtoend/runs/`, which is gitignored;
the published report is force-added as evidence.
