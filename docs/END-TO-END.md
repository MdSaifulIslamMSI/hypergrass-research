# End-to-end query workload

## Question

The archive-size benchmark answers how small a stored pangenome is. It does not
answer whether the index is useful for a query-heavy workload. This benchmark
serves 13 deterministic 2,000-base slices from the 121.8 MB pangenome, each in
a fresh process, and compares HyperGraSS with sequential gzip and xz
extraction.

## Protocol

- Input: the same 24-record, 121,825,884-byte *E. coli* pangenome used by the
  scale test, SHA-256
  `f03f0dd10c541b99ddf63d18053e3a9c6009910a81700ba18f15e24dcb7d2285`.
- Queries: every other record plus the last (records 0, 2, 4, ..., 22, 23),
  each 2,000 bases starting at 37% of the record length. The slice length is
  reduced only if a record is shorter than 2,000 bases.
- Each query runs in the locked Linux container with 2 CPUs and the profile's
  memory limit, in a fresh process. The 900 s timeout is never approached.
- HyperGraSS answers with `query` from its block/page index. gzip and xz use
  `benchmarks/endtoend/extract.py`, which decompresses sequentially and stops
  once the requested slice is complete; that is the honest cost of a
  non-indexed stream.
- Every profile's archive is the exact artifact from its scale report, checked
  by size and SHA-256 before querying. Full-restore wall time and peak RSS are
  taken from those scale reports (one full decompress per profile), not
  remeasured here.
- Output bytes are compared with the source FASTA. Process wall time, child
  CPU time, and peak RSS are measured inside the container. JVM startup is
  included for HyperGraSS; Python startup is included for gzip and xz.

<!-- E2E-TABLES:START -->
### Workload

| Item | Value |
| --- | --- |
| Input | `benchmarks/pangenome/data/pangenome.fasta`, 121,825,884 bytes, 24 records |
| Query | 13 slices of 2,000 bases, one fresh process each |
| Records queried | 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 23 |
| Restore | one full decompress per profile, from the scale report |

### Results

| Profile | Archive bytes | Query median (s) | Query p95 (s) | Query max (s) | Query total (s) | Query peak RSS (MiB) | Restore (s) | Restore peak RSS (MiB) | Exact |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| hypergrass | 9,215,764 | 8.02 | 12.54 | 12.54 | 101.1 | 130 | 43.0 | 491 | 13/13 |
| hypergrass-pangenome | 4,885,851 | 1.35 | 2.16 | 2.16 | 18.2 | 74 | 16.4 | 419 | 13/13 |
| hypergrass-lean | 28,910,764 | 0.73 | 1.19 | 1.19 | 10.3 | 44 | 43.4 | 48 | 13/13 |
| gzip-9 | 35,859,298 | 1.18 | 2.37 | 2.37 | 16.1 | 10 | 2.8 | 10 | 13/13 |
| xz-9 | 8,394,900 | 1.36 | 2.23 | 2.23 | 17.7 | 189 | 7.1 | 189 | 13/13 |

With 13 samples, the p95 column equals the maximum. Query time varies strongly by record and profile; the median is not a constant-time claim.

### Query speedup

Baseline query median divided by HyperGraSS query median. Above 1.0 means HyperGraSS answered the same slice faster.

| Comparison | Speedup |
| --- | ---: |
| hypergrass vs gzip-9 | 0.15x |
| hypergrass vs xz-9 | 0.17x |
| hypergrass-pangenome vs gzip-9 | 0.87x |
| hypergrass-pangenome vs xz-9 | 1.01x |
| hypergrass-lean vs gzip-9 | 1.62x |
| hypergrass-lean vs xz-9 | 1.87x |

### Break-even query count

Queries needed for `one restore + Q queries` to favour HyperGraSS on wall time: `ceil((HGS restore - baseline restore) / (baseline query - HGS query))`. `never` means the baseline has no per-query advantage to recover the restore cost.

| Comparison | Break-even queries |
| --- | ---: |
| hypergrass vs gzip-9 | never |
| hypergrass vs xz-9 | never |
| hypergrass-pangenome vs gzip-9 | never |
| hypergrass-pangenome vs xz-9 | 859 |
| hypergrass-lean vs gzip-9 | 91 |
| hypergrass-lean vs xz-9 | 58 |

### Pareto frontier

Non-dominated profiles across archive bytes, query median, restore time, query peak RSS: **`hypergrass-pangenome`, `hypergrass-lean`, `gzip-9`, `xz-9`**.

<!-- E2E-TABLES:END -->

## Interpretation

**The index is not automatically faster than sequential extraction.** On this
workload the tuned `--pangenome` profile answered a slice in a median of 1.35 s,
against 1.36 s for xz and 1.18 s for gzip. That is a tie with xz, not a query
win. xz starts from an 8.39 MB archive (72% more bytes) but restores in 7.1 s
instead of 16.4 s, so a workload with one restore and a handful of queries
finishes sooner with xz. The break-even calculation says the tuned profile
needs about 859 queries per restore to recover its slower restore through its
negligible per-query edge.

**The tuned profile's measured advantage is storage, not query speed.** It has
the smallest archive (4,885,851 bytes, 41.8% below xz), but it also has a
slower restore than both baselines. That is the right trade for a store-once,
read-rarely pangenome archive, and the wrong trade for a workload dominated by
full restores.

**The default profile is poor for queries at this scale.** Its median slice
took 8.0 s and its total for 13 slices was 101 s, 5.6x the tuned total. It is
dominated by the tuned profile on all four Pareto axes measured here: archive
bytes, query median, restore wall time, and query peak RSS. The default
configuration should not be used for query-heavy pangenome access.

**`--lean` is the query-speed option, at a large size cost.** It answered
slices fastest (0.73 s median, 1.6x gzip and 1.9x xz) with 44 MiB query RSS,
but its archive is 28,910,764 bytes and its restore is 43.4 s. It needs roughly
58 queries against xz, or 91 against gzip, to recover the restore cost. That is
plausible for a long-lived server that restores once and answers many queries,
not for a short batch.

**Pareto interpretation.** The non-dominated set contains four profiles, each
for a different reason: `hypergrass-pangenome` for archive size,
`hypergrass-lean` for query latency, `gzip-9` for restore time and query
memory, and `xz-9` for the best combination of restore time and moderate
archive size. No single profile wins all four axes. The default HyperGraSS
profile is not on the frontier.

The honest summary is: HyperGraSS's indexed query matches xz's query latency
on this input while storing 42% fewer bytes, but it does not beat xz
end-to-end for a small number of queries because its restore is slower. The
query interface is a real capability, not a universal speedup.

## Limits

- One pangenome, one host, one run, 13 queries, and one offset per queried
  record. No repeated-query cache, concurrency, or mixed read/write workload
  is measured.
- With 13 samples, the reported p95 equals the maximum. Query latency varies
  strongly by record; the median is not evidence of constant-time access.
- JVM or Python startup is included in every query, so absolute latencies are
  higher than the in-process work. The comparison is like-for-like only in the
  sense that every profile pays one fresh process per query.
- Full-restore numbers are reused from the scale reports rather than remeasured
  in this run; archive hashes were checked before reuse.
- The result is a workload study on one input, not a claim that indexed access
  is always faster or slower than sequential decompression.

## Reproduce

```sh
python -B benchmarks/endtoend/run.py --run-id endtoend-repro
python -B benchmarks/endtoend/render.py \
  --report benchmarks/endtoend/runs/endtoend-repro/report.json
```

Reports are written under `benchmarks/endtoend/runs/`, which is gitignored;
the published `endtoend-1` report is force-added as evidence.
