# Measured Benchmarks

Status: local, reproducible measurements. Not an independent third-party audit.
All numbers in the tables are rendered from machine-readable JSON by
`benchmarks/summarize.py`; none are hand-copied.

Two size-focused measurement sets are reported here; the end-to-end query
workload is reported separately in [`END-TO-END.md`](END-TO-END.md).

1. **Frozen held-out corpus** — 10 public datasets, 10 compressor profiles,
   sequence-only input contract, pinned Linux container, three repetitions
   except the explicitly bounded GeCo3 level-9 profile.
2. **100 MB+ pangenome scale test** — a 121.8 MB, 24-record, 120.1 Mbp
   multi-strain *E. coli* RefSeq pangenome, complete multi-FASTA input
   contract, including `--lean` under `-Xmx8m`.

Replication on two additional pangenomes is reported separately in
[`INDEPENDENT.md`](INDEPENDENT.md); the specialized-compressor scale check is
in [`SPECIALIZED.md`](SPECIALIZED.md). Those documents have their own protocols
and are not folded into the tables below.

The corpus is a held-out benchmark conducted in this workspace. The datasets
were selected for size and record structure before these measurements, not for
observed compressibility. They are not a random sample of all genomes, and the
result is not a claim of worldwide superiority.

## Corpus protocol

| Item | Value |
| --- | --- |
| Container | `hypergrass-corpus:locked`, Linux/amd64, pinned by image digest |
| Image digest | `sha256:d6ca62a74e3724172a5c49720e28ab2066eb31b3c965dadfb6220e112dbc3bb3` |
| Tool commits | JARVIS3 `4dcb37b`, GeCo3 `69642a0`; gzip 1.14, xz 5.8.4, OpenJDK 21.0.13 |
| CPU / memory | 2 CPUs, 2 GB container memory |
| Java heap | 1 GB default, 8 MB for `--lean`; `ActiveProcessorCount=1` |
| Operation timeout | 180 s per compression or decompression |
| Repetitions | 3, except `geco3-l9` once per dataset (see deviations) |
| Order | Tool order rotated by dataset and repetition |
| Warmup | None; every operation is a new process |
| Input | Exact `.seq` bytes: FASTA headers and line endings removed, no case or ambiguity changes |
| Verification | SHA-256 and byte length of decompressed output must equal the input; three slice queries per HyperGraSS archive |
| Timing | Child wall/CPU/peak RSS measured inside the container, Docker startup excluded |

Profiles: `hypergrass`, `hypergrass --optimize`, `hypergrass --lean`,
`gzip -9`, `xz -9`, `JARVIS3 -l 7`, `JARVIS3 -l 26`, `GeCo3 -l 9`,
plus supplementary `GeCo3 -l 6` and `GeCo3 -l 3`.

JARVIS3 and GeCo3 are DNA-specialized and run only when every sequence byte is
A/C/G/T. gzip, xz, and HyperGraSS run on every sequence input. Every tool sees
the same bytes for a given dataset.

### Datasets

| ID | Accessions | Sequence bytes | Raw A/C/G/T only |
| --- | --- | ---: | --- |
| `human-mito` | NC_012920.1 | 16,569 | no |
| `fly-mito` | NC_001709.1 | 19,517 | yes |
| `yeast-mito` | NC_001224.1 | 85,779 | yes |
| `plant-plastid` | NC_000932.1 | 154,478 | yes |
| `yeast-chr1` | NC_001133.9 | 230,218 | yes |
| `archaeal-reference` | NC_000909.1 | 1,664,970 | no |
| `ecoli-k12` | NC_000913.3 | 4,641,652 | yes |
| `bsubtilis-168` | NC_000964.3 | 4,215,606 | yes |
| `yeast-three-chr` | NC_001133.9, NC_001134.8, NC_001135.5 | 1,360,022 | yes |
| `ecoli-k12-pair` | NC_000913.3, NC_007779.1 | 9,287,984 | yes |

The two multi-record datasets reuse accessions that also appear alone, so the
ten inputs are not ten statistically independent genomes.

## Deviations and limits

- **The corpus run was stopped early.** It contains 179 of 280 planned results
  and 6 of 10 complete datasets; `ecoli-k12` is partial and `bsubtilis-168`,
  `yeast-three-chr`, and `ecoli-k12-pair` were never measured. The operator
  stopped it deliberately to free the host for the scale test. The report file
  was left in its live `running` state rather than rewritten, and
  `audit.py --allow-incomplete` verifies every present result without
  certifying the missing cells. This partial matrix is not a full benchmark.
- The frozen `geco3-l9` profile is measured **once** per dataset, not three
  times. A 1,048,576-byte pilot exceeded 120 s under a 180 s cap, so repeating
  an operation that cannot finish adds no information. Timeouts and
  container-memory kills are retained in the report and marked, not silently
  dropped.
- `jarvis3-l26` is the other high-effort profile. It is attempted three times,
  but on MB-scale inputs it can exceed the 180 s cap or the 2 GB container
  limit. Those cells are marked `infeasible` and excluded from baseline minima;
  `jarvis3-l7` remains the feasible JARVIS3 baseline on those datasets. No size
  is inferred from a partial compressed file.
- `geco3-l3` and `geco3-l6` are **supplementary** profiles added after the
  level-9 pilot to provide feasible GeCo3 baselines. They were not part of the
  frozen protocol and are labeled `origin=supplementary`.
- The first scale run (`pangenome-scale-1`) exposed a footer record-count bug:
  blocks that start inside a long record were excluded from the total, so a
  24-record FASTA was reported as 22 records and queries for the last two
  records were rejected. Compression, archive size, and round trip were
  unaffected. The fix changes only the stored count and is covered by the
  `Records with Mid-Record Block Starts` suite. The published scale numbers are
  from `pangenome-scale-2`, which uses the fixed build.
- Two scale runs produced byte-identical archives but very different wall
  times: 116.9 s and 155.6 s in the first run versus 45.9 s and 47.5 s in the
  second for the same two HyperGraSS profiles. The host or cache state changed
  between runs. Scale timings are therefore indicative, not precise; sizes and
  round-trip results were stable.
- **The default profile was not changed.** The published default remains the
  4 MiB bank / 2 MiB block configuration, which is 9.8% larger than xz on this
  pangenome. The `--pangenome` profile is opt-in: 16 MiB bank / 8 MiB blocks,
  41.8% smaller than xz. One pangenome, one host, and two runs are not enough
  evidence to replace the default for every input.
- The anchor-bank ablation measures compression and decompression in one
  process, so its peak-RSS values are not comparable with the per-operation
  scale numbers. See [`ABLATION.md`](ABLATION.md).
- The corpus compares **sequence bytes only**. FASTA headers, wrapping, case,
  and line endings are excluded from corpus sizes. The scale test uses complete
  multi-FASTA bytes, so those two size ratios are not directly comparable.
- The measurements come from one Windows host running a Linux container. They
  are not a cross-platform performance study.
- GeCo3 level 9 is included for completeness but is not a practical baseline at
  this timeout; its failures describe the configured cap, not the tool's
  best achievable size.

## Audit gates

1. **Automated report audit** — `benchmarks/corpus/audit.py` re-hashes every
   artifact and restored file on disk, checks the result matrix against the
   plan, recomputes every summary row from raw results, verifies round trips and
   query slices, and confirms the protocol and image digests. On a finished run
   it also requires the full planned matrix; `--allow-incomplete` audits the
   completed cells of a deliberately stopped run without certifying coverage.
2. **Lossless verification gate** — `test.ps1` / `test.sh` runs the 14
   edge-case suites, donor-index suites, anchor-codec suites, byte-accounting
   suites, and donor-cost suites. The published run passed all of them.

## Corpus results

Run status: `running` — the run was stopped before the runner wrote a final status. It contains 179 of 280 planned results and 6 of 10 complete datasets, with 3 infeasible high-effort cells. This is a partial matrix, not a finished benchmark. The report JSON is committed at `benchmarks/corpus/runs/corpus-full-1/report.json`; inputs and artifacts are not.

> Coverage: 179 of 280 planned results, 6 of 10 datasets complete. The report status is `running` because the runner was terminated before writing a final status. The tables below are the completed evidence only; they are not a full 10-dataset matrix.

| Cell with fewer than planned repetitions | Results | Planned |
| --- | ---: | ---: |
| ecoli-k12 / hypergrass | 1 | 3 |
| ecoli-k12 / hypergrass-optimize | 1 | 3 |
| ecoli-k12 / hypergrass-lean | 1 | 3 |
| ecoli-k12 / jarvis3-l7 | 1 | 3 |
| ecoli-k12 / jarvis3-l26 | 1 | 3 |
| ecoli-k12 / gzip-9 | 1 | 3 |
| ecoli-k12 / xz-9 | 1 | 3 |
| ecoli-k12 / geco3-l6 | 2 | 3 |
| ecoli-k12 / geco3-l3 | 1 | 3 |

### Head-to-head, frozen corpus

| Dataset | Input bytes | HyperGraSS | `--optimize` | `--lean` | Best completed baseline | Its bytes | Δ smaller than best (+ = smaller) |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| human-mito | 16,569 | 3,959 | 3,959 | 3,975 | xz-9 | 4,800 | +17.52% |
| fly-mito | 19,517 | 3,567 | 3,567 | 3,673 | jarvis3-l26 | 3,739 | +4.60% |
| yeast-mito | 85,779 | 15,338 | 15,338 | 15,733 | geco3-l9 | 15,425 | +0.56% |
| plant-plastid | 154,478 | 30,419 | 30,419 | 36,661 | geco3-l9 | 30,839 | +1.36% |
| yeast-chr1 | 230,218 | 51,864 | 51,864 | 54,521 | jarvis3-l26 | 52,331 | +0.89% |
| archaeal-reference | 1,664,970 | 370,196 | 370,196 | 376,881 | xz-9 | 410,584 | +9.84% |
| ecoli-k12 | 4,641,652 | 1,087,650 | 1,087,650 | 1,118,253 | geco3-l6 | 1,094,264 | +0.60% |

### Effective bits per base, frozen corpus

| Dataset | hypergrass | hypergrass-optimize | hypergrass-lean | gzip-9 | xz-9 | jarvis3-l7 | jarvis3-l26 | geco3-l9 | geco3-l6 | geco3-l3 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| human-mito | 1.912 | 1.912 | 1.919 | 2.436 | 2.318 | n/a | n/a | n/a | n/a | n/a |
| fly-mito | 1.462 | 1.462 | 1.506 | 1.927 | 1.828 | 1.767 | 1.533 | 1.549 | 1.560 | 1.565 |
| yeast-mito | 1.430 | 1.430 | 1.467 | 1.836 | 1.685 | 1.819 | 1.460 | 1.439 | 1.444 | 1.447 |
| plant-plastid | 1.575 | 1.575 | 1.899 | 2.283 | 2.154 | 1.665 | 1.597 | 1.597 | 1.600 | 1.601 |
| yeast-chr1 | 1.802 | 1.802 | 1.895 | 2.250 | 2.122 | 1.868 | 1.818 | infeasible | 1.823 | 1.824 |
| archaeal-reference | 1.779 | 1.779 | 1.811 | 2.150 | 1.973 | n/a | n/a | n/a | n/a | n/a |
| ecoli-k12 | 1.875 | 1.875 | 1.927 | 2.240 | 2.044 | 1.987 | infeasible | infeasible | 1.886 | 1.898 |

`infeasible` marks the high-effort `geco3-l9` and `jarvis3-l26` cells that exceeded the 180 s operation cap (compression or decompression) or the container memory limit; they are not treated as zero-byte results.

### Median of per-dataset medians: compression time and peak RSS

| Profile | Median wall (s) | Median CPU (s) | Median peak RSS (MiB) | Successful datasets |
| --- | ---: | ---: | ---: | ---: |
| hypergrass | 0.94 | 1.28 | 68.3 | 7 |
| hypergrass-optimize | 0.84 | 1.10 | 69.3 | 7 |
| hypergrass-lean | 0.59 | 0.70 | 47.0 | 7 |
| gzip-9 | 0.21 | 0.14 | 9.6 | 7 |
| xz-9 | 0.12 | 0.05 | 9.5 | 7 |
| jarvis3-l7 | 0.39 | 0.22 | 193.4 | 5 |
| jarvis3-l26 | 0.94 | 0.75 | 495.9 | 4 |
| geco3-l9 | 24.95 | 14.67 | 2042.2 | 3 |
| geco3-l6 | 0.76 | 0.62 | 348.0 | 5 |
| geco3-l3 | 0.49 | 0.32 | 123.8 | 5 |

Cells with fewer than the planned repetitions contribute their available samples; the coverage note above lists them.

### Non-success cells

70 dataset/profile cells: 57 successful, 10 not applicable, 3 infeasible or failed.

| Dataset | Profile | Status | Attempts | Successful |
| --- | --- | --- | ---: | ---: |
| yeast-chr1 | geco3-l9 | failed | 1 | 0 |
| ecoli-k12 | jarvis3-l26 | timeout | 1 | 0 |
| ecoli-k12 | geco3-l9 | timeout | 1 | 0 |


## 100 MB+ pangenome scale test

The pangenome is built by `benchmarks/pangenome/fetch.py` from NCBI Nucleotide
RefSeq complete *E. coli* records. The fetcher records accession versions,
lengths, and SHA-256 in `benchmarks/pangenome/manifest.json`; the FASTA itself
is not committed. The selected records are 24 distinct accessions totaling
120,108,223 bases (121,825,884 FASTA bytes). The runner verifies the manifest
hash before measuring, compresses the complete multi-FASTA, decompresses it,
compares SHA-256, and queries slices from the first, middle, and last records.
The `pangenome-scale-2` report JSON is committed at
`benchmarks/pangenome/runs/pangenome-scale-2/report.json` and the tuned
`pangenome-scale-3` report at
`benchmarks/pangenome/runs/pangenome-scale-3/report.json`; the FASTA and
archives are not. The opt-in `--pangenome` profile (16 MiB anchor bank, 8 MiB
blocks) was measured in `pangenome-scale-3`; the table below merges both
reports after checking that their input manifests are identical.

<!-- SCALE-TABLES:START -->
### Scale results

Input: `data/pangenome.fasta`, 121,825,884 bytes, 24 records, 120,108,223 bases. SHA-256 `f03f0dd10c541b99ddf63d18053e3a9c6009910a81700ba18f15e24dcb7d2285`.

| Profile | Archive bytes | Ratio | Bits/byte | Compress wall (s) | Compress CPU (s) | Peak RSS (MiB) | Round trip | Queries |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| hypergrass | 9,215,764 | 13.22x | 0.605 | 45.9 | 80.0 | 1083.0 | yes | 3/3 exact |
| hypergrass-pangenome | 4,885,851 | 24.93x | 0.321 | 30.6 | 47.1 | 1007.9 | yes | 3/3 exact |
| hypergrass-lean | 28,910,764 | 4.21x | 1.898 | 47.5 | 44.4 | 52.0 | yes | 3/3 exact |
| gzip-9 | 35,859,298 | 3.40x | 2.355 | 93.8 | 89.8 | 9.5 | yes | — |
| xz-9 | 8,394,900 | 14.51x | 0.551 | 186.3 | 186.6 | 737.9 | yes | — |

<!-- SCALE-TABLES:END -->

A second run with byte-identical archives measured 116.9 s and 155.6 s for the
two default HyperGraSS profiles instead of 45.9 s and 47.5 s. Host or cache
state changed between runs, so compression times here are indicative, not
precise; archive sizes and round trips were stable. The tuned profile was also
measured twice: 30.7 s in the ablation and 30.6 s in the scale harness, with a
byte-identical 4,885,851-byte archive both times. Its decompression measured
16.4 s at 419 MiB peak RSS, against 7.1 s at 189 MiB for xz. The xz archive is
a single stream with no index: it cannot answer record queries without
decompression, so its query cell is `—`, not a failed check.

## Interpretation

The results support a bounded claim, not a superiority claim.

**On the completed sequence-only corpus cells, the default archive was the
smallest of every successful profile on 7 of 7 datasets.** The margin over the
best completed baseline has a median of **+1.36%** and a range of
**+0.56% to +17.52%**. Four of the seven margins are below 1.4% (yeast
mitochondria +0.56%, chloroplast +1.36%, yeast chromosome I +0.89%, *E. coli*
+0.60%). Those are narrow wins against strong specialized compressors, not a
decisive gap. Only the 16.6 kB mitochondrial case (+17.52% versus xz) and the
1.66 MB archaeal genome (+9.84% versus xz) exceed 4%.

**On the 121.8 MB pangenome, the default profile loses to xz -9.** xz produced
8,394,900 bytes against 9,215,764 bytes for default HyperGraSS, a 9.8% larger
archive for HyperGraSS. The default was faster to compress (45.9 s versus
186.3 s) but used more memory (1,083 MiB versus 738 MiB). Any claim that the
default configuration supersedes general-purpose compressors is false on this
input.

**The opt-in `--pangenome` profile wins on size.** With a 16 MiB anchor bank and
8 MiB blocks, HyperGraSS produced 4,885,851 bytes: 41.8% smaller than xz -9 and
47.0% smaller than the default, at 24.93x ratio and 0.321 bits per input byte.
It compressed in 30.6 s versus 186.3 s for xz and answered all three record
queries exactly. It did not win every axis: xz decompressed in 7.1 s versus
16.4 s for HyperGraSS, and used 189 MiB versus 419 MiB for that operation. The
archive was byte-identical across two runs, so the size result is stable even
though timings are not.

This is a measured tuning win on one pangenome, not a new compression law. The
ablation shows the improvement comes from removing the 4 MiB anchor-bank cap,
not from a universal block-size rule; see [`ABLATION.md`](ABLATION.md).

**The query advantage is workload-dependent.** On 13 exact 2,000-base slices,
the tuned profile answered in a median of 1.35 s against 1.36 s for xz: a tie,
not a speedup. xz restored the archive in 7.1 s versus 16.4 s for the tuned
profile, so a small query batch finishes sooner with xz. The tuned profile
needs about 859 queries per restore to offset that difference. The default
profile was much slower at queries (8.0 s median) and is dominated by the tuned
profile on all four measured axes. Details are in
[`END-TO-END.md`](END-TO-END.md).

**`--lean` is a memory result, not a size result.** It round-tripped the
121.8 MB FASTA and answered all three queries with a 52 MiB peak RSS, but its
archive is 28,910,764 bytes: 3.14x the default archive, 3.44x the xz archive,
and 19% smaller than gzip. It is the option for a hard memory bound, at a large
size cost.

**Time.** On the corpus, the default profile has a median compression time of
0.94 s versus 0.21 s for gzip and 0.12 s for xz, so it is roughly 4-8x slower.
At the high-effort end, `geco3-l9` has a 24.95 s median and is infeasible on
larger inputs; `jarvis3-l26` and `geco3-l9` both exceed the 180 s cap on the
4.6 MB *E. coli* genome. HyperGraSS finishes those cases, but that is a
practicality result, not a size win against configurations that completed.

**Measured cause of the default's scale loss.** The default anchor bank is
capped at 4,194,304 bases (64 pages), while each pangenome record is about
4.7-5.4 Mbp. The bank cannot represent even one full record, let alone
variation across 24 strains, and the default encoder uses 2 MB blocks at this
input size. xz -9 has a 64 MiB dictionary and sees a contiguous 121.8 MB
stream, so it can reference cross-strain repeats that the default anchor bank
cannot. The ablation tested that explanation directly: at the same 2 MiB block
size, raising the bank to 8 MiB reduced the archive from 9,215,764 to
5,173,803 bytes (43.9% smaller); 16 MiB with 8 MiB blocks reached 4,885,851
bytes. Block size alone was not the lever: at a 4 MiB bank, 4 MiB blocks made
the archive larger. The 16 MiB cap is not fully consumed either; the winning
archive holds 144 anchor pages, about 9.4 Mbp. The cap was a real ceiling, but
a larger cap does not automatically add anchors.

**What the system does provide.** A pure-Java, zero-dependency implementation
with a documented archive format, byte-exact reconstruction, record and range
queries on the archive, and a bounded-memory mode. Those properties are not
measured by archive size alone, and xz has no equivalent query interface.

The honest summary is: competitive at small and medium sequence-only scales,
narrowly first on the completed corpus cells; **not** the smallest at 100 MB+
pangenome scale with the default profile, but **41.8% smaller than xz** on that
input with the opt-in `--pangenome` profile. It is differentiated mainly by
native record and range queries, a documented archive format, and a pure-Java,
zero-dependency implementation, not by universal compression superiority.

## Reproduce

```sh
# Corpus (requires the locked image and prepared data)
python -B benchmarks/corpus/prepare.py --verify
python -B benchmarks/corpus/run.py --run-id corpus-repro-1
python -B benchmarks/corpus/audit.py benchmarks/corpus/runs/corpus-repro-1
python -B benchmarks/summarize.py \
  --corpus benchmarks/corpus/runs/corpus-repro-1/report.json

# Scale test (downloads public data, then builds a 121.8 MB FASTA)
python -B benchmarks/pangenome/fetch.py
python -B benchmarks/pangenome/run.py --run-id pangenome-scale-repro
python -B benchmarks/summarize.py \
  --corpus benchmarks/corpus/runs/corpus-repro-1/report.json \
  --scale benchmarks/pangenome/runs/pangenome-scale-repro/report.json
```

The corpus runner refuses an existing run directory. Use a new `--run-id` or
delete the old run directory deliberately. Downloaded data and work artifacts
are gitignored; the published `report.json` evidence files are committed, and
individual reports are force-added when they live under a gitignored `runs/`
directory. If a run is stopped deliberately, add `--allow-incomplete` to the
audit command to verify the completed cells without certifying coverage.
