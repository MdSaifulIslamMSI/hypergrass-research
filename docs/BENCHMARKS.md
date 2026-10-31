# Measured Benchmarks

Status: local, reproducible measurements. Not an independent third-party audit.
All numbers in the tables are rendered from machine-readable JSON by
`benchmarks/summarize.py`; none are hand-copied.

Two measurement sets are reported:

1. **Frozen held-out corpus** — 10 public datasets, 10 compressor profiles,
   sequence-only input contract, pinned Linux container, three repetitions
   except the explicitly bounded GeCo3 level-9 profile.
2. **100 MB+ pangenome scale test** — a 121.8 MB, 24-record, 120.1 Mbp
   multi-strain *E. coli* RefSeq pangenome, complete multi-FASTA input
   contract, including `--lean` under `-Xmx8m`.

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

- The frozen `geco3-l9` profile is measured **once** per dataset, not three
  times. A 1,048,576-byte pilot exceeded 120 s under a 180 s cap, so repeating
  an operation that cannot finish adds no information. Timeouts and
  container-memory kills are retained in the report and marked, not silently
  dropped.
- `geco3-l3` and `geco3-l6` are **supplementary** profiles added after the
  level-9 pilot to provide feasible GeCo3 baselines. They were not part of the
  frozen protocol and are labeled `origin=supplementary`.
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
   artifact and restored file on disk, checks the expected result matrix,
   recomputes every summary row from raw results, verifies round trips and
   query slices, and confirms the protocol and image digests.
2. **Lossless verification gate** — `test.ps1` / `test.sh` runs the 13
   edge-case suites, donor-index suites, anchor-codec suites, byte-accounting
   suites, and donor-cost suites. The published run passed all of them.

## Corpus results

Run status: <!-- RUN-STATUS -->

<!-- CORPUS-TABLES -->

## 100 MB+ pangenome scale test

The pangenome is built by `benchmarks/pangenome/fetch.py` from NCBI Nucleotide
RefSeq complete *E. coli* records. The fetcher records accession versions,
lengths, and SHA-256 in `benchmarks/pangenome/manifest.json`; the FASTA itself
is not committed. The selected records are 24 distinct accessions totaling
120,108,223 bases (121,825,884 FASTA bytes). The runner verifies the manifest
hash before measuring, compresses the complete multi-FASTA, decompresses it,
compares SHA-256, and queries slices from the first, middle, and last records.

<!-- SCALE-TABLES -->

## Interpretation

<!-- INTERPRETATION -->

## Reproduce

```sh
# Corpus (requires the locked image and prepared data)
python -B benchmarks/corpus/prepare.py --verify
python -B benchmarks/corpus/run.py --run-id corpus-full-1
python -B benchmarks/corpus/audit.py benchmarks/corpus/runs/corpus-full-1
python -B benchmarks/summarize.py \
  --corpus benchmarks/corpus/runs/corpus-full-1/report.json

# Scale test (downloads public data, then builds a 121.8 MB FASTA)
python -B benchmarks/pangenome/fetch.py
python -B benchmarks/pangenome/run.py --run-id pangenome-scale-1
python -B benchmarks/summarize.py \
  --corpus benchmarks/corpus/runs/corpus-full-1/report.json \
  --scale benchmarks/pangenome/runs/pangenome-scale-1/report.json
```

The corpus runner refuses an existing run directory. Use a new `--run-id` or
delete the old run directory deliberately. `runs/` and downloaded data are
gitignored.
