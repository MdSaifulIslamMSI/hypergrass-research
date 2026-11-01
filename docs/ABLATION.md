# Anchor-bank and block-size ablation

## Question

The published 121.8 MB scale result left an explicit hypothesis: the default
4 MiB anchor bank cannot represent a full ~5 Mbp record, and 2 MiB blocks cannot
use the same long-range matches as xz's 64 MiB dictionary. If that ceiling is
real, increasing the bank or block size should reduce archive size. If it is
not, the 8.9% loss to xz needs a different explanation.

This ablation measures that question directly.

## Protocol

- Input: the 24-record, 121,825,884-byte *E. coli* pangenome from
  `pangenome-scale-2`, verified against manifest SHA-256
  `f03f0dd1...d2285`.
- Configurations: anchor-bank caps of 4, 8, and 16 MiB crossed with block sizes
  of 2, 4, and 8 MiB.
- Every configuration runs in a fresh JVM inside the locked Linux image with
  2 CPUs and a 2 GiB container limit.
- Each run compresses, decompresses, and checks the restored file byte-for-byte
  with `Files.mismatch`.
- Archive size, internal compression/decompression time, process wall/CPU time,
  and peak RSS are recorded.
- The `xz -9` baseline is the published 8,394,900-byte archive from
  `pangenome-scale-2`. xz was not rerun because the input and tool version are
  unchanged.

`Delta vs xz` is `(archive - 8,394,900) / 8,394,900`. Negative means smaller
than xz. `Delta vs default` uses the published 9,215,764-byte default archive.
Peak RSS covers both compression and decompression in the same process; the
scale report measures each operation separately, so these RSS numbers are not
directly interchangeable.

## Run 1: full grid

| Bank | Block | Archive | Ratio | Bits/byte | Δ vs xz | Δ vs default | Compress | Decompress | Wall | Peak RSS | Exact |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 4 MiB | 2 MiB | 9,215,764 B | 13.22x | 0.605 | +9.78% | 0.00% | 46.2 s | 40.1 s | 99.3 s | 1,263 MiB | yes |
| 4 MiB | 4 MiB | 9,535,419 B | 12.78x | 0.626 | +13.59% | +3.47% | 47.3 s | 79.9 s | 200.5 s | 1,162 MiB | yes |
| 4 MiB | 8 MiB | 7,669,270 B | 15.88x | 0.504 | −8.64% | −16.78% | 128.6 s | 84.0 s | 230.2 s | 1,275 MiB | yes |
| 8 MiB | 2 MiB | 5,173,803 B | 23.55x | 0.340 | −38.37% | −43.86% | 32.3 s | 13.2 s | 61.3 s | 1,179 MiB | yes |
| 8 MiB | 4 MiB | 5,173,643 B | 23.55x | 0.340 | −38.37% | −43.86% | 29.3 s | 13.3 s | 55.6 s | 1,192 MiB | yes |
| 8 MiB | 8 MiB | 4,974,040 B | 24.49x | 0.327 | −40.75% | −46.03% | 46.1 s | 19.5 s | 83.9 s | 1,201 MiB | yes |
| 16 MiB | 2 MiB | 5,061,060 B | 24.07x | 0.332 | −39.71% | −45.08% | 26.2 s | 17.9 s | 56.1 s | 1,200 MiB | yes |
| 16 MiB | 4 MiB | 5,055,576 B | 24.10x | 0.332 | −39.78% | −45.14% | 35.7 s | 9.8 s | 64.2 s | 1,205 MiB | yes |
| **16 MiB** | **8 MiB** | **4,885,851 B** | **24.93x** | **0.321** | **−41.80%** | **−46.98%** | **30.7 s** | **10.3 s** | **50.8 s** | **1,188 MiB** | **yes** |

## Run 2: repeat of the two best configurations

| Bank | Block | Archive | Δ vs xz | Compress | Decompress | Wall | Peak RSS | Exact |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 8 MiB | 8 MiB | 4,974,040 B | −40.75% | 37.8 s | 18.5 s | 72.0 s | 1,188 MiB | yes |
| **16 MiB** | **8 MiB** | **4,885,851 B** | **−41.80%** | **31.2 s** | **11.0 s** | **64.1 s** | **1,197 MiB** | **yes** |

Both archives are byte-identical to run 1. Wall times differ by 15–30%, which
is consistent with the timing variance already documented for the scale test.

## Run 3: 32 MiB bank check

The 16 MiB bank stored only 9.43 Mbp of anchors, so a larger cap should not
change the archive unless selection was still hitting the cap. This run tested
that prediction directly.

| Bank | Block | Archive | Δ vs xz | Δ vs default | Compress | Decompress | Wall | Peak RSS | Exact |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 32 MiB | 8 MiB | 4,885,851 B | −41.80% | −46.98% | 29.6 s | 10.2 s | 55.8 s | 1,194 MiB | yes |

The archive is byte-identical to the 16 MiB / 8 MiB result: the same 12 blocks,
144 anchor pages, and 9,434,480 anchor bases. Doubling the cap bought nothing
and did not reduce RSS. This is a negative result for the "more bank is better"
hypothesis and confirms that anchor selection, not the cap, is the binding
constraint on this input.

## Findings

1. **The 4 MiB bank was the ceiling.** Moving to 8 MiB reduced the archive from
   9,215,764 to 5,173,803 bytes at the same 2 MiB block size: **43.9% smaller**
   and **38.4% smaller than xz**.
2. **Block size alone is not the lever.** At a 4 MiB bank, 4 MiB blocks made the
   archive larger; 8 MiB blocks helped, but remained worse than any 8 MiB bank
   configuration.
3. **Bank and block size interact.** With a 16 MiB cap, 8 MiB blocks reached
   4,885,851 bytes, **41.8% smaller than xz** and **47.0% smaller than the old
   default**.
4. **The 16 MiB cap is not fully consumed.** The winning archive contains 144
   anchor pages, about 9.4 Mbp of anchors, not the full 16 MiB. The bank
   selection stops before the cap. Raising the cap to 32 MiB produced a
   byte-identical archive, so more bank space alone does not buy more anchors.
5. **The result is not marginal.** The tuned archive is 3.5 MB smaller than xz
   on this input, with faster compression than xz and faster compression and
   decompression than the 4 MiB default. xz still restores faster.

## Query verification

Slice queries were checked on the 8/8, 16/2, and 16/8 archives. Records 0, 12,
and 23 returned exact sequence bytes in every case, including the last record.
A query for an out-of-range coordinate returns an empty slice; that is a
coordinate error, not an archive failure.

## Limits

- One pangenome from one species, one host, one full grid, and one repeat of the
  two best configurations. The grid was not repeated on another species; the
  independent runs in [`INDEPENDENT.md`](INDEPENDENT.md) test only the tuned
  profile, not the full grid.
- xz was not rerun; its baseline comes from the same input and pinned tool.
- Peak RSS includes compression and decompression together.
- The ablation does not test whether the tuned profile generalizes to other
  pangenomes, read data, or much larger collections.
- The result is a tuning win, not a new compression law. It removes an arbitrary
  cap; it does not establish universal superiority.

## Reproduce

```sh
javac -cp bin -d bin benchmarks/ablation/AnchorBankAblation.java
python -B benchmarks/ablation/run.py --run-id anchor-ablation-1
python -B benchmarks/ablation/run.py --run-id anchor-ablation-2 --banks 8,16 --blocks 8
python -B benchmarks/ablation/run.py --run-id anchor-ablation-32 --banks 32 --blocks 8
```

The harness and protocol are in [`benchmarks/ablation`](../benchmarks/ablation).
Reports are written under `benchmarks/ablation/runs/`. The three published
`report.json` files are force-added as evidence; the work directories and
archives remain gitignored.
