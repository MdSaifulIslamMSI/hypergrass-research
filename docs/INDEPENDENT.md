# Independent pangenome replication

## Question

The 41.8% win over xz on the *E. coli* pangenome was measured after the
16 MiB / 8 MiB profile was selected on that same pangenome. That is a tuning
result on one input. Does it survive on pangenomes that played no part in
choosing the profile?

## Protocol

- Two additional NCBI RefSeq complete-genome pangenomes:
  *Salmonella enterica* (25 records, 120,527,038 bp, 122,251,418-byte FASTA)
  and *Klebsiella pneumoniae* (23 records, 122,700,601 bp, 124,455,339-byte
  FASTA).
- Both use the same selection rule as the *E. coli* set: NCBI relevance order,
  distinct accessions, records between 3 and 7.5 Mbp, stop once the target is
  reached. Accessions, lengths, and SHA-256 hashes are in
  `benchmarks/independent/data/<name>/manifest.json`.
- The locked Linux image, 2 CPUs, 2 GiB container limit, 1 GiB Java heap,
  fresh process per operation, and SHA-256 round-trip check are the same as the
  published scale test.
- Profiles: the published default (`hypergrass`), the scale-tuned
  `--pangenome` profile (16 MiB bank / 8 MiB blocks), `gzip -9`, and `xz -9`.
  No parameter was changed after seeing these results.
- Each HyperGraSS archive is additionally queried for records 0, the middle
  record, and the last record; every slice must match the source bytes.
- One full run per dataset. Times are indicative, as in the original scale
  test; archive sizes and round trips are the stable measurements.

<!-- INDEPENDENT-TABLES:START -->
| Dataset | Profile | Archive | Ratio | Bits/base | Wall | Peak RSS | Exact | Δ vs xz |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: | ---: |
| klebsiella | hypergrass | 8,571,498 B | 14.52x | 0.559 | 65.0 s | 1147 MiB | yes | +33.12% |
| klebsiella | hypergrass-pangenome | 3,464,423 B | 35.92x | 0.226 | 10.6 s | 912 MiB | yes | -46.19% |
| klebsiella | gzip-9 | 36,388,919 B | 3.42x | 2.373 | 201.1 s | 10 MiB | yes | +465.16% |
| klebsiella | xz-9 | 6,438,744 B | 19.33x | 0.420 | 521.3 s | 736 MiB | yes | — |
| salmonella | hypergrass | 7,317,099 B | 16.71x | 0.486 | 60.9 s | 922 MiB | yes | +19.56% |
| salmonella | hypergrass-pangenome | 3,443,600 B | 35.50x | 0.229 | 28.6 s | 732 MiB | yes | -43.73% |
| salmonella | gzip-9 | 35,969,250 B | 3.40x | 2.387 | 97.6 s | 10 MiB | yes | +487.71% |
| salmonella | xz-9 | 6,120,260 B | 19.97x | 0.406 | 184.6 s | 734 MiB | yes | — |

Tuned profile vs xz: median -44.96% across 2 independent pangenomes.
<!-- INDEPENDENT-TABLES:END -->

## Findings

The tuned profile is smaller than xz on both independent pangenomes:

- *Salmonella enterica*: 3,443,600 B vs 6,120,260 B, **43.7% smaller**.
- *Klebsiella pneumoniae*: 3,464,423 B vs 6,438,744 B, **46.2% smaller**.

The original *E. coli* pangenome was 41.8% smaller (4,885,851 B vs
8,394,900 B). That measurement is not part of this independent count because
the profile was selected on it. The two independent margins have a median of
44.96% and a range of 43.73–46.19%; the direction and rough magnitude of the
original result therefore replicate.

The default profile loses to xz on both new pangenomes, as it did on *E. coli*:
+19.56% on Salmonella and +33.12% on Klebsiella, versus +9.78% on *E. coli*.
The failure of the 4 MiB bank is also not an *E. coli* accident.

Compression time moved in the same direction as size on these runs: the tuned
profile compressed Salmonella in 28.6 s versus 184.6 s for xz, and Klebsiella
in 10.6 s versus 521.3 s. Restore time was mixed: 9.5 s versus 6.9 s on
Salmonella, and 9.9 s versus 11.4 s on Klebsiella. Do not read these times too
precisely. The gzip baseline took 97.6 s on Salmonella and 201.1 s on
Klebsiella for almost the same input size, which shows how much host variance
is present. Archive sizes and round trips are the stable measurements.

Two details matter for interpretation:

1. **This is still a small, convenience sample.** Three pangenomes from one
   family of Gram-negative bacteria, each built from an NCBI relevance-ordered
   selection, do not represent all pangenomes. The result supports
   "replicates on two additional *Enterobacteriaceae* pangenomes", not
   "wins on pangenomes in general".
2. **The margin depends on the input, not just the profile.** The tuned
   profile's archive size, compression time, and restore time differ between
   the three inputs. The comparison that is stable is the direction: tuned
   < xz < default on every pangenome measured here, with the default losing to
   xz.

## Three-codec comparison (EIMDL / JARVIS3 / original GraSS)

The same two locked panels were compressed with the unified `HGS2` / `EIMDL`
engine and the two external baselines. JARVIS3 v3.7 rejects FASTA input, and
original GraSS's raw path is headerless, so all three codecs received the
identical canonical flat sequence (headers removed, uppercased, non-ACGT
masked to `A`; both panels contained zero non-ACGT bases). Every archive was
restored in a fresh container and checked byte-for-byte.

| Dataset | EIMDL Mode C | JARVIS3 `-l 7` | original GraSS `-e2` | EIMDL bpb | vs JARVIS3 | vs GraSS |
|---|---:|---:|---:|---:|---:|---:|
| salmonella | **2,390,106 B** | 2,709,545 B | 19,991,646 B | **0.1586** | **-11.79%** | **-88.04%** |
| klebsiella | **2,523,991 B** | 2,948,926 B | 21,191,566 B | **0.1646** | **-14.41%** | **-88.09%** |

EIMDL is also 30.59% (Salmonella) and 27.15% (Klebsiella) smaller than the
earlier tuned HyperGraSS v1 archives on the same DNA-base denominator. The
full protocol, timings, pinned tool commits, and SHA-256 verdicts are in
[`FUSION-BENCHMARKS.md`](FUSION-BENCHMARKS.md) §3.3; the machine-readable
reports are under `benchmarks/independent/runs/<dataset>/trio-1/report.json`.

## Mode S / Mode R 25% subset with JARVIS3, GraSS, and RSS

The per-genome mode campaign was cut to 25%: six evenly spaced genomes per
panel. Mode S is like-for-like against JARVIS3 `-l 7` and original GraSS `-e2`
on identical flat bytes; Mode R is EIMDL-only leave-one-out because neither
baseline has a reference mode.

| Panel | EIMDL S | JARVIS3 | GraSS | EIMDL R | S vs JARVIS3 | S vs GraSS | R vs S |
|---|---:|---:|---:|---:|---:|---:|---:|
| salmonella (6) | **6,629,508 B** | 7,162,390 B | 7,029,710 B | **179,303 B** | **-7.44%** | **-5.69%** | **-97.30%** |
| klebsiella (6) | **7,238,897 B** | 8,041,937 B | 7,786,274 B | **254,195 B** | **-9.99%** | **-7.03%** | **-96.49%** |
| both (12) | **13,868,405 B** | 15,204,327 B | 14,815,984 B | **433,498 B** | **-8.79%** | **-6.40%** | **-96.87%** |

- Mode S uses **108 / 88 MiB** median cgroup RSS (compress / restore) versus
  **266 / 266 MiB** for JARVIS3 and **265 / 132 MiB** for GraSS.
- Mode S is **3.65x slower** to compress and **3.16x slower** to restore than
  JARVIS3; Mode R restores **5.74x faster** than Mode S but uses **3.18x more**
  compress RSS.
- All 48 runs round-tripped byte-exactly; 0 OOM kills; peak cgroup 379 MiB
  against a 2,200 MiB cap.

Full speed/memory table, per-target Mode R range, and protocol:
[`FUSION-BENCHMARKS.md`](FUSION-BENCHMARKS.md) §3.4. Reports:
`benchmarks/independent/runs/<dataset>/modes-subset-1/report.json`.

### Tuned JARVIS3 control (`-l 31`)

The published baseline above uses JARVIS3's default level 7. To check whether the
margin is a tuning artifact, the same 12 targets were recompressed at level 31
(4 context models, 2 repeat-class models, NN mixer) under the same 2 GiB cap and
the same fresh-container verification.

| Panel (targets) | EIMDL S | JARVIS3 `-l 31` | JARVIS3 `-l 7` | EIMDL vs L31 | L31 vs L7 |
|---|---:|---:|---:|---:|---:|
| salmonella (6) | **6,629,508 B** | 6,748,736 B | 7,162,390 B | **-1.77%** | -5.78% |
| klebsiella (6) | **7,238,897 B** | 7,448,657 B | 8,041,937 B | **-2.82%** | -7.38% |
| both (12) | **13,868,405 B** | 14,197,393 B | 15,204,327 B | **-2.32%** | -6.62% |

- All 12 targets verified byte-exactly; 0 OOM; L31 peak cgroup 1,158-1,159 MiB.
- L31 uses **10.7x / 13.2x** more RSS than EIMDL Mode S. Combined, EIMDL is
  **1.26x faster to compress** and **1.63x faster to restore**, but JARVIS3 is
  faster on Salmonella (100/106 s vs 153/131 s) and much slower on Klebsiella
  (318/346 s vs 180/147 s).
- Tuning narrows the Mode S size margin from **8.79% to 2.32%**.
- A full-panel level-31 run was started and deliberately stopped: it needed
  ~40 min for the 120.5 Mb Salmonella stream (~3x the subset per-base cost).
  No full-panel L31 result exists; Mode C claims remain versus level 7 and GraSS.

Reports: `benchmarks/independent/runs/<dataset>/jarvis-l31-subset-1/report.json`.

## Limits

- Two independent datasets, one run each, one host.
- The Mode S / Mode R comparison is a 25% subset: 6 of 25 Salmonella and 6 of 23 Klebsiella genomes, evenly spaced. It is not the full panels.
- The tuned JARVIS3 control is also subset-only (`-l 31`, 12 targets). The full-panel level-31 run was stopped after ~40 min of compression, so Mode C is not compared against a tuned JARVIS3.
- The datasets were fetched on 2026-11-01. NCBI relevance order can change;
  the manifests lock the exact accessions and hashes that were measured, not a
  permanent snapshot of the query.
- The FASTA inputs include headers and wrapping; the corpus comparison does
  not. These numbers are comparable to the published 121.8 MB scale test, not
  to the sequence-only corpus.
- Query probes are correctness checks, not a query-latency workload. The
  end-to-end workload is measured separately on the *E. coli* pangenome.

## Reproduce

```sh
# Re-fetch the locked accessions and verify every size and SHA-256 against its manifest.
python -B benchmarks/independent/fetch_locked.py --dataset salmonella
python -B benchmarks/independent/fetch_locked.py --dataset klebsiella

# EIMDL Mode C vs JARVIS3 -l 7 vs original GraSS -e2 on the canonical flat sequence.
python -B benchmarks/independent/run_trio.py --datasets salmonella klebsiella --run-id trio-1

# 25% subset: EIMDL Mode S vs JARVIS3/GraSS (like-for-like) plus Mode R LOO.
python -B benchmarks/independent/run_modes.py --datasets salmonella klebsiella --genomes 6 --run-id modes-subset-1

# Tuned JARVIS3 control on the same subset (level 31).
python -B benchmarks/independent/run_jarvis_level.py --scope subset --level 31 --genomes 6 --run-id jarvis-l31-subset-1

# Legacy HyperGraSS v1 baselines quoted in the tables above.
python -B benchmarks/independent/run.py --dataset salmonella --run-id salmonella-scale-1
python -B benchmarks/independent/run.py --dataset klebsiella --run-id klebsiella-scale-1
python -B benchmarks/independent/summarize.py
```

All runners use the pinned `hypergrass-corpus:jdk27` image and run strictly
sequentially; each archive is decompressed in a fresh container and checked by
`benchmarks/verify/verify.py`.

FASTAs are gitignored; the manifests and `report.json` files are the published
evidence.

## Artifact storage

Run artifacts are local evidence, not repository content: `.gitignore` excludes
`benchmarks/independent/runs/` along with the locked FASTA inputs. After the
`trio-1`, `modes-subset-1`, and `jarvis-l31-subset-1` campaigns the run tree
held **1.85 GiB in 1,151 files** (Salmonella 952 MiB, Klebsiella 943 MiB) on
2026-11-07. Each run directory contains `report.json`, per-operation
`spec.json` and logs, the archive, and the restored sequence.

Delete or archive the tree when disk is tight; it is reproducible from the
locked manifests with the commands above. Copy out any `report.json` you want to
keep, because the tree itself is not committed.
