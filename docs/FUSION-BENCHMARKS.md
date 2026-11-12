# HyperGraSS + RC-MDL Unified Engine (`HGS2` / `EIMDL`) — Verified Benchmarks

All measurements in this report are produced by [`benchmarks/fusion/run_fusion_benchmarks.py`](../benchmarks/fusion/run_fusion_benchmarks.py) and committed as machine-readable JSON under [`benchmarks/fusion/results/`](../benchmarks/fusion/results/). Every single archive reported below was decompressed in a fresh OS process and verified against the original input's **byte length and SHA-256** (`verified_sha256: true`).

In strict adherence to [`docs/FOUNDATIONS.md`](FOUNDATIONS.md), the three operational compression regimes are **never conflated**:
- **Mode S (`mode=0`):** Self-contained single genome / sequence ($\mathcal{R} = \emptyset$).
- **Mode R (`mode=1`):** Reference-conditioned delta archive relative to an external reference panel pinned by a 32-byte Merkle SHA-256 manifest digest ($\mathcal{R} \neq \emptyset$).
- **Mode C (`mode=2`):** Self-contained pangenome collection archive over $N$ genomes ($\mathcal{R} = \emptyset$, causal phylogenetic Prim DAG + unified pan-genome novelty superstream).

---

## 1. Mode S — Self-Contained Single Genomes (`mode=0`)

Source JSON: [`benchmarks/fusion/results/mode_s_and_synthetic.json`](../benchmarks/fusion/results/mode_s_and_synthetic.json)

| Dataset | Raw Bytes | DNA Bases | `HyperGraSS v1 --single` | `RC-MDL S` | **`HGS2` / `EIMDL` (Default Parallel)** | **`EIMDL` bpb** | Ratio | Compress / Restore | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| *E. coli* `00.fa` | 4,856,085 | 4,856,080 | 1,136,083 B *(host: 1,136,069 B)* | 1,192,494 B | **1,122,777 B** | **1.8497** | **4.33x** | 6.75 s / 6.57 s | PASS |
| *E. coli* `00.seq` | 4,856,080 | 4,856,080 | — | — | **1,122,749 B** | **1.8496** | **4.33x** | 8.53 s / 7.71 s | PASS |
| *M. tuberculosis* `00.seq` | 4,417,521 | 4,417,521 | — | — | **969,382 B** | **1.7555** | **4.56x** | 7.49 s / 6.25 s | PASS |
| *S. aureus* `00.seq` | 2,845,025 | 2,845,025 | — | — | **632,236 B** | **1.7778** | **4.50x** | 5.51 s / 5.67 s | PASS |

- Powered by the **6-State Causal Open-Reading-Frame HMM (`6F-ORF-HMM`, Law IIb)** and **Hierarchical Pitman-Yor Backoff + Two-Adic Tapered Adaptation (Laws VIa–VIc)** in [`FastNoveltyCoder.java`](../src/hypergrass/FastNoveltyCoder.java), `HGS2` Mode S beats `HyperGraSS v1 --single` by **`-13,306 B` (`-1.17%`)** vs. the Docker control (**`-13,292 B` [`-1.17%`]** vs. the local host control) and beats `RC-MDL S` by **`-69,717 B` (`-5.85%`)** on `00.fa`.

---

## 2. Mode R — Reference-Conditioned Leave-One-Out (LOO) Panels (`mode=1`)

Every `HGS2` Mode R archive embeds the compressed donor manifest plus a mandatory **32-byte cryptographic SHA-256 panel digest** (`L_pin(R)`) verified prior to decoding. All 73/73 genomes across all three species run with **one universal default configuration** (zero per-species or per-target tuning).

### 2.1 Multi-Species Leave-One-Out Summary

Sources: [`loo_ecoli.json`](../benchmarks/fusion/results/loo_ecoli.json), [`loo_mtb.json`](../benchmarks/fusion/results/loo_mtb.json), [`loo_saureus.json`](../benchmarks/fusion/results/loo_saureus.json)

| Species Panel | Genomes | Total Bases | `AGC 3.2.5` (Best Single) | `AGC 3.2.5` (13-Config Oracle) | `RC-MDL` (v2) | **`HGS2` / `EIMDL` (1 Universal Config)** | **`EIMDL` Total bpb** | **Median bpb** | Wins vs `RC-MDL` | Wins vs `AGC` (Single / Oracle) | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|:---:|
| ***E. coli*** | 24 | 120,108,223 | 8,384,117 B | — | 1,780,761 B | **1,174,113 B** (`-34.07%` vs RC-MDL; **`7.14x`** < AGC) | **0.07820** | **0.07162** | **24 / 24** | **24 / 24** *(—)* | 24/24 PASS |
| ***M. tuberculosis*** | 24 | 105,745,015 | 292,170 B | 96,580 B | 122,998 B | **47,339 B** (`-61.51%` vs RC-MDL; **`2.04x` < AGC Oracle**) | **0.00358** | **0.00367** | **24 / 24** | **24 / 24** *(19 / 24)* | 24/24 PASS |
| ***S. aureus*** | 25 | 71,043,114 | 2,979,078 B | 1,617,085 B | 267,412 B | **174,217 B** (`-34.85%` vs RC-MDL; **`9.28x` < AGC Oracle**) | **0.01962** | **0.00193** | **25 / 25** | **25 / 25** *(24 / 25)* | 25/25 PASS |

### 2.2 Key Scientific Milestones Resolved in Mode R
1. ***M. tuberculosis* Oracle Barrier Broken (`47,339 B` vs `96,580 B`):** Previously (`HANDOFF.md`), `AGC 3.2.5`'s 13-configuration per-target oracle (`96,580 B`) beat `RC-MDL` (`122,998 B`) by `21.5%`. With **Horizon-Bounded Two-Gear Coalescent Haplotype Copying (Law IIIb)** and **Laws VIa–VIb**, a single universal `HGS2` configuration compresses all 24 *M. tuberculosis* genomes (`105.7 Mbp`) into **`47,339 B` (`0.00358 bpb`, `2,233.78x` aggregate ratio)** — **`50.98%` (`2.04x`) smaller than AGC's 13-config per-target oracle**, **`6.17x` smaller than AGC's best single config**, and winning `19/24` targets against the oracle (yielding only slightly on 5 ultra-close genomes: `03.seq` `256 vs 243 B`, `08.seq` `287 vs 224 B`, `09.seq` `607 vs 368 B`, `16.seq` `1,150 vs 1,072 B`, and `19.seq` `536 vs 512 B`).
2. ***E. coli* 24/24 Clean Sweep Over AGC (Including Twin Pair `13.seq` & `14.seq`):** Previously, `AGC` beat `RC-MDL` on near-identical twins `13.seq` (`1,680 B` vs `2,032 B`) and `14.seq` (`1,624 B` vs `2,266 B`). `HGS2` compresses `13.seq` (`5.185 Mbp`) to **`633 B` (`0.00098 bpb`, `8,191.2x`)** and `14.seq` (`5.182 Mbp`) to **`632 B` (`0.00098 bpb`, `8,200.0x`)**, achieving a **24/24 sweep**.
3. **Conditional Monotonicity Verified on Clonal Outbreaks:**
   - *M. tuberculosis* `03.seq` (`4,367,416 bp`): **`256 B` (`0.000469 bpb`, `17,060x`, `18` substitutions vs `RC-MDL`'s `261` mismatches and `1,635 B`)**.
   - *S. aureus* `20.seq` (`2,879,034 bp`): **`134 B` (`0.000372 bpb`, `21,485x`, `1` syntenic token, `2` substitutions, including 32-byte SHA-256 manifest digest)**, while winning `24/25` against `AGC`'s oracle (`04.seq` `683` vs `617 B`).

### 2.3 Leave-Cluster-Out Unseen-Lineage Controls

Sources: [`cluster_loo_ecoli.json`](../benchmarks/fusion/results/cluster_loo_ecoli.json), [`cluster_loo_mtb.json`](../benchmarks/fusion/results/cluster_loo_mtb.json), [`cluster_loo_saureus.json`](../benchmarks/fusion/results/cluster_loo_saureus.json)

All members of the target's Mash/Jaccard cluster are excluded from the reference panel to measure generalization to unseen lineages:

| Species | Jaccard Cutoff | Clusters (Max Size) | `AGC 3.2.5 (-a)` | `RC-MDL` Cluster-LOO | **`HGS2` / `EIMDL` Cluster-LOO** | **`EIMDL` bpb** | Delta vs `RC-MDL` | Factor vs `AGC (-a)` | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| ***E. coli*** | 0.98 | 22 (max 3) | 9,203,994 B | 2,331,696 B | **1,619,296 B** | **0.10786** | **-30.55%** | **5.68x** | 24/24 PASS |
| ***M. tuberculosis*** | 0.99 | 13 (max 6) | 1,188,909 B | 153,182 B | **78,603 B** | **0.00595** | **-48.69%** | **15.13x** | 24/24 PASS |
| ***S. aureus*** | 0.98 | 14 (max 8) | 4,132,771 B | 741,650 B | **524,535 B** | **0.05907** | **-29.27%** | **7.88x** | 25/25 PASS |

---

## 3. Mode C — Self-Contained Pangenome Collections (`mode=2`)

Source JSON: [`benchmarks/fusion/results/mode_c_collections.json`](../benchmarks/fusion/results/mode_c_collections.json)

### 3.1 Canonical 24-Genome *E. coli* Flat Sequence (`120,108,223 B`, SHA-256 `51cb3ecb...`)

| Compressor | Archive Bytes | Bits / Base | Ratio | Compress Time | Restore Time | SHA-256 |
|---|---:|---:|---:|---:|---:|:---:|
| `GraSS (-e2)` | 24,777,592 B | 1.6503 | 4.85x | 224.5 s | — | PASS |
| `AGC 3.2.5` *(multi-FASTA)* | 8,540,341 B | 0.5688 | 14.06x | — | — | PASS |
| `xz -9` *(multi-FASTA)* | 8,394,900 B | 0.5592 | 14.51x | 186.3 s | — | PASS |
| `RC-MDL Prim` | 3,972,975 B | 0.2646 | 30.23x | — | — | PASS |
| `JARVIS3 (level 7)` | 3,896,361 B | 0.2595 | 30.83x | 95.4 s | 84.2 s | PASS |
| `HyperGraSS v1` *(tuned + LZ gap)* | 3,810,409 B | 0.2538 | 31.52x | 887.3 s | 82.4 s | PASS |
| **`HGS2` / `EIMDL` (`--threads 1`)** | **3,371,489 B** | **0.2246** | **35.62x** | **40.31 s** | **56.78 s** | **PASS** |
| **`HGS2` / `EIMDL` (Default Parallel `8T` / `12T`)** | **3,388,639 B** | **0.2257** | **35.44x** | **19.65 s** *(12T)* / **33.78 s** *(default)* | **10.40 s** *(12T, 8.66 s @ 8T)* | **PASS** |
| **`HGS2` / `EIMDL` (`--bounded -Xmx256m`)** | **3,422,314 B** | **0.2279** | **35.10x** | **22.38 s** *(24.03 s @ 224m)* | **10.61 s** *(13.89 s @ 224m)* | **PASS** |
| **`HGS2` / `EIMDL` (`--bounded -Xmx192m`)** | **3,452,703 B** | **0.2300** | **34.79x** | **62.68 s** | **27.46 s** | **PASS** |
| **`HGS2` / `EIMDL` (`--bounded -Xmx128m`)** | **3,452,703 B** | **0.2300** | **34.79x** | **25.59 s** | **11.67 s** | **PASS** |

- `HGS2` Mode C default parallel (`3,388,639 B`) beats `HyperGraSS v1` (`3,810,409 B`) by **`-421,770 B` (`-11.07%`)** at **`26.3x–45.2x` faster compression speed** (`33.78 s` default / `19.65 s` at `--threads 12` vs `887.3 s`), beats `JARVIS3 l7` (`3,896,361 B`) by **`-507,722 B` (`-13.03%`)**, and beats `RC-MDL Prim` (`3,972,975 B`) by **`-584,336 B` (`-14.71%`)**.

### 3.2 Multi-FASTA Full Pangenome Collections (6-Channel Exact File Reconstruction)

| Collection | Raw File Bytes | DNA Bases | Prior Baseline | **`HGS2` / `EIMDL` Archive** | **DNA bpb** | **File Ratio** | Compress / Restore | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| *E. coli* (24 genomes, multi-FASTA) | 121,825,884 | 120,108,223 | `HyperGraSS --pangenome`: 4,885,851 B | **3,387,485 B** (`-30.67%`) | **0.2256** | **35.96x** | 32.5 s / 16.5 s | PASS |
| *M. tuberculosis* (24 genomes, multi-FASTA) | 107,257,786 | 105,745,015 | `xz -9`: ~2.3 MB | **1,049,266 B** | **0.0794** | **102.22x** | 41.4 s / 24.8 s | PASS |
| *S. aureus* (25 genomes, multi-FASTA) | 72,060,360 | 71,043,114 | `RC-MDL Prim`: 1,233,433 B / `AGC`: 2,145,694 B | **1,086,129 B** (`-11.94%` vs RC-MDL; **`1.98x`** < AGC) | **0.1223** | **66.35x** | 31.8 s / 16.1 s | PASS |

### 3.3 Independent Locked Pangenomes (*Salmonella* + *Klebsiella*)

Source JSON: [`salmonella`](../benchmarks/independent/runs/salmonella/trio-1/report.json), [`klebsiella`](../benchmarks/independent/runs/klebsiella/trio-1/report.json). Harness: [`benchmarks/independent/run_trio.py`](../benchmarks/independent/run_trio.py).

Both panels were re-fetched by exact accession from their locked manifests and accepted only after matching the recorded byte counts and SHA-256. JARVIS3 v3.7 accepts raw DNA only (`currently JARVIS3 only supports DNA format`), and original GraSS's raw path is headerless, so all three codecs received the **identical canonical flat sequence** (headers removed, uppercased, non-ACGT masked to `A`; both panels contained zero non-ACGT bases). Archives therefore compare like with like, with no metadata asymmetry.

| Dataset (genomes) | DNA Bases | `JARVIS3 -l 7` | original `GraSS -e2` | **`HGS2` / `EIMDL`** (`--collection --threads 1`) | **`EIMDL` bpb** | **Ratio** | vs `JARVIS3` | vs `GraSS` | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| *Salmonella enterica* (25) | 120,527,038 | 2,709,545 B | 19,991,646 B | **2,390,106 B** | **0.1586** | **50.43x** | **`-11.79%`** (`-319,439 B`) | **`-88.04%`** (`-17,601,540 B`) | PASS |
| *Klebsiella pneumoniae* (23) | 122,700,601 | 2,948,926 B | 21,191,566 B | **2,523,991 B** | **0.1646** | **48.61x** | **`-14.41%`** (`-424,935 B`) | **`-88.09%`** (`-18,667,575 B`) | PASS |

- **Verified round trips:** each archive was decompressed in a fresh container and checked by [`benchmarks/verify/verify.py`](../benchmarks/verify/verify.py) (streaming byte comparison plus SHA-256). Both panels: `PASS`.
- **Environment:** `hypergrass-corpus:jdk27`, pinned by digest `sha256:580e58b6ee16d43f4dd99a2871b1c920548b7c6e9b0c695a01f521c17701a4ae`, 2 CPUs, strictly sequential runs; EIMDL `-Xmx2200m` under a 2800 MiB cap, JARVIS3 and GraSS under 2 GiB.
- **Baselines:** JARVIS3 v3.7 at level 7; original GraSS commit `6d93700f5a2de8841f7cf2c20064413c80df2c4a` with libbsc `1dc676540cb40a41562af6f2db04bd1adce97f90` (`bsc` 3.4.0, SHA-256 `84f7960dc211b5b146f45e0b6186569dd73db96ee97e36ecb534728c63e64640`), default `-e2` path.
- **Times (compress / restore):** EIMDL 100.8 / 53.4 s and 32.7 / 16.2 s; JARVIS3 90.0 / 85.9 s and 55.3 / 52.0 s; GraSS 212.3 / 94.5 s and 87.2 / 75.0 s. Host load varied between runs; archive sizes and round trips are the stable measurements.
- **Versus prior generation:** the earlier tuned `HyperGraSS v1` results on these panels were 3,443,600 B (0.2286 bpb) and 3,464,423 B (0.2259 bpb); `EIMDL` is **30.59%** and **27.15%** smaller on the same DNA-base denominator.

### 3.4 Mode S / Mode R 25% Subset with JARVIS3, GraSS, and RSS

Source JSON: [`salmonella`](../benchmarks/independent/runs/salmonella/modes-subset-1/report.json), [`klebsiella`](../benchmarks/independent/runs/klebsiella/modes-subset-1/report.json). Harness: [`benchmarks/independent/run_modes.py`](../benchmarks/independent/run_modes.py).

The full per-genome mode campaign was cut to **25%**: six evenly spaced genomes per panel, always including the first and last record (6 of 25 *Salmonella*, 6 of 23 *Klebsiella*). Every record was written as its canonical flat sequence, so Mode S compares like-for-like with the two reference-free baselines on identical bytes. Mode R is `EIMDL`-only leave-one-out against the full 24/22-genome panel; JARVIS3 v3.7 and original GraSS expose no reference mode, so their columns are context for Mode R, not a like-for-like baseline.

**Archive size**

| Panel (targets) | `EIMDL` S | `JARVIS3 -l 7` | original `GraSS -e2` | `EIMDL` R | S bpb | R bpb | S vs `JARVIS3` | S vs `GraSS` | R vs S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *Salmonella* (6) | **6,629,508 B** | 7,162,390 B | 7,029,710 B | **179,303 B** | **1.8540** | **0.0501** | **`-7.44%`** | **`-5.69%`** | **`-97.30%`** |
| *Klebsiella* (6) | **7,238,897 B** | 8,041,937 B | 7,786,274 B | **254,195 B** | **1.7983** | **0.0632** | **`-9.99%`** | **`-7.03%`** | **`-96.49%`** |
| **Both (12)** | **13,868,405 B** | 15,204,327 B | 14,815,984 B | **433,498 B** | **1.8245** | **0.0570** | **`-8.79%`** | **`-6.40%`** | **`-96.87%`** |

**Speed and memory (12 targets, 60,810,308 bp)**

| Regime | Compress wall | Restore wall | Compress CPU | Restore CPU | Compress MiB/s | Restore MiB/s | cgroup C (med–max) | cgroup R (med–max) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **`EIMDL` S** | 332.5 s | 277.9 s | 368.8 s | 299.3 s | 0.17 | 0.21 | **108–113 MiB** | **88–89 MiB** |
| `JARVIS3` S | **91.0 s** | 87.9 s | **52.3 s** | **49.3 s** | **0.64** | 0.66 | 266–267 MiB | 266–267 MiB |
| `GraSS` S | 102.6 s | **81.2 s** | 83.6 s | 63.7 s | 0.57 | 0.71 | 265–310 MiB | 132–135 MiB |
| **`EIMDL` R** | 398.2 s | **48.4 s** | 451.9 s | 62.2 s | 0.15 | **1.20** | 343–379 MiB | 222–234 MiB |

- **Size:** `EIMDL` Mode S is **8.79% smaller than JARVIS3** (−1,335,922 B) and **6.40% smaller than original GraSS** (−947,579 B) across the 12 genomes. The margin is larger on *Klebsiella* (−9.99% / −7.03%) than on *Salmonella* (−7.44% / −5.69%).
- **Mode R:** 433,498 B total at **0.0570 bpb**, **31.99x smaller** than `EIMDL` Mode S. Per-target wins range from **−90.78%** to **−99.98%**. Near-identical genomes collapse to **271–1,311 B**; the most divergent targets still need 41.8–112.5 KB.
- **Speed:** Mode S is **3.65x slower to compress** and **3.16x slower to restore** than JARVIS3 (3.24x / 3.42x versus GraSS). Mode R is **1.20x slower to compress** than Mode S (it loads the full panel) but **5.74x faster to restore**; against JARVIS3, Mode R compresses 4.38x slower and restores 1.82x faster. By CPU time the Mode S gap is wider (**7.05x** compress / **6.07x** restore versus JARVIS3); Mode R restore uses **4.81x less CPU** than Mode S, though its 62.2 CPU-s is 1.26x JARVIS3's 49.3 CPU-s because JARVIS3 is wall-clock I/O-bound.
- **Memory:** Mode S uses **2.46x less compress RSS** and **3.02x less restore RSS** than JARVIS3 (2.45x / 1.50x versus GraSS). Mode R uses **3.18x more compress RSS** than Mode S because of the panel, but **1.20x less restore RSS** than JARVIS3. Peak cgroup was **379 MiB** (Mode R compress) against a 2,200 MiB cap; no OOM kills across all 48 runs (96 compress/restore phases).
- **Protocol:** `hypergrass-corpus:jdk27` at digest `sha256:580e58b6ee16d43f4dd99a2871b1c920548b7c6e9b0c695a01f521c17701a4ae`, 2 CPUs, `--threads 1` for `EIMDL`. Heaps/caps: Mode S `-Xmx1200m` / 1800 MiB, Mode R `-Xmx1600m` / 2200 MiB, JARVIS3 and GraSS 2 GiB. Every archive was decompressed in a fresh container and byte/SHA-256 verified against its target.
- **Measurement:** wall is the in-container wrapper duration (Docker creation/startup excluded); CPU and RSS are container-side `os.wait4` rusage of the tool process (including descendants it reaped); cgroup is v2 `memory.peak`. Speed is host-dependent, so treat it as a tool-process ranking, not a storage-system benchmark.
- **Baseline caveat:** the earlier *E. coli* Mode S numbers in §1 were measured with default parallelism; these runs are single-threaded, so absolute times are not comparable. A smoke rerun reproduced `EIMDL` and JARVIS3 archives byte-identically but GraSS moved by 6 B on a 1.14 MB archive (5 ppm), so small GraSS deltas are noise.

**Tuned baseline control: JARVIS3 `-l 31`**

JARVIS3's default is level 7. To test whether the Mode S margin above is a tuning artifact, the same 12 canonical flat targets were recompressed at level 31 (4 context models, 2 repeat-class models, NN mixer, `-hs 48`) under the identical 2 GiB cap, with a fresh-container restore and byte/SHA-256 verification. Harness: [`run_jarvis_level.py`](../benchmarks/independent/run_jarvis_level.py).

| Panel (targets) | `EIMDL` S | `JARVIS3 -l 31` | `JARVIS3 -l 7` | EIMDL vs L31 | L31 vs L7 | L31 compress/restore wall | L31 cgroup C/R |
|---|---:|---:|---:|---:|---:|---:|---:|
| *Salmonella* (6) | **6,629,508 B** | 6,748,736 B | 7,162,390 B | **`-1.77%`** | `-5.78%` | 100.0 / 106.4 s | 1,158 / 1,157 MiB |
| *Klebsiella* (6) | **7,238,897 B** | 7,448,657 B | 8,041,937 B | **`-2.82%`** | `-7.38%` | 317.5 / 345.7 s | 1,159 / 1,159 MiB |
| **Both (12)** | **13,868,405 B** | 14,197,393 B | 15,204,327 B | **`-2.32%`** | `-6.62%` | 417.5 / 452.1 s | 1,158–1,159 MiB |

- **Verified:** 12/12 targets round-tripped byte-exactly, 0 OOM kills, 1.16 GiB peak against the 2 GiB cap.
- **Memory:** level 31 uses **10.7x / 13.2x** more cgroup RSS than `EIMDL` Mode S (1,158–1,159 MiB vs 108 / 88 MiB).
- **Speed:** combined, `EIMDL` is **1.26x faster to compress** and **1.63x faster to restore** (417.5 / 452.1 s versus 332.5 / 277.9 s); compress CPU is essentially tied (363.0 vs 368.8 s) and restore CPU favors `EIMDL` by 1.30x. Per panel the ranking flips: JARVIS3 L31 is faster on *Salmonella* (100.0 / 106.4 s vs 153.1 / 130.9 s) but 1.77x / 2.35x slower on *Klebsiella*.
- **Effect on the claim:** tuning JARVIS3 narrows the Mode S size margin from **8.79% to 2.32%**, while costing ~11x the memory. Level 7 numbers earlier in this section are still the published default-setting results; the L31 row is the tuned control.
- **Mode C caveat:** a full-panel level-31 run was started but **not completed**. JARVIS3 consumed ~40 min to compress the 120.5 Mb *Salmonella* stream (~3x the per-base cost measured on the subset, i.e. superlinear in stream length) and was stopped deliberately; no full-panel L31 archive exists. Mode C comparisons in §3.3 therefore remain versus `-l 7` and `GraSS` only.

---

## 4. Amdahl Multi-Core Scaling & Strict Bounded-Heap Execution

Source JSON: [`benchmarks/fusion/results/amdahl_and_memory.json`](../benchmarks/fusion/results/amdahl_and_memory.json)

### 4.1 Multi-Core Scaling on `120,108,223 B` *E. coli* Collection

| Threads (`--threads`) | Archive Bytes | Size Delta vs `1T` | Encode Time | Encode Speedup vs `v1` (`887.3 s`) | Decode Time | Decode Speedup vs `1T` | SHA-256 |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| `1` | **3,371,489 B** | `0.00%` | 40.31 s | **22.01x** | 56.78 s | 1.00x | PASS |
| `2` | **3,373,077 B** | `+0.05%` | 37.44 s | **23.70x** | 12.58 s | **4.51x** | PASS |
| `4` | **3,383,225 B** | `+0.35%` | 25.47 s | **34.83x** | 12.58 s | **4.51x** | PASS |
| `8` | **3,388,639 B** | `+0.51%` | 21.32 s | **41.61x** | **8.66 s** | **6.56x** | PASS |
| `12` | **3,388,639 B** | `+0.51%` | **19.65 s** | **45.16x** | 10.40 s | **5.46x** | PASS |

### 4.2 Strict Bounded-Heap Execution on Pure JVM (`120,108,223 B` *E. coli* Collection)

Previously (`docs/FOUNDATIONS.md` §6), every pure-JVM bounded profile of `HyperGraSS v1` OOM-killed under a strict 256 MiB container cap (`0/7` natural collections passed), requiring a GraalVM native image (`3,864,464 B` at `-Xmx192m`) to beat `JARVIS3` (`3,896,361 B`). With **Laws V and VII** (streaming 64 KiB flat-file parser, `0.25 B/base` 2-bit bank, 2-pass zero-materialization CSR index, 32-bit `int[]` per-worker syncmer hashes, strict Stage-1 vs. Stage-2 lifecycle release, and adaptive `ultraTight`/`extremeTight` thread throttling down to $P \le 2$), pure JVM scales all the way to **`-Xmx128m`**:

| Pure JVM Heap Cap (Encode & Decode) | Subsampling (`1 / 2^s`) | Archive Bytes | Bits / Base | Delta vs Unbounded (`3,388,639 B`) | Margin vs `JARVIS3` (`3,896,361 B`) | Margin vs `HyperGraSS v1` (`3,810,409 B`) | Encode / Decode | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| **`-Xmx256m`** | $s = 1$ (`1/2`) | **3,422,314 B** | **0.2279** | **`+0.99%`** *(passes `<2%` gate)* | **`-12.17%`** (`-474,047 B`) | **`-10.19%`** (`-388,095 B`) | 22.38 s / 10.61 s | **PASS** |
| **`-Xmx224m`** | $s = 1$ (`1/2`) | **3,422,314 B** | **0.2279** | **`+0.99%`** *(passes `<2%` gate)* | **`-12.17%`** (`-474,047 B`) | **`-10.19%`** (`-388,095 B`) | 24.03 s / 13.89 s | **PASS** |
| **`-Xmx192m`** | $s = 2$ (`1/4`) | **3,452,703 B** | **0.2300** | **`+1.89%`** *(passes `<2%` gate)* | **`-11.39%`** (`-443,658 B`) | **`-9.39%`** (`-357,706 B`) | 62.68 s / 27.46 s | **PASS** |
| **`-Xmx160m`** | $s = 2$ (`1/4`) | **3,452,703 B** | **0.2300** | **`+1.89%`** *(passes `<2%` gate)* | **`-11.39%`** (`-443,658 B`) | **`-9.39%`** (`-357,706 B`) | 55.90 s / 14.34 s | **PASS** |
| **`-Xmx128m`** | $s = 2$ (`1/4`, `extremeTight`) | **3,452,703 B** | **0.2300** | **`+1.89%`** *(passes `<2%` gate)* | **`-11.39%`** (`-443,658 B`) | **`-9.39%`** (`-357,706 B`) | 25.59 s / 11.67 s | **PASS** |

---

## 5. Non-Existing, Xenobiological, Synthetic & Adversarial Genomes

Source JSON: [`benchmarks/fusion/results/mode_s_and_synthetic.json`](../benchmarks/fusion/results/mode_s_and_synthetic.json)

| Corpus ID | Regime | Biological / Adversarial Stress Design | Raw Bytes | `gzip -9` | **`HGS2` / `EIMDL`** | **DNA bpb** | SHA-256 |
|---|---|---|---:|---:|---:|---:|:---:|
| `01_euk_centromere` | Mode S | Higher-order alpha-satellite repeats + HOR inversions (`485 kbp`) | 491,604 B | 83,797 B | **42,075 B** | **0.6934** | PASS |
| `02_viral_rna_acgu` | Mode C | Pure RNA (`A/C/G/U`) segmented quasispecies swarm (`299 kbp`) | 304,013 B | 42,923 B | **11,995 B** | **0.3207** | PASS |
| `03_chimeric_iupac_crlf` | Mode S | Mixed DNA+RNA (`T`+`U`), all 11 IUPAC codes, soft-masking, `N`-scaffolds, ragged `CRLF` | 263,817 B | 86,566 B | **51,536 B** | **1.6893** *(1.563 raw)* | PASS |
| `04_xeno_at88` | Mode S | Extreme `88% AT` xenobiological / *Carsonella*-like genome (`500 kbp`) | 506,294 B | 123,139 B | **94,248 B** | **1.5080** | PASS |
| `05_xeno_gc86` | Mode S | Extreme `86% GC` actinomycete / synthetic genome (`500 kbp`) | 506,291 B | 126,076 B | **97,525 B** | **1.5604** | PASS |
| `06_denovo_grammar` | Mode S | Non-existing synthetic regulatory grammar + palindromic hairpins (`500 kbp`) | 506,287 B | 141,540 B | **41,159 B** | **0.6585** | PASS |
| `07_hypermutator_pangenome` | Mode C | 6-strain `6%` hypermutator + heavy indel pangenome (`1.495 Mbp`) | 1,513,506 B | 454,203 B | **102,503 B** | **0.5486** | PASS |
| `07_hypermutator_strain05` | Mode R | Single `6%` hypermutator strain vs 5-donor panel (`249 kbp`) | 249,436 B | — | **7,809 B** | **0.2505** | PASS |
| `08_crypto_white_noise` | Mode S | Cryptographic `SHA-256` counter-mode white noise (`500,000 bp`, $H_0 = 2.0000$) | 500,000 B | 144,017 B (`2.304 bpb`) | **125,075 B** | **2.0012** *($\le 2.005$ gate)* | PASS |

---

## 6. Real-Genome Law-by-Law Ablations

Source JSON: [`benchmarks/fusion/results/ablations_real_genomes.json`](../benchmarks/fusion/results/ablations_real_genomes.json)

| Ablation Variant | *E. coli* `00.seq` vs Distant Single Donor `20.seq` (`54,124` SNPs, `619` splices) | *E. coli* `00.seq` LOO (23 Donors) | *M. tuberculosis* `00.seq` LOO (23 Donors) |
|---|---:|---:|---:|
| **`full_eimdl` (All Laws I–VII Active)** | **197,163 B** (`0.3248 bpb`) | **40,588 B** (`0.0669 bpb`) | **1,240 B** (`0.00225 bpb`) |
| `no_kimura_v4_law_I` (4-symbol allele fallback) | 203,340 B (`+6,177 B` / `+3.13%`) | 41,178 B (`+590 B` / `+1.45%`) | 1,261 B (`+21 B` / `+1.69%`) |
| `no_codon_wobble_law_II` (flat gap coding) | 199,024 B (`+1,861 B` / `+0.94%`) | 40,697 B (`+109 B` / `+0.27%`) | 1,237 B (`-3 B`, only 251 SNPs) |
| `no_indel_splice_law_III` (no micro-indel splicing) | 196,936 B (`-227 B`, `+428` broken tokens `1,945` vs `1,517`) | 40,469 B (`-119 B` on dense 23-donor panel) | 1,240 B (`0 B`, `0` splices) |
| `no_v4_no_wobble_no_splice` (all 3 disabled) | 205,071 B (`+7,908 B` / `+4.01%`) | 41,175 B (`+587 B` / `+1.45%`) | 1,258 B (`+18 B` / `+1.45%`) |
| `exact_matches_only_no_approx` (exact copy only) | 238,550 B (`+41,387 B` / **`+20.99%`**) | 43,357 B (`+2,769 B` / **`+6.82%`**) | 1,267 B (`+27 B` / **`+2.18%`**) |
| Prior `RC-MDL` v2 Baseline | 249,491 B (`+52,328 B` / **`+26.54%`**) | 73,194 B (`+32,606 B` / **`+80.33%`**) | 4,012 B (`+2,772 B` / **`+223.5%`**) |
