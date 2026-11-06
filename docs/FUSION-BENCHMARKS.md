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
| *E. coli* `00.fa` | 4,856,085 | 4,856,080 | 1,136,083 B *(host: 1,136,069 B)* | 1,192,494 B | **1,123,893 B** *(`1T`: **1,122,240 B**)* | **1.8515** *(1.8488)* | **4.32x** | 6.18 s / 6.35 s | PASS |
| *E. coli* `00.seq` | 4,856,080 | 4,856,080 | — | — | **1,123,865 B** | **1.8515** | **4.32x** | 6.72 s / 5.18 s | PASS |
| *M. tuberculosis* `00.seq` | 4,417,521 | 4,417,521 | — | — | **978,281 B** | **1.7716** | **4.52x** | 7.54 s / 5.90 s | PASS |
| *S. aureus* `00.seq` | 2,845,025 | 2,845,025 | — | — | **629,323 B** | **1.7696** | **4.52x** | 4.91 s / 3.61 s | PASS |

- Powered by the **6-State Causal Open-Reading-Frame HMM (`6F-ORF-HMM`, Law IIb)** in [`FastNoveltyCoder.java`](../src/hypergrass/FastNoveltyCoder.java), `HGS2` Mode S beats `HyperGraSS v1 --single` by **`-12,190 B` (`-1.07%`)** in parallel default mode (**`-13,843 B` [`-1.22%`]** single-threaded) and beats `RC-MDL S` by **`-68,601 B` (`-5.75%`)** on `00.fa`.

---

## 2. Mode R — Reference-Conditioned Leave-One-Out (LOO) Panels (`mode=1`)

Every `HGS2` Mode R archive embeds the compressed donor manifest plus a mandatory **32-byte cryptographic SHA-256 panel digest** (`L_pin(R)`) verified prior to decoding. All 73/73 genomes across all three species run with **one universal default configuration** (zero per-species or per-target tuning).

### 2.1 Multi-Species Leave-One-Out Summary

Sources: [`loo_ecoli.json`](../benchmarks/fusion/results/loo_ecoli.json), [`loo_mtb.json`](../benchmarks/fusion/results/loo_mtb.json), [`loo_saureus.json`](../benchmarks/fusion/results/loo_saureus.json)

| Species Panel | Genomes | Total Bases | `AGC 3.2.5` (Best Single) | `AGC 3.2.5` (13-Config Oracle) | `RC-MDL` (v2) | **`HGS2` / `EIMDL` (1 Universal Config)** | **`EIMDL` Total bpb** | **Median bpb** | Wins vs `RC-MDL` | Wins vs `AGC` | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|:---:|
| ***E. coli*** | 24 | 120,108,223 | 8,384,117 B | — | 1,780,761 B | **1,177,778 B** (`-33.86%` vs RC-MDL; **`7.12x`** < AGC) | **0.07845** | **0.07185** | **24 / 24** | **24 / 24** | 24/24 PASS |
| ***M. tuberculosis*** | 24 | 105,745,015 | 292,170 B | 96,580 B | 122,998 B | **47,305 B** (`-61.54%` vs RC-MDL; **`2.04x` < AGC Oracle**) | **0.00358** | **0.00366** | **24 / 24** | **24 / 24** | 24/24 PASS |
| ***S. aureus*** | 25 | 71,043,114 | 2,979,078 B | 1,617,085 B | 267,412 B | **174,639 B** (`-34.69%` vs RC-MDL; **`9.26x` < AGC Oracle**) | **0.01967** | **0.00192** | **25 / 25** | **25 / 25** | 25/25 PASS |

### 2.2 Key Scientific Milestones Resolved in Mode R
1. ***M. tuberculosis* Oracle Barrier Broken (`47,305 B` vs `96,580 B`):** Previously (`HANDOFF.md`), `AGC 3.2.5`'s 13-configuration per-target oracle (`96,580 B`) beat `RC-MDL` (`122,998 B`) by `21.5%`. With **Horizon-Bounded Two-Gear Coalescent Haplotype Copying (Law IIIb)**, a single universal `HGS2` configuration compresses all 24 *M. tuberculosis* genomes (`105.7 Mbp`) into **`47,305 B` (`0.00358 bpb`, `2,235x`)** — **`51.02%` (`2.04x`) smaller than AGC's 13-config per-target oracle** and **`6.18x` smaller than AGC's best single config**.
2. ***E. coli* 24/24 Clean Sweep Over AGC (Including Twin Pair `13.seq` & `14.seq`):** Previously, `AGC` beat `RC-MDL` on near-identical twins `13.seq` (`1,680 B` vs `2,032 B`) and `14.seq` (`1,624 B` vs `2,266 B`). `HGS2` compresses `13.seq` (`5.185 Mbp`) to **`635 B` (`0.00098 bpb`, `8,165x`)** and `14.seq` (`5.182 Mbp`) to **`633 B` (`0.00098 bpb`, `8,187x`)**, achieving a **24/24 sweep**.
3. **Conditional Monotonicity Verified on Clonal Outbreaks:**
   - *M. tuberculosis* `03.seq` (`4,367,416 bp`): **`256 B` (`0.000469 bpb`, `17,060x`, `18` substitutions vs `RC-MDL`'s `261` mismatches and `1,635 B`)**.
   - *S. aureus* `20.seq` (`2,879,034 bp`): **`134 B` (`0.000372 bpb`, `21,485x`, `1` syntenic token, `2` substitutions, including 32-byte SHA-256 manifest digest)**.

### 2.3 Leave-Cluster-Out Unseen-Lineage Controls

Sources: [`cluster_loo_ecoli.json`](../benchmarks/fusion/results/cluster_loo_ecoli.json), [`cluster_loo_mtb.json`](../benchmarks/fusion/results/cluster_loo_mtb.json), [`cluster_loo_saureus.json`](../benchmarks/fusion/results/cluster_loo_saureus.json)

All members of the target's Mash/Jaccard cluster are excluded from the reference panel to measure generalization to unseen lineages:

| Species | Jaccard Cutoff | Clusters (Max Size) | `AGC 3.2.5 (-a)` | `RC-MDL` Cluster-LOO | **`HGS2` / `EIMDL` Cluster-LOO** | **`EIMDL` bpb** | Delta vs `RC-MDL` | Factor vs `AGC (-a)` | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| ***E. coli*** | 0.98 | 22 (max 3) | 9,203,994 B | 2,331,696 B | **1,625,424 B** | **0.10826** | **-30.29%** | **5.66x** | 24/24 PASS |
| ***M. tuberculosis*** | 0.99 | 13 (max 6) | 1,188,909 B | 153,182 B | **78,582 B** | **0.00595** | **-48.70%** | **15.13x** | 24/24 PASS |
| ***S. aureus*** | 0.98 | 14 (max 8) | 4,132,771 B | 741,650 B | **525,828 B** | **0.05921** | **-29.10%** | **7.86x** | 25/25 PASS |

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
| **`HGS2` / `EIMDL` (`--threads 1`)** | **3,377,802 B** | **0.2250** | **35.56x** | **96.9 s** | **55.4 s** | **PASS** |
| **`HGS2` / `EIMDL` (Default Parallel)** | **3,393,302 B** | **0.2260** | **35.40x** | **44.0 s** *(27.9 s cold)* | **25.0 s** *(14.0 s cold)* | **PASS** |
| **`HGS2` / `EIMDL` (`--bounded -Xmx256m`)** | **3,423,050 B** | **0.2280** | **35.09x** | **59.5 s** *(37.5 s cold)* | **29.0 s** *(16.3 s cold)* | **PASS** |
| **`HGS2` / `EIMDL` (`--bounded -Xmx192m`)** | **3,455,169 B** | **0.2301** | **34.76x** | **57.6 s** *(34.9 s cold)* | **27.2 s** *(15.6 s cold)* | **PASS** |

- `HGS2` Mode C default parallel (`3,393,302 B`) beats `HyperGraSS v1` (`3,810,409 B`) by **`-417,107 B` (`-10.95%`)** at **`20.2x–31.8x` faster compression speed**, beats `JARVIS3 l7` (`3,896,361 B`) by **`-503,059 B` (`-12.91%`)**, and beats `RC-MDL Prim` (`3,972,975 B`) by **`-579,673 B` (`-14.59%`)**.

### 3.2 Multi-FASTA Full Pangenome Collections (6-Channel Exact File Reconstruction)

| Collection | Raw File Bytes | DNA Bases | Prior Baseline | **`HGS2` / `EIMDL` Archive** | **DNA bpb** | **File Ratio** | Compress / Restore | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| *E. coli* (24 genomes, multi-FASTA) | 121,825,884 | 120,108,223 | `HyperGraSS --pangenome`: 4,885,851 B | **3,392,869 B** (`-30.56%`) | **0.2260** | **35.91x** | 63.8 s / 57.1 s | PASS |
| *M. tuberculosis* (24 genomes, multi-FASTA) | 107,257,786 | 105,745,015 | `xz -9`: ~2.3 MB | **1,041,596 B** | **0.0788** | **102.97x** | 40.0 s / 49.9 s | PASS |
| *S. aureus* (25 genomes, multi-FASTA) | 72,060,360 | 71,043,114 | `RC-MDL Prim`: 1,233,433 B / `AGC`: 2,145,694 B | **1,082,866 B** (`-12.21%` vs RC-MDL; **`1.98x`** < AGC) | **0.1219** | **66.55x** | 31.4 s / 35.0 s | PASS |

---

## 4. Amdahl Multi-Core Scaling & Strict Bounded-Heap Execution

Source JSON: [`benchmarks/fusion/results/amdahl_and_memory.json`](../benchmarks/fusion/results/amdahl_and_memory.json)

### 4.1 Multi-Core Scaling on `120,108,223 B` *E. coli* Collection

| Threads (`--threads`) | Archive Bytes | Size Delta vs `1T` | Encode Time | Encode Speedup vs `v1` (`887.3 s`) | Decode Time | Decode Speedup vs `1T` | SHA-256 |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| `1` | **3,377,802 B** | `0.00%` | 96.93 s | **9.15x** | 55.44 s | 1.00x | PASS |
| `2` | **3,380,134 B** | `+0.07%` | 62.06 s | **14.30x** | 39.83 s | **1.39x** | PASS |
| `4` | **3,388,767 B** | `+0.32%` | 51.94 s | **17.08x** | 26.02 s | **2.13x** | PASS |
| `8` | **3,393,302 B** | `+0.46%` | 51.74 s | **17.15x** | 23.27 s | **2.38x** | PASS |
| `12` | **3,393,302 B** | `+0.46%` | **44.04 s** | **20.15x** | **18.36 s** | **3.02x** | PASS |

### 4.2 Strict Bounded-Heap Execution on Pure JVM (`120,108,223 B` *E. coli* Collection)

Previously (`docs/FOUNDATIONS.md` §6), every pure-JVM bounded profile of `HyperGraSS v1` OOM-killed under a strict 256 MiB container cap (`0/7` natural collections passed), requiring a GraalVM native image (`3,864,464 B` at `-Xmx192m`) to beat `JARVIS3` (`3,896,361 B`). With `HGS2`'s streaming 64 KiB flat-file parser, incremental `concatenateAndConsume` 2-bit bank (`30.0 MB`), and 2-pass zero-materialization subsampled CSR index across **all 20 units**:

| Pure JVM Heap Cap (Encode & Decode) | Subsampling (`1 / 2^s`) | Archive Bytes | Bits / Base | Delta vs Unbounded (`3,393,302 B`) | Margin vs `JARVIS3` (`3,896,361 B`) | Margin vs `HyperGraSS v1` (`3,810,409 B`) | SHA-256 |
|---|---:|---:|---:|---:|---:|---:|:---:|
| **`-Xmx256m`** | $s = 1$ (`1/2`) | **3,423,050 B** | **0.2280** | **`+0.88%`** *(passes `<2%` gate)* | **`-12.15%`** (`-473,311 B`) | **`-10.17%`** (`-387,359 B`) | **PASS** |
| **`-Xmx224m`** | $s = 1$ (`1/2`) | **3,423,050 B** | **0.2280** | **`+0.88%`** *(passes `<2%` gate)* | **`-12.15%`** (`-473,311 B`) | **`-10.17%`** (`-387,359 B`) | **PASS** |
| **`-Xmx192m`** | $s = 2$ (`1/4`) | **3,455,169 B** | **0.2301** | **`+1.82%`** *(passes `<2%` gate)* | **`-11.32%`** (`-441,192 B`) | **`-9.32%`** (`-355,240 B`) | **PASS** |

---

## 5. Non-Existing, Xenobiological, Synthetic & Adversarial Genomes

Source JSON: [`benchmarks/fusion/results/mode_s_and_synthetic.json`](../benchmarks/fusion/results/mode_s_and_synthetic.json)

| Corpus ID | Regime | Biological / Adversarial Stress Design | Raw Bytes | `gzip -9` | **`HGS2` / `EIMDL`** | **DNA bpb** | SHA-256 |
|---|---|---|---:|---:|---:|---:|:---:|
| `01_euk_centromere` | Mode S | Higher-order alpha-satellite repeats + HOR inversions (`485 kbp`) | 491,604 B | 83,797 B | **42,193 B** | **0.6954** | PASS |
| `02_viral_rna_acgu` | Mode C | Pure RNA (`A/C/G/U`) segmented quasispecies swarm (`299 kbp`) | 304,013 B | 42,923 B | **12,122 B** | **0.3241** | PASS |
| `03_chimeric_iupac_crlf` | Mode S | Mixed DNA+RNA (`T`+`U`), all 11 IUPAC codes, soft-masking, `N`-scaffolds, ragged `CRLF` | 263,817 B | 86,566 B | **51,518 B** | **1.6887** *(1.562 raw)* | PASS |
| `04_xeno_at88` | Mode S | Extreme `88% AT` xenobiological / * Carsonella*-like genome (`500 kbp`) | 506,294 B | 123,139 B | **94,222 B** | **1.5076** | PASS |
| `05_xeno_gc86` | Mode S | Extreme `86% GC` actinomycete / synthetic genome (`500 kbp`) | 506,291 B | 126,076 B | **97,524 B** | **1.5604** | PASS |
| `06_denovo_grammar` | Mode S | Non-existing synthetic regulatory grammar + palindromic hairpins (`500 kbp`) | 506,287 B | 141,540 B | **41,096 B** | **0.6575** | PASS |
| `07_hypermutator_pangenome` | Mode C | 6-strain `6%` hypermutator + heavy indel pangenome (`1.495 Mbp`) | 1,513,506 B | 454,203 B | **101,566 B** | **0.5436** | PASS |
| `07_hypermutator_strain05` | Mode R | Single `6%` hypermutator strain vs 5-donor panel (`249 kbp`) | 249,436 B | — | **7,514 B** | **0.2410** | PASS |
| `08_crypto_white_noise` | Mode S | Cryptographic `SHA-256` counter-mode white noise (`500,000 bp`, $H_0 = 2.0000$) | 500,000 B | 144,017 B (`2.304 bpb`) | **125,075 B** | **2.0012** *($\le 2.005$ gate)* | PASS |

---

## 6. Real-Genome Law-by-Law Ablations

Source JSON: [`benchmarks/fusion/results/ablations_real_genomes.json`](../benchmarks/fusion/results/ablations_real_genomes.json)

| Ablation Variant | *E. coli* `00.seq` vs Distant Single Donor `20.seq` (`54,124` SNPs, `619` splices) | *E. coli* `00.seq` LOO (23 Donors) | *M. tuberculosis* `00.seq` LOO (23 Donors) |
|---|---:|---:|---:|
| **`full_eimdl` (All Laws I–IV Active)** | **197,112 B** (`0.3247 bpb`) | **40,735 B** (`0.0671 bpb`) | **1,237 B** (`0.00224 bpb`) |
| `no_kimura_v4_law_I` (4-symbol allele fallback) | 203,140 B (`+6,028 B` / `+3.06%`) | 41,335 B (`+600 B` / `+1.47%`) | 1,258 B (`+21 B` / `+1.70%`) |
| `no_codon_wobble_law_II` (flat gap coding) | 198,966 B (`+1,854 B` / `+0.94%`) | 40,856 B (`+121 B` / `+0.30%`) | 1,235 B (`-2 B`, only 251 SNPs) |
| `no_indel_splice_law_III` (no micro-indel splicing) | 197,205 B (`+93 B`, `+428` broken tokens) | 40,607 B (`-128 B` on dense 23-donor panel) | 1,236 B (`-1 B`, `0` splices) |
| `no_v4_no_wobble_no_splice` (all 3 disabled) | 204,912 B (`+7,800 B` / `+3.96%`) | 41,330 B (`+595 B` / `+1.46%`) | 1,255 B (`+18 B` / `+1.46%`) |
| `exact_matches_only_no_approx` (exact copy only) | 239,083 B (`+41,971 B` / **`+21.29%`**) | 43,596 B (`+2,861 B` / **`+7.02%`**) | 1,275 B (`+38 B` / **`+3.07%`**) |
| Prior `RC-MDL` v2 Baseline | 249,491 B (`+52,379 B` / **`+26.57%`**) | 73,194 B (`+32,459 B` / **`+79.68%`**) | 4,012 B (`+2,775 B` / **`+224.3%`**) |
