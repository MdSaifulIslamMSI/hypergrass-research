# HyperGraSS / RC-MDL / Unified `HGS2` (`EIMDL`) Research Handoff

Snapshot: **2026-11-07**
Repository: `https://github.com/MdSaifulIslamMSI/hypergrass-research`
Release asset: `hypergrass-handoff-v1.tar.gz`
Research Monograph (22-Page PDF): [`paper/HGS2_EIMDL_Research_Paper.pdf`](paper/HGS2_EIMDL_Research_Paper.pdf) (`0` overfull boxes, `1,274 / 1,274` automated numerical audit checks PASS via [`paper/audit_paper_numbers.py`](paper/audit_paper_numbers.py))
Research Monograph (Markdown): [`docs/PAPER-DRAFT.md`](docs/PAPER-DRAFT.md)

This package is the compact, byte-exact research state for the three-generation lossless genomic, reference-conditioned, and pangenomic sequence compressor project (**Generation 1 `HyperGraSS v1`** $\to$ **Generation 2 `RC-MDL v2`** $\to$ **Generation 3 `HGS2` / `EIMDL`**). It contains all Java 21+ source files, documentation, result JSONs, manifests, publication figures, and the LaTeX/PDF monograph. The ~29 GB of raw benchmark corpus and run output is regenerable from the locked manifests and scripts; the machine-readable evidence (`benchmarks/fusion/results/*.json`, `benchmarks/independent/runs/*/*.json`, `prototype/rcmdl/work/*/*.json`, `benchmarks/corpus/runs/*/*.json`) is committed in the repository.

---

## Canonical Input & Baseline Hashes

| File | Bytes | DNA Bases | SHA-256 |
|---|---:|---:|---|
| `benchmarks/specialized/runs/full-all-methods-v1/input.seq` (masked *E. coli* 24 genomes) | 120,108,223 | 120,108,223 | `51cb3ecb2e9f532833caa0b60345f37d0342474e9f4a239801486f9bc9219acb` |
| `benchmarks/pangenome/data/pangenome.fasta` (raw *E. coli* 24-genome multi-FASTA) | 121,825,884 | 120,108,223 | `f03f0dd10c541b99ddf63d18053e3a9c6009910a81700ba18f15e24dcb7d2285` |
| `prototype/rcmdl/work/mtb/panel.fasta` (*M. tuberculosis* 24-genome multi-FASTA) | 107,257,786 | 105,745,015 | `4a62d9cf8563b01fe926d5a779dc5c44cfb3c04fd6e5ea4d07939056f905618f` |
| `prototype/rcmdl/work/saureus/panel.fasta` (*S. aureus* 25-genome multi-FASTA) | 72,060,360 | 71,043,114 | `9413289ee51e253b6adf6c389d04f20765ad1d5ac815b33d2552d9b7732adbd5` |
| `benchmarks/independent/data/salmonella/pangenome.fasta` (*S. enterica* 25 genomes, locked) | 122,251,418 | 120,527,038 | `c54a499e6172f3474dba6c363c9f05ebc764d73976ed272d97ddd01a387f26a2` |
| `benchmarks/independent/data/klebsiella/pangenome.fasta` (*K. pneumoniae* 23 genomes, locked) | 124,455,339 | 122,700,601 | `fa5d322d96e585dbc779665341c346f74eeedfc1fe888befe5e70e27da9390ae` |
| `prototype/rcmdl/vendor/agc.zip` (`AGC 3.2.5` binary archive) | 802,120 | — | `4467b4183b37affe60b0444dee255761f085c060935c2b2e29d1b7c4cd683acd` |

---

## Unified `HGS2` / `EIMDL` Final Results (2026-11-07)

`HyperGraSS` and `RC-MDL` are unified in [`HyperGraSSFusionEngine.java`](src/hypergrass/HyperGraSSFusionEngine.java) (`.hgs2` container magic `"HGS2"`, CLI [`EimdlCli.java`](src/hypergrass/EimdlCli.java), mathematical laws in [`docs/COMPRESSION-LAW.md`](docs/COMPRESSION-LAW.md), full benchmark tables in [`docs/FUSION-BENCHMARKS.md`](docs/FUSION-BENCHMARKS.md) and [`docs/INDEPENDENT.md`](docs/INDEPENDENT.md)):

| Regime & Benchmark | Previous Best (`HyperGraSS v1` / `RC-MDL` / `AGC` / `JARVIS3`) | **`HGS2` / `EIMDL` Result** | Verification & Status |
|---|---|---|---|
| **Mode S (Host):** *E. coli* `00.fa` (`4,856,085 B`) | `HyperGraSS v1`: `1,136,083 B` vs `RC-MDL S`: `1,192,494 B` | **1,123,893 B (`1.8515 bpb`)** parallel / **1,122,240 B (`1.8488 bpb`)** `1T` in `6.18 s` | Beats both (`SHA-256` exact) |
| **Mode S (Independent 12-Genome Subset):** *Salmonella* (6) + *Klebsiella* (6) (`60,810,308 bp`) | `JARVIS3 -l 7`: `15,204,327 B`<br>**`JARVIS3 -l 31`**: `14,197,393 B` (`1,159 MiB` RSS)<br>`GraSS -e2`: `14,815,984 B` | **13,868,405 B (`1.8245 bpb`)** (`6,629,508 B` + `7,238,897 B`) at **`108 MiB` C / `88 MiB` R RSS** | **12/12 wins** (`-8.79%` vs L7; **`-2.32%` vs tuned L31** at **`10.7x` less RAM** & **`1.26x` faster**) |
| **Mode R LOO:** *E. coli* (24 genomes, `120.1 Mbp`) | `RC-MDL`: `1,780,761 B` (`0.1186 bpb`); `AGC`: `8,384,117 B` (won `13, 14`) | **1,177,778 B (`0.07845 bpb`)**; `13.seq` = **635 B**, `14.seq` = **633 B** | **24/24 clean sweep** vs RC-MDL (`-33.86%`) & AGC (`7.12x`) |
| **Mode R LOO:** *M. tuberculosis* (24 genomes, `105.7 Mbp`) | `AGC` 13-config oracle: `96,580 B` beat `RC-MDL`: `122,998 B` | **47,305 B (`0.00358 bpb`)** (`03.seq` = **256 B**, `18` SNPs) | **24/24 clean sweep**: **`2.04x` smaller than AGC 13-config oracle** (`-51.02%`) with 1 universal config |
| **Mode R LOO:** *S. aureus* (25 genomes, `71.0 Mbp`) | `RC-MDL`: `267,412 B` vs `AGC` oracle: `1,617,085 B` | **174,639 B (`0.01967 bpb`)** (`20.seq` = **134 B**, `21,485x`) | **25/25 clean sweep** (`-34.69%` vs RC-MDL; **`9.26x` < AGC oracle**) |
| **Mode R LOO (Independent 12-Genome Subset):** *Salmonella* + *Klebsiella* | `EIMDL` Mode S: `13,868,405 B` | **433,498 B (`0.0570 bpb`)** (`179,303 B` + `254,195 B`) | **`31.99x` smaller** than Mode S; restores in `48.4 s` (`1.20 MiB/s`) |
| **Leave-Cluster-Out:** *E. coli* / *MTB* / *S. aureus* | `RC-MDL`: `2,331,696` / `153,182` / `741,650 B` | **1,625,424** (`-30.29%`) / **78,582** (`-48.70%`) / **525,828 B** (`-29.10%`) | 73/73 byte-exact (`5.66x` / `15.13x` / `7.86x` < `AGC -a`) |
| **Mode C Collection:** *E. coli* `input.seq` (`120.1 Mbp`) | `HyperGraSS v1`: `3,810,409 B` (`887.3 s`); `JARVIS3 l7`: `3,896,361 B`; `RC-MDL`: `3,972,975 B` | **3,393,302 B (`0.2260 bpb`)** in `44.0 s` (`27.9 s` cold); `1T` = **3,377,802 B** (`0.2250 bpb`) | **`-10.95%` vs v1** at **`20.2x–31.8x` speedup**; **`-12.91%` vs JARVIS3** |
| **Mode C (Independent Locked Trio):** *Salmonella* (25) & *Klebsiella* (23) | `JARVIS3 -l 7`: `2,709,545 B` & `2,948,926 B`<br>`GraSS -e2`: `19,991,646 B` & `21,191,566 B` | ***Salmonella***: **2,390,106 B (`0.1586 bpb`, `50.43x`)**<br>***Klebsiella***: **2,523,991 B (`0.1646 bpb`, `48.61x`)** | **`-11.79%`** & **`-14.41%` vs JARVIS3**; **`-88.04%`** & **`-88.09%` vs GraSS** |
| **Bounded Heap (`-Xmx256m` / `-Xmx224m` / `-Xmx192m`):** | `v1` pure JVM OOMed `0/7` at 256 MiB | **3,423,050 B** (`-Xmx256m`/`-Xmx224m`) & **3,455,169 B** (`-Xmx192m`) | Pure JVM `-Xmx192m` beats unbounded `v1`, native `v1`, `JARVIS3`, & `RC-MDL` |
| **Law IV Cryptographic Safety:** `08_crypto_white_noise.seq` (`500,000 B`) | `gzip -9`: `144,017 B` (`2.3043 bpb`) | **125,075 B (`2.0012 bpb`)** | Sub-block 2-bit bypass (`+0.0012 bpb` over $H_0 = 2.0000$) |

---

## Historical Pre-Fusion Baseline (2026-11-04)

| Experiment | Result | Verification |
|---|---|---|
| RC-MDL LOO, *E. coli* (24 genomes, 120,108,223 bases) | **1,780,761 B / 0.1186 bits/base** | 24/24 byte-exact |
| AGC 3.2.5 LOO, same panel (`-s 200000`) | 8,384,117 B / 0.5584 bits/base | 24/24 byte-exact |
| RC-MDL vs AGC, *E. coli* (target-00-selected setting) | **4.71x smaller** | median 4.09x; AGC wins only targets 13, 14 |
| RC-MDL LOO, *M. tuberculosis* (24 genomes, 105,745,015 bases) | **122,998 B / 0.009305 bits/base** | 24/24 byte-exact |
| AGC 3.2.5 LOO, *M. tuberculosis* (full 13-config sweep) | 292,170 B best single / 96,580 B oracle | 2.38x larger / **0.79x** (AGC oracle 21.5% smaller than Gen 2 RC-MDL) |
| Self-contained single genome (*E. coli* 00) | HyperGraSS 1,136,083 B vs RC-MDL 1,192,494 B | HyperGraSS v1 won by 4.7% |
| Progressive collection (24 *E. coli*) | HyperGraSS 3,810,409 B vs RC-MDL Prim **3,972,975 B** | HyperGraSS v1 won by 4.3% |
| Panel-size ablation (target 00) | 1/2/4/8/16/23 genomes $\to$ 249,491 / 185,439 / 146,758 / 102,931 / 87,254 / 73,194 B | all byte-exact |
| *S. aureus* LOO (25 genomes, 71,043,355 bases) | **267,412 B / 0.03011 bits/base** | 25/25 byte-exact |
| AGC 3.2.5 LOO, same *S. aureus* panel | 3,539,729 B (`-a`) / 2,979,078 B (best single) / 1,617,085 B (oracle) | 13.24x / 11.14x / **6.05x** larger than RC-MDL |
| Leave-cluster-out, *E. coli* (cluster `13,14,15` excluded) | **2,331,696 B / 0.15531 bits/base** vs AGC 9,203,994 B (`-a`) | **3.95x**; RC wins 24/24 |
| Leave-cluster-out, *M. tuberculosis* (13 clusters, largest 6) | **153,182 B / 0.01159 bits/base** vs AGC 1,188,909 B (`-a`) | **7.76x**; RC wins 16/24 |
| Leave-cluster-out, *S. aureus* (14 clusters, largest 8) | **741,650 B / 0.08352 bits/base** vs AGC 4,132,771 B (`-a`) | **5.57x**; RC wins 23/25 |
| Collection order: cluster-aware Prim (24 *E. coli*) | 3,982,813 B | +0.25% vs Prim; still +4.5% vs HyperGraSS v1 |
| Collection order: greedy measured marginal | 3,990,249 B | +0.43% vs Prim |
| Collection order: fixed anchors 4/8/12/16/20 | 6,621,279 / 5,818,621 / 5,494,235 / 4,921,615 / 4,115,368 B | all worse than Prim |
| Collection, *S. aureus* (25 genomes) | RC-MDL Prim **1,233,433 B** vs AGC full archive 2,145,694 B | 1.74x smaller |

---

## Restore from the Release Asset & Verify Everything

```powershell
git clone https://github.com/MdSaifulIslamMSI/hypergrass-research.git
cd hypergrass-research
gh release download handoff-v1 --pattern 'hypergrass-handoff-v1.tar.gz'
tar -xzf hypergrass-handoff-v1.tar.gz
```

Verify every file against `MANIFEST.sha256`:

```powershell
Get-Content MANIFEST.sha256 | ForEach-Object {
  $hash, $path = $_ -split '  ', 2
  $actual = (Get-FileHash $path -Algorithm SHA256).Hash.ToLower()
  if ($actual -ne $hash) { throw "mismatch: $path" }
}
```

Rebuild the per-genome panel directories:

```powershell
python prototype/rcmdl/split_panel.py --fasta benchmarks/pangenome/data/pangenome.fasta --out prototype/rcmdl/work/genomes
python prototype/rcmdl/split_input.py --input benchmarks/specialized/runs/full-all-methods-v1/input.seq --manifest benchmarks/pangenome/manifest.json --out prototype/rcmdl/work/genomes_masked
python prototype/rcmdl/split_panel.py --fasta prototype/rcmdl/work/mtb/panel.fasta --out prototype/rcmdl/work/mtb/genomes
python prototype/rcmdl/split_panel.py --fasta prototype/rcmdl/work/saureus/panel.fasta --out prototype/rcmdl/work/saureus/genomes
```

Compile the unified `HGS2` / `EIMDL` engine, run the verification suites, and execute the `1,274`-assertion numerical audit:

```powershell
New-Item -ItemType Directory -Force -Path bin | Out-Null
javac --release 21 -d bin src/hypergrass/*.java
java -cp bin hypergrass.VerifyCompressionLaws
java -cp bin hypergrass.VerifyFusionEngine
python paper/audit_paper_numbers.py
```

---

## Regenerating Bulk Benchmark Panels

- ***E. coli* panel:** `python benchmarks/pangenome/fetch.py` (`benchmarks/pangenome/manifest.json`).
- ***M. tuberculosis* & *S. aureus* panels:** `python prototype/rcmdl/fetch_panel.py`.
- **Independent locked *Salmonella* & *Klebsiella* panels:**
  ```powershell
  python -B benchmarks/independent/fetch_locked.py --dataset salmonella
  python -B benchmarks/independent/fetch_locked.py --dataset klebsiella
  python -B benchmarks/independent/run_trio.py --datasets salmonella klebsiella --run-id trio-1
  python -B benchmarks/independent/run_modes.py --datasets salmonella klebsiella --genomes 6 --run-id modes-subset-1
  python -B benchmarks/independent/run_jarvis_level.py --scope subset --level 31 --genomes 6 --run-id jarvis-l31-subset-1
  ```
- **Unified `HGS2` / `EIMDL` benchmark suite:** `python benchmarks/fusion/run_fusion_benchmarks.py`.
