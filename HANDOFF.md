# HyperGraSS / RC-MDL / Unified `HGS2` (`EIMDL`) Research Handoff

Snapshot: **2026-11-14**
Repository: `https://github.com/MdSaifulIslamMSI/hypergrass-research`
Release assets: [`v2.0.0`](https://github.com/MdSaifulIslamMSI/hypergrass-research/releases/tag/v2.0.0) (`hypergrass-hgs2-v2.0.0.jar`, `HGS2_EIMDL_Research_Paper.pdf`, `hypergrass-research-v2.0.0.tar.gz`; prior `handoff-v1`: `hypergrass-handoff-v1.tar.gz`)
Research Monograph (33-Page PDF): [`paper/HGS2_EIMDL_Research_Paper.pdf`](paper/HGS2_EIMDL_Research_Paper.pdf) (`0` overfull boxes, `1,859 / 1,859` automated numerical & figure audit checks PASS via [`paper/audit_paper_numbers.py`](paper/audit_paper_numbers.py))
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

## Unified `HGS2` / `EIMDL` Final Results (2026-11-13)

`HyperGraSS` and `RC-MDL` are unified in [`HyperGraSSFusionEngine.java`](src/hypergrass/HyperGraSSFusionEngine.java) (`.hgs2` container magic `"HGS2"`, CLI [`EimdlCli.java`](src/hypergrass/EimdlCli.java), Laws I–VII + `openai2026eth` §1.1.3 OVC boundary in [`docs/COMPRESSION-LAW.md`](docs/COMPRESSION-LAW.md), full benchmark tables in [`docs/FUSION-BENCHMARKS.md`](docs/FUSION-BENCHMARKS.md) and [`docs/INDEPENDENT.md`](docs/INDEPENDENT.md)):

| Regime & Benchmark | Previous Best (`HyperGraSS v1` / `RC-MDL` / `AGC` / `JARVIS3`) | **`HGS2` / `EIMDL` Result** | Verification & Status |
|---|---|---|---|
| **Mode S (Host):** *E. coli* `00.fa` (`4,856,085 B`) | `HyperGraSS v1`: `1,136,069 B` vs `RC-MDL S`: `1,192,494 B` | **1,122,777 B (`1.8497 bpb`)** `00.fa` in **`6.75 s`** (`6.57 s` dec); `00.seq` = **1,122,749 B** | Beats both (`-1.17%` vs `v1`, `SHA-256` exact) |
| **Mode S (Independent 12-Genome Subset):** *Salmonella* (6) + *Klebsiella* (6) (`60,810,308 bp`) | `JARVIS3 -l 7`: `15,204,327 B`<br>**`JARVIS3 -l 31`**: `14,197,393 B` (`1,159 MiB` RSS)<br>`GraSS -e2`: `14,815,984 B` | **13,868,405 B (`1.8245 bpb`)** (`6,629,508 B` + `7,238,897 B`) at **`108 MiB` C / `88 MiB` R RSS** | **12/12 wins** (`-8.79%` vs L7; **`-2.32%` vs tuned L31** at **`10.7x` less RAM** & **`1.26x` faster**) |
| **Mode R LOO:** *E. coli* (24 genomes, `120.1 Mbp`) | `RC-MDL` v1/v2: `1,780,761` / `1,636,367 B`; `AGC`: `8,384,117 B` (won `13, 14`) | **1,174,113 B (`0.07820 bpb`)**; `13.seq` = **633 B**, `14.seq` = **632 B** | **24/24 clean sweep** vs RC-MDL (`-34.07%` v1 / `-28.25%` v2) & AGC (`7.14x`) |
| **Mode R LOO:** *M. tuberculosis* (24 genomes, `105.7 Mbp`) | `AGC` 13-cfg oracle: `96,580 B`; `RC-MDL` v1/v2: `122,998` / `92,506 B` | **47,339 B (`0.00358 bpb`)** (`03.seq` = **256 B**, `18` SNPs) | **24/24 vs Single & RC-MDL**; **`2.04x` < AGC 13-cfg oracle** (`19/24` wins; 5 misses `13–239 B`) |
| **Mode R LOO:** *S. aureus* (25 genomes, `71.0 Mbp`) | `RC-MDL` v1/v2: `267,412` / `243,807 B` vs `AGC` oracle: `1,617,085 B` | **174,217 B (`0.01962 bpb`)** (`20.seq` = **134 B**, `21,485x`) | **25/25 vs Single & RC-MDL**; **`9.28x` < AGC oracle** (`24/25` wins; `04.seq` `683` vs `617 B`) |
| **Mode R LOO (Independent 12-Genome Subset):** *Salmonella* + *Klebsiella* | `EIMDL` Mode S: `13,868,405 B` | **433,498 B (`0.0570 bpb`)** (`179,303 B` + `254,195 B`) | **`31.99x` smaller** than Mode S; restores in `48.4 s` (`1.20 MiB/s`) |
| **Leave-Cluster-Out:** *E. coli* / *MTB* / *S. aureus* | `RC-MDL`: `2,331,696` / `153,182` / `741,650 B` | **1,619,296** (`-30.55%`) / **78,603** (`-48.69%`) / **524,535 B** (`-29.27%`) | 73/73 byte-exact (`5.68x` / `15.13x` / `7.88x` < `AGC -a`) |
| **Mode C Collection:** *E. coli* `input.seq` (`120.1 Mbp`) | `HyperGraSS v1`: `3,810,409 B` (`887.3 s`); `JARVIS3 l7`: `3,896,361 B`; `RC-MDL`: `3,972,975 B` | **3,388,639 B (`0.2257 bpb`)** in `19.65 s` (`12T` Amdahl; `33.78 s` Mode C); `1T` = **3,371,489 B** (`0.2246 bpb`) | **`-11.52%` (`1T`) / `-11.07%` (`12T`) vs v1** at **`45.16x` speedup**; **`-13.47%` (`1T`) vs JARVIS3** |
| **Mode C (Independent Locked Trio):** *Salmonella* (25) & *Klebsiella* (23) | `JARVIS3 -l 7`: `2,709,545 B` & `2,948,926 B`<br>`GraSS -e2`: `19,991,646 B` & `21,191,566 B` | ***Salmonella***: **2,390,106 B (`0.1586 bpb`, `50.43x`)**<br>***Klebsiella***: **2,523,991 B (`0.1646 bpb`, `48.61x`)** | **`-11.79%`** & **`-14.41%` vs JARVIS3**; **`-88.04%`** & **`-88.09%` vs GraSS** |
| **Bounded Heap (`-Xmx256m` to `-Xmx128m`):** | `v1` pure JVM OOMed `0/7` at 256 MiB | **3,422,314 B** (`-Xmx256m`/`-Xmx224m`) & **3,452,703 B** (`-Xmx192m`/`-Xmx160m`/`-Xmx128m`; `25.59 s` enc / `11.67 s` dec at `128m`) | Pure JVM **`-Xmx128m`** beats unbounded `v1`, native `v1`, `JARVIS3`, & `RC-MDL` |
| **Law IV Cryptographic Safety:** `08_crypto_white_noise.seq` (`500,000 B`) | `gzip -9`: `144,017 B` (`2.3043 bpb`) | **125,075 B (`2.0012 bpb`)** | Sub-block 2-bit bypass (`+0.0012 bpb` over $H_0 = 2.0000$) |

---

## Kaggle E2E Cloud Verification & 15-Codec World Shootout (`benchmarks/kaggle/results/`)

Mandatory Kaggle-only testing policy (`test.ps1` / `test.sh` disable local execution outside Kaggle/CI; local checks compile with `javac --release 21 -Xlint:all` only). All 7 Java `Verify*` suites passed across all three cloud workers:

| Kaggle Artifact / Regime | Scope & Hardware | Key Verified `HGS2` / `EIMDL` Results | Honest External & Boundary Notes |
|---|---|---|---|
| **`world_shootout_v2_results.json`** | `19,088.08 s`, 4 EPYC vCPUs, 15 codecs | Mode C FA (*E. coli* / *Salm.* / *Kleb.*): **`3,479,510`** / **`2,417,055`** / **`2,608,659 B`**; Mode R: **`41,192`** / **`45,446`** / **`42,164 B`**; Law IV `10M` noise: **`2,500,074 B`** (`2.0001 bpb`); Pathology `15.6M`: **`2,673 B`** (`-Xmx4g` & `-Xmx128m`); Diploid `Chr21`/`Chr22`/`Chr1`: **`8,878,088`** / **`8,859,549`** / **`50,463,121 B`** | `MFCompress` GCC link fail (`unavailable`); `MBGC`/`NAF` `raw_sha256_verified=False` (`NAF` fails seq SHA on pathology; `MBGC -m3` 900s TO on `Diploid_Chr22`); `GraSS-e2` crashes (`AIOOB`) on multi-FASTA & corrupts pathology; `NAF -l 22` (`1,862,725 B`) / `MBGC -m2` (`1,966,421 B`) edge `HGS2` FA (`2,152,306 B`) on `covid_mega`; unpaired `chr1@-Xmx256m` & raw mega FA OOM |
| **`gpu_v3_largescale_results.json`** | `2x Tesla T4` GPUs (`cuda`), Human & COVID | Sketch CUDA speedup: **`43.2x–50.77x`** (`0.0` err); T4 neural probes: `949,438 B` (`1.9019 bpb`) / `926,939 B` (`1.8569 bpb`); `GRCh38_chr1_8`: **`314,875,226 B`** (`1,562,097,512 B`, `4.96x` vs. `zstd-19` `384,873,207 B`); `Diploid_Chr1_to_4`: **`189,424,047 B`** (`1,791,158,556 B`, `9.46x` vs. `AGC` `213,613,422 B` & `zstd-19` `449,440,145 B`); `covid_mega .seq@-Xmx512m`: **`1,657,005 B`** (`1,409,028,147 B`, `850.35x`, `508.74 MB` RSS) | `ChrX@-Xmx256m`, `Chr1@-Xmx256m`, `Chr1..8@-Xmx1g` OOM; `AGC` encodes `Diploid_Chr1_to_4` faster (`74.8 s` vs `671.6 s`) while `HGS2` is `-11.3%` smaller |
| **`tpu_v3_scale_results.json`** | All-24 `T2T-CHM13v2.0 \| GRCh38` & `13/13` Heap Matrix | All-24 Mode R: **`3,156,242,851 B` (`3,117,275,501 bp`) $\to$ `35,277,041 B` (`0.0894 bpb`, `89.47x`)**, `24/24` SHA-256 (`chr21`: **`938,046 B`**, `chr22`: **`1,213,472 B`**, `chr1`: **`2,589,930 B`** vs. `GeCo3-ref` `8,261,548` / `8,308,207 B`); `13/13` sweep (`100M@-Xmx128m` **`145,974 B`** .. `1000M@-Xmx512m` **`1,176,101 B`** .. `1000M@-Xmx16g` **`1,053,431 B`**); Amdahl `covid_group_22`: **`467,386 B`** (`56.26 s` 1T $\to$ `25.06 s` 4T) | Kaggle TPU batch runtime exposed 4 Intel Xeon CPUs where `jax 0.11.1` fell back to `cpu:0` (`0.27x–0.31x` vs NumPy) |
| **`world_shootout_v3_results.json`** | `29,168.35 s` (`8.1 hr`), 4 EPYC vCPUs, Java 25 + Native `eimdl`, `.hgs3` & CD | **`8/8`** `Verify*` incl. `VerifyParadigm` (`1.083 s`); Native vs Java: Bacterial Mode C (**`3,522,612`** / **`2,453,829`** / **`2,605,617 B`**), `GRCh38 chr21/22/1` (**`7,857,257`** / **`7,555,256`** / **`47,511,069 B`** @ `-61%` RSS), `Diploid_Chr21/22` (`8,846,955` / `8,808,668 B`), `10M` Law IV (`2,500,090 B`, `91.46 MB`), `1000M` Nat bounded (`2,016,887 B`), `CHM13\|GRCh38` Nat/Java (`1,018,276`/`945,989`; `1,299,810`/`1,234,617`; `3,207,153`/`2,627,685 B`); `.hgs3` *E. coli* FA **`58,987 B`** (`27,176` ORFs, `31,334` SNPs, `33,293` VCF); `CHM13 chr21` **`1,024,592 B`** (`252,280` ORFs, `285,603` SNPs, **`172,088` VCFv4.2**, `66.9 MB` RSS); Illumina FASTQ `HGS3-Java` **`6,056,142 B`** / `Native` **`6,226,185 B`** (`-17.4%` vs `NoPanel` `7,536,092 B`; `libbsc` `3,697,004 B`) | **`5/5`** Java$\leftrightarrow$Native CD parity = True; `MFCompress` GCC link `unavailable`; v3.0 FASTQ `qual=4,890,843 B` (`78.55%`), `exact_zero_target_spans=0`, and `GRCh38_chr1@-Xmx256m` OOM fixed in EIMDL v3.1 source (`ANS_MODE_RANS4_QUAL_BWT = 4`, lazy `kind >= 2`, `Packed2BitSequence` novelty streaming) |

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

**Mandatory Kaggle-Only Testing Policy:** Local execution of verification suites and benchmarks is disabled (`test.ps1` / `test.sh` delegate to Kaggle cloud / GitHub CI). Locally verify Java syntax compilation only, run the 7 `Verify*` suites on Kaggle, and run the automated numerical audit:

```powershell
javac --release 21 -Xlint:all -d build_check src/hypergrass/*.java; Remove-Item -Recurse -Force build_check
python benchmarks/kaggle/run_kaggle_e2e.py --target all
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
