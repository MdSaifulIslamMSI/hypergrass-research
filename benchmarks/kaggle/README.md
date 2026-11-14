# Kaggle-Only End-to-End Verification & World-Leading Shootouts

**Mandatory Testing Policy:** Local unit tests, verification suites, and compression benchmarks are permanently disabled (`test.ps1` / `test.sh` enforce Kaggle-only execution outside GitHub CI). Run all verification suites (`VerifyAllEdgeCases`, `VerifyDonorIndex`, `VerifyAnchorCodecs`, `VerifyArchiveAccounting`, `VerifyDonorCosts`, `VerifyCompressionLaws`, `VerifyFusionEngine`), bounded-heap stress tests (`-Xmx128m .. -Xmx1g`), and world-leading compressor shootouts end-to-end on Kaggle:

```bash
python benchmarks/kaggle/run_kaggle_e2e.py --target shootout
python benchmarks/kaggle/run_kaggle_e2e.py --target all
```

## Cloud Kernels (`mdsaifulislammsi/*`)

- **World Shootout (`world_shootout/hypergrass-eimdl-world-shootout.py`, CPU-only, 0 GPU/TPU quota):**
  Runs all 7 `Verify*` suites + adversarial pathology/Law IV noise/Law VII bounded-heap stress (`100 Mbp@-Xmx128m .. 1.0 Gbp@-Xmx512m`) + full Mode S / Mode R / Mode C shootouts comparing `HGS2 / EIMDL` and `HyperGraSS v1` against **15 open-source world leaders**: Original `GraSS` (`-e2`), `AGC 3.2.5` (`default` & `-a`), `MBGC` (`-m3`/`-m2`), `NAF` (`ennaf --level 22 --long 31`), `MFCompress` (`-3`), `GeCo3` (`-l 5` Mode S & `-r` Mode R), `JARVIS3` (`-l 7` & `-l 31`), `libbsc` (`BWT+QLFC`), `zpaq -m5`, `zstd --ultra -22 --long=31`, `zstd -19`, `xz -9e`, `brotli -q 11`, `bzip2 -9`, and `gzip -9`.
- **GPU Large Scale (`gpu_largescale/hypergrass-eimdl-gpu-scale.py`, `NvidiaTeslaT4`):**
  Single chromosomes (`GRCh38 Chr21/22/X/1`), multi-chromosome `GRCh38 Chr1..8` (`~1.48 GB`), `T2T-CHM13v2.0 | GRCh38` Mode R & Diploid Mode C (`Chr21/22/1/1..4`), SARS-CoV-2 scaling up to `1.45 GB`, and PyTorch CUDA pairwise genomic sketch similarity + causal GRU+Conv neural DNA entropy baseline.
- **TPU Scale (`tpu_scale/hypergrass-eimdl-tpu-scale.py`, `Tpu1VmV38`):**
  Amdahl scaling sweep, complete `100 Mbp .. 1.0 Gbp` `-Xmx128m .. -Xmx16g` heap sweep, all-24-chromosome `CHM13v2.0 | GRCh38` Mode R (`~6.2 Gbp` cumulative), and JAX XLA pairwise sketch matrix benchmarks.

Every run verifies lossless round-trips via SHA-256 and captures Peak RSS via `/usr/bin/time -v`.

## Committed Canonical Results (`benchmarks/kaggle/results/`)

- [`results/world_shootout_v2_results.json`](results/world_shootout_v2_results.json) — 15-codec CPU shootout (`19,088.08 s`, 4 EPYC vCPUs): all 7 `Verify*` suites PASS; `HGS2` Mode C wins bacterial multi-FASTA (*E. coli* `3,479,510 B`, *S. enterica* `2,417,055 B`, *K. pneumoniae* `2,608,659 B`) and Mode R (`41,192` / `45,446` / `42,164 B`); wins Law IV `10M` crypto noise (`2,500,074 B`, `2.0001 bpb`), `15.6M` pathology FASTA (`2,673 B` at `-Xmx4g` & `-Xmx128m`), and Diploid `Chr21`/`22`/`1` (`8,878,088` / `8,859,549` / `50,463,121 B`).
- [`results/gpu_v3_largescale_results.json`](results/gpu_v3_largescale_results.json) — `2x Tesla T4` GPU & multi-chromosome human scale: `43.2x–50.77x` CUDA sketch speedup (`0.0` max err); `GRCh38 Chr1..8` `1,562,097,512 B` $\to$ `314,875,226 B` (`4.96x` vs. `zstd-19` `384,873,207 B`); `Diploid_Chr1_to_4` `1,791,158,556 B` $\to$ `189,424,047 B` (`9.46x` vs. `AGC` `213,613,422 B`); `covid_mega .seq@-Xmx512m` `1,409,028,147 B` $\to$ `1,657,005 B` (`850.35x`, `508.74 MB` RSS).
- [`results/tpu_v3_scale_results.json`](results/tpu_v3_scale_results.json) — All-24 `T2T-CHM13v2.0 | GRCh38` Mode R (`3,156,242,851 B` / `3,117,275,501 bp` $\to$ `35,277,041 B`, `0.0894 bpb`, `89.47x`, `24/24` SHA-256) and `13/13` Law VII bounded-heap sweep (`100M@-Xmx128m` `145,974 B` .. `1000M@-Xmx512m` `1,176,101 B` .. `1000M@-Xmx16g` `1,053,431 B`).

> **Kaggle Runtime Error Visibility:** `kaggle kernels status <slug>` prints only the worker state enum (`RUNNING` / `COMPLETE` / `ERROR`). Live cell `stdout`/`stderr` streams in the Kaggle web notebook viewer during execution, while `kaggle kernels output <slug>` downloads the JSON-lines `.log` and JSON artifacts after completion or failure. Per-codec failures inside the Python benchmark harnesses (e.g., `MFCompress` GCC link failure, `Original GraSS` multi-FASTA `ArrayIndexOutOfBoundsException`, `NAF`/`GraSS` pathology SHA-256 mismatches, `MBGC` 900s timeout on `Diploid_Chr22`, and bounded human OOMs) are caught individually and recorded in JSON so the batch worker completes cleanly.
