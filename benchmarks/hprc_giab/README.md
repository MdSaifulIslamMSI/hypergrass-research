# HPRC Pangenome & GIAB HG002 Benchmark Harness (`benchmarks/hprc_giab/`)

Reproducible benchmark suite for evaluating **Next-Gen EIMDL v3** (`eimdl` Native + `HGS3`) against community standards on public **Human Pangenome Reference Consortium (HPRC)** multi-haplotype assemblies and **Genome in a Bottle (`GIAB` HG002 / Ashkenazi Son)** variant & Illumina FASTQ truth sets.

## 1. Tracks Measured
1. **HPRC Multi-Haplotype Assembly Mode C (`eimdl` vs. `AGC 3.2.5`, `MBGC`, `zstd-19`):**
   - Archives multi-haplotype chromosome FASTA shards (`Chr21`, `Chr22`, or full cohorts) while recording exact raw `SHA-256` round-trip verification, wall-clock encode/decode time, and `/usr/bin/time -v` Peak RSS.
2. **GIAB HG002 Illumina FASTQ Shootout + Direct Compressed-Domain VCF Concordance:**
   - Compresses `HG002` Illumina reads with `--panel <GRCh38+CHM13>` and `--cds` (`HGS3`), comparing `.hgs3` archive bytes against `libbsc`, `zstd`, `xz`, `gzip`, and `SPRING` (when present).
   - Streams `eimdl variants --archive ... --out calls.vcf` directly without decompressing FASTQ/FASTA and evaluates SNV/indel F1 concordance (using `hap.py` / `rtg vcfeval` when installed, or normalized exact-locus POS/REF/ALT matching inside GIAB high-confidence BED regions).

## 2. Usage
```bash
# Dry-run configuration & tool discovery check
python benchmarks/hprc_giab/run_hprc_giab_bench.py --dry-run

# Run against local HPRC haplotypes and/or GIAB HG002 inputs
python benchmarks/hprc_giab/run_hprc_giab_bench.py \
  --eimdl ./native/eimdl-rs/target/release/eimdl \
  --hprc-fasta /data/hprc_chr21_47haps.fa \
  --panel /data/grch38_chr21/ \
  --giab-fastq /data/hg002_chr21_reads.fastq \
  --giab-truth-vcf /data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --giab-truth-bed /data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --min-support 2 \
  --out benchmarks/hprc_giab/hprc_giab_results.json

# Offline end-to-end synthetic smoke check (runs in CI)
python benchmarks/hprc_giab/test_hprc_giab_smoke.py
```

## 3. Verified Real Public HPRC + NIST GIAB HG002 Results (`v3.2.0`)

Executed end-to-end on real public human pangenome and Illumina read data on Kaggle (`4` Intel Xeon @ `2.20 GHz` vCPUs, `32 GiB` RAM, archived in [`benchmarks/kaggle/results/hprc_giab_v3_2_results.json`](../kaggle/results/hprc_giab_v3_2_results.json)):

### HPRC `chr21` (14 Real Human Haplotypes, `470,323,902 B` Raw FASTA)

| Codec | Archive (B) | Rate (bpb) | Ratio | Enc / Dec (s) | Enc / Dec Peak RSS (MB) | SHA-256 Verified |
|---|---:|---:|---:|---:|---:|:---:|
| **`EIMDL_Native_collection` (`v3.2.0`)** | **9,959,347** | **0.169404** | **47.2244x** | 213.44 / 42.20 | **1,653.80** / **567.38** | PASS |
| `MBGC-m3` (`fef6376`) | 16,801,144 | 0.285780 | 27.9936x | 510.97 / 4.86 | 2,997.21 / 617.06 | PASS |
| `AGC` (`v3.2.5`) | 105,239,515 | 1.790077 | 4.4691x | 53.41 / — | 4,160.93 / — | — |
| `zstd-19` | 112,536,914 | 1.914203 | 4.1793x | 322.69 / 1.34 | 593.02 / 12.09 | PASS |

`EIMDL` is **`1.687x` smaller**, **`2.39x` faster to encode**, and **`1.81x` lower RSS** than `MBGC-m3`, and **`10.57x` smaller** at **`2.52x` lower RSS** than `AGC v3.2.5` (which encodes `4.00x` faster).

### NIST GIAB `HG002` `2x250bp` Illumina (`chr21:20000000-23000000`, `457,245,491 B`, `842,463` Reads, `5,753` NIST v4.2.1 Truth Variants)

| Track / Metric | Codec / Filter | Archive / Counts | Precision | Recall / Sens. | F1 | Notes / Resources |
|---|---|---|---:|---:|---:|---|
| **FASTQ Compression** | `EIMDL_Native_HGS3_CDS` | **72,240,847 B** (`1.2639 bpb`, `6.3295x`) | — | — | — | `85.38 s` enc / `19.18 s` dec (`3,488.5 MB` RSS), beats `xz-9` (`74,608,508 B`) & `zstd-19` (`77,236,258 B`), `1,011,194` ORFs |
| **FASTQ Baseline Lead** | `libbsc` | **66,748,702 B** (`1.1678 bpb`, `6.8503x`) | — | — | — | `17.16 s` enc / `8.73 s` dec (`1,284.4 MB` RSS); leads flat FASTQ size/speed |
| **CD `analyze`** | `eimdl analyze` | `2.04 s`, `406.19 MB` RSS | — | — | — | `97,275,131` panel-match bases, `1,178,003` SNPs (`Ti/Tv=0.4979`) |
| **Direct VCF (`rtg vcfeval 3.13`)** | `ms2_q10_af0.10` (`QUAL >= 18.0`) | `tp_base=5,387`, `tp_call=5,561`, `fp=617`, `fn=319` | **0.9001** | **0.9441** | **0.9216** | `20.14 s`, `1,264.14 MB` RSS (`8,696` total variants; unfiltered `F1=0.9184`) |
| **Direct VCF (BED exact)** | `ms2_q10_af0.10` | `TP=5,527`, `FP=807`, `FN=226` (`6,334` in BED) | 0.8726 | **0.9607** | 0.9145 | Zero-FASTA Phred diploid calling + GRCh38 `N`-gap `GapCoordMap` + `<= 4%` repeat guard |
| **Direct VCF (BED exact)** | `ms3_q15_af0.15` | `TP=5,472`, `FP=564`, `FN=281` (`6,036` in BED) | **0.9066** | 0.9512 | **0.9283** | `20.26 s` (`7,736` total calls: `2,616` SNV, `537` INS, `827` DEL, `3,756` complex) |

> **Honest Tool Notes:** `SPRING` failed C++ compilation from upstream source under GCC 14 (`max_readlen`), and `hap.py` was absent in the Kaggle container (`rtg vcfeval 3.13` performed normalized haplotype-aware VCF evaluation). `MFCompress` also failed GCC linking (`inline` definitions).
