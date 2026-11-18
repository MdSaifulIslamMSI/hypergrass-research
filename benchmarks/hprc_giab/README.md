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
  --out benchmarks/hprc_giab/hprc_giab_results.json
```
