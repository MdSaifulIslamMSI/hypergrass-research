# EIMDL Research Launch & Adoption Playbook

This step-by-step launch guide turns **HyperGraSS / RC-MDL / EIMDL** (`HGS2` / `.hgs3` + Lean 4 proofs) from a verified open repository into a cited, package-installed community standard across **bioRxiv**, **Zenodo DOI**, **Bioconda**, **Nextflow/Snakemake**, and **computational genomics / formal-methods communities**.

### Current Live Launch Status (`2026-11-18`)

| Track | Identifier / URL | Status |
|---|---|---|
| **Author ORCID** | [`0009-0006-0570-9934`](https://orcid.org/0009-0006-0570-9934) (`Academy of Technology`) | Linked in repo, papers, Zenodo & bioRxiv |
| **Zenodo Permanent Archive** | Concept DOI [`10.5281/zenodo.24344825`](https://doi.org/10.5281/zenodo.24344825) (`v3.1.1`: [`10.5281/zenodo.24344826`](https://doi.org/10.5281/zenodo.24344826)) | Live at CERN/Zenodo + `CITATION.cff` / `README.md` / PDFs |
| **Bioconda Recipe PR** | [`bioconda/bioconda-recipes#71846`](https://github.com/bioconda/bioconda-recipes/pull/71846) (`eimdl 3.1.0`) | `Lint`, `Linux`, `OSX-64`, `ARM` PASS · `please review & merge` |
| **bioRxiv Preprint** | `BIORXIV/2026/767583` (`Version 1`, submitted `2026-11-18`) | Submitted (`7p` main + `33p` supp PDF); awaiting screening & DOI |

---

## 1. Pre-Launch Checklist (Verify Green CI First)

1. Confirm `master` is green across all GitHub Actions jobs in [`.github/workflows/ci.yml`](../.github/workflows/ci.yml):
   - Java verification matrix (`8/8` suites including `VerifyCompressionLaws`)
   - Native Linux `x86_64` & macOS Apple Silicon (`arm64`) C23/C++23/Rust 2024 checks (`cargo test --offline`)
   - Offline Bioconda package verification (`packaging/bioconda/eimdl/build.sh` + `run_test.sh`)
   - **Formal Lean 4 Proof Verification (`proofs/lean/` — zero `sorry`)**
   - Paper number & PDF audit (`paper/audit_paper_numbers.py`)
2. Inspect all three compiled manuscripts before submission:
   - Concise main paper (~7 pages, journal/bioRxiv target): [`paper/EIMDL_Main_Manuscript.pdf`](../paper/EIMDL_Main_Manuscript.pdf)
   - Standalone Lean 4 formal-methods pearl (~4 pages, CPP/ITP/AFP target): [`paper/EIMDL_Lean4_Paper.pdf`](../paper/EIMDL_Lean4_Paper.pdf)
   - Exhaustive technical monograph & 73-genome appendices (~33 pages): [`paper/HGS2_EIMDL_Research_Paper.pdf`](../paper/HGS2_EIMDL_Research_Paper.pdf)

---

## 2. Mint Permanent Archival DOI via Zenodo (Do First!)

1. Sign in to [Zenodo](https://zenodo.org) with your GitHub account (`MdSaifulIslamMSI`).
2. Go to **GitHub -> Enable repository** and toggle ON **`MdSaifulIslamMSI/hypergrass-research`**.
3. Wait for the **Kaggle Version #7 (`v3.1.0`)** tag or click **Create a new release** on GitHub. Zenodo automatically reads [`.zenodo.json`](../.zenodo.json) and mints a permanent versioned + concept DOI (e.g., `10.5281/zenodo.XXXXXXX`).
4. Paste the minted DOI into [`CITATION.cff`](../CITATION.cff), [`README.md`](../README.md), and the bioRxiv submission header.

---

## 3. Submit Preprint to bioRxiv (Subject: *Bioinformatics*)

Submit [`paper/EIMDL_Main_Manuscript.pdf`](../paper/EIMDL_Main_Manuscript.pdf) as the **Main Manuscript** and [`paper/HGS2_EIMDL_Research_Paper.pdf`](../paper/HGS2_EIMDL_Research_Paper.pdf) as **Supplementary Data / Comprehensive Monograph**.

**Ready-to-paste bioRxiv Metadata & Abstract:**
- **Subject Area:** Bioinformatics (Secondary: Genomics)
- **Title:** *Evolutionary Information Minimum Description Length (EIMDL): Formally Verified Evolutionary Coding, Sub-128 MiB Pangenome Archiving, and Compressed-Domain Genomics*
- **License:** CC BY 4.0 (Preprint) / Apache-2.0 (Software & Lean Proofs)
- **Code & Formal Proofs Links:**
  - Repository: `https://github.com/MdSaifulIslamMSI/hypergrass-research`
  - Machine-Checked Lean 4 Laws: `https://github.com/MdSaifulIslamMSI/hypergrass-research/tree/master/proofs/lean`
  - Archival Release & Binaries: `https://github.com/MdSaifulIslamMSI/hypergrass-research/releases`

**Suggested Target Peer-Reviewed Journals (Post-Preprint):**
1. **Primary:** *Bioinformatics* (Oxford Academic — Original Paper or Applications Note; fits concise main manuscript + formal theorems + CLI tool) or *Nucleic Acids Research* (*NAR* / *NAR Genomics and Bioinformatics*).
2. **Methods / Pangenomics Alternative:** *Genome Research* or *Algorithms for Molecular Biology* (strong fit for formal proofs + finite-state evolutionary transducers).
3. **Formal Verification Companion Note:** *Archive of Formal Proofs* / *CPP* / *ITP* short pearl showcasing SAT-certified 64-bit SWAR DNA word involutions (`BitVec 64`) and range-coder invariant proofs (`[11, 4085]`) in Lean 4 (`proofs/lean/`).

---

## 4. Open Upstream `bioconda-recipes` Pull Request

Once the release tarball (`hypergrass-research-v<TAG>.tar.gz` or GitHub source archive) is tagged:

1. Compute the exact SHA-256 of the GitHub release source archive:
   ```bash
   curl -sL https://github.com/MdSaifulIslamMSI/hypergrass-research/archive/refs/tags/v3.1.0.tar.gz | sha256sum
   ```
2. Update `packaging/bioconda/eimdl/meta.yaml` (`source.url` + `source.sha256` pointing to the GitHub archive instead of local path `../../..`).
3. Fork `bioconda/bioconda-recipes`, copy `packaging/bioconda/eimdl/` to `recipes/eimdl/`, run `bioconda-utils lint recipes/eimdl`, and submit the PR following [`packaging/bioconda/README.md`](../packaging/bioconda/README.md).
4. Mention that the identical `build.sh` and `run_test.sh` already pass end-to-end under GitHub Actions Linux CI.

---

## 5. Run Multi-Haplotype Human Benchmarks (HPRC & GIAB HG002)

On a Linux cluster node (`>= 16 GB RAM`, optional `agc`/`mbgc`/`spring` installed via Bioconda):

```bash
conda install -c bioconda -c conda-forge eimdl agc mbgc spring zstd
python benchmarks/hprc_giab/run_hprc_giab_bench.py \
    --hprc-dir /path/to/hprc_chr21_haplotypes \
    --fastq /path/to/HG002_R1.fastq \
    --ref /path/to/GRCh38_chr21.fa \
    --vcf /path/to/HG002_GRCh38.benchmark.vcf.gz \
    --out benchmarks/hprc_giab/hprc_hg002_report.json
```

See [`benchmarks/hprc_giab/README.md`](../benchmarks/hprc_giab/README.md) for VCF concordance metrics (precision / recall / F1 against GIAB truth sets).

---

## 6. Community Announcement & Author Outreach Templates

### A. Show HN / Reddit `r/bioinformatics` / Biostars / Bluesky & X Thread
> **Headline:** Show HN: EIMDL – Lossless genome & FASTQ compressor with Lean 4–verified compression laws, sub-128 MiB JVM streaming, and native Rust/C++23/C23 CLI (`eimdl`)
>
> We unified biological coding-frame symmetries (Klein $V_4$ strand-invariant SNPs, codon-wobble bijection, 6-frame ORF-HMM, two-adic $\tau=20$ range adaptation, and two-gear coalescent DAG copying) into an MIT/Apache open engine spanning pure Java 25 (`.hgs2`/`.hgs3`), a native `C23/C++23/Rust 2024` CLI (`eimdl`), and zero-`sorry` **Lean 4 formal proofs** verifying the group algebra, 64-bit SWAR bit-tricks, DAG acyclicity, and range-coder bounds.
>
> **Results on 73 SHA-256 verified bacterial genomes + human chr21/chr1:**
> - Mode R LOO vs `AGC 3.2.5`: **7.14x smaller** on *E. coli* (`1,174,113 B` vs `8.38 MB`), **6.17x smaller** than single-cfg `AGC` (`2.04x` smaller than a 13-cfg oracle) on *M. tuberculosis*, **17.10x smaller** on *S. aureus*.
> - Runs pangenome Mode C (`120.1 MB` *E. coli*) inside `-Xmx128m` (`3.45 MB` archive) and adds compressed-domain `eimdl variants` / `eimdl analyze` (`VCF v4.2` directly from compressed `.hgs3` without writing FASTA to disk).
> - Preprint + 33-page monograph + Lean 4 CI proofs + Bioconda/Nextflow/Snakemake/Docker wrappers all included!

### B. Concise Benchmark Verification Note to Baseline Authors (AGC / MBGC / PanMAN / JARVIS3 / GraSS)
Send *only after* Zenodo DOI + bioRxiv preprint link + CI badges are live:

> **Subject:** Reproducible benchmark + formal Lean 4 laws for unaligned pangenome compression (`EIMDL` / `HGS2` / `.hgs3`)
>
> Dear Prof./Dr. [Name],
>
> Building directly on your insights in [AGC / JARVIS3 / MBGC / GraSS / PanMAN], we recently released **EIMDL (`HGS2` / `.hgs3`)**, an open-source Apache-2.0 multi-regime evolutionary archiver and compressed-domain caller (`https://github.com/MdSaifulIslamMSI/hypergrass-research`, Zenodo DOI: `10.5281/zenodo.24344825`, bioRxiv DOI: `[insert live 10.1101/... DOI from BIORXIV/2026/767583]`).
>
> We benchmarked against your tool across 73 Mode R Leave-One-Out bacterial genomes plus locked Docker *Salmonella*/*Klebsiella* and Kaggle *Homo sapiens* chr1/chr21 panels, preserving full reproducible scripts, SHA-256 checks, and an explicit *Honest Boundaries* section noting where existing tools remain faster or smaller (e.g., ultra-clonal viral FASTA cohorts and single-genome `JARVIS3 -l 31` rates). We also mechanized the core algebraic and range-coder invariants in **Lean 4** (`proofs/lean/`).
>
> If you spot any configuration we should re-evaluate for your tool or wish to suggest additional panels (e.g., HPRC haplotypes via `benchmarks/hprc_giab/`), we would be grateful to incorporate them and run an updated comparison.
>
> Best regards,
> Md Saiful Islam
