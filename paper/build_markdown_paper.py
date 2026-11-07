#!/usr/bin/env python3
"""
Generate the complete, audited end-to-end Markdown research monograph:
  1. docs/PAPER-DRAFT.md (Repository Markdown monograph)
  2. paper/artifact_preview.md (Source for User-facing artifact HGS2_EIMDL_Research_Paper.md)
Covers 100% of Generation 1 (HyperGraSS v1), Generation 2 (RC-MDL v2), and Generation 3 (HGS2 / EIMDL).
"""

import os
from table_builders import build_all_tables

BASE_DIR = r"C:\Users\mdsai\.gemini\antigravity\scratch\hypergrass"
DOCS_DIR = os.path.join(BASE_DIR, "docs")
ARTIFACT_DIR = r"C:\Users\mdsai\.gemini\antigravity\brain\8fd6ada9-0a0f-403f-95e4-ee038b5c4fa3"

MD_PART1 = r"""# Evolutionary Information Minimum Description Length (EIMDL): A Unified Five-Law Theory and Bounded-Memory Architecture for Single-Genome, Reference-Conditioned, and Pangenome Compression

**Author:** Md Saiful Islam  
**Engine & Container Specification:** HyperGraSS v2 (`HGS2` / `EIMDL`, `.hgs2` v2.0) & HyperGraSS v1 (`.hgs` v1.0)  
**Repository:** [https://github.com/MdSaifulIslamMSI/hypergrass-research](https://github.com/MdSaifulIslamMSI/hypergrass-research)  
**Compiled 22-Page LaTeX Monograph (PDF):** [HGS2_EIMDL_Research_Paper.pdf](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/paper/HGS2_EIMDL_Research_Paper.pdf) · **LaTeX Source:** [HGS2_EIMDL_Research_Paper.tex](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/paper/HGS2_EIMDL_Research_Paper.tex)

---

## Abstract

Lossless genomic sequence compression has historically fractured into three incompatible algorithmic silos:
1. **Alignment- and anchor-dependent pangenome archivers** such as `AGC 3.2.5`, which require conserved $k$-mer anchor segmentation and suffer up to $20\times$ archive inflation when structural rearrangements, mobile elements, or elevated substitution rates disrupt anchor colinearity;
2. **Exact-match Lempel–Ziv and Burrows–Wheeler codecs** (`xz -9`, `bzip2 -9`, `zstd -19`, `gzip -9`, and the original `GraSS -e2`), which shatter syntenic homologies whenever synonymous third-codon point mutations interrupt exact substrings; and
3. **Pure finite-context Markov and weighted-repeat models** (`GeCo3`, `JARVIS3`), which model local nucleotide statistics at substantial memory cost (up to `1,158 MiB` for `JARVIS3 -l 31`) without capturing multi-donor coalescent haplotype switching, affine syntenic micro-indels, or 6-channel FASTA container formatting.

This monograph presents the complete, end-to-end theoretical and experimental record of the **HyperGraSS** research program across all three engineering generations—**Generation 1 (`HyperGraSS v1`)**, **Generation 2 (`RC-MDL v2`)**, and **Generation 3 (`HGS2` / `EIMDL`)** ([HyperGraSSFusionEngine.java](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/HyperGraSSFusionEngine.java), [FastNoveltyCoder.java](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/FastNoveltyCoder.java), [FastaStreamFactorizer.java](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/FastaStreamFactorizer.java), [EimdlCli.java](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/EimdlCli.java), [Main.java](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/Main.java))—governed by five mathematically proved, experimentally falsifiable laws:
- **Law I (Klein Four-Group $V_4 = (\Sigma_4, \oplus)$ Strand-Invariant Kimura–CpG Differential Substitution Coding):** Proving that bitwise XOR differentials under $\phi(\texttt{A})=0, \phi(\texttt{C})=1, \phi(\texttt{T})=2, \phi(\texttt{G})=3$ commute identically with Watson–Crick reverse-complement inversion ($\bar{b} = b \oplus 2$), algebraically separate transitions from transversions on both strands simultaneously, and reduce substitution entropy from $\log_2 3 \approx 1.5850\text{ bits}$ to $1.10\text{--}1.32\text{ bits/SNP}$.
- **Law II (Bijective Codon-Wobble Gap Decomposition $\psi(\text{gap}) = (\lfloor \text{gap}/3 \rfloor, (\text{gap}+1)\bmod 3)$ & 6-State Causal Open-Reading-Frame HMM):** Proving zero description-length overhead on non-coding DNA while capturing synonymous third-codon periodicity ($\omega = d_N/d_S \ll 1$) in aligned matches and causal 6-frame stop-codon depletion in unaligned novelty sequences (`6F-ORF-HMM`).
- **Law III (Conditional Monotonicity, Horizon-Bounded Two-Gear Coalescent Haplotype Copying & Affine Micro-Indel Splicing):** Combining a strict clonal copying gear (`+0.008` bits/bp, `-11.2` bits/SNP, $\ge 128\text{ bp}$ tail guard) with a horizon-clamped diverged syntenic gear (`findClonalReentryHorizon` at $\ge 192\text{ bp}$ lookahead) and 2D affine micro-indel splicing ($\Delta_{\text{ref}} \in [1,256]$, $\Delta_{\text{tgt}} \in [1,64]$) to eliminate recombination bulldozing (`29,623` false SNPs / `56,813 B` $\to$ `5,706` true SNPs / `40,735 B` on *E. coli* `00.seq`).
- **Law IV (Universal Five-Regime Existing & Non-Existing Sequence Completeness):** Guaranteeing a strict worst-case safety ceiling of $\le 2.0012\text{ bits/base}$ ($+0.06\%$ over raw 2-bit packing) on cryptographic white noise via per-sub-block raw 2-bit bypass (`raw2BitBypass`).
- **Law V (Generalized Amdahl-Memory Streaming Scaling):** Uniting a $0.25\text{ B/base}$ packed 2-bit DNA bank with a two-pass zero-materialization canonical open-syncmer Compressed Sparse Row (`CSR`) index ($k=16, s=10, w=6, t=3$) of exact size $N_{\text{panel}}(0.25 + 2^{-s}) + 4\text{ MiB}$, enabling full `120.11 Mbp` pangenome compression inside a strict **`-Xmx192m` (`192 MiB`) pure-JVM heap cap** (`3,455,169 B`) and scaling across 12 threads (`96.93 s` $\to$ `44.04 s` encode [$2.20\times$]; `55.44 s` $\to$ `18.36 s` decode [$3.02\times$]).

![Figure 1: Multi-Regime Description Length and Information Decomposition](%%FIG_PREFIX%%/fig1_regime_overview.png)

---

## 1. Introduction, Three-Generation Evolution & Scientific Contract

### 1.1 Executive Empirical Summary Across All Canonical Regimes

Every number in this paper is backed by a committed machine-readable JSON report in [`benchmarks/`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks) or [`prototype/rcmdl/work/`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/prototype/rcmdl/work) and verified by [`paper/audit_paper_numbers.py`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/paper/audit_paper_numbers.py).

| Regime | Benchmark Corpus | Total Bases (`bp`) | `HGS2` / `EIMDL` Archive (`B`) | Rate (`bpb`) | Ratio | Strongest External Baseline (`B`) | `HGS2` Advantage & Win Record |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| **Mode R (LOO)** | *E. coli* (24 genomes) | `120,108,223` | **`1,177,778`** | **`0.07845`** | **`101.98x`** | `AGC 3.2.5 -s 200k`: `8,384,117` | **`7.12x` smaller (`24/24` wins)** |
| **Mode R (LOO)** | *M. tuberculosis* (24 genomes) | `105,745,015` | **`47,305`** | **`0.00358`** | **`2,235.39x`** | `AGC Single`: `292,170` · `Oracle`: `96,580` | **`6.18x` vs Single (`24/24`); `2.04x` vs 13-Cfg Oracle (`19/24`)** |
| **Mode R (LOO)** | *S. aureus* (25 genomes) | `71,043,114` | **`174,639`** | **`0.01967`** | **`406.80x`** | `AGC Single`: `2,979,078` · `Oracle`: `1,617,085` | **`17.06x` vs Single (`25/25`); `9.26x` vs 13-Cfg Oracle (`24/25`)** |
| **Mode R (LOO)** | *S. enterica* + *K. pneumoniae* (12g subset) | `60,810,308` | **`433,498`** | **`0.05703`** | **`140.28x`** | `EIMDL Mode S`: `13,868,405` | **`31.99x` smaller than reference-free Mode S** |
| **Mode R (Cluster-LOO)** | *E. coli* ($J=0.98$, 22 clusters) | `120,108,223` | **`1,625,424`** | **`0.10826`** | **`73.89x`** | `AGC 3.2.5 -a`: `9,193,975` | **`5.66x` smaller on unseen lineages** |
| **Mode R (Cluster-LOO)** | *M. tuberculosis* ($J=0.99$, 13 clusters) | `105,745,015` | **`78,582`** | **`0.00595`** | **`1,345.66x`** | `AGC 3.2.5 -a`: `1,189,036` | **`15.13x` smaller on unseen lineages** |
| **Mode R (Cluster-LOO)** | *S. aureus* ($J=0.98$, 14 clusters) | `71,043,114` | **`525,828`** | **`0.05921`** | **`135.11x`** | `AGC 3.2.5 -a`: `4,132,331` | **`7.86x` smaller on unseen lineages** |
| **Mode C (Collection)** | *E. coli* 24g Flat (`--threads 1` / `12T`) | `120,108,223` | **`3,377,802` / `3,393,302`** | **`0.2250` / `0.2260`** | **`35.56x` / `35.40x`** | `JARVIS3 -l 7`: `3,896,361` · `GraSS`: `24,777,592` | **`-13.31%` vs `JARVIS3`; `7.34x` vs `GraSS`** |
| **Mode C (Collection)** | *E. coli* 24g Bounded Heap (`-Xmx192m`) | `120,108,223` | **`3,455,169`** | **`0.2301`** | **`34.76x`** | `HyperGraSS v1`: `3,810,409` (`1,452 MB` RSS) | **`-9.32%` smaller inside `192 MiB` pure-JVM cap** |
| **Mode C (Collection)** | *E. coli* 24g Multi-FASTA | `120,108,223` | **`3,392,869`** | **`0.2260`** | **`35.91x`** | `AGC 3.2.5 -a`: `8,540,341` | **`2.52x` smaller (6-channel exact FASTA)** |
| **Mode C (Collection)** | *M. tuberculosis* 24g Multi-FASTA | `105,745,015` | **`1,041,596`** | **`0.0788`** | **`102.97x`** | `RC-MDL v2`: `1,148,486` · `AGC`: `3,318,508` | **`-9.31%` vs `RC-MDL`; `3.19x` vs `AGC` (exact FASTA)** |
| **Mode C (Collection)** | *S. aureus* 25g Multi-FASTA | `71,043,114` | **`1,082,866`** | **`0.1219`** | **`66.55x`** | `RC-MDL v2`: `1,253,650` · `AGC`: `5,637,705` | **`-13.62%` vs `RC-MDL`; `5.21x` vs `AGC` (exact FASTA)** |
| **Mode C (Locked Trio)** | *Salmonella enterica* (25 genomes) | `120,527,038` | **`2,390,106`** | **`0.1586`** | **`50.43x`** | `JARVIS3 -l 7`: `2,709,545` · `GraSS`: `19,991,646` | **`-11.79%` vs `JARVIS3`; `-88.04%` (`8.36x`) vs `GraSS`** |
| **Mode C (Locked Trio)** | *Klebsiella pneumoniae* (23 genomes) | `122,695,431` | **`2,523,991`** | **`0.1646`** | **`48.61x`** | `JARVIS3 -l 7`: `2,948,926` · `GraSS`: `21,191,566` | **`-14.41%` vs `JARVIS3`; `-88.09%` (`8.40x`) vs `GraSS`** |
| **Mode S (Single)** | *S. enterica* + *K. pneumoniae* (12g subset) | `60,810,308` | **`13,868,405`** | **`1.8245`** | **`4.38x`** | `JARVIS3 -l 31`: `14,197,393` · `-l 7`: `15,204,327` | **`-2.32%` vs `-l 31` (`10.7x` less RAM, `12/12` wins)** |
| **Law IV Safety** | `08_crypto_white_noise.seq` (SHA-256 CTR) | `500,000` | **`125,075`** | **`2.0012`** | **`4.00x`** | `gzip -9`: `144,017` (`2.3043 bpb`) | **`+0.06%` above `2.0000 bpb` Shannon limit** |

### 1.2 Three-Generation Architectural Evolution & Resolution of Naming Debt

To maintain complete historical and technical clarity across the codebase and documentation:
1. **Original `GraSS` (`AnshuOishik/GraSS -e2`):** External two-stage LZ-style match finder coupled with `libbsc` block-sorting. On our canonical `120,108,223 bp` 24-genome *E. coli* collection, `GraSS -e2` produces `24,777,592 B` (`1.6503 bpb`, `4.85x`).
2. **Generation 1 — `HyperGraSS v1` ([`Main.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/Main.java), `.hgs` v1.0 container, [`docs/REPORT.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/REPORT.md)):** Our first-generation multi-block archiver combining open-syncmer LZ77, divsufsort Burrows–Wheeler Transform (`BWT`) + `rANS` entropy coding (`--fast`/`--lean`), an anchor-bank pangenome differencing engine (`--pangenome`) with 4-donor anchor bank + bit-packed mismatch mask + order-6 Context-Repeat Mixer (`CRM`), and a single-genome quadratic CRM (`--single`, [`SingleGenomeCoder.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/SingleGenomeCoder.java)). Generation 1 reduced the 24-genome *E. coli* archive from `24,777,592 B` to `3,810,409 B` (`0.2538 bpb`, `31.52x`), passed a strict `256 MiB` Linux cgroup v2 gate (`3,864,464 B` at `251.1 MiB` peak via GraalVM 25 static-musl native compilation), scaled to `2.077 GiB` (`2,230,000,000 bp`) via 64-bit `SequenceSource` streaming (`30,171,539 B` at `235 MB` RSS), and mathematically proved via an 8-channel information audit ([`docs/ENTROPY-FLOOR.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/ENTROPY-FLOOR.md)) why fixed-anchor block architectures hit a hard structural wall above `1.20 MB` (`100:1`).
3. **Generation 2 — `RC-MDL v2` ([`prototype/rcmdl`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/prototype/rcmdl), [`prototype/rcmdl/RESULTS.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/prototype/rcmdl/RESULTS.md)):** Our standalone Reference-Conditioned Minimum Description Length prototype, introducing Mode R (unaligned donor-panel conditioning) and Mode C (Containment/Prim tree-ordered pangenome coding without fixed anchor replication). `RC-MDL v2` broke the Generation 1 anchor floor in Mode R (`1,636,367 B` on 24g *E. coli* LOO; `92,506 B` on 24g *M. tuberculosis* LOO; `243,807 B` on 25g *S. aureus* LOO) and reached `3,972,975 B` in Mode C, while executing a 13-configuration parameter sweep of `AGC 3.2.5`.
4. **Generation 3 — Unified `HGS2` / `EIMDL` ([`HyperGraSSFusionEngine.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/HyperGraSSFusionEngine.java), [`FastNoveltyCoder.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/FastNoveltyCoder.java), [`FastaStreamFactorizer.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/FastaStreamFactorizer.java), [`EimdlCli.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/EimdlCli.java), `.hgs2` v2.0 container):** The unified five-law architecture that replaces Generation 1's slow `CRM` (`887.45 s`) and Generation 2's hash lists (`1,812 MB` RSS) with the five fundamental laws of EIMDL, achieving **`1,177,778 B` (`101.98x`)** on 24g *E. coli* Mode R LOO, **`47,305 B` (`2,235.39x`)** on 24g *M. tuberculosis* LOO, **`174,639 B` (`406.80x`)** on 25g *S. aureus* LOO, **`3,377,802 B` (`35.56x`)** on 24g *E. coli* Mode C (`44.04 s` encode / `18.36 s` decode at 12T), and **`3,455,169 B` (`34.76x`)** inside a **`192 MiB` (`-Xmx192m`) pure-JVM heap**.

### 1.3 Four-Part Scientific Contract ([`docs/SCIENTIFIC-CONTRACT.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/SCIENTIFIC-CONTRACT.md))
1. **Regime Hygiene:** Never compare reference-conditioned Mode R byte counts against self-contained Mode C or reference-free Mode S tools without explicit labeling; every Mode R archive embeds a mandatory **32-byte SHA-256 digest** of the canonical donor panel and verifies it prior to decompression.
2. **Baseline Fairness:** Tune external baselines across their full parameter spaces (13 configurations of `AGC 3.2.5` plus a 13-configuration per-target oracle; both default `-l 7` and maximum high-memory `-l 31` for `JARVIS3 v3.7`; both `-l 3` and `-l 6` for `GeCo3`).
3. **Container Integrity:** Preserve $100\%$ byte-for-byte SHA-256 identity across raw `.seq`, single-FASTA, and multi-FASTA files (headers, line widths, `LF`/`CRLF`, lowercase soft-masking, IUPAC codes, and RNA `U`/`u`).
4. **Falsifiability:** Pair every mathematical law with an automated falsification test ([`LawsFalsificationSuite.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/LawsFalsificationSuite.java)) and a real-genome ablation study ([`ablations_real_genomes.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/fusion/results/ablations_real_genomes.json)), and document every rejected hypothesis.

### 1.4 Complete Record of Rejected Hypotheses & Negative Levers Across Generations 1–3

| Generation | Tested Lever / Hypothesis | Observed Empirical Result | Root Cause & Disposition |
| :--- | :--- | :--- | :--- |
| **Gen 1 (`v1`)** | 1-bit identity stream + separate gap stream (`docs/ABLATION.md`) | `+252,188 B` (`+6.53%`) regression on 24g *E. coli* | Redundant with sparse mismatch mask (`0.088 b/bp`) + gap coding |
| **Gen 1 (`v1`)** | Match-length conditioning of post-gap literal | `+519 B` regression (`3,864,983 B` vs `3,864,464 B`) | Context dilution across sparse gap-length buckets |
| **Gen 1 (`v1`)** | Order-2 arithmetic coding of mismatch mask bytes | `+32,740 B` (`+0.85%`) regression on 24g *E. coli* | $256^2 = 65,536$ state space too sparse per `16 MiB` block |
| **Gen 1 (`v1`)** | Greedy novelty-selected anchor banks (`docs/ABLATION.md`) | `+10,948 B` (`16M`) to `+181,690 B` (`8M`) worse | Fragmented long syntenic runs; chronological stride anchors won |
| **Gen 1 (`v1`)** | Direct BWT on `100 Mbp` human DNA (`docs/SPECIALIZED.md`) | `2.189 bpb` (`27.36 MB`), **worse than 2-bit packing** | High-entropy repeat variation scrambles BWT contexts without CRM |
| **Gen 2 (`RC-MDL`)** | Fixed small anchor subsets ($N=4..20$) in Mode C (`RESULTS.md`) | `6,621,279 B` ($N=4$) $\to$ `4,115,368 B` ($N=20$) vs `3,972,975 B` | Limiting donor pool misses clade-specific accessory plasmids/phages |
| **Gen 2 (`RC-MDL`)** | Naive file-order or cluster-first Mode C ordering | `4,056,419 B` (`+2.10%`) / `4,081,636 B` (`+2.73%`) | Suboptimal tree jumps across distant phylogroups vs Containment/Prim |
| **Gen 3 (`HGS2`)** | Unclamped two-gear approximate matching (`clonalGear` w/o horizon) | `56,813 B` (`29,623` SNPs) on *E. coli* `00.seq` LOO | Bulldozed across mosaic recombination breakpoints; fixed by `findClonalReentryHorizon` (`40,735 B`, `5,706` SNPs) |
| **Gen 3 (`HGS2`)** | Micro-indel splicing on dense 23-donor *E. coli* LOO panel | `40,735 B` (ON) vs `40,607 B` (OFF) (`+128 B`) | Alternative donor bridges indel directly; retained because it saves `1,131 B` (`-15.56%`) on distant/indel-rich targets |

---

## 2. Formal Problem Statement, Six-Channel Bijective Factorization & Container Specifications

### 2.1 Definition 1 (Six-Channel Exact Container Bijection $\Phi$)
Implemented in [`FastaStreamFactorizer.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/FastaStreamFactorizer.java), let $\Sigma_{\text{raw}} = \{0, \dots, 255\}$ denote raw input bytes and $\Sigma_4 = \{0, 1, 2, 3\}$ denote the canonical 2-bit nucleotide alphabet under the complement-symmetric encoding:
$$\phi(\texttt{A}) = 0\;(00_2), \qquad \phi(\texttt{C}) = 1\;(01_2), \qquad \phi(\texttt{T/U}) = 2\;(10_2), \qquad \phi(\texttt{G}) = 3\;(11_2).$$
The forward bijection $\Phi: \Sigma_{\text{raw}}^* \to \mathcal{M} \times \Sigma_4^n$ decomposes any input byte stream $X \in \Sigma_{\text{raw}}^L$ into a canonical 2-bit DNA sequence $S \in \Sigma_4^n$ and a formatting metadata tuple $\mathcal{M} = (\mathcal{H}, \mathcal{W}, \mathcal{C}, \mathcal{E}, \tau)$:
1. **Header Channel $\mathcal{H}$:** Ordered sequence of FASTA header strings $\{h_1, \dots, h_m\}$ and their exact 0-based nucleotide start offsets $\{p_1, \dots, p_m\}$.
2. **Line-Width Run-Length Channel $\mathcal{W}$:** Run-length encoded line widths $(w_j, r_j)$ and line-ending type ($\text{LF}=0$, $\text{CRLF}=1$).
3. **Case-Mask Toggle Channel $\mathcal{C}$:** Delta-coded run lengths of uppercase/lowercase character transitions.
4. **Non-ACGT IUPAC Exception Channel $\mathcal{E}$:** Sparse delta-coded offsets $\{e_1, \dots, e_u\}$ and raw ASCII bytes for any non-$\{\texttt{A,C,G,T/U}\}$ base (substituted by $0$ [`A`] in $S$ so biological codon phase is preserved without coordinate shift).
5. **Terminal Newline & Alphabet Flags $\tau$:** Bit-flags recording trailing newline presence, $\text{CRLF}$ mode (`flags & 0x02`), and RNA $\texttt{T}\leftrightarrow\texttt{U}$ translation (`flags & 0x04`).

**Proposition 1 (Strict Lossless Identity):** For every raw input $X \in \Sigma_{\text{raw}}^L$, $\Phi^{-1}(\Phi(X)) = X$ identically. When $X$ is a headerless canonical `.seq` file, `flags & 0x01` is set and $\mathcal{M}$ occupies `0 bytes`.

### 2.2 Binary Container Specifications: `.hgs` v1.0 vs. `.hgs2` v2.0 ([`docs/HGS-FORMAT.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/HGS-FORMAT.md))

| Container Field / Feature | `.hgs` v1.0 Compact (`--pangenome` / `--single`) | `.hgs` v1.0 Multi-Block Indexed (`--fast` / `--lean`) | `.hgs2` v2.0 Unified (`Mode S` / `Mode R` / `Mode C`) |
| :--- | :--- | :--- | :--- |
| **Magic & Version** | `0x48475301` (`"HGS\x01"`, `4 B`) | `0x48475301` (`"HGS\x01"`, `4 B`) | `0x48475332` (`"HGS2"`, `4 B`) + `ver=2` (`1 B`) |
| **Mode / Flags Header** | `CompactFlag = 1` (`1 B`) + `origLen` (`4 B`) | `CompactFlag = 0` (`1 B`) + `BlockCount` (`4 B`) | `mode` (`0=S, 1=R, 2=C`, `1 B`) + `flags` (`1 B`) |
| **Length & Integrity Fields** | 32-bit `origLen` + `payloadLen` (`8 B`) | Per-block `origLen` (`4 B`) + `CRC-32` (`4 B`) | 64-bit `rawLen` (`8 B`) + `dnaLen` (`4 B`) + `donorCount` (`4 B`) |
| **Mode R Reference Pin** | N/A (self-contained only) | N/A (self-contained only) | **Mandatory 32-byte SHA-256 donor-panel digest** |
| **FASTA Formatting Metadata** | Coupled inside stream | Coupled inside stream | **6-channel decoupled $\mathcal{M} = (\mathcal{H},\mathcal{W},\mathcal{C},\mathcal{E},\tau)$** (`0 B` on `.seq`) |
| **Indexed Random Access** | Block-level (`8 MiB` sub-blocks) | **`HGF1` (`0x48474631`) footer** (`24 B/block` + `16 B` trailer) | Parallel `4 Mbp` Stage-2 `FastCRM` sub-blocks |
| **Worst-Case Safety Bypass** | Block mode `0` (`RAW`) | Block mode `0` (`RAW`) | Sub-block bit-31 raw 2-bit bypass (`len \| 0x80000000`, `2.0012 bpb`) |

---

## 3. Mathematical Formulation & Proofs of the Five Fundamental Laws

All five laws are formally specified in [`docs/COMPRESSION-LAW.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/COMPRESSION-LAW.md) and verified by [`LawsFalsificationSuite.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/LawsFalsificationSuite.java).

### 3.1 Law I: Klein Four-Group $V_4 = (\Sigma_4, \oplus)$ Strand-Invariant Kimura–CpG Substitution Law

Under $\phi(\texttt{A})=0, \phi(\texttt{C})=1, \phi(\texttt{T})=2, \phi(\texttt{G})=3$, Watson–Crick complementation is the bitwise XOR involution $\bar{b} = b \oplus 2$ ($0 \leftrightarrow 2$, $1 \leftrightarrow 3$). Equipping $\Sigma_4$ with bitwise XOR $\oplus$ forms the abelian Klein four-group $V_4 \cong \mathbb{Z}_2 \times \mathbb{Z}_2$.

#### Theorem 1 (Reverse-Complement Equivariance & Kimura–CpG Entropy Reduction)
Let $r \in \Sigma_4$ be a donor reference base and $t \in \Sigma_4$ ($t \neq r$) be a mutated target base, with non-zero group differential $\delta(r, t) = r \oplus t \in \{1, 2, 3\}$. Then:
1. **Exact Strand Invariance:** For all $r, t \in \Sigma_4$,
   $$\delta(\bar{r}, \bar{t}) = (r \oplus 2) \oplus (t \oplus 2) = r \oplus t \oplus (2 \oplus 2) = r \oplus t = \delta(r, t).$$
   Thus forward-strand and reverse-complement matches emit the **exact same** differential symbol without strand-conditional alphabet permutation.
2. **Algebraic Partition of Kimura's Two-Parameter Model:** The three non-zero cosets of $V_4 \setminus \{0\}$ partition all $4 \times 3 = 12$ directional nucleotide substitutions into three strand-invariant biochemical classes:
   - $r \oplus t = 3$ ($\delta_{\text{K}} = 1$ in [`HyperGraSSFusionEngine.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/HyperGraSSFusionEngine.java)): **Biochemical Transitions** ($\texttt{A}(0)\leftrightarrow\texttt{G}(3)$ and $\texttt{C}(1)\leftrightarrow\texttt{T}(2)$);
   - $r \oplus t = 2$ ($\delta_{\text{K}} = 2$): **Complement Transversions** ($\texttt{A}(0)\leftrightarrow\texttt{T}(2)$ and $\texttt{C}(1)\leftrightarrow\texttt{G}(3)$);
   - $r \oplus t = 1$ ($\delta_{\text{K}} = 3$): **Cross Transversions** ($\texttt{A}(0)\leftrightarrow\texttt{C}(1)$ and $\texttt{G}(3)\leftrightarrow\texttt{T}(2)$).
3. **Strict Entropy Reduction:** Conditioning on the donor's flanking CpG dinucleotide context $C_{\text{CpG}} \in \{0,1,2,3\}$ and preceding differential $\Delta_{\text{prev}}$ yields:
   $$H(\Delta_{\text{K}} \mid C_{\text{CpG}}, \Delta_{\text{prev}}) \le H(\Delta_{\text{K}}) = -\alpha\log_2\alpha - \beta_1\log_2\beta_1 - \beta_2\log_2\beta_2 < \log_2 3 \approx 1.5850\text{ bits} < 2.0000\text{ bits}.$$
   Empirically, Law I encodes substitutions in **$1.10\text{--}1.32\text{ bits/SNP}$**, saving **`6,028 B` (`3.06%`)** on distant single-donor *E. coli* (`00.seq | 20.seq`) and **`600 B` (`1.47%`)** on 23-donor *E. coli* LOO.

### 3.2 Law II: Bijective Codon-Wobble Gap Decomposition & Six-State Causal ORF-HMM

#### Theorem 2 (Zero Non-Coding Penalty & Positive Coding Mutual Information)
Define the bijective inter-SNP gap mapping $\psi: \mathbb{N}_0 \to \mathbb{N}_0 \times \{0, 1, 2\}$ by:
$$\psi(\text{gap}) = (q, \phi) = \left(\left\lfloor \frac{\text{gap}}{3} \right\rfloor,\; (\text{gap} + 1) \bmod 3\right), \qquad \psi^{-1}(q, \phi) = 3q + ((\phi + 2) \bmod 3).$$
By Euclidean division and the chain rule of Shannon entropy, $H(\text{gap}) = H(Q, \Phi) = H(Q) + H(\Phi \mid Q)$.
- **In non-coding DNA:** $\Phi$ is uniform on $\{0, 1, 2\}$ ($\log_2 3\text{ bits}$) and independent of $Q$, so decomposing $\text{gap}$ into $(q, \phi)$ incurs **zero asymptotic bit overhead** relative to direct integer coding of $\text{gap}$.
- **In protein-coding CDS ($\approx 85\text{--}90\%$ of bacterial chromosomes):** Purifying selection ($\omega = d_N/d_S \ll 1$) concentrates mutations at synonymous third-base wobble positions, creating a strong resonance at $\text{gap} \equiv 2 \pmod 3 \iff \phi = 0$. Consequently, $\Pr(\Phi = 0 \mid \text{CDS}) \gg 1/3$, strictly reducing per-gap code length by $I(\Phi; \text{CDS}) = \log_2 3 - H(\Phi \mid \text{CDS}) > 0\text{ bits/SNP}$ (saving **`1,854 B`** on *E. coli* `00.seq | 20.seq`).

#### Six-Frame Causal Open-Reading-Frame HMM (`6F-ORF-HMM` in [`FastNoveltyCoder.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/FastNoveltyCoder.java))
For unaligned novelty regions, `FastNoveltyCoder` tracks causal reading-frame survival across all 6 biological frames (3 forward phases $c = i \bmod 3 \in \{0,1,2\}$ on both $+$ and $-$ strands) in $O(1)$ time per base, testing 6-bit triplet history $t_3 \in \{0,\dots,63\}$ against forward stop codons $\{16\;(\texttt{TAA}),\; 18\;(\texttt{TAG}),\; 24\;(\texttt{TGA})\}$ and reverse-complement stop codons $\{36\;(\texttt{TTA}),\; 48\;(\texttt{CTA}),\; 50\;(\texttt{TCA})\}$, and classifying each base into $\text{orfActive} \in \{0, 1, 2\}$ to condition the phase-locked `codonOrder5Table` ($3 \times 2^{14}$ entries), four 24-byte Markov tables (`order10`, `order8`, `order6`, `order4`), and the `repeatTable` inside an 8-regime logistic mixer + 64-bin Secondary Symbol Estimator (`SSE`).

### 3.3 Law III: Conditional Monotonicity, Horizon-Bounded Two-Gear Coalescent Haplotype Copying & Affine Micro-Indel Splicing

#### Theorem 3 (Monotonic Conditional Description Length & Coalescent Scaling)
Under prefix-free coding of donor switch tokens, expanding a reference panel $D_{1:N} \to D_{1:N+1}$ satisfies $K(X \mid D_{1:N+1}) \le K(X \mid D_{1:N}) + O(1)$, with expected description length scaling inversely with panel size along the Li–Stephens terminal coalescent branch length $\mathbb{E}[\tau_N] \propto 1/N$.

#### Horizon-Bounded Two-Gear State Machine (`clonalGear` in [`HyperGraSSFusionEngine.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/HyperGraSSFusionEngine.java))
- **Gear 1 (Clonal Gear, `clonalGear = true`):** Credits `+0.008` bits/bp and penalizes each SNP by `-11.2` bits ($\approx 1$ SNP per $1,400\text{ bp}$ steady-state tolerance), enforcing a strict $\ge 128\text{ bp}$ exact 64-bit word-XOR tail guard (`countExactMatch64`) and halting whenever a 32-bp window exceeds $2$ SNPs.
- **Gear 2 (Diverged Syntenic Gear, `clonalGear = false`):** Engaged upon a local SNP burst ($>2$ SNPs in $32\text{ bp}$) **if and only if** the $O(1)$ syncmer lookahead probe `findClonalReentryHorizon` finds no alternative $\ge 192\text{ bp}$ exact clonal match in the donor panel within $4 \times 96 = 384\text{ bp}$. Gear 2 credits `+0.135` bits/bp and penalizes SNPs by `-10.2` bits ($\ge 18\text{ bp}$ tail guard), and immediately yields back to Gear 1 as soon as `findClonalReentryHorizon` detects a clonal donor reentry.
- **Affine Syntenic Micro-Indel Splicing (`attemptMicroIndelSplice`):** Searches $(\Delta_{\text{ref}}, \Delta_{\text{tgt}}) \in [1, 256] \times [1, 64]$ at match termination for a downstream continuation sharing $\ge 24\text{ bp}$ exact identity, replacing a $\sim 65\text{-bit}$ full token restart with a $\sim 12\text{-bit}$ affine splice descriptor.

### 3.4 Law IV: Universal Five-Regime Completeness Law

#### Theorem 4 (Strict Worst-Case Redundancy Bound)
Across any input sequence $S \in \Sigma_4^n$, the 5-tier grammar hierarchy—(1) Panel Syntenic Homology, (2) Bidirectional Self-Homology, (3) Tandem Periodicity Attractor ($1 \le P \le 64$), (4) `FastCRM` Context-Repeat Mixer, and (5) Sub-Block Raw 2-Bit Safety Bypass (`len | 0x80000000`)—guarantees:
$$L_{\text{HGS2}}(S) \le 2n + O(\log n)\text{ bits} \qquad \left(\le 2.0012\text{ bpb on } 500,000\text{ bp cryptographic white noise}\right).$$

### 3.5 Law V: Generalized Amdahl-Memory Streaming Scaling Law

#### Theorem 5 (Sub-Linear CSR Index Footprint & Two-Stage Memory Decoupling)
Using canonical open syncmers ($k=16, s_{\text{mer}}=10, w=6, t=3$) with subsampling shift $s \in \{0, 1, 2\}$, the expected number of indexed syncmers over $N_{\text{panel}}$ donor bases is $M = \frac{N_{\text{panel}}}{6 \cdot 2^s}$. Constructed via a two-pass zero-materialization Compressed Sparse Row (`CSR`) index (`int[1<<20]` heads = $4\text{ MiB}$, `int[M]` positions = $4M\text{ B}$, `short[M]` tags = $2M\text{ B}$) alongside the $0.25 N_{\text{panel}}$ packed 2-bit DNA bank, the exact Stage-1 core memory footprint is:
$$\mathcal{M}_{\text{Stage-1}}(N_{\text{panel}}, s) = N_{\text{panel}}\left(0.25 + 2^{-s}\right) + 4\text{ MiB},$$
yielding **`1.25 B/bp` raw (`1.58 B/bp` with JVM metadata) at $s=0$**, **`0.75 B/bp` raw (`0.93 B/bp` total) at $s=1$**, and **`0.50 B/bp` raw (`0.59 B/bp` total) at $s=2$**. Releasing the CSR arrays immediately after Stage-1 token parsing ensures Stage-2 parallel `FastCRM` coding ($6.8\text{ MB/thread}$ L3-cache-resident tables) executes comfortably inside `-Xmx256m`, `-Xmx224m`, and `-Xmx192m` pure-JVM heaps.

---

## 4. Generation 1 (`HyperGraSS v1`) Complete Experimental Record

Before developing `RC-MDL` and `HGS2`, we conducted an exhaustive experimental campaign with **Generation 1 (`HyperGraSS v1`, [`Main.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/Main.java))** across held-out synthetic and real datasets ([`docs/REFERENCE.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/REFERENCE.md)), anchor/block parameter sweeps ([`docs/SWEEP.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/SWEEP.md)), component ablations ([`docs/ABLATION.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/ABLATION.md)), single-genome `QCRM` and `100 Mbp` human chromosomes ([`docs/SPECIALIZED.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/SPECIALIZED.md)), cross-taxon bacterial collections ([`docs/BROAD.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/BROAD.md)), `256 MiB` native cgroup sweeps ([`docs/MEMORY-SWEEP.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/MEMORY-SWEEP.md)), `2.077 GiB` 64-bit streaming ([`docs/SCALING.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/SCALING.md)), the 8-channel `100:1` entropy floor breakdown ([`docs/ENTROPY-FLOOR.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/ENTROPY-FLOOR.md)), and 13-slice indexed random-access queries ([`docs/END-TO-END.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/END-TO-END.md)).

![Figure 7: Generation 1 (HyperGraSS v1) and Generation 2 (RC-MDL) Foundational Studies](%%FIG_PREFIX%%/fig7_gen1_gen2_foundations.png)

### 4.1 Frozen Held-Out Reference Corpus Across Seven Datasets ([`docs/REFERENCE.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/REFERENCE.md))

| Dataset | Input (`B`) | `HyperGraSS v1` (`B`, `bpb`) | `--optimize` (`B`, `bpb`) | `--lean` (`B`, `bpb`) | `gzip -9` (`B`, `bpb`) | `xz -9` (`B`, `bpb`) | `jarvis3-l7` (`B`, `bpb`) | `jarvis3-l26` (`B`, `bpb`) | `geco3-l6` (`B`, `bpb`) | `geco3-l9` (`B`, `bpb`) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `human-mito` | `16,569` | **`3,959`** (`1.912`) | **`3,959`** (`1.912`) | `3,975` (`1.919`) | `5,180` (`2.501`) | `4,800` (`2.318`) | `12,183` (`5.882`) | `4,877` (`2.355`) | `4,884` (`2.358`) | `4,884` (`2.358`) |
| `fly-mito` | `19,517` | **`3,567`** (`1.462`) | **`3,567`** (`1.462`) | `3,673` (`1.506`) | `4,947` (`2.028`) | `4,680` (`1.918`) | `4,394` (`1.801`) | `3,739` (`1.533`) | `3,760` (`1.541`) | `3,759` (`1.541`) |
| `yeast-mito` | `85,779` | **`15,338`** (`1.430`) | **`15,338`** (`1.430`) | `15,733` (`1.467`) | `21,390` (`1.995`) | `18,984` (`1.771`) | `16,116` (`1.503`) | `15,517` (`1.447`) | `15,433` (`1.439`) | `15,425` (`1.439`) |
| `plant-plastid` | `154,478` | **`30,419`** (`1.575`) | **`30,419`** (`1.575`) | `36,661` (`1.899`) | `39,192` (`2.030`) | `34,212` (`1.772`) | `30,623` (`1.586`) | `36,193` (`1.874`) | `30,847` (`1.597`) | `30,839` (`1.597`) |
| `yeast-chr1` | `230,218` | **`51,864`** (`1.802`) | **`51,864`** (`1.802`) | `54,521` (`1.895`) | `66,474` (`2.310`) | `56,884` (`1.977`) | `53,372` (`1.855`) | `52,331` (`1.818`) | `52,557` (`1.826`) | `52,541` (`1.826`) |
| `archaeal-ref` | `1,664,970` | **`370,196`** (`1.779`) | **`370,196`** (`1.779`) | `376,881` (`1.811`) | `484,565` (`2.328`) | `410,584` (`1.973`) | `414,492` (`1.992`) | *infeasible* | `413,569` (`1.987`) | *timeout* |
| `ecoli-k12` | `4,641,652` | **`1,087,650`** (`1.875`) | **`1,087,650`** (`1.875`) | `1,118,253` (`1.927`) | `1,341,697` (`2.312`) | `1,149,544` (`1.981`) | `1,100,664` (`1.897`) | `1,103,759` (`1.902`) | `1,094,264` (`1.886`) | *timeout* |

### 4.2 $3\times 3$ Anchor-Bank $\times$ Block-Size Grid & 10 MB Component Ablation ([`docs/ABLATION.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/ABLATION.md), [`docs/SPECIALIZED.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/SPECIALIZED.md))

| Anchor Cap | Block Cap | Blocks | Anchor Bases | Archive (`B`) | Bits/Base | $\Delta$ vs Default | $\Delta$ vs `xz -9` | In-Proc Time | Peak RSS |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `4 MiB` (Default) | `2 MiB` | `58` | `4,194,304` | `9,215,764` | `0.6138` | `0.00%` | `+9.78%` | `54.02 s` | `982.4 MB` |
| `4 MiB` | `4 MiB` | `29` | `4,194,304` | `8,370,127` | `0.5575` | `-9.18%` | `-0.30%` | `64.08 s` | `996.4 MB` |
| `4 MiB` | `8 MiB` | `15` | `4,194,304` | `7,264,536` | `0.4839` | `-21.17%` | `-13.46%` | `77.72 s` | `1,071.5 MB` |
| `8 MiB` | `2 MiB` | `58` | `8,388,608` | `6,026,650` | `0.4014` | `-34.60%` | `-28.21%` | `84.56 s` | `993.8 MB` |
| `8 MiB` | `4 MiB` | `29` | `8,388,608` | `5,863,560` | `0.3906` | `-36.37%` | `-30.15%` | `93.80 s` | `1,040.2 MB` |
| `8 MiB` | `8 MiB` | `15` | `8,388,608` | `5,551,039` | `0.3697` | `-39.77%` | `-33.88%` | `108.13 s` | `1,288.2 MB` |
| `16 MiB` | `2 MiB` | `58` | `9,405,717` | `5,180,683` | `0.3451` | `-43.78%` | `-38.29%` | `89.81 s` | `1,035.0 MB` |
| `16 MiB` | `4 MiB` | `29` | `9,405,717` | `5,060,967` | `0.3371` | `-45.08%` | `-39.71%` | `98.71 s` | `1,148.5 MB` |
| **`16 MiB` (`--pangenome`)** | **`8 MiB`** | **`15`** | **`9,405,717`** | **`4,885,851`** | **`0.3254`** | **`-46.98%`** | **`-41.80%`** | **`108.43 s`** | **`1,385.1 MB`** |
| `32 MiB` (Control) | `8 MiB` | `15` | `9,405,717` | `4,885,851` | `0.3254` | `-46.98%` | `-41.80%` | `110.22 s` | `1,297.2 MB` |
| **`16 MiB` + LZ77 Gap (`v1 Final`)** | **`8 MiB`** | **`15`** | **`9,405,717`** | **`3,810,409`** | **`0.2538`** | **`-58.65%`** | **`-54.61%`** | **`887.45 s`** | **`1,452.0 MB`** |

- **10 MB Component Ablation (`codec-ablation-10m-v2`, `9,859,515 bp`, [`docs/SPECIALIZED.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/SPECIALIZED.md)):** `QCRM` forced (no anchor bank) = `2,406,156 B` (`1.952 bpb`); **`Auto / AC-HMC` (Default) = `1,450,787 B` (`1.177 bpb`)**; `BC-LZG` forced = `1,454,172 B` (`1.180 bpb`); Anchor CM pages (parallel) = `1,463,055 B`; Aux without LZ77 = `1,450,787 B`; Aux without structured planes = `1,451,603 B`; Aux basic Deflate = `1,459,406 B`; Original `GraSS -e2` = `2,341,794 B` (`1.900 bpb`).

### 4.3 Single-Genome `QCRM` & `100 Mbp` Human `hg38` Chromosome Slices ([`docs/SPECIALIZED.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/SPECIALIZED.md))

| Dataset | Bases (`bp`) | `Orig. GraSS -e2` (`B`, `bpb`) | `v1 Default` (`B`, `bpb`) | `v1 --single` (`B`, `bpb`) | `v1 --single -t 2` (`B`, `bpb`) | $\Delta$ (`--single` vs `GraSS`) | Peak RSS (`--single` vs `GraSS`) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| *E. coli* K-12 MG1655 | `4,641,652` | `1,141,608` (`1.968`) | `1,086,715` (`1.873`) | **`1,084,172` (`1.869`)** | `1,090,807` (`1.880`) | **`-5.03%` (`-57,436 B`)** | **`88 MiB`** vs `278 MiB` |
| *B. subtilis* 168 | `4,215,606` | `1,034,078` (`1.962`) | `995,744` (`1.890`) | **`994,834` (`1.888`)** | `998,693` (`1.895`) | **`-3.80%` (`-39,244 B`)** | **`87 MiB`** vs `248 MiB` |
| *S. cerevisiae* (3 chr) | `1,360,022` | `333,798` (`1.963`) | `315,396` (`1.855`) | **`315,289` (`1.855`)** | `317,353` (`1.867`) | **`-5.55%` (`-18,509 B`)** | **`80 MiB`** vs `178 MiB` |
| *P. aeruginosa* NY4383 | `7,316,070` | `1,690,428` (`1.848`) | `1,533,558` (`1.677`) | **`1,529,536` (`1.673`)** | `1,548,347` (`1.693`) | **`-9.52%` (`-160,892 B`)** | **`95 MiB`** vs `296 MiB` |
| Human `hg38 chr1` (`100 Mbp`) | `100,000,000` | `22,264,066` (`1.781`) | `21,056,843` (`1.684`) | **`20,626,877` (`1.650`)** | `20,635,964` (`1.651`, `70.1 s`) | **`-7.35%` (`-1,637,189 B`)** | **`396/396 MiB`** vs `893/1,553 MiB` |
| Human `hg38 chr14` (`100 Mbp`) | `100,000,000` | `18,538,762` (`1.483`) | — | **`17,084,921` (`1.367`)** | — | **`-7.84%` (`-1,453,841 B`)** | **`380/389 MiB`** vs `960/1,524 MiB` |

### 4.4 Six Cross-Taxon Natural Collections & $3\times 3$ Mash Strata ([`docs/BROAD.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/docs/BROAD.md))

| Category / Species | Collection / Stratum | Genomes | Input Bytes (`B`) | `HyperGraSS v1` (`B`, `bpb`) | `JARVIS3 -l 7` (`B`, `bpb`) | `xz -9` (`B`, `bpb`) | `Orig. GraSS` (`B`) | $\Delta$ vs Best Baseline |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| Natural Collection | *Campylobacter jejuni* (`campy`) | `72` | `122,419,352` | **`1,535,296` (`0.100`)** | `2,017,046` (`0.132`) | `4,211,348` (`0.275`) | `24,598,042` | **`+23.88%`** vs `JARVIS3` |
| Natural Collection | *Chlamydia trachomatis* (`chlamydia`) | `66` | `69,754,471` | **`333,507` (`0.038`)** | `472,224` (`0.054`) | `862,388` (`0.099`) | `11,163,682` | **`+29.38%`** vs `JARVIS3` |
| Natural Collection | *Mycobacterium tuberculosis* (`mtb`) | `28` | `123,736,606` | **`1,139,535` (`0.074`)** | `1,418,155` (`0.092`) | `2,505,964` (`0.162`) | `26,389,531` | **`+19.65%`** vs `JARVIS3` |
| Natural Collection | *Streptococcus pneumoniae* (`pneumo`) | `57` | `121,007,978` | **`1,754,114` (`0.116`)** | `1,899,711` (`0.126`) | `5,960,712` (`0.394`) | `25,947,810` | **`+7.66%`** vs `JARVIS3` |
| Natural Collection | *Pseudomonas aeruginosa* (`pseudo`) | `19` | `123,241,439` | **`3,083,542` (`0.200`)** | `3,350,844` (`0.218`) | `10,847,460` (`0.704`) | `25,794,452` | **`+7.98%`** vs `JARVIS3` |
| Natural Collection | *Staphylococcus aureus* (`staph`) | `43` | `123,028,922` | **`1,807,399` (`0.118`)** | `1,895,866` (`0.123`) | `5,672,772` (`0.369`) | `25,470,552` | **`+4.67%`** vs `JARVIS3` |
| Diversity Stratum | `campy` — Low / Med / High | `24 x 3` | `~40.8 MB` ea. | **`444,619` / `586,251` / `961,772`** | `563,914` / `734,882` / `1,148,711` | `973K` / `1.45M` / `2.66M` | — | **`+21.2%` / `+20.2%` / `+16.3%`** |
| Diversity Stratum | `ecoli` — Low / Med / High | `20 x 3` | `~101.6 MB` ea. | **`1,829,036`** / `3,381,136` / **`4,373,930`** | `1,981,969` / **`3,375,860`** / `4,404,956` | `3.74M` / `8.07M` / `11.49M` | — | **`+7.7%`** / `-0.16%` / **`+0.7%`** |
| Diversity Stratum | `staph` — Low / Med / High | `14 x 3` | `~40.0 MB` ea. | **`345,519` / `583,848` / `1,037,833`** | `397,298` / `618,750` / `1,049,092` | `694K` / `1.48M` / `3.29M` | — | **`+13.0%` / `+5.6%` / `+1.1%`** |

### 4.5 `256 MiB` Native cgroup Gate, `2.077 GiB` 64-Bit Streaming, 8-Channel `100:1` Floor & 13-Slice Query Workload

| Study / Experiment | Profile / Channel | Bytes / Metric | Share / Rate | Compress / Query | Peak Memory | Verified Outcome |
| :--- | :--- | ---: | ---: | :--- | :--- | :--- |
| **`256 MiB` cgroup Gate** (`MEMORY-SWEEP.md`) | OpenJDK 25 JVM (`-Xmx136m..168m`) | OOM Killed | `0/3` pass | — | `>256.0 MiB` cgroup | Native JIT/malloc overhead exceeds `256 MiB` |
| | **GraalVM 25 Native (`-R:MaxHeapSize=180m`)** | **`3,864,464 B`** | **`0.257 bpb`** | `681.9 s` / `674.8 s` | **`240.1` / `225.5 MiB`** | **`3/3 PASS` (`251.1 MiB` cgroup peak)** |
| **`2.077 GiB` Streaming** (`SCALING.md`) | 64-bit `SequenceSource` (`132` blocks) | **`30,171,539 B`** | **`73.90x`** | `798.0 s` / `516.0 s` | **`235.0 MB` (`-Xmx160m`)** | **`PASS` (`CRC-32 0xEE111322` on `2.23 GB`)** |
| **8-Channel `100:1` Floor** (`ENTROPY-FLOOR.md`, `3,810,409 B` total) | Ch. 1: Anchor Bank (`9,405,717 bp` CM) | `2,257,234 B` | `59.24%` (`1.920 bpb`) | — | — | Alone $= 1.88\times$ the `1.20 MB` `100:1` budget |
| | Ch. 2: Variant Gaps (LZ77+range) | `882,215 B` | `23.15%` (`2.21 b/gap`) | — | — | `3,190,112` gaps across 15 blocks |
| | Ch. 3–6: Strand + Diff + Geom + Mask | `624,633 B` | `16.39%` | — | — | `239.8K + 160.2K + 143.8K + 80.8K B` |
| | Ch. 7–8: Literals + Donor/Mode Headers | `46,327 B` | `1.22%` | — | — | `30,647 B` literals (`128,808 bp`) + `15,680 B` hdr |
| **13-Slice Query** (`END-TO-END.md`, `13 x 2,000 bp`) | **HyperGraSS v1 `--lean`** | `8,952,202 B` | `0.596 bpb` | `31.33 s` / **`0.73 s`** | `804.6` / **`43.9 MiB`** | **`15.2x` faster query** vs `xz -9` (break-even `58`) |
| | **HyperGraSS v1 `--pangenome`** | **`4,909,412 B`** | **`0.327 bpb`** | `910.4 s` / **`1.35 s`** | `1,452` / **`74.1 MiB`** | **`8.2x` faster query** vs `xz -9` (break-even `859`) |
| | `xz -9` / `gzip -9` (Full Stream) | `8.39M` / `34.24M` | `0.559` / `2.280` | `48.0 s` / `11.12–34.14 s` | `674` / `32.9 MiB` | Requires full-archive decompression per query |

---

## 5. Generation 2 (`RC-MDL v2` Prototype) Complete Experimental Record

Generation 2 ([`prototype/rcmdl`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/prototype/rcmdl), [`prototype/rcmdl/RESULTS.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/prototype/rcmdl/RESULTS.md)) bridged the gap between Generation 1's block-anchored architecture and Generation 3's five-law engine by testing unaligned reference-conditioned MDL coding (**Mode R**) and tree-ordered pangenome coding (**Mode C**).

### 5.1 Complete 13-Configuration `AGC 3.2.5` Parameter Sweep ([`prototype/rcmdl/RESULTS.md`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/prototype/rcmdl/RESULTS.md))

| `AGC 3.2.5` Config ID | CLI Flags | *E. coli* `00.seq` (`4.86 Mbp`, `B`) | *E. coli* `00` `bpb` | *M. tuberculosis* 24g LOO (`B`) | *M. tb* `bpb` | *S. aureus* 25g LOO (`B`) | *S. aureus* `bpb` |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `default` | `-a` | `316,868` | `0.52201` | `443,669` | `0.03357` | `3,539,729` | `0.39860` |
| `fallback` | `-a -f 0.05` | `322,005` | `0.53048` | `517,362` | `0.03914` | `3,541,988` | `0.39885` |
| `min15` | `-a -l 15` | `329,505` | `0.54283` | `446,925` | `0.03381` | `3,560,914` | `0.40099` |
| `min15-fallback` | `-a -l 15 -f 0.05` | `334,729` | `0.55144` | `520,536` | `0.03938` | `3,596,087` | `0.40495` |
| `seg50k` | `-a -s 50000` | `307,143` | `0.50599` | `1,024,001` | `0.07747` | **`2,979,078`** | **`0.33547`** |
| `seg100k` | `-a -s 100000` | `316,369` | `0.52119` | **`292,170`** | **`0.02210`** | `3,737,764` | `0.42090` |
| `seg200k` | `-a -s 200000` | **`189,411`** | **`0.31204`** | `1,390,783` | `0.10522` | `4,188,783` | `0.47169` |
| `plain-seg200k` | `-s 200000` | **`189,411`** (24g: **`8,384,117`**) | **`0.31204`** (`0.5584`) | `1,390,783` | `0.10522` | `4,188,783` | `0.47169` |
| `seg200k-min15-fallback` | `-a -s 200000 -l 15 -f 0.05` | `189,900` | `0.31284` | `1,506,569` | `0.11398` | `4,235,561` | `0.47696` |
| `seg500k` | `-a -s 500000` | `224,332` | `0.36957` | `766,757` | `0.05801` | `4,458,726` | `0.50209` |
| `seg1m` | `-a -s 1000000` | `220,122` | `0.36263` | `384,994` | `0.02913` | `3,849,024` | `0.43343` |
| `seg1m-min15` | `-a -s 1000000 -l 15` | `229,178` | `0.37755` | `391,248` | `0.02960` | `3,956,000` | `0.44547` |
| `k25-seg200k` | `-a -k 25 -s 200000` | `394,600` | `0.65007` | `1,359,004` | `0.10281` | `3,744,704` | `0.42168` |
| **13-Config Per-Target `AGC` Oracle** | $\min_{c \in \mathcal{C}_{13}}$ per target | `189,411` (`8,384,117` 24g) | `0.31204` | **`96,580`** | **`0.00731`** | **`1,617,085`** | **`0.18210`** |
| **Generation 2 (`RC-MDL v1` / `v2`)** | Single default config | `73,194` / **`57,540`** (`1,636,367` 24g) | `0.1206` / **`0.0948`** | `122,998` / **`92,506`** | `0.00931` / **`0.00700`** | `267,412` / **`243,807`** | `0.03011` / **`0.02745`** |
| **Generation 3 (`HGS2` / `EIMDL`)** | **Five-Law Mode R** | **`40,735`** (**`1,177,778`** 24g) | **`0.06711`** (**`0.07845`**) | **`47,305`** | **`0.00358`** | **`174,639`** | **`0.01967`** |

### 5.2 Mode C Ordering Topologies, Fixed Anchor Sets & Mode R Donor Panel Scaling (`N=1..23`)

| Mode C Ordering / Anchor Topology (`ecoli24`) | Archive (`B`) | Bits/Base | $\Delta$ vs `RC-MDL` Best | Mode R Panel Size (*E. coli* `00.seq`) | `RC-MDL v1` (`B`, `bpb`) | `RC-MDL v2` (`B`, `bpb`) | `HGS2 / EIMDL` (`B`, `bpb`) |
| :--- | ---: | ---: | :--- | :--- | ---: | ---: | ---: |
| **Gen 3 `HGS2` (`Containment/Prim`, 1T)** | **`3,377,802`** | **`0.2250`** | **`-14.98%`** | $N=1$ donor | `249,491` (`0.4110`) | `273,654` (`0.4508`) | **`197,112` (`0.3247`)** |
| **Gen 1 `HyperGraSS v1` (`BWT+CRM`)** | **`3,810,409`** | **`0.2538`** | **`-4.09%`** | $N=2$ donors | `185,439` (`0.3055`) | `197,390` (`0.3252`) | — |
| **Gen 2 `RC-MDL v2` (`Containment/Prim`)** | **`3,972,975`** | **`0.2646`** | **Baseline** | $N=4$ donors | `146,758` (`0.2418`) | `185,481` (`0.3056`) | — |
| Gen 2 `RC-MDL` (`Cluster-aware Prim`) | `3,982,813` | `0.2653` | `+0.25%` | $N=8$ donors | `102,931` (`0.1696`) | `79,595` (`0.1311`) | — |
| Gen 2 `RC-MDL` (`Greedy Marginal`) | `3,990,249` | `0.26578` | `+0.43%` | $N=16$ donors | `87,254` (`0.1437`) | `65,671` (`0.1082`) | — |
| Gen 2 `RC-MDL` (`File Order / Cluster-First`) | `4,056,419` / `4,081,636` | `0.2702` / `0.2719` | `+2.10%` / `+2.73%` | $N=23$ donors (Full LOO) | `73,194` (`0.1206`) | `57,540` (`0.0948`) | **`40,735` (`0.0671`)** |
| Gen 2 Fixed Anchors ($N=20 / 16 / 12$) | `4,115,368` / `4,921,615` / `5,494,235` | `0.2741–0.3660` | `+3.58%` to `+38.29%` | **8-Donor Share of Gain** | **`83.1%` of $1\to 23$** | **`89.8%` of $1\to 23$** | — |
| Gen 2 Fixed Anchors ($N=8 / 4$) | `5,818,621` / `6,621,279` | `0.3876` / `0.4410` | `+46.46%` / `+66.66%` | **Total Gain ($1 \to 23$)** | **`-70.66%` (`3.41x`)** | **`-78.97%` (`4.76x`)** | **`-79.33%` (`4.84x`)** |
"""

MD_PART2 = r"""
---

## 6. Generation 3 (`HGS2` / `EIMDL`) Unified Multi-Regime Evaluation

### 6.1 Regime 1: Mode S (Reference-Free Single-Genome Compression)

#### Table 6A: Host Single-Genome Benchmarks ([`mode_s_and_synthetic.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/fusion/results/mode_s_and_synthetic.json))
> [!NOTE]
> `JARVIS3` was not run on these four host files (`docs/FUSION-BENCHMARKS.md` §1); see **Table 6B** directly below for the locked 12-genome independent Mode S evaluation against both `JARVIS3 -l 7` and `JARVIS3 -l 31`.

| Target File | Raw Bytes (`B`) | Canonical Bases (`bp`) | `HGS2` Archive (`B`) | Rate (`bpb`) | Ratio | Encode / Decode (`s`) | SHA-256 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| *E. coli* `00.fa` (FASTA) | `4,856,085` | `4,856,080` | **`1,123,893`** | `1.8515` | `4.32x` | `6.18` / `6.35` | `PASS` |
| *E. coli* `00.seq` (12T default) | `4,856,080` | `4,856,080` | **`1,123,865`** | `1.8515` | `4.32x` | `6.72` / `5.18` | `PASS` |
| *E. coli* `00.seq` (1T exact) | `4,856,080` | `4,856,080` | **`1,100,453`** | `1.8129` | `4.41x` | `14.82` / `14.51` | `PASS` |
| *E. coli* `00.seq` (`HyperGraSS v1 --single`) | `4,856,080` | `4,856,080` | `1,136,069` | `1.8716` | `4.27x` | `6.45` / `6.20` | `PASS` |
| *M. tuberculosis* `00.seq` | `4,417,521` | `4,417,521` | **`978,281`** | `1.7716` | `4.52x` | `7.54` / `5.90` | `PASS` |
| *S. aureus* `00.seq` | `2,845,025` | `2,845,025` | **`629,323`** | `1.7696` | `4.52x` | `4.91` / `3.61` | `PASS` |

#### Table 6B: Locked 12-Genome Independent *Salmonella* & *Klebsiella* Mode S and Mode R Evaluation ([`modes-subset-1`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/independent/runs/salmonella/modes-subset-1/report.json), [`jarvis-l31-subset-1`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/independent/runs/salmonella/jarvis-l31-subset-1/report.json))

| Species | Target Accession | Bases (`bp`) | `EIMDL Mode S` (`B`, `bpb`) | `JARVIS3 -l 31` (`B`, `bpb`) | `JARVIS3 -l 7` (`B`, `bpb`) | `Orig. GraSS -e2` (`B`, `bpb`) | `EIMDL Mode R` (`B`, `bpb`) |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
%%IND12_MD%%
| **Wall Time (Enc / Dec)** | **All 12 Genomes** | — | `332.5 s` / `277.9 s` | `417.5 s` / `452.1 s` | **`91.0 s` / `87.9 s`** | `127.0 s` / `100.6 s` | `325.1 s` / `120.4 s` |
| **Max cgroup Peak RSS** | **All 12 Genomes** | — | **`108 MiB` / `88 MiB`** | `1,158 MiB` / `1,159 MiB` | `266 MiB` / `266 MiB` | `1,680 MiB` / `618 MiB` | `526 MiB` / `98 MiB` |

![Figure 6: Independent 12-Genome Rate vs. Memory Pareto Frontier and Per-Target Archive Sizes](%%FIG_PREFIX%%/fig6_independent_modes_pareto.png)

---

### 6.2 Regime 2: Mode R (Reference-Conditioned Leave-One-Out & Leave-Cluster-Out Across 73 Genomes)

![Figure 2: Mode R Reference-Conditioned Leave-One-Out Across 73 Bacterial Genomes and Leave-Cluster-Out Controls](%%FIG_PREFIX%%/fig2_mode_r_loo_73genomes.png)

#### Table 7A: Complete 24-Genome *Escherichia coli* Mode R Leave-One-Out Benchmark ([`loo_ecoli.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/fusion/results/loo_ecoli.json))
- **Aggregate:** `120,108,223 bp` $\to$ **`1,177,778 B` (`0.07845 bpb`, `101.98x`)** vs. `RC-MDL v2` `1,636,367 B` (`-28.02%`) and `AGC 3.2.5` Best Single (`-s 200000`) `8,384,117 B` (**`7.12x` smaller, `24/24` per-target wins**). Total encode: `199.99 s`; total decode: `39.32 s`.

| Target | Raw (`B`) | `HGS2` (`B`) | `bpb` | Ratio | `RC-MDL v2` (`B`) | `AGC -s 200k` (`B`) | `AGC / HGS2` | Panel Match (`bp`) | Self Match (`bp`) | Novelty (`bp`) | SNPs | Splices | SHA-256 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
%%ECOLI_LOO_MD%%

#### Table 7B: Complete 24-Genome *Mycobacterium tuberculosis* Mode R Leave-One-Out Benchmark ([`loo_mtb.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/fusion/results/loo_mtb.json), [`agc-sweep-mtb-full/report.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/prototype/rcmdl/work/agc-sweep-mtb-full/report.json))
- **Aggregate:** `105,745,015 bp` $\to$ **`47,305 B` (`0.00358 bpb`, `2,235.39x`)** vs. `RC-MDL v2` `92,506 B` (`-48.86%`), `AGC 3.2.5` Best Single (`-s 100000`) `292,170 B` (**`6.18x` smaller, `24/24` wins**), and the **13-Configuration Per-Target `AGC` Oracle** `96,580 B` (**`2.04x` [`-51.02%`] smaller in aggregate, winning `19/24` targets**). Total encode: `157.57 s`; total decode: `23.84 s`.

| Target | Raw (`B`) | `HGS2` (`B`) | `bpb` | Ratio | `RC-MDL v2` (`B`) | `AGC 13-Cfg Oracle` (`B`, Flag) | Panel Match (`bp`) | Self Match (`bp`) | Novelty (`bp`) | SNPs | Splices | SHA-256 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
%%MTB_LOO_MD%%

#### Table 7C: Complete 25-Genome *Staphylococcus aureus* Mode R Leave-One-Out Benchmark ([`loo_saureus.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/fusion/results/loo_saureus.json), [`agc-verify-saureus-oracle/results.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/prototype/rcmdl/work/agc-verify-saureus-oracle/results.json))
- **Aggregate:** `71,043,114 bp` $\to$ **`174,639 B` (`0.01967 bpb`, `406.80x`)** vs. `RC-MDL v2` `243,807 B` (`-28.37%`), `AGC 3.2.5` Best Single (`-s 50000`) `2,979,078 B` (**`17.06x` smaller, `25/25` wins**; `20.27x` smaller than `AGC -a` `3,539,729 B`), and the **13-Configuration Per-Target `AGC` Oracle** `1,617,085 B` (**`9.26x` [`-89.20%`] smaller in aggregate, winning `24/25` targets**). Total encode: `189.48 s`; total decode: `25.41 s`.

| Target | Raw (`B`) | `HGS2` (`B`) | `bpb` | Ratio | `RC-MDL v2` (`B`) | `AGC 13-Cfg Oracle` (`B`, Flag) | Panel Match (`bp`) | Self Match (`bp`) | Novelty (`bp`) | SNPs | Splices | SHA-256 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
%%SAUREUS_LOO_MD%%

---

### 6.3 Regime 3: Mode C (Self-Contained Pangenome Collections & Independent Locked Trio)

![Figure 3: Mode C Self-Contained Pangenome Collection Benchmarks Across 5 Bacterial Species](%%FIG_PREFIX%%/fig3_mode_c_pangenomes.png)

- **Canonical 24-Genome *E. coli* Flat & Multi-FASTA Collection (`120,108,223 bp`):**
  - `EIMDL / HGS2 (--threads 1)`: **`3,377,802 B` (`0.2250 bpb`, `35.56x`)** (`96.93 s` encode, `55.44 s` decode)
  - `EIMDL / HGS2 (--threads 2)`: **`3,380,134 B` (`0.2251 bpb`, `35.53x`)** (`62.06 s` encode, `39.83 s` decode)
  - `EIMDL / HGS2 (--threads 4)`: **`3,388,767 B` (`0.2257 bpb`, `35.44x`)** (`51.94 s` encode, `26.02 s` decode)
  - `EIMDL / HGS2 (--threads 8 / Default 12T)`: **`3,393,302 B` (`0.2260 bpb`, `35.40x`)** (`51.74 s` / `44.04 s` encode, `23.27 s` / `18.36 s` decode)
  - `EIMDL / HGS2 (6-Channel Multi-FASTA)`: **`3,392,869 B` (`0.2260 bpb`, `35.91x`)** (`63.78 s` encode, `57.07 s` decode)
  - `EIMDL / HGS2 (-Xmx256m / -Xmx224m)`: **`3,423,050 B` (`0.2280 bpb`, `35.09x`)**
  - `EIMDL / HGS2 (-Xmx192m)`: **`3,455,169 B` (`0.2301 bpb`, `34.76x`)** (`57.64 s` encode, `27.19 s` decode)
  - `HyperGraSS v1 (BWT+CRM)`: `3,810,409 B` (`0.2538 bpb`, `31.52x`, `887.45 s` encode, `1,452 MB` RSS)
  - `JARVIS3 v3.7 (-l 7)`: `3,896,361 B` (`0.2595 bpb`, `30.83x`, `95.36 s` encode, `549 MB` RSS)
  - `RC-MDL v2 (Prim Tree)`: `3,972,975 B` (`0.2646 bpb`, `30.23x`, `164.59 s` encode, `1,812 MB` RSS)
  - `xz -9 (Flat / Multi-FASTA)`: `5,663,640 B` (`0.3772 bpb`) / `8,394,900 B` (`0.5592 bpb`)
  - `GeCo3 (-l 6 / -l 3)`: `8,409,358 B` (`0.5601 bpb`) / `14,389,043 B` (`0.9584 bpb`)
  - `AGC 3.2.5 (-a -t 1)`: `8,540,341 B` (`0.5688 bpb`)
  - `Original GraSS (-e2)`: `24,777,592 B` (`1.6503 bpb`)

- **Locked Independent *Enterobacteriaceae* Trio (`benchmarks/independent/runs/{salmonella,klebsiella}/trio-1/report.json`):**
  - ***Salmonella enterica* (25 genomes, `120,527,038 bp`):** `EIMDL Mode C` achieves **`2,390,106 B` (`0.1586 bpb`, `50.43x`)** vs. `JARVIS3 -l 7` `2,709,545 B` (`0.1798 bpb`, `-11.79%`) and `Original GraSS -e2` `19,991,646 B` (`1.3269 bpb`, `-88.04%`, `8.36x` smaller), with `53.36 s` decode (`1.61x` faster than `JARVIS3` at `125.1 MiB` vs. `294.4 MiB` cgroup RSS).
  - ***Klebsiella pneumoniae* (23 genomes, `122,695,431 bp`):** `EIMDL Mode C` achieves **`2,523,991 B` (`0.1646 bpb`, `48.61x`)** vs. `JARVIS3 -l 7` `2,948,926 B` (`0.1923 bpb`, `-14.41%`) and `Original GraSS -e2` `21,191,566 B` (`1.3817 bpb`, `-88.09%`, `8.40x` smaller), with `32.69 s` encode / `16.19 s` decode (`1.69x` faster encode and `3.21x` faster decode than `JARVIS3 -l 7` at `177.3 MiB` vs. `295.5 MiB` decode cgroup RSS).

---

### 6.4 Verification of Law V (Amdahl Multi-Core Scaling & `-Xmx192m` Bounded Memory) and Law IV (Synthetic & Adversarial Genomes)

![Figure 4: Law V Amdahl Multi-Core Scaling and Pure-JVM Bounded-Heap Sweep](%%FIG_PREFIX%%/fig4_amdahl_and_bounded_memory.png)

![Figure 5: Law-by-Law Real-Genome Ablations, Donor Panel Scaling, Entropy Floor Progression, and Synthetic/Adversarial Genomes](%%FIG_PREFIX%%/fig5_laws_ablations_and_synthetic.png)

#### Table 8: Law IV Evaluation Across Nine Non-Existing, Xenobiological, Synthetic & Adversarial Genomes ([`mode_s_and_synthetic.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/fusion/results/mode_s_and_synthetic.json))

| Corpus ID | Target File | Regime | Raw Bytes (`B`) | Bases (`bp`) | `HGS2` (`B`) | `bpb` (Bases) | `bpb` (Raw) | Ratio | `gzip -9` (`B`) | SHA-256 |
| :--- | :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| `01_euk_centromere` | `01_euk_centromere.fasta` | Mode S | `491,604` | `485,418` | **`42,193`** | `0.6954` | `0.6866` | `11.65x` | `83,797` | `PASS` |
| `02_viral_rna_acgu` | `02_viral_rna_acgu.fasta` | Mode C | `304,013` | `299,253` | **`12,122`** | `0.3241` | `0.3190` | `25.08x` | `42,923` | `PASS` |
| `03_chimeric_iupac_crlf` | `03_chimeric_iupac_crlf.fa` | Mode S | `263,817` | `244,065` | **`51,518`** | `1.6887` | `1.5622` | `5.12x` | `86,566` | `PASS` |
| `04_xeno_at88` | `04_xeno_at88.fasta` | Mode S | `506,294` | `500,000` | **`94,222`** | `1.5076` | `1.4888` | `5.37x` | `123,139` | `PASS` |
| `05_xeno_gc86` | `05_xeno_gc86.fasta` | Mode S | `506,291` | `500,000` | **`97,524`** | `1.5604` | `1.5410` | `5.19x` | `126,076` | `PASS` |
| `06_denovo_grammar` | `06_denovo_grammar.fasta` | Mode S | `506,287` | `500,000` | **`41,096`** | `0.6575` | `0.6494` | `12.32x` | `141,540` | `PASS` |
| `07_hypermutator_pangenome` | `07_hypermutator_pangenome.fasta` | Mode C | `1,513,506` | `1,494,677` | **`101,566`** | `0.5436` | `0.5369` | `14.90x` | `454,203` | `PASS` |
| `07_hypermutator_strain05_loo` | `strain_05.seq` (vs 5-strain panel) | Mode R | `249,436` | `249,436` | **`7,514`** | `0.2410` | `0.2410` | `33.20x` | — | `PASS` |
| `08_crypto_white_noise` | `08_crypto_white_noise.seq` | Mode S | `500,000` | `500,000` | **`125,075`** | **`2.0012`** | **`2.0012`** | `4.00x` | `144,017` | `PASS` |

---

### 6.5 Controlled Real-Genome Law-by-Law Ablations ([`ablations_real_genomes.json`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/benchmarks/fusion/results/ablations_real_genomes.json))

| Scenario | Ablation Configuration | Archive (`B`) | Rate (`bpb`) | $\Delta$ vs. Full `EIMDL` | Panel Match (`bp`) | Panel Tokens | SNPs | Splices | SHA-256 |
| :--- | :--- | ---: | ---: | :--- | ---: | ---: | ---: | ---: | :---: |
| **(A) *E. coli* `00.seq` \| `20.seq`** | **Full `EIMDL` (All Laws ON)** | **`197,112`** | **`0.32473`** | **Baseline** | `4,258,080` | `1,517` | `54,124` | `619` | `PASS` |
| (1 distant donor, `4,856,080 bp`) | `no_kimura_v4_law_I` (Flat SNP) | `203,140` | `0.33466` | **`+6,028 B` (`+3.06%`)** | `4,257,425` | `1,490` | `53,987` | `618` | `PASS` |
| | `no_codon_wobble_law_II` (Flat gap) | `198,966` | `0.32778` | **`+1,854 B` (`+0.94%`)** | `4,258,080` | `1,517` | `54,124` | `619` | `PASS` |
| | `no_indel_splice_law_III` (No splice) | `197,205` | `0.32488` | **`+93 B` (`+428` broken tokens)** | `4,259,145` | `1,945` | `54,157` | `0` | `PASS` |
| | `no_v4_no_wobble_no_splice` | `204,912` | `0.33758` | **`+7,800 B` (`+3.96%`)** | `4,258,398` | `1,918` | `54,006` | `0` | `PASS` |
| | `exact_matches_only_no_approx` | `239,083` | `0.39387` | **`+41,971 B` (`+21.29%`)** | `4,028,211` | `30,701` | `0` | `0` | `PASS` |
| **(B) *E. coli* `00.seq` 23-Donor LOO** | **Full `EIMDL` (All Laws ON)** | **`40,735`** | **`0.06711`** | **Baseline** | `4,766,688` | `5,531` | `5,706` | `99` | `PASS` |
| (23 donors, `4,856,080 bp`) | `no_kimura_v4_law_I` (Flat SNP) | `41,335` | `0.06810` | **`+600 B` (`+1.47%`)** | `4,766,538` | `5,524` | `5,680` | `99` | `PASS` |
| | `no_codon_wobble_law_II` (Flat gap) | `40,856` | `0.06731` | **`+121 B` (`+0.30%`)** | `4,766,688` | `5,531` | `5,706` | `99` | `PASS` |
| | `no_indel_splice_law_III` (No splice) | **`40,607`** | **`0.06690`** | **`-128 B` (dense 23-donor panel)** | `4,766,614` | `5,590` | `5,685` | `0` | `PASS` |
| | `no_v4_no_wobble_no_splice` | `41,330` | `0.06809` | **`+595 B` (`+1.46%`)** | `4,766,417` | `5,583` | `5,656` | `0` | `PASS` |
| | `exact_matches_only_no_approx` | `43,596` | `0.07182` | **`+2,861 B` (`+7.02%`)** | `4,747,677` | `7,791` | `0` | `0` | `PASS` |
| **(C) *M. tuberculosis* `00.seq` LOO** | **Full `EIMDL` (All Laws ON)** | **`1,237`** | **`0.00224`** | **Baseline** | `4,412,383` | `116` | `251` | `0` | `PASS` |
| (23 donors, `4,417,521 bp`) | `no_kimura_v4_law_I` (Flat SNP) | `1,258` | `0.00228` | **`+21 B` (`+1.70%`)** | `4,412,383` | `116` | `251` | `0` | `PASS` |
| | `no_codon_wobble_law_II` (Flat gap) | `1,235` | `0.00224` | `-2 B` (`251` SNPs total) | `4,412,383` | `116` | `251` | `0` | `PASS` |
| | `no_indel_splice_law_III` (No splice) | `1,236` | `0.00224` | `-1 B` (`0` splices) | `4,412,383` | `116` | `251` | `0` | `PASS` |
| | `no_v4_no_wobble_no_splice` | `1,255` | `0.00227` | **`+18 B` (`+1.46%`)** | `4,412,383` | `116` | `251` | `0` | `PASS` |
| | `exact_matches_only_no_approx` | `1,275` | `0.00231` | **`+38 B` (`+3.07%`)** | `4,412,761` | `325` | `0` | `0` | `PASS` |

---

### 6.6 Programmatic Falsification Suites (`171,933` Automated Assertions in [`LawsFalsificationSuite.java`](file:///C:/Users/mdsai/.gemini/antigravity/scratch/hypergrass/src/hypergrass/LawsFalsificationSuite.java))
- **Suite 1 (Law I $V_4$ Reverse-Complement Invariance & Kimura Gain):** Exhaustively verifies $(r \oplus 2) \oplus (t \oplus 2) = r \oplus t$ across all $16$ pairs and confirms `31,599 B` ($V_4$ ON) vs. `32,126 B` (OFF), saving **`527 B`** on synthetic Kimura mutated genomes.
- **Suite 2 (Law II Codon-Wobble Bijection & 3-Periodicity Gain):** Exhaustively verifies $\psi^{-1}(\psi(\text{gap})) = \text{gap}$ for all $\text{gap} \in [0, 100{,}000)$ (`100,000` assertions) and confirms `22,521 B` (Wobble ON) vs. `24,782 B` (OFF), saving **`2,261 B` (`-9.12%`)**.
- **Suite 3 (Law III Monotonicity & Affine Micro-Indel Splicing):** Confirms monotonic description length ($D_1$: `32,202 B` $\ge D_{1:2}$: `15,879 B` $\ge D_{1:4}$: `187 B`) and **`1,131 B` (`-15.56%`)** savings (`6,136 B` vs. `7,267 B`) from affine micro-indel splicing.
- **Suite 4 (Law IV Six-Channel Container & `2.0012 bpb` Safety Ceiling):** Verifies exact SHA-256 reconstruction across multi-record FASTA with `CRLF`, variable widths, soft-masking, IUPAC codes, and RNA `U/u`, plus `250,183 B` (`2.0015 bpb`) on `1,000,000 bp` white noise.
- **Suite 5 (Law V Sub-Linear CSR Index Scaling & Thread Determinism):** Verifies `1.58 B/bp` ($s=0$), `0.93 B/bp` ($s=1$), and `0.59 B/bp` ($s=2$), plus byte-identical archives across `1`, `2`, `4`, and `8` threads.
- **Suite 6 (`FastCRM` 6-Frame Causal ORF-HMM):** Confirms `34,475 B` (`1.3790 bpb`) on stop-free coding reading frames vs. `35,983 B` (`1.4393 bpb`) on stop-rich control (**`-1,508 B`, `-0.0603 bpb`**).
- **Suite 7 (Mode R Cryptographic Donor-Panel Integrity Pin):** Confirms immediate rejection (`IllegalArgumentException: Mode R donor panel SHA-256 mismatch`) upon a single-base donor mutation.

---

## 7. Honest Scientific Boundaries, Negative Results & Threats to Validity

1. **Dense-Panel vs. Distant-Donor Micro-Indel Splicing Trade-Off:** As shown in Section 6.5, affine micro-indel splicing (Law IIIb) saves `1,131 B` (`-15.56%`) on synthetic indel benchmarks and eliminates `428` broken tokens (`-93 B`) on distant single-donor *E. coli* (`00.seq | 20.seq`), but costs `+128 B` (`40,735 B` vs. `40,607 B`) on the dense 23-donor *E. coli* LOO panel where an alternative donor frequently bridges the indel directly without needing a splice.
2. **Near-Clonal Per-Target Differences vs. 13-Config `AGC` Oracle:** While `HGS2` wins **`73/73` targets** against any single deployable `AGC 3.2.5` configuration and beats the 13-configuration per-target `AGC` oracle by **`2.04x`** on *M. tuberculosis* (`47,305 B` vs. `96,580 B`) and **`9.26x`** on *S. aureus* (`174,639 B` vs. `1,617,085 B`), the per-target oracle is smaller by `13–236 B` on 5/24 near-clonal *M. tuberculosis* targets (`03, 08, 09, 16, 19.seq`) and by `64 B` on 1/25 *S. aureus* target (`04.seq`), where `HGS2`'s mandatory 32-byte SHA-256 donor-panel manifest pin represents a non-trivial fraction of the archive.
3. **Mode S CPU Time vs. Lightweight `JARVIS3 -l 7`:** In Mode S single-genome compression (Table 6B), `EIMDL Mode S` (`--threads 1`) is `8.79%` smaller than default `JARVIS3 -l 7` (`13,868,405 B` vs. `15,204,327 B`) with `2.46x` lower memory (`108 MiB` vs. `266 MiB`), and `2.32%` smaller and `1.26x` faster than tuned `JARVIS3 -l 31` (`14,197,393 B`, `1,158 MiB` RAM), winning `12/12` targets against both. However, against default `JARVIS3 -l 7`, `EIMDL Mode S` is `3.65x` slower to compress (`332.5 s` vs. `91.0 s`) and `3.16x` slower to restore (`277.9 s` vs. `87.9 s`).
4. **32-Bit Signed Integer CSR Addressing Horizon (`~2.14 Gbp`):** `DonorPanelIndex` stores concatenated donor coordinates in 32-bit signed `int[]` arrays, bounding a single unsharded reference panel to $2^{31}-1 \approx 2.14\text{ Gbp}$ (whereas Generation 1's 64-bit `SequenceSource` streaming path in `Main.java` handles $>2\text{ GiB}$ block-by-block streams directly).
"""


def write_all_markdown():
    tables = build_all_tables()
    body = (
        (MD_PART1 + MD_PART2)
        .replace("%%IND12_MD%%", tables["ind12_md"])
        .replace("%%ECOLI_LOO_MD%%", tables["ecoli_loo_md"])
        .replace("%%MTB_LOO_MD%%", tables["mtb_loo_md"])
        .replace("%%SAUREUS_LOO_MD%%", tables["saureus_loo_md"])
    )

    repo_md = body.replace("%%FIG_PREFIX%%", "../paper/figures")
    repo_path = os.path.join(DOCS_DIR, "PAPER-DRAFT.md")
    with open(repo_path, "w", encoding="utf-8") as f:
        f.write(repo_md)
    print(f"[OK] Wrote {repo_path} ({len(repo_md):,} bytes)")

    art_md = body.replace("%%FIG_PREFIX%%", ARTIFACT_DIR.replace("\\", "/"))
    art_temp = os.path.join(BASE_DIR, "paper", "artifact_preview.md")
    with open(art_temp, "w", encoding="utf-8") as f:
        f.write(art_md)
    print(f"[OK] Wrote {art_temp} ({len(art_md):,} bytes)")


if __name__ == "__main__":
    write_all_markdown()
