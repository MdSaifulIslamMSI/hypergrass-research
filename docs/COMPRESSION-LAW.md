# The Evolutionary-Information MDL (EIMDL) & Amdahl-Memory Laws

This document states the mathematical, biological, and architectural laws governing the **HyperGraSS + RC-MDL Unified Engine (`HGS2` / `EIMDL`)**. In accordance with [FOUNDATIONS.md](FOUNDATIONS.md), every law states an explicit **invariant**, **mechanism**, **predicted effect**, and **executable falsification test** (implemented in [`VerifyCompressionLaws.java`](../src/hypergrass/VerifyCompressionLaws.java)).

---

## 1. Unified Problem Statement Across All Three Contracts

Let $\Sigma_{\text{raw}}$ be the arbitrary byte alphabet of real-world sequence files (headers, line widths, `LF`/`CRLF`, lowercase masking, `N`-runs, RNA `U`, and IUPAC ambiguity codes), and let $\Phi: \Sigma_{\text{raw}}^* \leftrightarrow \mathcal{M} \times \Sigma_4^n$ be the invertible 6-channel bijection implemented by [`FastaStreamFactorizer`](../src/hypergrass/FastaStreamFactorizer.java), where $\Sigma_4 = \{0(\text{A}), 1(\text{C}), 2(\text{G}), 3(\text{T/U})\}$.

For any target sequence $X \in \Sigma_4^n$ conditioned on a (possibly empty) reference state $\mathcal{R} = (R_1, \dots, R_m)$:

$$L_{\text{total}}(\text{file} \mid \mathcal{R}) = L_{\text{meta}}(\mathcal{M}) + L_{\text{pin}}(\mathcal{R}) + L_{\text{EIMDL}}(X \mid \mathcal{R})$$

All three operational compression regimes are specializations of this single functional:
- **Mode S (Universal Single Genome, $\mathcal{R} = \emptyset$):** Self-contained single sequence or assembly (existing biological or non-existing/synthetic/de novo).
- **Mode R (Reference-Conditioned Delta, $\mathcal{R} \neq \emptyset$ external):** Target encoded relative to an external panel pinned by a cryptographic 32-byte Merkle SHA-256 manifest digest $L_{\text{pin}}(\mathcal{R})$.
- **Mode C (Self-Contained Pangenome Collection, $X = (G_1, \dots, G_N), \mathcal{R} = \emptyset$):** Joint archive over $N$ genomes factored along a bounded-depth phylogenetic reference DAG $\mathcal{T} = (V, E)$ ($\text{depth}(\mathcal{T}) \le D$) coupled to a unified pan-genome novelty superstream $\mathcal{L}$.

---

## 2. Law I — Klein Four-Group ($V_4$) Strand-Invariant Kimura-CpG Substitution Law

### Mathematical & Biological Invariant
Equip $\Sigma_4 = \{0(\text{A}), 1(\text{C}), 2(\text{G}), 3(\text{T})\}$ with bitwise XOR $\oplus$, forming the abelian Klein four-group $V_4 \cong \mathbb{Z}_2 \times \mathbb{Z}_2$:
- Watson-Crick complementation is translation by $3 = 11_2$: $\overline{x} = x \oplus 3$.
- At any homologous substitution locus where target base $y \neq x$ (with donor base $x$ known to the decoder), the identity element $\delta = 0$ has **exact probability zero** ($P(y = x \mid y \neq x) = 0$).
- The non-zero group elements $\delta = y \oplus x \in \{1, 2, 3\}$ partition all 12 possible nucleotide substitutions into **three strand-invariant evolutionary operators**:
  1. $\delta = 2$ ($10_2$): **Kimura Transition** ($\text{A} \leftrightarrow \text{G}$ purine-purine, $\text{C} \leftrightarrow \text{T}$ pyrimidine-pyrimidine).
  2. $\delta = 1$ ($01_2$): **Keto-Amino Transversion** ($\text{A} \leftrightarrow \text{C}$, $\text{G} \leftrightarrow \text{T}$).
  3. $\delta = 3$ ($11_2$): **Watson-Crick Complement Transversion** ($\text{A} \leftrightarrow \text{T}$, $\text{C} \leftrightarrow \text{G}$).

Because $(\overline{y}) \oplus (\overline{x}) = (y \oplus 3) \oplus (x \oplus 3) = y \oplus x = \delta$, every operator $\delta \in \{1, 2, 3\}$ commutes identically with reverse-complementation:

$$\mathcal{RC}(\delta(x)) = \delta(\mathcal{RC}(x)) \quad \forall x \in \Sigma_4, \; \delta \in V_4 \setminus \{0\}$$

### Entropy Bound
Let $p_{\text{ti}}(c) = P(\delta = 2 \mid c, y \neq x)$ denote the transition probability in local reference trinucleotide context $c = (x_{-1}, x_0, x_{+1})$ (capturing 5-methylcytosine $\text{CpG} \to \text{TpG}/\text{CpA}$ hypermutability) and divergence regime $d$, with transversion split $p_1(c) + p_3(c) = 1 - p_{\text{ti}}(c)$. The per-substitution allele code length satisfies:

$$H_{V_4}(Y \mid X, c, Y \neq X) = h_2(p_{\text{ti}}(c)) + (1 - p_{\text{ti}}(c)) \, h_2\!\left(\frac{p_1(c)}{1 - p_{\text{ti}}(c)}\right) \le \log_2 3 \approx 1.5850 \text{ bits}$$

strictly below unconditioned 4-symbol coding ($2.0000$ bits). Under typical Kimura transition bias ($\kappa \in [3, 6] \implies p_{\text{ti}} \in [0.60, 0.75]$), $H_{V_4} \in [1.10, 1.32]$ bits/substitution, saving $0.68\text{--}0.90$ bits per SNP.

### Falsification Test
1. Verify algebraic involution invariance $(y \oplus 3) \oplus (x \oplus 3) = y \oplus x$ across all $4 \times 3 = 12$ mutant pairs.
2. On synthetic Kimura sequences ($\kappa = 5.0$) and real *E. coli* homolog pairs, ablating $V_4$ conditioning back to 4-symbol allele coding must strictly increase substitution bitstream size; on uniform synthetic transversion+transition noise ($\kappa = 1$), $V_4$ coding must converge to $\le 1.60$ bits/SNP (never exceeding $\log_2 3 + \epsilon$).

---

## 3. Law II — Bijective Codon-Wobble Gap & 6-State Causal ORF-HMM Phase Law

### Mathematical & Biological Invariant
In protein-coding sequence (~85–90% of prokaryotic genomes and coding exons), the genetic code imposes two distinct 3-periodic invariants:
1. **Homologous Substitution Phase Law (Law IIa):** Purifying selection ($\omega = d_N / d_S \ll 1$) concentrates homologous SNPs at 3rd-codon synonymous ("wobble") positions, inducing strong 3-periodicity $\Delta p = \text{gap} + 1 \equiv 0 \pmod 3$ when $\text{gap} \ge 2$, whereas multinucleotide/repeat mutations concentrate at $\text{gap} = 0$ ($\Delta p = 1$). Simultaneously, horizontal gene transfer and recombination create piecewise-constant coalescent TMRCA branch lengths.
2. **De Novo 6-Frame Coding Phase & Strand Law (Law IIb):** Across real chromosomes, thousands of open reading frames (ORFs) sit on **both strands** ($+$ and $-$) separated by intergenic spacers whose lengths are *not* multiples of $3$. Consequently, global coordinate phase $t \bmod 3$ decorrelates from biological codon position $\{0, 1, 2\}$ across gene boundaries. However, inside every sense ($+$) ORF of offset $f \in \{0, 1, 2\}$, internal stop triplets $\{\text{TAA}, \text{TAG}, \text{TGA}\}$ are depleted to $\approx 0$ at $c = (t - f) \bmod 3 = 2$ and dicodon usage is strongly non-uniform; symmetrically, inside every antisense ($-$) ORF of offset $f \in \{0, 1, 2\}$, reverse-complement stop triplets $\{\text{TTA}, \text{CTA}, \text{TCA}\}$ terminate at the **exact same forward phase** $c = (t - f) \bmod 3 = 2$, while the antisense wobble base maps to $c = 0$ ($\text{bioPhase} = 2 - c$).

### Mechanism
1. **Exact Gap Bijection ($\psi: \mathbb{Z}_{\ge 0} \leftrightarrow \mathbb{Z}_{\ge 0} \times \{0, 1, 2\}$):**

$$q = \lfloor \text{gap} / 3 \rfloor, \qquad r = (\text{gap} + 1) \bmod 3, \qquad \text{gap} = 3q + ((r + 2) \bmod 3)$$

   By the chain rule of Shannon entropy, $H(\text{gap} \mid \text{regime}) = H(q \mid \text{regime}) + H(r \mid q, \text{regime})$. When $r$ is uniform (non-coding or i.i.d. sequence), reducing $\text{gap} \mapsto q$ saves $\log_2 3$ bits while $r$ costs $\log_2 3$ bits (**zero non-coding penalty**); when purifying selection skews $r = 0$, $\Delta H = \log_2 3 - H(r \mid q, \text{regime}) > 0$ bits/SNP is pure net gain.
2. **6-State Causal Open-Reading-Frame Hidden Markov Model (`6F-ORF-HMM` in [`FastNoveltyCoder`](../src/hypergrass/FastNoveltyCoder.java)):**
   Maintains 6 causal integer fixed-point log-posteriors $L_t(s)$ over hidden reading frames $s \in \{+0, +1, +2, -0, -1, -2\}$ with a gene-boundary switching floor $\max_{s'} L_t(s') - L_t(s) \le \Lambda_{\text{switch}}$ ($\approx 6.5$ nats, corresponding to $p_{\text{switch}} \approx 1.5 \times 10^{-3}$ at ORF boundaries), seeded only by universal stop-triplet depletion at $c = 2$. Codon and dicodon tables (`orfPlus2/3/5`, `orfMinus2/3/5`, grouped into 24-byte L1 cache lines `[((ctx * 3) + c) << 2]`) marginalize predictions across all 6 states via `EXP_NEG_LUT`, gate table updates by MAP posterior confidence ($\ge 4/9$), and route the inferred biological wobble indicator ($\text{bioPhase} = c$ on $+$, $2 - c$ on $-$) into the 8-regime logistic mixer and SSE calibrator.

### Falsification Test
1. Prove bijection $\psi^{-1}(\psi(g)) = g$ for all $g \in [0, 2 \times 10^6]$.
2. Verify codon-periodic homologous pairs compress strictly smaller with $\psi$ than flat-gap ablation.
3. Construct a multi-gene synthetic chromosome with random non-multiple-of-3 intergenic spacers and alternating $+/-$ strand ORFs: `FastNoveltyCoder` (`6F-ORF-HMM`) must compress the multi-frame coding chromosome significantly smaller than a frame-scrambled control with identical dinucleotide composition.

---

## 4. Law III — Conditional Monotonicity, Horizon-Bounded Two-Gear Coalescent Copying & Micro-Indel Continuity Law

### Mathematical & Biological Invariant
1. **Conditional Monotonicity (Law IIIa):** For nested reference sets $\mathcal{R}_1 \subseteq \mathcal{R}_2$, ideal conditional description length obeys monotonicity up to panel-index addressing:

$$L^*(X \mid \mathcal{R}_2) \le L^*(X \mid \mathcal{R}_1) + \Delta L_{\text{addr}}(|\mathcal{R}_2|, |\mathcal{R}_1|)$$

   Bounded-bucket FIFO/LIFO seed tables with small candidate caps violate conditional monotonicity whenever conserved $k$-mers in later panel genomes evict a closer twin genome in $\mathcal{R}_1$ (observed empirically in `Rcmdl` v2 on *M. tuberculosis* `03.seq`: `261` mismatches against 23 genomes vs. `59` mismatches against twin `08.seq` alone; `EIMDL` reduces `03.seq` against all 23 genomes to **`18` substitutions and `256 B`**).
2. **Horizon-Bounded Two-Gear Coalescent Spectral Gap (Law IIIb):** Bacterial and recombining population-genetic genealogies exhibit an ~20x divergence spectral gap between:
   - **Vertical clonal inheritance:** $\theta_{\text{clonal}} \le 5 \times 10^{-4}$ ($\le 1$ SNP per $2,000\text{--}50,000\text{ bp}$), and
   - **Horizontal homologous recombination / diverged homolog imports:** $\theta_{\text{recomb}} \in [3 \times 10^{-3}, 1.2 \times 10^{-1}]$ ($\ge 1$ SNP per $10\text{--}300\text{ bp}$).

   A naive single-gear approximate matcher calibrated against the single-genome novelty entropy ($R_{\text{lit}} \approx 1.865\text{ bpb}$, tolerating up to $\sim 15\%$ SNP divergence) commits catastrophic **recombination bulldozing** in multi-donor panels: when donor $d_A$ hits a recombination scar or mobile-element boundary while donor $d_B$ carries an exact $0$-SNP clonal haplotype $200\text{ bp}$ downstream, extending $d_A$ across the boundary at $5\text{--}12\%$ SNP density beats $1.865\text{ bpb}$ locally (`~0.85 bpb`) yet destroys $d_B$'s downstream $0.001\text{ bpb}$ clonal run (`29,623` false SNPs on *E. coli* `00.seq` before horizon bounding vs. `5,706` true SNPs with Law IIIb, cutting archive size from `56,813 B` to **`40,588 B`**).
3. **Affine Syntenic Micro-Indel Continuity (Law IIIc):** When a localized insertion or deletion $(\Delta_{\text{tgt}}, \Delta_{\text{ref}})$ occurs inside a syntenic block on donor $(g^*, o^*)$, breaking the match pays full global addressing $L_{\text{break}} \approx L(\text{event}) + L(g^*, o^*) + L(\text{pos}) + L(\text{len}) \approx 42\text{--}60$ bits, whereas an **in-match syntenic micro-indel splice** pays only:

$$L_{\text{splice}}(\Delta_{\text{tgt}}, \Delta_{\text{ref}}) = L(\text{spliceFlag} \mid \text{ctx}_8) + L(\Delta_{\text{tgt}}, \Delta_{\text{ref}}) + R_{\text{lit}} \cdot \Delta_{\text{tgt}}$$

   conditioned on an 8-state `(kind, spanBin)` `hasSplice` context prior so clonal exact chains pay $<0.02\text{ bits}$ overhead when unspliced while spliced homologs maintain megabase syntenic donor lock.

### Mechanism
1. **Strand-Symmetric Canonical Open Syncmer CSR Index** $(k=16, s=10, w=6, t=3)$ built via a **2-pass zero-materialization SWAR builder** with sticky active-donor collinear/indel probing (`mtf[0..1]`), distinct-donor bucket preservation, and 32-bit exact signature tags (`sig16` + CSR bucket).
2. **Two-Gear Coalescent Parser + O(1)-Amortized Clonal Re-Entry Horizon (`findClonalReentryHorizon` in [`EimdlMosaicCodec`](../src/hypergrass/EimdlMosaicCodec.java)):**
   - **Gear 1 (`clonalGear`):** Prices extension against the multi-donor clonal baseline (`baseGain = 0.008` bits/bp, `snpPenalty = 11.2` bits/SNP, `>= 128 bp` exact tail guard), immediately locking onto exact/near-clonal donor haplotypes and terminating cleanly at recombination boundaries.
   - **Gear 2 (`!clonalGear`, Horizon-Bounded Homolog Gear):** Invoked only when no clonal donor covers `pos`; first scans ahead via `horizonCache` to locate the earliest coordinate `gear2Limit > pos` where any panel donor resumes an exact clonal haplotype of length $\ge 192\text{ bp}$, strictly clamping Gear 2 (`[pos .. gear2Limit)`) so diverged homologs ($V_4$ Kimura/CpG + wobble + micro-indel splices) compress recombinant islands without ever bulldozing downstream clonal blocks.
3. **Bounded-window Syntenic Micro-Indel Lookahead** re-anchoring $(g^*, o^*)$ across pure deletions ($\Delta_{\text{tgt}} = 0, 1 \le \Delta_{\text{ref}} \le 256$), pure insertions ($1 \le \Delta_{\text{tgt}} \le 64, \Delta_{\text{ref}} = 0$), and compound micro-indels under exact MDL pricing.

### Falsification Test
1. **Conditional Monotonicity Test:** Adding distractors to a panel containing a twin donor must never degrade compressed payload by more than $\log_2(|\mathcal{R}_2|/|\mathcal{R}_1|)$ bits per token (and on *M. tuberculosis* `03.seq`, 23-panel LOO substitution count `18` is strictly below single-twin `08.seq` mismatch count `59`).
2. **Clonal Horizon & Distant-Homolog Ablation Test:** On real 23-donor *E. coli* `00.seq` LOO and single-distant-donor (`00.seq` vs `20.seq`, `54,124` SNPs, `619` indel splices), ablating Kimura $V_4$, codon wobble, or approximate homology increases archive bytes (`197,163 B` full EIMDL vs. `203,340 B` no-$V_4$, `199,024 B` no-wobble, `205,071 B` no-$V_4$/no-wobble/no-splice, and `238,550 B` exact-only); disabling micro-indel splicing alone fragments synteny from `1,517` to `1,945` panel tokens (`196,936 B`, `-227 B`), reflecting the honest Law VIa/VIb trade-off where cheaper short-gap and literal coding offsets long-tail splice signaling on this strain.
3. **Indel Continuity Test:** Injecting random 1–20 bp insertions and deletions into a homologous chromosome reduces compressed bytes by $>30\%$ compared to splice-disabled parsing (`6,818 B` vs `10,067 B`, `-32.3%` in `VerifyCompressionLaws`).

---

## 5. Law IV — Universal 5-Regime Existing & Non-Existing Sequence Completeness Law

To compress **every existing and non-existing genome** without catastrophic failure on any input distribution, the EIMDL parser evaluates five competing generative mechanisms under a strict MDL bit-budget:

$$L_{\text{EIMDL}}(X \mid \mathcal{R}) = \min \Big\{ 2|X| + c_0, \; \min_{\mathcal{P}} \sum_{p \in \mathcal{P}} L_{\text{phrase}}(p) \Big\}$$

where phrases $p \in \mathcal{P}$ draw from:
1. **Panel Syntenic Homology** ($\pm$ strand, exact + $V_4$ Kimura/CpG + wobble + micro-indel splices),
2. **Bidirectional Self-Homology** ($\pm$ strand, exact + approximate diverged paralogs/inverted repeats),
3. **Tandem Periodicity Attractor** (overlapping period $1 \le P \le 64$ run-exponent copies for satellites, telomeres, STRs, homopolymers),
4. **Multi-Order Cache-Resident Context + Repeat Mixer (`FastCRM`)** for novel coding, non-coding, Markov, or xenobiological sequence,
5. **Raw 2-Bit Safety Bypass** guaranteeing $L_{\text{DNA}}(X) \le \lceil 2|X|/8 \rceil + O(1)$ bytes even on cryptographically random sequences ($H_0 = 2.0000$ bpb).

---

## 6. Law V — Generalized Amdahl-Memory Streaming Scaling Law

### Bottleneck Diagnosis in v1
In [`HyperGraSS.compressFileBounded`](../src/hypergrass/HyperGraSS.java) ([SCALING.md](SCALING.md)), four sequential file passes plus a single-threaded 9.43 Mbp `QCRM` anchor page created a serial fraction $f_{\text{serial}} \approx 0.62$, capping multi-core speedup at $S_{\max} = 1 / f_{\text{serial}} \approx 1.61\times$ (`1.15x` measured at 2 threads), while open-addressed `SeedHashTable` consumed $7.1\text{ B/base}$.

### Fused Parallel & Memory Law
By:
1. **Single-pass streaming factorization (`tryStreamFlatFile` + incremental `concatenateAndConsume`)** directly into one contiguous 2-bit word-packed DNA bank ($0.25\text{ B/base}$ with zero duplicate copies),
2. **Zero-materialization 2-pass strand-symmetric canonical open-syncmer CSR indexing** ($\approx \frac{1}{w \cdot 2^s}(4 + 2) + 0.25\text{ B/base}$: $1.58\text{ B/base}$ full at $s=0$, $0.93\text{ B/base}$ at $s=1$ inside `-Xmx256m`/`-Xmx224m`, and $0.59\text{ B/base}$ at $s=2$ inside `-Xmx192m`/`-Xmx160m`/`-Xmx128m`, keeping **100% of donors indexed** with $0\text{ bytes}$ of intermediate per-donor syncmer arrays),
3. **Rank-ordered lock-free parallel DAG mosaic parsing** (every genome/unit $r \in \{0, \dots, N-1\}$ slices its target view on demand from the immutable 2-bit bank and parses concurrently across $P$ worker threads against strictly causal ancestors $< r$), and
4. **Warm-primed L3-cache-sized parallel `FastCRM` novelty stripes** (encoding a single shared warm-up prefix once, cloning the adapted Markov/6F-ORF-HMM snapshot, and right-sizing per-worker repeat rings to fit CPU L3 cache across concurrent workers),

the end-to-end wall time and peak live heap obey:

$$T(P) = T_{\text{stream+csr}} + \frac{T_{\text{mosaic}} + T_{\text{novelty}}}{P} \implies S(P) = \frac{1}{f_{\text{serial}} + (1 - f_{\text{serial}})/P}$$

$$M_{\text{peak}}(L_{\text{total}}, s, P) = L_{\text{total}} \cdot \left(0.25 + \frac{6}{w \cdot 2^s}\right)\text{ B} + P \cdot B_{\text{worker}}$$

where Stage 1 releases the entire CSR index before Stage 2 (`FastCRM`) begins, preventing cross-stage heap stacking.

### Falsification Test
1. Parallel determinism & scaling: `compress` across `--threads 1..12` must produce byte-exact SHA-256 decompressed reconstructions with $\le 0.5\%$ multi-stripe size delta and $>35\times$ end-to-end speedup over `HyperGraSS v1` (`887.3 s`).
2. Strict bounded-heap gate: Compressing and restoring the full 24-genome *E. coli* collection (`120,108,223 B`) must complete and pass fresh-process SHA-256 verification under `-Xmx256m`, `-Xmx224m`, `-Xmx192m`, `-Xmx160m`, and **`-Xmx128m`** JVM heap caps while still beating unbounded `HyperGraSS v1` (`3,810,409 B`), `JARVIS3` (`3,896,361 B`), and `RC-MDL` (`3,972,975 B`).

---

## 7. Law VI — Hierarchical Pitman-Yor Backoff & Two-Adic Tapered Adaptation Law

### Mathematical & Algorithmic Invariant
1. **Two-Adic Tapered Binary Range Adaptation (Law VIa):** In shift-based binary arithmetic range coding over probability scale $M = 2^{12} = 4096$ with shift $k = 5$, standard integer updates $p_{t+1} = p_t - (p_t \gg 5)$ stall at the integer truncation fixed point $p_{\min}^{\text{classic}} = 2^5 - 1 = 31 / 4096 \approx 0.757\%$, lower-bounding the asymptotic code length of high-certainty structural flags (`isMtf0`, `useCollinearPred`, `hasSplice`, `hasSnps`) at $-\log_2(1 - 31/4096) \approx 0.0110\text{ bits/symbol}$. Introducing a symmetric two-adic rounding taper $\tau = 20$ ($p_{t+1} = p_t - ((p_t + \tau) \gg 5)$ on $b = 0$; $p_{t+1} = p_t + ((M - p_t + \tau) \gg 5)$ on $b = 1$) lowers the asymptotic fixed point to:

$$p_{\min}(\tau) = 31 - \tau = 11 / 4096 \approx 0.269\%, \qquad L_{\min}(\tau) = -\log_2\!\left(1 - \frac{11}{4096}\right) \approx 0.00388 \text{ bits/symbol}$$

   cutting steady-state binary flag entropy by **$2.8\times$** at zero extra clock cycles, while joint Law I $\times$ Law II substitution conditioning (`(kimuraCtx << 1) | prevMod3Zero`) and slot-isolated low-bit contexts (`lowBitProbs[8][128]`, `posDeltaSlots[2][64]`) eliminate cross-stream probability pollution.
2. **Hierarchical Pitman-Yor Suffix-Tree Backoff & Post-SNP Homology Recovery (Law VIb):** Because genomic $k$-mer frequencies follow a heavy-tailed Pitman-Yor / Zipf distribution, high-order Markov and codon models ($k \in \{1, \dots, 9, \text{Sparse-}12, \text{ORF-3}, \text{ORF-5}\}$) spend most updates in the low-count regime. Instead of backing off sparse counts to the flat uniform baseline $1/4$ ($2.0\text{ bpb}$), each order-$k$ estimator interpolates directly with its already-computed order-$(k-1)$ posterior $P_{k-1}(a \mid c_{k-1})$:

$$P_k(a \mid c_k) = \frac{n_k(c_k, a) + \alpha_k \, P_{k-1}(a \mid c_{k-1})}{n_k(c_k) + \alpha_k}$$

   recovering $(2.0 - H_{k-1})\text{ bpb}$ on every novel or rare context transition, while `MicroRepeatHead` conditions post-mismatch confidence on cumulative homologous match length (`totalMatched`) once a single SNP is cleared (`misses == 0`).
3. **Reciprocal-Multiplication Mixer Normalization (Law VIc):** Replacing variable-denominator 64-bit hardware integer divisions (`idiv`) inside `FastCRM`'s per-symbol prediction loop with an immutable 16,385-entry reciprocal lookup table $R(d) = \lfloor (2^{32} + d - 1)/d \rfloor$ (`(num * R(d)) >>> 32`) eliminates hardware division latency across all 17 experts.

### Falsification Test
1. Verify algebraically that two-adic tapered adaptation ($\tau = 20$) converges to exact fixed point $p_{\min} = 11 / 4096$ vs. $31 / 4096$ for classic shift-5 updates ($>2.8\times$ lower residual probability floor) and compresses coding+repeat novelty test streams to `1.0897 bpb` in `VerifyCompressionLaws` (`Novelty=1.0897 bpb`).
2. Verify in fresh-process Mode S execution (`mode_s_and_synthetic.json`) that `FastNoveltyCoder` compresses real *E. coli* `00.fa` (`4,856,085 B`) down to **`1,122,777 B`** (`1.8497 bpb`) in **`6.75 s`** in-process encode (`6.92 s` wall) / **`6.57 s`** in-process decode (`6.79 s` wall) and *M. tuberculosis* `00.seq` (`4,417,521 B`) down to **`969,382 B`** (`1.7555 bpb`) in **`7.49 s`** encode / **`6.25 s`** decode.

---

## 8. Law VII — Disjoint Index/Superstream Lifecycle & Sub-Linear Bounded-Heap Law (`-Xmx128m`)

### Mathematical & Architectural Invariant
In two-stage pangenome collection compression (Stage 1: rank-ordered parallel DAG mosaic parsing over CSR syncmer index $\mathcal{I}_{\text{CSR}}$; Stage 2: unified pan-genome novelty superstream coding over $\mathcal{L}_{\text{super}}$), the memory lifespans of $\mathcal{I}_{\text{CSR}}$ and $\mathcal{L}_{\text{super}}$ are **strictly disjoint** once per-rank mosaic bitstreams are emitted:

$$\text{Live}(\text{Stage 1}) = \{\mathcal{B}_{2\text{bit}}, \, \mathcal{I}_{\text{CSR}}, \, \{\mathcal{S}_r\}\}, \qquad \text{Live}(\text{Stage 2}) = \{\mathcal{M}_{\text{out}}, \, \mathcal{L}_{\text{super}}, \, \{\mathcal{W}_p\}_{p=1}^{P}\}$$

$$\implies M_{\text{peak}} = \max\!\big(M_{\text{Stage 1}}, \, M_{\text{Stage 2}}\big) \ll M_{\text{Stage 1}} + M_{\text{Stage 2}}$$

By:
1. Explicitly releasing `sharedDagIdx`, `mstParentRank`, and `rankDonorSimilarityOrder` immediately after Stage 1 mosaic emission and incrementally nulling `rankMosaics[r]` and `unitNovSlices[r]` as $\mathcal{L}_{\text{super}}$ is assembled,
2. Narrowing per-worker target syncmer canonical hashes from 64-bit `long[]` to 32-bit `int[]` (`TargetSyncmers.canonKmers`, `-38%` per-worker memory),
3. Decoupling `FastCRM` warm-primed stripe count (`workerCount`) from `ForkJoinPool` concurrency (`poolThreads`) and dynamically throttling `buildThreads`/`workerThreads` $\le 2$ under `extremeTight` ($\le 145\text{ MiB}$) while enabling $s=2$ (`ultraTight` $\le 205\text{ MiB}$) so **full `panSuper` (`Order-8`, `Order-9`, `Sparse-12`) context modeling remains active even when `--bounded` is enabled**, and
4. Replacing 3.43 million per-line kernel writes and JNI CRC updates in `FastaStreamFactorizer.reconstructTo` with a 64 KB staging buffer and parallelizing Multi-FASTA chunk factorization across CPU cores,

pure-JVM Mode C scales down to **`-Xmx160m` and `-Xmx128m`** ($1.06\text{ B/base}$ total JVM heap limit on the $120,108,223\text{ B}$ 24-genome *E. coli* collection), hitting **`3,422,314 B`** (`0.2279 bpb`, `+0.99%`) at `-Xmx256m`/`-Xmx224m` and **`3,452,703 B`** (`0.2300 bpb`, `+1.89%`) at `-Xmx192m`/`-Xmx160m`/`-Xmx128m` (all within the `<2%` unbounded delta gate).

### Falsification Test
1. Verify in `VerifyCompressionLaws` (`Unbounded=27257 B, Bounded=27257 B, delta=0.000%`, part of **`171,941` passing assertions**) that Mode C unbounded and `--bounded` archives round-trip with identical `SHA-256` and $<2\%$ archive size difference.
2. Verify in fresh-process CLI execution (`run_fusion_benchmarks.py`) that the $120,108,223\text{ B}$ 24-genome *E. coli* collection compresses and decompresses with `SHA-256 = True` at `-Xmx256m`, `-Xmx224m`, `-Xmx192m`, **`-Xmx160m`**, and **`-Xmx128m`**.

---

## 9. Architectural Transfer from OpenAI's October 2026 ETH Refutation Preprint

Three structural proof techniques from OpenAI's 67-page ETH refutation (*"The Exponential Time Hypothesis is False"*, 17 Oct 2026 / rev. 25 Oct 2026; `openai2026eth` in [`paper/references.bib`](../paper/references.bib)) directly inspired Laws VI–VII, subject to an explicit asymptotic scope boundary:

1. **Honest Asymptotic Boundary (§1.1.3 OVC Caveat):** Theorem 1.1 solves `#k-SAT` in $(n+m)^{O_k(1)} 2^{O_k(n \log\log n / \sqrt{\log n})}$ time and polynomial space, but Section 1.1.3 explicitly cautions that polynomial-time fine-grained complexity lower bounds tied to the Orthogonal Vectors Conjecture (OVC) — including quadratic-time edit distance and Longest Common Subsequence (LCS) bottlenecks — remain intact. Accordingly, subexponential `#k-SAT` algebraic reductions do **not** bypass polynomial-time string-matching complexity in biological sequence compression.
2. **Shared Linear Summaries (§2.4 Def. 2.1, §2.7 eqs. 2.14–2.16, §4.2 Lemma 4.1) $\to$ Law VIb:** Computing shared hierarchical summaries rather than independent sparse tables maps directly to `FastNoveltyCoder`'s incremental Hierarchical Pitman-Yor suffix backoff chain (`ord0 -> ... -> ord9 -> sparse12` and `orf2 -> orf3 -> orf5`).
3. **Two-Adic Valuation Tapering (Appendix D.4.3–D.4.4, Lemma D.5 eq. D.10) $\to$ Law VIa:** Replacing uniform integer degree/step floors with two-adic multiplicity tapering inspired the symmetric `tau = 20` rounding taper in `BitRangeEnc`/`BitRangeDec` (`(p + 20) >>> 5` and `(4096 - p + 20) >>> 5`), lowering the stalling floor from `31/4096` to `11/4096`.
4. **Disjoint Partitioning & Sequential Polynomial-Space Accumulation (Lemma 3.2, Prop. 3.4, §5.2 $\mathbb{Q}[V]/(V^\ell)$, Prop. 6.4, App. A & §D.4.7) $\to$ Law VII:** Decomposing global state into disjoint partitions accumulated sequentially in polynomial space motivated strict Stage-1 (`sharedDagIdx`, `mstParentRank`, `rankDonorSimilarityOrder`) vs. Stage-2 (`rankMosaics`, `unitNovSlices`, pooled `FastCRM`) buffer deallocation and `ultraTight`/`extremeTight` scaling down to `-Xmx128m`.

