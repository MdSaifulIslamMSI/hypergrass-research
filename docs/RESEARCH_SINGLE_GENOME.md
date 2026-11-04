# Beyond HyperGraSS: lossless single-genome compression from genome and hardware architecture

Status: research note. No code changes and no new benchmark runs are reported here.
The HyperGraSS profile campaign is closed. This note asks a different question: if
HyperGraSS is retained as a pangenome/reference engine, what algorithm class should
replace or extend it for the single-genome case?

Evidence tags used below:

- **[M]** measured in this project and independently verified (SHA-256 plus a separate
  byte comparison).
- **[R]** reported by the cited paper or project, not independently reproduced here.
- **[D]** derived arithmetic from reported numbers.
- **[P]** projection or hypothesis; explicitly not measured.

The short version:

1. Self-contained single-genome compression is clustered around **1.38–1.9 bits/base**
   on human-scale DNA. The best published human result found in this investigation is
   JARVIS3 at **1.381 bits/base** (538,155,679 B for 3,117,292,120 bases) [R].
   HyperGraSS measures **1.650 bits/base** on the project's 100 Mb chr1 input [M].
2. A 100:1 archive is **0.08 bits/base**. That is roughly **17x below** the best
   self-contained rate and no published self-contained codec comes close. The 100:1
   regime is a **reference/pangenome-conditioned** regime, where RLZ, ERGC, AGC and
   MBGC report 0.1, 0.019, 0.014–0.040 and 0.006 bits/base respectively [R]; all but
   RLZ are at or beyond 100:1.
3. The strongest next design is a **dual-mode MDL mosaic**: parse the target into
   reference copies, self-copies, grammar expansions and literals; price every phrase
   with the same entropy model that codes it; and code only the residual with an
   expensive context mixer. Most bases take the fast copy path; the slow model is
   reserved for the divergent minority.
4. The expected envelope is **5–15% better than JARVIS3 in self-contained mode** and
   **10–100x better with a panel**, at **10–100x the speed of a full context model**
   because most bases are copied rather than modeled [P]. The projections are
   hypotheses with explicit falsification gates in §9.
5. Build order: **RC-MDL** prototype (reference-conditioned parser), **GAH**
   self-contained grammar/attractor track, **HNP** hardware-native rANS backend,
   **NCM** neural copy model as a research bet, **PODA** for collections.

---

## 1. What we are trying to beat

### 1.1 The measured HyperGraSS baseline

The 100 Mb chr1 input is the project's single-genome reference point. The input
SHA-256 is `0c1aeb4eaf5c5ca33c3bdba9d1f7236b480e9c17415a39e82057904bbcddee5e`;
305,111 non-ACGT bases were masked to `A`. All rows below are byte-exact [M].

| Profile | Archive (B) | bits/base | Compress (s) | Peak RSS (MiB) | Exact |
|---|---:|---:|---:|---:|:---:|
| `--single` | 20,626,877 | 1.650 | 355.0 | 396 | yes |
| `--threads 2` | 20,635,964 | 1.651 | 70.1 | 586 | yes |
| `--threads 4` | 20,649,693 | 1.652 | 124.9 | 617 | yes |
| `--threads 8` | 20,677,028 | 1.654 | 225.9 | 697 | yes |
| default | 21,056,843 | 1.685 | 152.9 | 967 | yes |
| `--lean` | 22,788,010 | 1.823 | 40.5 | 94 | yes |
| `--pangenome` | 22,735,705 | 1.819 | 744.2 | 962 | yes |
| GraSS | 22,264,066 | 1.781 | 137.5 | 893 | yes |
| `--bclzg` | 21,120,171 | 1.690 | 181.1 | 991 | yes |
| `--achmc` | 21,151,440 | 1.692 | 632.6 | 1,438 | yes |
| `--qcrm` | 22,788,010 | 1.823 | 40.9 | 94 | yes |
| `--auto` | 22,735,705 | 1.819 | 208.7 | 981 | yes |

`--pangenome` and `--auto` are byte-identical; `--qcrm` and `--lean` are byte-identical.
The best serial result is **7.35% below GraSS** on the same input. The best practical
trade is `--threads 2`: 20,635,964 B in 70.1 s. The `--bclzg` ablation is the strongest
forced-payload variant at 21,120,171 B, still 2.39% above serial `--single`. The
`--achmc` ablation is dominated: 0.148% fewer bytes than `--bclzg` for 3.49x the
compress time and 1.45x the RSS.

The 120,108,223-byte, 24-strain *E. coli* pangenome tells a different story. HyperGraSS
is **2.21% smaller than JARVIS3 level 7** there (3,810,409 B vs 3,896,361 B), and the
opt-in GraalVM native image passes a strict 256 MiB gate at 3,864,464 B [M]. That is
the pangenome strength the project already has. The single-genome case is the weak
case: 1.650 bits/base is a GraSS-class result, not a JARVIS3-class result.

### 1.2 What "single genome" means

Three modes must be kept separate, because they have different information budgets:

- **Mode S — self-contained.** The archive and decoder are sufficient; no external
  data. This is the only mode that is unconditionally lossless end-to-end.
- **Mode R — reference-conditioned.** The archive plus a pinned reference panel
  reconstruct the target. Lossless is conditional on the panel. The panel digest and
  version must be in the archive; if the panel is unavailable, decoding is undefined
  by contract.
- **Mode C — collection.** Several genomes are compressed jointly. The unit is the
  collection, not one genome. AGC, MBGC and phylogenetic compression live here.

The user's requirement is end-to-end lossless. Mode S satisfies it without caveats.
Mode R can satisfy it if the reference is bundled or pinned by digest; otherwise it is
"lossless given the reference," which is a different contract. The recommended format
therefore carries an explicit `reference_required` flag and a `reference_digest`; a
self-contained fallback stores the reference when the application demands it.

### 1.3 What 100:1 actually costs

For a genome of `n` bases, 100:1 means `n/100` bytes, or **0.08 bits/base**. The
arithmetic is unforgiving:

| Input | 100:1 archive | bits/base | Best archive found (mode) | Gap |
|---|---:|---:|---:|---:|
| 3.117 Gb human | 31.2 MB | 0.080 | 1.381 (JARVIS3, S) | 17.3x |
| 120.1 Mb E. coli panel | 1.20 MB | 0.080 | 0.254 (HyperGraSS, C) | 3.2x |
| 100 Mb chr1 | 1.00 MB | 0.080 | 1.650 (HyperGraSS, S) | 20.6x |

`S` is a self-contained single sequence; `C` is a jointly coded collection. The *E. coli*
row is a collection result, which is why it comes closest to 100:1. The table compares
the target archive with the best archive found for the same input, not with a universal
codec ranking.

The E. coli pangenome is the closest case because the 24 strains share most of their
sequence. Even there, the measured entropy decomposition in
[`ENTROPY-FLOOR.md`](ENTROPY-FLOOR.md) rejects 100:1: the anchor channel alone is
2,257,234 B (1.88x the target), and the non-anchor channels alone are 1,553,175 B
(1.29x). The 100:1 target is not a tuning problem on that corpus; it is a modeling
problem.

For a single human genome, 100:1 is **reference territory**. AGC's 5.4 MB incremental
cost per human assembly is 0.014 bits/base, about 5.7x better than 100:1; AGC's 1.45 GB
for 96 assemblies is 0.040 bits/base, about 2x better than 100:1 [R][D]. So 100:1 is
already demonstrated on human collections, but only with a panel. A self-contained
codec would need a 17x improvement over JARVIS3. Nothing in the literature suggests
that is available from a better context model alone.

---

## 2. Method of the investigation

This is a targeted literature synthesis, not a formal systematic review. The working
backends were Europe PMC REST, arXiv, OpenAlex, Crossref, publisher pages, and the
Cobilab `HumanGenome` benchmark repository. Queries covered:

- reference-based and relative Lempel-Ziv genome compression;
- grammar-based compression, RePair, string attractors, smallest grammars;
- BWT, run-length BWT, r-index, prefix-free parsing, graph BWT;
- pangenome graphs, variation graphs, GBWT, AGC, MBGC/MBGC2;
- finite-context models, context mixing, substitution-tolerant models, neural DNA
  compressors, genomic language models;
- GPU/FPGA/SIMD acceleration, rANS, parallel decompression, compressed-resident
  decoding, pangenome serving;
- DNA storage codes and fountain codes where they bear on constrained sequence coding.

Inclusion criteria were: lossless DNA/genome compression; a described mechanism; at
least one quantitative result or a reusable data structure; and enough detail to
compare ratio, speed or memory. Lossy methods were excluded. Read-compression methods
(FASTQ/CRAM) were included only when they illustrated a reference or hardware
mechanism that transfers to assembled genomes. Two prior reviews [47][48] were used as
entry points and cross-checks.

Three caveats apply throughout:

- Author-reported numbers use different hardware, thread counts, preprocessing and
  input lengths. Bits/base is the least unstable comparison, and even that is only
  valid on the same input.
- Several 2025–2026 results are preprints or technical reports. They are marked [R]
  and should be re-run before being treated as established.
- The HyperGraSS numbers are the only numbers in this note that were independently
  verified in this project.

---

## 3. The mathematical object

### 3.1 Shannon, MDL, and the one-sample problem

For a source `X` with distribution `p`, an arithmetic code needs about
`-log2 p(x)` bits for a particular sequence `x`. The expected length is bounded below
by the entropy `H(X)` [50]. For a genome we do not know `p`, and we usually have one sample.
The practical substitute is the **minimum description length (MDL)** two-part code:

```text
L(x) = min over models M of [ L(M) + L(x | M) ].
```

`L(M)` is the cost of describing the model; `L(x | M)` is the cost of the data under
the model. A model that is too simple underfits; a model that is too complex spends
too many bits on `L(M)`. This is the correct frame for single-genome compression
because it makes model cost explicit. It also explains why "just use a bigger model"
is not automatically better: the model has to pay for itself in the archive.

### 3.2 Finite-context entropy and why it saturates

An order-`k` finite-context model estimates `p(x_i | x_{i-k..i-1})`. The plug-in
entropy `H_k` falls as `k` grows, but the number of possible contexts grows as `4^k`.
At `k = 12` there are 16.7 million contexts; at `k = 16`, 4.3 billion. A 3 Gb genome
has only 3 billion positions, so high-order contexts are sparsely sampled and the
plug-in estimate is optimistically biased. Smoothing, backoff and context mixing are
not optional; they are the difference between a model and a lookup table.

The sequence-space data make the point directly. In the 2008 analysis of seven
genomes, **more than 98% of all 12-mers appear in vertebrate genomes**, while **less
than 2% of 19-mers do** [46]. High-order k-mer contexts are saturated: many 12-mers are
seen, but their continuations are not reliably estimated from one genome. This is why
JARVIS3 and GeCo3 use multiple context orders, substitution-tolerant contexts and
mixers rather than a single high-order table.

### 3.3 Repeats, grammar, and attractors

A finite-context model cannot see a 10 kb repeat that starts 1 Mb away. Long-range
redundancy needs a different mechanism:

- **LZ77 / relative LZ.** A phrase copies a previous substring (or a substring of a
  reference). The code is the distance and length. This is the mechanism behind
  HyperGraSS's anchors, RLZ, AGC and MBGC.
- **Grammar.** A straight-line grammar `G` generates `x` by repeated substitution.
  Rules capture nested and overlapping repeats, tandem arrays and copy-number
  structure. The smallest grammar is NP-hard; RePair and its variants are the practical
  heuristics.
- **String attractor.** A set of positions `Γ` such that every distinct substring of
  `x` has an occurrence crossing `Γ`. The minimum attractor size `γ*` lower-bounds
  many dictionary compressors. Computing `γ*` is NP-hard; LZ, run-length BWT, the
  smallest grammar and macro schemes all induce attractors within approximation
  factors [R].
- **BWT / r-index.** The Burrows-Wheeler transform sorts rotations; equal contexts
  become adjacent, and runs in the BWT capture repeats. The number of runs `r` is a
  repetitiveness measure; the r-index supports pattern matching in `O(r)` space [R].

These are not interchangeable. A grammar is good at nested tandem repeats; LZ is good
at long exact copies; a context model is good at short-range statistics. A serious
single-genome codec needs at least two of the three.

### 3.4 Reference-conditioned rate

If the target `x` is within edit distance `d` of a reference `R`, a relative LZ parse
can encode `x` as a sequence of copies and an edit script. The information cost is
roughly `O(d log(n/d))` bits plus the cost of the copy geometry. For human genomes,
`d` is on the order of millions of substitutions and indels, not billions; that is why
reference-based methods reach 0.01–0.1 bits/base while self-contained methods stay
near 1.4.

The practical lesson is not "references are magic." It is that **the reference supplies
the model that a single genome cannot estimate from itself**. The archive still has to
encode the differences, the copy geometry, and the parts of the target absent from the
panel. A good reference codec is therefore a good **delta codec**, not a lookup table.

### 3.5 Implications

1. **Self-contained compression has an empirical plateau near 1.3–1.4 bits/base**
   with current model classes. The exact floor is unknown, but no published method is
   an order of magnitude below JARVIS3.
2. **The large gains are in the prior, not the entropy coder.** A 5% better arithmetic
   coder is useful; a 10x better prior is transformative.
3. **The right architecture is hybrid.** Use copies for the long-range redundancy and
   a context model for the residual. A single model class leaves bits on the table.
4. **The parser and the coder must be co-designed.** A match that looks cheap in
   length but expensive in edit script is not cheap. The parser should minimize the
   actual code length, not the match count.

---

## 4. Taxonomy of lossless DNA compression

### 4.1 Statistical and finite-context models

**Mechanism.** Predict the next base from a bounded context; entropy-code the residual.
Enhancements include multiple context orders, substitution-tolerant contexts (allow
mismatches in the context), reverse-complement contexts, codon-phase models, and a
mixer that combines expert probabilities.

**Representatives.** JARVIS3 [3], GeCo3 [4], DNA-COMPACT [16], MFCompress [49], XM,
and the context-mixing family descended from PAQ.

**Strengths.** Strong on short-range composition, CpG/GC structure, codon bias and
moderately repetitive sequence. Constant memory in GeCo3; JARVIS3 has 40 presets and
fits under 4 GB except the maximal profiles [R].

**Weaknesses.** Blind to long repeats beyond the model order; slow because every base
is modeled; model cost is significant; high-order contexts are sparse.

**Reported human result.** JARVIS3 1.381 bits/base; GeCo3 1.425 bits/base [R]. On the
ArTh benchmark, JARVIS3 averages 1.930 bits/base across nine datasets, while PMKLC and
AgentGC-CP average 1.854 and 1.850 [R]. The ranking is dataset-dependent: JARVIS3 wins
ArTh at 1.855, but loses the nine-dataset average to the learning-based methods.

### 4.2 Dictionary and LZ methods

**Mechanism.** Replace repeated substrings with references. The reference can be the
same sequence (LZ77), a relative reference (RLZ), or a collection (AGC/MBGC).

**Representatives.** NAF [6], FFC [10], xz/LZMA2, zstd, RLZ [11], RLZAP [12], GDC2
[21], AGC [9], MBGC [7], MBGC2 [8].

**Strengths.** Fast decode; excellent for collections; naturally supports random
access if block-structured. NAF is 30–80x faster to decode than DELIMINATE/MFCompress
at a comparable ratio [R]. FFC is 4.7x/11.4x faster to compress than zstd/NAF and 3.5x/
2.7x faster to decompress, with a 23% better ratio than pzstd [R]. MBGC reports 1,265x
on 168,311 bacterial genomes (587 GB) at ~1,580 MB/s compression and ~780 MB/s
decompression with 8 threads [R]. AGC compresses 96 HPRC human assemblies (293.2 GB)
to 1.45 GB in ~12 minutes with 32 threads, about 200x [R].

**Weaknesses.** Match geometry and literals are often coded with fixed or weak models;
divergent regions and short repeats are left to a generic coder. Collection methods
require a collection; a single genome with no panel only has self-references.

### 4.3 Grammar and attractor methods

**Mechanism.** Build a grammar or a dictionary of repeated substrings; encode the
grammar, not the sequence. RePair repeatedly replaces the most frequent pair with a
nonterminal. Attractor-based methods aim at a small set of positions that "covers" all
distinct substrings.

**Representatives.** GraSS [5], RePair/MR-RePair, DNACompress, grammar-based
compression, and the string-attractor literature [24][25].

**Strengths.** Captures nested and tandem repeats, copy-number structure, and some
reverse-complement symmetry if the grammar allows it. GraSS reports weighted average
compression ratios of 4.5 for DNA and 19.6 for RNA, with a 246.8 s total DNA corpus
time [R].

**Weaknesses.** Grammar construction is expensive and memory-hungry; RePair is a
heuristic; grammar size is not the same as code length; GraSS is weak on human-scale
repeats (its 100 Mb chr1 archive is 22,264,066 B, 1.781 bits/base [M]).

**Opportunity.** A modern grammar codec that (a) allows reverse-complement rules,
(b) allows run-length/exponentiation rules for tandem arrays, (c) scores rules by MDL
gain including the entropy coder, and (d) encodes the grammar with a context model has
not been built for human-scale DNA. That is the GAH track.

### 4.4 BWT, r-index and suffix-based methods

**Mechanism.** Sort rotations, run-length encode the BWT, and encode the runs. The
r-index stores the run-length BWT plus samples to support pattern matching in space
proportional to the number of runs `r`.

**Representatives.** bzip2, Cox et al.'s BWT-based genome compression [23], the r-index
[26], prefix-free parsing for BWT construction [27], and graph BWT/GBWT for pangenomes
[28][30].

**Strengths.** Mature theory; `r` is a principled repetitiveness measure; supports
compressed indexing and random access. The r-index is optimal in `r` for locating
occurrences [R]. Prefix-free parsing makes BWT construction practical on large and
repetitive collections [R].

**Weaknesses.** BWT construction is memory-heavy (about 11.9 bytes/character with the
GPU prefix-free parser, 2,108 MiB input on a 24 GB RTX 3090 [R]); BWT alone is not the
best ratio for DNA; the r-index is an index, not a complete compressor. For a single
3 Gb genome, a full BWT is a 10–60 GB data structure unless it is blockwise or
external-memory.

**Opportunity.** BWT is best used as a **repetition analyzer and index**, not as the
whole codec. A blockwise BWT can identify repeat families and feed a grammar/parser;
the r-index can provide random access to the resulting archive.

### 4.5 Reference-based and delta methods

**Mechanism.** Encode the target as a set of matches to one or more references plus
literals and an edit script. Variants differ in match finding (exact vs approximate),
reference selection, edit encoding, and random access.

**Representatives.** RLZ [11], RLZAP [12], iDoComp [13], ERGC [14], GReEn [15],
DNA-COMPACT [16], HRCM [17], NRGC [18], SparkGC [19], LMSRGC [20], GDC2 [21], DNAzip
[22], AGC [9], MBGC2 [8].

**Strengths.** The only route to 0.01–0.1 bits/base. RLZ represents individual human
genomes at around 0.1 bits/base with rapid access [R]. iDoComp compresses 2,987 MB of
human data to ~64 MB against hg18 (0.171 bits/base); ERGC reduces the same pair to
~7 MB (0.019 bits/base), 9x better than iDoComp [R][D]. AGC reaches 0.040 bits/base on
96 human assemblies and 0.014 bits/base incremental per assembly [R][D].

**Weaknesses.** Requires a reference; performance is sensitive to divergence, assembly
quality and structural variation; approximate-match encoding can be lossy if not
verified; random access is often sacrificed. Many methods are chromosome-scale and do
not handle contigs, inversions or copy-number variation well.

**Opportunity.** The delta codec is where the ratio is. The open problem is not "find
matches" but "price matches, edits and literals with one model and parse optimally."

### 4.6 Pangenome graphs and graph indexes

**Mechanism.** Represent a population as a graph whose paths are haplotypes. Encode a
target as a path plus edits, or index the graph with a graph BWT/GBWT for search.

**Representatives.** vg [29], GBWT/haplotype-aware graph indexes [28], minigraph [54],
Minigraph-Cactus [56], ODGI [55], Movi [53], ropebwt3 [43], panvram [44].

**Strengths.** Captures SNPs, indels, inversions and structural variants in one
structure; supports graph-aware search and genotyping; can serve pangenome windows
from GPU memory. The HPRC panel [2] is the canonical human collection for this class.
panvram stores 558 haplotype assemblies in 16.642 GB (default q4k) or
11.634 GB compact, and serves ~740k windows/s [R]. ropebwt3 indexed 320 human genomes
in 65 hours with 170 GB peak RAM [R].

**Weaknesses.** Graph construction and indexing are expensive; graph formats and
coordinate systems are complex; compression is usually a by-product of the index, not
the primary objective; a graph is only useful if the target is represented in it.

**Opportunity.** A graph prior is the most powerful reference prior for human
pangenomes. The missing piece is a **rate-optimal path encoder**: given a graph and a
target, choose the path and edits that minimize the actual code length, not just the
alignment score.

### 4.7 Neural and language-model methods

**Mechanism.** Train an autoregressive model `p_θ(x_i | x_{<i})`; arithmetic-code the
sequence with the model's probabilities. Variants use CNNs, RNNs, transformers,
state-space models, BPE tokenization, and retrieval/memory.

**Representatives.** GeCo3's shallow MLP mixer [4], DeepGeCo [33], PMKLC [34],
AgentGC [35], DNAGPT2 [36], NNCP, and DNABERT-family models.

**Strengths.** Learn motifs, codon structure, regulatory grammar and feature-specific
predictability without hand engineering. DNAGPT2 reaches 1.470 bits/base on T2T, only
about 6% behind JARVIS3, and produces an information-content map where exons, introns,
intergenic regions and Alu repeats separate cleanly (1.875, 1.814, 1.767 and 0.631
bits/base respectively) [R].

**Weaknesses.** Model side information is large. DNAGPT2's 86M-parameter bf16 model is
about 170 MB; on a 3.1 Gb genome that is ~0.45 bits/base if not amortized, which would
make it worse than JARVIS3 [D]. It is also slow: 728 minutes on human vs 22 minutes for
MFCompress [R]. Long context does not automatically help; DNAGPT2's 32-token BPE model
outperformed longer-context genomic LMs in its benchmark, and the authors conclude that
practical DNA compression will likely require hybridization with explicit copy/repeat
mechanisms [R].

**Opportunity.** Neural models are best used as **probability experts** inside a hybrid,
especially for the literal/residual channel and for annotation-aware region models, not
as standalone single-genome codecs. A small, quantized, shared model can be amortized
across a collection; a 170 MB model cannot be amortized across one genome.

### 4.8 DNA storage codes

**Mechanism.** Encode bits into DNA with constraints (GC balance, no long homopolymers,
no repeats), error-correcting codes, and fountain codes for robustness.

**Representatives.** DNA Fountain [37], constrained codes, deletion/insertion codes.

**Strengths.** DNA Fountain approaches the information capacity per nucleotide and
achieves robust retrieval from sequencing [R]. The coding-theory ideas (erasure/
insertion correction, constrained coding) are mature.

**Weaknesses.** This is a different channel: synthesis and sequencing errors, not disk
storage. The codec is not a substitute for genome compression; it is relevant only if
the target storage medium is DNA.

**Opportunity.** Constrained coding and fountain codes are a source of ideas for
**robust archival formats**, not for the ratio problem on disk. Keep them out of the
critical path unless DNA storage is the deployment target.

### 4.9 Comparison

| Family | Long repeats | Short context | Reference prior | Speed | Single-genome ceiling |
|---|---|---|---|---|---|
| Finite-context (JARVIS3/GeCo3) | weak/medium | strong | no | very slow | ~1.38–1.43 bpb |
| LZ/dictionary (NAF/FFC) | strong | weak | no | fast | ~1.8–2.0 bpb |
| Grammar (GraSS/RePair) | strong | medium | no | medium | ~1.7–1.9 bpb |
| BWT/r-index | strong | medium | no | medium | ~1.7–1.9 bpb |
| Reference/delta (RLZ/AGC) | strong | weak | yes | fast | ~0.006–0.1 bpb |
| Pangenome graph | strong | medium | yes | medium | ~0.01–0.1 bpb |
| Neural/LM | medium | strong | optional | very slow | ~1.47 bpb |
| Hybrid (this proposal) | strong | strong | optional | fast copy path | target 1.3–1.6 / 0.02–0.08 |

The table is a synthesis, not a controlled benchmark. The "ceiling" column is the best
reported or measured point in that family, not a guarantee for every dataset.

---

## 5. Where the bits are: genome architecture

### 5.1 Repeat composition of the human genome

The T2T-CHM13 assembly is the best current reference for repeat architecture. Its
summary statistics include [1]:

| Feature | GRCh38 | T2T-CHM13 | Change |
|---|---:|---:|---:|
| Segmental duplications (Mbp) | 151.71 | 201.93 | +33.1% |
| Segmental duplications (% of assembly) | 5.00% | 6.61% | +1.61 pp |
| RepeatMasker repeats (% of assembly) | 51.89% | 53.94% | +2.05 pp |
| LINE (Mbp) | 626.33 | 631.64 | +0.8% |
| SINE (Mbp) | 386.48 | 390.27 | +1.0% |
| LTR (Mbp) | 267.52 | 269.91 | +0.9% |
| Satellite (Mbp) | 76.51 | 150.42 | +96.6% |
| Simple repeat (Mbp) | 36.5 | 77.69 | +112.9% |
| Low complexity (Mbp) | 6.16 | 6.44 | +4.6% |

More than half of the human genome is annotated repeat. The T2T assembly added
~74 Mbp of satellite and ~41 Mbp of simple repeat relative to GRCh38, exactly the
regions that were missing from earlier references and that are hardest for short-read
assemblies. A single-genome codec that does not model satellites, tandem arrays and
segmental duplications is leaving a large fraction of the genome on the table.

### 5.2 The repeat hierarchy

A useful taxonomy for compression, distinct from the biological taxonomy:

1. **Homopolymers and short tandem repeats (STRs).** Run-length and copy-number
   models. Cheap to encode, but a context model can waste bits on them.
2. **Tandem arrays and satellites.** Monomer + higher-order repeat structure +
   copy-number variation. Grammar with exponentiation rules, or a run-length grammar.
3. **Interspersed repeats (Alu, LINE, LTR, DNA transposons).** Family consensus +
   variants + insertion site. A typed grammar with family IDs, orientation and
   divergence.
4. **Segmental duplications.** Long paralogous blocks with divergence. A copy graph
   with match length, identity and orientation.
5. **Gene families and conserved domains.** Coding-sequence similarity; a context
   model with codon phase and a protein-level prior.
6. **Unique sequence.** The residual; a high-order context model with backoff.

The compression opportunity is to assign each region to the cheapest model, and to
encode the assignment compactly. This is a **latent region-type model**: a hidden Markov
model or change-point model over region classes, each with its own expert. The class
sequence is small (a few bits per change point), and the per-class models are much
better specified than one global model.

### 5.3 Sequence-space saturation and context order

The 2008 sequence-space analysis is still the clearest warning against high-order
k-mer tables [R]:

- More than 98% of all 12-mers appear in vertebrate genomes.
- Less than 2% of 19-mers appear.
- The number of distinct k-mers grows rapidly with k, while the number of observations
  per k-mer falls.

For a single human genome, a 16-mer table has 4.3 billion possible contexts and only
3.1 billion positions. Even a 12-mer table has 16.7 million contexts, many seen once.
The practical consequence is that a codec must **back off** from high-order contexts,
**mix** experts with different orders, and **use long matches** for the repeats that
context models cannot see. This is exactly what JARVIS3, GeCo3 and HyperGraSS do; the
next step is to make the long-match path explicit and rate-optimal.

### 5.4 Information content by feature

DNAGPT2's information-content map gives a model-relative view of where the bits are
[R]:

| Feature | Median information content (bits/base) |
|---|---:|
| Exons | 1.875 |
| Introns | 1.814 |
| Intergenic | 1.767 |
| Alu repeats | 0.631 |

These are not absolute entropies; they are the surprise under one model. But the
ordering is informative: Alu repeats are much more predictable than exons, and the
lowest-information bins concentrate in centromeric tandem repeats. A codec that models
repeats explicitly should beat a generic context model on exactly these regions.

### 5.5 What not to claim: 3D chromatin

It is tempting to say that chromatin loops, TADs or 3D genome architecture provide a
compression prior. They do not, unless the architecture is available to the decoder and
is not counted in the archive. A 3D contact map is a **phenotype**, not a free codebook:
if you store it, you pay for it; if you predict it from sequence, you have not added
information about the sequence. The same caution applies to annotations: gene models,
repeat masks and epigenetic tracks can help if they are available and cheap to encode,
but they are side information, not a free lunch.

The legitimate uses of annotation are:

- choosing a region type/model (and encoding the type sequence compactly);
- seeding a repeat-family grammar (and encoding the family ID);
- conditioning a context model on codon phase or CpG status (features derived from the
  sequence itself, so no side information).

Everything else must be paid for in the archive.

### 5.6 Genome-architecture implications

1. **Model repeats explicitly.** A context-only codec leaves satellite/SD/transposon
   redundancy on the table.
2. **Use a region-type mixture.** One global model is misspecified across isochores,
   repeats and genes.
3. **Exploit reverse complement.** DNA is double-stranded; a repeat and its reverse
   complement are the same opportunity. Canonicalization roughly halves the search
   space.
4. **Separate streams.** Match geometry, edit scripts, literals, region types and
   quality scores have different statistics. Mixing them in one byte stream wastes bits.
5. **Encode the model choice.** The parser's decisions (which region uses which model)
   are part of the archive. A good codec makes them cheap.

---

## 6. Where the cycles are: hardware architecture

### 6.1 The memory wall is the index

A genome compressor spends most of its time looking up contexts and matches. The index
footprint decides whether that lookup is a cache hit or a DRAM miss.

HyperGraSS's seed index is about **67 MB for a saturated 9.43 Mbp anchor** at
`SEED_STEP=4`, both strands, 70% load: **7.1 bytes per base** [M]. A human-scale index
at that density would be ~22 GB. That is larger than the working set of most CPUs and
far larger than a GPU cache. It is also why HyperGraSS's bounded mode struggles on
diverse pangenomes: the index, not the input, sets the memory bound.

The alternatives are:

- **Sampled indexes.** Syncmers/minimizers at density 1/8–1/16 reduce the index by an
  order of magnitude at a small ratio cost.
- **Compressed indexes.** SSHash, r-index, GBWT, Movi and prefix-free parsing store the
  index in space proportional to a repetitiveness measure rather than the text.
- **Blockwise indexes.** Index each block independently; accept a small ratio loss for
  bounded memory and random access.
- **External-memory indexes.** Keep the index on SSD/NVMe and use compressed-resident
  decode or near-memory processing.

The target for a human panel should be **1–2 bytes/base**, not 7. That is the difference
between a 3–6 GB index (fit on a GPU or a large server) and a 22 GB index (fit nowhere
convenient).

### 6.2 Parallelism and entropy coding

Adaptive arithmetic coding is serial: the probability state depends on every previous
symbol. There are three ways to parallelize it:

1. **Independent blocks.** Split the input into blocks, reset the model per block, and
   encode blocks in parallel. This is simple and fast, but the ratio loss grows as the
   block size shrinks. HyperGraSS's block experiments show the trade-off directly: on
   the 120 MB pangenome, 16 MiB blocks with `--seed-step 8` give 3,965,049 B, while
   32 MiB blocks give 3,824,166 B [M].
2. **Shared static model.** Train a model on a panel or a first pass, freeze it, and
   encode blocks in parallel with rANS. The model is stored once; the blocks are
   independent. This is the best fit for a reference-conditioned codec.
3. **Interleaved streams.** Encode `k` independent streams in one pass, interleaving
   their arithmetic coders. This hides the serial latency in SIMD/GPU pipelines. The
   streams can share a model if the model update is commutative or if the model is
   static [52][42].

rANS is the natural coder for options 2 and 3 [51]. It is table-driven, SIMD-friendly, and
parallelizable across streams. The cost is a small amount of ratio loss from the
quantized CDF and the independent streams. For DNA, the alphabet is four symbols, so
the tables are tiny and fit in L1.

### 6.3 GPU

The GPU is strongest at **match finding, BWT construction, and batch entropy coding**,
not at adaptive model updates. Reported results:

| Technique | Reported result | Relevance |
|---|---|---|
| GPU LZ77 decode (ACEAPEX) [38] | up to 260 GB/s on FASTQ; 165.7 GB/s on a 50 GB genome; 0.362 ms random read | copy-path decode; random access |
| DietGPU ANS on H100 | 592 GB/s decode | entropy-coding backend |
| GPU BWT via prefix-free parsing [39] | 11.9 bytes/char; 2,108 MiB input on a 24 GB RTX 3090 | index construction |
| panvram [44] | 558 haplotypes in 16.642 GB; ~740k windows/s | pangenome serving |
| ropebwt3 [43] | 320 human genomes in 65 h; 170 GB peak | pangenome index |

The GPU LZ77 result is particularly relevant: it shows that a copy-heavy format can
decode at memory-class speeds and support random access. The bottleneck moves to the
literal/residual channel, which is exactly where the context model lives. A hybrid that
copies 95% of the genome and models 5% can be both fast and compact.

### 6.4 FPGA and ASIC

FPGAs are attractive for **rANS/arithmetic coding** and for **streaming match
filtering** because they can sustain high throughput at low energy. Reported results:

- **Bancroft** [40] (Alveo U50, 8 GB HBM) provides the illusion of infinite on-device memory
  by compressing genomic data over PCIe; it reaches over 30% of HBM bandwidth, an order
  of magnitude above PCIe-limited architectures, and >6x over the conventional
  PCIe-attached architecture in a pre-alignment filtering application [R].
- **RAS** [41] integrates rANS with a probabilistic generator, stores distributions in BF16,
  and reports 121.2x encode and 70.9x decode speedups over a Python rANS baseline, with
  55% fewer decoder binary-search steps [R].
- **EXMA** [45] accelerates FM-index exact matching with a multi-task-learning index and
  CHAIN compression, improving search throughput 4.9x over state-of-the-art FM-index
  PIMs [R].

The caveat is development cost and portability. An FPGA rANS core is a good offload for
a production codec; it is not a research prototype's first priority. The CPU/GPU path
should come first.

### 6.5 Near-memory and compressed-resident decoding

The cheapest byte is the one you do not move. If the archive is served repeatedly
(e.g., a pangenome reference), keep it compressed in GPU/HBM and decode on demand.
panvram demonstrates the pattern: 558 human haplotypes in 16.6 GB, ~740k windows/s,
with random access [R]. ACEAPEX demonstrates position-invariant random access in a
GPU-resident LZ77 format, with a read-to-block index 6.3x smaller than a `.fai` [R].

For a single genome, the same idea applies at a smaller scale: store the archive in
blocks with checkpoints, decode only the block containing the query, and avoid
materializing the whole genome. This is already a HyperGraSS feature (page-indexed
random access); the next step is to make the block format compatible with GPU/FPGA
decode and with a reference panel.

### 6.6 Hardware-architecture implications

1. **Design for the copy path first.** Most bases should be copied, not modeled.
2. **Use rANS for the literal/residual channel.** It is parallel, SIMD-friendly and
   GPU/FPGA-friendly.
3. **Shrink the index.** A 7 bytes/base index is a scalability blocker; target 1–2.
4. **Block for memory, not just for parallelism.** Bounded memory and random access are
   product features, not afterthoughts.
5. **Measure end-to-end, not just the kernel.** A 592 GB/s ANS decoder is irrelevant if
   the model update is 10 MB/s.
6. **Keep a CPU fallback.** GPU/FPGA availability is not guaranteed; the archive format
   must decode on a CPU.

---

## 7. Candidate algorithms

Five candidates are worth considering. They are not mutually exclusive; the
recommended architecture combines RC-MDL, GAH and HNP.

### 7.1 RC-MDL: reference-conditioned MDL parser

**Idea.** Parse the target into phrases drawn from a reference panel, from the already
decoded prefix, from a grammar, or as literals. Price every phrase with the same
entropy model that codes it. Minimize the total code length, not the number of matches.

**Why it can win.** AGC and MBGC2 find matches and encode them with a separate literal
coder. RC-MDL chooses matches based on their actual bit cost, including the edit script
and the literal model. It also uses a stronger literal model (context mixing) and
reverse-complement canonicalization. The expected gain is 5–15% on divergent genomes,
larger on repeat-rich regions [P].

**Expected envelope.** 0.02–0.08 bits/base with a human panel; 1.4–1.6 bits/base in
self-contained mode; 10–100 MB/s per core on the copy path, slower on the literal path;
2–16 GB for a human panel, 256 MiB–4 GiB bounded blocks.

**Main risk.** The parser is a shortest-path problem with a context-dependent cost. The
cost depends on the parse, so the optimization is a fixed point. Beam search and
alternating minimization are approximations; the quality of the approximation decides
whether the gain is real.

### 7.2 GAH: grammar–attractor hybrid

**Idea.** Build a grammar with four rule types: concatenation, terminal, tandem
exponentiation, and reverse-complement copy. Score rules by MDL gain under the actual
entropy coder. Use a string-attractor approximation to choose anchors and bound the
grammar search. Encode the grammar with a context model.

**Why it can win.** Grammar captures nested and tandem repeats that a finite-context
model cannot see. GraSS is an early grammar codec but is weak on human-scale repeats;
a modern grammar with reverse-complement and exponentiation rules, MDL scoring and a
context-coded grammar stream has not been built for human DNA. The expected gain over
JARVIS3 is 3–10% on repeat-rich genomes [P].

**Expected envelope.** 1.30–1.60 bits/base on T2T; 1–10 MB/s; 4–32 GB; high
implementation risk.

**Main risk.** Grammar construction is expensive and easy to overfit. The grammar must
pay for its own description; a grammar that is large but not entropy-efficient loses to
a context model.

### 7.3 NCM: neural copy model

**Idea.** A small transformer or state-space model predicts the next base from a long
context, with a retrieval memory over a reference panel. The model is quantized and
shared; arithmetic coding uses its probabilities. The copy mechanism is explicit: the
model can attend to a matching reference segment and predict the next base from it.

**Why it might win.** Neural models learn motifs and feature-specific predictability
that hand-engineered models miss. DNAGPT2 is only 6% behind JARVIS3 despite no explicit
repeat module; a hybrid with a copy mechanism could close the gap and surpass it [P].

**Expected envelope.** 1.30–1.50 bits/base self-contained; 0.05–0.20 bits/base with
retrieval; 0.1–1 MB/s on CPU, 10–100 MB/s on GPU; 8–80 GB memory.

**Main risk.** Model side information. A 170 MB model costs ~0.45 bits/base on one human
genome. NCM is only viable if the model is small (<10 MB), heavily quantized, or
amortized across a collection. It is a research bet, not the first production track.

### 7.4 PODA: phylogeny-ordered delta archive

**Idea.** Build a phylogeny of the collection, order genomes by a DFS traversal, and
delta-code each genome against its parent/nearest neighbor with RC-MDL. Encode the tree
and the deltas.

**Why it can win.** Phylogenetic compression [31][32] reports 27.8x (xz) to 38.9x
(xz + MBGC) on 661k bacterial assemblies, and 18.1 bytes/genome on 591k SARS-CoV-2
genomes [R]. AGC and MBGC already exploit collection redundancy; a phylogeny-aware
order plus a stronger delta coder can improve the ratio further.

**Expected envelope.** 0.005–0.05 bits/base on collections; 10–100 MB/s; 4–32 GB.

**Main risk.** This is Mode C, not single-genome. It is the natural extension for
pangenomes and repositories, but it does not answer the self-contained single-genome
question.

### 7.5 HNP: hardware-native rANS pipeline

**Idea.** A block-parallel rANS backend with SIMD/GPU/FPGA acceleration, static or
slowly adapting CDFs, and a compressed-resident block format with checkpoints.

**Why it can win.** It makes the literal/residual channel fast enough that the hybrid
is practical. GPU LZ77 decode at 165–260 GB/s and DietGPU ANS at 592 GB/s show the
hardware can sustain memory-class throughput [R]. The bottleneck is the model, not the
coder.

**Expected envelope.** Same ratio as the underlying model; 5–10x speedup at <1% ratio
loss; 1–8 GB memory.

**Main risk.** It is an enabler, not a ratio breakthrough. Do not expect HNP alone to
beat JARVIS3.

### 7.6 Trade-off summary

| Candidate | Mode | Expected bits/base | Speed | Memory | Risk | Recommendation |
|---|---|---|---:|---:|---|---|
| RC-MDL | R/S | 0.02–0.08 / 1.4–1.6 | 10–100 MB/s | 2–16 GB | medium | **primary** |
| GAH | S | 1.30–1.60 | 1–10 MB/s | 4–32 GB | high | **self-contained track** |
| NCM | S/R | 1.30–1.50 / 0.05–0.20 | 0.1–100 MB/s | 8–80 GB | high | research bet |
| PODA | C | 0.005–0.05 | 10–100 MB/s | 4–32 GB | medium | collection extension |
| HNP | backend | same as model | 5–10x | 1–8 GB | low–medium | **enabler** |

The projections are deliberately wide. They are hypotheses to be tested, not promises.

---

## 8. Recommended architecture: MDL-Mosaic

The recommended design is **MDL-Mosaic**, with RC-MDL as its reference-conditioned
instantiation and GAH as its self-contained long-repeat component. The name is
descriptive: the genome is a mosaic of copied, repeated and novel segments, and the
codec chooses the mosaic that minimizes the description length.

### 8.1 Formal model

Let `x ∈ {A,C,G,T}^n`. Let `R = {R_1, ..., R_m}` be an optional reference panel; `R`
may be empty. A parse

```text
P = ((t_1, p_1), ..., (t_k, p_k))
```

covers `x`, where each phrase type `t_i` is one of:

- `REF(r, a, l, o)`: copy length `l` from reference `r` at `a`, orientation `o`.
- `SELF(a, l, o)`: copy from the already decoded prefix of `x`.
- `GRAM(A)`: expand grammar rule `A`.
- `LIT(s)`: literal string `s`.
- `APPROX(r, a, l, o, E)`: approximate copy with edit script `E`.

The code length is

```text
L(x) = min over P, θ of
       [ L(θ) + Σ_i -log2 p(t_i, p_i | c_i, θ)
              + Σ_{i: LIT} -log2 p(s_i | c_i, θ) ].
```

`c_i` is the context available to the decoder at phrase `i`: the previous phrases, the
reference sequence (if any), the target prefix, and the current region type. `θ` is the
model: context experts, copy-length/distance distributions, edit distributions, grammar
probabilities and region-type transitions.

This is a two-part MDL code. The first term pays for the model; the second pays for the
data given the model. The parser and the entropy coder use the same `p`, so there is no
mismatch between "match cost" and "actual bits." For an adaptive model, `θ` is not
stored: it is a deterministic function of the decoded prefix, so `L(θ)` is the cost of
the initial state and quantization. For a static model, `L(θ)` is the table cost. Because
the coder uses the same `p`, the realized code length is at most the parser's cost plus
the arithmetic-coder redundancy (at most 2 bits per block for a binary coder).

### 8.2 The parser

The parser is a shortest-path search in a DAG:

- **Nodes:** positions `0..n` in `x`.
- **Edges:** candidate phrases starting at a position.
- **Edge cost:** `-log2 p(t, p | c, θ)` plus the literal cost.
- **Goal:** a path from `0` to `n` with minimum total cost.

Candidate generation:

1. Build a syncmer/minimizer index of the reference panel (and of the decoded prefix
   for self mode).
2. For each position, look up seeds and extend exact matches. Keep the longest matches
   per reference and orientation.
3. For approximate matches, use a seed-and-extend aligner with a banded edit distance;
   keep candidates whose estimated edit cost is below the literal cost.
4. For grammar phrases, look up rules whose yield matches the current position (using a
   grammar index or a hash of rule prefixes).
5. For literals, the edge is a single base or a short run; the context model prices it.

Search:

- A full dynamic program over all candidates is too expensive.
- Use a beam search with width `B` (e.g., 64–256) and a lookahead of `W` bases (e.g.,
  1–4 kb).
- Use A* with an admissible heuristic based on the minimum bits per base; the heuristic
  can be precomputed from the context model's average rate.
- Re-parse after model updates (alternating minimization): 2–5 iterations.

The parser is where the ratio gain lives. A greedy longest-match parser is fast but
chooses matches that may be expensive to encode; a cost-aware parser chooses fewer,
better matches.

### 8.3 The entropy model

The model is a mixture of experts:

```text
p(x_i | c) = Σ_k w_k(c) p_k(x_i | c)
```

with experts:

- **Order-k context** for `k = 1..24`, with Krichevsky–Trofimov or KT-style smoothing.
- **Reverse-complement context**: the same order-k model on the reverse complement of
  the context; useful for double-stranded symmetry.
- **Substitution-tolerant context**: a context pattern that allows mismatches, as in
  GeCo3/JARVIS3.
- **Reference-predicted base**: if the position is aligned to a reference, the reference
  base and its flanking k-mers are features.
- **Copy expert**: the match length, distance, orientation and local identity.
- **Run-length/homopolymer expert**: for low-complexity regions.
- **Region-type expert**: a latent class (repeat family, exon, intron, intergenic) with
  its own context model.
- **Neural expert (optional)**: a small quantized model, used only if it pays for
  itself.

The mixing weights `w_k(c)` are a logistic function of features: context order, match
length, local GC, region type, recent error rate. They are updated online by gradient
descent or Bayesian mixing. The model is integer-only at encode/decode time to avoid
cross-platform drift; training can use floating point offline.

### 8.4 The entropy coder

Two coder modes:

- **Adaptive range coder** for self-contained mode and for the final residual. It
  follows the model exactly; it is serial but simple.
- **Static rANS** for reference-conditioned mode and for block-parallel decode. The CDF
  is trained on a first pass or on the panel, quantized to 12–16 bits, and stored once
  per block. Blocks are independent; decode is SIMD/GPU-friendly.

The archive records the coder mode, the CDF tables (or a hash), and the block
checkpoints. The decoder must reproduce the exact integer probabilities; any floating
point in the coder is a bug.

### 8.5 The archive format

A block-structured format:

```text
header:
  magic, version, mode (S/R/C), alphabet, input length, input SHA-256
  reference_required, reference_digest, reference_lengths
  block_size, coder_mode, model_hash
per block:
  block_offset, block_length, raw_length, crc32
  phrase stream:
    type, params, payload
  checkpoint:
    model_state_hash, output_offset
footer:
  block index, total_length, archive_sha256
```

The block size is a tunable trade-off: 16–64 MiB is a good starting point. Each block
is independently decodable given the model state at its checkpoint; the checkpoint is
the bridge between random access and compression ratio.

### 8.6 Self-contained fallback

If no reference panel is available:

- Set `R = ∅`.
- Allow `SELF` phrases from the decoded prefix (LZ77 with a large window).
- Enable `GRAM` phrases for tandem and nested repeats.
- Use the adaptive range coder for the residual.

The same parser becomes a fuzzy LZ77 + grammar + context codec. This is the GAH track.
It will be slower than the reference mode because more bases take the literal path, but
it is unconditionally lossless and does not depend on an external panel.

### 8.7 Complexity

Let `n` be the input length, `h` the average number of seed hits, `B` the beam width,
and `W` the lookahead window.

- Candidate generation: `O(n·h)` time, `O(index)` memory.
- Beam search: `O(n·B·W)` time in the worst case; in practice `W` is small and most
  positions have few candidates.
- Model update: `O(n)` time, `O(model)` memory.
- Entropy coding: `O(n)` time.

With a 16–64 MiB block, the working set is bounded. The index is the dominant memory
term; a sampled/compressed index is essential for human-scale panels.

### 8.8 Why this is mathematically better

1. **It optimizes the actual code length.** The parser's edge costs are the coder's
   probabilities. There is no heuristic "match score" that diverges from the bits.
2. **It combines the two redundancy classes.** Short-range context and long-range
   repeats are modeled by different experts, and the mixture chooses between them per
   position.
3. **It makes the prior explicit.** A reference panel is a model, not a magic oracle.
   The archive records the panel digest and pays for the differences.
4. **It amortizes the expensive model.** Most bases are copied; the context model runs
   on the residual. This is the only way to get both a low rate and a practical speed.
5. **It is falsifiable.** Each component has an ablation and a go/no-go gate (§9).

### 8.9 Lossless verification

The format must be lossless by construction, not by hope:

- **Reference digest.** The archive stores the SHA-256 of every reference sequence and
  the panel version. A mismatch is a hard error, not a silent fallback.
- **Match validation.** Every exact match is verified byte-for-byte at encode time; the
  decoder also checks bounds. Approximate matches carry an explicit edit script; the
  decoder applies it and verifies the resulting phrase.
- **Per-block CRC-32 and full SHA-256.** CRC catches corruption; SHA-256 catches
  logical errors.
- **Independent verifier.** A separate process reconstructs the input and compares
  SHA-256 and bytes, as in `benchmarks/verify/verify.py`.
- **Fuzzing and property tests.** Random sequences, all-IUPAC, all-N, empty records,
  single-base records, huge records, lowercase/uppercase, CRLF/LF, RNA `U`, chimeric
  DNA/RNA, and reference permutations.
- **Deterministic encoder.** The encoder must be deterministic for a given input,
  model and parameters; ties in the parser are broken by a fixed rule.
- **Cross-implementation tests.** At least one independent decoder (e.g., Python) for
  the format, even if slow, to catch spec ambiguities.
- **Model hash.** For neural or learned models, store the model hash, quantization
  scheme and kernel version. Floating-point non-determinism is not acceptable in the
  decoder.

The verification contract for Mode R is explicit: "lossless given the pinned panel."
For Mode S, the contract is unconditional. The archive header states which contract
applies.

---

## 9. Experimental plan

### 9.1 Questions and hypotheses

| ID | Question | Hypothesis | Falsification |
|---|---|---|---|
| Q1 | How low can a self-contained grammar+context codec go on T2T? | 1.30–1.50 bits/base | If GAH is not ≥3% below JARVIS3, drop the grammar track |
| Q2 | How much does a human panel buy over JARVIS3? | 10–30x | If RC-MDL is not ≥5x below JARVIS3 at same input, the panel is not paying |
| Q3 | Does cost-aware parsing beat AGC/MBGC2? | 5–15% fewer bytes | If not ≥5% at same panel/speed, the parser is not worth it |
| Q4 | What is the speed/ratio Pareto of block-parallel rANS? | 5–10x speed at <1% ratio loss | If <3x at >1% loss, the backend is not worth it |
| Q5 | Does a neural expert beat explicit repeat models on repeats? | 5–10% on repeat-rich regions | If not ≥3%, keep neural as analysis only |
| Q6 | How small can the panel be? | 10–50 references for most of the gain | If 100+ references are needed for 5%, the deployment cost is too high |

### 9.2 Datasets

- **T2T-CHM13v2.0** (3,117,292,120 bases) [1]: the human self-contained benchmark.
- **100 Mb chr1** (the project's existing input, SHA-256 `0c1a...`): continuity with
  HyperGraSS.
- **HPRC 96/558 haplotypes** [2]: reference-conditioned human panel.
- **24-strain E. coli pangenome** (120,108,223 bases): the project's pangenome anchor.
- **ArTh / Arabidopsis thaliana**: cross-species plant benchmark; compare to the
  AgentGC/JARVIS3 table.
- **Random sequence**: negative control; a good codec should not expand it much.
- **Satellite/low-complexity synthetic**: stress the repeat model.
- **IUPAC/N/low-complexity edge cases**: lossless verification.
- **Ancient DNA / metagenome** (optional): divergence and contamination stress.

Preprocessing must be identical across codecs: uppercase, non-ACGT handling, line
wrapping, headers. Report the exact input hash.

### 9.3 Baselines

Self-contained: JARVIS3, GeCo3, GraSS, NAF, FFC, xz -9, zstd -19, HyperGraSS.

Reference/collection: AGC, MBGC2, RLZ/RLZAP, ERGC, iDoComp, phylogenetic compression.

Neural: DNAGPT2, PMKLC, AgentGC, DeepGeCo (if runnable).

Every baseline must report archive bytes, bits/base, wall time, CPU time, peak RSS,
thread count, and exactness. If a baseline cannot run under the same memory cap, record
the asymmetry.

### 9.4 Metrics

- **Ratio:** archive bytes, bits/base, compression factor.
- **Speed:** compress and decompress MB/s; random-access latency.
- **Memory:** peak RSS, GPU memory, index size.
- **Energy:** joules per GB compressed/decompressed (optional but valuable).
- **Exactness:** SHA-256 and byte comparison, independent verifier.
- **Robustness:** edge cases, malformed input, reference mismatch.
- **Pareto:** ratio vs speed vs memory; no single-number winner.

### 9.5 Phases and gates

| Phase | Work | Gate |
|---|---|---|
| 0 | Freeze the HyperGraSS baseline and the exact input set | Reproduce the 100 Mb and 120 Mb numbers |
| 1 | Build the verification/benchmark harness for the new format | Independent verifier passes on all edge cases |
| 2 | RC-MDL prototype on 10–100 Mb with a small panel | ≥5% below AGC/MBGC2 at same panel, or stop |
| 3 | GAH prototype on 10–100 Mb chr1 | ≥3% below JARVIS3, or stop |
| 4 | Full human T2T run for the winning track | Byte-exact; ratio and speed recorded |
| 5 | HNP backend (SIMD/rANS, optional GPU) | ≥5x speedup at <1% ratio loss |
| 6 | Integration, docs, release candidate | All gates pass; no regression on pangenome |

Phase 2 is the decisive experiment. If a cost-aware parser cannot beat AGC/MBGC2 on a
small panel, the central thesis is wrong and the project should pivot to GAH or to a
pure collection codec.

### 9.6 Ablations

- Exact matches only vs exact + approximate.
- Cost-aware parser vs greedy longest-match.
- Context model vs fixed literal costs.
- Reverse-complement canonicalization on/off.
- Region-type mixture on/off.
- Grammar rules on/off (self mode).
- Panel size: 1, 5, 10, 50, 100, 558.
- Block size: 1, 16, 64, 256 MiB.
- Adaptive range coder vs static rANS.
- SIMD/GPU on/off.

Each ablation reports bytes, bits/base, speed and memory. The point is attribution, not
just a headline number.

### 9.7 Statistical rigor

- Use at least three interleaved repeats for speed/RSS; report the worst case, not the
  best.
- Report archive bytes exactly; they are deterministic.
- Separate host-noise-sensitive metrics (wall time, RSS) from stable metrics (bytes,
  hashes).
- Do not compare bits/base across different inputs without saying so.
- Keep the input manifests and hashes.

---

## 10. Failure modes and mitigations

| Failure mode | Consequence | Mitigation |
|---|---|---|
| Reference mismatch | Wrong or corrupt decode | Store panel SHA-256; hard error on mismatch |
| Approximate-match false positive | Corrupt output | Verify every phrase byte-for-byte; never trust a hash-only match |
| Parser/model mismatch | Ratio worse than expected | Price phrases with the actual coder; ablate |
| Grammar blowup | Archive larger than input | Fallback to literal/context; cap grammar size |
| Model side information | Neural codec loses to JARVIS3 | Count weights; quantize; amortize across collection |
| Block independence | Ratio loss | Larger blocks; shared static model; checkpoints |
| Index memory | OOM on human panel | Sampled/compressed index; blockwise; external memory |
| GPU/FPGA unavailable | No decode | CPU fallback; portable format |
| Non-deterministic float | Decode drift | Integer-only coder; fixed kernels; model hash |
| IUPAC/N/lowercase | Expansion or corruption | Separate alphabet streams; explicit edge-case tests |
| Low-complexity runs | Context model wastes bits | Run-length/homopolymer expert |
| Contamination/adapters | Benchmark distortion | Preprocess identically; report input hash |
| Overfitting to T2T | False generalization | Cross-species panel; ArTh; random controls |

The two failure modes that kill projects are **reference mismatch** and **unverified
approximate matches**. Both are solved by explicit contracts and byte-level validation,
not by trusting the aligner.

---

## 11. Open problems

1. **Joint parse/model optimization.** The parser's cost depends on the model, and the
   model's parameters depend on the parse. Alternating minimization is a heuristic; a
   convex or game-theoretic formulation is open.
2. **Optimal reference panel selection.** Given a budget of `k` references, which ones
   maximize compression of an unseen target? This is a submodular optimization problem
   with a non-additive objective.
3. **Streaming grammar construction.** RePair on 3 Gb is memory-heavy. A bounded-memory,
   blockwise grammar with a global dictionary is an open engineering problem.
4. **Reverse-complement-aware attractors.** The string-attractor theory is oriented to
   forward strings; a double-stranded attractor (closed under reverse complement) is
   not fully developed.
5. **Neural priors that pay for themselves.** A model that is small enough to amortize
   on one genome and strong enough to beat JARVIS3 is not yet demonstrated.
6. **Random access vs ratio.** Checkpoints cost bits; the optimal checkpoint density
   for a given query distribution is an open trade-off.
7. **Benchmark standardization.** Bits/base is not enough; the field needs a standard
   protocol that includes model side information, preprocessing, thread count and
   exactness. The Cobilab HumanGenome benchmark is a good start.
8. **Privacy and security.** Reference panels can leak; compression side channels can
   reveal membership. A production codec needs a threat model.
9. **Hardware portability.** A GPU/FPGA codec that cannot decode on a CPU is a niche
   product. Bit-exact cross-platform decoding is a requirement.
10. **The theoretical floor.** The Kolmogorov complexity of a human genome is
    uncomputable; the best practical floor is model-dependent. A principled lower bound
    for the grammar+context class would tell us when to stop.

---

## 12. References

1. Nurk S. et al. The complete sequence of a human genome. *Science* 376, 44–53 (2022).
   doi:10.1126/science.abj6987.
2. Liao W.-W. et al. A draft human pangenome reference. *Nature* 617, 312–324 (2023).
   doi:10.1038/s41586-023-05896-x.
3. Sousa A., Pratas D., Pinho A. J. JARVIS3: an efficient encoder for genomic data.
   *Bioinformatics* 40, btae725 (2024). doi:10.1093/bioinformatics/btae725.
4. Silva M., Pratas D., Pinho A. J. Efficient DNA sequence compression with neural
   networks. *GigaScience* 9, giaa119 (2020). doi:10.1093/gigascience/giaa119.
5. Roy A., Maity S., Mukhopadhyay A. A lossless reference-free sequence compression
   algorithm leveraging grammatical, statistical, and substitution rules.
   *Briefings in Functional Genomics* 24, elae050 (2025). doi:10.1093/bfgp/elae050.
6. Kryukov K. et al. Nucleotide Archival Format (NAF) enables efficient lossless
   reference-free compression of DNA sequences. *Bioinformatics* 35, 3826–3828 (2019).
   doi:10.1093/bioinformatics/btz144.
7. Grabowski S., Kowalski T. M. MBGC: Multiple Bacteria Genome Compressor.
   *GigaScience* 11, giab099 (2022). doi:10.1093/gigascience/giab099.
8. Kowalski T. M. MBGC2: Boosting compression via efficient encoding of approximate
   matches in genome collections. *GigaScience* 15, giag008 (2026).
   doi:10.1093/gigascience/giag008.
9. Deorowicz S., Danek A., Niemiec M. AGC: compact representation of assembled
   genomes with fast queries and updates. *Bioinformatics* 39, btad097 (2023).
   doi:10.1093/bioinformatics/btad097.
10. FFC: a scalable FASTA compressor. *Bioinformatics* 42, btag132 (2026).
    doi:10.1093/bioinformatics/btag132.
11. Kuruppu S., Puglisi S. J., Zobel J. Relative Lempel-Ziv compression of genomes for
    large-scale storage and retrieval. *SPIRE 2010*, LNCS 6393, 201–212.
    doi:10.1007/978-3-642-16321-0_20.
12. Kuruppu S., Puglisi S. J., Zobel J. Relative Lempel-Ziv with constant-time random
    access. *SPIRE 2014*, LNCS 8799. doi:10.1007/978-3-319-11918-2_2. See also RLZAP:
    doi:10.1007/978-3-319-46049-9_1.
13. Ochoa I., Hernaez M., Weissman T. iDoComp: a compression scheme for assembled
    genomes. *Bioinformatics* 31, 626–633 (2015). doi:10.1093/bioinformatics/btu698.
14. Saha S., Rajasekaran S. ERGC: an efficient referential genome compression
    algorithm. *Bioinformatics* 31, 3468–3475 (2015). doi:10.1093/bioinformatics/btv399.
15. Pinho A. J., Pratas D., Garcia S. P. GReEn: a tool for efficient compression of
    genome resequencing data. *Nucleic Acids Research* 40, e27 (2012).
    doi:10.1093/nar/gkr1124.
16. Zhu Z. et al. DNA-COMPACT: DNA COMpression based on a pattern-aware contextual
    modeling technique. *PLoS ONE* 8, e80377 (2013).
    doi:10.1371/journal.pone.0080377.
17. Yao J. et al. HRCM: an efficient hybrid referential compression method for genomic
    big data. *BioMed Research International* 2019, 3108950.
    doi:10.1155/2019/3108950.
18. Saha S., Rajasekaran S. NRGC: a novel referential genome compression algorithm.
    *Bioinformatics* 32, 3405–3412 (2016). doi:10.1093/bioinformatics/btw505.
19. Wang S. et al. SparkGC: Spark based genome compression for large collections of
    genomes. *BMC Bioinformatics* 23, 297 (2022).
    doi:10.1186/s12859-022-04825-5.
20. Reference-based genome compression using the longest matched substrings with
    parallelization consideration. *BMC Bioinformatics* 24, 371 (2023).
    doi:10.1186/s12859-023-05500-z.
21. Deorowicz S., Danek A., Grabowski S. GDC 2: compression of large collections of
    genomes. *Scientific Reports* 5, 11565 (2015). doi:10.1038/srep11565.
22. Christley S., Lu Y., Li C., Xie X. Human genomes as email attachments.
    *Bioinformatics* 25, 274–275 (2009). doi:10.1093/bioinformatics/btn582.
23. Cox A. J., Bauer M. J., Jakobi T., Rosone G. Large-scale compression of genomic
    sequence databases with the Burrows–Wheeler transform. *Bioinformatics* 28,
    1415–1419 (2012). doi:10.1093/bioinformatics/bts173.
24. Prezza N. String attractors. arXiv:1709.05314 (2017).
25. Kempa D., Prezza N. At the roots of dictionary compression: string attractors.
    *STOC 2018*, 827–840. doi:10.1145/3188745.3188814.
26. Gagie T., Navarro G., Prezza N. Optimal-time text indexing in BWT-runs bounded
    space. *SODA 2018*, 1459–1477. doi:10.1137/1.9781611975031.96.
27. Boucher C. et al. Prefix-free parsing for building big BWTs. *Algorithms for
    Molecular Biology* 14, 13 (2019). doi:10.1186/s13015-019-0148-5.
28. Sirén J., Garrison E., Novak A. M., Paten B., Durbin R. Haplotype-aware graph
    indexes. *Bioinformatics* 36, 400–407 (2020). doi:10.1093/bioinformatics/btz575.
29. Garrison E. et al. Variation graph toolkit improves read mapping by representing
    genetic variation in the reference. *Nature Biotechnology* 36, 875–879 (2018).
    doi:10.1038/nbt.4227.
30. Gagie T., Manzini G., Sirén J. Wheeler graphs: a framework for BWT-based data
    structures. *Theoretical Computer Science* 698, 67–78 (2017).
    doi:10.1016/j.tcs.2017.06.016.
31. Břinda K. et al. Efficient and robust search of microbial genomes via phylogenetic
    compression. *Nature Methods* 22, 1231–1241 (2025).
    doi:10.1038/s41592-025-02625-2.
32. Břinda K. et al. Efficient and robust search of microbial genomes via phylogenetic
    compression. bioRxiv 2023.04.15.536996. doi:10.1101/2023.04.15.536996.
33. DeepGeCo: Genomics data lossless compression with (s,k)-mer encoding and deep
    neural networks. *AAAI 2025*. doi:10.1609/aaai.v39i12.33371.
34. PMKLC: Parallel multi-knowledge learning-based lossless compression for large-scale
    genomics database. *CIKM 2025*. doi:10.1145/3711896.3737083. arXiv:2507.12805.
35. AgentGC: Evolutionary learning-based lossless compression for genomics data with
    LLM-driven multiple agent. arXiv:2601.13559 (2026).
36. Macala V., Simecek P. DNA compression with genomic language models: tokenization,
    benchmarking, and an information-content map. bioRxiv 2026.06.10.731316.
    doi:10.64898/2026.06.10.731316.
37. Erlich Y., Zielinski D. DNA Fountain enables a robust and efficient storage
    architecture. *Science* 355, 950–954 (2017). doi:10.1126/science.aaj2038.
38. Compressed-resident genomics: full-pipeline device-resident GPU LZ77 decode with
    position-invariant random access. arXiv:2606.18900 (2026).
39. Engineering lightweight Burrows-Wheeler transform construction on the GPU via
    prefix-free parsing. arXiv:2610.33889 (2026).
40. Bancroft: genomics acceleration beyond on-device memory. arXiv:2502.16470 (2025).
41. RAS: a bit-exact rANS accelerator for high-performance neural lossless compression.
    arXiv:2511.04684 (2025).
42. Recoil: parallel rANS decoding with decoder-adaptive scalability. arXiv:2306.12141
    (2023).
43. Li H. BWT construction and search at the terabase scale. *Bioinformatics* 40,
    btae717 (2024). doi:10.1093/bioinformatics/btae717. arXiv:2409.00613.
44. panvram: GPU-resident random access to a human pangenome cohort. arXiv:2610.25149
    (2026).
45. EXMA: a genomics accelerator for exact-matching. arXiv:2101.05314 (2021).
46. Sequence space coverage, entropy of genomes and the potential to detect non-human
    DNA in human samples. *BMC Genomics* 9, 509 (2008).
    doi:10.1186/1471-2164-9-509.
47. Vertical lossless genomic data compression tools for assembled genomes: a
    systematic literature review. *PLoS ONE* 15, e0232942 (2020).
    doi:10.1371/journal.pone.0232942.
48. A survey on gene sequence compression algorithms based on reference sequences.
    *Frontiers of Data and Computing* 6, 59–76 (2024).
    doi:10.11871/jfdc.issn.2096-742X.2024.04.005.
49. Pinho A. J., Pratas D. MFCompress: a compression tool for FASTA and multi-FASTA
    data. *Bioinformatics* 30, 117–118 (2014). doi:10.1093/bioinformatics/btt594.
50. Witten I. H., Neal R. M., Cleary J. G. Arithmetic coding for data compression.
    *Communications of the ACM* 30, 520–540 (1987). doi:10.1145/214762.214771.
51. Duda J. Asymmetric numeral systems: entropy coding combining speed of Huffman
    coding with compression rate of arithmetic coding. arXiv:1311.2540 (2013).
52. Massively-parallel lossless data decompression. arXiv:1606.00519 (2016).
53. Movi: a fast and cache-efficient full-text pangenome index. *iScience* 27, 111464
    (2024). doi:10.1016/j.isci.2024.111464.
54. Li H. The design and construction of reference pangenome graphs with minigraph.
    *Genome Biology* 21, 265 (2020). doi:10.1186/s13059-020-02168-z.
55. Guarracino A. et al. ODGI: understanding pangenome graphs. *Bioinformatics* 38,
    3319–3323 (2022). doi:10.1093/bioinformatics/btac308.
56. Hickey G. et al. Pangenome graph construction from genome alignments with
    Minigraph-Cactus. *Nature Biotechnology* 42, 663–673 (2024).
    doi:10.1038/s41587-023-01793-w.

---

## Appendix A: measured HyperGraSS baseline (for continuity)

The full reports are in `benchmarks/specialized/runs/`. The rows below are the stable
measurements; wall time and RSS are host-sensitive.

| Run | Archive (B) | bits/base | Compress (s) | Peak RSS (MiB) | Restore (s) | Exact |
|---|---:|---:|---:|---:|---:|:---:|
| 100 Mb `--single` | 20,626,877 | 1.650 | 355.0 | 396 | — | yes |
| 100 Mb `--threads 2` | 20,635,964 | 1.651 | 70.1 | 586 | — | yes |
| 100 Mb `--threads 4` | 20,649,693 | 1.652 | 124.9 | 617 | — | yes |
| 100 Mb `--threads 8` | 20,677,028 | 1.654 | 225.9 | 697 | — | yes |
| 100 Mb default | 21,056,843 | 1.685 | 152.9 | 967 | — | yes |
| 100 Mb `--lean` | 22,788,010 | 1.823 | 40.5 | 94 | — | yes |
| 100 Mb `--pangenome` | 22,735,705 | 1.819 | 744.2 | 962 | — | yes |
| 100 Mb GraSS | 22,264,066 | 1.781 | 137.5 | 893 | — | yes |
| 100 Mb `--bclzg` | 21,120,171 | 1.690 | 181.1 | 991 | 167.0 | yes |
| 100 Mb `--achmc` | 21,151,440 | 1.692 | 632.6 | 1,438 | 568.9 | yes |
| 120 Mb tuned | 3,810,409 | 0.254 | 887.3 | 999 | — | yes |
| 120 Mb JARVIS3 l7 | 3,896,361 | 0.260 | 95.4 | 286 | — | yes |
| 120 Mb native bounded | 3,864,464 | 0.257 | 54.8 | 240.1 | 25.8 | yes |

The 120 Mb native bounded row passes the strict 256 MiB gate: 240.1 MiB compress,
225.5 MiB restore, size ratio 1.014, wall ratio 0.940. It is the only bounded
configuration that passes all four gate checks.

## Appendix B: glossary

- **bits/base (bpb):** archive bits per input nucleotide. 2.0 is no compression; 1.0 is
  8x; 0.08 is 100x.
- **MDL:** minimum description length; the two-part code `model + data`.
- **Parse:** a segmentation of the target into copied, repeated and literal phrases.
- **Relative LZ:** LZ77 where the dictionary is a reference, not the sequence itself.
- **String attractor:** a set of positions that covers all distinct substrings.
- **r-index:** a BWT-run-based index supporting pattern matching in `O(r)` space.
- **GBWT:** graph Burrows-Wheeler transform; a graph index for haplotype path sets.
- **rANS:** range asymmetric numeral systems; a table-driven entropy coder.
- **Mode S/R/C:** self-contained / reference-conditioned / collection.

---

## Appendix C: one-page decision

**Decision.** Do not build another self-contained context-model-only codec. Build
MDL-Mosaic with RC-MDL as the primary track and GAH as the self-contained fallback.

**Why.** The evidence says the ratio is in the prior, not the entropy coder. Self-
contained methods are clustered at 1.38–1.9 bits/base; reference methods reach
0.006–0.1 bits/base. A 5% better context model is useful but not transformative. A
cost-aware parser plus a panel is the only path to 100:1.

**What to build first.** A 10–100 Mb RC-MDL prototype on a small panel, with a byte-
exact verifier, compared against AGC and MBGC2 at the same panel and input. If it does
not win by ≥5%, stop and pivot to GAH.

**What not to do.** Do not claim self-contained 100:1. Do not compare bits/base across
different inputs. Do not count a neural model's bits without its weights. Do not trust
an approximate match without byte-level verification. Do not resume the HyperGraSS
profile sweep; it answered its question.

**What would change the decision.** A self-contained grammar+context codec that beats
JARVIS3 by ≥10% on T2T, or a neural model that beats JARVIS3 with <10 MB of side
information, would reopen the self-contained track. Neither exists today.
