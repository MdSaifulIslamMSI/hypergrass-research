# HyperGraSS: Lossless Genomic and Pangenomic Compression with a Stride-Sampled 16-mer Seed Index and Context Mixing

**Working paper draft — not peer reviewed.** This document describes the
implementation and the measurement protocol. Numeric results are in
[`BENCHMARKS.md`](BENCHMARKS.md) and are generated from the reports in this
repository; they are not repeated as independent claims here.

## Abstract

HyperGraSS is a pure-Java 27 lossless compressor for FASTA, multi-FASTA, raw
DNA, and RNA sequence data. It factorizes input into six orthogonal channels,
encodes the DNA stream with a stride-sampled 16-mer referential model and an
integer-only context-mixing range coder, and stores a page-indexed footer that
supports record and subsequence queries without decompressing the whole
archive. A bounded-memory streaming mode is intended for inputs larger than the
Java heap. We evaluate the implementation on a frozen held-out corpus that was
stopped after 179 of 280 planned results (6 of 10 datasets complete) and on a
121.8 MB real multi-strain *E. coli* pangenome. On the completed sequence-only
cells the default archive is smallest by a median margin of 1.36%; on the
121.8 MB complete-FASTA pangenome, the default profile is 9.8% larger than
xz -9, while an opt-in 16 MiB bank / 8 MiB block profile is 41.8% smaller and
answers exact record queries. The tuned profile replicated on two additional
~120 Mbp RefSeq pangenomes, *Salmonella enterica* (43.7% smaller than xz) and
*Klebsiella pneumoniae* (46.2% smaller). A measured ablation attributes that
result to the default's 4 MiB anchor-bank cap. An LZ77 plus range-coder gap
channel then changes the specialized comparison: on the identical 120.1 Mbp
sequence the tuned profile reaches 3,810,409 bytes, 2.21% below JARVIS3 level 7
(3,896,361 bytes) at 9.3x its compression time, so the margin against a
specialized DNA compressor is real but narrow. On six additional natural
collections spanning five phyla, the tuned profile was smaller than JARVIS3
level 7 on all six, by a median of 13.8% (range 4.7–29.4%), while the published
GraSS pipeline was 7.3–11.7x larger. The JVM 256 MiB bounded variant does not generalize: it stays under the RSS
cap on one of the six collections, reaches 271–273 MiB on four, and exhausts
its heap on *P. aeruginosa*; its archives are 6.6–19.8% larger than the
unbounded profile. An opt-in GraalVM native image of the same sources, with
split-tail blocking, completes 3/3 under a strict 256 MiB no-swap cap on the
*E. coli* anchor at 3,864,464 bytes (1.42% above unbounded, 0.82% below
JARVIS3) in a 54.8 s median compress and 25.8 s restore; that mode was not run
across the cross-taxon panel. The contribution is a working
system plus a reproducible measurement harness with byte-exact verification and
explicit deviations; the results do not establish universal superiority over
specialized compressors.

## 1. Introduction

Lossless genomic compression must preserve not only the four canonical bases
but also headers, line wrapping, case, ambiguous symbols, and scaffold gaps.
Pangenome archives add a second requirement: repeated and near-repeated
haplotypes should share storage, and users should be able to retrieve a locus
without expanding the entire cohort.

HyperGraSS addresses these requirements with four mechanisms:

1. a six-channel factorization that separates DNA from metadata;
2. a stride-sampled 16-mer anchor bank that aligns records in both orientations
   and stores substitutions by anchor coordinate;
3. an integer-only context-mixing coder for residual and non-referential
   sequence; and
4. a footer directory with per-block record offsets for direct queries.

The system is deliberately dependency-free: it compiles with `javac` and runs
on a stock JDK. The measured tables report the margins, including every cell
where HyperGraSS is not the smallest; they do not establish universal
superiority over specialized compressors.

## 2. System design

### 2.1 Six-channel factorization

`FastaStreamFactorizer` splits a block into:

| Channel | Content |
| --- | --- |
| 0 | DNA bases (A/C/G/T/U), 2-bit or model-coded |
| 1 | FASTA headers, front-coded against the previous header |
| 2 | Sequence line widths and line-ending types, run-length encoded |
| 3 | Lowercase (soft-masked) intervals |
| 4 | Scaffold `N` runs |
| 5 | Non-ACGTUN exception symbols |

Channels 1–5 are serialized with a raw/Deflate adaptive wrapper. The decoder
replays them in order, so headers, wrapping, case, and exceptions are restored
exactly. The format is specified in [`HGS-FORMAT.md`](HGS-FORMAT.md).

### 2.2 Stride-sampled 16-mer anchor bank

`CanonicalSyncmerEngine` indexes every 16-mer at a fixed stride (default 4) on
both strands, storing the forward and reverse-complement orientations with a
strand bit. The stride guarantees that any exact island of at least
`16 + stride - 1` bases contains an indexed seed (19 bp at the default). It
builds a bounded anchor bank (default at most 4 Mbp), and partitions the bank into
65,536-base pages. A target block is aligned to the bank; matches, strands,
literal runs, and SNPs are encoded as geometry plus residuals. The anchor bank
is stored once per archive and can be shared by all blocks. Global donor
catalogs let a block reuse variants observed in other clades.

The encoder chooses among four DNA payload modes: empty, pure context-mixer,
seed matches with SNP residuals, and anchor-coordinate haplotype matrices.
It falls back to the pure context mixer when referential coverage is low, and
otherwise compares the serialized mode 2 and mode 3 payloads rather than
choosing by a label.

### 2.3 Integer-only context mixing

`QuantumContextRepeatMixer` combines dense Markov orders, sparse hashed
contexts, spaced seeds, and codon-phase experts with a logistic/Bayesian mixer.
All arithmetic is integer, avoiding cross-CPU floating-point drift and making
archives deterministic for a fixed input and option set.

### 2.4 Container, footer, and queries

The `.hgs` container is documented in [`HGS-FORMAT.md`](HGS-FORMAT.md). A
multi-block archive stores a footer directory with each block's file offset,
compressed length, first record index, record count, and first record base
offset. A query reads the directory, locates the containing block or blocks,
decodes only those blocks and the anchor pages they reference, and returns the
requested slice. The footer is Deflate-compressed when that is smaller.

### 2.5 Bounded-memory mode

`--lean` uses the streaming block path, disables the global anchor bank, and
keeps per-block buffers bounded. The scale test measures the actual heap and
peak RSS; the implementation does not claim a specific bound beyond what that
test reports.

`--bounded` is a two-pass path for the anchor-bank format: the input is
factorized one block at a time, the anchor is built and polished by streaming
the file again, and the seed index is rebuilt only after the vote buffer is
released. The bounded encoder packs the four consensus counters per anchor
position into 4 bits and carries a 256 KiB gap history across block boundaries
(aux mode 7), which lets it run at a 160 MiB heap. On the 120.1 Mbp pangenome
with 16 MiB blocks it produces 3,997,232 bytes, 4.90% above the 64 MiB-block
profile. Its first run peaked at 252.7 MiB compressing and 221.2 MiB restoring
and appeared to meet the 256 MiB RSS gate. That pass did not survive repetition:
in the final three-repeat run every stride-4 attempt was OOM-killed at the
256 MiB cgroup cap (cgroup `oom_kill=1` each; wrapper RSS 269.3–269.4 MiB). It
therefore misses the RSS gate as well as the 2% size gate (110,615 bytes over).
The 4-bit counters saturate at 15 votes, which costs 2,218 bytes against exact
8-bit counts; the unbounded path keeps the exact form. Stride 8 at 16 MiB blocks
reduces the archive to 3,965,049 bytes, remains 1.76% larger than JARVIS3, and
is OOM-killed 3/3 by the same hard cap (269.3–269.4 MiB). Stride 8 at 32 MiB
blocks reaches 3,824,166 bytes (1.85% below JARVIS3) but needs at least 225-232
MiB of decoder heap; at 232 MiB it peaks at 375.4 / 286.6 MiB, and at 240 MiB
at 391.0 / 298.5 MiB. It is OOM-killed at 5.9 s under a strict 256 MiB cap (and
was OOM-killed 3/3 at 332.7-332.9 MiB under an earlier 320 MiB cap), so it is
not a 256 MiB JVM mode either. The anchor
header has a separate speed trade: on the 40.5 MB test, one context-mixing
anchor page takes 11.70 s, 64 KiB context-mixing pages (`--anchor-cm-pages`)
take 8.85 s for +0.46% archive, and 2-bit pages take 6.85 s for +2.33%.
Parallel context-mixing pages were slower because each page allocates its own
model, so that mode stays serial. At the 160 MiB bounded heap on a 120.1 Mbp
input, the paged mode cut compress wall from 144-158 s to 112-123 s and
decompress wall from 66-71 s to 40-53 s for +0.49% archive, with the same
272 MiB peak RSS in both modes (over the 256 MiB cap).

The bound is set by anchor saturation, not input size. The bank costs about
1 byte per base and its seed index 7.1 bytes per base (`SEED_STEP=4`, both
strands, 70% load), so a saturated 9 MiB anchor carries 67 MB of index before a
block is factorized. On the 122.9 Mbp *E. coli* medium stratum a 160 MiB heap
raises `OutOfMemoryError` in the factorizer; the 120.5 Mbp high stratum
completes the same profile at 272.2 MiB peak RSS, byte-exact but over the
256 MiB cap. A 2.08 GiB scale probe built from 55 copies of one sequence passes
because the bank never saturates. The JVM 256 MiB claim is therefore conditional
on pangenome novelty: at 120 MB of diverse input the 9 MiB anchor and the
256 MiB RSS cap cannot both hold in the JVM, and either the anchor or the heap
has to give.

A native image changes the runtime term without changing the anchor term. The
opt-in GraalVM 25 static-musl build (`build-native.ps1`) compiles the same
sources and keeps the archive format; `--split-tail-block` and the streaming
reconstruct path remove a 56 MiB tail block and a third block-sized output
array. With a 192 MiB heap, stride 8, 32 MiB blocks, and a strict 256 MiB
no-swap cgroup cap, the profile completed 3/3 at 3,864,464 bytes (1.42% above
the unbounded archive, 0.82% below JARVIS3), 54.8 s median compress at
240.1 MiB, 25.8 s restore at 225.5 MiB, every round trip byte-exact. The cgroup
compress peak leaves only 4.9 MiB of headroom. The native mode has not been run
on the six-collection panel, so the general claim remains open.

## 3. Evaluation protocol

### 3.1 Frozen corpus

Ten public datasets (mitochondria, plastid, yeast chromosomes, archaeal
reference, *E. coli*, *B. subtilis*, and a two-record *E. coli* pair) are
downloaded once and locked by SHA-256. The exact `.seq` bytes are the input
contract: FASTA headers and line endings are removed, but case, ambiguity
symbols, and record order are preserved. DNA-specialized tools run only when
every byte is A/C/G/T.

Every operation runs in a fresh process inside a pinned Linux/amd64 container
with 2 CPUs, 2 GB memory, and a 180 s cap. Compression and decompression are
timed separately; peak RSS is the child's `ru_maxrss`. Three repetitions are
used, except for the frozen GeCo3 level-9 profile, which is measured once
because a 1 MiB pilot exceeded 120 s. Supplementary GeCo3 levels 3 and 6 are
reported separately.

The output must match the input SHA-256 and length exactly. Three slice queries
per HyperGraSS archive must match the corresponding sequence bytes. The report
is then re-hashed and recomputed by an independent audit script.

The reported run was stopped deliberately after 179 of 280 planned results:
six datasets are complete, `ecoli-k12` is partial, and three datasets were not
measured. The report preserves its live `running` status, the audit is run with
`--allow-incomplete` to verify every completed result without certifying the
missing cells, and every table states its coverage. The corpus numbers are
therefore evidence about the completed cells, not a full 10-dataset benchmark.

### 3.2 Scale test

A 121.8 MB multi-FASTA pangenome is built from 24 complete RefSeq *E. coli*
records (120.1 Mbp). The manifest records accessions, lengths, and SHA-256. The
test compresses the complete FASTA (headers included), decompresses it, checks
SHA-256, and queries the first, middle, and last records. `--lean` runs with
`-Xmx8m` and a 256 MB container limit. A separate run measures the opt-in
`--pangenome` profile (16 MiB anchor bank, 8 MiB blocks). The anchor-bank and
block-size ablation that selected those values is reported in
[`ABLATION.md`](ABLATION.md).

### 3.3 Cross-taxon and diversity panel

The scale test and its two replications all use *Enterobacteriaceae*. To test
whether the pangenome margin survives unrelated phyla, six additional natural
collections were fetched with the same NCBI complete-genome rule: *M.
tuberculosis* (Actinomycetota, ~66% GC), *S. aureus* and *S. pneumoniae*
(Bacillota), *P. aeruginosa* (Pseudomonadota), *C. jejuni*
(Campylobacterota), and *C. trachomatis* (Chlamydiota). Each pool holds a
fixed-seed random sample of up to 72 eligible records. A natural collection is
then taken in farthest-point order from the pool medoid until it reaches about
120 MB, so the strain count varies with genome size.

Diversity is measured, not assumed. Each genome is sketched as the bottom
20,000 hashes of all canonical 21-mers; Jaccard similarity gives a Mash
distance and an approximate ANI. For *E. coli*, *S. aureus*, and *C. jejuni*,
three disjoint 24-genome strata are drawn from the same pool: low is the most
central members of the tightest cluster, high is a farthest-point sample of the
remainder, and medium is what remains. Every stratum records its mean, median,
and range of pairwise distance, cluster coverage, and ANI, and monotonicity of
the low/medium/high ordering is checked rather than assumed.

Natural collections and the nine diversity strata both use the 64 MiB-block
maximum-compression profile, so the gradient is not confounded with a profile
change. Both panels include JARVIS3 level 7, `xz -9`, and the bounded 256 MiB
profile; GeCo3 level 6 runs on three natural collections and the published
GraSS implementation on two. All inputs are normalized to one flat A/C/G/T
sequence, and every archive is verified against the input SHA-256 and a
streaming byte comparison.

### 3.4 Metrics

- archive bytes, compression ratio, and bits per input byte;
- wall time and child CPU time;
- peak child RSS;
- round-trip SHA-256 equality;
- query exactness;
- status, including `not_applicable` and `timeout`.

Medians are used for repeated operations. Failed and timed-out cells remain in
the tables. No failed cell is replaced by a synthetic value.

## 4. Results

All numbers below are rendered from the reports; the full tables, coverage
notes, and interpretation are in [`BENCHMARKS.md`](BENCHMARKS.md). The corpus
run was stopped at 179 of 280 planned results, so it supports statements about
six complete datasets plus a partial *E. coli* dataset, not about all ten.

**Sequence-only corpus.** On the seven datasets with a completed default
archive, the HyperGraSS archive was the smallest of every successful profile.
The margin over the best completed baseline has a median of +1.36% and a range
of +0.56% to +17.52%. Four margins are below 1.4%; the two largest are the
16.6 kB mitochondrial case (+17.52% versus xz) and the 1.66 MB archaeal genome
(+9.84% versus xz). The *E. coli* genome is a single-repetition cell:
HyperGraSS 1,087,650 bytes at 1.875 bits/base versus GeCo3 level 6 at
1,094,264 bytes and 1.886 bits/base. JARVIS3 level 26 and GeCo3 level 9 both
exceed the 180 s cap on that input. The default profile's median compression
time is 0.94 s versus 0.21 s for gzip and 0.12 s for xz.

**121.8 MB pangenome.** The default profile produced 9,215,764 bytes (0.605
bits per input byte) against 8,394,900 bytes for xz -9 (0.551 bits per byte):
the default is 9.8% larger. It was faster to compress (45.9 s versus 186.3 s)
and used more peak RSS (1,083 MiB versus 738 MiB). The opt-in `--pangenome`
profile produced 4,885,851 bytes (0.321 bits per byte), 41.8% smaller than xz,
with 30.6 s compression and 16.4 s decompression at 419 MiB peak RSS; its round
trip and all three record queries were exact. xz decompressed faster (7.1 s at
189 MiB) and has no query interface, so it cannot return a record slice without
full decompression. `--lean` completed the round trip and all three queries
with 52 MiB peak RSS, but its archive is 28,910,764 bytes, 3.14x the default.

**Independent replication.** Two additional ~120 Mbp pangenomes were fetched
with the same NCBI selection rule after the tuned profile was fixed. On
*Salmonella enterica* (25 records, 120.5 Mbp), the default produced 7,317,099
bytes (+19.56% versus xz at 6,120,260) and the tuned profile produced 3,443,600
bytes (43.73% smaller than xz). On *Klebsiella pneumoniae* (23 records,
122.7 Mbp), the default produced 8,571,498 bytes (+33.12% versus xz at
6,438,744) and the tuned profile produced 3,464,423 bytes (46.19% smaller).
Every round trip and all three slice queries per HyperGraSS archive were exact.
The median independent margin is 44.96%, so the direction and rough magnitude
of the original 41.8% result replicate. The default loses to xz on all three
pangenomes. Compression and restore times moved in the tuned profile's favor
on these runs, but host variance was large (the gzip baseline took 97.6 s on
Salmonella and 201.1 s on Klebsiella for nearly the same input size), so sizes
and round trips are the stable measurements.

**Cross-taxon replication.** To test whether the pangenome result is an
*Enterobacteriaceae* artifact, six additional natural collections were fetched
with the same rule, one each from Actinomycetota (*M. tuberculosis*), Bacillota
(*S. aureus*, *S. pneumoniae*), Pseudomonadota (*P. aeruginosa*),
Campylobacterota (*C. jejuni*), and Chlamydiota (*C. trachomatis*). Each targets
about 120 MB by farthest-point selection from the pool medoid; *C. trachomatis*
has only 66 eligible records and reaches 68.9 MB. On the identical flat
sequence, the 64 MiB HyperGraSS profile was smaller than JARVIS3 level 7 on all
six collections, with a median margin of 13.8% (range 4.7–29.4%), and every
round trip was byte-exact. The margin is not ordered by Mash diversity: the
largest are *C. jejuni* (d = 0.0190, 23.9%) and *C. trachomatis* (d = 0.0045,
29.4%), while *S. aureus* (d = 0.0125) has the smallest at 4.7%. The published
GraSS pipeline is 7.3–11.7x larger than HyperGraSS on the two collections where
it was run. This establishes that the advantage is not confined to one
bacterial family; it does not establish a monotone diversity law, which the
gradient below tests directly.

**Within-species diversity gradient.** Nine strata (three species, 24 genomes
each) test whether the margin depends on diversity. HyperGraSS was smaller than
JARVIS3 on eight of nine. *S. aureus* fell monotonically as mean Mash distance
rose: +11.14% (d = 0.00481), +5.07% (d = 0.01140), and +1.25% (d = 0.01519).
*E. coli* crossed zero at the medium stratum: +7.02% (d = 0.01188), -0.16%
(d = 0.02533), and +2.18% (d = 0.02595); its medium and high strata are nearly
identical in diversity, so that panel is effectively binary. *C. jejuni* was
non-monotonic: +17.57% (d = 0.01076), +19.67% (d = 0.01162), and +16.07%
(d = 0.01901), with low and medium close in distance. Against `xz -9` the
HyperGraSS archive stayed smaller on all nine strata (+9.3% to +50.4%). The
data do not support a single monotone law: the *S. aureus* decline is clean,
the other two are not, and Mash distance, genome size, and anchor coverage per
base move together in this design. The plausible mechanism is anchor coverage:
a fixed 9 MiB bank covers a smaller fraction of 5 Mb *E. coli* genomes than of
2.8 Mb *S. aureus* genomes. This panel does not separate that from genome size
or variable-site count; an anchor-size sweep on one species would.

**Specialized compressors.** The corpus comparison against JARVIS3 and GeCo3
uses single genomes. To test whether the pangenome margin survives against
tools designed for DNA, the same 120.1 Mbp sequence was extracted from the
pangenome (headers, line breaks, case, and 2,357 non-ACGT symbols removed or
normalized) and given to JARVIS3 level 7, GeCo3 levels 3 and 6, and xz -9. This
input is slightly smaller than the complete multi-FASTA given to HyperGraSS, so
the normalization favors the specialized tools. JARVIS3 level 7 produced
3,896,361 bytes (0.260 bits/base) and restored the sequence exactly in 84.2 s,
using 286 MiB peak RSS. The tuned HyperGraSS profile, extended with an LZ77
plus range-coder channel for variant gaps, produced **3,810,409 bytes** on the
same input, **2.21% smaller**, and restored it exactly in 82.4 s. The cost is
compression time and memory: 887.3 s versus 95.4 s and 999 MiB versus 286 MiB.
The scale result therefore shows a narrow size win over a specialized DNA
compressor at this scale, paid for with 9.3x the compression time, not the
earlier 20.3% deficit. On the *E. coli* anchor, the JVM 256 MiB bounded variant
does **not** meet the memory gate once measured repeatedly: 16 MiB blocks, a
160 MiB heap, packed 4-bit vote counters, and a 256 KiB cross-block gap history
produce 3,997,232 bytes (+4.90%), but all three repeats were OOM-killed at the
256 MiB cgroup cap (wrapper RSS 269.3–269.4 MiB). The earlier 252.7 MiB single
sample had swap headroom and is not a gate pass. It also misses the 2% size gate
by 110,615 bytes. Stride 8 at 16 MiB blocks gives 3,965,049 bytes and was also
OOM-killed 3/3 at 269.3–269.4 MiB; 32 MiB blocks reach 3,824,166 bytes (1.85%
below JARVIS3) but need at least 225-232 MiB of decoder heap; at 232 MiB they
peak at 375.4 / 286.6 MiB, and at 240 MiB at 391.0 / 298.5 MiB. They are
OOM-killed at 5.9 s under a strict 256 MiB cap (and were OOM-killed 3/3 at
332.7–332.9 MiB under an earlier 320 MiB cap). No JVM bounded profile completes
in a strict no-swap 256 MiB container.

The native image supplies the missing mode. Interleaved three-repeat against
JARVIS3, with the candidate under a strict 256 MiB no-swap cap, produced
**3,864,464 bytes** (1.42% above the unbounded archive, 0.82% below JARVIS3) in
a **54.8 s median compress at 240.1 MiB** and a **25.8 s restore at 225.5 MiB**.
All three round trips were byte-exact; the gate passes status, both RSS checks,
the 1.014 size ratio against 3,810,409 bytes, and the 0.940 wall ratio against
the JARVIS3 median. The cgroup compress margin is 4.9 MiB.

The published GraSS pipeline, the ancestor of HyperGraSS, was also run
unmodified on the identical sequence. It produced 24,777,592 bytes (1.650
bits/base), 6.50x the HyperGraSS archive and 6.36x JARVIS3, with 224.5 s
compress at 1,254.5 MiB peak RSS and 165.5 s restore at 1,019.3 MiB. JARVIS3
is smaller, faster in both directions, and lower-memory, so the original
grammar/substitution plus BSC design is not a competitive pangenome
compressor. The gain comes from the HyperGraSS anchor and variant-gap model,
not from the shared GraSS lineage.

**Interpretation.** HyperGraSS is competitive at small and medium sequence-only
scales and narrowly first on the completed corpus cells. With the default
profile it is not the smallest at 100 MB+ pangenome scale; the measured
ablation shows the 4 MiB anchor-bank cap was the ceiling. Raising the bank to
8 MiB at the same block size cut the archive by 43.9%, and 16 MiB with 8 MiB
blocks reached 4,885,851 bytes. The winning archive holds 144 anchor pages,
about 9.4 Mbp, so the 16 MiB cap is not fully consumed. The tuned profile also
beat JARVIS3 level 7 on six collections from five additional phyla, with a
median margin of 13.8%. This is a tuning result that replicates on nine
pangenomes spanning five phyla, not a new compression law; the opt-in profile
has not been tested on read data or outside bacteria. The 256 MiB bounded mode
holds its RSS cap on the *E. coli* anchor and on *C. trachomatis* but misses
the size gate on every natural collection and cannot complete on *P.
aeruginosa*, so it remains an *E. coli*-calibrated prototype rather than a
general mode. The demonstrated advantages are pure-Java implementation,
byte-exact round trips, an indexed query interface, and a consistent size win
over JARVIS3 on the tested pangenomes, not universal compression superiority.

**Query workload.** A separate workload served 13 deterministic 2,000-base
slices from the same pangenome, each in a fresh process. The tuned profile
answered in a median of 1.35 s against 1.36 s for xz and 1.18 s for gzip, so
the index did not beat sequential extraction on this input. It did tie xz while
storing 42% fewer bytes, but xz restored in 7.1 s versus 16.4 s, and roughly
859 queries per restore would be needed to offset that difference. `--lean`
was the fastest query profile (0.73 s median) at a 3.4x size cost versus xz.
The default profile was much slower (8.0 s median). Indexed access is a
capability with workload-dependent value, not an automatic speedup.

## 5. Limitations

- The corpus is locally conducted and not an independent third-party audit.
- The datasets are few and were chosen for size and record structure, not
  sampled from a broader distribution.
- The corpus excludes FASTA metadata; the scale test includes it. Ratios from
  the two sets are not directly comparable.
- GeCo3 level 9 is not a practical baseline under the 180 s cap, and the
  reported timeouts describe the cap rather than the tool's best size.
- HyperGraSS with the default profile does not beat xz on any of the three
  120 Mbp pangenomes; the anchor bank is capped below one record's length, so
  cross-record redundancy is only partly captured. The opt-in `--pangenome`
  profile does beat xz on all three, but they are all from one bacterial
  family, and each independent dataset has one measured run.
- The natural-collection panel covers nine pangenomes spanning five phyla, but
  each has one measured run on one host. The six cross-taxon collections are
  random RefSeq samples, not phylogenetically balanced panels, and the 120 MB
  target yields 19–72 records depending on genome size. Mash distance is an
  alignment-free approximation. The two best ablation configurations were
  repeated with byte-identical archives, while wall times varied by up to 2.5x
  across runs.
- The JVM 256 MiB bounded mode is not a general mode. It misses the 2% size gate
  on all seven natural collections (including the *E. coli* anchor), exceeds the
  RSS cap on four, and exhausts its heap on *P. aeruginosa* and the *E. coli*
  medium stratum. The heap limit is anchor saturation, not input size: a
  saturated 9 MiB bank carries a 67 MB seed index, so the 160 MiB heap sits on
  the OOM boundary, and the 120 MB high stratum completes at 272 MiB RSS. The
  absolute overhead is roughly constant, so the percentage penalty grows as the
  unbounded archive shrinks.
- The native 256 MiB mode is measured only on the *E. coli* anchor and only on
  one host. It passes 3/3 with 240.1 / 225.5 MiB worst RSS and 4.9 MiB of cgroup
  compress headroom, but it has not been run on the cross-taxon panel or on the
  large scale probe. The binary is a gitignored x86-64 Linux static-musl build
  made with GraalVM 25 while the default JVM build targets `--release 27`; the
  report records the binary size and SHA-256 for each run.
- The diversity gradient holds up to 24 strains per stratum and uses the same
  64 MiB profile as the natural collections. The strata differ in total size
  with genome size, and three species with three strata each support an
  ordering, not a significance test. The gradients are also uneven:
  *S. aureus* separates cleanly, but *C. jejuni* low/medium and *E. coli*
  medium/high are close in Mash distance.
- The specialized comparison is one run on one sequence. JARVIS3 level 26 and
  GeCo3 level 9 were infeasible at this scale under the 2 GiB container and
  180 s corpus conditions, so it is not a comparison against either tool's
  maximum setting. The specialized tools received header-free sequence input,
  which slightly favors them.
- GraSS was run unmodified with `-Xmx1500m` and its fixed `-e2`/25 MB BSC
  settings; the backend is bsc 3.4.0, newer than the revision used in the
  GraSS paper. Its tar container embeds F1/F2 modification times, so the
  archive size varies by a few bytes between runs.
- The query workload is one run of 13 slices on one host. It does not measure
  repeated queries, concurrency, or a long-lived server; JVM/Python startup is
  included in every query.
- Timing is from one Windows host with a Linux container; cross-platform
  behavior is not characterized.
- Earlier development claims (a 250-genome Ebola cohort and a synthetic
  2,000-strain outbreak) are not backed by retained artifacts in this
  repository and are not part of the results.
- The format is implementation-defined and versioned by a magic byte; it is not
  yet an external standard.

## 6. Related work

gzip and xz are general-purpose baselines. JARVIS3 and GeCo3 are specialized
DNA compressors with context-mixing models. PanMAN targets pangenome archives
and preserves additional alignment and tree information. HyperGraSS is
sequence-storage oriented; the PanMAN comparison in this repository is a
fixture-level archive-size result, not feature parity.

## 7. Reproducibility

The repository contains the pinned container definition, image digest, data
locks, corpus orchestrator, audit script, scale-test fetcher and runner, and
the report summarizer. Commands are listed in [`BENCHMARKS.md`](BENCHMARKS.md).
Downloaded sequences and run directories are gitignored; the locks and scripts
are committed.

## References

1. Roy, Maity & Mukhopadhyay. A lossless reference-free sequence compression
   algorithm leveraging grammatical, statistical, and substitution rules.
   *Briefings in Functional Genomics*, 2025, `elae050`.
2. JARVIS3, pinned commit `4dcb37b05b1a1482c2bf864d4968244817b0ad0a`.
3. GeCo3, pinned commit `69642a0a913e543599f35aee372e300601e4cbbe`.
4. PanMAN 0.1.4, as used by the pinned fixture runner.
5. GraSS implementation, commit `6d93700f5a2de8841f7cf2c20064413c80df2c4a`,
   with libbsc commit `1dc676540cb40a41562af6f2db04bd1adce97f90`.
