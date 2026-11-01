# HyperGraSS: Lossless Genomic and Pangenomic Compression with Canonical Syncmers and Context Mixing

**Working paper draft — not peer reviewed.** This document describes the
implementation and the measurement protocol. Numeric results are in
[`BENCHMARKS.md`](BENCHMARKS.md) and are generated from the reports in this
repository; they are not repeated as independent claims here.

## Abstract

HyperGraSS is a pure-Java 21 lossless compressor for FASTA, multi-FASTA, raw
DNA, and RNA sequence data. It factorizes input into six orthogonal channels,
encodes the DNA stream with a canonical syncmer referential model and an
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
result to the default's 4 MiB anchor-bank cap. On the identical sequence,
however, JARVIS3 level 7 produced a 3,896,361-byte archive, 20.3% smaller than
the tuned profile, so the xz margin is not a state-of-the-art DNA compression
result. The contribution is a working system plus a reproducible measurement
harness with byte-exact verification and explicit deviations; the results do
not establish universal superiority over specialized compressors.

## 1. Introduction

Lossless genomic compression must preserve not only the four canonical bases
but also headers, line wrapping, case, ambiguous symbols, and scaffold gaps.
Pangenome archives add a second requirement: repeated and near-repeated
haplotypes should share storage, and users should be able to retrieve a locus
without expanding the entire cohort.

HyperGraSS addresses these requirements with four mechanisms:

1. a six-channel factorization that separates DNA from metadata;
2. a canonical syncmer anchor bank that aligns records in both orientations
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

### 2.2 Canonical syncmer anchor bank

`CanonicalSyncmerEngine` selects canonical 16-mers on both strands, builds a
bounded anchor bank (default at most 4 Mbp), and partitions the bank into
65,536-base pages. A target block is aligned to the bank; matches, strands,
literal runs, and SNPs are encoded as geometry plus residuals. The anchor bank
is stored once per archive and can be shared by all blocks. Global donor
catalogs let a block reuse variants observed in other clades.

The encoder chooses among four DNA payload modes: empty, pure context-mixer,
syncmer matches with SNP residuals, and anchor-coordinate haplotype matrices.
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

### 3.3 Metrics

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

**Specialized compressors.** The corpus comparison against JARVIS3 and GeCo3
uses single genomes. To test whether the pangenome margin survives against
tools designed for DNA, the same 120.1 Mbp sequence was extracted from the
pangenome (headers, line breaks, case, and 2,357 non-ACGT symbols removed or
normalized) and given to JARVIS3 level 7, GeCo3 levels 3 and 6, and xz -9. This
input is slightly smaller than the complete multi-FASTA given to HyperGraSS, so
the normalization favors the specialized tools. JARVIS3 level 7 produced
3,896,361 bytes (0.260 bits/base) and restored the sequence exactly in 84.2 s,
using 286 MiB peak RSS. That is 20.3% smaller than the tuned HyperGraSS archive,
although HyperGraSS compressed faster (30.6 s versus 95.4 s) and restored faster
(16.4 s versus 84.2 s). The scale result therefore does not show that
HyperGraSS compresses pangenomes better than specialized DNA compressors; it
shows that the tuned profile beats xz, a general-purpose baseline, by a wide
margin, and that the margin over specialized tools is negative at this scale.

**Interpretation.** HyperGraSS is competitive at small and medium sequence-only
scales and narrowly first on the completed corpus cells. With the default
profile it is not the smallest at 100 MB+ pangenome scale; the measured
ablation shows the 4 MiB anchor-bank cap was the ceiling. Raising the bank to
8 MiB at the same block size cut the archive by 43.9%, and 16 MiB with 8 MiB
blocks reached 4,885,851 bytes. The winning archive holds 144 anchor pages,
about 9.4 Mbp, so the 16 MiB cap is not fully consumed. This is a tuning result
that replicates on three *Enterobacteriaceae* pangenomes, not a new compression
law; the opt-in profile has not been tested on read data or outside this
bacterial family. The demonstrated advantages are pure-Java implementation,
byte-exact round trips, an indexed query interface, and a bounded-memory mode,
not universal compression superiority.

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
- The scale test covers three pangenomes from one bacterial family. Most
  profiles have one measured run; the two best ablation configurations were
  repeated with byte-identical archives, while wall times varied by up to 2.5x
  across runs.
- The specialized comparison is one run on one sequence. JARVIS3 level 26 and
  GeCo3 level 9 were infeasible at this scale under the 2 GiB container and
  180 s corpus conditions, so it is not a comparison against either tool's
  maximum setting. The specialized tools received header-free sequence input,
  which slightly favors them.
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

1. Roy, Maity & Mukhopadhyay. GraSS: genomic sequence compression.
   *Briefings in Functional Genomics*, 2025, `elae050`.
2. JARVIS3, pinned commit `4dcb37b05b1a1482c2bf864d4968244817b0ad0a`.
3. GeCo3, pinned commit `69642a0a913e543599f35aee372e300601e4cbbe`.
4. PanMAN 0.1.4, as used by the pinned fixture runner.
