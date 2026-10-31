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
121.8 MB complete-FASTA pangenome, xz -9 is 8.9% smaller, while HyperGraSS
answers exact record queries and finishes faster. The contribution is a working
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
`-Xmx8m` and a 256 MB container limit.

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

**121.8 MB pangenome.** xz -9 produced 8,394,900 bytes (0.551 bits per input
byte) against 9,215,764 bytes for HyperGraSS (0.605 bits per byte): xz is 8.9%
smaller. HyperGraSS was faster (45.9 s versus 186.3 s) and used more peak RSS
(1,083 MiB versus 738 MiB). `--lean` completed the round trip and all three
record queries with 52 MiB peak RSS, but its archive is 28,910,764 bytes, 3.14x
the default. The three HyperGraSS queries were exact. The xz archive has no
query interface, so it cannot return a record slice without full decompression.

**Interpretation.** HyperGraSS is competitive at small and medium sequence-only
scales and narrowly first on the completed corpus cells. It is not the smallest
at 100 MB+ pangenome scale. A likely cause is that the anchor bank is capped at
4.19 Mbp while each record is about 5 Mbp, whereas xz -9 can reference repeats
across a 64 MiB window of the 121.8 MB stream; this is a hypothesis, not a
measured ablation. The demonstrated advantages are pure-Java implementation,
byte-exact round trips, an indexed query interface, and a bounded-memory mode,
not universal compression superiority.

## 5. Limitations

- The corpus is locally conducted and not an independent third-party audit.
- The datasets are few and were chosen for size and record structure, not
  sampled from a broader distribution.
- The corpus excludes FASTA metadata; the scale test includes it. Ratios from
  the two sets are not directly comparable.
- GeCo3 level 9 is not a practical baseline under the 180 s cap, and the
  reported timeouts describe the cap rather than the tool's best size.
- HyperGraSS does not beat xz on the 121.8 MB pangenome; the anchor bank is
  capped below one record's length, so cross-record redundancy is only partly
  captured.
- The scale test is one pangenome from one species and one repetition per
  profile. Two runs produced identical archives but wall times differing by up
  to 2.5x, so its timing data is indicative.
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
