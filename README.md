# HyperGraSS

**Next-Generation Lossless Genomic & Pangenomic Sequence Compression Engine (Pure Java 27, Zero External Dependencies)**

HyperGraSS extends the **GraSS** pipeline (*Roy, Maity & Mukhopadhyay, Briefings in Functional Genomics, 2025, `elae050`*) by replacing ad-hoc ternary prefix/digram substitution and monolithic `libbsc` shells with a unified architecture grounded in biological sequence symmetry, information theory, and self-indexed pangenomic graph structures:

1. **Universal 6-Channel Orthogonal Lossless Factorization (`FastaStreamFactorizer.java`):** Separates metadata from the DNA stream instead of applying one substitution alphabet to the whole file, avoiding the error propagation that a single Kimura transition can cause in a shared digram code. The verification suites cover raw bytes, FASTA, multi-FASTA, RNA (`U`), chimeric DNA+RNA (`T+U`), IUPAC ambiguity codes, scaffold `N`-runs, soft-masked lowercase intervals, ragged line widths, and CRLF/LF line endings, and require byte-exact `SHA-256` identity on round trip.
2. **Bidirectional Canonical 16-mer Referential Engine with Anchor-Coordinate Haplotype Sets (`CanonicalSyncmerEngine.java`, v2):** Indexes both strands, bridges homologous SNPs, predicts circular-genome seams, and stores substitutions by absolute anchor coordinate so circularly permuted and reverse-strand genomes share clade variant sets. Earlier development notes report 190x on a 250-genome GenBank Ebola cohort and 2,087x on a synthetic 2,000-strain outbreak; those artifacts are **not retained in this repository**, so they are not counted as reproduced evidence here.
3. **Integer-Only Context-Mixing Range Coder (`QuantumContextRepeatMixer.java`, QCRM):** Dense Markov orders, sparse hashed k=12/15, spaced-seed and codon-phase experts, bidirectional repeat heads, and a logistic plus Bayesian mixer, all in integer arithmetic (no float/double, so no cross-CPU drift). The reproduced held-out comparison against gzip, xz, JARVIS3, and GeCo3 is in [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md); an earlier "11 of 14 genomes" development claim is not independently reproduced here.
4. **Page-Indexed Random Access & Bounded `--lean` Streaming (`HyperGraSS.java`):** Supports record and subsequence queries through a page index that decodes only the containing block, multi-core parallel block execution, and bounded-memory streaming. The `--lean` heap and 100 MB+ pangenome measurements are reported in [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md), not asserted here.

---

## Project Structure

```text
hypergrass/
├── src/hypergrass/
│   ├── HyperGraSS.java                # .hgs container engine, parallel pipeline, --lean streaming, page-indexed query CLI
│   ├── FastaStreamFactorizer.java     # 6-channel lossless FASTA/RNA/DNA factorizer & selective slice reconstructor
│   ├── CanonicalSyncmerEngine.java    # BC-LZG bidirectional (+/-) k-mer + Kimura SNP-bridging referential engine
│   ├── HaplotypeDonorIndex.java       # Exact inverted variant index + 64-lane bit-sliced donor counting
│   ├── CostAwareDonorSelector.java    # Bounded pre-Deflate byte-cost model for optional size search
│   ├── ArchiveByteAccounting.java     # Read-only exact encoded-byte breakdown; no payload decoding
│   ├── QuantumContextRepeatMixer.java # Integer-only context-mixing range coder
│   ├── BayesianMarkovRangeCoder.java  # Multi-order Markov (0..9) + codon-phase + Chargaff RC 64-bit range coder
│   ├── VarIntCodec.java               # LEB128/ZigZag varints, 2-bit DNA packing, adaptive raw/Deflate aux streams
│   ├── VerifyAllEdgeCases.java        # 14-suite SHA-256 byte-exact lossless verification harness
│   ├── VerifyDonorIndex.java          # Differential donor-count tests + optional synthetic microbenchmark
│   ├── VerifyDonorArchives.java       # Cross-version archive-byte checks + optional synthetic pipeline timings
│   ├── VerifyAnchorCodecs.java        # Bounded codec selection, size non-regression, previous-decoder checks
│   ├── VerifyArchiveAccounting.java  # Format coverage, exact totals, malformed-framing checks
│   ├── VerifyDonorCosts.java          # Independent raw-column serialization oracle for donor costs
│   ├── VerifyCostAwareEncoding.java  # Saved-build size, default-byte, query, and decoder checks
│   └── BenchmarkSuite.java            # Reproducible 5-experiment scientific benchmark suite
├── benchmarks/panman/                # Optional pinned comparison runner, reports, inputs, and command logs
├── benchmarks/corpus/                # Frozen held-out corpus, locked container, orchestrator, audit gate
├── benchmarks/pangenome/             # Reproducible >100 MB real E. coli scale test
├── benchmarks/ablation/              # Anchor-bank and block-size ablation harness
├── benchmarks/endtoend/              # Query-heavy end-to-end workload
├── benchmarks/cost-aware/            # Byte profiles, size/CPU results, and compatibility logs
├── docs/HGS-FORMAT.md                # Byte-level `.hgs` version-1 container specification
├── docs/BENCHMARKS.md                # Measured corpus and scale results (generated from JSON)
├── docs/ABLATION.md                  # Anchor-bank/block-size ablation results
├── docs/FOUNDATIONS.md               # Invariant/mechanism/evidence/falsification and knob policy
├── docs/END-TO-END.md                # Query-heavy end-to-end workload results
├── docs/RESEARCH_SINGLE_GENOME.md    # Algorithm research: self-contained vs reference/pangenome compression
└── README.md
```

---

## Build, Test, and Version

```powershell
# Windows
powershell -NoProfile -ExecutionPolicy Bypass -File .\build.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File .\test.ps1
java -cp bin hypergrass.HyperGraSS --version
```

```sh
# Linux / macOS
./build.sh
./test.sh
java -cp bin hypergrass.HyperGraSS --version
```

The build uses only `javac` and JDK 27. `pom.xml` is an optional Maven packaging
path; CI uses the scripts. The test gate runs the lossless edge-case,
donor-index, anchor-codec, byte-accounting, and donor-cost suites. It does not
run cross-version tests that require a saved prior build.

---

## Quick Start (JDK 27+)

### 1. Compile

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\build.ps1
# Linux/macOS: ./build.sh
```

### 2. Compress, Inspect, Query, and Decompress

```powershell
# High-compression parallel mode
java -cp bin hypergrass.HyperGraSS compress input.fna archive.hgs

# Optional bounded extra size search for small archives (more CPU; never retains a larger archive)
java -cp bin hypergrass.HyperGraSS compress --optimize input.fna archive_optimized.hgs

# Scale-tuned pangenome profile: 16 MiB anchor bank and 8 MiB blocks
java -cp bin hypergrass.HyperGraSS compress --pangenome input.fna archive_pangenome.hgs

# Single genome / one assembly: one QCRM superblock, no anchor probe
java -cp bin hypergrass.HyperGraSS compress --single ecoli.fna ecoli.hgs

# Same, split into N independently coded QCRM segments for parallel encode/decode
java -cp bin hypergrass.HyperGraSS compress --single --threads 4 ecoli.fna ecoli-4t.hgs

# Bounded-memory streaming mode; the measured scale result is in docs/BENCHMARKS.md
java -Xmx8m -cp bin hypergrass.HyperGraSS compress --lean input.fna archive_lean.hgs

# Bounded anchored mode above 2 GiB; one block resident, archive streamed to disk
java -Xmx160m -cp bin hypergrass.HyperGraSS compress --bounded --anchor-mib 9 --block-mib 16 big.fna big.hgs

# Same, encoding independent blocks on N threads (drops the shared gap channel)
java -Xmx4g -cp bin hypergrass.HyperGraSS compress --bounded --threads 8 --anchor-mib 64 --block-mib 64 big.fna big.hgs

# Native 256 MiB mode (opt-in): build the static-musl image, then run it
powershell -NoProfile -ExecutionPolicy Bypass -File build-native.ps1
bin/native/hypergrass -Xmx192m -Dhypergrass.matchMin=104 -Dhypergrass.rangeModel=count compress --pangenome --bounded --single-thread --anchor-single-page --seed-step 8 --block-mib 32 --anchor-mib 9 --split-tail-block big.fna big.hgs

# Inspect archive metadata in O(1)
java -cp bin hypergrass.HyperGraSS inspect archive.hgs

# Exact encoded-byte breakdown (reads the archive into memory; not an O(1) inspection)
java -cp bin hypergrass.HyperGraSS profile archive.hgs

# Page-indexed random-access query: extract record #1337, locus [1000, 1500)
java -cp bin hypergrass.HyperGraSS query archive.hgs 1337 1000 1500

# Full lossless decompression (with automatic CRC-32 integrity verification)
java -cp bin hypergrass.HyperGraSS decompress archive.hgs restored.fna
```

### 3. Run Verification & Reproducible Benchmarks

```powershell
java -cp bin hypergrass.VerifyAllEdgeCases
java -cp bin hypergrass.VerifyDonorIndex
java -Xmx512m -cp bin hypergrass.VerifyAnchorCodecs
java -Xmx512m -cp bin hypergrass.VerifyArchiveAccounting
java -Xmx512m -cp bin hypergrass.VerifyDonorCosts
java -Xms512m -Xmx512m -cp bin hypergrass.VerifyDonorIndex --bench
java -cp bin hypergrass.BenchmarkSuite
java -Xmx3g -cp bin hypergrass.RealCorpusBench <dir-with-genomes>
```

---

## Verified PanMAN Archive-Size Result

On PanMAN's official **20-record demo**, complete archive sizes were:

| Encoder | Archive bytes |
| --- | ---: |
| Released PanMAN 0.1.4 | 13,816 |
| HyperGraSS before anchor-codec selection | 10,665 |
| HyperGraSS default | 10,419 |
| HyperGraSS `--optimize` | **10,257** |

Default HyperGraSS is **24.59% smaller than PanMAN on this fixture**. Optional `--optimize` reaches **25.76% smaller than PanMAN**, or **1.55% smaller than the current default**. The earlier anchor-codec selection change reduced the original HGS size by 2.31%. All 20 identifiers, full headers, and exact sequence bytes matched PanMAN's extraction. HGS also passed byte-exact FASTA/SHA-256 reconstruction and 60 slice queries.

The change is bounded **anchor-codec selection**, not a claim of a fundamental algorithmic breakthrough: for banks fitting in one 65,536-base page, try both existing anchor codecs and keep the smaller one. The format and decoder are unchanged; larger banks and `--lean` retain their original paths. Extra context mixing trades CPU time for archive size and can slow decoding/queries when it replaces packed anchors.

**Scope:** PanMAN preserves additional alignment/tree/evolutionary information, so this is a sequence-storage size result, not feature parity or worldwide superiority. A larger published archive could not be compared because the pinned PanMAN loader exceeded its message traversal limit; that failure is recorded, not counted as a win. Speed and peak-memory superiority have not been established.

See [`benchmarks/panman/README.md`](benchmarks/panman/README.md) for the pinned input/image hashes, reproducible Python/Docker runner, machine-readable report, compatibility tests, and failed larger-data attempt. The engine itself still has no external dependencies.

---

## Byte Accounting and Optional Cost-Aware Search

`profile` partitions every encoded archive byte into container/directory, anchor, metadata, literal, match, variant, and donor categories, with separate DNA-mode counts. It skips opaque compressed streams without inflating or decoding them. Totals must equal the file size. This validates framing, **not** compressed contents, directory semantics, or CRC integrity; `decompress` remains the reconstruction check. The API is `ArchiveByteAccounting.report(Path)` or `report(byte[])`.

`compress --optimize` tries cost-aware donor choices and an archive without the global donor dictionary. Donor costs count actual varint widths, rounded mask bytes, and added-position gaps/differences before Deflate. Because those estimates do not predict final Deflate sizes exactly, **complete serialized archives, including compressed footers, are compared against the default archive**. The default wins ties.

The extra search applies only to non-lean inputs up to 2 MiB, at most 16 blocks, one anchor page, and a nonempty global donor dictionary. Individual donor trials are bounded to 8,192 differences, 256 matches, and 1,048,576 modeled scan units; exhausted trials are discarded. This budget does not include intersection/index construction or all encoder work. Larger/excluded inputs and `--lean` retain the default path. Default compression and the public four-argument block encoder keep their existing archive bytes.

The gain is modest: the demo improves by 162 bytes; only one of twelve additional synthetic regression cases improves (18 bytes), with none larger. A separate four-case synthetic size/timing run improved two archive sizes but took substantially more CPU, including cases with no size improvement. **That is why the search is opt-in, not a default performance regression.** The format, decoder, query dependency depth, and engine dependencies are unchanged.

See [`benchmarks/cost-aware/README.md`](benchmarks/cost-aware/README.md) for measurements, scope, logs, and saved-build regression commands. Representative-reference selection, cost-aware match parsing, and bounded-memory global sharing remain future work; this is the first measured stage, not the entire research roadmap.

---

## Exact Bit-Sliced Haplotype Donor Search

The measurements below describe the earlier index-only optimization, not the optional cost-aware search.

`HaplotypeDonorIndex` accelerates both global clade harvesting and Mode 3 donor selection without changing the archive format, donor scoring, insertion order, or first-donor tie-breaking. It is an encoder optimization, not a new compression format or a compression-ratio improvement.

Instead of merging every donor's entire variant list with each sweep, it indexes the **full packed variant `(anchorPos << 2) | xorDiff`**. Each posting holds a 64-donor membership mask. Queries add these masks into bit-sliced binary counters using XOR/AND carry propagation, counting up to 64 donor intersections in parallel. Counts remain exact; different alleles at the same coordinate remain distinct. Only donors sharing queried variants contribute to the counters.

The index is lazy and local to an encoder operation, so parallel blocks do not share mutable state. Small pools (fewer than 16 donors), tiny initial queries (fewer than 4 variants), and pools exceeding 131,072 stored variant memberships use the original sorted-merge approach. The membership budget limits additional index storage; index construction must still be amortized across queries. Decoding and `--lean` compression are unchanged.

### Measured Results, Not Worldwide Superiority Claims

One local run on Windows with Oracle JDK 21.0.6 and `-Xms512m -Xmx512m` produced these **synthetic donor-count microbenchmark** medians. Each batch performs 64 queries against 512 donors containing 65,536 stored variants; both implementations consume every donor count and produce identical checksums. There are 4 warmup batches and 7 measured batches, with execution order alternating.

| Synthetic workload | Sorted-merge scan | Bit-sliced index | Warmed speedup |
| --- | ---: | ---: | ---: |
| Sparse-related | 4.858 ms | 0.209 ms | 23.2x |
| Dense-shared | 4.374 ms | 0.247 ms | 17.7x |
| Disjoint | 3.487 ms | 0.064 ms | 54.5x |

These timings **exclude index construction**, which the harness reports separately. First queries including lazy construction took 2.8–5.8 ms in this run. These are not whole-compressor speedups and do not establish superiority over other compressors.

A separate original-versus-updated encoder comparison used four deterministic 1,559,186-byte synthetic pangenomes, with reverse complements and circular rotations. Serial full-archive compression medians (4 warmups, 7 measured runs) were:

| Synthetic cohort | Original encoder | Updated encoder |
| --- | ---: | ---: |
| Private variants | 45.0 ms | 49.9 ms |
| 30 shared variants plus private variants | 60.3 ms | 40.9 ms |
| 90 shared variants plus subclades/private variants | 112.2 ms | 77.8 ms |
| 90 shared variants plus few private variants | 49.7 ms | 37.6 ms |

All four archives were **byte-for-byte identical**, and lossless reconstruction and 60 slice queries passed. The private-variant case was slower: end-to-end performance is workload-dependent, and these small synthetic timings are not a substitute for a real-corpus benchmark. All 14 lossless verification suites also passed.

### Reproduce the Checks

`VerifyDonorIndex` needs no external data. Its 8 deterministic suites compare 1,928,153 donor intersection counts and 10,039 selections against a sorted-merge oracle, including lazy additions, memory-budget fallback, allele distinctions, carry boundaries, signed mask lanes, and scoring ties.

For cross-version comparison, retain the **original** compiled encoder classes in `original-bin/`. On Windows, placing that directory first lets the same test harness run against the old implementation:

```powershell
# Requires an actual saved original build in original-bin/; do not copy the updated build there.
java -Xms512m -Xmx512m -cp "original-bin;bin" hypergrass.VerifyDonorArchives baseline donor-baseline --bench
java -Xms512m -Xmx512m -cp bin hypergrass.VerifyDonorArchives verify donor-baseline --bench
```

Use `:` instead of `;` in the classpath on Linux/macOS. Omit `--bench` to run archive-byte, round-trip, and query checks without repeated timing runs. `baseline` writes the four generated archives; `verify` compares against them and fails if even one archive byte differs.

---

## Reproduced Benchmarks

The frozen corpus and the 100 MB+ pangenome scale test are documented in
[`docs/BENCHMARKS.md`](docs/BENCHMARKS.md); independent replication on two
additional pangenomes is in
[`docs/INDEPENDENT.md`](docs/INDEPENDENT.md), the anchor-bank ablation is in
[`docs/ABLATION.md`](docs/ABLATION.md), the query workload is in
[`docs/END-TO-END.md`](docs/END-TO-END.md), and the specialized-compressor
scale check is in [`docs/SPECIALIZED.md`](docs/SPECIALIZED.md). Every table in
`BENCHMARKS.md` and `END-TO-END.md` is rendered from a machine-readable
`report.json`; the independent, specialized, and ablation tables are generated
from or transcribed from the committed reports. The corpus report is also
re-hashed and recomputed by `benchmarks/corpus/audit.py` before publication.
The corpus is a locally conducted held-out benchmark against pinned tool
versions, not an independent third-party audit. Deviations, timeouts, and
non-applicable cells are reported rather than dropped.

Measured headlines, stated without rounding up:

| Result | Measurement |
| --- | --- |
| Corpus coverage | 179 of 280 planned results; 6 of 10 datasets complete |
| Smallest archive on completed corpus cells | 7 of 7 datasets |
| Median margin over the best completed baseline | +1.36% (range +0.56% to +17.52%) |
| 121.8 MB pangenome, default HyperGraSS | 9,215,764 B, 13.22x, 45.9 s, 1,083 MiB peak RSS |
| 121.8 MB pangenome, `--pangenome` (tuned, multi-FASTA) | **4,885,851 B, 24.93x**, 30.6 s, 1,008 MiB peak RSS |
| 121.8 MB pangenome, xz -9 | 8,394,900 B, 14.51x, 186.3 s, 738 MiB peak RSS |
| 122.3 MB Salmonella pangenome, tuned | 3,443,600 B, 35.50x, 28.6 s, **43.7% smaller than xz** |
| 124.5 MB Klebsiella pangenome, tuned | 3,464,423 B, 35.92x, 10.6 s, **46.2% smaller than xz** |
| Cross-taxon panel (6 collections, 5 phyla) | **6/6 smaller than JARVIS3 l7**, median 13.8% (range 4.7–29.4%), byte-exact |
| Cross-taxon panel, published GraSS | 7.3–11.7x larger than HyperGraSS on the two collections measured |
| Cross-taxon panel, 256 MiB bounded mode | 0/6 pass the size gate; 4/6 exceed 256 MiB RSS; 1 exhausts the heap |
| Within-species gradient (3 species × 3 strata) | HyperGraSS smaller than JARVIS3 on **8/9**; *E. coli* medium is the first loss (−0.16%), byte-exact |
| Single genome, *E. coli* K-12 (4.6 Mb), `--single` | **1,084,172 B, 1.869 bits/base**, 17.3 s, 88 MiB peak RSS; **5.03% smaller than GraSS**, 3.2x lower RSS |
| Single genome, *B. subtilis* 168 (4.2 Mb), `--single` | **994,834 B, 1.888 bits/base**, 14.6 s, 87 MiB peak RSS; **3.80% smaller than GraSS**, 2.8x lower RSS |
| Single genome, *P. aeruginosa* NY4383 (7.3 Mb), `--single` | **1,529,536 B, 1.673 bits/base**, 20.5 s, 95 MiB peak RSS; **9.52% smaller than GraSS**, 3.1x lower RSS |
| Single genome, human chr1 first 100 Mb (0.3% Ns), `--single` | **20,626,877 B, 1.650 bits/base**, 355.0 s, 396 MiB peak RSS; **7.35% smaller than GraSS**, 2.25x lower RSS; 2.58x slower |
| Single genome, human chr1 first 100 Mb, `--single --threads 2` | **20,635,964 B, 1.651 bits/base**, 70.1 s, 586 MiB peak RSS; 0.044% above serial but **5.07x faster**; 1.96x faster than GraSS with a 7.31% smaller archive |
| Single genome, human chr1 first 100 Mb, `--single --threads 4 / 8` | 20,649,693 / 20,677,028 B, 1.652 / 1.654 bits/base, 124.9 / 225.9 s; both slower than `--threads 2` on a 2-CPU container |
| Single genome, human chr1 first 100 Mb, default (no flag) | 21,056,843 B, 1.684 bits/base, 152.9 s, 967 MiB peak RSS; 2.08% larger but 2.32x faster than `--single`; 5.42% smaller than GraSS |
| Single genome, human chr1 first 100 Mb, `--lean` | 22,788,010 B, 1.823 bits/base, **40.5 s, 94 MiB peak RSS**; 2.35% larger than GraSS but **3.40x faster at 9.5x lower RSS** |
| Single genome, human chr14 first 100 Mb (16.3% Ns), `--single` | **17,084,921 B, 1.367 bits/base**, 115.6 s, 380 MiB peak RSS; **7.84% smaller than GraSS**, 2.5x lower RSS; 1.41x slower |
| Single genomes, `--single --threads 2` | 0.39–1.23% larger than serial `--single`; 1.3–1.9x faster on the bacterial genomes, byte-exact |
| Same 120.1 Mbp sequence, tuned + LZ gap coder | **3,810,409 B, 0.254 bits/base, 31.52x**, 887.3 s, 999 MiB peak RSS |
| Same 120.1 Mbp sequence, JARVIS3 l7 | 3,896,361 B, 0.260 bits/base, 95.4 s, 286 MiB peak RSS |
| Same sequence, GraSS (public implementation, `-e2`) | 24,777,592 B, 1.650 bits/base, 4.85x, 224.5 s, 1,255 MiB compress RSS |
| Same sequence, bounded 256 MiB attempt | 4,131,860 B, 272.5 MiB, 149.5 s; missed size and RSS gates |
| Same sequence, bounded tuned (stride 4, 16 MiB blocks) | 3,997,232 B, 57.1 s; OOM-killed 3/3 at the 256 MiB cap (269.4 MiB wrapper RSS); misses size and RSS gates (+4.90%) |
| Same sequence, bounded stride 8 (16 MiB blocks) | 3,965,049 B, 51.2 s; OOM-killed 3/3 (269.4 MiB); 1.76% larger than JARVIS3, misses both gates |
| Same sequence, bounded stride 8 (32 MiB blocks) | **3,824,166 B, 1.85% below JARVIS3**, 52.3 s; needs 225-232 MiB decoder heap, peaks 375.4 / 286.6 MiB, OOM-killed at 5.9 s under a strict 256 MiB cap |
| Same sequence, bounded native (GraalVM 25, 192 MiB heap) | **3,864,464 B**, 1.42% above unbounded and 0.82% below JARVIS3; 54.8 s / 240.1 MiB compress, 25.8 s / 225.5 MiB restore, **3/3 under a strict 256 MiB cap** |
| 121.8 MB pangenome, `--lean` (`-Xmx8m`) | 28,910,764 B, 4.21x, 52 MiB peak RSS |
| Verification | every round trip byte-exact; 3/3 record queries exact per HyperGraSS profile |

On a single assembly the opt-in `--single` profile is the intended choice: it
beat GraSS by 3.80–9.52% on four genomes from 1.36 Mb to 7.32 Mb, used 2.2–3.2x
less peak RSS, and stayed byte-exact on every repeat. Serial compression is
about twice GraSS's time; `--threads 2` recovers most of that on the bacterial
genomes for 0.39–1.23% more bytes. On 100 Mb human slices the ratio and memory
margins hold (**7.35% and 7.84% smaller**, 2.25x and 2.5x lower RSS), but serial
compression is 2.58x and 1.41x GraSS's time; the tuned path is not a serial
speed win on repeat-rich sequence of that size. The un-flagged default is the
speed-oriented choice at 100 Mb: it used 47 blocks with a 4.19 Mbp anchor bank,
produced an archive 2.08% larger than `--single` but 2.32x faster, and stayed
5.42% below GraSS. `--single --threads 2` is the practical 100 Mb choice: 0.044%
more bytes than serial for 5.07x the speed, and 1.96x faster than GraSS with a
7.31% smaller archive. More segments oversubscribe the 2-CPU container:
`--threads 4` and `--threads 8` took 124.9 s and 225.9 s for +0.111% and
+0.243% bytes. `--lean` is the speed/memory extreme: 40.5 s and 94 MiB, at
2.35% above GraSS but 3.40x faster. The pangenome profile is not a substitute for a single
record: on the 100 Mb chr1 slice it is 10.22% larger and 2.10x slower than
`--single`, and 2.12% larger than GraSS.

So the **default** archive is narrowly smallest on the completed small/medium
sequence-only cells, but **xz -9 is 9.8% smaller than the default on the
100 MB+ pangenome**. The opt-in `--pangenome` profile reverses that result:
4,885,851 B is **41.8% smaller than xz -9** on the original pangenome, and the
same profile was **43.7% and 46.2% smaller** on two independent *Salmonella*
and *Klebsiella* pangenomes (median 44.96%). It also compressed faster than xz
in those runs and answered exact record queries. xz decompressed faster on the
original pangenome; restore time was mixed on the independent sets. The profile
remains opt-in: nine pangenomes from five phyla, one host, and one run each are
a sample rather than a guarantee for every input. The bounded-memory mode is not
a general 256 MiB mode on the JVM. On the *E. coli* anchor its first single run
peaked at 252.7 MiB compressing and 221.2 MiB restoring, but that pass had swap
headroom. The final three-repeat JVM run OOM-killed every bounded attempt at the
cgroup cap: stride 4 and stride 8 with 16 MiB blocks died at 256 MiB
(269.3–269.4 MiB wrapper RSS). It also gives back 4.90% in size. Stride 8 at
16 MiB blocks recovers 0.80% but remains 1.76% larger than JARVIS3 and is
OOM-killed 3/3. The only JVM profile that beats JARVIS3 (3,824,166 B, 1.85%
below) uses 32 MiB blocks, needs 225-232 MiB of decoder heap, and peaks at
375.4 / 286.6 MiB (worst observed 391.0 / 298.5 MiB). It is OOM-killed at 5.9 s
under a strict 256 MiB cap. No JVM bounded profile completes in a strict no-swap
256 MiB container, so there is no ≤256 MiB JVM mode that beats JARVIS3.

The opt-in native image changes that specific verdict. With GraalVM 25, a
192 MiB heap, stride 8, 32 MiB blocks, and `--split-tail-block`, HyperGraSS
completed 3/3 under a strict 256 MiB no-swap cap on the same input and image:
**3,864,464 B** (1.42% above the unbounded archive, 0.82% below JARVIS3),
54.8 s median compress at 240.1 MiB, 25.8 s restore at 225.5 MiB, every round
trip byte-exact. The cgroup compress peak leaves only 4.9 MiB of headroom. The
JVM path remains unable to meet the cap, and the native mode was not rerun on
the cross-taxon panel; across those six collections the JVM bounded mode misses
the 2% size gate on all six, exceeds 256 MiB RSS on four, and runs out of heap
on *P. aeruginosa*.

**The margin over the best specialized tool is real but narrow.** On the
identical 120.1 Mbp sequence, the tuned HyperGraSS archive is **2.21% smaller
than JARVIS3 level 7** (3,810,409 versus 3,896,361 B), and both round-trip
exactly. The cost is 9.3x longer compression (887.3 versus 95.4 s) and 3.5x the
compression RSS (999 versus 286 MiB); restore time is equal (82.4 versus
84.2 s). That is a size win, not a speed win, and it is not a compression
breakthrough: 2.21% on one pangenome is a narrow margin, and the 41.8% xz
margin is measured against a general-purpose baseline rather than the state of
the art for DNA. The published GraSS implementation, the ancestor pipeline, is
**6.50x the size of HyperGraSS and 6.36x the size of JARVIS3** on the same
input, and JARVIS3 also beats it on compression time, restore time, and memory,
so the specialized frontier remains JARVIS3. See
[`docs/SPECIALIZED.md`](docs/SPECIALIZED.md).

**The advantage is not an *E. coli* artifact.** Six unrelated natural
collections were fetched with the same rule, one each from Actinomycetota,
Bacillota, Pseudomonadota, Campylobacterota, and Chlamydiota, and given the
identical flat sequence. HyperGraSS was smaller than JARVIS3 level 7 on all
six, with a median margin of 13.8% (range 4.7–29.4%) and byte-exact round
trips. The margin is not ordered by Mash diversity: *C. jejuni* (d = 0.0190)
and *C. trachomatis* (d = 0.0045) have the largest margins, while *S. aureus*
(d = 0.0125) has the smallest. The published GraSS pipeline is 7.3–11.7x
larger than HyperGraSS on the two collections where it was measured.

**The advantage shrinks, but does not vanish, as within-species diversity
rises.** A nine-run gradient drew low/medium/high strata from the same
72-genome pool for three species, with identical flat sequences. HyperGraSS
beat JARVIS3 on eight of nine strata. *S. aureus* fell monotonically from
+11.14% to +5.07% to +1.25%; *E. coli* crossed zero at the medium stratum
(−0.16%) and recovered to +2.18% at high; *C. jejuni* stayed between +16% and
+20%. The gradient is not one curve: Mash distance, genome size, and anchor
coverage per base move together in this panel, so it does not separate them.
Protocol, per-stratum tables, and limits are in [`docs/BROAD.md`](docs/BROAD.md).

On the query workload (13 exact 2,000-base slices), the tuned profile tied xz's
query latency (1.35 s versus 1.36 s median) while storing 42% fewer bytes, but
xz restored faster; about 859 queries per restore would be needed to offset
that difference. The default profile was much slower at queries (8.0 s median).
See [`docs/ABLATION.md`](docs/ABLATION.md) and
[`docs/END-TO-END.md`](docs/END-TO-END.md).

The published report is committed at
`benchmarks/corpus/runs/corpus-full-1/report.json`; inputs and artifacts are
not. It was stopped early, so it still says `running`.
`python -B benchmarks/corpus/audit.py <run-dir> --allow-incomplete` re-verifies
every completed result without claiming the missing cells.

To regenerate the corpus report from scratch (requires Docker and downloads):

```sh
python -B benchmarks/corpus/prepare.py --verify
python -B benchmarks/corpus/run.py --run-id corpus-repro-1
python -B benchmarks/corpus/audit.py benchmarks/corpus/runs/corpus-repro-1
python -B benchmarks/summarize.py --corpus benchmarks/corpus/runs/corpus-repro-1/report.json
```

To rebuild and run the scale test:

```sh
python -B benchmarks/pangenome/fetch.py
python -B benchmarks/pangenome/run.py --run-id pangenome-scale-repro
python -B benchmarks/pangenome/run.py --run-id pangenome-tuned-repro --profiles hypergrass-pangenome
python -B benchmarks/render_benchmarks.py \
  --corpus benchmarks/corpus/runs/corpus-full-1/report.json \
  --scale benchmarks/pangenome/runs/pangenome-scale-repro/report.json \
  --scale benchmarks/pangenome/runs/pangenome-tuned-repro/report.json
python -B benchmarks/endtoend/run.py --run-id endtoend-repro
```

The published scale reports are committed at
`benchmarks/pangenome/runs/pangenome-scale-2/report.json` (default, lean, gzip,
xz) and `benchmarks/pangenome/runs/pangenome-scale-3/report.json` (tuned
profile). The anchor-bank ablation is documented in
[`docs/ABLATION.md`](docs/ABLATION.md).

`prepare.py` and `fetch.py` download public data; the repository commits the
locks, not the sequences. The Docker image is pinned by digest in
`benchmarks/corpus/image.lock.json`.

## Scaling Beyond 2 GiB

The bounded path previously rejected inputs above 2 GiB because block offsets were
32-bit and the archive was assembled in one `byte[]`. Both are fixed: block offsets are
`long`, and the writer streams header, blocks and footer to disk.

Measured on a 2.08 GiB synthetic input (`-Xmx160m`, 16 MiB blocks, 9 MiB anchor):

| Metric | Value |
|---|---|
| Input / archive | 2,229,710,670 B / 30,171,539 B |
| Blocks | 132 |
| Compress | 864.5 s, peak sampled working set 235 MB |
| Restore | 238.2 s, CRC-32 OK, SHA-256 identical |

Archive bytes are unchanged for inputs below 2 GiB: the campy-high bounded archive
still reproduces as `1,125,629 B` with SHA-256 `14DD7E94…E1C5D6`.

Throughput, not memory, is now the limit: 2.58 MB/s single-threaded is about 65 minutes
for 10 GB and 4.5 days for 1 TB. `--threads N` parallelises the encode pass, but the
bounded path makes four passes over the input and only one is parallel, so threads
alone cap out well below linear. The anchor header has a separate trade:
`--anchor-cm-pages` uses 64 KiB context-mixing pages and costs 0.46% in archive size
for 1.32x in wall on the 40.5 MB test; at the 160 MiB bounded heap on a 120 Mbp input
it costs 0.49% and cuts compress wall by 1.2-1.3x. Fast 2-bit pages are the fastest
mode at +2.33% and are already the default when `--anchor-single-page` is omitted. CM pages are
encoded serially because parallel pages were slower. The remaining work - threads,
sharded anchors, and resumability - is sequenced in [`docs/SCALING.md`](docs/SCALING.md).
A native image has since been measured for the 256 MiB gate (above); it is a
deployment choice that removes JVM runtime overhead, not an algorithm change.

## Format Specification

[`docs/HGS-FORMAT.md`](docs/HGS-FORMAT.md) documents the version-1 `.hgs`
container, footer directory, six-channel metadata layout, DNA modes, and query
semantics. The format is implementation-defined, not yet an external standard.

## License

Apache License 2.0. See [`LICENSE`](LICENSE). The license choice is a default
for this packaging pass, not legal advice; change it before publication if your
institution requires a different terms.
