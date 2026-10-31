# HyperGraSS

**Next-Generation Lossless Genomic & Pangenomic Sequence Compression Engine (Pure Java 21, Zero External Dependencies)**

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
├── benchmarks/cost-aware/            # Byte profiles, size/CPU results, and compatibility logs
├── docs/HGS-FORMAT.md                # Byte-level `.hgs` version-1 container specification
├── docs/BENCHMARKS.md                # Measured corpus and scale results (generated from JSON)
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

The build uses only `javac` and JDK 21. `pom.xml` is an optional Maven packaging
path; CI uses the scripts. The test gate runs the lossless edge-case,
donor-index, anchor-codec, byte-accounting, and donor-cost suites. It does not
run cross-version tests that require a saved prior build.

---

## Quick Start (JDK 21+)

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

# Bounded-memory streaming mode; the measured scale result is in docs/BENCHMARKS.md
java -Xmx8m -cp bin hypergrass.HyperGraSS compress --lean input.fna archive_lean.hgs

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
[`docs/BENCHMARKS.md`](docs/BENCHMARKS.md). Every table in that document is
rendered from a machine-readable `report.json`; the corpus report is also
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
| 121.8 MB pangenome, HyperGraSS | 9,215,764 B, 13.22x, 45.9 s, 1,083 MiB peak RSS |
| 121.8 MB pangenome, xz -9 | **8,394,900 B, 14.51x**, 186.3 s, 738 MiB peak RSS |
| 121.8 MB pangenome, `--lean` (`-Xmx8m`) | 28,910,764 B, 4.21x, 52 MiB peak RSS |
| Verification | every round trip byte-exact; 3/3 record queries exact per HyperGraSS profile |

So the default archive is narrowly smallest on the completed small/medium
sequence-only cells, but **xz -9 is 8.9% smaller on the 100 MB+ pangenome**.
HyperGraSS is faster there and supports indexed record queries; xz is not
queryable. The bounded-memory mode fits a small heap at a large size cost.

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
```

The published scale report is committed at
`benchmarks/pangenome/runs/pangenome-scale-2/report.json`.

`prepare.py` and `fetch.py` download public data; the repository commits the
locks, not the sequences. The Docker image is pinned by digest in
`benchmarks/corpus/image.lock.json`.

## Format Specification

[`docs/HGS-FORMAT.md`](docs/HGS-FORMAT.md) documents the version-1 `.hgs`
container, footer directory, six-channel metadata layout, DNA modes, and query
semantics. The format is implementation-defined, not yet an external standard.

## License

Apache License 2.0. See [`LICENSE`](LICENSE). The license choice is a default
for this packaging pass, not legal advice; change it before publication if your
institution requires a different terms.
