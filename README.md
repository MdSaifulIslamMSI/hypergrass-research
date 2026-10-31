# HyperGraSS

**Next-Generation Lossless Genomic & Pangenomic Sequence Compression Engine (Pure Java 21, Zero External Dependencies)**

HyperGraSS supersedes the **GraSS** pipeline (*Roy, Maity & Mukhopadhyay, Briefings in Functional Genomics, 2025, `elae050`*) by replacing ad-hoc ternary prefix/digram substitution and monolithic `libbsc` shells with a unified mathematical architecture grounded in biological sequence symmetry, information theory, and self-indexed pangenomic graph structures:

1. **Universal 6-Channel Orthogonal Lossless Factorization (`FastaStreamFactorizer.java`):** Eliminates GraSS's **Parity-Avalanche Catastrophe** (where a single Kimura transition SNP corrupts >92% of downstream `F2.SSG` digram bytes) and achieves 100% `SHA-256` lossless bijectivity over arbitrary Raw, FASTA, Multi-FASTA, RNA (`U`), chimeric DNA+RNA (`T+U`), IUPAC ambiguity codes, scaffold `N`-runs, soft-masked lowercase intervals, ragged line widths, and CRLF/LF line endings.
2. **Bidirectional Canonical 16-mer Referential Engine with Anchor-Coordinate Haplotype Sets (`CanonicalSyncmerEngine.java`, v2):** Indexes both strands, bridges homologous SNPs, predicts circular-genome seams, and stores substitutions by absolute anchor coordinate so circularly permuted and reverse-strand genomes share clade variant sets (real 250-genome GenBank Ebola cohort: 190x; synthetic 2,000-strain outbreak: 2,087x).
3. **Integer-Only Context-Mixing Range Coder (`QuantumContextRepeatMixer.java`, QCRM):** Dense Markov orders, sparse hashed k=12/15, spaced-seed and codon-phase experts, bidirectional repeat heads, and a logistic plus Bayesian mixer, all in integer arithmetic (no float/double, so no cross-CPU drift). Beats the best published specialized compressor on 11 of 14 real genomes; loses by 0.4-1.1% on PlFa, ScPo and C. parvum. See the report artifact for exact numbers.
4. **1-Hop Page-Indexed Random Access & Bounded `--lean` Streaming (`HyperGraSS.java`):** Supports sub-millisecond `O(1)` subsequence and record queries (`query` CLI) directly via `FileChannel`, multi-core parallel block execution, and strictly bounded `< 2 MB` live-heap streaming (`--lean`, verified inside `-Xmx8m` on >100 MB pangenomes where GraSS crashes with `OutOfMemoryError` below `-Xmx64m` even on 18 KB inputs).

---

## Project Structure

```text
hypergrass/
├── src/hypergrass/
│   ├── HyperGraSS.java                # .hgs container engine, parallel pipeline, --lean streaming, O(1) query CLI
│   ├── FastaStreamFactorizer.java     # 6-channel lossless FASTA/RNA/DNA factorizer & selective slice reconstructor
│   ├── CanonicalSyncmerEngine.java    # BC-LZG bidirectional (+/-) k-mer + Kimura SNP-bridging referential engine
│   ├── HaplotypeDonorIndex.java       # Exact inverted variant index + 64-lane bit-sliced donor counting
│   ├── CostAwareDonorSelector.java    # Bounded pre-Deflate byte-cost model for optional size search
│   ├── ArchiveByteAccounting.java     # Read-only exact encoded-byte breakdown; no payload decoding
│   ├── QuantumContextRepeatMixer.java # Integer-only context-mixing range coder
│   ├── BayesianMarkovRangeCoder.java  # Multi-order Markov (0..9) + codon-phase + Chargaff RC 64-bit range coder
│   ├── VarIntCodec.java               # LEB128/ZigZag varints, 2-bit DNA packing, adaptive raw/Deflate aux streams
│   ├── VerifyAllEdgeCases.java        # 13-suite SHA-256 byte-exact lossless verification harness
│   ├── VerifyDonorIndex.java          # Differential donor-count tests + optional synthetic microbenchmark
│   ├── VerifyDonorArchives.java       # Cross-version archive-byte checks + optional synthetic pipeline timings
│   ├── VerifyAnchorCodecs.java        # Bounded codec selection, size non-regression, previous-decoder checks
│   ├── VerifyArchiveAccounting.java  # Format coverage, exact totals, malformed-framing checks
│   ├── VerifyDonorCosts.java          # Independent raw-column serialization oracle for donor costs
│   ├── VerifyCostAwareEncoding.java  # Saved-build size, default-byte, query, and decoder checks
│   └── BenchmarkSuite.java            # Reproducible 5-experiment scientific benchmark suite
├── benchmarks/panman/                # Optional pinned comparison runner, reports, inputs, and command logs
├── benchmarks/cost-aware/            # Byte profiles, size/CPU results, and compatibility logs
└── README.md
```

---

## Quick Start (JDK 21+)

### 1. Compile

```powershell
mkdir -Force bin
javac -d bin src/hypergrass/*.java
```

### 2. Compress, Inspect, Query (`O(1)`), and Decompress

```powershell
# High-compression parallel mode
java -cp bin hypergrass.HyperGraSS compress input.fna archive.hgs

# Optional bounded extra size search for small archives (more CPU; never retains a larger archive)
java -cp bin hypergrass.HyperGraSS compress --optimize input.fna archive_optimized.hgs

# Ultra-lean bounded-memory streaming mode (runs safely even under -Xmx8m on 100+ MB files!)
java -Xmx8m -cp bin hypergrass.HyperGraSS compress --lean input.fna archive_lean.hgs

# Inspect archive metadata in O(1)
java -cp bin hypergrass.HyperGraSS inspect archive.hgs

# Exact encoded-byte breakdown (reads the archive into memory; not an O(1) inspection)
java -cp bin hypergrass.HyperGraSS profile archive.hgs

# Sub-millisecond O(1) random-access query: extract record #1337, genomic locus [1000, 1500)
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

All four archives were **byte-for-byte identical**, and lossless reconstruction and 60 slice queries passed. The private-variant case was slower: end-to-end performance is workload-dependent, and these small synthetic timings are not a substitute for a real-corpus benchmark. All 13 existing lossless verification suites also passed.

### Reproduce the Checks

`VerifyDonorIndex` needs no external data. Its 8 deterministic suites compare 1,928,153 donor intersection counts and 10,039 selections against a sorted-merge oracle, including lazy additions, memory-budget fallback, allele distinctions, carry boundaries, signed mask lanes, and scoring ties.

For cross-version comparison, retain the **original** compiled encoder classes in `original-bin/`. On Windows, placing that directory first lets the same test harness run against the old implementation:

```powershell
# Requires an actual saved original build in original-bin/; do not copy the updated build there.
java -Xms512m -Xmx512m -cp "original-bin;bin" hypergrass.VerifyDonorArchives baseline donor-baseline --bench
java -Xms512m -Xmx512m -cp bin hypergrass.VerifyDonorArchives verify donor-baseline --bench
```

Use `:` instead of `;` in the classpath on Linux/macOS. Omit `--bench` to run archive-byte, round-trip, and query checks without repeated timing runs. `baseline` writes the four generated archives; `verify` compares against them and fails if even one archive byte differs.
