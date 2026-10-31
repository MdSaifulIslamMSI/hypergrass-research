# Byte accounting and bounded cost-aware encoding

This is the first measured stage of the compression roadmap: **exact encoded-byte accounting**, a **donor byte-cost model**, and **whole-archive dictionary admission**. It does not implement representative-reference selection, cost-aware match parsing, or bounded-memory global sharing. It is not a claimed algorithmic breakthrough.

## Use

```sh
javac --release 21 -Xlint:all -d bin src/hypergrass/*.java
java -cp bin hypergrass.HyperGraSS compress input.fa default.hgs
java -cp bin hypergrass.HyperGraSS compress --optimize input.fa optimized.hgs
java -cp bin hypergrass.HyperGraSS profile optimized.hgs
```

The default encoder remains byte-identical to the saved pre-change build on all twelve new synthetic regression cases and the four existing donor-archive cases. The public four-argument `encodeBlockDna` also keeps its default behavior. **Extra size search is opt-in** because the observed gains do not justify slowing default compression.

`--optimize` compares three complete representations when eligible:

1. The unchanged default archive.
2. Cost-aware donor selection with the existing global dictionary.
3. Cost-aware donor selection without that dictionary.

Only a strictly smaller **complete archive**, including its compressed footer/directory, replaces the current best. A smaller block alone is not sufficient: changing block lengths can change footer compression. Ties retain the earlier representation. There is no format migration, decoder change, external engine dependency, or additional query dependency chain.

### Bounds and cost model

Archive trials require all of: non-lean mode, raw input at most **2 MiB**, at most **16 blocks**, at most **65,536 anchor bases**, and an existing nonempty global donor dictionary. Excluded archives retain the default encoding path.

Each cost-aware Mode 3 block trial requires at most **8,192 differences** and **256 matches**. Its selector shares a **1,048,576-unit** scan budget across calls; exceeding it discards the block trial. The budget charges the scratch-query length and each evaluated query/donor length pair. It is **not** a cap on total encoder instructions: intersection counting, index construction, matching, entropy coding, and other encoding work are outside it.

The model counts pre-Deflate stream bytes: sweep mode, donor-ID/count varints, byte-rounded retained masks, added-position gap varints, and one-byte differences. Added-gap prediction resets for each list, including when variants are inherited. Distinct alleles at the same position do not match. Scratch wins equal estimated costs; otherwise the first minimum-cost donor wins ties.

These are not exact final compressed costs because the auxiliary streams are compressed together by channel. Actual candidate output sizes, and ultimately complete archive sizes, decide whether a trial survives. Larger dictionaries and longer dependency chains are not introduced. `--lean` ignores the option; lean CLI outputs with and without it matched under `-Xmx8m` on the demo.

## Measured results

### Official demo

| Mode | Complete archive bytes |
| --- | ---: |
| Saved pre-change/default HyperGraSS | 10,419 |
| HyperGraSS `--optimize` | **10,257** |
| Released PanMAN 0.1.4 | 13,816 |

The new search saves **162 bytes (1.55%)** versus default on this one 609,255-byte input. The matched PanMAN sequence-storage difference is **25.76%**, not a feature-equivalent representation comparison. PanMAN also retains alignment/tree/evolutionary information. Exact FASTA/SHA-256 reconstruction and 60 queries passed in the [pinned comparison](../panman/README.md).

The initial donor-only experiment reached 10,412 bytes. Comparing dictionary overhead against its collection-wide benefit reached 10,257. [Before](profile-before.txt) and [after](profile-after.txt) profiles show the tradeoff: dictionary bytes fall from **673 to 4**, while variant/donor streams grow. The anchor payload remains **7,051 bytes**; all categories sum exactly to their archive length.

### Four development workloads: size and CPU tradeoff

The existing deterministic `VerifyDonorArchives` corpus contains four 1,559,186-byte synthetic collections with rotations and reverse complements. One local Windows/Java 21.0.6 run used two active processors, serial compression, `-Xms512m -Xmx512m`, four warmups, and medians of seven measured runs:

| Synthetic case | Default bytes | Optimized bytes | Saved bytes | Saved-build median | Optimized median |
| --- | ---: | ---: | ---: | ---: | ---: |
| Private variants | 9,027 | 8,958 | 69 | 95.890 ms | 168.073 ms |
| Shared plus private variants | 9,126 | 9,126 | 0 | 48.601 ms | 245.392 ms |
| Shared variants plus subgroups | 9,761 | 9,756 | 5 | 110.002 ms | 152.158 ms |
| Shared variants plus few private variants | 4,635 | 4,635 | 0 | 88.322 ms | 112.848 ms |

See [default timing log](timing-before.log) and [optimized timing log](timing-after.log). These are separate JVM forks, run sequentially, not interleaved repetitions with statistical confidence intervals. Timings vary; they establish no speed superiority. The extra search costs substantially more CPU, even when it saves no bytes. Peak memory was not measured; the JVM heap setting is only a limit.

### Twelve additional synthetic regression inputs

`VerifyCostAwareEncoding` covers private/shared/dense patterns, subgroups, rotations/reverse complements, a small collection, formatting/case/ambiguities, identical and unrelated records, lean encoding, a >2 MiB collection, and >16 blocks.

- Default outputs were byte-identical to the actual saved build on **all twelve**.
- Optimized output was no larger on any case. **Only the private-pattern case improved: 5,926 to 5,908 bytes.** The other eleven were unchanged.
- Serial and parallel optimized archives were identical.
- Exact reconstruction and nine slices per archive passed for old and optimized archives in the current decoder.
- The **actual saved previous decoder** read all twelve optimized archives, including their slices.
- Large/excluded and lean paths retained byte-identical output.

Logs: [before](compatibility-before.log), [after](compatibility-after.log), [previous decoder](previous-decoder.log). These are deterministic synthetic regressions, **not an independent real-data evaluation**. The overall evidence supports a modest optional improvement, not a broad compression leadership claim.

## Reproduce the saved-build checks

The genuine pre-change Java 21 build is retained in [`baseline-bin/`](baseline-bin/): 43 class files totaling 206,226 bytes. Its per-class hashes are recorded in [`../panman/report.anchor-selection.json`](../panman/report.anchor-selection.json); that historical report names the temporary path from which the same files were copied. Do not replace this build with updated classes and call it an old-build comparison. Compile the new harness in `bin/`. From the project root on Windows:

```sh
java -Xmx512m -cp "benchmarks/cost-aware/baseline-bin;bin" hypergrass.VerifyCostAwareEncoding baseline cost-baseline
java -Xmx512m -cp bin hypergrass.VerifyCostAwareEncoding verify cost-baseline
java -Xmx512m -cp "benchmarks/cost-aware/baseline-bin;bin" hypergrass.VerifyCostAwareEncoding decode cost-baseline
```

Use `:` instead of `;` on Linux/macOS. `baseline` does not enable the new option and can run with the older options class. `verify` tests both default-byte identity and optimized size non-regression, then writes `.updated.hgs` files; `decode` reads those with the previous engine. Optional `--bench` adds four warmups and seven measured samples.

The four-workload timing commands were:

```sh
java -XX:ActiveProcessorCount=2 -Xms512m -Xmx512m -cp "benchmarks/cost-aware/baseline-bin;bin" hypergrass.VerifyDonorArchives baseline timing-before --bench
java -XX:ActiveProcessorCount=2 -Xms512m -Xmx512m -cp bin hypergrass.VerifyDonorArchives baseline timing-after --bench --optimize
```

`VerifyDonorArchives verify` intentionally still requires exact archive bytes. Use it for default/index compatibility, not for comparison against an intentionally different optimized archive.

## Accounting API and verification

```java
ArchiveByteAccounting.Report report = ArchiveByteAccounting.report(archivePath);
long total = report.totalBytes();
var categories = report.categoryBytes();
var modeCounts = report.modeCounts();
```

The `Path` API reads the archive into memory; `report(byte[])` accepts an existing array. Category and mode maps are immutable and ordered. Categories partition the encoded representation, including mode/length prefixes as documented in the class. Compressed footer/directory bytes remain one opaque category.

Accounting skips encoded streams by their lengths without inflating auxiliary channels or decoding primary/anchor data. It checks framing and bounds but **does not verify compressed contents, directory entries, CRCs, or reconstruction**. Use full decompression for integrity verification. The parser is not a replacement for decoder hardening or fuzzing the decompressor.

```sh
java -Xmx512m -cp bin hypergrass.VerifyArchiveAccounting
java -Xmx512m -cp bin hypergrass.VerifyDonorCosts
java -Xmx1g -cp bin hypergrass.VerifyAllEdgeCases
java -Xmx512m -cp bin hypergrass.VerifyAnchorCodecs
java -Xmx512m -cp bin hypergrass.VerifyDonorIndex
python -m unittest discover -s benchmarks/panman -p test_run.py -v
```

Validation passed: **34 accounting suites** (3,326 checks, 1,721 malformed-framing rejections); **nine donor-cost suites** (19,804 selections, 398,060 independently serialized candidates, 2,880 randomized cases); the twelve saved-build cases; all original 13 lossless suites; nine anchor-codec suites; eight donor-index suites; four default donor-archive comparisons; and six Python benchmark tests. Java 21 `-Xlint:all` compilation and editor diagnostics were clean. See the [validation log](validation.log); the retained baseline class hashes were also verified against the historical report.
