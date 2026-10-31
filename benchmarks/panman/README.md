# Verified PanMAN sequence-storage comparison

## Measured result

On the official **20-record demo**, using the same record identifiers and exact sequence bytes:

| Representation | Complete archive bytes | FASTA / archive |
| --- | ---: | ---: |
| PanMAN 0.1.4 | 13,816 | 44.10x |
| HyperGraSS before anchor selection change | 10,665 | 57.13x |
| HyperGraSS default, with bounded anchor selection | 10,419 | 58.48x |
| HyperGraSS `--optimize` | **10,257** | **59.40x** |

The input FASTA is **609,255 bytes**. Default HGS is **24.59% smaller than PanMAN** on this fixture; `--optimize` is **25.76% smaller than PanMAN**, or **1.55% smaller than the current default**. The earlier anchor selection change saved 2.31% over the original HGS archive. HyperGraSS was already smaller before these changes; do not attribute the entire difference to either new selection rule.

This is a measured result on one small dataset, **not worldwide superiority, a new compression-theory breakthrough, or a replacement for PanMAN's full representation**. PanMAN also preserves alignment/tree/evolutionary information that HyperGraSS does not. These are complete archive sizes for **sequence storage**, not feature-equivalent representations.

## What changed

Previously, anchor encoding was selected solely by reuse: at least 12 input bases per anchor base selected 2-bit packing plus adaptive raw/Deflate; otherwise it selected QCRM. Reuse does not determine which codec is smaller.

For a non-lean anchor bank fitting in **one 65,536-base page**, the encoder now tries both existing codecs and retains the smaller payload. Ties retain the previous choice. The mode byte, page directory, decoder, and random-access format are unchanged. Larger anchor banks and `--lean` keep their original paths. This changes archive bytes when the alternative wins; it does not require a format migration.

On the demo's 29,857-base bank, packed/Deflate uses 7,297 bytes and QCRM uses 7,051 bytes, saving 246 bytes. This is a **size/CPU tradeoff**: the extra QCRM trial costs encoding time, and choosing QCRM instead of packed anchors can also slow decompression and slice queries. No speed or memory superiority is claimed.

Optional `--optimize` additionally compares the default archive against cost-aware donor choices with and without the global donor dictionary. Comparisons include the complete serialized footer. On this fixture it removes a 673-byte dictionary (leaving a 4-byte empty-dictionary header), at the cost of larger variant/donor streams, for a net saving of **162 bytes**. The search is bounded and opt-in because it costs substantial extra CPU for modest size gains. See [`../cost-aware/README.md`](../cost-aware/README.md).

## Verification contract

The successful **`--optimize`** run in [`report.json`](report.json) verifies the following. A separately rerun saved pre-change build produced the default result in [`report.anchor-selection.json`](report.anchor-selection.json).

- PanMAN extracts **all 20 unique identifiers**, with the same complete headers and **case-sensitive sequence bytes**. Record order, FASTA wrapping, line endings, and inter-record blank separators may differ.
- No trimming, ambiguity replacement, gap removal, case folding, reverse complementation, or circular rotation is permitted to make records match.
- HGS decompression reproduces the **entire original FASTA byte-for-byte**, including its ordering and formatting, with matching SHA-256.
- **60 HGS slice queries** cover the beginning, middle, and end of every record.
- The pinned input hashes and compiled Java class files are checked; command outputs, versions, image identity, and wall times are logged.

Separately, `VerifyAnchorCodecs` passed nine deterministic synthetic cases, including both codec choices, low/high reuse, unpadded lengths, exactly 65,536 anchor bases, a multi-page bank, mixed case/line endings, and no-anchor lean encoding. Serial and parallel archives matched. Updated archives never exceeded the saved pre-change archive sizes; multi-page and lean archives remained byte-identical. Both the updated decoder reading old archives and the **actual saved previous decoder reading new archives** passed round-trip and slice checks. The previous decoder also reconstructed the demo exactly.

The original 13 lossless suites, eight donor-index suites, and four **default-mode** cross-version donor archive checks passed. Those four default archives remained byte-identical. The additional cost-aware regression harness checked twelve synthetic inputs: default bytes were unchanged, optimized archives were no larger, and both current and saved previous decoders passed exact reconstruction and slice queries. Thirty-four byte-accounting suites, nine donor-cost suites, and six Python benchmark tests also passed.

## Reproduce

The engine still needs only **JDK 21**. This optional comparison additionally needs **Python 3** and a running **Linux Docker daemon** with the image below available. Pulling the image requires network access and roughly 140 MB of compressed layers. The runner never installs tools or pulls images itself.

From the HyperGraSS project root:

```sh
javac --release 21 -Xlint:all -d bin src/hypergrass/*.java
docker pull quay.io/biocontainers/panman@sha256:ae2e1ab9f4ef57d3bef248e5580f95a615a5aff9bb53be05cac03a0c56ea0a27
python benchmarks/panman/run.py --optimize
python -m unittest discover -s benchmarks/panman -p test_run.py -v
java -Xmx512m -cp bin hypergrass.VerifyAnchorCodecs
```

The runner downloads only missing, hash-pinned demo inputs, rejects modified cached inputs, reconstructs both formats, runs the comparisons, and exits nonzero on failure. It replaces `report.json` and the documented demo output files; each run retains a separate command-log directory. Omit `--optimize` to benchmark default compression; the report records the selected options. `--classpath` selects an already compiled build; `--timeout` bounds each command. A timed-out Docker command triggers cleanup of only its uniquely named benchmark container.

Docker gets **only this benchmark directory**, two CPUs, a 2 GiB memory limit, and no network. The Java CLI runs on the host. These resource limits **are not peak-memory measurements**. Command wall times include JVM/process/container startup, and PanMAN construction starts from supplied graph JSON/Newick inputs: generating those inputs is excluded. This is **not an end-to-end raw-sequence speed comparison**. Peak memory has not been measured.

For cross-version regression, preserve a genuine pre-change engine build in `original-bin/` before editing. On Windows:

```sh
java -Xmx512m -cp "original-bin;bin" hypergrass.VerifyAnchorCodecs baseline anchor-baseline
java -Xmx512m -cp bin hypergrass.VerifyAnchorCodecs verify anchor-baseline
java -Xmx512m -cp "original-bin;bin" hypergrass.VerifyAnchorCodecs decode anchor-baseline
```

Use `:` instead of `;` in the classpath on Linux/macOS. The new test harness must be in `bin/`; the engine classes from `original-bin/` take precedence in the first and last commands. Do not copy the updated build into `original-bin/` and call it a previous-version check.

## Pinned provenance

- Demo source: [PanMAN repository](https://github.com/TurakhiaLab/panman/tree/1c3c6204134af94c10a447fe7c6182636c19e31a/test), commit `1c3c6204134af94c10a447fe7c6182636c19e31a`.
- Comparator: released **PanMAN 0.1.4**, BioContainers image `quay.io/biocontainers/panman@sha256:ae2e1ab9f4ef57d3bef248e5580f95a615a5aff9bb53be05cac03a0c56ea0a27`.
- Paper: [PanMAN, Nature Genetics](https://www.nature.com/articles/s41588-025-02478-7).

| File | SHA-256 |
| --- | --- |
| `data/sars_20.fa` | `7a5cf9dffcd22d03ab6f50e02e2750ae9f58bacb328440cedf4b068c013cbefe` |
| `data/sars_20.json` | `6d1390e22d73b81d007d0368acd8a812bcf7c05c35470595385c93be0c605e5d` |
| `data/sars_20.nwk` | `f6fcf515751f0018988fc25673f12c86ac96b9c02b620eca055184974126d660` |
| `data/sars_20.before.hgs` | `b88cc495bebd5138418e16c7ee3a19911817d42c5bff6e07bdde88c1d833bbca` |
| `data/sars_20.anchor-selection.hgs` (default) | `82b11b26de452206edb1cbb43b42847d6b2df386d8aeda0b7703688581bab72c` |
| `data/sars_20.hgs` (`--optimize`) | `97c38c0e29f7d32535253a168400d25202a8ab5766206d8b8ba923869fe9c6e9` |

## Larger published archive: failed extraction, no comparison

A separate attempt used the published **725,748-byte `rsv_4000.panman`** from [Zenodo record 17781629](https://doi.org/10.5281/zenodo.17781629) (CC-BY-4.0). Its published MD5 and downloaded SHA-256 were verified. The same pinned PanMAN binary failed before producing FASTA:

```text
capnp/arena.c++:153: failed: Exceeded message traversal limit. See capnp::ReaderOptions.
```

See [`rsv_4000_failure.json`](rsv_4000_failure.json) and [`logs/rsv-extract.log`](logs/rsv-extract.log). **No valid HGS-versus-PanMAN comparison was completed for that dataset.** This extraction failure is not evidence of a compression win or a general PanMAN limitation. The successful runner deliberately stays scoped to the validated 20-record demo.
