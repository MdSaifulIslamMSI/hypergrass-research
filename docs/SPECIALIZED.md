# Specialized DNA compressors at pangenome scale

## Question

The frozen corpus compares HyperGraSS with JARVIS3 and GeCo3 on individual
genomes, where the margins are narrow. Does either specialized DNA compressor
close the gap on a 121.8 MB multi-strain pangenome, where most of the redundant
information is between genomes rather than within one genome?

## Protocol

- Input: the same 24-record *E. coli* pangenome as the scale test, converted to
  one flat sequence: headers removed, line breaks removed, letters uppercased,
  and the 2,357 non-A/C/G/T symbols (0.002%) replaced with `A`. The result is
  the identical 120,108,223-byte `input.seq` (SHA-256
  `51cb3ecb2e9f532833caa0b60345f37d0342474e9f4a239801486f9bc9219acb`) for
  every profile below.
- The normalization is deliberately conservative for HyperGraSS: the
  specialized tools receive slightly less information (no headers and no rare
  ambiguity symbols) than the multi-FASTA archive.
- Same locked Linux image (`d6ca62a74e37`), 2 CPUs, 2 GiB container limit,
  fresh process per operation.
- Baseline profiles: JARVIS3 level 7, GeCo3 levels 3 and 6, and `xz -9`.
  GeCo3 level 9 and JARVIS3 level 26 are excluded at this scale: in the frozen
  corpus their 4.6 MB *E. coli* cells reached about 2.14 GB RSS and the 180 s
  cap, so they are not feasible on 120 MB under the same 2 GiB container.
- HyperGraSS profiles: `--pangenome` with 64 MiB blocks (the tuned unbounded
  profile) and `--bounded` with 16 MiB blocks and a 192 MiB heap in a 256 MiB
  container (the bounded profile).
- Bounded-mode gates, fixed before the run: at most 256 MiB peak RSS for
  compress and restore, at most 2% archive growth, and at most 2x wall time
  against the tuned unbounded profile.
- Each decompressed sequence is verified in a separate process by SHA-256 and
  a streaming byte comparison against the original, independent of the codec.

<!-- SPECIALIZED-TABLES:START -->
| Profile | Archive | Bits/base | Compress | Decompress | Peak RSS (c / d) | Exact |
| --- | ---: | ---: | ---: | ---: | ---: | :---: |
| JARVIS3 l7 | 3,896,361 B | 0.260 | 95.4 s | 84.2 s | 286 / 286 MiB | yes |
| GeCo3 l3 | 14,389,043 B | 0.958 | 506.8 s | 488.5 s | 129 / 129 MiB | yes |
| GeCo3 l6 | 8,409,358 B | 0.560 | 534.4 s | 357.7 s | 513 / 513 MiB | yes |
| xz -9 | 5,663,640 B | 0.377 | 144.3 s | 6.6 s | 730 / 185 MiB | yes |
| HyperGraSS unbounded (64 MiB blocks) | **3,810,409 B** | **0.254** | 887.3 s | 82.4 s | 999 / 535 MiB | yes |
| HyperGraSS bounded (16 MiB blocks) | 4,131,860 B | 0.275 | 149.5 s | 73.4 s | 273 / 237 MiB | yes |
<!-- SPECIALIZED-TABLES:END -->

## Findings

HyperGraSS wins this comparison, narrowly, and pays for it in compression
time. On the identical 120.1 Mbp sequence the tuned profile produced
**3,810,409 bytes (0.254 bits/base)**, **2.21% smaller than JARVIS3 level 7**
(3,896,361 B) and 31.5x smaller than the input. Both restores are byte-exact
under independent verification. Against general-purpose compressors the margin
is larger: 32.7% below `xz -9`, 54.7% below GeCo3 level 6, and 73.5% below
GeCo3 level 3.

The trade is compression time and memory, not restore time. HyperGraSS needed
887.3 s to compress (9.3x JARVIS3's 95.4 s) and peaked at 998.9 MiB
(3.5x JARVIS3's 285.8 MiB). Restore was 82.4 s versus 84.2 s, at 534.9 MiB
versus 285.8 MiB. The earlier 4,885,851-byte tuned archive was 20.3% *larger*
than JARVIS3; the LZ77 plus range-coder gap channel described in
[`ENTROPY-FLOOR.md`](ENTROPY-FLOOR.md) removed that deficit and then some.

The bounded profile does not meet its gates. At 16 MiB blocks and a 192 MiB
heap it produces 4,131,860 B (0.275 bits/base), **8.44% larger** than the
unbounded archive, and peaks at **272.5 MiB**, 16.5 MiB over the 256 MiB cap.
It does compress 5.9x faster (149.5 s) and restores at 237.3 MiB, so the speed
gate passes comfortably and the restore side of the memory gate passes; the
archive-size gate and the compress-side memory gate fail. The 256 MiB cgroup
did not kill the process: `ru_maxrss` counts resident pages, including shared
file-backed pages that cgroup accounting does not charge to the container, so
the cgroup limit and the RSS metric are not the same quantity. The gate is
defined on measured peak RSS, and it fails.

The size regression is the variant-gap channel. Seven of the eight bounded
blocks already encode 1,130,591 B of gaps against **882,215 B** for the single
64 MiB-block archive, and the 2.7 MB remainder adds roughly 30-40 KB more.
Over those seven blocks the other channels are close to their global totals:
mask is 50,680 B *smaller*, diff is 19,173 B larger, geom 33,691 B larger,
strands 61,037 B larger, and lits 4,170 B larger. Splitting the sequence forces the factorizer to emit more
gap entries and shortens the match history the LZ coder can see.

Larger blocks fix the ratio but not the memory. Measured on the same host with
the same input: 32 MiB blocks at a 192 MiB heap and 64 MiB blocks at a 192 MiB
heap both die with `OutOfMemoryError` in `FastaStreamFactorizer.factorize`
before finishing the first block, and 16 MiB blocks at a 160 MiB heap also die.
16 MiB is the largest block that completes at a 192 MiB heap, and that
configuration misses the size gate by 6.4 percentage points. A gap model that
carries history across block boundaries is the remaining lever; I did not
build one, so the gate is open rather than proven impossible.

This does **not** overturn the frozen-corpus result. On single genomes, JARVIS3
and GeCo3 are close to HyperGraSS, and on some cells they are better. What it
overturns is the earlier interpretation of the pangenome result: the anchor
bank alone does not make HyperGraSS the smallest archive when many
near-identical genomes are stored together. The LZ77 gap coder is what makes
the difference at this scale.

JARVIS3 also has no indexed query interface. Its restore is a full 84.2 s
decompression, so it is not a competitor for the slice-query workload in
[`END-TO-END.md`](END-TO-END.md); it is a better choice for store-once,
restore-rarely archives, which is the workload where HyperGraSS's size claim
was strongest.

## Limits

- One pangenome, one run, one host.
- The sequence-only normalization removes headers, wrapping, case, and rare
  ambiguity symbols. The HyperGraSS archive includes those bytes. The
  asymmetry favors the specialized tools.
- JARVIS3 level 26 and GeCo3 level 9 are not measured at this scale, so this is
  not a comparison against either tool's maximum-compression setting.
- The specialized tools have no equivalent of the record-aware anchor bank or
  the page index; this measures their out-of-the-box sequence mode.
- The bounded profile misses two of its three gates. The failure is measured,
  not inferred: it is the 16 MiB factorization granularity, which the 256 MiB
  cap cannot raise.
- Times are indicative; the original scale test already documents large
  wall-time variance on this host.

## Reproduce

```sh
python -B benchmarks/specialized/run.py \
  --run-id specialized-scale-1 \
  --profiles jarvis3-l7,geco3-l3,geco3-l6,xz-9 \
  --timeout 7200
python -B benchmarks/specialized/run.py \
  --run-id hypergrass-lz-1 \
  --profiles hypergrass-pangenome,hypergrass-bounded-256m \
  --timeout 7200
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-1
```

The sequence extract and archives are gitignored; `report.json` is the
published evidence, and `verify.json` in each operation directory is the
independent codec-agnostic round-trip verdict.
