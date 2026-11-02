# Specialized DNA compressors at pangenome scale

## Question

The frozen corpus compares HyperGraSS with JARVIS3 and GeCo3 on individual
genomes, where the margins are narrow. Does either specialized DNA compressor
close the gap on a 121.8 MB multi-strain pangenome, where most of the redundant
information is between genomes rather than within one genome? The comparison
also includes the published GraSS implementation, the ancestor of HyperGraSS,
so the effect of the HyperGraSS redesign is measured rather than assumed.

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
  profile) and `--bounded` with 16 MiB blocks in a 256 MiB container. The
  bounded profile carries a 256 KiB gap history across block boundaries (aux
  mode 7), packs consensus votes into 4 bits per counter, and runs at a 160 MiB
  heap with SerialGC, a 32 MiB code cache, and a 64 MiB metaspace cap. The
  earlier bounded attempt used the same blocks with a 192 MiB heap but no
  cross-block history; it is kept in the table as the failure baseline.
- GraSS (public implementation, commit
  `6d93700f5a2de8841f7cf2c20064413c80df2c4a`) runs unmodified on the same flat
  sequence through its headerless raw path. Its `bsc` backend is built from
  libbsc commit `1dc676540cb40a41562af6f2db04bd1adce97f90` (bsc 3.4.0) with
  the commands in `benchmarks/specialized/vendor/README.md`. The code fixes
  `-e2`, 25 MB BSC blocks, and LZP preprocessing; the JVM gets `-Xmx1500m`
  because the default heap in this 2 GiB container is 512 MiB and the input is
  one 120M-character line. Its native archive is `FinalBsc.bsc`, a
  BSC-compressed tar of `F1.SSG` and `F2.SSG`.
- Bounded-mode gates, fixed before the run: at most 256 MiB peak RSS for
  compress and restore, at most 2% archive growth, and at most 2x wall time
  against the tuned unbounded profile.
- Each decompressed sequence is verified in a separate process by SHA-256 and
  a streaming byte comparison against the original, independent of the codec.

<!-- SPECIALIZED-TABLES:START -->
| Profile | Archive | Bits/base | Compress | Decompress | Peak RSS (c / d) | Exact |
| --- | ---: | ---: | ---: | ---: | ---: | :---: |
| JARVIS3 l7 | 3,896,361 B | 0.260 | 95.4 s | 84.2 s | 286 / 286 MiB | yes |
| GraSS (public implementation, `-e2`) | 24,777,592 B | 1.650 | 224.5 s | 165.5 s | 1,255 / 1,019 MiB | yes |
| GeCo3 l3 | 14,389,043 B | 0.958 | 506.8 s | 488.5 s | 129 / 129 MiB | yes |
| GeCo3 l6 | 8,409,358 B | 0.560 | 534.4 s | 357.7 s | 513 / 513 MiB | yes |
| xz -9 | 5,663,640 B | 0.377 | 144.3 s | 6.6 s | 730 / 185 MiB | yes |
| HyperGraSS unbounded (64 MiB blocks) | **3,810,409 B** | **0.254** | 887.3 s | 82.4 s | 999 / 535 MiB | yes |
| HyperGraSS bounded, first attempt (16 MiB blocks, 192 MiB heap) | 4,131,860 B | 0.275 | 149.5 s | 73.4 s | 273 / 237 MiB | yes |
| HyperGraSS bounded, tuned (16 MiB blocks, 160 MiB heap, shared gap) | 3,997,232 B | 0.266 | 149.5 s | 67.6 s | **253 / 221 MiB** | yes |
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

The published GraSS implementation, the ancestor of HyperGraSS, is not
competitive at this scale. Run unmodified on the identical 120,108,223-byte
sequence, it produced **24,777,592 B (1.650 bits/base, 4.85x)**, which is
**6.50x the HyperGraSS archive** and **6.36x JARVIS3 level 7**. It compressed
in 224.5 s (121.7 CPU s) at 1,254.5 MiB peak RSS and restored in 165.5 s at
1,019.3 MiB. JARVIS3 therefore dominates GraSS on every measured axis: 6.36x
smaller, 2.35x faster to compress, 1.97x faster to restore, and 4.4x lower
compress RSS. HyperGraSS trades 3.95x more compress time for a 6.50x smaller
archive. The gap is structural: GraSS's grammar/substitution stage plus BSC has
no cross-genome anchor bank or variant-gap channel, so the 24-strain redundancy
reaches the backend only as generic BWT/LZP matches. This does not challenge the
GraSS paper's single-genome results; it measures the published pipeline out of
the box on a pangenome, with the same input and container as the other
profiles.

The tuned bounded profile meets two of its three gates. At 16 MiB blocks, a
160 MiB heap, and a 256 KiB gap history shared across blocks it produces
**3,997,232 B (0.266 bits/base)** and peaks at **252.7 MiB** while compressing
and **221.2 MiB** while restoring, both under the 256 MiB cap; the
compress-side margin is only 3.3 MiB, so this is a thin pass rather than a
comfortable one. Compression
takes 149.5 s, 5.9x faster than the unbounded profile and well inside the 2x
speed gate. The archive is still **4.90% larger** than the unbounded one
(+186,823 B), so the 2% size gate fails by 110,615 B. The first bounded
attempt (4,131,860 B) and the intermediate 8-bit-vote build (3,995,014 B)
both restore byte-exactly as well; they are kept as committed evidence for the
failure and the ablation.

Two changes produced the improvement. First, the gap channel now shares a
256 KiB LZ history across block boundaries (aux mode 7): over the first seven
blocks the gap channel falls from 1,130,591 B to 993,743 B, and the archive,
still using exact 8-bit vote counters, drops by 136,846 B to 3,995,014 B.
Second, the
bounded encoder packs the four vote counters per anchor position into 4 bits,
cutting the 37.7 MB vote buffer to 18.9 MB and letting the run complete at a
160 MiB heap; the 8-bit buffer exhausts that heap. The packing is not free:
counters saturate at 15, the exact 8-bit path reaches 351 votes at one
position, and 1,891 positions end with two bases tied at the cap, so some
anchor calls differ. That costs 2,218 B (3,995,014 to 3,997,232 B). The
unbounded path keeps exact 8-bit counts, so its archive is unchanged.

The size regression is the variant-gap channel plus the block seams of the
other channels. Measured against the single 64 MiB-block archive: gap
+111,528 B, strands +61,037 B, geom +33,691 B, diff +19,173 B, lits +4,170 B,
cnt +946 B, donor -2,814 B, mask -50,680 B, modes -48 B, framing +7,602 B,
and +2,218 B from vote saturation. Mask and donor improve because smaller
blocks change which anchors and donors are selected, but the 119,017 B of extra
seam cost across diff/geom/strands/lits/cnt is structural: eight block
boundaries force the factorizer to emit more gap entries and shorten the match
history the LZ coder can see. A predictor ensemble does not remove that cost.

Larger blocks fix the ratio but not the memory. On the same host and input,
32 MiB blocks at a 224 MiB heap complete but peak at 272 MiB RSS, over the
cap, and 64 MiB blocks at a 192 MiB heap exhaust the heap before the first
block completes. 16 MiB is the largest block size that meets the RSS gate. The
remaining 110,615 B to the 2% target would need a longer cross-block gap
history or a better seam model; the plausible levers are estimated at tens of
KB, so the gate is probably unreachable at this block size, but it is not
proven impossible.

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
- GraSS was given `-Xmx1500m`; the container default is 512 MiB. The code has
  no compression-level switch, so `-e2`, 25 MB BSC blocks, and LZP are fixed.
  It was built with bsc 3.4.0 (libbsc commit
  `1dc676540cb40a41562af6f2db04bd1adce97f90`), which is newer than the libbsc
  revision used in the 2025 paper; archive size may differ slightly from the
  published figures. The GraSS source was not modified. Its tar container
  embeds the F1/F2 modification times, so the archive bytes and size vary by a
  few bytes between runs; two 2 MB reruns produced identical F1/F2 streams and
  archives differing by 4 bytes. The 24,777,592 B figure is this run's size.
- The bounded profile misses one of its three gates. It meets the 256 MiB RSS
  cap (252.7 MiB compress, 221.2 MiB restore) and the 2x speed gate, but its
  archive is 4.90% larger than the unbounded one, 110,615 B over the 2% target.
  The residual is measured, not inferred: it is the 16 MiB factorization
  granularity, which the RSS cap cannot raise.
- Times are indicative; the original scale test already documents large
  wall-time variance on this host. The same 8-bit ablation archive took 51.4 s
  in the lz-2 run and 119.9 s in a later local rerun with byte-identical
  output; the host reports a 1.5 GHz nominal clock, so cross-session wall and
  CPU times move by more than 2x. Archive sizes and within-run ratios are the
  stable comparisons.

## Reproduce

```sh
powershell -NoProfile -ExecutionPolicy Bypass \
  -File benchmarks/specialized/vendor/fetch_vendor.ps1
python -B benchmarks/specialized/run.py \
  --run-id specialized-scale-1 \
  --profiles jarvis3-l7,geco3-l3,geco3-l6,xz-9 \
  --timeout 7200
python -B benchmarks/specialized/run.py \
  --run-id grass-1 \
  --profiles grass \
  --timeout 7200
python -B benchmarks/specialized/run.py \
  --run-id hypergrass-lz-1 \
  --profiles hypergrass-pangenome \
  --timeout 7200
# The 192 MiB profile now selects the exact 8-bit vote counters for the ablation.
python -B benchmarks/specialized/run.py \
  --run-id hypergrass-lz-2 \
  --profiles hypergrass-bounded-256m \
  --timeout 3600
python -B benchmarks/specialized/run.py \
  --run-id hypergrass-lz-4 \
  --profiles hypergrass-bounded-256m-tuned \
  --timeout 3600
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-1
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-2
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-4
python -B benchmarks/specialized/summarize.py --run grass-1
```

The GraSS run uses the native `bsc` binary built by `fetch_vendor.ps1`; the
harness copies it into each operation directory because GraSS invokes
`./bsc` relative to its working directory.

The lz-1 report also carries the first bounded attempt (4,131,860 B); that
entry predates the shared gap history and the packed vote array and is kept as
failure evidence. The current `hypergrass-bounded-256m` profile selects the
exact 8-bit vote counters and reproduces the lz-2 entry.

The sequence extract and archives are gitignored; `report.json` is the
published evidence, and `verification.json` in each operation directory is the
independent codec-agnostic round-trip verdict.
