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
- The original table below was measured on the frozen JDK 21 image
  (`d6ca62a74e37`) and used one observation per profile. The JDK 27 replication
  uses `hypergrass-corpus:jdk27` (`sha256:580e58b6…01a4ae`), the same 2 CPUs and
  2 GiB container limit. Its four-profile table is also one observation per
  profile; the later repeat harness interleaves observations and reports median
  wall with worst-case RSS. A single sample is never treated as a gate pass.
- Bounded-mode gates, fixed before the run: at most 256 MiB peak RSS for
  compress and restore, at most 2% archive growth, and at most 2x wall time
  against the tuned unbounded profile.
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
  `--anchor-cm-pages` is an opt-in speed variant, not part of the published
  comparison: it encodes 64 KiB context-mixing anchor pages instead of one
  whole-anchor page. The tuned profile above keeps `--anchor-single-page` so
  the 256 MiB gate and these archives stay comparable.
- The JDK 27 replication adds `--seed-step 8` variants: 16 MiB blocks at a
  160 MiB heap, and 32 MiB blocks at a 232-240 MiB heap. Stride 8 halves the
  seed index but is opt-in; the published default remains stride 4.
- The native bounded gate compiles the same sources with GraalVM 25 into a
  static-musl image (`build-native.ps1`, `build-native.sh`) because NMT showed
  that JVM runtime overhead, not codec data, caused the 256 MiB overrun. The
  native profile adds `--split-tail-block`, so the 120 MB single-line record
  becomes four blocks (32 + 32 + 32 + 18.5 MiB) instead of three
  (32 + 32 + 56 MiB). The default build remains `--release 27`; `--release 25`
  is a documented exception forced by the newest available native-image
  release, and the native archive is verified against the JVM archive.
- GraSS (public implementation, commit
  `6d93700f5a2de8841f7cf2c20064413c80df2c4a`) runs unmodified on the same flat
  sequence through its headerless raw path. Its `bsc` backend is built from
  libbsc commit `1dc676540cb40a41562af6f2db04bd1adce97f90` (bsc 3.4.0) with
  the commands in `benchmarks/specialized/vendor/README.md`. The code fixes
  `-e2`, 25 MB BSC blocks, and LZP preprocessing; the JVM gets `-Xmx1500m`
  because the default heap in this 2 GiB container is 512 MiB and the input is
  one 120M-character line. Its native archive is `FinalBsc.bsc`, a
  BSC-compressed tar of `F1.SSG` and `F2.SSG`.
- Each decompressed sequence is verified in a separate process by SHA-256 and
  a streaming byte comparison against the original, independent of the codec.

<!-- SPECIALIZED-TABLES:START -->
First run, JDK 21 image, one observation per profile:

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

**Correction.** The 253 MiB bounded entry is one observation, not a robust gate
pass. Repeated runs of the byte-identical archive reached 274.4 MiB; see
*Repeatability and the 256 MiB gate* below.

### JDK 27 replication

The migration was checked by rebuilding the identical archives, not by assuming
that a newer JVM is neutral. The two byte identities were reproduced from a
fresh build: campy-high bounded stride 4 is `1,125,629 B` with SHA-256
`14DD7E94…E1C5D6`, and the 120.1 Mbp stride-8 bounded archive is `3,965,049 B`
with SHA-256 `9A61A029…4F33207E`. Both match the pre-migration archives.

The rows below come from two sessions, with one observation per profile.
Percentages are against JARVIS3 l7; RSS is the peak of that observation.

| Profile | Archive | vs JARVIS3 | Compress | Restore | Peak RSS (c / d) | Exact |
| --- | ---: | ---: | ---: | ---: | ---: | :---: |
| JARVIS3 l7 | 3,896,361 B | — | 62.4 s | 63.9 s | 286.0 / 285.6 MiB | yes |
| Bounded stride 8, 16 MiB blocks, 160 MiB heap | 3,965,049 B | +1.76% | 47.5 s | 20.7 s | 274.9 / 214.6 MiB | yes |
| Bounded stride 8, 32 MiB blocks, 232 MiB heap | 3,824,166 B | **-1.85%** | 52.3 s | 36.6 s | 375.4 / 286.6 MiB | yes |
| Bounded stride 8, 32 MiB blocks, paged CM | 3,849,995 B | −1.19% | 39.3 s | 14.0 s | 339.5 / 292.3 MiB | yes |
| Control: bounded stride 4, 16 MiB blocks, 160 MiB heap | 3,997,232 B | +2.59% | 57.1 s | 24.6 s | 274.4 / 216.4 MiB | yes |

The archive bytes are stable across runs; wall time on this host is not (the
same stride-8 restore took 20.7 s in one session and 66.5 s in another). Only
size and exactness are compared across sessions.

### Full-profile sweep on JDK 27

`full-all-methods-v1` ran the complete published profile set, one observation per
profile, on the same `hypergrass-corpus:jdk27` image and 2 GiB container limit.
Four profiles completed and restored byte-exactly. Six strict-cap profiles were
OOM-killed before writing an archive and contribute no size: all four 32 MiB-block
variants in the 256 MiB container (stride 4, tuned, paged CM, stride 8) at
268.8–270.3 MiB wrapper RSS, and both wide variants at the 320 MiB cap
(332.8–332.9 MiB).

| Profile | Archive | Bits/base | Compress | Restore | Peak RSS | Exact |
| --- | ---: | ---: | ---: | ---: | ---: | :---: |
| HyperGraSS unbounded (64 MiB blocks) | **3,810,409 B** | **0.254** | 753.8 s | 70.9 s | 581 MiB | yes |
| HyperGraSS bounded native (256 MiB limit) | 3,864,464 B | 0.257 | 158.6 s | 75.3 s | 240 MiB | yes |
| HyperGraSS bounded stride 8 (16 MiB blocks, 144 MiB heap) | 3,965,049 B | 0.264 | 129.8 s | 65.4 s | 262 MiB | yes |
| GraSS (public implementation, `-e2`) | 24,777,548 B | 1.650 | 166.9 s | 116.8 s | 1,322 MiB | yes |

The archive bytes reproduce the earlier sessions; the single-observation wall
times do not (the native profile took 158.6 s here against a 54.8 s median in
`native-gate-3`). Only size and exactness are compared across sessions. This
sweep also confirms the strict-cap result: no JVM profile fits 256 MiB, while
the native build completes at 240 MiB.

### Codec ablation on a 10 MB slice

`codec-ablation-10m-v2` ran a 10,000,000-base prefix through the same harness to
separate the payload models from the anchor and aux layers. Two diagnostic JVM
properties make the ablation possible without changing defaults:
`-Dhypergrass.dnaMode=1|2|3` forces QCRM, BC-LZG, or AC-HMC, and
`-Dhypergrass.noCompact=true` bypasses the ≤32 MiB single-genome path so the
anchor pipeline runs. Neither property is part of the archive format.

| Profile | Archive | Bits/base | Compress | Restore | Peak RSS | Exact |
| --- | ---: | ---: | ---: | ---: | ---: | :---: |
| QCRM (lean) | 2,406,156 B | 1.925 | 13.5 s | 13.4 s | 88 MiB | yes |
| Auto (all codecs) | **1,450,787 B** | **1.161** | 38.0 s | 29.8 s | 370 MiB | yes |
| BC-LZG forced | 1,454,172 B | 1.163 | 35.3 s | 50.1 s | 224 MiB | yes |
| AC-HMC forced | 1,450,787 B | 1.161 | 42.7 s | 34.0 s | 299 MiB | yes |
| Anchor CM pages | 1,463,055 B | 1.170 | 19.5 s | 18.3 s | 269 MiB | yes |
| Aux without LZ77 | 1,450,787 B | 1.161 | 37.4 s | 31.6 s | 219 MiB | yes |
| Aux without planes | 1,451,603 B | 1.161 | 33.4 s | 31.8 s | 227 MiB | yes |
| Aux basic (no optional aux) | 1,459,406 B | 1.168 | 34.5 s | 31.1 s | 212 MiB | yes |
| GraSS | 2,341,794 B | 1.873 | 14.0 s | 14.4 s | 294 MiB | yes |

Forced AC-HMC reproduces the auto archive byte-for-byte, so AC-HMC is selected
for every block that matters on this slice. BC-LZG is 3,385 B behind; QCRM is
955,369 B behind and 2.75% larger than GraSS on its own. The HyperGraSS margin
therefore comes from the anchor bank plus AC-HMC, not from QCRM in isolation.
Anchor CM pages cost 12,268 B (0.85%), which buys header speed rather than size.
Removing the LZ77 aux coder changes nothing on this slice, the plane coder is
worth 816 B, and all optional aux coders together are worth 8,619 B. Every
HyperGraSS variant except pure QCRM is 37.5–38.1% smaller than GraSS here.

### Single-genome fast path (`--single`)

The pangenome profiles spend time on work that cannot pay off for one assembly:
the anchor-bank probe, the block directory, and experts for cross-strain
variation. `--single` factorizes the input once and encodes a single QCRM
superblock. It is opt-in; every archive above still uses the default path and
keeps its published bytes. The tuned compact path applies to inputs up to
128 MiB; larger inputs fall back to the blocked path.

Tuning was measured on the exact 4,641,652-base *E. coli* K-12 sequence with the
direct QCRM harness. Each variant is the second of two observations in one
interleaved session, so compare rows within this table rather than across
sessions.

| QCRM variant | QCRM payload | Δ vs published | Encode | Decode |
| --- | ---: | ---: | ---: | ---: |
| Published (dense orders 0–10, sparse 12/15) | 1,086,677 B | — | 7.95 s | 8.03 s |
| Dense orders capped at 7 | 1,084,392 B | −2,285 B | 6.16 s | 6.13 s |
| **Dense ≤7 + no sparse hash tables** | **1,084,133 B** | **−2,544 B** | **4.19 s** | **4.22 s** |
| Dense ≤7 + no sparse + 4 heads + no verify | 1,084,864 B | −1,813 B | 3.99 s | 4.03 s |
| Dense ≤6 + all prunes | 1,085,534 B | −1,143 B | 4.42 s | 4.12 s |
| Dense ≤7, repeat heads disabled (control) | 1,096,410 B | +9,733 B | 5.21 s | 5.07 s |

These are QCRM payload bytes. The compact container adds 39 bytes of header and
metadata, so the tuned *E. coli* archive is 1,084,172 B.

Two structural results. High-order dense tables (orders 9–10) are worse than
useless on a single genome: dropping them saves 2,285 B and 23% of the time.
The two sparse hash experts (k=12, k=15) cost roughly 45% of the remaining
model time and 16 MiB per model while making the archive 2,544 B *larger*;
removing both is the single largest win in the table. Repeat heads are not
optional: disabling them costs 9,733 B. The chosen profile keeps six heads and
exact verification; a variant with four heads and exact verification skipped
measured about 0.2 s faster in a separate batch but gave back 731 B, so it was
rejected.

Allocation dominated the original profile. JFR showed `RepeatHead.predict`
allocating a four-int array per base per active head — 71.7% of all allocation.
Reusing one scratch array per head cut the harness's used heap from about
226 MiB to 83 MiB without changing archive bytes. Dropping the sparse tables
removes a further 16 MiB per model, so the serial tuned archive completes under
a small heap and a segmented run stays proportional to the thread count.

`--threads N` splits the DNA channel into N independently coded QCRM segments.
Each segment resets the model, so bytes grow; the trade is explicit, measured
on the same genome with three interleaved observations per point (medians):

| Segments | Archive | Δ vs serial | Encode | Speedup | Peak RSS |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | **1,084,172 B** | — | 17.84 s | 1.00x | 90 MiB |
| 2 | 1,090,807 B | +0.61% | 10.58 s | 1.69x | 117 MiB |
| 4 | 1,097,843 B | +1.26% | 7.32 s | 2.44x | 161 MiB |
| 6 | 1,102,073 B | +1.65% | 7.24 s | 2.46x | 206 MiB |
| 8 | 1,104,694 B | +1.89% | 8.69 s | 2.05x | 247 MiB |

Decode scales slightly better than encode on the same archives: 1.94x at two
segments, 2.67x at four, 2.75x at six, and 2.44x at eight, with peak RSS
tracking the encode column to within a few MiB. Wall times move with machine
load, so only the ratios inside this table should be compared; the archive
bytes are reproducible across sessions.

This host is an Intel i3-1220P: two performance cores, eight efficiency cores,
and a 28 W package. The encode curve flattens after four segments because the
serial model work is no longer the bottleneck and the efficiency-core clock and
package power take over; eight segments regresses. Under a strict two-CPU
container quota the picture is cleaner: 20.5 s serial to 10.3 s at two
segments, a 1.99x speedup at 99% efficiency, with every restore byte-exact
(the t4/t8 container runs in that sweep oversubscribed the two-CPU quota and
are not scaling evidence). The segment count is a property of the archive:
decode uses the archived count, not the machine's processor count, so a
restore is reproducible on any host.

### End-to-end single-genome results

Four genomes, three interleaved repeats per profile, run in the same 2-CPU
container used for the pangenome comparisons. Every HyperGraSS profile produced
the same archive size on all three repeats and restored byte-exactly (3/3).
GraSS archives move by a few bytes because its tar container records mtimes;
those runs are marked `partial`, and every GraSS restore was also exact.

| Genome | Bases | Profile | Archive | bits/base | Compress | Peak RSS |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| *E. coli* K-12 | 4,641,652 | GraSS `-e2` | 1,141,608 B | 1.968 | 7.6 s | 278 MiB |
| | | default auto-single | 1,086,715 B | 1.873 | 28.7 s | 204 MiB |
| | | **`--single`** | **1,084,172 B** | **1.869** | 17.3 s | **88 MiB** |
| | | `--single --threads 2` | 1,090,807 B | 1.880 | 8.9 s | 121 MiB |
| *B. subtilis* 168 | 4,215,606 | GraSS `-e2` | 1,034,078 B | 1.962 | 6.8 s | 248 MiB |
| | | default auto-single | 995,744 B | 1.890 | 28.6 s | 198 MiB |
| | | **`--single`** | **994,834 B** | **1.888** | 14.6 s | **87 MiB** |
| | | `--single --threads 2` | 998,693 B | 1.895 | 8.4 s | 113 MiB |
| Yeast, 3 chromosomes | 1,360,022 | GraSS `-e2` | 333,798 B | 1.963 | 3.0 s | 178 MiB |
| | | default auto-single | 315,396 B | 1.855 | 9.2 s | 133 MiB |
| | | **`--single`** | **315,289 B** | **1.855** | 5.9 s | **80 MiB** |
| | | `--single --threads 2` | 317,353 B | 1.867 | 5.9 s | 104 MiB |
| *P. aeruginosa* NY4383 | 7,316,070 | GraSS `-e2` | 1,690,428 B | 1.848 | 9.9 s | 296 MiB |
| | | default auto-single | 1,533,558 B | 1.677 | 41.4 s | 217 MiB |
| | | **`--single`** | **1,529,536 B** | **1.673** | 20.5 s | **95 MiB** |
| | | `--single --threads 2` | 1,548,347 B | 1.693 | 15.2 s | 129 MiB |

GraSS values are medians over the three repeats; its archive varied by up to
22 B between runs. HyperGraSS sizes were identical in every repeat. The
*P. aeruginosa* row is the longest record (7,316,070 bases) extracted from
`benchmarks/broad/data/pseudo/strata/high.fasta`.

The un-flagged default path detects a low-redundancy single assembly and takes
its compact branch, but still pays the anchor-bank probe; that is the
`default auto-single` row. With no flag it is 1,086,715 B on *E. coli*.

The `--single` profile is smaller than GraSS on all four genomes, by 5.03%,
3.80%, 5.55%, and **9.52%** respectively; the margin is widest on the largest
genome. It also shrinks the un-flagged default archive and beats its time on
every input (2,543 B and 1.66x on *E. coli*; 910 B and 1.96x on *B. subtilis*;
107 B and 1.56x on yeast; 4,022 B and 2.02x on *P. aeruginosa*).

Serial `--single` compression is 1.97-2.28x slower than GraSS but needs 2.2-3.2x
less peak RSS. Two segments bring the three bacterial genomes within 1.2-1.5x
of GraSS time for 0.39-1.23% more bytes; the 1.36 Mb yeast genome is too small
for the reset cost to pay off. Every archive stores its segment count, so a
two-segment file restores reproducibly on a single-core host.

### 100 Mbp single sequence

The compact single-genome cap was raised from 64 MiB to 128 MiB so `--single`
engages on a 100,000,000-base record. Two hg38 slices were measured: the first
100,000,000 bases of chr1 (305,111 non-ACGT bases, 0.305%) and of chr14
(16,315,569 non-ACGT bases, 16.3%, a gap-heavy centromeric slice). The harness
masks non-ACGT bases to A in `input.seq`, so both tools receive the identical
100,000,000-byte sequence. One repeat each; HyperGraSS used a 2 GiB heap and
GraSS kept its 1.5 GiB heap, both in a 3 GiB container. Peak RSS columns are
compress / restore.

| Input | Profile | Archive | bits/base | Compress | Peak RSS (c / r) | Exact |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| chr1, 100 Mb | GraSS `-e2` | 22,264,066 B | 1.781 | 137.5 s | 893 / 1,553 MiB | 1/1 |
| | default (no flag) | 21,056,843 B | 1.684 | 152.9 s | 967 / 278 MiB | 1/1 |
| | **`--single`** | **20,626,877 B** | **1.650** | 355.0 s | **396 / 396 MiB** | 1/1 |
| | `--single --threads 2` | 20,635,964 B | 1.651 | 70.1 s | 586 / 536 MiB | 1/1 |
| | `--single --threads 4` | 20,649,693 B | 1.652 | 124.9 s | 617 / 575 MiB | 1/1 |
| | `--single --threads 8` | 20,677,028 B | 1.654 | 225.9 s | 697 / 584 MiB | 1/1 |
| | `--lean` | 22,788,010 B | 1.823 | 40.5 s | 94 / 94 MiB | 1/1 |
| | `--pangenome` | 22,735,705 B | 1.819 | 744.2 s | 962 / 384 MiB | 1/1 |
| chr14, 100 Mb | GraSS `-e2` | 18,538,762 B | 1.483 | 81.9 s | 960 / 1,524 MiB | 1/1 |
| | **`--single`** | **17,084,921 B** | **1.367** | 115.6 s | **380 / 389 MiB** | 1/1 |

`--single` is **7.35% smaller on chr1** and **7.84% smaller on chr14**, at
2.25x and 2.5x lower compression RSS. The pangenome profile is the wrong tool
for one record: on chr1 it is 10.22% larger and 2.10x slower than `--single`,
and 2.12% larger than GraSS, because its 9.4 Mbp anchor page has no
cross-genome redundancy to exploit. The un-flagged default is the
speed-oriented choice at this size: it selected 47 blocks with a 4.19 Mbp
anchor bank, not the compact single branch, and produced 21,056,843 B. That is
2.08% larger than `--single` but 2.32x faster, and 5.42% smaller than GraSS at
1.11x GraSS's time. Its compress RSS is 967 MiB (2.4x `--single`), but its
restore RSS is 278 MiB. Two segments are nearly free at this size: 20,635,964 B
is only 9,087 B (0.044%) above serial, while compress falls to 70.1 s (5.07x)
and restore to 71.9 s (5.04x). That is 1.96x faster than GraSS compressing and
1.46x faster restoring, with a 7.31% smaller archive and 1.52x lower compress
RSS. More segments do not pay on a 2-CPU container: `--threads 4` grows to
20,649,693 B (+0.111%) and is slower than `--threads 2` at 124.9 s compress /
126.3 s restore, because four 25 Mb model resets outweigh the parallelism.
`--threads 8` is worse again: 20,677,028 B (+0.243%), 225.9 s compress / 228.7 s
restore, and 444 s of CPU time, 3.35x the CPU of two segments. Two segments is
the sweet spot at this size and CPU budget.
`--lean` is the speed and memory extreme: 40.5 s / 94
MiB compressing and 36.2 s / 94 MiB restoring, at 22,788,010 B. It is 10.48%
larger than `--single` but 8.77x faster with 4.21x lower compress RSS. Against
GraSS it is 2.35% larger but 3.40x faster compressing and 2.91x faster
restoring, at 9.5x lower compress RSS and 16.5x lower restore RSS. The serial
speed trade is real at this
size: `--single` is 2.58x slower on chr1 and 1.41x slower on chr14. chr1 is the harder
case because 99.7% of its bases are real repeat-rich sequence, versus 83.7% for
the chr14 slice; the per-base model cost on chr1 is roughly 3x the chr14 cost.
The 100 Mb result is therefore a ratio and memory win, not a serial speed win.
`--threads N` is the speed lever, since it splits the record into independently
coded segments. At 100 Mb on this 2-CPU container, two segments is the sweet
spot; more segments add model resets and oversubscribe the CPUs. It is measured
on smaller genomes above as well.

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
1,019.3 MiB. A same-image JDK 27 rerun produced **24,777,548 B**, 44 B smaller;
the tar container embeds modification times, so a few bytes move between runs
without changing the 6.50x ratio. JARVIS3 therefore dominates GraSS on every
measured axis: 6.36x
smaller, 2.35x faster to compress, 1.97x faster to restore, and 4.4x lower
compress RSS. HyperGraSS trades 3.95x more compress time for a 6.50x smaller
archive. The gap is structural: GraSS's grammar/substitution stage plus BSC has
no cross-genome anchor bank or variant-gap channel, so the 24-strain redundancy
reaches the backend only as generic BWT/LZP matches. This does not challenge the
GraSS paper's single-genome results; it measures the published pipeline out of
the box on a pangenome, with the same input and container as the other
profiles.

The tuned bounded profile fails two of its three gates once RSS is measured
more than once. At 16 MiB blocks, a 160 MiB heap, and a 256 KiB gap history
shared across blocks it produces **3,997,232 B (0.266 bits/base)**. Its first
run peaked at **252.7 MiB** compress / **221.2 MiB** restore, which looked like
a 3.3 MiB pass. A repeat of the byte-identical archive on JDK 27 peaked at
**274.4 MiB**; the stride-8 archive measured **274.9–275.1 MiB** in three
earlier runs. The final three-repeat run with `--memory-swap=memory` was
OOM-killed 3/3 at the 256 MiB cap (269.3–269.4 MiB wrapper RSS). Peak RSS on
this host moves by about **22 MiB for identical archive bytes**, so a single
observation cannot establish a 256 MiB pass; the final evidence is a kill, not
a measured overage. Compression is 57.1 s
in the same-session control and 149.5 s in the original run, both well inside
the 2x speed gate. The archive is still **4.90% larger** than the unbounded one
(+186,823 B), so the 2% size gate fails by 110,615 B. The first bounded attempt
(4,131,860 B) and the intermediate 8-bit-vote build (3,995,014 B) both restore
byte-exactly as well; they are kept as committed evidence for the failure and
the ablation.

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

Larger blocks fix the ratio but not the memory. Stride 8 at 16 MiB blocks
produces **3,965,049 B**, 0.80% smaller than the stride-4 bounded archive and
**1.76% larger than JARVIS3**. That still misses the 2% size gate: it is
154,640 B (4.06%) larger than the unbounded archive, 78,432 B beyond the
allowance. It does not fix RSS either: the
160 MiB heap peaked at 274.9–275.1 MiB. Lowering the heap to 144 MiB brings
the same archive to **258.3 MiB**, still 2.3 MiB over the cap, and stretches
compression to 134 s; 128 MiB exhausts the heap in
`FastaStreamFactorizer.factorize`. No tested heap between 128 and 160 MiB meets
both gates.

Stride 8 at **32 MiB blocks** reaches **3,824,166 B**, which is **1.85% below
JARVIS3** and only 0.36% above the unbounded archive, so it passes the size
gate for the first time. It needs at least 225-232 MiB of decoder heap. At the
smallest completing heap (232 MiB) it peaks at **375.4 MiB compress / 286.6 MiB
restore**; a 240 MiB heap peaks at 391.0 / 298.5 MiB. An earlier 256 MiB-heap
sample recorded 339.1 / 302.7 MiB, but that is not the worst case. The paged-CM
variant archived 3,849,995 B (1.19% below JARVIS3) in a separate sample and was
not part of this sweep. The size gate is therefore reachable; the RSS cap is the
binding constraint, and the 16 MiB profile that overruns least does not beat
JARVIS3.

The JVM consequence is simple: **no JVM bounded profile completes in a strict
256 MiB container.** If the cap is relaxed just enough for the 16 MiB profile
to finish, its smallest archive is 3,965,049 B (stride 8), still 68,688 B
(1.76%) larger than JARVIS3. The only JVM profile that beats JARVIS3 is
3,824,166 B at 32 MiB blocks, which peaks at 375.4 / 286.6 MiB at its smallest
completing heap (worst observed 391.0 / 298.5 MiB). The native runtime removes
that overhead and does meet the cap; *Native runtime and the 256 MiB gate*
below records the result.

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

### Repeatability and the 256 MiB gate

The original harness recorded one observation per profile. That is enough to
compare deterministic archive bytes, but not enough to certify a memory gate:
byte-identical stride-4 archives measured 252.7 MiB in the first run and
274.4 MiB in a later one, and byte-identical stride-8 archives measured
274.9–275.1 MiB in three runs. The spread is host noise, not a code change.

The hardened harness also sets `--memory-swap` equal to `--memory`, removing
the swap headroom that earlier runs inherited from Docker's default 2x swap.
The final three-repeat run (`specialized-jdk27-repeats-1`) settles the question:
every bounded attempt was **OOM-killed**, with the cgroup reporting
`oom_kill=1` and a peak equal to its hard cap. JARVIS3 completed 3/3 at the same
3,896,361 bytes. The earlier 253 MiB "pass" therefore depended on swap and on a
single sample. A bounded mode that cannot finish in its own 256 MiB container
does not meet the gate.

| Profile (3 interleaved repeats) | Result | Archive | Attempt wall | Worst RSS | cgroup peak | OOM |
|---|---|---:|---:|---:|---:|---:|
| JARVIS3 l7 | success 3/3 | 3,896,361 B | 84.2 s median | 285.9 MiB | 294.9 MiB | 0/3 |
| Bounded stride 4, 16 MiB blocks, 160m heap | **OOM 3/3** | — (partial output only) | 116.3 s median | 269.4 MiB | 256.0 MiB | 3/3 |
| Bounded stride 8, 16 MiB blocks, 160m heap | **OOM 3/3** | — (partial output only) | 110.5 s median | 269.4 MiB | 256.0 MiB | 3/3 |
| Bounded stride 8, 32 MiB blocks, 256m heap | **OOM 3/3** | — (partial output only) | 16.9 s median | 332.7 MiB | 320.0 MiB | 3/3 |
| Bounded stride 8 + CM pages, 32 MiB blocks, 256m heap | **OOM 3/3** | — (partial output only) | 10.9 s median | 332.9 MiB | 320.0 MiB | 3/3 |

The gate reports 20 failed checks across the four candidates (status,
compress RSS, restore RSS, size ratio, and wall ratio for each) and exits 3.
The partial archives written before a kill are truncated and are not reported as
sizes; only completed round trips contribute archive bytes.

A later one-observation sweep under the same strict 256 MiB cap killed 224, 240,
and 248 MiB heaps at 5.8-6.4 s (cgroup 256.0 MiB, wrapper 247.1-247.9 MiB).
The 320 MiB cap rows above are the earlier three-repeat run; both caps kill
before the natural peak. See [`MEMORY-SWEEP.md`](MEMORY-SWEEP.md) for the
natural-heap and NMT measurements.

The harness now expands each profile into interleaved repeats
(`r0/A, r0/B, r1/A, r1/B, …`), reports median wall and worst-case RSS, and
keeps the per-repeat observations in the report. Memory comes from two
independent sources: `ru_maxrss` for the child and the container cgroup's
high-water mark (`memory.peak`, or v1 `memory.max_usage_in_bytes`). A JVM
`OutOfMemoryError` is classified separately from a cgroup OOM kill; a raw
SIGKILL or exit 137 without a cgroup `oom_kill` counter is reported as
`killed`, not `oom-killed`, and an unreadable counter is reported as unknown
rather than zero. Reports record the
git commit and dirty state, image tag and digest, Docker version, Python
version, platform, repeats, and timeouts; writes are atomic and the run
directory is never reused. A cross-process lock enforces one HyperGraSS
container at a time. The gate is evaluated on the aggregated profile and exits
3 when a bounded candidate misses it, so a single lucky sample cannot pass.

### Native runtime and the 256 MiB gate

The JVM measurements above led to a native test, not a language rewrite. The
codec, archive format, and option parser are unchanged; the same sources are
compiled ahead of time into a static-musl binary so the JIT compiler arenas,
malloc transients, and shared class metadata disappear. The opt-in streaming
decoder (`FastaStreamFactorizer.reconstructTo`) removes a third block-sized
output array, and `--split-tail-block` cuts the 120 MB single-line record into
four blocks instead of leaving a 56 MiB tail. Neither change affects the default
archive path.

`native-gate-3` interleaves three JARVIS3 l7 runs with three native runs on the
identical input and image. The native profile runs under a strict
`--memory=256m --memory-swap=256m`; JARVIS3 keeps the usual 2 GiB limit because
it cannot complete under 256 MiB.

| Profile (3 interleaved repeats) | Archive | Compress median | Worst compress RSS | Worst restore RSS | Worst cgroup peak | Exact |
|---|---:|---:|---:|---:|---:|:---:|
| JARVIS3 l7 (2 GiB limit) | 3,896,361 B | 58.3 s | 286.0 MiB | 285.9 MiB | 294.8 MiB | 3/3 |
| Native bounded (256 MiB limit) | **3,864,464 B** | **54.8 s** | **240.1 MiB** | **225.5 MiB** | **251.1 MiB** | **3/3** |

The archive is 1.42% above the unbounded 3,810,409 B (limit 3,886,617 B) and
0.82% below JARVIS3. Every fixed gate check passes on the worst repeat:

| Check | Value | Limit | Verdict |
|---|---:|---:|:---:|
| status | success / success | success | pass |
| compress peak RSS | 240.1 MiB | 256.0 MiB | pass |
| restore peak RSS | 225.5 MiB | 256.0 MiB | pass |
| size ratio vs unbounded | 1.014 | 1.020 | pass |
| wall ratio vs JARVIS3 median | 0.940 | 2.000 | pass |

All three archives are byte-identical at
`8f351911f318da52c368196a47d63dc4cf0ced3791b749e0ca1bd9d7e13c3705`; each of the
three restores was compared over all 120,108,223 bytes by the separate
verifier process. Compression used 50.2 / 51.3 / 48.0 CPU s and restore 24.6 /
25.1 / 23.1 CPU s. The cgroup compress headroom is only 4.9 MiB (251.1 of
256.0 MiB) even though `ru_maxrss` leaves 15.9 MiB, so the gate is real but
thin. The native binary is 13,698,160 B, SHA-256
`508f4d53d5d81563f408dd59f8ffdbb79ae60fd86ab4ff1eb14dfdbc3530b657`.

Default-path bytes are unchanged: recompressing the 40.5 MB campy-high
collection with the current `bin/` at 160 MiB and the original bounded flags
reproduced `1,125,629 B`, SHA-256
`14DD7E9413C10B8414E90EE831666B420D81A9A3BA2A9F8DBBC1BED149E1C5D6`, exactly.

## Cross-taxon extension

The *E. coli* comparison above is the anchor for a broader panel: six unrelated
natural collections (*M. tuberculosis*, *S. aureus*, *S. pneumoniae*,
*P. aeruginosa*, *C. jejuni*, *C. trachomatis*), each about 120 MB except
*C. trachomatis* at 68.9 MB, run with the same flat-sequence protocol and locked
image. The 64 MiB HyperGraSS profile was smaller than JARVIS3 level 7 on all
six, with a median margin of 13.8% (range 4.7–29.4%) and byte-exact round
trips. The published GraSS pipeline was 7.3x larger on *P. aeruginosa* and
11.7x larger on *M. tuberculosis*. The JVM bounded profile did not generalize:
it missed the 2% size gate on every natural collection, exceeded 256 MiB peak
RSS on four, and exhausted its heap on *P. aeruginosa*.

This does not change the *E. coli* verdict: the unbounded profile wins narrowly
there (2.21%) and by more on several unrelated collections, while the JVM
bounded mode remains an *E. coli*-calibrated prototype rather than a general
256 MiB mode. The native mode was gated only on the *E. coli* anchor; the
cross-taxon panel was not rerun with it, so no general 256 MiB claim follows.
The protocol, diversity strata, and complete tables are in
[`BROAD.md`](BROAD.md) and `benchmarks/broad/summary.md`.

## Limits

- One pangenome and one host. Archive bytes are deterministic and were
  reproduced from a fresh JDK 27 build; wall time is not stable across sessions.
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
  archives differing by 4 bytes. The 24,777,592 B figure is the JDK 21 run's
  size; the JDK 27 sweep produced 24,777,548 B.
- The JVM bounded profile misses two of its three gates. Its archive is 4.90%
  larger than the unbounded one at stride 4 (110,615 B over the 2% target;
  4.06% and 78,432 B over at stride 8), and its 256 MiB RSS gate fails in the
  final three-repeat run: all three stride-4 attempts and all three stride-8
  attempts were OOM-killed at the 256 MiB cgroup cap. The earlier 252.7 MiB
  figure was a single-sample pass with swap headroom; later repeats of the same
  bytes measured 269.4 MiB before the hard cap killed them. The 144 MiB heap
  still peaks at 258.3 MiB and triples wall time, and 128 MiB exhausts the heap.
- The 32 MiB stride-8 JVM profile passes the size gate in a swap-enabled run
  (0.36% growth, 1.85% below JARVIS3) but needs 225-232 MiB of decoder heap.
  At 232 MiB it peaks at 375.4 / 286.6 MiB; the worst observed is 391.0 / 298.5
  MiB at 240 MiB. Under a strict 256 MiB cap it is OOM-killed at 5.9 s (247.8
  MiB wrapper RSS, cgroup at cap). The earlier 320 MiB cap killed the 256 MiB
  run 3/3 at 332.7 MiB; the paged-CM variant at 332.9 MiB. Neither is a 256 MiB
  mode, and no JVM profile that meets 256 MiB beats JARVIS3.
- The native bounded profile passes every fixed gate on the *E. coli* anchor:
  **3,864,464 B** (1.42% above unbounded, 0.82% below JARVIS3), 54.8 s median
  compress, 240.1 / 225.5 MiB worst RSS, and 3/3 byte-exact round trips under
  a no-swap 256 MiB cgroup cap. Its cgroup compress peak leaves only 4.9 MiB
  of headroom, and the cross-taxon panel was not rerun natively, so the result
  is an *E. coli* and host-specific 256 MiB mode, not a universal one.
- The native build is an opt-in GraalVM 25 static-musl image. The default JVM
  build remains `--release 27`; `--release 25` is a forced exception because
  GraalVM 25 is the newest published native-image. The binary is gitignored
  (13,698,160 B, SHA-256 `508f4d53…0b657`); each run records its size and hash.
  It is x86-64 Linux only, and the native archive is byte-identical to the JVM
  path so archives remain portable between the two runtimes.
- Peak RSS has a measured ~22 MiB spread for byte-identical output. Gate
  verdicts are therefore based on the worst observation across interleaved
  repeats, with cgroup and `ru_maxrss` evidence, not on a single run.
- The JVM bounded profile's memory is set by anchor saturation, not input size.
  At the 9 MiB cap the stride-4 seed index alone is 67 MB, so the 160 MB heap
  sits on the OOM boundary: the 123 MB *E. coli* medium stratum fails in the
  factorizer, while the 120 MB high stratum completes at 272 MiB peak RSS. The
  native gate uses `--seed-step 8` (33.6 MB index) and a 192 MiB heap, which is
  why it clears the cap; the 2.08 GiB repetitive probe passes because its bank
  never saturates. See [`SCALING.md`](SCALING.md).
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
# JDK 27: three interleaved observations per profile. The gate compares the
# bounded candidates with JARVIS3 (rss cap is absolute; size/wall are ratios).
python -B benchmarks/specialized/run.py \
  --run-id specialized-jdk27-repeats-1 \
  --profiles jarvis3-l7,hypergrass-bounded-256m-tuned,hypergrass-bounded-256m-stride8,hypergrass-bounded-wide,hypergrass-bounded-wide-cm \
  --repeats 3 --timeout 3600 \
  --gate-baseline jarvis3-l7 \
  --gate-profiles hypergrass-bounded-256m-tuned,hypergrass-bounded-256m-stride8,hypergrass-bounded-wide,hypergrass-bounded-wide-cm
# Native 256 MiB gate. The binary is gitignored; the report records its hash.
powershell -NoProfile -ExecutionPolicy Bypass -File build-native.ps1
python -B benchmarks/specialized/run.py \
  --run-id native-gate-3 \
  --profiles jarvis3-l7,hypergrass-bounded-native \
  --repeats 3 --timeout 1800 \
  --gate-baseline jarvis3-l7 \
  --gate-profiles hypergrass-bounded-native \
  --gate-size-baseline-bytes 3810409
# Full published profile set plus GraSS on the identical input.
python -B benchmarks/specialized/run.py \
  --run-id full-all-methods-v1 \
  --profiles hypergrass-pangenome,hypergrass-bounded-256m,hypergrass-bounded-256m-tuned,hypergrass-bounded-256m-cm-pages,hypergrass-bounded-256m-stride8,hypergrass-bounded-144m-stride8,hypergrass-bounded-wide,hypergrass-bounded-wide-cm,hypergrass-bounded-native,grass \
  --repeats 1 --timeout 7200
# Codec ablation on a 10 MB prefix; forced modes are diagnostic only.
python -B benchmarks/specialized/run.py \
  --run-id codec-ablation-10m-v2 --slice-bases 10000000 \
  --profiles hypergrass-qcrm,hypergrass-auto,hypergrass-bclzg,hypergrass-achmc,hypergrass-anchor-cm,hypergrass-aux-no-lz,hypergrass-aux-no-plane,hypergrass-aux-basic,grass \
  --repeats 1 --timeout 3600
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-1
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-2
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-4
python -B benchmarks/specialized/summarize.py --run grass-1
python -B benchmarks/specialized/summarize.py --run full-all-methods-v1
python -B benchmarks/specialized/summarize.py --run codec-ablation-10m-v2
```

The gate's ratio baseline is explicit. The JDK 27 command asks whether a
bounded candidate is within 2% of JARVIS3 while staying under 256 MiB. The
native command instead passes `--gate-size-baseline-bytes 3810409`, so the same
candidate is measured against the tuned unbounded archive (2% = 3,886,617 B)
while the wall gate still uses the measured JARVIS3 median. Both questions
matter, and the report records which baseline was used. A failed profile exits
2; a failed gate exits 3.

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
