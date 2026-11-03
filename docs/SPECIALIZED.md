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
  160 MiB heap, and 32 MiB blocks at a 256 MiB heap. Stride 8 halves the seed
  index but is opt-in; the published default remains stride 4.
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
| Bounded stride 8, 32 MiB blocks, 256 MiB heap | 3,824,166 B | **−1.85%** | 47.6 s | 23.5 s | 339.1 / 302.7 MiB | yes |
| Bounded stride 8, 32 MiB blocks, paged CM | 3,849,995 B | −1.19% | 39.3 s | 14.0 s | 339.5 / 292.3 MiB | yes |
| Control: bounded stride 4, 16 MiB blocks, 160 MiB heap | 3,997,232 B | +2.59% | 57.1 s | 24.6 s | 274.4 / 216.4 MiB | yes |

The archive bytes are stable across runs; wall time on this host is not (the
same stride-8 restore took 20.7 s in one session and 66.5 s in another). Only
size and exactness are compared across sessions.

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
gate for the first time. It needs a 256 MiB heap and peaks at **339.1 MiB
compress / 302.7 MiB restore**, 83 and 47 MiB over the cap. The paged-CM
variant is 3,849,995 B (1.19% below JARVIS3) at 339.5 / 292.3 MiB. The size gate
is therefore reachable; the RSS cap is the binding constraint, and the 16 MiB
profile that overruns least does not beat JARVIS3.

The consequence is simple: **no bounded profile completes in a strict 256 MiB
container, and no profile that beats JARVIS3 fits 256 MiB.** If the cap is
relaxed just enough for the 16 MiB profile to finish, its smallest archive is
3,965,049 B (stride 8), still 68,688 B (1.76%) larger than JARVIS3. The only
profile that beats JARVIS3 is 3,824,166 B at 32 MiB blocks, which needs
339.1 / 302.7 MiB. The published 110,615 B size gap is therefore not the only
failure; the sub-cap profile is behind the specialized baseline even after the
full stride-8 recovery.

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

## Cross-taxon extension

The *E. coli* comparison above is the anchor for a broader panel: six unrelated
natural collections (*M. tuberculosis*, *S. aureus*, *S. pneumoniae*,
*P. aeruginosa*, *C. jejuni*, *C. trachomatis*), each about 120 MB except
*C. trachomatis* at 68.9 MB, run with the same flat-sequence protocol and locked
image. The 64 MiB HyperGraSS profile was smaller than JARVIS3 level 7 on all
six, with a median margin of 13.8% (range 4.7–29.4%) and byte-exact round
trips. The published GraSS pipeline was 7.3x larger on *P. aeruginosa* and
11.7x larger on *M. tuberculosis*. The bounded profile did not generalize: it
missed the 2% size gate on every natural collection, exceeded 256 MiB peak RSS
on four, and exhausted its heap on *P. aeruginosa*.

This does not change the *E. coli* verdict: the unbounded profile wins narrowly
there (2.21%) and by more on several unrelated collections, while the bounded
mode remains an *E. coli*-calibrated prototype rather than a general 256 MiB
mode. The protocol, diversity strata, and complete tables are in
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
  archives differing by 4 bytes. The 24,777,592 B figure is this run's size.
- The bounded profile misses two of its three gates. Its archive is 4.90%
  larger than the unbounded one at stride 4 (110,615 B over the 2% target;
  4.06% and 78,432 B over at stride 8), and its 256 MiB RSS gate fails in the
  final three-repeat run: all three stride-4 attempts and all three stride-8
  attempts were OOM-killed at the 256 MiB cgroup cap. The earlier 252.7 MiB
  figure was a single-sample pass with swap headroom; later repeats of the same
  bytes measured 269.4 MiB before the hard cap killed them. The 144 MiB heap
  still peaks at 258.3 MiB and triples wall time, and 128 MiB exhausts the heap.
- The 32 MiB stride-8 profile passes the size gate in a swap-enabled run
  (0.36% growth, 1.85% below JARVIS3) but needs a 256 MiB heap. Under the strict
  320 MiB cap it was OOM-killed 3/3 (332.7 MiB wrapper RSS, cgroup at cap); the
  paged-CM variant was OOM-killed 3/3 at 332.9 MiB. Neither is a 256 MiB mode,
  and no profile that meets 256 MiB beats JARVIS3.
- Peak RSS has a measured ~22 MiB spread for byte-identical output. Gate
  verdicts are therefore based on the worst observation across interleaved
  repeats, with cgroup and `ru_maxrss` evidence, not on a single run.
- The bounded profile's memory is set by anchor saturation, not input size. At
  the 9 MiB cap the seed index alone is 67 MB, so the 160 MB heap sits on the
  OOM boundary: the 123 MB *E. coli* medium stratum fails in the factorizer,
  while the 120 MB high stratum completes at 272 MiB peak RSS. The 2.08 GiB
  repetitive probe passes because its bank never saturates. See
  [`SCALING.md`](SCALING.md).
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
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-1
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-2
python -B benchmarks/specialized/summarize.py --run hypergrass-lz-4
python -B benchmarks/specialized/summarize.py --run grass-1
```

The gate's ratio baseline is explicit. The published *bounded versus
unbounded* gate uses `hypergrass-pangenome` as the baseline and compares the
2% size growth with 3,810,409 B; the command above instead asks whether a
bounded candidate is within 2% of JARVIS3 while staying under 256 MiB. Both
questions matter, and the report records which baseline was used. A failed
profile exits 2; a failed gate exits 3.

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
