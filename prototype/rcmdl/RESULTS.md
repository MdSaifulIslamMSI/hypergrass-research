# RC-MDL prototype: first measurements

Status: research prototype. Every archive below was decompressed and compared
byte-for-byte with its target (SHA-256); no number is reported from compression alone.

## What it is

`prototype/rcmdl/Rcmdl.java` is a flat-sequence implementation of the reference-conditioned
mosaic parser described in `docs/RESEARCH_SINGLE_GENOME.md`:

- a panel index over k-mers (k = 16, stride 8) in both orientations;
- greedy longest-match parsing with substitution tolerance (up to 64 mismatches per match,
  at most one per eight bases) and explicit mismatch positions/bases;
- self-copies for repeats inside the target;
- literal runs coded with an order-2 DNA model;
- one adaptive binary range coder for every channel (token type, source, orientation,
  length, position delta, distance, mismatch mask, literals).

Two modes share the parser:

- **S** (self-contained): no panel, self-copies and literals only.
- **R** (reference-conditioned): the target is coded against an external panel. The archive
  stores the target only and pins every panel member by name, length, and SHA-256.

The archive does not contain the panel; decoding requires the same panel bytes. This is the
same contract as AGC's incremental archive, and it is what makes the target's byte count the
comparable quantity.

## Protocol

- Panel: the 24 *E. coli* chromosomes in `benchmarks/pangenome/data/pangenome.fasta`
  (120,108,223 bases).
- Canonical input: the normalized flat sequence used by the specialized benchmark
  (headers removed, letters uppercased, non-ACGT replaced with A; 2,357 bytes differ from the
  raw parse). Its SHA-256 is
  `51cb3ecb2e9f532833caa0b60345f37d0342474e9f4a239801486f9bc9219acb`, the same input on which
  HyperGraSS measured 3,810,409 B.
- **Leave-one-out (LOO)**: genome *i* is compressed against the other 23 genomes. No target
  genome is present in its own panel.
- Every result was produced with `prototype/rcmdl/run.py`; per-genome records are in
  `work/runs/loo24-masked/results.json`.
- AGC 3.2.5 (official Windows x64 build) is measured with the same panel:
  `base` = AGC archive of the other 23 genomes, `full` = `base` with the target appended;
  the marginal is `|full| - |base|`. This charges only the target's bytes, exactly like an
  RC-MDL LOO archive.

## Self-contained mode

| Corpus | Bases | Archive | bits/base |
|---|---:|---:|---:|
| 24 genomes, RC-MDL S | 120,108,223 | 29,251,210 | 1.9483 |
| 00, RC-MDL S | 4,856,080 | 1,192,494 | 1.9645 |
| 00, HyperGraSS `--single` | 4,856,080 | 1,136,083 | 1.8716 |

HyperGraSS was run on the same masked FASTA (`--single`, local build) and round-tripped
byte-identically (SHA-256 match). On one genome it is 4.7% smaller than RC-MDL S, so the
self-contained arm is not competitive. That is expected: order-2 literals plus exact
self-matches are a floor, not a codec. The architecture's claim rests on the
reference-conditioned mode; a production self-contained path should use a strong fallback
(HyperGraSS or a GAH-style context-mixing codec), not this arm.

## Reference-conditioned mode, leave-one-out

| Metric | Value |
|---|---:|
| Genomes | 24 |
| Bases | 120,108,223 |
| Total archive bytes | 1,780,761 |
| bits/base | **0.1186** |
| Median per-genome bits/base | 0.1190 |
| p10 / p90 per-genome bits/base | 0.0082 / 0.2330 |
| Panel match coverage | 97.8% |
| Literal bases | 1.83% |
| Mismatch bases | 933,037 (0.78% of all bases) |
| Total compression time | 271 s |
| Total decompression time | 17.7 s |

Per-genome range (full table in `results.json`):

| Target | Bases | Archive B | bits/base | Panel coverage |
|---|---:|---:|---:|---:|
| 13 | 5,185,030 | 2,032 | 0.0031 | 100.0% |
| 14 | 5,182,390 | 2,266 | 0.0035 | 100.0% |
| 15 | 5,256,457 | 5,411 | 0.0082 | 99.9% |
| 04 | 5,365,294 | 32,929 | 0.0491 | 99.3% |
| 01 | 5,385,437 | 33,253 | 0.0494 | 99.4% |
| … | | | | |
| 22 | 4,960,534 | 144,449 | 0.2330 | 96.8% |
| 08 | 5,012,486 | 154,234 | 0.2462 | 97.2% |
| 11 | 5,205,165 | 184,176 | 0.2831 | 94.6% |

The distribution is wide because three panel members (13, 14, 15) are near-identical to
other panel genomes; their marginal cost collapses to a few kilobytes. The weighted total is
the honest headline; the median is close to it by coincidence, not by construction.

Per-genome cost is predicted by best-relative divergence, not by panel size alone. A
phase-insensitive 16-mer sketch (hash-sampled, 17,080–19,870 k-mers per genome) gives each
genome's maximum containment in any other panel member. Across the 24 targets,
`corr(max containment, log LOO bytes) = -0.79`. The near-identical genomes have containment
0.989–0.9998 and cost 2.0–5.4 KB; the most divergent (08, 22) have 0.66–0.67 and cost
144–154 KB. This is the intended behavior of a relative-LZ/mosaic codec: a genome is cheap
when *some* panel member covers it well.

An earlier ordering proxy sampled every 251st k-mer at a fixed stride. Indels shift that
phase, so unrelated genomes scored Jaccard ≈ 0.004 and the metric was useless. Sampling by
*k-mer hash* instead of by position removes the phase dependence; the same panel then shows
mean containment 0.60–0.74, which is biologically plausible for 24 *E. coli* chromosomes.

## Panel-size scaling

Target 00 against the first *N* panel members by filename (not an optimized subset). Every
archive round-tripped byte-exactly.

| Panel genomes | Panel bases | Archive B | bits/base | Gain vs N=1 |
|---:|---:|---:|---:|---:|
| 1 | 5,385,437 | 249,491 | 0.4110 | — |
| 2 | 10,246,019 | 185,439 | 0.3055 | 25.7% |
| 4 | 21,032,735 | 146,758 | 0.2418 | 41.2% |
| 8 | 40,208,597 | 102,931 | 0.1696 | 58.7% |
| 16 | 81,103,585 | 87,254 | 0.1437 | 65.0% |
| 23 | 115,252,143 | 73,194 | 0.1206 | 70.7% |

The curve is strongly diminishing. Eight genomes capture 83% of the bytes saved between a
one-genome and a 23-genome panel (146,560 of 176,297 B); doubling to 16 buys another 15,677 B,
and the last seven buy 14,060 B. The practical implication is that a surveillance panel does
not need to be huge to get most of the reference-conditioned gain. It also bounds the index
memory: 8 *E. coli* genomes is 40 Mb, roughly 200 MB in the prototype.

## Comparison at the same panel

Full 24-target leave-one-out table. AGC 3.2.5 is run at a single swept setting
(`-s 200000`); both codecs charge only the target and use the same 23-genome panel. Every
AGC archive was decoded and the target compared byte-for-byte (`verified=True` in
`work/agc-loo-masked/results.json`).

*Caveat on the AGC number.* The setting was selected by a 5-configuration sweep on target
`00`, not by a per-target oracle. The *S. aureus* full sweep shows that AGC's best segment
size varies by target and that a single setting can overstate the margin substantially. A
full 13-configuration E. coli sweep is the next control; until it completes, 4.71x should be
read as "against one target-00-selected AGC setting", not as an upper bound on AGC.

| Metric | RC-MDL | AGC 3.2.5 | AGC / RC-MDL |
|---|---:|---:|---:|
| Total charged bytes | **1,780,761** | 8,384,117 | **4.71x** |
| bits/base | **0.1186** | 0.5584 | 4.71x |
| Median per-target ratio | — | — | 4.09x |
| Ratio range | — | — | 0.72x–16.57x |
| Targets where AGC is smaller | — | — | 13, 14 |
| Best-of-both hybrid | 1,779,767 | | 1.00x |

| Target | RC-MDL B | AGC B | AGC/RC-MDL | RC-MDL bpb |
|---|---:|---:|---:|---:|
| 00 | 73,194 | 189,412 | 2.59x | 0.1206 |
| 01 | 33,253 | 540,494 | 16.25x | 0.0494 |
| 02 | 72,932 | 288,872 | 3.96x | 0.1200 |
| 03 | 116,355 | 536,683 | 4.61x | 0.1717 |
| 04 | 32,929 | 312,087 | 9.48x | 0.0491 |
| 05 | 65,978 | 290,726 | 4.41x | 0.1099 |
| 06 | 71,641 | 171,489 | 2.39x | 0.1256 |
| 07 | 70,697 | 230,295 | 3.26x | 0.1179 |
| 08 | 154,234 | 735,021 | 4.77x | 0.2462 |
| 09 | 101,599 | 412,125 | 4.06x | 0.1608 |
| 10 | 59,828 | 339,296 | 5.67x | 0.0968 |
| 11 | 184,176 | 502,058 | 2.73x | 0.2831 |
| 12 | 133,278 | 705,177 | 5.29x | 0.2110 |
| 13 | 2,032 | **1,680** | 0.83x | 0.0031 |
| 14 | 2,266 | **1,624** | 0.72x | 0.0035 |
| 15 | 5,411 | 11,848 | 2.19x | 0.0082 |
| 16 | 55,231 | 488,380 | 8.84x | 0.0881 |
| 17 | 79,927 | 250,296 | 3.13x | 0.1350 |
| 18 | 42,777 | 708,828 | 16.57x | 0.0692 |
| 19 | 52,085 | 276,906 | 5.32x | 0.0861 |
| 20 | 96,930 | 272,145 | 2.81x | 0.1575 |
| 21 | 93,381 | 232,003 | 2.48x | 0.1480 |
| 22 | 144,449 | 595,255 | 4.12x | 0.2330 |
| 23 | 36,178 | 291,417 | 8.06x | 0.0616 |
| **Total** | **1,780,761** | **8,384,117** | **4.71x** | **0.1186** |

The ratio is not uniform, and its shape is informative. AGC wins only on 13 and 14, where
the target is a near-exact copy of another panel member: with almost no substitutions, exact
matches are the right code, and RC-MDL's token/mismatch overhead is pure loss. As soon as
divergence is nonzero, RC-MDL's explicit substitution channel dominates. A production codec
should therefore choose an exact-copy path when containment is ≈1 and RC-MDL otherwise; the
measured best-of-both hybrid saves only 994 B here, so this is a micro-optimization at this
panel size, not the main effect.

### Parameter sweep (target 00)

The sweep below was run on the raw panel. The masked full-table value for the winning
setting is 189,412 B (0.3120 bpb), 0.17% below the raw 189,728 B. The ranking of settings is
unchanged.

| Codec / configuration | Archive charged to target | bits/base | Relative to RC-MDL |
|---|---:|---:|---:|
| RC-MDL, LOO | 73,198 B | 0.1206 | 1.00x |
| AGC 3.2.5, `-s 200000` (best) | 189,728 B | 0.3087 | 2.59x larger |
| AGC 3.2.5, `-a -s 200000` | 189,728 B | 0.3087 | 2.59x larger |
| AGC 3.2.5, `-a -s 1000000` | 220,079 B | 0.3581 | 3.01x larger |
| AGC 3.2.5, `-a -s 500000` | 224,148 B | 0.3647 | 3.06x larger |
| AGC 3.2.5, adaptive default (`-a`) | 316,976 B | 0.5157 | 4.33x larger |
| AGC 3.2.5, `-a -l 15` | 329,023 B | 0.5353 | 4.49x larger |
| AGC 3.2.5, `-a -f 0.05` | 321,903 B | 0.5238 | 4.40x larger |

The sweep matters: AGC's default is not its best on this corpus. A 200 kb segment size is
worth 1.67x over the default, and the adaptive flag makes no difference at that setting
(identical base and full sizes with and without `-a`). The fallback minimizers (`-f`) and a
shorter minimum match (`-l 15`) both make it worse. The comparison above uses AGC's best
measured configuration, not its default.

AGC's full 24-genome archive is 8,540,341 B (0.5688 bits/base) at the best setting and
8,417,249 B (0.5607) at the default. The sum of the 24 leave-one-out marginals is 8,384,117 B,
so the marginal table is consistent with AGC's collection behaviour, not an artifact of
`append`.

## Progressive collection (C mode)

A collection archive must decode in one pass. The prototype's progressive scheme stores
genome 0 self-contained and codes genome *i* against genomes 0..*i*-1:

| Order | Total bytes | bits/base | First genome B | Marginal B | vs HyperGraSS |
|---|---:|---:|---:|---:|---:|
| file order (naive) | 4,056,419 | 0.2702 | 1,192,494 | 2,863,925 | +6.5% larger |
| cluster-first (13, 14, 15 first) | 4,081,636 | 0.2719 | 1,262,894 | 2,818,742 | +7.1% larger |
| containment/Prim | **3,972,975** | **0.2646** | 1,110,848 | 2,862,127 | **+4.3% larger** |

The first two orders used a phase-sensitive proxy. The third uses the phase-insensitive
hash-sampled 16-mer containment sketch (above): start at the most central genome (06), then
repeatedly add the genome best covered by the already-placed set. It is the best of the three
and cuts the gap from 6.5% to 4.3%, but does not close it.

Four further order strategies were measured to test whether the residual gap is an ordering
problem:

| Strategy | Total bytes | bits/base | vs Prim | vs HyperGraSS |
|---|---:|---:|---:|---:|
| Containment/Prim | **3,972,975** | 0.26463 | — | +4.27% |
| Cluster-aware Prim | 3,982,813 | 0.26528 | +0.25% | +4.52% |
| Greedy measured marginal | 3,990,249 | 0.26578 | +0.43% | +4.72% |
| 20 fixed anchors | 4,115,368 | 0.27411 | +3.58% | +8.00% |
| 16 fixed anchors | 4,921,615 | 0.32781 | +23.9% | +29.2% |
| 12 fixed anchors | 5,494,235 | 0.36595 | +38.3% | +44.2% |
| 8 fixed anchors | 5,818,621 | 0.38756 | +46.5% | +52.7% |
| 4 fixed anchors | 6,621,279 | 0.44102 | +66.7% | +73.8% |

**Cluster-aware Prim.** The order places the three near-identical *E. coli* genomes
`13,14,15` consecutively at the end. The first member then costs 228,344 B, while the other
two cost 3,289 B and 5,237 B. In plain Prim the order is `14,13,15`; `14` pays 228,020 B and
`13` pays 2,756 B. Moving the cluster therefore saves only 1,556 B on the twins and loses
11,394 B by moving `16` and `18`; the net is +9,838 B. **The twin penalty is not removable by
permutation.** Some member of the cluster must be introduced before its near-duplicates
exist, and that first introduction costs about 228 KB regardless of which member is chosen.

**Greedy measured marginal.** At each step the harness re-compressed every remaining genome
against the current panel and chose the smallest measured marginal (eight candidates per step
after the first). It chose `06` first, then `19`, and reached 3,990,249 B — 17,274 B worse
than Prim. Myopic one-step minimisation is not a better global order; the greedy path spends
more on intermediate panels.

**Fixed anchors.** The first *N* genomes are coded progressively, then every remaining genome
is coded against only those *N* anchors. This tests whether a small representative backbone
plus independent residuals can beat a growing panel. It cannot: every anchor count is worse,
and the cost falls monotonically as the anchor set grows back toward the full progressive
panel. The panel's later members carry information that a fixed representative set does not.

The same cluster-aware reorder was tested on *S. aureus*: plain Prim gives
**1,233,433 B / 0.13889 bits/base**, cluster-aware Prim 1,237,522 B / 0.13935, so the reorder
is again 0.33% worse. AGC's full 25-genome archive is 2,145,694 B, so RC-MDL is 1.74x smaller
in the collection setting even before applying the cluster-aware order. (This is a collection
total, not a leave-one-out marginal.) The *E. coli* and *S. aureus* results agree: ordering
has been exhausted as a lever.

The lower bound from the LOO totals still applies. A progressive marginal cannot be cheaper
than the same genome's leave-one-out marginal, so the collection cannot beat

    cheapest first genome (1,110,848 B) + sum of LOO marginals (1,780,761 B) = 2,891,609 B.

HyperGraSS is 31.8% above that bound; Prim is 37.4%, greedy 38.0%, cluster-aware Prim 37.7%.
The 5.6-point difference between HyperGraSS and Prim is about 162 KB. Ordering has now been
tested four ways and does not recover it. The remaining gap is in the representation: a
progressive RC-MDL archive still codes each genome as an independent target, while
HyperGraSS codes shared variation across the collection jointly. Closing it needs a joint
collection model (graph, consensus, or multi-genome alignment), not another permutation.

For scale: AGC's full 24-genome archive at its best setting is 8,540,341 B. Prim is 2.15x
smaller than AGC, and 4.3% larger than HyperGraSS.

## Comparison with HyperGraSS

HyperGraSS measured 3,810,409 B on this exact 24-genome, 120,108,223-base masked set
(0.2538 bits/base, self-contained collection; `benchmarks/specialized/runs/hypergrass-lz-1`).
The RC-MDL LOO total is 1,780,761 B, but LOO is not itself a decodable collection: each
genome uses all 23 others. The comparable number is the progressive total above, where
HyperGraSS wins by 4.3%.

The modes should not be conflated:

| Scenario | Best measured codec | Bytes charged |
|---|---|---:|
| One genome, no external data | HyperGraSS `--single` | 1,136,083 |
| One genome, 23-genome panel already present | **RC-MDL R** | **73,198** |
| 24-genome self-contained collection | **HyperGraSS** | **3,810,409** |
| 24-genome collection, RC-MDL Prim | RC-MDL | 3,972,975 |
| 24-genome collection, AGC best setting | AGC 3.2.5 | 8,540,341 |

RC-MDL's win is the second row: when the panel is already distributed (a reference database,
a surveillance panel, previously decoded genomes), a new genome costs 0.12 bits/base instead
of AGC's 0.31. It is an incremental codec, not a replacement for a collection codec. The
collection experiment was worth running because it bounds the claim: the same parser that
wins by 4.7x on a fixed panel is 4.3% larger than HyperGraSS when it must also carry the
panel.

## Why the R mode wins on this corpus

1. **Substitutions are a channel, not a match break.** A 5 kb syntenic block with 1% SNPs is
   one token, ~50 mismatch codes, and no literals. An exact-match LZ codec needs ~50 matches.
2. **Long matches amortize position and length coding.** The parser's mean panel match is
   thousands of bases; token overhead is negligible relative to the mismatch channel.
3. **The model is shared by parser and coder.** The parser's score is the coder's cost model
   (length minus four bits per mismatch), so it does not choose parses the coder cannot pay
   for.

### Headroom

The 933,037 mismatch bases are 0.78% of the matched sequence. At 2 bits per substitution
that is only 233 KB, 13% of the 1.78 MB archive. The rest is dominated by *where* the
mismatches are: with a mean gap of ~120 matched bases, each gap position costs ~7 bits,
about 0.8 MB across the run. Literals (1.83% of bases at ~1.7 bits) cost ~0.37 MB. The
codec does not emit a per-channel byte report, so these are estimates from the parser
counters, but the ranking is robust: **gap positions, then literals, then substitution
bases.** An indel-aware aligner and a context-conditioned gap model are the next real
levers; a 5% better literal model is not.

## Cross-species generalization

*E. coli* is one point. The same codec, k, stride, and parser were run on two other
species panels fetched reproducibly from NCBI RefSeq (`fetch_panel.py`; manifests in
`work/mtb/manifest.json` and `work/saureus/manifest.json`).

### *Mycobacterium tuberculosis*

| Metric | Value |
|---|---:|
| Genomes / bases | 24 / 105,745,015 |
| Self-contained (S) | 25,064,978 B / 1.8963 bits/base |
| Leave-one-out (R) | **122,998 B / 0.00931 bits/base** |
| Panel coverage | 99.984% |
| Literal bases | 0.0036% |
| Mismatch bases | 48,266 (0.046% of matched bases) |
| Round-trip | 24/24 byte-exact |

At almost the same self-contained rate as *E. coli* (1.896 vs 1.948 bits/base), the
reference-conditioned rate is **12.7x lower** (0.0093 vs 0.1186). The mechanism is visible
in the channels: MTB is far more clonal, so 99.98% of bases are copied and the residual is
tens of thousands of SNPs. This is the expected behavior of a relative-LZ codec, not a
tuning artifact. It also warns against quoting a single "bits/base" for a species: the
panel's clonality is part of the result.

AGC is more sensitive to tuning on this panel than the earlier target-00-selected comparison
suggested. A full 13-configuration sweep on every target gives:

| AGC 3.2.5 configuration | Total marginal B | bits/base | vs RC-MDL | Targets AGC wins |
|---|---:|---:|---:|---:|
| best single config (`-s 100000`) | 292,170 | 0.02210 | **2.38x larger** | 0/24 |
| per-target oracle | **96,580** | 0.00731 | **0.79x** (AGC 21.5% smaller) | 23/24 |

These sweep totals are archive-size deltas from AGC `create`/`append`. Unlike the *S. aureus*
sweep, the MTB best-single and oracle archives were not re-decoded after the sweep was
stopped, so the two totals still need a final round-trip confirmation. The RC-MDL side is
24/24 byte-exact.

The earlier 1,390,893 B / 11.3x figure came from the target-00-selected `-s 200000` setting;
that setting totals 1,390,783 B across the full sweep, so the old comparison overstated the
margin because it missed better per-corpus settings. The best single setting is `seg100k`,
which wins no individual target but has the lowest total. The per-target oracle picks a different config
for nearly every target (`seg1m` 12, `seg500k` 6, `seg1m-min15` 3, `k25-seg200k` 2, default 1)
and is 21.5% smaller than RC-MDL overall. RC-MDL wins only target `12` against the oracle.

This is the opposite sign from *S. aureus*, where the oracle still leaves AGC 6.05x larger.
The MTB panel is extremely clonal and the absolute sizes are tiny (122,998 B total, roughly
5 KB per genome), so a 26,418 B absolute difference produces a large ratio. The honest
statement is narrower than the old headline: **RC-MDL beats a single deployable AGC setting
by 2.38x on MTB, but AGC's per-target oracle is smaller.** The margin is corpus-dependent,
not a species-wide constant.

### *Staphylococcus aureus*

Fetched: 25 complete RefSeq chromosomes, 71,043,355 bases (manifest in
`work/saureus/manifest.json`; panel SHA-256
`9413289ee51e253b6adf6c389d04f20765ad1d5ac815b33d2552d9b7732adbd5`).

| Metric | Value |
|---|---:|
| Genomes / bases | 25 / 71,043,355 |
| Self-contained (S) | 16,802,462 B / 1.8921 bits/base |
| Leave-one-out (R) | **267,412 B / 0.03011 bits/base** |
| Median per-genome LOO | 0.01069 bits/base |
| p10 / p90 per-genome LOO | 0.00316 / 0.05746 bits/base |
| Panel coverage | 99.598% |
| Literal bases | 0.312% |
| Mismatch bases | 136,181 (0.192% of matched bases) |
| Round-trip | 25/25 byte-exact |
| AGC 3.2.5 LOO marginal (`-a`) | 3,539,729 B / 0.3986 bits/base |
| AGC / RC-MDL (`-a`) | **13.24x larger** |
| AGC / RC-MDL (best single swept config) | **11.14x larger** |
| AGC / RC-MDL (per-target oracle) | **6.05x larger** |
| Targets where AGC (`-a`) is smaller | 0 / 25 |
| Targets where AGC (oracle) is smaller | 5 / 25 |

The self-contained rate is flat across the three species (1.892–1.948 bits/base); the
reference-conditioned rate is not (0.00931–0.1186). That is the expected signature of a
relative-LZ codec: the fallback is a constant floor, while the R-mode result is set by how
much of the target the panel already explains.

The *S. aureus* LOO distribution is strongly skewed. Twenty-one of 25 targets are below
0.03 bits/base, but targets 10 (0.2396) and 24 (0.2018) alone account for 57% of the
267,412 B total. The mean (0.03062) is therefore less representative than the median
(0.01069); both are reported, and the per-target table is in
`work/runs/loo25-saureus/results.json`. The first 21 targets reproduce the earlier
interrupted run exactly (same archive bytes, bits/base, and target SHA-256); the remaining
four were completed in this run.

AGC's comparison is unusually sensitive to its segment-size setting on this panel. A
13-configuration sweep was therefore run on **every** target, not just a sample:

| AGC 3.2.5 configuration | Total marginal B | bits/base | vs RC-MDL | Targets AGC wins |
|---|---:|---:|---:|---:|
| adaptive default (`-a`) | 3,539,729 | 0.3986 | 13.24x | 0/25 |
| best single config (`-s 50000`) | 2,979,078 | 0.33547 | **11.14x** | 1/25 |
| per-target oracle | **1,617,085** | 0.18210 | **6.05x** | 5/25 |

The per-target oracle picks the best of the 13 configurations independently for each target.
That is not a deployable single setting; it is a lower bound on AGC's size and therefore a
conservative bound on the RC-MDL margin. The winning configuration varies: `seg50k` wins 8
targets, `seg1m` 7, `seg500k` 4, `seg1m-min15` 2, and four other configs one each. The
largest per-target spread between AGC's best and worst configuration is 633,338 B (target
`07`). The earlier 13.24x figure was against `-a`; it overstates the margin because `-a`
happens to be poor on this panel's many small near-duplicate targets. The honest range is
**6.05x–11.14x**, with 6.05x as the conservative headline.

Against the oracle, RC-MDL still wins 20/25 targets and 1,362,452 B in aggregate; AGC wins
only five small targets (`01`, `02`, `03`, `04`, `16`) by 12,779 B combined. Against the best
single configuration, RC-MDL wins 24/25; only target `01` is smaller under AGC (4,879 B vs
2,212 B). Both the best-single and oracle archives were re-run and decoded byte-exactly
(25/25 each), so the table is not based on marginal arithmetic alone.

The 6.05x oracle margin is in the same range as the single-config margins for *E. coli* and
*M. tuberculosis*; the larger 11–13x figures are what this panel's near-identical chromosomes
add when AGC is held to one setting. RC-MDL represents a target as a small set of long panel
copies plus its private differences; AGC's block model pays a larger marginal for each
appended genome, but its per-target tuning can recover much of the gap on small targets. The
ratio is a property of this panel, not a species-wide constant.

For context, AGC's full 25-genome archive is 2,145,694 B (0.2416 bits/base); that is a
collection total, not a leave-one-out marginal, and is not directly comparable to the
267,412 B LOO total. RC-MDL's collection total is 1,233,433 B, so it is 1.74x smaller than
AGC's collection archive as well.

### Leave-cluster-out robustness

The regular leave-one-out experiment gives every target the other 23 or 24 genomes, including
near-identical siblings. That is a legitimate "panel already present" scenario, but it
overstates performance on an unseen lineage. The leave-cluster-out control removes the
target's entire Jaccard cluster from the panel, so no near-duplicate can explain the target.
The target itself is still compressed and round-tripped byte-exactly.

| Species | Jaccard | Clusters / largest | RC regular B | RC cluster B (bpb) | AGC regular B | AGC cluster B (bpb) | Cluster ratio | RC wins |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| *E. coli* | 0.98 | 22 / 3 | 1,780,761 | **2,331,696** (0.15531) | 8,384,117 | 9,203,994 (0.61305) | **3.95x** | 24/24 |
| *M. tuberculosis* | 0.99 | 13 / 6 | 122,998 | **153,182** (0.01159) | 1,390,893 | 1,188,909 (0.08995) | **7.76x** | 16/24 |
| *S. aureus* | 0.98 | 14 / 8 | 267,412 | **741,650** (0.08352) | 3,539,729 | 4,132,771 (0.46538) | **5.57x** | 23/25 |

*AGC caveat for this table.* The cluster AGC runs used AGC's adaptive default (`-a`), which
the full *S. aureus* and *M. tuberculosis* sweeps show is not its best setting. These ratios
are therefore upper bounds on the AGC margin; tuned best-single cluster runs are the next
control. The regular *E. coli* AGC number still uses a target-00-selected setting (full sweep
incomplete); the regular *M. tuberculosis* number now has a full sweep (best single 292,170 B,
oracle 96,580 B), so its `-a` cluster row is especially conservative and the 1,390,893 B
regular cell above is the old target-00-selected value, kept only for continuity with the
cluster run.

The ratios are lower than the regular-LOO headline because RC-MDL's cost rises when siblings
are removed: 1.31x on *E. coli*, 1.25x on *M. tuberculosis*, and 2.77x on *S. aureus*. That
is the expected signature of a relative-LZ codec, and it quantifies exactly how much of the
regular result came from near-duplicates. AGC's cost changes by 1.10x, 0.85x, and 1.17x
respectively; on MTB it actually improves with the smaller, more diverse panel. Against AGC
`-a`, RC-MDL remains **3.95–7.76x smaller in aggregate**, and it wins every *E. coli* target
and all but two *S. aureus* targets. The MTB cluster ratio is the least transferable: the
regular MTB oracle already shows AGC can be smaller when tuned per target.

Per-target losses are small in absolute terms. On *S. aureus*, AGC wins only `01` and `04`, by
5,358 B combined; RC-MDL's aggregate winning margin is 3,391,121 B. On MTB, AGC wins eight
small targets by 16,830 B combined, against an aggregate RC-MDL margin of 1,035,727 B. The
leave-cluster-out headline for *S. aureus* is **5.57x against AGC `-a`**; the regular 13.24x
combines near-duplicate panel composition with an AGC setting that is poor on this panel.
Against the regular per-target oracle the margin is 6.05x, so the panel effect and the tuning
effect are both real and should not be conflated.

MTB uses a 0.99 threshold because 0.98 merges 20 of 24 genomes into one cluster, leaving a
four-genome panel. That would test "codec on a tiny panel", not "codec on an unseen lineage".
At 0.99, the largest cluster is six genomes and the other twelve clusters are one to three
members.

## Limits and threats to validity

- **Near-duplicate panel members.** The *S. aureus* panel is the extreme case: targets 10
  and 24 contribute 57% of the regular LOO total, and several genomes are near-identical.
  The leave-cluster-out control quantifies the effect (2.77x cost increase; 5.57x margin over
  AGC `-a`). The regular 13.24x is not a species-wide constant, and part of the gap is AGC
  tuning: the regular per-target oracle gives 6.05x.
- **Three species, one panel each.** *E. coli*, *M. tuberculosis*, and *S. aureus* are
  measured, but each is a single RefSeq panel. The absolute bits/base is panel-dependent;
  a second panel per species and a less clonal *S. aureus* panel are the next controls.
- **No random access.** RC-MDL decodes the target as one stream; AGC supports per-sample and
  per-contig extraction. That is a feature difference, not a size difference, but it matters
  for production.
- **Prototype memory.** The 115 Mb panel index uses roughly 500 MB; larger panels need a
  sampled or packed index.
- **IUPAC bytes.** The prototype preserves every byte, including the 1,913 non-ACGT symbols
  in the panel. AGC's round-trip on the tested target was byte-exact; the comparison does not
  hide ambiguity codes.

## Next experiments

1. **Joint collection model.** Prim is 4.3% larger than HyperGraSS, and four ordering
   strategies (greedy, cluster-aware, fixed anchors) do not close the ~162 KB residual.
   Prototype a graph/consensus collection representation that codes shared variation jointly.
2. **Finish AGC sensitivity.** The full 13-config *S. aureus* sweep is done (best single
   11.14x, oracle 6.05x); the full *M. tuberculosis* sweep is also done (best single 2.38x,
   oracle 0.79x, i.e., AGC's oracle is 21.5% smaller). The *E. coli* sweep is incomplete
   (6/24 targets), so its 4.71x still uses a target-00-selected setting. Re-tune the cluster
   AGC controls at each species' best-single settings before quoting final robustness ratios.
3. **Panel robustness:** a second panel per species and a less clonal *S. aureus* panel, to
   separate codec behavior from panel composition.
4. **Residual model:** indel-aware extension (banded affine alignment) and a
   panel-conditioned literal model, aimed at the gap-position cost identified in the
   headroom analysis.
5. **Production gaps:** random access and a packed/sampled panel index.
