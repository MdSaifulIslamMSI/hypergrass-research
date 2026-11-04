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

Full 24-target leave-one-out table. AGC 3.2.5 is run at its best swept setting
(`-s 200000`); both codecs charge only the target and use the same 23-genome panel. Every
AGC archive was decoded and the target compared byte-for-byte (`verified=True` in
`work/agc-loo-masked/results.json`).

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

Ordering helps for two reasons. The first genome is the cheapest self-contained one
(1,110,848 B for 06, vs 1,192,494 for 00), and the early panel is built from central genomes:
the first ten marginals are 68–171 KB, whereas file order pays 349 KB for genome 01 and
190 KB for 02. It also puts 13/14/15 adjacent at the end, where 13 and 15 cost only 2,756 and
7,650 B. But 14 costs 228,020 B. That is the **twin penalty**, and it is structural: in
leave-one-out each twin sees the other (2,032 and 2,266 B), but a progressive order must
encode one of them before its near-duplicate exists. No permutation avoids paying one twin
against a panel that lacks its twin.

A lower bound makes the residual gap concrete. A progressive marginal can never be cheaper
than the same genome's LOO marginal (the prefix panel is a subset of the 23-genome LOO
panel), so the collection cannot beat

    cheapest first genome (1,110,848 B) + sum of LOO marginals (1,780,761 B) = 2,891,609 B.

HyperGraSS is 31.8% above that bound; Prim is 37.4% above it, naive 40.3%, cluster-first
41.2%. The roughly six-point difference between HyperGraSS and Prim is the value of coding
shared variation jointly: HyperGraSS stores the near-duplicate cluster once, while a
progressive RC-MDL archive stores one full twin and then an edit script. Closing the last
4.3% needs a collection representation (a graph or joint alignment), not a better
permutation of independent per-genome archives.

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
| AGC 3.2.5 LOO marginal | 3,539,729 B / 0.3986 bits/base |
| AGC / RC-MDL | **13.24x larger** |
| Targets where AGC is smaller | 0 / 25 |

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

AGC was swept on targets 0, 10, and 24 across the same 13 configurations used for *E. coli*
and *M. tuberculosis*. Its adaptive default (`-a`) had the smallest total marginal on the
three-target sample (541,013 B), so the full table uses it. The segment-size settings that
won on the other two species are worse here, which is itself a caution: "AGC at its best"
must be swept per corpus. For context, AGC's full 25-genome archive is 2,145,694 B
(0.2416 bits/base); that is a collection total, not a leave-one-out marginal, and is not
directly comparable to the 267,412 B LOO total.

The margin over AGC is larger than for *E. coli* (4.7x) or *M. tuberculosis* (11.3x)
because this panel contains several near-identical chromosomes. RC-MDL represents a target
as a small set of long panel copies plus its private differences; AGC's block model still
pays a larger marginal for each appended genome. The ratio is a property of this panel, not
a species-wide constant.

## Limits and threats to validity

- **Near-duplicate panel members.** The *S. aureus* panel is the extreme case: targets 10
  and 24 contribute 57% of the LOO total, and several genomes are near-identical. Report
  medians and per-genome ranges, and repeat on a less clonal panel.
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

1. Progressive collection total versus HyperGraSS's 3,810,409 B (valid decode order).
2. AGC sensitivity: the setting is selected on three representative targets per species;
   verify it on every target, and report the best-per-target oracle as an upper bound on AGC.
3. Ablations: k, stride, `--no-approx`, mismatch cap, panel size scaling (1, 2, 4, 8, 16, 23).
4. Panel robustness: a second panel per species and a less clonal *S. aureus* panel, to
   separate codec behavior from panel composition.
5. If the margin holds: indel-aware extension (banded affine alignment) and a
   panel-conditioned literal model.
