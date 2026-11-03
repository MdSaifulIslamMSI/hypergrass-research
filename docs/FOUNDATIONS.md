# Foundations: invariants, mechanisms, and falsification

This note is the contract between the compressor and the biology it claims to
exploit. Every behaviour that is described as scientific must fit one line:

> **invariant → mechanism → measured effect → falsification test**

A feature that cannot fill all four cells is either engineering (kept, but not
claimed as biology) or speculation (removed). Naming is not evidence. The
measured numbers live in [`SPECIALIZED.md`](SPECIALIZED.md),
[`BROAD.md`](BROAD.md), [`ABLATION.md`](ABLATION.md),
[`ENTROPY-FLOOR.md`](ENTROPY-FLOOR.md), and
[`SCALING.md`](SCALING.md); this note cites them rather than restating the
protocols.

## 1. Exactness is the first invariant

**Invariant.** Lossless compression is a bijection on the normalized input
bytes. Nothing downstream is meaningful if reconstruction is not byte-identical.

**Mechanism.** A mode byte selects the block codec; the anchor and footer are
self-describing; CRC-32 covers the restored stream; every test compares SHA-256
and byte length in a separate process.

**Measured effect.** `VerifyAllEdgeCases` (15 suites), `VerifyAnchorCodecs`
(9), and `VerifyArchiveAccounting` (34 suites, 3,326 checks, 1,721 malformed
archives rejected) pass on JDK 27. The broad panel round-trips every row, and
every 120.1 Mbp specialized profile that completes is verified byte-for-byte;
a profile that runs out of memory is recorded as a failure, not a pass.

**Falsification.** Any SHA-256 or length mismatch; any malformed archive that
is accepted; any "successful" restore that differs by one byte.

## 2. The signal is homology, not nucleotide composition

**Invariant.** A pangenome is redundant because genomes descend from common
ancestors and recombine: homologous segments stay collinear over long spans,
while substitutions, indels, and rearrangements are sparse relative to the
genome length.

**Mechanism.** One anchor bank represents shared sequence. Each target block is
aligned to the anchor; exact and near-exact spans become matches, and the
residual is coded in separate gap, SNP, geometry, mask, and literal channels.
Matches shorter than 22 bp are rejected, so the coder spends bits only where
collinearity is long enough to amortize a reference.

**Measured effect.** On the identical 120,108,223-byte *E. coli* sequence the
tuned profile produces 3,810,409 B (0.254 bits/base), 2.21% below JARVIS3 l7
(3,896,361 B) and 32.7% below `xz -9` (5,663,640 B). Across six unrelated
natural collections HyperGraSS is smallest on 6/6, median margin 13.8% over the
best baseline (range 4.7–29.4%). The published GraSS implementation is 6.5x
larger on *E. coli*, 7.3x on *P. aeruginosa*, and 11.7x on *M. tuberculosis*.

**Falsification.** Destroy collinearity while preserving composition: shuffle
bases within each record (or reverse the order of records and blocks) and
compress the same collection. If the margin over JARVIS3 survives, the gain is
not from homology. A second control is per-genome compression summed over
genomes; if that matches the pangenome archive, cross-genome reference adds
nothing.

**Known boundary.** The advantage is not universal. On the *E. coli* medium
diversity stratum HyperGraSS is 0.16% *larger* than JARVIS3, the first loss in
the broad panel. At high diversity the shared-anchor fraction falls and the
referential model runs out of amortizable matches.

## 3. Both strands are indexed, and the strand is recorded

**Invariant.** DNA is double-stranded. A homologous segment can align to the
anchor in either orientation, and the reverse complement of a functional
segment is the same physical homology.

**Mechanism.** At every seed position the encoder inserts the forward 16-mer
and its reverse complement into the seed table, with a strand bit in the value.
A match carries that bit, so the decoder reconstructs the target in the correct
orientation. The two orientations are indexed separately rather than collapsed
to a `min(fwd, rev)` key, because the decoder needs the orientation, not just
the equivalence class.

**Measured effect.** The mechanism is exercised by every pangenome row and by
the inverted-repeat cases in the edge suite; no reverse-complement-specific
failure has appeared in the 58 suites. Its isolated size contribution has not
been ablated.

**Falsification.** Disable reverse-complement insertion and measure archive
growth on a panel containing inverted repeats; or restore a synthetic inverted
repeat and compare the strand field against the input.

## 4. Seed stride has a provable recall guarantee

**Mathematical invariant.** With k-mer length `K = 16` and stride `s`, a
contiguous exact island of length `L` is guaranteed to contain at least one
indexed seed start iff `L >= K + s - 1`. The default `s = 4` therefore
guarantees every island of 19 bp or longer; the opt-in `s = 8` guarantees 23 bp.

**Mechanism.** Seeds are dense 16-mers sampled every `s` positions on the
anchor, both strands. The index is an open-addressed table of 32-bit keys and
32-bit values (8 bytes per slot) at a 0.7 load factor. Its measured footprint
at a saturated 9.43 Mbp anchor is 67 MB at `s = 4` and 33.6 MB at `s = 8`,
about 7.1 and 3.6 bytes per anchor base.

**Measured effect.** On the 120.1 Mbp input, 16 MiB blocks:

| Configuration | Archive | Heap | Peak RSS (c / d) | Exact |
|---|---:|---:|---:|:---:|
| stride 4, 160 MiB heap (published bounded) | 3,997,232 B | 160m | 274.4 / 216.4 MiB | yes |
| stride 8, 160 MiB heap | 3,965,049 B | 160m | 275.1 / 216.4 MiB | yes |
| stride 8, 32 MiB blocks | 3,824,166 B | 232m | 375.4 / 286.6 MiB | yes |

Stride 8 is **0.80% smaller** at 16 MiB blocks and **1.85% smaller than
JARVIS3** at 32 MiB blocks. The direction is not obvious a priori: a sparser
index can *improve* size by removing short, expensive candidate matches while
still finding every long one. It does **not** make the 16 MiB profile fit a
256 MiB cap: peak RSS is unchanged within noise because the index is only part
of the live set. The 16 MiB RSS rows above are earlier swap-enabled
observations; under the strict no-swap cap the final run OOM-killed both 16 MiB
profiles 3/3 at 256 MiB. The 32 MiB archive passes the size gate but needs at
least 225-232 MiB of decoder heap and peaks at 375.4 / 286.6 MiB at the smallest
completing heap; 240 MiB peaks at 391.0 / 298.5 MiB.

**Cost, stated honestly.** `MATCH_ACCEPT_MIN` is 22 bp, but `s = 8` only
guarantees islands of 23 bp. A 22-bp exact island can be missed. The measured
net effect on this input is favourable; the guarantee is a lower bound on
recall, not a claim that all accepted matches are found. At a 144 MiB heap the
same stride-8 archive peaks at 258.3 MiB (still over the cap) and compression
rises to 134 s; 128 MiB exhausts the heap. No tested heap between 128 and
160 MiB satisfies both RSS and wall.

**Falsification.** Build a synthetic ladder of exact islands at 19–24 bp and
measure hit rate at `s = 4` and `s = 8`; the observed hit sets must match the
`K + s - 1` bound. On real panels, if `s = 8` ever grows an archive by more than
0.5%, the stride is not free and the default stays at 4.

**Decision.** The default remains `s = 4`, so every published archive stays
byte-identical. `--seed-step 8` is opt-in and documented as a memory-for-recall
trade.

## 5. The anchor is a model, and a model must pay rent

**Invariant.** Under a symmetric substitution model with equal base priors, the
majority base at an aligned column is the maximum-likelihood ancestral state;
the polished anchor is the median genome in Hamming distance. Coding against a
median reference minimizes the expected number of differences.

**Mechanism.** One pass accumulates per-position votes from every aligned
chunk; the encoder applies the majority base and stores the polished bases.
The bounded path packs four counters per byte so the vote buffer fits the
256 MiB envelope, then rebuilds the seed index after the buffer dies.

**Measured effect.** The consensus pass costs 9% of wall time on the 40.5 MB
campy-high input. Its *size* benefit has not been isolated, so it is the one
mechanism in this note without a measured effect. That is a debt, not a
justification.

**Falsification.** Run the same input with
`-Dhypergrass.consensusAnchor=false` and compare archives. If the polished
anchor is never smaller, the pass is not paying rent and should be removed from
the default profile.

## 6. Memory is a format constraint, not a runtime detail

**Invariant.** Peak memory is set by the live anchor index plus one block, not
by the input length. For a fixed anchor, doubling the input does not double
RSS. For a fixed input, increasing pangenome novelty saturates the anchor and
does increase RSS.

**Mechanism.** The bounded encoder streams the input twice, holds one block at
a time, carries a 256 KiB gap history across block boundaries, and builds the
large seed index only after the vote buffer is dead. The decoder streams one
block at a time. A native image is a deployment mechanism that removes JVM
compiler and metadata overhead without changing the codec or archive bytes.

**Measured effect.** The first single observation of the published bounded
profile restored the 120.1 Mbp input at 253 MiB compress / 221 MiB restore,
byte-exact. That was not a robust pass. In the final three-repeat run with
`--memory-swap` equal to `--memory`, every stride-4 and stride-8 bounded attempt
was OOM-killed at the 256 MiB cgroup cap (`returncode -9`, cgroup `oom_kill=1`;
wrapper `ru_maxrss` 269.3–269.4 MiB), while JARVIS3 completed 3/3 at 3,896,361
bytes. The spread for byte-identical output is about 22 MiB (252.7 MiB in the
single swap-enabled sample versus 269.4–274.4 MiB in later runs), so only the
worst case counts; the earlier pass depended on swap headroom. A 144 MiB heap
lowers the peak to 258.3 MiB but still misses the cap and triples wall time;
128 MiB exhausts the heap. No JVM bounded profile completes in a no-swap 256 MiB
container. The 32 MiB block profile that beats JARVIS3 needs at least 225-232
MiB of decoder heap; at 232 MiB it peaks at 375.4 / 286.6 MiB, and at 240 MiB
at 391.0 / 298.5 MiB. Under a strict 256 MiB cap it is OOM-killed at 5.9 s,
before its natural peak. NMT attributes the gap to JVM native memory: at a
240 MiB heap the encoder exits with 262.6 MiB committed but peaks at 389.3 MiB
RSS, including a 109.0 MiB `malloc` transient and a 56.0 MiB compiler arena.

The invariant is about the codec's live set, not the runtime. A GraalVM 25
native image of the same sources removes that tax: with a 192 MiB heap,
`--seed-step 8`, 32 MiB blocks, and `--split-tail-block`, it completed 3/3 under
the same strict 256 MiB no-swap cap, peaking at **240.1 MiB compress / 225.5 MiB
restore** on **3,864,464 B**. That is 1.42% above the unbounded archive (inside
the 2% gate) and 0.82% below JARVIS3, with every round trip byte-exact. The
cgroup compress peak leaves only 4.9 MiB of headroom, so the native mode is
still an envelope boundary rather than a comfortable default. A 2.08 GiB scale
probe passes because 55 copies of one sequence never saturate the bank.

**Known failure.** The JVM bounded profile fails the 256 MiB gate on 0/7 natural
collections: six exceed the 2% size gate, four exceed 256 MiB RSS, and
*P. aeruginosa* exhausts the heap before producing an archive. The native image
was not run across that panel, so no general native claim follows. The JVM
mechanism is the saturated seed index, not the input size: at stride 4 the
9.43 Mbp anchor index is 67 MB, leaving too little heap for factorization.
Halving the index with stride 8 is the direct response, and is why that knob
exists.

**Falsification.** Any input where peak RSS grows with input bytes at fixed
anchor saturation; any round trip that OOMs at the documented heap; any archive
that cannot be restored at the same heap that encoded it.

## 7. Codec choice is measured in bits, not asserted

**Invariant.** For a page of anchor bases, the smaller of the candidate codecs
is the better description under a two-part code; the choice itself is part of
the archive, so the decoder never repeats the trial.

**Mechanism.** The encoder runs the context-mixing and packed/Deflate paths on
an anchor page, keeps the smaller payload, and records the codec in the header.
`--anchor-cm-pages` applies the same idea to 64 KiB pages, trading a small size
loss for a large wall-time gain.

**Measured effect.** The whole-anchor context-mixing page costs 4.94 s versus
0.25 s for the fast 2-bit path on campy-high, but wins enough bits to be the
published default. Paged context mixing costs +0.489% on the 120.1 Mbp input
(4,636,907 → 4,659,584 B), cuts the campy-high header from 5.19 s to 2.38 s
(2.18x) and total wall by 1.32x; restore falls from 66.3–71.4 s to 40.3–53.3 s
at unchanged RSS (272–273 MiB). The decision is
deterministic: three paged runs produced identical archive SHA-256.

**Falsification.** If the chosen codec is never the smaller one on a panel, the
trial is dead weight; if paged mode ever changes restored bytes, the page
boundary is wrong.

## Knob-selection policy

Knobs are chosen from measurable properties of the collection, before the
archive is measured, and the choice is recorded in the run manifest. This is
the MDL discipline applied to configuration: pick the model with the smallest
`L(model) + L(data | model)`, where block size buys longer matches (lower data
term) at the cost of a larger working set (higher model term), and the anchor
cap buys reference coverage at the cost of index memory.

**Preflight statistics** (one streaming pass, no compression):

- records, total bases, N50, maximum record length, GC, non-ACGT fraction;
- a bottom-20,000 canonical 21-mer sketch per record, Mash distance, and
  approximate ANI (the estimator already used in [`BROAD.md`](BROAD.md));
- anchor saturation estimate: distinct shared 16-mers under the chosen cap.

**Rules, in order:**

1. **Single genome or ANI below 0.95.** Do not promise a pangenome gain. Use
   the single-genome profile; enable bounded mode only after a preflight
   round-trip at the target heap.
2. **ANI between 0.95 and 0.98.** Run a 10% pilot before committing the full
   collection. The measured boundary is not clean: the *E. coli* medium stratum
   at ANI 0.9747 lost 0.16%, while *C. jejuni* at 0.9810 won 23.9%, so ANI
   alone does not decide the margin. If the pilot margin is under 2%, report
   the collection as a tie rather than a win.
3. **Anchor cap.** Set the cap to `min(9 MiB, estimated saturation)`; use a
   4 MiB floor only when the estimate supports it. Never raise the cap above
   measured saturation: 32 MiB bought exactly 0 bytes on the *E. coli* panel
   ([`ABLATION.md`](ABLATION.md)).
4. **Seed stride.** Start at 4. Switch to 8 only when the projected stride-4
   index (7.1 bytes per anchor base) plus the block working set exceeds 85% of
   the heap. Record that the 22-bp gap in the recall guarantee was accepted.
5. **Block size.** Choose the largest block from the measured envelope, then
   preflight a round trip at the target heap. Do not extrapolate: 16 MiB blocks
   fit a 160 MiB heap, 24 MiB fails at 176 MiB, and 32 MiB needs 225-232 MiB
   for decode. The envelope is a lookup table, not a formula.
6. **Anchor paging.** Enable `--anchor-cm-pages` only when header wall is the
   bottleneck and +0.5% size is acceptable. It is not part of the published
   bounded profile.

**Escalation.** If a rule fails on two or more datasets, the rule is wrong, not
the datasets. Record the failure and change the rule before the next plan.

## Rejected levers

| Lever | Why rejected | What would change the verdict |
|---|---|---|
| Quantum-assisted coding | No known lossless-compression speedup; state preparation and QRAM cost dominate; the bottleneck here is memory bandwidth and serial passes. `QuantumContextRepeatMixer` is a name, not physics. | A reproducible lossless quantum-assisted encoder with a measured end-to-end wall or RSS win on this input. |
| Larger anchor bank alone | 32 MiB produced a byte-identical archive to 16 MiB: selection saturates at 9.43 Mbp. | A collection where saturation exceeds the cap, demonstrated by a strictly smaller archive. |
| More threads on few blocks | 2 threads on 2 blocks gave 1.15x and +0.79% size; the serial stages cap Amdahl. | A sharded design that removes the serial anchor build, measured on 2 CPUs. |
| 100:1 compression | The anchor alone needs 1.018 bits/base; the best measured DNA rate is 1.806, and non-anchor channels already exceed the whole budget. | An anchor coder below 1.018 bits/base on the distinct anchor content. |

## Naming debt

`CanonicalSyncmerEngine` implements **fixed-stride 16-mer seeds on both
strands**, not syncmers. There is no syncmer selection in the code path: no
s-mer minimum, no closed/open syncmer predicate. The paper draft and format
note previously called this a "canonical syncmer" model; they have been
corrected to "stride-sampled 16-mer seed index". The class name remains for
source compatibility. This matters because syncmers and fixed-stride seeds have
different density and recall guarantees; anyone reading the class name should
not infer an s-mer selection scheme that is not there. Either the code should
gain a real syncmer predicate (with density and guarantee measured), or the
class should be renamed in a deliberate, separately verified change.

`QuantumContextRepeatMixer` is a context-mixing range coder; the "quantum" is
naming only and is not a physics claim.

## Adding a feature

1. Name the invariant in one sentence. If it is only an implementation
   convenience, call it engineering and keep it out of the scientific claims.
2. State the mechanism in terms of the archive format, not the JVM.
3. Measure the effect on a named input with the published protocol; record
   archive bytes, heap, wall, and byte-exactness.
4. Write the falsification test before running the measurement.
5. If the effect is within noise, say so and leave the knob off by default.
6. If the feature changes archive bytes, make it opt-in and keep the published
   default byte-identical.
