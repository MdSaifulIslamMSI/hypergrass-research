# Broad pangenome and diversity experiment

## Question

The *E. coli* result and the two independent *Enterobacteriaceae* replications
all sit in one family of Gram-negative bacteria. Does the HyperGraSS advantage
survive on pangenomes from unrelated phyla, and how does it change as the
genetic diversity of the collection increases?

This experiment separates two axes that were confounded in the earlier work:

1. **Taxonomic breadth** — one natural collection per species, spanning
   Actinomycetota, Bacillota, Pseudomonadota, Campylobacterota, and
   Chlamydiota, with GC content from about 30% to 66% and genome sizes from
   about 1.0 to 6.6 Mbp.
2. **Within-species diversity** — for three species, low/medium/high strata
   drawn from the same 72-genome pool, with a measured Mash distance and
   approximate ANI for every stratum.

## Protocol

### Data acquisition

- For each species, NCBI Nucleotide is queried for RefSeq complete genomes
  (`"<organism>"[Organism] AND "complete genome"[Title] AND refseq[Filter]`).
- Up to 300 candidate records are retrieved; records outside 0.8–9.0 Mbp are
  discarded to remove plasmids and incomplete assemblies.
- A uniform random sample of 72 eligible records is fetched with a fixed seed
  (`20261102`) and sorted by accession. The sample, query, per-record lengths,
  GC, ambiguity counts, sequence hashes, and file offsets are written to
  `benchmarks/broad/data/<name>/manifest.json`.
- The pool is not filtered for phylogenetic balance. It is a random sample of
  the RefSeq complete-genome catalogue returned by NCBI, and the diversity
  strata are then selected from it.
- A k-mer con-specificity QC step drops records whose median Mash distance to
  the rest of the pool is more than `max(0.05, 3 × pool median)`, whose length
  is outside 0.6–1.4× the pool median, or whose ambiguity fraction exceeds 5%.
  The relaxed fallback is documented in `diversity.json` if it is used.

### Diversity metric and strata

- Each genome is sketched as the bottom 20,000 hashes of **all** canonical
  21-mers, using a splitmix64 hash of the 2-bit k-mer code. The earlier
  position-strided sketch was rejected because homologous k-mers offset by
  indels were missed; this is recorded here so the mistake is not repeated.
- Jaccard similarity is converted to a Mash distance
  `D = -ln(2J/(1+J))/k` and an approximate ANI `1-D`.
- Genomes are clustered with average-linkage hierarchical clustering on the
  Mash distance matrix (`k = min(4, max(2, n//18))`, small clusters merged).
- The three strata are disjoint and fixed-size (`n = min(24, pool//3)`):
  - **low** — the `n` most central members of the tightest cluster;
  - **high** — a farthest-point sample of `n` genomes from the remainder;
  - **medium** — the remaining `n` genomes.
  The measured mean pairwise distance, median, range, cluster coverage, and
  approximate ANI are stored for every stratum. Monotonicity is checked and
  reported; it is not assumed from the labels.
- A **natural** collection is also written for each species: genomes are taken
  in farthest-point order from the pool medoid until the collection reaches
  about 120 MB (or the whole pool, if smaller).

### Compression protocol

- The primary comparison is the same sequence-only protocol as
  [`SPECIALIZED.md`](SPECIALIZED.md): headers removed, letters uppercased, and
  non-A/C/G/T symbols replaced with `A`. Every profile receives the identical
  flat sequence and is verified against its SHA-256 and a streaming byte
  comparison in a separate process.
- Locked Linux image `hypergrass-corpus:locked` (`d6ca62a74e37`), 2 CPUs,
  2 GiB container limit, fresh process per operation, no network.
- Profiles:
  - HyperGraSS tuned (natural collections and diversity strata):
    `--pangenome --single-thread --anchor-single-page --block-mib 64
    --anchor-mib 9`, 1536 MiB heap. The same profile is used for both panels so
    the diversity gradient is not confounded with a profile change.
  - HyperGraSS bounded: `--pangenome --bounded --single-thread
    --anchor-single-page --block-mib 16 --anchor-mib 9`, 160 MiB heap,
    256 MiB container, shared gap history and packed votes.
  - JARVIS3 level 7.
  - `xz -9`.
  - GeCo3 level 6 on three taxonomically spread natural collections.
  - GraSS (published implementation, `-e2`) on one low-diversity and one
    high-diversity natural collection.
- The 64 MiB HyperGraSS profile is the maximum-compression profile used in the
  *E. coli* specialized comparison. The bounded profile is the 256 MiB gate
  configuration. No parameter was selected on the new datasets.
- Times are indicative and host variance is large; archive sizes, bits/base,
  round-trip verdicts, and within-run ratios are the stable measurements.

### Matrix

- **Cross-taxon natural collections (6):** *Mycobacterium tuberculosis*,
  *Staphylococcus aureus*, *Streptococcus pneumoniae*, *Pseudomonas
  aeruginosa*, *Campylobacter jejuni*, *Chlamydia trachomatis*. The original
  24-strain *E. coli* collection remains the anchor for this comparison and is
  not re-fetched; its result is reported in
  [`SPECIALIZED.md`](SPECIALIZED.md).
- **Diversity gradient (3 species × 3 strata):** *E. coli*, *S. aureus*, and
  *C. jejuni*. The *E. coli* gradient uses a separate 72-genome pool fetched
  for this experiment. Each stratum has 24 genomes (`min(24, kept//3)`), drawn
  from the same pool.
- The exact run matrix, including per-run profile lists and diversity metadata,
  is `benchmarks/broad/matrix.json`.

## Results

<!-- BROAD-TABLES:START -->
### Cross-taxon natural collections

HyperGraSS produced the smallest archive on all six collections, with a median
margin of 13.8% over the best baseline (range 4.7–29.4%) and byte-exact round
trips on every row. JARVIS3 level 7 is the best baseline throughout. `xz -9` is
within 1.1% of JARVIS3 on *M. tuberculosis* and *C. jejuni* but 23–79% larger
on the other four collections. GeCo3 level 6, measured on three collections, is
3.4–4.8x larger than JARVIS3.

The margin is not monotone in the measured Mash diversity: the largest margins
are *C. jejuni* (d = 0.0190, +23.9%) and *C. trachomatis* (d = 0.0045, +29.4%),
while *S. aureus* (d = 0.0125) has the smallest (+4.7%). These collections
differ in genome size, GC content, repeat structure, and clade composition at
the same time, so the panel establishes taxonomic breadth but does not isolate
a diversity effect. The within-species gradient below controls the species and
varies only the sampled diversity.

| Dataset | Organism | Records | Mb | Mash d | ANI | HGS | JARVIS3 | xz -9 | GeCo3 l6 | Best baseline | HGS vs best | Exact |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| campy | *Campylobacter jejuni* | 72 | 120 | 0.01897 | 0.9810 | 1,535,296 | 2,017,046 | 2,019,208 | 7,673,951 | jarvis3-l7 | +23.88% | yes |
| chlamydia | *Chlamydia trachomatis* | 66 | 69 | 0.00450 | 0.9955 | 333,507 | 472,224 | 634,356 | — | jarvis3-l7 | +29.38% | yes |
| mtb | *Mycobacterium tuberculosis* | 28 | 124 | 0.00083 | 0.9992 | 1,139,535 | 1,418,155 | 1,433,772 | 6,772,719 | jarvis3-l7 | +19.65% | yes |
| pneumo | *Streptococcus pneumoniae* | 57 | 120 | 0.01223 | 0.9878 | 1,754,114 | 1,899,711 | 2,328,208 | — | jarvis3-l7 | +7.66% | yes |
| pseudo | *Pseudomonas aeruginosa* | 19 | 126 | 0.01133 | 0.9887 | 3,083,542 | 3,350,844 | 6,005,620 | — | jarvis3-l7 | +7.98% | yes |
| staph | *Staphylococcus aureus* | 43 | 121 | 0.01246 | 0.9875 | 1,807,399 | 1,895,866 | 2,741,684 | 6,427,702 | jarvis3-l7 | +4.67% | yes |

The published GraSS pipeline is 7.3x larger than HyperGraSS on *P. aeruginosa*
and 11.7x larger on *M. tuberculosis*, with exact round trips on both. That is
the same direction and rough magnitude as the 6.5x gap on the *E. coli* anchor,
so the redesign rather than a single-dataset artifact accounts for the
difference.

### Bounded-memory gate

The bounded profile does not generalize at the 256 MiB gate. It fails on all
seven natural collections: six produce an archive over the 2% size gate, and
the seventh (*P. aeruginosa*) exhausts its 160 MiB heap before producing one.
It exceeds 256 MiB peak RSS on four of the six that complete. The six new
collections are tabulated below; the *E. coli* anchor is the seventh and is
reported in [`SPECIALIZED.md`](SPECIALIZED.md).
The absolute overhead is roughly constant at +22 kB to +251 kB, so the
percentage penalty is worst where the unbounded archive is smallest. Only
*E. coli* (252.7 MiB) and *C. trachomatis* (205 MiB) stay under the memory cap,
and both still miss the size gate. Archive sizes and round trips are stable
across sessions, so this is a systemic result, not host noise.

| Run | HGS archive | Bounded archive | Δ bytes | Growth | Bounded peak RSS | Wall ratio | Verdict |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| broad-campy-natural-1 | 1,535,296 B | 1,785,856 B | +250,560 B | +16.32% | 272 MiB | 0.13x | fail |
| broad-chlamydia-natural-1 | 333,507 B | 355,433 B | +21,926 B | +6.57% | 205 MiB | 1.81x | fail |
| broad-mtb-natural-1 | 1,139,535 B | 1,364,788 B | +225,253 B | +19.77% | 273 MiB | 0.12x | fail |
| broad-pneumo-natural-1 | 1,754,114 B | 1,905,683 B | +151,569 B | +8.64% | 271 MiB | 0.13x | fail |
| broad-pseudo-natural-1 | 3,083,542 B | — | — | — | 208 MiB | — | fail (nonzero) |
| broad-staph-natural-1 | 1,807,399 B | 1,990,707 B | +183,308 B | +10.14% | 271 MiB | 0.13x | fail |

The *P. aeruginosa* entry is an `OutOfMemoryError` in the factorizer, not a
size or RSS overrun: the profile cannot complete on that input at all.

The *E. coli* diversity panel isolates the same limit. At a saturated 9 MiB
anchor the seed index is 67 MB (7.1 bytes per base: `SEED_STEP=4`, both
strands, 70% load), so the 160 MB heap is within a few MB of full once the
bank saturates. Low and medium both exhaust the heap in `factorize` — after
6.8 s at 203.7 MiB and 18.7 s at 194.7 MiB — while high completes in 49.7 s
at 272.2 MiB, byte-exact but over the 256 MiB cap. The low stratum dies early,
not on a tail block, so block size does not decide the outcome; the margin
between the live anchor index and the heap does, and GC timing selects the
block that fails. The 2.08 GiB scale probe passes because 55 copies of one
sequence never saturate the bank. The bound follows pangenome novelty, not
input bytes; the mechanism and the `--anchor-mib 4` control are in
[`SCALING.md`](SCALING.md).

### Within-species diversity gradient

<!-- DIVERSITY-TABLES:START -->
| Dataset | Stratum | Mean d | ANI | HGS | JARVIS3 | xz -9 | HGS vs JARVIS3 | HGS vs xz | Exact |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| campy | low | 0.01076 | 0.9892 | 606,441 | 735,705 | 668,728 | +17.57% | +9.31% | yes |
| campy | medium | 0.01162 | 0.9884 | 786,859 | 979,551 | 1,164,136 | +19.67% | +32.41% | yes |
| campy | high | 0.01901 | 0.9810 | 1,136,495 | 1,354,165 | 1,546,744 | +16.07% | +26.52% | yes |
| ecoli | low | 0.01188 | 0.9881 | 2,513,762 | 2,703,678 | 5,071,132 | +7.02% | +50.43% | yes |
| ecoli | medium | 0.02533 | 0.9747 | 3,381,136 | 3,375,860 | 5,484,492 | -0.16% | +38.35% | yes |
| ecoli | high | 0.02595 | 0.9740 | 4,477,584 | 4,577,435 | 6,484,272 | +2.18% | +30.95% | yes |
| staph | low | 0.00481 | 0.9952 | 891,829 | 1,003,593 | 1,704,184 | +11.14% | +47.67% | yes |
| staph | medium | 0.01140 | 0.9886 | 1,145,736 | 1,206,866 | 1,957,676 | +5.07% | +41.47% | yes |
| staph | high | 0.01519 | 0.9848 | 1,744,870 | 1,766,912 | 2,641,104 | +1.25% | +33.93% | yes |

### Trend

Advantage is `(baseline - HGS) / baseline`; positive means HyperGraSS is
smaller. Spearman is over three strata only and is indicative, not a
significance test.

| Dataset | Low | Medium | High | Spearman(d, advantage; n=3) |
| --- | ---: | ---: | ---: | ---: |
| campy | +17.57% | +19.67% | +16.07% | -0.50 |
| ecoli | +7.02% | -0.16% | +2.18% | -0.50 |
| staph | +11.14% | +5.07% | +1.25% | -1.00 |

Three patterns, only one of them clean:

- *S. aureus* loses advantage monotonically as diversity rises: +11.14%,
  +5.07%, +1.25% (Spearman -1.00).
- *E. coli* crosses zero at the medium stratum (-0.16%), the first
  HyperGraSS loss on any run. Its medium and high strata are nearly the same
  diversity (d = 0.02533 vs 0.02595), so this panel is effectively binary and
  the medium/high ordering is not a reliable trend.
- *C. jejuni* is non-monotonic: +17.57%, +19.67%, +16.07% (Spearman -0.50).
  Its low and medium strata are nearly the same diversity (d = 0.01076 vs
  0.01162), so with three points the non-monotonic middle is weak evidence.

The common factor is not Mash distance by itself. HyperGraSS wins when the
fixed 9 MiB anchor covers a large fraction of the input; the 5 Mb *E. coli*
genomes and their higher within-species structure leave less relative
coverage than the 2.8 Mb *S. aureus* genomes, which is why the crossover
appears there first. That is a hypothesis consistent with these three
species, not a measured law: this panel does not separate anchor coverage
from genome size or from the number of variable sites. An anchor-size sweep
on one species would test it directly.
<!-- DIVERSITY-TABLES:END -->

The full profile detail, including wall/CPU times and peak RSS, is in
`benchmarks/broad/summary.md`.
<!-- BROAD-TABLES:END -->

## Limits

- One host, one run per dataset. The host's wall and CPU times vary by more
  than 2× across sessions; archive sizes and round trips are the stable
  comparisons.
- The 72-genome pools are random samples of the NCBI RefSeq complete-genome
  catalogue returned by the query at fetch time. NCBI relevance order and the
  catalogue can change; the manifests lock the exact accessions and hashes,
  not a permanent snapshot.
- The natural collections target a fixed total size (about 120 MB), not a fixed
  strain count. Small genomes therefore contribute more records than large
  ones (for example, roughly 28 *M. tuberculosis* chromosomes versus 57
  *S. pneumoniae* or all 66 *C. trachomatis* records). Bits per base is the
  primary size metric; the diversity strata fix the strain count at up to 24
  to control this confound for the three gradient species.
- The Mash distance is an alignment-free approximation. It is adequate for
  ordering within-species diversity but is not a substitute for a pairwise ANI
  pipeline with alignment or an established ANI tool.
- Clustering on a distance matrix is sensitive to the pool composition. The
  strata are defined by the measured distances and the reported selection
  rule, not by a claim that they are universal diversity classes.
- The three-level gradients are uneven. *S. aureus* separates cleanly
  (d = 0.0048/0.0114/0.0152), but *C. jejuni* low and medium are close
  (0.0108/0.0116) and *E. coli* medium and high are close (0.0253/0.0260).
  A three-point Spearman coefficient is descriptive, not a significance test.
- Strata hold up to 24 strains, not a fixed number of bases. Total stratum
  size therefore follows genome size (about 40 MB for *C. jejuni*, 67 MB for
  *S. aureus*, and 120 MB for *E. coli*). Both panels use the same 64 MiB
  profile, but absolute archive sizes still follow input size and are not
  compared across species.
- *C. trachomatis* has only 66 eligible records in the pool, so its natural
  collection is about 69 MB rather than 120 MB. All other natural collections
  are near 120 MB.
- GeCo3 level 6 and GraSS are not run on every collection; the matrix records
  exactly where they are measured.

## Reproduce

```sh
python -B benchmarks/broad/fetch.py --name staph --organism "Staphylococcus aureus"
python -B benchmarks/broad/strata.py --dataset staph
python -B benchmarks/broad/plan.py
python -B benchmarks/broad/run.py --dry-run
python -B benchmarks/broad/run.py --panel cross-taxon
python -B benchmarks/broad/run.py --panel diversity
python -B benchmarks/broad/summarize.py
```

Pool FASTAs, subset FASTAs, and archives are gitignored. The manifests,
`diversity.json`, `matrix.json`, `report.json`, `verification.json`, and the
summary are the published evidence.
