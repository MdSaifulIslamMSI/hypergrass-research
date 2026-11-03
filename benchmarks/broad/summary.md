# Broad pangenome experiment

Diversity is the mean pairwise Mash distance from canonical 21-mer Jaccard; ANI is the Mash approximation. Lower archive and higher diversity are independent axes.

## Cross-taxon natural collections

Positive `HGS vs best` means the HyperGraSS archive is smaller than the smallest baseline archive on that row.

| Dataset | Organism | Records | Mb | Mash d | ANI | HGS | JARVIS3 | xz -9 | GeCo3 l6 | Best baseline | HGS vs best | Exact |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| campy | Campylobacter jejuni | 72 | 120 | 0.01897 | 0.9810 | 1,535,296 | 2,017,046 | 2,019,208 | 7,673,951 | jarvis3-l7 | +23.88% | yes |
| chlamydia | Chlamydia trachomatis | 66 | 69 | 0.00450 | 0.9955 | 333,507 | 472,224 | 634,356 | — | jarvis3-l7 | +29.38% | yes |
| mtb | Mycobacterium tuberculosis | 28 | 124 | 0.00083 | 0.9992 | 1,139,535 | 1,418,155 | 1,433,772 | 6,772,719 | jarvis3-l7 | +19.65% | yes |
| pneumo | Streptococcus pneumoniae | 57 | 120 | 0.01223 | 0.9878 | 1,754,114 | 1,899,711 | 2,328,208 | — | jarvis3-l7 | +7.66% | yes |
| pseudo | Pseudomonas aeruginosa | 19 | 126 | 0.01133 | 0.9887 | 3,083,542 | 3,350,844 | 6,005,620 | — | jarvis3-l7 | +7.98% | yes |
| staph | Staphylococcus aureus | 43 | 121 | 0.01246 | 0.9875 | 1,807,399 | 1,895,866 | 2,741,684 | 6,427,702 | jarvis3-l7 | +4.67% | yes |

### Profile detail

| Dataset | Profile | Archive | bits/base | Ratio | Wall | CPU | Peak RSS | Exact | Δ vs JARVIS3 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: |
| campy | hypergrass-pangenome | 1,535,296 B | 0.102 | 78.19x | 212.9 s | 213.5 s | 897 MiB | yes | -23.88% |
| campy | hypergrass-bounded-256m-tuned | 1,785,856 B | 0.119 | 67.22x | 26.9 s | 26.5 s | 272 MiB | yes | -11.46% |
| campy | jarvis3-l7 | 2,017,046 B | 0.134 | 59.52x | 65.9 s | 12.5 s | 284 MiB | yes | +0.00% |
| campy | geco3-l6 | 7,673,951 B | 0.511 | 15.64x | 168.5 s | 159.0 s | 506 MiB | yes | +280.45% |
| campy | xz-9 | 2,019,208 B | 0.135 | 59.45x | 64.5 s | 64.8 s | 726 MiB | yes | +0.11% |
| chlamydia | hypergrass-pangenome | 333,507 B | 0.039 | 206.44x | 6.0 s | 6.8 s | 371 MiB | yes | -29.38% |
| chlamydia | hypergrass-bounded-256m-tuned | 355,433 B | 0.041 | 193.71x | 10.9 s | 10.0 s | 183 MiB | yes | -24.73% |
| chlamydia | jarvis3-l7 | 472,224 B | 0.055 | 145.80x | 35.7 s | 7.3 s | 274 MiB | yes | +0.00% |
| chlamydia | xz-9 | 634,356 B | 0.074 | 108.54x | 24.5 s | 24.7 s | 645 MiB | yes | +34.33% |
| mtb | hypergrass-pangenome | 1,139,535 B | 0.074 | 108.39x | 283.8 s | 285.2 s | 898 MiB | yes | -19.65% |
| mtb | hypergrass-bounded-256m-tuned | 1,364,788 B | 0.088 | 90.50x | 34.1 s | 35.0 s | 273 MiB | yes | -3.76% |
| mtb | jarvis3-l7 | 1,418,155 B | 0.092 | 87.09x | 130.4 s | 54.8 s | 285 MiB | yes | +0.00% |
| mtb | geco3-l6 | 6,772,719 B | 0.439 | 18.24x | 361.0 s | 343.9 s | 509 MiB | yes | +377.57% |
| mtb | xz-9 | 1,433,772 B | 0.093 | 86.14x | 131.5 s | 133.0 s | 729 MiB | yes | +1.10% |
| mtb | grass | 13,327,718 B | 0.863 | 9.27x | 109.1 s | 33.1 s | 1252 MiB | yes | +839.79% |
| pneumo | hypergrass-pangenome | 1,754,114 B | 0.117 | 68.47x | 714.8 s | 718.8 s | 923 MiB | yes | -7.66% |
| pneumo | hypergrass-bounded-256m-tuned | 1,905,683 B | 0.127 | 63.02x | 91.1 s | 93.6 s | 271 MiB | yes | +0.31% |
| pneumo | jarvis3-l7 | 1,899,711 B | 0.127 | 63.22x | 74.4 s | 19.8 s | 286 MiB | yes | +0.00% |
| pneumo | xz-9 | 2,328,208 B | 0.155 | 51.58x | 73.0 s | 73.2 s | 726 MiB | yes | +22.56% |
| pseudo | hypergrass-pangenome | 3,083,542 B | 0.196 | 40.92x | 666.5 s | 671.9 s | 1027 MiB | yes | -7.98% |
| pseudo | jarvis3-l7 | 3,350,844 B | 0.212 | 37.66x | 114.8 s | 46.7 s | 287 MiB | yes | +0.00% |
| pseudo | xz-9 | 6,005,620 B | 0.381 | 21.01x | 386.8 s | 386.8 s | 736 MiB | yes | +79.23% |
| pseudo | grass | 22,616,774 B | 1.434 | 5.58x | 115.7 s | 34.1 s | 1260 MiB | yes | +574.96% |
| staph | hypergrass-pangenome | 1,807,399 B | 0.120 | 66.92x | 284.4 s | 287.9 s | 891 MiB | yes | -4.67% |
| staph | hypergrass-bounded-256m-tuned | 1,990,707 B | 0.132 | 60.76x | 37.0 s | 38.1 s | 271 MiB | yes | +5.00% |
| staph | jarvis3-l7 | 1,895,866 B | 0.125 | 63.80x | 136.3 s | 60.4 s | 286 MiB | yes | +0.00% |
| staph | geco3-l6 | 6,427,702 B | 0.425 | 18.82x | 337.4 s | 325.7 s | 512 MiB | yes | +239.04% |
| staph | xz-9 | 2,741,684 B | 0.181 | 44.11x | 134.7 s | 135.4 s | 728 MiB | yes | +44.61% |

### Bounded-memory gate (≤256 MiB, ≤2% growth, ≤2× wall)

Reference is the unbounded `hypergrass-pangenome` archive for the same input.

| Run | HGS archive | Bounded archive | Δ bytes | Growth | Bounded peak RSS | Wall ratio | Verdict |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| broad-campy-high-1 | 1,136,495 B | 1,125,629 B | -10,866 B | -0.96% | 226 MiB | 1.09x | pass |
| broad-campy-low-1 | 606,441 B | 630,254 B | +23,813 B | +3.93% | 215 MiB | 1.25x | fail |
| broad-campy-medium-1 | 786,859 B | 799,633 B | +12,774 B | +1.62% | 222 MiB | 1.29x | pass |
| broad-campy-natural-1 | 1,535,296 B | 1,785,856 B | +250,560 B | +16.32% | 272 MiB | 0.13x | fail |
| broad-chlamydia-natural-1 | 333,507 B | 355,433 B | +21,926 B | +6.57% | 205 MiB | 1.81x | fail |
| broad-ecoli-high-1 | 4,477,584 B | 4,636,907 B | +159,323 B | +3.56% | 272 MiB | 0.18x | fail |
| broad-ecoli-low-1 | 2,513,762 B | — | — | — | 204 MiB | — | fail (nonzero) |
| broad-ecoli-medium-1 | 3,381,136 B | — | — | — | 195 MiB | — | fail (nonzero) |
| broad-mtb-natural-1 | 1,139,535 B | 1,364,788 B | +225,253 B | +19.77% | 273 MiB | 0.12x | fail |
| broad-pneumo-natural-1 | 1,754,114 B | 1,905,683 B | +151,569 B | +8.64% | 271 MiB | 0.13x | fail |
| broad-pseudo-natural-1 | 3,083,542 B | — | — | — | 208 MiB | — | fail (nonzero) |
| broad-staph-high-1 | 1,744,870 B | 1,820,635 B | +75,765 B | +4.34% | 272 MiB | 1.07x | fail |
| broad-staph-low-1 | 891,829 B | 955,373 B | +63,544 B | +7.13% | 224 MiB | 1.39x | fail |
| broad-staph-medium-1 | 1,145,736 B | 1,155,136 B | +9,400 B | +0.82% | 238 MiB | 1.24x | pass |
| broad-staph-natural-1 | 1,807,399 B | 1,990,707 B | +183,308 B | +10.14% | 271 MiB | 0.13x | fail |

## Diversity gradient

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

Advantage is `(baseline - HGS) / baseline`; positive means HyperGraSS is smaller. Spearman is over three strata only and is indicative, not a significance test.

| Dataset | Low | Medium | High | Spearman(d, advantage; n=3) |
| --- | ---: | ---: | ---: | ---: |
| campy | +17.57% | +19.67% | +16.07% | -0.50 |
| ecoli | +7.02% | -0.16% | +2.18% | -0.50 |
| staph | +11.14% | +5.07% | +1.25% | -1.00 |

## Coverage

- Runs in matrix: 15
- Runs with a report: 15
- Failed profile/run entries: 3

```
broad-pseudo-natural-1/hypergrass-bounded-256m-tuned: status=nonzero exact=False verification=None error=None
broad-ecoli-low-1/hypergrass-bounded-256m-tuned: status=nonzero exact=False verification=None error=None
broad-ecoli-medium-1/hypergrass-bounded-256m-tuned: status=nonzero exact=False verification=None error=None
```
