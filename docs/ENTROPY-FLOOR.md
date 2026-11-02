# Entropy floor for the 100:1 target

Input: 120,108,223 bases, 24 genomes. Target: 1,200,000 B (100:1).

## Measured decomposition of the current archive

The tuned archive is **3,810,409 B** (0.254 bits/base), 2.21% smaller than
JARVIS3 level 7 (3,896,361 B). The bounded-memory variant is 3,997,232 B, so
holding 256 MiB only widens the gap; the argument below uses the smaller
unbounded archive. It splits into two parts:

| Part | Bytes | Basis |
|---|---:|---|
| Anchor (9,434,480 bases) + container | 2,257,234 | 1.914 bits/base |
| Gap (variant distances) | 882,215 | LZ77 + range coder |
| Strand bits | 239,808 | packed |
| Diff (SNP alleles) | 160,236 | range coder |
| Geom (match geometry) | 143,771 | digit planes |
| Mask | 80,818 | range coder |
| Literals | 30,647 | QCRM residuals |
| Donor, counts, modes | 15,680 | varints |
| **Total** | **3,810,409** | |

Two measured facts bound the target:

1. **The non-anchor channels alone are 1,553,175 B**, which is 1.29x the entire
   1.2 MB budget. Reaching 100:1 would require coding all variant metadata for
   free *and* shrinking it by 23% on top.
2. **The anchor is 1.88x the target by itself.** Fitting the 9,434,480 distinct
   anchor bases in 1.2 MB requires 1.018 bits/base. The best measured DNA rate in
   this project is QCRM at 1.806 bits/base on a single *E. coli* genome, and
   JARVIS3 reaches 2.035 bits/base on this anchor. No measured DNA coder is
   within 1.7x of 1.018 bits/base.

## Earlier estimate, corrected

An earlier version of this note summed an anchor term (1.806 bits/base), an
independence-model SNP-position term, one bit per SNP allele, and a literal
term to 3,933,118 B and called it a floor. That was an *estimate*, not a bound,
and the finished archive beats it by 122,709 B: the independence model charged
1,514,899 B for SNP positions, while the measured gap + diff + geom + mask
channels total 1,267,040 B. It should not have been labelled a floor.

## Verdict

**infeasible**

- 100:1 needs 1.018 bits/base on the distinct anchor content alone, 1.77x better
  than the best DNA rate measured here and below every specialized tool measured
  on this input.
- Even with a free anchor, the measured variant metadata is already 29% over the
  whole budget.
- The gap from the current archive to 1.2 MB is 3.18x. Closing it would need a
  coding gain larger than the entire gain HyperGraSS achieved over xz (5,663,640
  -> 3,810,409 B, 1.49x), repeated more than four times.

These are measured bounds from four coders and a channel decomposition, not an
information-theoretic proof. A radically different representation (for example,
one that reconstructs genomes from a generative model) could in principle do
better, but no measured method is close.
