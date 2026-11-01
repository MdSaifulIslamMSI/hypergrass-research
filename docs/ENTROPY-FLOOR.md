# Entropy floor for the 100:1 target

Input: 120,108,223 bases, 24 genomes.

| Term | Bytes | Basis |
|---|---:|---|
| Anchor content (9,434,480 bases) | 2,129,748 | 1.806 bits/base |
| SNP positions (1,577,229 instances) | 1,514,899 | independence model |
| SNP alleles (>=1 bit each) | 197,154 | lower bound |
| Literal bases (404,521) | 91,317 | measured DNA rate |
| **Case A total** | **3,933,118** | anchor + edits |
| Case B total (one genome + 23 edited) | 2,869,971 | reference + edits |

Target: 1,200,000 B. Current: 4,153,891 B. JARVIS3 l7: 3,896,361 B.

The SNP-position term alone is 1.262x the entire target budget.

## Verdict

**infeasible**

- QCRM's measured 1.806 bits/base on a single 5.07 Mbp E. coli genome is the best DNA rate measured in this project and is used as the floor for distinct sequence; no known method compresses E. coli DNA below ~1.7 bits/base.
- SNP positions are bounded by the independence model, which gives shared variants no credit; phylogenetic sharing can only lower this term.
- The anchor bank (9,434,480 bases) is the distinct non-redundant content of the 24 genomes; any faithful archive must convey it or convey a reference genome plus the same amount of variation.
- Literal bases are the unmatched (indel/accessory) bases measured in the current archive; they are in addition to the SNP instances.
