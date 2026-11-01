"""Entropy-floor analysis for the 100:1 archive claim on the 24-genome pangenome.

The 100:1 target means an archive of <= 1.2 MB for the 120,108,223-base flat
A/C/G/T input. This script computes information-theoretic lower bounds from the
measured corpus statistics and from the best rates this project has measured,
then checks whether the target is reachable under any plausible design.

Two independent designs are bounded:

  A. Store the pangenome's distinct content (the 9.43 Mbp anchor bank) and encode
     each genome as edits.
  B. Store one genome as a reference and encode the other 23 as edits.

Both are lower bounds in the sense that every term uses the *best measured* rate
for that kind of data (QCRM 1.806 bits/base on a single E. coli genome), and the
variant-position term uses the independence model, which ignores any help from
shared phylogeny. Real compressors do worse on both counts (JARVIS3 l7 lands at
3,896,361 B against a computed floor of ~3.9 MB).

Usage:
    python -B benchmarks/entropy_floor.py [--json OUT]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

# Measured on benchmarks/specialized/runs/specialized-scale-1/input.seq.
BASES = 120_108_223
GENOMES = 24
GENOME_BASES = BASES / GENOMES

# Anchor bank saturates at 9,434,480 bases (anchor-mib >= 9 does not grow it).
ANCHOR_BASES = 9_434_480
# QCRM single-genome rate: 1,144,863 B for 5,071,579 bases.
QCRM_BITS_PER_BASE = 1_144_863 * 8 / 5_071_579

# Post-consensus SNP instances measured by STATS on the 3-block archive
# (359,880 + 470,075 + 747,274). Pre-consensus was 2,007,000.
SNPS = 359_880 + 470_075 + 747_274
SNPS_PRE_CONSENSUS = 2_007_000
# Literal (unmatched) bases that the archive must spell out.
LITERAL_BASES = 61_385 + 123_245 + 219_891

# Allele cost: a variant is one of the three other bases. A hard lower bound is
# 1 bit; uniform substitutions cost log2(3).
ALLELE_FLOOR_BITS = 1.0
ALLELE_UNIFORM_BITS = math.log2(3)

TARGET_BYTES = 1_200_000
BEST_SPECIALIZED_BYTES = 3_896_361   # JARVIS3 -l 7, identical input
CURRENT_BYTES = 4_153_891            # HyperGraSS, consensus anchor + range aux


def log2_choose(n: int, k: int) -> float:
    """log2 of the binomial coefficient, computed through lgamma."""
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2)


def bits_to_bytes(bits: float) -> float:
    return bits / 8.0


def analyse() -> dict:
    per_genome_snps = SNPS / GENOMES
    position_bits = GENOMES * log2_choose(round(GENOME_BASES), round(per_genome_snps))

    anchor_bits = ANCHOR_BASES * QCRM_BITS_PER_BASE
    literal_bits = LITERAL_BASES * QCRM_BITS_PER_BASE
    allele_low_bits = SNPS * ALLELE_FLOOR_BITS
    allele_uniform_bits = SNPS * ALLELE_UNIFORM_BITS

    case_a = bits_to_bytes(anchor_bits + position_bits + allele_low_bits + literal_bits)
    case_a_uniform = bits_to_bytes(anchor_bits + position_bits + allele_uniform_bits + literal_bits)

    # Case B omits the anchor: one genome at the QCRM rate, the other 23 as edits
    # against it. The position term is the same per-genome combinatorial cost.
    one_genome_bits = GENOME_BASES * QCRM_BITS_PER_BASE
    case_b = bits_to_bytes(
        one_genome_bits + (GENOMES - 1) * log2_choose(round(GENOME_BASES), round(per_genome_snps))
        + allele_low_bits + literal_bits
    )

    # Sensitivity: what would it take to fit 1.2 MB? Solve for the anchor rate
    # after granting a zero-cost variation layer, then repeat with the position
    # floor alone.
    anchor_rate_to_fit = TARGET_BYTES * 8 / ANCHOR_BASES
    position_only = bits_to_bytes(position_bits)

    return {
        "schema_version": 1,
        "kind": "entropy-floor-analysis",
        "input": {
            "bases": BASES,
            "genomes": GENOMES,
            "genome_bases_mean": GENOME_BASES,
            "flat_input_sha256": None,
        },
        "measured": {
            "qcrms_bits_per_base_single_genome": round(QCRM_BITS_PER_BASE, 6),
            "anchor_bases": ANCHOR_BASES,
            "snps_post_consensus": SNPS,
            "snps_pre_consensus": SNPS_PRE_CONSENSUS,
            "literal_bases": LITERAL_BASES,
            "hypergrass_archive_bytes": CURRENT_BYTES,
            "jarvis3_l7_archive_bytes": BEST_SPECIALIZED_BYTES,
            "target_100x_bytes": TARGET_BYTES,
        },
        "floor": {
            "anchor_bytes": round(bits_to_bytes(anchor_bits)),
            "snp_position_bytes": round(position_only),
            "snp_allele_bytes_lower_bound": round(bits_to_bytes(allele_low_bits)),
            "snp_allele_bytes_uniform": round(bits_to_bytes(allele_uniform_bits)),
            "literal_bytes": round(bits_to_bytes(literal_bits)),
            "case_a_total_bytes": round(case_a),
            "case_a_total_bytes_uniform_alleles": round(case_a_uniform),
            "case_b_total_bytes": round(case_b),
        },
        "sensitivity": {
            "anchor_rate_to_fit_target_bits_per_base": round(anchor_rate_to_fit, 4),
            "position_floor_exceeds_target": position_only > TARGET_BYTES,
            "position_over_target_ratio": round(position_only / TARGET_BYTES, 3),
            "case_a_over_target_ratio": round(case_a / TARGET_BYTES, 2),
            "case_b_over_target_ratio": round(case_b / TARGET_BYTES, 2),
        },
        "verdict": (
            "infeasible" if case_a > TARGET_BYTES and case_b > TARGET_BYTES else "plausible"
        ),
        "assumptions": [
            "QCRM's measured 1.806 bits/base on a single 5.07 Mbp E. coli genome is the "
            "best DNA rate measured in this project and is used as the floor for distinct "
            "sequence; no known method compresses E. coli DNA below ~1.7 bits/base.",
            "SNP positions are bounded by the independence model, which gives shared "
            "variants no credit; phylogenetic sharing can only lower this term.",
            "The anchor bank (9,434,480 bases) is the distinct non-redundant content of the "
            "24 genomes; any faithful archive must convey it or convey a reference genome "
            "plus the same amount of variation.",
            "Literal bases are the unmatched (indel/accessory) bases measured in the current "
            "archive; they are in addition to the SNP instances.",
        ],
    }


def render(result: dict) -> str:
    f = result["floor"]
    m = result["measured"]
    lines = [
        "# Entropy floor for the 100:1 target",
        "",
        f"Input: {m['bases'] if 'bases' in m else result['input']['bases']:,} bases, "
        f"{result['input']['genomes']} genomes.",
        "",
        "| Term | Bytes | Basis |",
        "|---|---:|---|",
        f"| Anchor content ({m['anchor_bases']:,} bases) | {f['anchor_bytes']:,} | "
        f"{result['measured']['qcrms_bits_per_base_single_genome']:.3f} bits/base |",
        f"| SNP positions ({m['snps_post_consensus']:,} instances) | {f['snp_position_bytes']:,} | "
        "independence model |",
        f"| SNP alleles (>=1 bit each) | {f['snp_allele_bytes_lower_bound']:,} | lower bound |",
        f"| Literal bases ({m['literal_bases']:,}) | {f['literal_bytes']:,} | measured DNA rate |",
        f"| **Case A total** | **{f['case_a_total_bytes']:,}** | anchor + edits |",
        f"| Case B total (one genome + 23 edited) | {f['case_b_total_bytes']:,} | reference + edits |",
        "",
        f"Target: {m['target_100x_bytes']:,} B. Current: {m['hypergrass_archive_bytes']:,} B. "
        f"JARVIS3 l7: {m['jarvis3_l7_archive_bytes']:,} B.",
        "",
        "The SNP-position term alone is "
        f"{result['sensitivity']['position_over_target_ratio']}x the entire target budget.",
        "",
        "## Verdict",
        "",
        f"**{result['verdict']}**",
        "",
        *[f"- {item}" for item in result["assumptions"]],
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--markdown", type=Path, default=None)
    args = parser.parse_args(argv)

    result = analyse()
    payload = json.dumps(result, indent=2) + "\n"
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload, encoding="utf-8")
    text = render(result)
    if args.markdown is not None:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(text, encoding="utf-8")
    print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
