"""Measured feasibility analysis for the 100:1 archive claim.

The 100:1 target is <= 1,200,000 bytes for the 120,108,223-base flat A/C/G/T
input. This script records the measured channel decomposition of the finished
3,810,409-byte archive and checks whether the target is reachable given the DNA
rates this project has actually measured.

The argument is deliberately labelled *measured evidence*, not an
information-theoretic proof. Two independent measured components already exceed
the whole budget:

  * the anchor channel alone is 2,257,234 B, 1.88x the target, and would need
    1.018 bits/base versus the best measured rate of 1.806 bits/base;
  * the non-anchor channels alone are 1,553,175 B, 1.29x the target.

An earlier version summed an anchor rate, an independence-model SNP-position
term, one bit per allele, and a literal term to 3,933,118 B and called it a
floor. The finished archive is smaller than that number, so the sum is retained
below as an *estimate*, not a bound.

Usage:
    python -B benchmarks/entropy_floor.py [--json OUT] [--markdown OUT]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

BASES = 120_108_223
GENOMES = 24
GENOME_BASES = BASES / GENOMES
FLAT_INPUT_SHA256 = "51cb3ecb2e9f532833caa0b60345f37d0342474e9f4a239801486f9bc9219acb"

# Saturated anchor bank measured on this input (anchor-mib >= 9 does not grow).
ANCHOR_BASES = 9_434_480

# QCRM single-genome rate: 1,144,863 B for 5,071,579 bases.
QCRM_BITS_PER_BASE = 1_144_863 * 8 / 5_071_579
# JARVIS3 l7 anchor rate on this input: 2.035 bits/base (measured).
JARVIS3_ANCHOR_BITS_PER_BASE = 2.035

TARGET_BYTES = 1_200_000
JARVIS3_BYTES = 3_896_361
CURRENT_BYTES = 3_810_409

# Measured channel decomposition of the tuned archive. These are the
# ArchiveByteAccounting figures reported in docs/ENTROPY-FLOOR.md; they sum
# exactly to CURRENT_BYTES.
CHANNEL_BYTES = {
    "anchor_and_container": 2_257_234,
    "gap": 882_215,
    "strand": 239_808,
    "diff": 160_236,
    "geom": 143_771,
    "mask": 80_818,
    "literals": 30_647,
    "donor_counts_modes": 15_680,
}

# Earlier independence-model estimate. This is not a lower bound: the finished
# archive beats it, because the coder exploits cross-genome sharing that the
# independence model gives no credit.
SNPS_POST_CONSENSUS = 359_880 + 470_075 + 747_274
SNPS_PRE_CONSENSUS = 2_007_000
LITERAL_BASES = 61_385 + 123_245 + 219_891


def bits_to_bytes(bits: float) -> float:
    return bits / 8.0


def analyse() -> dict:
    non_anchor_bytes = sum(
        value for key, value in CHANNEL_BYTES.items() if key != "anchor_and_container"
    )
    anchor_bytes = CHANNEL_BYTES["anchor_and_container"]
    if sum(CHANNEL_BYTES.values()) != CURRENT_BYTES:
        raise ValueError("channel decomposition does not sum to the measured archive")

    anchor_bits_per_base = anchor_bytes * 8 / ANCHOR_BASES
    needed_anchor_bits_per_base = TARGET_BYTES * 8 / ANCHOR_BASES
    rate_gap = QCRM_BITS_PER_BASE / needed_anchor_bits_per_base

    # Retained only to show why the earlier estimate was withdrawn.
    position_bytes = 1_514_899
    independence_estimate = 3_933_118

    return {
        "schema_version": 2,
        "kind": "entropy-feasibility-analysis",
        "input": {
            "bases": BASES,
            "genomes": GENOMES,
            "genome_bases_mean": GENOME_BASES,
            "flat_input_sha256": FLAT_INPUT_SHA256,
            "target_100x_bytes": TARGET_BYTES,
        },
        "measured": {
            "hypergrass_archive_bytes": CURRENT_BYTES,
            "hypergrass_bits_per_base": round(CURRENT_BYTES * 8 / BASES, 6),
            "jarvis3_l7_archive_bytes": JARVIS3_BYTES,
            "anchor_bases": ANCHOR_BASES,
            "channel_bytes": CHANNEL_BYTES,
            "anchor_bytes": anchor_bytes,
            "non_anchor_bytes": non_anchor_bytes,
            "anchor_bits_per_base": round(anchor_bits_per_base, 6),
            "qcrms_bits_per_base_single_genome": round(QCRM_BITS_PER_BASE, 6),
            "jarvis3_anchor_bits_per_base": JARVIS3_ANCHOR_BITS_PER_BASE,
        },
        "sensitivity": {
            "anchor_rate_needed_for_target_bits_per_base": round(
                needed_anchor_bits_per_base, 4
            ),
            "best_measured_over_needed_anchor_rate": round(rate_gap, 3),
            "anchor_over_target_ratio": round(anchor_bytes / TARGET_BYTES, 3),
            "non_anchor_over_target_ratio": round(non_anchor_bytes / TARGET_BYTES, 3),
            "current_over_target_ratio": round(CURRENT_BYTES / TARGET_BYTES, 3),
        },
        "withdrawn_independence_estimate": {
            "snp_position_bytes": position_bytes,
            "case_a_total_bytes": independence_estimate,
            "note": (
                "Not a lower bound. The finished 3,810,409-byte archive is "
                "122,709 B smaller, because the independence model gives no "
                "credit for cross-genome sharing."
            ),
        },
        "verdict": "infeasible",
        "verdict_basis": (
            "The measured anchor channel alone is 2,257,234 B (1.88x the target) "
            "and would need 1.018 bits/base, 1.77x better than the best DNA rate "
            "measured here. Independently, the measured non-anchor channels alone "
            "are 1,553,175 B (1.29x the target). Both components exceed the whole "
            "budget. This is measured evidence, not an information-theoretic proof."
        ),
        "assumptions": [
            "The channel bytes are the measured ArchiveByteAccounting decomposition "
            "of the finished tuned archive and sum exactly to 3,810,409 B.",
            "The best measured DNA rate is QCRM at 1.806 bits/base on a single "
            "5.07 Mbp E. coli genome; JARVIS3 reaches 2.035 bits/base on this "
            "anchor, and the finished archive reaches 1.914 bits/base on it.",
            "The anchor bank (9,434,480 bases) is the distinct non-redundant "
            "content of the 24 genomes; any faithful archive must convey it or a "
            "reference plus the same amount of variation.",
            "A radically different generative representation could in principle "
            "beat these measured bounds, but no measured method is close.",
        ],
    }


def render(result: dict) -> str:
    m = result["measured"]
    s = result["sensitivity"]
    lines = [
        "# Measured feasibility analysis for the 100:1 target",
        "",
        f"Input: {result['input']['bases']:,} bases, {result['input']['genomes']} genomes, "
        f"SHA-256 `{result['input']['flat_input_sha256']}`.",
        f"Target: {result['input']['target_100x_bytes']:,} B. "
        f"Current: {m['hypergrass_archive_bytes']:,} B. "
        f"JARVIS3 l7: {m['jarvis3_l7_archive_bytes']:,} B.",
        "",
        "| Measured component | Bytes | Ratio to target |",
        "|---|---:|---:|",
        f"| Anchor ({m['anchor_bases']:,} bases, {m['anchor_bits_per_base']:.3f} bits/base) "
        f"| {m['anchor_bytes']:,} | {s['anchor_over_target_ratio']:.2f}x |",
        f"| Non-anchor channels | {m['non_anchor_bytes']:,} | "
        f"{s['non_anchor_over_target_ratio']:.2f}x |",
        "",
        f"Fitting the anchor in the whole budget needs "
        f"{s['anchor_rate_needed_for_target_bits_per_base']:.4f} bits/base, "
        f"which is {s['best_measured_over_needed_anchor_rate']:.2f}x below the best "
        "measured DNA rate.",
        "",
        "## Verdict",
        "",
        f"**{result['verdict']}** - {result['verdict_basis']}",
        "",
        "## Withdrawn estimate",
        "",
        result["withdrawn_independence_estimate"]["note"],
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
