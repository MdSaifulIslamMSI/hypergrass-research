"""Build the broad-experiment run matrix from the selected datasets.

Usage:
    python -B benchmarks/broad/plan.py
    python -B benchmarks/broad/plan.py --out benchmarks/broad/matrix.json
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent.parent
DATA = ROOT / "data"

CROSS_PROFILES = [
    "hypergrass-pangenome",
    "hypergrass-bounded-256m-tuned",
    "jarvis3-l7",
    "xz-9",
    "geco3-l6",
]
DIVERSITY_PROFILES = [
    "hypergrass-pangenome",
    "hypergrass-bounded-256m-tuned",
    "jarvis3-l7",
    "xz-9",
]

# GeCo3 level 6 is slow; keep it on three taxonomically spread collections.
GECO3_DATASETS = {"mtb", "staph", "campy"}
# The published GraSS ancestor is measured on a low- and a high-diversity set.
GRASS_DATASETS = {"mtb", "pseudo"}

DATASETS = [
    {
        "name": "mtb",
        "organism": "Mycobacterium tuberculosis",
        "panels": ["cross-taxon"],
        "natural": True,
        "strata": False,
    },
    {
        "name": "staph",
        "organism": "Staphylococcus aureus",
        "panels": ["cross-taxon", "diversity"],
        "natural": True,
        "strata": True,
    },
    {
        "name": "pneumo",
        "organism": "Streptococcus pneumoniae",
        "panels": ["cross-taxon"],
        "natural": True,
        "strata": False,
    },
    {
        "name": "pseudo",
        "organism": "Pseudomonas aeruginosa",
        "panels": ["cross-taxon"],
        "natural": True,
        "strata": False,
    },
    {
        "name": "campy",
        "organism": "Campylobacter jejuni",
        "panels": ["cross-taxon", "diversity"],
        "natural": True,
        "strata": True,
    },
    {
        "name": "chlamydia",
        "organism": "Chlamydia trachomatis",
        "panels": ["cross-taxon"],
        "natural": True,
        "strata": False,
    },
    {
        "name": "ecoli",
        "organism": "Escherichia coli",
        "panels": ["diversity"],
        "natural": False,
        "strata": True,
    },
]


def diversity_entry(path: Path, stratum: str) -> dict:
    if not path.exists():
        raise SystemExit(f"missing diversity selection: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    entry = payload["strata"][stratum]
    return {
        "mean_pairwise_distance": entry["mean_pairwise_distance"],
        "median_pairwise_distance": entry["median_pairwise_distance"],
        "mean_mash_ani": entry["mean_mash_ani"],
        "record_count": entry["record_count"],
        "total_bases": entry["total_bases"],
        "clusters_represented": entry["clusters_represented"],
        "file": entry["file"],
    }


def natural_entry(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing diversity selection: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    entry = payload["natural"]
    return {
        "mean_pairwise_distance": entry["mean_pairwise_distance"],
        "median_pairwise_distance": entry["median_pairwise_distance"],
        "mean_mash_ani": entry["mean_mash_ani"],
        "record_count": entry["record_count"],
        "total_bases": entry["total_bases"],
        "file": entry["file"],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(ROOT / "matrix.json"))
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    runs = []
    for dataset in DATASETS:
        name = dataset["name"]
        selection_path = DATA / name / "diversity.json"
        if dataset["natural"]:
            entry = natural_entry(selection_path)
            profiles = [value for value in CROSS_PROFILES if value != "geco3-l6" or name in GECO3_DATASETS]
            if name in GRASS_DATASETS:
                profiles.append("grass")
            runs.append(
                {
                    "run_id": f"broad-{name}-natural-1",
                    "panel": "cross-taxon",
                    "dataset": name,
                    "organism": dataset["organism"],
                    "stratum": "natural",
                    "source": (DATA / name / "natural.fasta").relative_to(PROJECT).as_posix(),
                    "profiles": profiles,
                    "diversity": entry,
                }
            )
        if dataset["strata"]:
            for stratum in ("low", "medium", "high"):
                entry = diversity_entry(selection_path, stratum)
                runs.append(
                    {
                        "run_id": f"broad-{name}-{stratum}-1",
                        "panel": "diversity",
                        "dataset": name,
                        "organism": dataset["organism"],
                        "stratum": stratum,
                        "source": (DATA / name / "strata" / f"{stratum}.fasta").relative_to(PROJECT).as_posix(),
                        "profiles": DIVERSITY_PROFILES,
                        "diversity": entry,
                    }
                )

    for run in runs:
        source = PROJECT / run["source"]
        if not source.exists():
            raise SystemExit(f"missing source FASTA: {source}")
        digest = hashlib.sha256()
        with source.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != run["diversity"]["file"]["sha256"]:
            raise SystemExit(f"source SHA-256 does not match selection: {run['run_id']}")
        run["diversity"]["file"]["path"] = run["source"]

    matrix = {
        "schema_version": 1,
        "kind": "broad-pangenome-matrix",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "profiles": {
            "cross-taxon": CROSS_PROFILES,
            "diversity": DIVERSITY_PROFILES,
            "geco3-datasets": sorted(GECO3_DATASETS),
            "grass-datasets": sorted(GRASS_DATASETS),
        },
        "runs": runs,
    }
    out = Path(args.out)
    out.write_text(json.dumps(matrix, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}: {len(runs)} runs")
    for run in runs:
        print(
            f"  {run['run_id']:32s} {run['diversity']['record_count']:3d} records "
            f"{run['diversity']['total_bases']/1e6:7.1f} Mb "
            f"d={run['diversity']['mean_pairwise_distance']:.5f} "
            f"{','.join(run['profiles'])}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
