"""Build the natural collection and diversity strata for one broad dataset.

Usage:
    python -B benchmarks/broad/strata.py --dataset staph
    python -B benchmarks/broad/strata.py --dataset staph --strata low,medium,high
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import sys
from pathlib import Path

import numpy as np

import diversity


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent.parent
DATA = ROOT / "data"


def portable(entry: dict, project: Path) -> dict:
    path = Path(entry["path"])
    entry["path"] = path.resolve().relative_to(project.resolve()).as_posix()
    return entry


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--target-bases", type=int, default=120_000_000)
    parser.add_argument("--strata", default="low,medium,high")
    parser.add_argument("--stratum-size", type=int, default=24)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--refresh-sketches", action="store_true")
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    dataset = DATA / args.dataset
    manifest_path = dataset / "manifest.json"
    if not manifest_path.exists():
        raise SystemExit(f"missing pool manifest: {manifest_path}")
    manifest = diversity.load_manifest(manifest_path)
    pool = dataset / manifest["fasta"]["path"]
    if not pool.exists():
        raise SystemExit(f"missing pool FASTA: {pool}")
    if pool.stat().st_size != manifest["fasta"]["size_bytes"]:
        raise SystemExit("pool size does not match manifest")

    output = dataset / "diversity.json"
    if output.exists() and not args.force:
        print(f"exists: {output}; use --force to replace")
        return 0

    print(f"sketching {len(manifest['records'])} records from {pool.name} ...", flush=True)
    sketches_all = diversity.sketch_manifest(
        pool, manifest, cache_path=dataset / "sketches.npz",
        refresh=args.refresh_sketches,
    )
    print("computing pairwise distances ...", flush=True)
    distances_all = diversity.distance_matrix(sketches_all)
    keep, qc = diversity.qc_filter(sketches_all, distances_all)
    sketches = [sketches_all[index] for index in keep]
    records = [manifest["records"][index] for index in keep]
    distances = distances_all[np.ix_(keep, keep)]
    print(
        f"qc: kept {qc['kept']}/{len(sketches_all)} records "
        f"(dropped {qc['dropped']}, relaxed={qc['relaxed']})"
    )

    natural = diversity.natural_selection(distances, sketches, args.target_bases)
    natural_path = dataset / "natural.fasta"
    natural_file = diversity.write_subset(
        pool, records, natural["indices"], natural_path
    )
    natural["file"] = portable(natural_file, PROJECT)
    print(
        f"natural: {natural['record_count']} records, "
        f"{natural['total_bases']:,} bases, "
        f"mean distance {natural['mean_pairwise_distance']:.6f}, "
        f"{natural_path.name}"
    )

    wanted = [name.strip() for name in args.strata.split(",") if name.strip()]
    result: dict = {
        "schema_version": 1,
        "kind": "broad-diversity-selection",
        "dataset": args.dataset,
        "organism": manifest["organism"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "pool": {
            "manifest": manifest_path.relative_to(PROJECT).as_posix(),
            "fasta": pool.relative_to(PROJECT).as_posix(),
            "records": len(manifest["records"]),
            "records_used": int(len(keep)),
            "total_bases": manifest["total_bases"],
            "sha256": manifest["fasta"]["sha256"],
        },
        "qc": qc,
        "metric": {
            "k": diversity.K,
            "sampling": "all canonical k-mers",
            "sketch_size": diversity.SKETCH_SIZE,
            "hash_seed": diversity.HASH_SEED,
            "distance": "Mash distance from canonical k-mer Jaccard",
        },
        "natural": natural,
        "strata": {},
    }

    if wanted:
        selection = diversity.select_strata(
            distances, sketches, target_n=args.stratum_size
        )
        print(
            f"clusters: {selection['k']}; strategy: {selection['strategy']}; "
            f"monotonic: {selection['monotonic']}"
        )
        for name in wanted:
            if name not in selection["strata"]:
                raise SystemExit(f"unknown stratum: {name}")
            entry = selection["strata"][name]
            path = dataset / "strata" / f"{name}.fasta"
            entry["file"] = portable(
                diversity.write_subset(pool, records, entry["indices"], path),
                PROJECT,
            )
            print(
                f"  {name:6s}: {entry['record_count']} records, "
                f"{entry['total_bases']:,} bases, "
                f"mean distance {entry['mean_pairwise_distance']:.6f}, "
                f"ANI~{entry['mean_mash_ani']:.4f}"
            )
        result["clusters"] = selection["clusters"]
        result["strategy"] = selection["strategy"]
        result["stratum_size"] = selection["stratum_size"]
        result["monotonic"] = selection["monotonic"]
        result["span"] = selection["span"]
        result["strata"] = selection["strata"]

    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
