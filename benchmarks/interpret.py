"""Compute the claims that the benchmark documents are allowed to make.

Prints a small set of counts and medians from the audited JSON. This is a
decision aid, not a source of prose: every number still comes from the report.
"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

HGS = ("hypergrass", "hypergrass-optimize", "hypergrass-lean")
FROZEN = ("gzip-9", "xz-9", "jarvis3-l7", "jarvis3-l26", "geco3-l9")
SUPPLEMENTARY = ("geco3-l6", "geco3-l3")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--scale")
    args = parser.parse_args(argv)
    report = json.loads(Path(args.corpus).read_text(encoding="utf-8"))
    summary = report["summary"]
    datasets = []
    for row in summary:
        if row["dataset"] not in datasets:
            datasets.append(row["dataset"])

    def get(dataset, profile):
        for row in summary:
            if row["dataset"] == dataset and row["profile"] == profile:
                return row
        return None

    def best(dataset, exclude_l9=False):
        candidates = []
        for row in summary:
            if row["dataset"] != dataset:
                continue
            if row["profile"] in HGS:
                continue
            if exclude_l9 and row["profile"] == "geco3-l9":
                continue
            if row["status"] == "success" and row.get("artifact_bytes"):
                candidates.append(row)
        return min(candidates, key=lambda r: r["artifact_bytes"]) if candidates else None

    wins_all = []
    wins_feasible = []
    losses = []
    deltas = []
    lean_overhead = []
    for dataset in datasets:
        hgs = get(dataset, "hypergrass")
        if not hgs or hgs["status"] != "success" or not hgs.get("artifact_bytes"):
            continue
        all_best = best(dataset)
        feasible_best = best(dataset, exclude_l9=True)
        if all_best:
            if hgs["artifact_bytes"] < all_best["artifact_bytes"]:
                wins_all.append((dataset, all_best))
            else:
                losses.append((dataset, all_best, (hgs["artifact_bytes"] - all_best["artifact_bytes"]) / all_best["artifact_bytes"] * 100))
        if feasible_best:
            delta = (feasible_best["artifact_bytes"] - hgs["artifact_bytes"]) / feasible_best["artifact_bytes"] * 100
            deltas.append(delta)
            if delta > 0:
                wins_feasible.append((dataset, feasible_best, delta))
        lean = get(dataset, "hypergrass-lean")
        if lean and lean["status"] == "success" and lean.get("artifact_bytes"):
            lean_overhead.append((lean["artifact_bytes"] - hgs["artifact_bytes"]) / hgs["artifact_bytes"] * 100)

    print(f"datasets with a successful default archive: {len(deltas)}")
    print(f"default smaller than best of all successful profiles: {len(wins_all)}/{len(deltas)}")
    print(f"default smaller than best feasible (excluding geco3-l9): {len(wins_feasible)}/{len(deltas)}")
    if deltas:
        print(f"median delta vs best feasible (+ = HGS smaller): {statistics.median(deltas):+.2f}%")
        print(f"range: {min(deltas):+.2f}% .. {max(deltas):+.2f}%")
    if lean_overhead:
        print(f"lean vs default archive overhead: median {statistics.median(lean_overhead):+.2f}%, max {max(lean_overhead):+.2f}%")
    print("losses:")
    for dataset, row, delta in losses:
        print(f"  {dataset}: {row['profile']} {row['artifact_bytes']:,} B; HGS larger by {delta:.2f}%")
    print("wins vs all successful profiles:")
    for dataset, row in wins_all:
        print(f"  {dataset}: best baseline {row['profile']} {row['artifact_bytes']:,} B")
    print("geco3-l9 cells:")
    for row in summary:
        if row["profile"] == "geco3-l9":
            print(f"  {row['dataset']}: {row['status']}, {row.get('artifact_bytes')}")

    if args.scale:
        scale = json.loads(Path(args.scale).read_text(encoding="utf-8"))
        print("scale test:")
        for entry in scale["results"]:
            print(
                f"  {entry['profile']}: status={entry['status']}, "
                f"bytes={entry.get('artifact_bytes')}, ratio={entry.get('ratio')}, "
                f"rss={entry.get('peak_rss_bytes')}, wall={entry.get('wall_seconds')}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
