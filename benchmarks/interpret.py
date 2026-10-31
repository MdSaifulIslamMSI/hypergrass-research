"""Compute the claims that the benchmark documents are allowed to make.

Prints a small set of counts and medians from the audited JSON. This is a
decision aid, not a source of prose: every number still comes from the report.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summarize import corpus_coverage  # noqa: E402

HGS = ("hypergrass", "hypergrass-optimize", "hypergrass-lean")
HIGH_EFFORT = ("geco3-l9", "jarvis3-l26")


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

    present, planned, complete, total = corpus_coverage(report)
    print(
        f"corpus coverage: {present}/{planned} results, {complete}/{total} datasets "
        f"complete, report status {report.get('status')!r}"
    )

    def get(dataset, profile):
        for row in summary:
            if row["dataset"] == dataset and row["profile"] == profile:
                return row
        return None

    def best(dataset, exclude_high_effort=False):
        candidates = []
        for row in summary:
            if row["dataset"] != dataset:
                continue
            if row["profile"] in HGS:
                continue
            if exclude_high_effort and row["profile"] in HIGH_EFFORT:
                continue
            if row["status"] == "success" and row.get("artifact_bytes"):
                candidates.append(row)
        return min(candidates, key=lambda r: r["artifact_bytes"]) if candidates else None

    wins = []
    losses = []
    deltas = []
    deltas_feasible = []
    lean_overhead = []
    for dataset in datasets:
        hgs = get(dataset, "hypergrass")
        if not hgs or hgs["status"] != "success" or not hgs.get("artifact_bytes"):
            continue
        best_all = best(dataset)
        best_feasible = best(dataset, exclude_high_effort=True)
        if best_all:
            delta = (best_all["artifact_bytes"] - hgs["artifact_bytes"]) / best_all["artifact_bytes"] * 100
            deltas.append(delta)
            if delta > 0:
                wins.append((dataset, best_all, delta))
            else:
                losses.append((dataset, best_all, delta))
        if best_feasible:
            deltas_feasible.append(
                (best_feasible["artifact_bytes"] - hgs["artifact_bytes"])
                / best_feasible["artifact_bytes"] * 100
            )
        lean = get(dataset, "hypergrass-lean")
        if lean and lean["status"] == "success" and lean.get("artifact_bytes"):
            lean_overhead.append((lean["artifact_bytes"] - hgs["artifact_bytes"]) / hgs["artifact_bytes"] * 100)

    print(f"datasets with a successful default archive: {len(deltas)}")
    print(f"default smaller than best of all successful profiles: {len(wins)}/{len(deltas)}")
    if deltas:
        print(f"median delta vs best successful baseline (+ = HGS smaller): {statistics.median(deltas):+.2f}%")
        print(f"range: {min(deltas):+.2f}% .. {max(deltas):+.2f}%")
    if deltas_feasible:
        print(
            "median delta vs best baseline excluding high-effort profiles: "
            f"{statistics.median(deltas_feasible):+.2f}%"
        )
    if lean_overhead:
        print(
            f"lean vs default archive overhead: median {statistics.median(lean_overhead):+.2f}%, "
            f"max {max(lean_overhead):+.2f}%"
        )
    print("head-to-head:")
    for dataset, row, delta in wins:
        print(f"  {dataset}: HGS smaller than {row['profile']} by {delta:.2f}%")
    for dataset, row, delta in losses:
        print(f"  {dataset}: HGS larger than {row['profile']} by {-delta:.2f}%")
    print("high-effort cells:")
    for row in summary:
        if row["profile"] in HIGH_EFFORT:
            print(f"  {row['dataset']} / {row['profile']}: {row['status']}")

    if args.scale:
        scale = json.loads(Path(args.scale).read_text(encoding="utf-8"))
        print("scale test:")
        for entry in scale["results"]:
            if entry["status"] != "success":
                print(f"  {entry['profile']}: {entry['status']} ({entry.get('error')})")
                continue
            queries = entry.get("queries", [])
            exact = sum(1 for query in queries if query["exact"])
            print(
                f"  {entry['profile']}: {entry['artifact_bytes']:,} B, "
                f"{entry['ratio']:.2f}x, {entry['bits_per_byte']:.3f} bits/byte, "
                f"{entry['wall_seconds']:.1f}s wall, "
                f"{entry['peak_rss_bytes'] / 1048576:.0f} MiB RSS, "
                f"roundtrip={entry.get('roundtrip_ok')}, queries={exact}/{len(queries)}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
