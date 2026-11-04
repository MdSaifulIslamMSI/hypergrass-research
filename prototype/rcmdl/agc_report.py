#!/usr/bin/env python3
"""Summarize an AGC per-target configuration sweep against RC-MDL LOO results.

The sweep directory contains one results-NN.json per target, each holding the
marginals for every swept configuration.  This script reports the best single
configuration (the fairest deployable AGC number) and the per-target oracle (a
lower bound on AGC's size), both against the same RC-MDL targets.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path


def rc_bytes(record: dict) -> int:
    if "loo" in record:
        return record["loo"]["archive_bytes"]
    return record["archive_bytes"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sweep-dir", type=Path, required=True)
    parser.add_argument("--rc-results", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    rc = json.loads(args.rc_results.read_text())
    sweeps = sorted(args.sweep_dir.glob("results-*.json"))
    if len(sweeps) != len(rc):
        raise SystemExit("sweep has %d targets, RC has %d" % (len(sweeps), len(rc)))

    totals: dict[str, int] = defaultdict(int)
    wins: Counter[str] = Counter()
    rows = []
    for record, path in zip(rc, sweeps):
        sweep = json.loads(path.read_text())
        best = min(sweep, key=lambda row: row["marginal_bytes"])
        for row in sweep:
            totals[row["config"]] += row["marginal_bytes"]
        wins[best["config"]] += 1
        charged = rc_bytes(record)
        rows.append(
            {
                "target": record["target"],
                "target_bytes": record["target_bytes"],
                "rc_bytes": charged,
                "best_config": best["config"],
                "best_bytes": best["marginal_bytes"],
                "ratio": best["marginal_bytes"] / charged,
            }
        )

    best_config = min(totals, key=totals.get)
    target_bases = sum(row["target_bytes"] for row in rows)
    rc_total = sum(row["rc_bytes"] for row in rows)
    oracle_total = sum(row["best_bytes"] for row in rows)
    report = {
        "targets": len(rows),
        "target_bases": target_bases,
        "rc_bytes": rc_total,
        "best_single": {
            "config": best_config,
            "bytes": totals[best_config],
            "bpb": totals[best_config] * 8.0 / target_bases,
            "ratio": totals[best_config] / rc_total,
        },
        "oracle": {
            "bytes": oracle_total,
            "bpb": oracle_total * 8.0 / target_bases,
            "ratio": oracle_total / rc_total,
        },
        "config_totals": dict(sorted(totals.items(), key=lambda item: item[1])),
        "config_wins": dict(wins.most_common()),
        "targets_detail": rows,
    }
    print(json.dumps({k: v for k, v in report.items() if k != "targets_detail"}, indent=2))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
