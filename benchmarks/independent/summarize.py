"""Render independent-pangenome reports as Markdown."""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path


RUNS = Path(__file__).resolve().parent / "runs"


def latest_report(dataset: str) -> Path | None:
    candidates = sorted(
        (RUNS / dataset).glob("*/report.json"),
        key=lambda path: path.stat().st_mtime,
    )
    return candidates[-1] if candidates else None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", action="append", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    datasets = args.dataset or sorted(
        path.name for path in RUNS.iterdir() if path.is_dir()
    )
    lines = [
        "| Dataset | Profile | Archive | Ratio | Bits/base | Wall | Peak RSS | Exact | Δ vs xz |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: | ---: |",
    ]
    deltas = []
    for dataset in datasets:
        report_path = latest_report(dataset)
        if report_path is None:
            continue
        report = json.loads(report_path.read_text(encoding="utf-8"))
        results = {
            entry["profile"]: entry
            for entry in report.get("results", [])
            if entry.get("status") == "success"
        }
        xz = results.get("xz-9")
        xz_bytes = xz.get("artifact_bytes") if xz else None
        for profile in ("hypergrass", "hypergrass-pangenome", "gzip-9", "xz-9"):
            entry = results.get(profile)
            if not entry:
                continue
            input_bytes = report["manifest"]["fasta"]["size_bytes"]
            input_bases = report["manifest"]["total_bases"]
            delta = ""
            if xz_bytes and profile != "xz-9":
                percent = (entry["artifact_bytes"] - xz_bytes) / xz_bytes * 100
                delta = f"{percent:+.2f}%"
                if profile == "hypergrass-pangenome":
                    deltas.append(percent)
            lines.append(
                "| {dataset} | {profile} | {size:,} B | {ratio:.2f}x | "
                "{bpb:.3f} | {wall:.1f} s | {rss:.0f} MiB | {exact} | {delta} |".format(
                    dataset=dataset,
                    profile=profile,
                    size=entry["artifact_bytes"],
                    ratio=input_bytes / entry["artifact_bytes"],
                    bpb=entry["artifact_bytes"] * 8 / input_bases,
                    wall=entry["wall_seconds"],
                    rss=entry["peak_rss_bytes"] / 1048576,
                    exact="yes" if entry.get("roundtrip_ok") else "no",
                    delta=delta or "—",
                )
            )
    if deltas:
        lines.extend(
            [
                "",
                f"Tuned profile vs xz: median {statistics.median(deltas):+.2f}% "
                f"across {len(deltas)} independent pangenomes.",
            ]
        )
    text = "\n".join(lines) + "\n"
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
