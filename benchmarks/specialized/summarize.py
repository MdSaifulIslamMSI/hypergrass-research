"""Render specialized-compressor reports as Markdown."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


RUNS = Path(__file__).resolve().parent / "runs"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args(argv)

    report_path = RUNS / args.run / "report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    sequence = report["input"]["sequence"]
    lines = [
        "| Profile | Archive | Bits/base | Compress | Decompress | Peak RSS | Exact |",
        "| --- | ---: | ---: | ---: | ---: | ---: | :---: |",
    ]
    for entry in report.get("results", []):
        if entry.get("status") != "success" or not entry.get("artifact_bytes"):
            lines.append(
                f"| {entry['profile']} | — | — | — | — | — | "
                f"failed: {entry.get('status')} |"
            )
            continue
        decompress = entry.get("decompress") or {}
        lines.append(
            "| {profile} | {size:,} B | {bpb:.3f} | {compress:.1f} s | "
            "{decompress:.1f} s | {rss:.0f} MiB | {exact} |".format(
                profile=entry["profile"],
                size=entry["artifact_bytes"],
                bpb=entry["artifact_bytes"] * 8 / sequence["bases"],
                compress=entry["wall_seconds"],
                decompress=decompress.get("wall_seconds") or 0.0,
                rss=(entry.get("peak_rss_bytes") or 0) / 1048576,
                exact="yes" if entry.get("roundtrip_ok") else "no",
            )
        )
    text = "\n".join(lines) + "\n"
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
