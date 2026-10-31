"""Render docs/END-TO-END.md tables from an end-to-end query report.

The prose is human-written; the tables come from the report JSON. Repeat runs
replace the `E2E-TABLES:START` / `E2E-TABLES:END` block.

Usage:
    python -B benchmarks/endtoend/render.py \
        --report benchmarks/endtoend/runs/endtoend-1/report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent.parent


def fmt_int(value) -> str:
    return f"{value:,}" if isinstance(value, int) else "?"


def fmt_seconds(value, digits: int = 2) -> str:
    return f"{value:.{digits}f}" if isinstance(value, (int, float)) else "?"


def fmt_mib(value) -> str:
    return f"{value / 1048576:.0f}" if isinstance(value, (int, float)) else "?"


def query_records(profiles: list[dict]) -> str:
    for profile in profiles:
        if profile.get("queries"):
            records = [query["record"] for query in profile["queries"]]
            return ", ".join(str(record) for record in records)
    return "?"


AXIS_LABELS = {
    "archive_bytes": "archive bytes",
    "query_wall_median_seconds": "query median",
    "restore_wall_seconds": "restore time",
    "query_peak_rss_bytes": "query peak RSS",
}


def tables(report: dict) -> str:
    workload = report["workload"]
    profiles = report["profiles"]
    lines = ["### Workload\n"]
    lines.append("| Item | Value |")
    lines.append("| --- | --- |")
    lines.append(
        "| Input | `{path}`, {size} bytes, {records} records |".format(
            path=report["input"]["path"].replace("\\", "/"),
            size=fmt_int(report["input"]["size_bytes"]),
            records=fmt_int(report["input"]["records"]),
        )
    )
    lines.append(
        "| Query | {count} slices of {length} bases, one fresh process each |".format(
            count=workload["query_count"], length=fmt_int(workload["query_length"])
        )
    )
    lines.append(f"| Records queried | {query_records(profiles)} |")
    lines.append("| Restore | one full decompress per profile, from the scale report |")
    lines.append("")

    lines.append("### Results\n")
    lines.append(
        "| Profile | Archive bytes | Query median (s) | Query p95 (s) | "
        "Query max (s) | Query total (s) | Query peak RSS (MiB) | Restore (s) | "
        "Restore peak RSS (MiB) | Exact |"
    )
    lines.append(
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"
    )
    for profile in profiles:
        exact = f"{profile['query_exact']}/{len(profile['queries'])}"
        lines.append(
            "| {id} | {archive} | {median} | {p95} | {max} | {total} | {query_rss} | "
            "{restore} | {restore_rss} | {exact} |".format(
                id=profile["id"],
                archive=fmt_int(profile["archive_bytes"]),
                median=fmt_seconds(profile["query_wall_median_seconds"]),
                p95=fmt_seconds(profile["query_wall_p95_seconds"]),
                max=fmt_seconds(profile["query_wall_max_seconds"]),
                total=fmt_seconds(profile["query_total_wall_seconds"], 1),
                query_rss=fmt_mib(profile["query_peak_rss_bytes"]),
                restore=fmt_seconds(profile["restore_wall_seconds"], 1),
                restore_rss=fmt_mib(profile["restore_peak_rss_bytes"]),
                exact=exact,
            )
        )
    lines.append("")
    lines.append(
        "With 13 samples, the p95 column equals the maximum. Query time varies "
        "strongly by record and profile; the median is not a constant-time "
        "claim.\n"
    )

    analysis = report["analysis"]
    lines.append("### Query speedup\n")
    lines.append(
        "Baseline query median divided by HyperGraSS query median. Above 1.0 "
        "means HyperGraSS answered the same slice faster.\n"
    )
    lines.append("| Comparison | Speedup |")
    lines.append("| --- | ---: |")
    for name, value in analysis["query_speedup"].items():
        lines.append(f"| {name} | {value:.2f}x |")
    lines.append("")

    lines.append("### Break-even query count\n")
    lines.append(
        "Queries needed for `one restore + Q queries` to favour HyperGraSS on "
        "wall time: "
        "`ceil((HGS restore - baseline restore) / "
        "(baseline query - HGS query))`. `never` means the baseline has no "
        "per-query advantage to recover the restore cost.\n"
    )
    lines.append("| Comparison | Break-even queries |")
    lines.append("| --- | ---: |")
    for name, value in analysis["breakeven_queries_one_restore"].items():
        if value is None:
            cell = "never"
        elif value == 0:
            cell = "0 (ahead before queries)"
        else:
            cell = f"{value:,}"
        lines.append(f"| {name} | {cell} |")
    lines.append("")

    axes = ", ".join(
        AXIS_LABELS.get(axis, axis.replace("_", " ")) for axis in analysis["pareto_axes"]
    )
    lines.append("### Pareto frontier\n")
    lines.append(
        "Non-dominated profiles across {axes}: **{frontier}**.\n".format(
            axes=axes,
            frontier=", ".join(f"`{name}`" for name in analysis["pareto_frontier"]),
        )
    )
    return "\n".join(lines)


def replace_generated(text: str, content: str) -> str:
    start = "<!-- E2E-TABLES:START -->"
    end = "<!-- E2E-TABLES:END -->"
    if start not in text or end not in text:
        raise SystemExit("missing E2E-TABLES markers in the document")
    head_end = text.index(start) + len(start)
    return text[:head_end] + "\n" + content + "\n" + text[text.index(end):]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--doc", default="docs/END-TO-END.md")
    args = parser.parse_args(argv)

    report = json.loads(Path(args.report).read_text(encoding="utf-8"))
    if report.get("status") != "success":
        raise SystemExit(f"report is not finished: {report.get('status')}")
    doc = Path(args.doc)
    text = doc.read_text(encoding="utf-8")
    doc.write_text(replace_generated(text, tables(report)), encoding="utf-8")
    print(f"rendered {doc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
