"""Render measured benchmark JSON as Markdown tables.

Reads a corpus report and/or a pangenome scale report and prints a Markdown
document. No numbers are typed by hand; every cell comes from the JSON.

Usage:
    python -B benchmarks/summarize.py --corpus benchmarks/corpus/runs/corpus-full-1/report.json
    python -B benchmarks/summarize.py --corpus ... --scale benchmarks/pangenome/runs/.../report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "corpus"))
from audit import EXPECTED_REPS as CORPUS_EXPECTED_REPS  # noqa: E402

CORPUS_BASELINES = ("gzip-9", "xz-9", "jarvis3-l7", "jarvis3-l26", "geco3-l9", "geco3-l6", "geco3-l3")
CORPUS_HGS = ("hypergrass", "hypergrass-optimize", "hypergrass-lean")
# `running` is what a terminated runner leaves behind: the results already
# appended are still evidence, but the matrix is incomplete.
CORPUS_ACCEPTED = ("success", "partial", "running", "stopped")


def corpus_coverage(report):
    """Return (present, planned, complete_datasets, total_datasets)."""
    planned = sum(CORPUS_EXPECTED_REPS.values()) * len(report["protocol"]["datasets"])
    present = len(report["results"])
    by_cell = {}
    for entry in report["results"]:
        by_cell.setdefault((entry["dataset"], entry["profile"]), []).append(entry)
    complete = 0
    for dataset in report["protocol"]["datasets"]:
        if all(
            len(by_cell.get((dataset["id"], profile), [])) == reps
            for profile, reps in CORPUS_EXPECTED_REPS.items()
        ):
            complete += 1
    return present, planned, complete, len(report["protocol"]["datasets"])


def corpus_gaps(report):
    """Cells that started but have fewer results than the frozen matrix."""
    counts = {}
    for entry in report["results"]:
        key = (entry["dataset"], entry["profile"])
        counts[key] = counts.get(key, 0) + 1
    gaps = []
    for dataset in report["protocol"]["datasets"]:
        for profile, reps in CORPUS_EXPECTED_REPS.items():
            got = counts.get((dataset["id"], profile), 0)
            if 0 < got < reps:
                gaps.append((dataset["id"], profile, got, reps))
    return gaps


def fmt_int(value):
    return f"{value:,}" if isinstance(value, int) else "—"


def fmt_float(value, digits=2):
    return f"{value:.{digits}f}" if isinstance(value, (int, float)) else "—"


def row_by(summary, dataset, profile):
    for row in summary:
        if row["dataset"] == dataset and row["profile"] == profile:
            return row
    return None


def corpus_tables(report: dict) -> list[str]:
    lines = []
    status = report.get("status")
    if status not in ("success", "partial"):
        present, planned, complete, total = corpus_coverage(report)
        lines.append(
            f"> Coverage: {present:,} of {planned:,} planned results, "
            f"{complete} of {total} datasets complete. The report status is "
            f"`{status}` because the runner was terminated before writing a "
            f"final status. The tables below are the completed evidence only; "
            f"they are not a full {total}-dataset matrix.\n"
        )
        gaps = corpus_gaps(report)
        if gaps:
            lines.append("| Cell with fewer than planned repetitions | Results | Planned |")
            lines.append("| --- | ---: | ---: |")
            for dataset, profile, got, want in gaps:
                lines.append(f"| {dataset} / {profile} | {got} | {want} |")
            lines.append("")
    summary = report["summary"]
    datasets = []
    for row in summary:
        if row["dataset"] not in datasets:
            datasets.append(row["dataset"])
    profiles = []
    for row in summary:
        if row["profile"] not in profiles:
            profiles.append(row["profile"])

    lines.append("### Head-to-head, frozen corpus\n")
    lines.append(
        "| Dataset | Input bytes | HyperGraSS | `--optimize` | `--lean` | Best completed baseline | Its bytes | Δ smaller than best (+ = smaller) |\n"
        "| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |"
    )
    for dataset in datasets:
        hgs = row_by(summary, dataset, "hypergrass")
        optimize = row_by(summary, dataset, "hypergrass-optimize")
        lean = row_by(summary, dataset, "hypergrass-lean")
        candidates = []
        for name in CORPUS_BASELINES:
            row = row_by(summary, dataset, name)
            if row and row["status"] == "success" and row.get("artifact_bytes"):
                candidates.append(row)
        best = min(candidates, key=lambda r: r["artifact_bytes"]) if candidates else None
        delta = "—"
        if best and hgs and hgs.get("artifact_bytes"):
            delta = f"{(best['artifact_bytes'] - hgs['artifact_bytes']) / best['artifact_bytes'] * 100:+.2f}%"
        lines.append(
            "| {ds} | {inp} | {hgs} | {opt} | {lean} | {bn} | {bb} | {delta} |".format(
                ds=dataset,
                inp=fmt_int(hgs["input_bytes"] if hgs else None),
                hgs=fmt_int(hgs.get("artifact_bytes") if hgs else None),
                opt=fmt_int(optimize.get("artifact_bytes") if optimize else None),
                lean=fmt_int(lean.get("artifact_bytes") if lean else None),
                bn=best["profile"] if best else "—",
                bb=fmt_int(best["artifact_bytes"]) if best else "—",
                delta=delta,
            )
        )
    lines.append("")

    lines.append("### Effective bits per base, frozen corpus\n")
    header = "| Dataset | " + " | ".join(profiles) + " |"
    divider = "| --- | " + " | ".join("---:" for _ in profiles) + " |"
    lines.append(header)
    lines.append(divider)
    for dataset in datasets:
        cells = []
        for profile in profiles:
            row = row_by(summary, dataset, profile)
            if row and row["status"] == "success" and row.get("bits_per_base") is not None:
                cells.append(f"{row['bits_per_base']:.3f}")
            elif row and row["status"] == "not_applicable":
                cells.append("n/a")
            elif row and row["profile"] in ("geco3-l9", "jarvis3-l26") and row["status"] in ("timeout", "failed"):
                cells.append("infeasible")
            elif row and row["status"] == "timeout":
                cells.append("timeout")
            elif row:
                cells.append("failed")
            else:
                cells.append("—")
        lines.append(f"| {dataset} | " + " | ".join(cells) + " |")
    lines.append("")
    lines.append(
        "`infeasible` marks the high-effort `geco3-l9` and `jarvis3-l26` cells that "
        "exceeded the 180 s operation cap (compression or decompression) or the "
        "container memory limit; they are not treated as zero-byte results.\n"
    )

    lines.append("### Median of per-dataset medians: compression time and peak RSS\n")
    lines.append(
        "| Profile | Median wall (s) | Median CPU (s) | Median peak RSS (MiB) | Successful datasets |\n"
        "| --- | ---: | ---: | ---: | ---: |"
    )
    for profile in profiles:
        rows = [
            row for row in summary
            if row["profile"] == profile and row["status"] == "success" and row.get("median_wall_seconds")
        ]
        if not rows:
            lines.append(f"| {profile} | — | — | — | 0 |")
            continue
        wall = sorted(row["median_wall_seconds"] for row in rows)
        cpu = sorted(row["median_cpu_seconds"] for row in rows if row.get("median_cpu_seconds") is not None)
        rss = sorted(row["median_peak_rss_bytes"] for row in rows if row.get("median_peak_rss_bytes") is not None)
        mid = lambda values: values[len(values) // 2] if values else None
        lines.append(
            "| {p} | {w} | {c} | {r} | {n} |".format(
                p=profile,
                w=fmt_float(mid(wall)),
                c=fmt_float(mid(cpu)),
                r=fmt_float(mid(rss) / 1048576 if mid(rss) else None, 1),
                n=len(rows),
            )
        )
    lines.append("")
    lines.append(
        "Cells with fewer than the planned repetitions contribute their "
        "available samples; the coverage note above lists them.\n"
    )

    lines.append("### Non-success cells\n")
    anomalies = [
        row for row in summary
        if row["status"] not in ("success", "not_applicable")
    ]
    na = sum(1 for row in summary if row["status"] == "not_applicable")
    lines.append(
        f"{len(summary)} dataset/profile cells: "
        f"{sum(1 for row in summary if row['status'] == 'success')} successful, "
        f"{na} not applicable, {len(anomalies)} infeasible or failed.\n"
    )
    if anomalies:
        lines.append("| Dataset | Profile | Status | Attempts | Successful |\n| --- | --- | --- | ---: | ---: |")
        for row in anomalies:
            lines.append(
                "| {ds} | {p} | {s} | {a} | {ok} |".format(
                    ds=row["dataset"],
                    p=row["profile"],
                    s=row["status"],
                    a=row["attempts"],
                    ok=row["successful"],
                )
            )
        lines.append("")
    return lines


def scale_tables(report: dict) -> list[str]:
    lines = []
    manifest = report["manifest"]
    lines.append("### Scale results\n")
    lines.append(
        "Input: `{path}`, {size} bytes, {records} records, {bases} bases. "
        "SHA-256 `{sha}`.\n".format(
            path=manifest["fasta"]["path"],
            size=fmt_int(manifest["fasta"]["size_bytes"]),
            records=len(manifest["records"]),
            bases=fmt_int(manifest["total_bases"]),
            sha=manifest["fasta"]["sha256"],
        )
    )
    lines.append(
        "| Profile | Archive bytes | Ratio | Bits/byte | Compress wall (s) | Compress CPU (s) | Peak RSS (MiB) | Round trip | Queries |\n"
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |"
    )
    for entry in report["results"]:
        if entry["status"] == "success":
            queries = entry.get("queries", [])
            query_cell = f"{sum(1 for q in queries if q['exact'])}/{len(queries)} exact" if queries else "—"
            lines.append(
                "| {p} | {b} | {r:.2f}x | {bpb:.3f} | {w:.1f} | {c:.1f} | {rss:.1f} | yes | {q} |".format(
                    p=entry["profile"],
                    b=fmt_int(entry["artifact_bytes"]),
                    r=entry["ratio"],
                    bpb=entry["bits_per_byte"],
                    w=entry["wall_seconds"],
                    c=entry["cpu_seconds"],
                    rss=entry["peak_rss_bytes"] / 1048576,
                    q=query_cell,
                )
            )
        else:
            lines.append(f"| {entry['profile']} | — | — | — | — | — | — | **{entry['status']}** | — |")
    lines.append("")
    return lines


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", required=True)
    parser.add_argument("--scale")
    parser.add_argument("--out")
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    corpus = json.loads(Path(args.corpus).read_text(encoding="utf-8"))
    lines = []
    if corpus.get("status") not in CORPUS_ACCEPTED:
        raise SystemExit(f"corpus report is not finished: {corpus.get('status')}")
    lines.extend(corpus_tables(corpus))
    if args.scale:
        scale = json.loads(Path(args.scale).read_text(encoding="utf-8"))
        if scale.get("status") not in ("success", "partial"):
            raise SystemExit(f"scale report is not finished: {scale.get('status')}")
        lines.extend(scale_tables(scale))
    text = "\n".join(lines)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    main()
