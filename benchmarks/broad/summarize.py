"""Summarize the broad cross-taxon and diversity-gradient runs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent.parent
RUNS = ROOT / "runs"
MATRIX = ROOT / "matrix.json"

HGS = "hypergrass-pangenome"
HGS_BOUNDED = "hypergrass-bounded-256m-tuned"
BASELINES = ("jarvis3-l7", "xz-9", "geco3-l6")
PROFILE_ORDER = (
    HGS,
    HGS_BOUNDED,
    "jarvis3-l7",
    "geco3-l6",
    "xz-9",
    "grass",
)


def load_report(run: dict) -> dict | None:
    path = RUNS / run["run_id"] / "report.json"
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def successful(report: dict | None) -> dict[str, dict]:
    if report is None:
        return {}
    return {
        entry["profile"]: entry
        for entry in report.get("results", [])
        if (
            entry.get("status") == "success"
            and entry.get("artifact_bytes")
            and entry.get("roundtrip_ok")
            and (entry.get("verification") or {}).get("verdict") == "pass"
            and (entry.get("verification") or {}).get("byte_identical") is True
        )
    }


def ratio(size: int, input_bytes: int) -> float:
    return input_bytes / size if size else 0.0


def delta(value: int, baseline: int) -> float:
    return (value - baseline) / baseline * 100.0


def exact(entry: dict) -> str:
    return "yes" if entry.get("roundtrip_ok") else "no"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(MATRIX))
    parser.add_argument("--out", default=str(ROOT / "summary.md"))
    parser.add_argument("--json", default=str(ROOT / "summary.json"))
    args = parser.parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    matrix = json.loads(Path(args.matrix).read_text(encoding="utf-8"))
    rows = []
    for run in matrix["runs"]:
        report = load_report(run)
        results = successful(report)
        input_bytes = (
            report.get("input", {}).get("sequence", {}).get("size_bytes")
            if report else None
        )
        input_bases = (
            report.get("input", {}).get("sequence", {}).get("bases")
            if report else None
        )
        rows.append(
            {
                "run": run,
                "report": report,
                "results": results,
                "input_bytes": input_bytes,
                "input_bases": input_bases,
            }
        )

    lines = ["# Broad pangenome experiment", ""]
    lines.extend(
        [
            "Diversity is the mean pairwise Mash distance from canonical 21-mer "
            "Jaccard; ANI is the Mash approximation. Lower archive and higher "
            "diversity are independent axes.",
            "",
        ]
    )

    cross = [row for row in rows if row["run"]["panel"] == "cross-taxon"]
    if cross:
        lines.extend(["## Cross-taxon natural collections", ""])
        lines.append(
            "Positive `HGS vs best` means the HyperGraSS archive is smaller than "
            "the smallest baseline archive on that row."
        )
        lines.append("")
        lines.append(
            "| Dataset | Organism | Records | Mb | Mash d | ANI | HGS | JARVIS3 | "
            "xz -9 | GeCo3 l6 | Best baseline | HGS vs best | Exact |"
        )
        lines.append(
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |"
        )
        for row in sorted(cross, key=lambda item: item["run"]["dataset"]):
            run = row["run"]
            results = row["results"]
            hgs = results.get(HGS)
            jarvis = results.get("jarvis3-l7")
            xz = results.get("xz-9")
            geco = results.get("geco3-l6")
            baselines = {
                name: results[name]["artifact_bytes"]
                for name in BASELINES if name in results
            }
            best_name = min(baselines, key=baselines.get) if baselines else None
            best_bytes = baselines[best_name] if best_name else None
            advantage = (
                (best_bytes - hgs["artifact_bytes"]) / best_bytes * 100.0
                if hgs and best_bytes else None
            )
            def size(entry):
                return f"{entry['artifact_bytes']:,}" if entry else "—"
            def bpb(entry):
                return f"{entry['bits_per_base']:.3f}" if entry else "—"
            lines.append(
                "| {dataset} | {organism} | {records} | {mb:.0f} | {distance:.5f} | "
                "{ani:.4f} | {hgs} | {jarvis} | {xz} | {geco} | {best} | {adv} | {exact} |".format(
                    dataset=run["dataset"],
                    organism=run["organism"],
                    records=run["diversity"]["record_count"],
                    mb=run["diversity"]["total_bases"] / 1e6,
                    distance=run["diversity"]["mean_pairwise_distance"],
                    ani=run["diversity"]["mean_mash_ani"],
                    hgs=size(hgs),
                    jarvis=size(jarvis),
                    xz=size(xz),
                    geco=size(geco),
                    best=best_name or "—",
                    adv=f"{advantage:+.2f}%" if advantage is not None else "—",
                    exact=exact(hgs) if hgs else "—",
                )
            )

        lines.extend(["", "### Profile detail", ""])
        lines.append(
            "| Dataset | Profile | Archive | bits/base | Ratio | Wall | CPU | Peak RSS | Exact | Δ vs JARVIS3 |"
        )
        lines.append(
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: |"
        )
        for row in sorted(cross, key=lambda item: item["run"]["dataset"]):
            run = row["run"]
            jarvis = row["results"].get("jarvis3-l7")
            for profile in PROFILE_ORDER:
                entry = row["results"].get(profile)
                if not entry:
                    continue
                change = (
                    delta(entry["artifact_bytes"], jarvis["artifact_bytes"])
                    if jarvis else None
                )
                lines.append(
                    "| {dataset} | {profile} | {size:,} B | {bpb:.3f} | {ratio:.2f}x | "
                    "{wall:.1f} s | {cpu:.1f} s | {rss:.0f} MiB | {exact} | {change} |".format(
                        dataset=run["dataset"],
                        profile=profile,
                        size=entry["artifact_bytes"],
                        bpb=entry["bits_per_base"],
                        ratio=ratio(entry["artifact_bytes"], row["input_bytes"]),
                        wall=entry["wall_seconds"],
                        cpu=entry.get("cpu_seconds") or 0.0,
                        rss=entry["peak_rss_bytes"] / 1048576,
                        exact=exact(entry),
                        change=f"{change:+.2f}%" if change is not None else "—",
                    )
                )

    diversity = [row for row in rows if row["run"]["panel"] == "diversity"]
    bounded_rows = []
    for row in sorted(rows, key=lambda item: item["run"]["run_id"]):
        reference = HGS
        hgs = row["results"].get(reference)
        if not hgs:
            continue
        raw = next(
            (
                entry
                for entry in (row["report"] or {}).get("results", [])
                if entry.get("profile") == HGS_BOUNDED
            ),
            None,
        )
        if raw is None:
            continue
        if raw.get("status") != "success" or not raw.get("roundtrip_ok"):
            bounded_rows.append(
                {
                    "run_id": row["run"]["run_id"],
                    "reference": reference,
                    "hgs": hgs,
                    "bounded": None,
                    "absolute_growth": None,
                    "growth": None,
                    "time_ratio": None,
                    "rss": (raw.get("peak_rss_bytes") or 0) / 1048576.0,
                    "passed": False,
                    "failure": raw.get("status") or "failed",
                }
            )
            continue
        bounded = row["results"][HGS_BOUNDED]
        growth = (bounded["artifact_bytes"] - hgs["artifact_bytes"]) / hgs["artifact_bytes"] * 100.0
        time_ratio = bounded["wall_seconds"] / hgs["wall_seconds"] if hgs["wall_seconds"] else float("inf")
        restore_rss = (bounded.get("decompress") or {}).get("peak_rss_bytes") or 0
        rss = max(bounded["peak_rss_bytes"], restore_rss) / 1048576.0
        bounded_rows.append(
            {
                "run_id": row["run"]["run_id"],
                "reference": reference,
                "hgs": hgs,
                "bounded": bounded,
                "absolute_growth": bounded["artifact_bytes"] - hgs["artifact_bytes"],
                "growth": growth,
                "time_ratio": time_ratio,
                "rss": rss,
                "passed": rss <= 256.0 and growth <= 2.0 and time_ratio <= 2.0,
                "failure": None,
            }
        )
    if bounded_rows:
        lines.extend(["", "### Bounded-memory gate (≤256 MiB, ≤2% growth, ≤2× wall)", ""])
        lines.append(
            "Reference is the unbounded `hypergrass-pangenome` archive for the "
            "same input."
        )
        lines.append("")
        lines.append("| Run | HGS archive | Bounded archive | Δ bytes | Growth | Bounded peak RSS | Wall ratio | Verdict |")
        lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |")
        for item in bounded_rows:
            if item["failure"]:
                lines.append(
                    "| {run_id} | {hgs:,} B | — | — | — | {rss:.0f} MiB | — | "
                    "fail ({failure}) |".format(
                        run_id=item["run_id"],
                        hgs=item["hgs"]["artifact_bytes"],
                        rss=item["rss"],
                        failure=item["failure"],
                    )
                )
                continue
            lines.append(
                "| {run_id} | {hgs:,} B | {bounded:,} B | {absolute:+,} B | {growth:+.2f}% | "
                "{rss:.0f} MiB | {ratio:.2f}x | {verdict} |".format(
                    run_id=item["run_id"],
                    hgs=item["hgs"]["artifact_bytes"],
                    bounded=item["bounded"]["artifact_bytes"],
                    absolute=item["absolute_growth"],
                    growth=item["growth"],
                    rss=item["rss"],
                    ratio=item["time_ratio"],
                    verdict="pass" if item["passed"] else "fail",
                )
            )

    if diversity:
        lines.extend(["", "## Diversity gradient", ""])
        lines.append(
            "| Dataset | Stratum | Mean d | ANI | HGS | JARVIS3 | xz -9 | "
            "HGS vs JARVIS3 | HGS vs xz | Exact |"
        )
        lines.append(
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |"
        )
        for row in sorted(
            diversity,
            key=lambda item: (item["run"]["dataset"], ("low", "medium", "high").index(item["run"]["stratum"])),
        ):
            run = row["run"]
            results = row["results"]
            hgs = results.get(HGS)
            jarvis = results.get("jarvis3-l7")
            xz = results.get("xz-9")
            def change(entry):
                if not hgs or not entry:
                    return "—"
                return f"{(entry['artifact_bytes'] - hgs['artifact_bytes']) / entry['artifact_bytes'] * 100.0:+.2f}%"
            lines.append(
                "| {dataset} | {stratum} | {distance:.5f} | {ani:.4f} | {hgs} | {jarvis} | {xz} | {dj} | {dx} | {exact} |".format(
                    dataset=run["dataset"],
                    stratum=run["stratum"],
                    distance=run["diversity"]["mean_pairwise_distance"],
                    ani=run["diversity"]["mean_mash_ani"],
                    hgs=f"{hgs['artifact_bytes']:,}" if hgs else "—",
                    jarvis=f"{jarvis['artifact_bytes']:,}" if jarvis else "—",
                    xz=f"{xz['artifact_bytes']:,}" if xz else "—",
                    dj=change(jarvis),
                    dx=change(xz),
                    exact=exact(hgs) if hgs else "—",
                )
            )

        lines.extend(["", "### Trend", ""])
        lines.append(
            "Advantage is `(baseline - HGS) / baseline`; positive means HyperGraSS "
            "is smaller. Spearman is over three strata only and is indicative, "
            "not a significance test."
        )
        lines.append("")
        lines.append("| Dataset | Low | Medium | High | Spearman(d, advantage; n=3) |")
        lines.append("| --- | ---: | ---: | ---: | ---: |")
        for dataset in sorted({row["run"]["dataset"] for row in diversity}):
            points = []
            advantages = {}
            for stratum in ("low", "medium", "high"):
                row = next(
                    (item for item in diversity
                     if item["run"]["dataset"] == dataset and item["run"]["stratum"] == stratum),
                    None,
                )
                if not row or HGS not in row["results"] or "jarvis3-l7" not in row["results"]:
                    continue
                hgs = row["results"][HGS]["artifact_bytes"]
                jarvis = row["results"]["jarvis3-l7"]["artifact_bytes"]
                advantage = (jarvis - hgs) / jarvis * 100.0
                points.append((row["run"]["diversity"]["mean_pairwise_distance"], advantage))
                advantages[stratum] = advantage
            correlation = "—"
            if len(points) >= 3:
                rho, _ = spearmanr([point[0] for point in points], [point[1] for point in points])
                correlation = f"{rho:+.2f}"
            lines.append(
                "| {dataset} | {low} | {medium} | {high} | {rho} |".format(
                    dataset=dataset,
                    low=f"{advantages['low']:+.2f}%" if "low" in advantages else "—",
                    medium=f"{advantages['medium']:+.2f}%" if "medium" in advantages else "—",
                    high=f"{advantages['high']:+.2f}%" if "high" in advantages else "—",
                    rho=correlation,
                )
            )

    failures = []
    for row in rows:
        if row["report"] is None:
            failures.append(f"{row['run']['run_id']}: no report")
            continue
        for entry in row["report"].get("results", []):
            verification = entry.get("verification") or {}
            if (
                entry.get("status") != "success"
                or not entry.get("roundtrip_ok")
                or verification.get("verdict") != "pass"
                or verification.get("byte_identical") is not True
            ):
                failures.append(
                    f"{row['run']['run_id']}/{entry.get('profile')}: "
                    f"status={entry.get('status')} exact={entry.get('roundtrip_ok')} "
                    f"verification={verification.get('verdict')} "
                    f"error={entry.get('error')}"
                )
    lines.extend(["", "## Coverage", ""])
    lines.append(f"- Runs in matrix: {len(rows)}")
    lines.append(f"- Runs with a report: {sum(1 for row in rows if row['report'])}")
    lines.append(f"- Failed profile/run entries: {len(failures)}")
    if failures:
        lines.extend(["", "```"])
        lines.extend(failures)
        lines.append("```")

    text = "\n".join(lines) + "\n"
    Path(args.out).write_text(text, encoding="utf-8")
    payload = {
        "schema_version": 1,
        "runs": [
            {
                "run_id": row["run"]["run_id"],
                "panel": row["run"]["panel"],
                "dataset": row["run"]["dataset"],
                "organism": row["run"]["organism"],
                "stratum": row["run"]["stratum"],
                "diversity": row["run"]["diversity"],
                "input_bytes": row["input_bytes"],
                "input_bases": row["input_bases"],
                "profiles": row["results"],
            }
            for row in rows
        ],
        "failures": failures,
    }
    Path(args.json).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(text, end="")
    print(f"wrote {args.out} and {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
