"""Run specialized DNA compressors on the >100 MB pangenome.

JARVIS3 and GeCo3 both document FASTA input and use only the A/C/G/T sequence,
so this harness feeds them the exact same multi-FASTA as the HyperGraSS scale
test. Their archives exclude metadata; HyperGraSS and xz archives include it.
That asymmetry slightly favors the specialized tools, which is the conservative
direction for a HyperGraSS comparison.

GraSS is the published ancestor of HyperGraSS. Its raw path accepts the flat
headerless sequence and restores it byte-for-byte, so it receives the identical
120,108,223-byte `input.seq` rather than a FASTA wrapper. Its archive is the
native `FinalBsc.bsc` produced by the bundled grammar/substitution pipeline and
`libbsc -e2`; the harness records the exact GraSS and libbsc commits.

A small prefix slice can be used as a pilot to estimate full-scale runtime
before committing hours to a 121.8 MB run.

Usage:
    python -B benchmarks/specialized/run.py \
        --run-id specialized-pilot-5m --slice-bases 5000000 --timeout 1800
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import sys


PROJECT = Path(__file__).resolve().parent.parent.parent
DATA = Path(__file__).resolve().parent / "data"
RUNS = Path(__file__).resolve().parent / "runs"
SOURCE = PROJECT / "benchmarks/pangenome/data/pangenome.fasta"

SPEC = importlib.util.spec_from_file_location(
    "pangenome_run", PROJECT / "benchmarks/pangenome/run.py"
)
pangenome_run = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(pangenome_run)

PROFILES = (
    {"id": "jarvis3-l7", "kind": "jarvis3", "level": "7"},
    {"id": "jarvis3-l26", "kind": "jarvis3", "level": "26"},
    {"id": "geco3-l3", "kind": "geco3", "level": "3"},
    {"id": "geco3-l6", "kind": "geco3", "level": "6"},
    {"id": "geco3-l9", "kind": "geco3", "level": "9"},
    {"id": "xz-9", "kind": "xz", "level": "9"},
    {
        "id": "hypergrass-pangenome",
        "kind": "hypergrass",
        "level": "consensus-anchor",
        "options": [
            "--pangenome", "--single-thread", "--anchor-single-page",
            "--block-mib", "64", "--anchor-mib", "9",
        ],
        "jvm": ["-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count"],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "hypergrass-bounded-256m",
        "kind": "hypergrass",
        "level": "bounded-256m",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-single-page",
            "--block-mib", "16", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-Dhypergrass.boundedVotes8=true",
        ],
        "heap": "192m",
        "memory": "256m",
    },
    {
        "id": "hypergrass-bounded-256m-tuned",
        "kind": "hypergrass",
        "level": "bounded-256m",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-single-page",
            "--block-mib", "16", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-XX:+UseSerialGC", "-XX:ReservedCodeCacheSize=32m",
            "-XX:MaxMetaspaceSize=64m", "-XX:CICompilerCount=2", "-Xss512k",
        ],
        "heap": "160m",
        "memory": "256m",
    },
    {
        "id": "hypergrass-bounded-256m-cm-pages",
        "kind": "hypergrass",
        "level": "bounded-256m-cm-pages",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-cm-pages",
            "--block-mib", "16", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-XX:+UseSerialGC", "-XX:ReservedCodeCacheSize=32m",
            "-XX:MaxMetaspaceSize=64m", "-XX:CICompilerCount=2", "-Xss512k",
        ],
        "heap": "160m",
        "memory": "256m",
    },
    {
        "id": "hypergrass-bounded-256m-stride8",
        "kind": "hypergrass",
        "level": "bounded-256m-stride8",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-single-page",
            "--seed-step", "8", "--block-mib", "16", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-XX:+UseSerialGC", "-XX:ReservedCodeCacheSize=32m",
            "-XX:MaxMetaspaceSize=64m", "-XX:CICompilerCount=2", "-Xss512k",
        ],
        "heap": "160m",
        "memory": "256m",
    },
    {
        "id": "hypergrass-bounded-144m-stride8",
        "kind": "hypergrass",
        "level": "bounded-144m-stride8",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-single-page",
            "--seed-step", "8", "--block-mib", "16", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-XX:+UseSerialGC", "-XX:ReservedCodeCacheSize=32m",
            "-XX:MaxMetaspaceSize=64m", "-XX:CICompilerCount=2", "-Xss512k",
        ],
        "heap": "144m",
        "memory": "256m",
    },
    {
        "id": "hypergrass-bounded-wide",
        "kind": "hypergrass",
        "level": "bounded-wide",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-single-page",
            "--seed-step", "8", "--block-mib", "32", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-XX:+UseSerialGC", "-XX:ReservedCodeCacheSize=32m",
            "-XX:MaxMetaspaceSize=64m", "-XX:CICompilerCount=2", "-Xss512k",
        ],
        "heap": "256m",
        "memory": "320m",
    },
    {
        "id": "hypergrass-bounded-wide-cm",
        "kind": "hypergrass",
        "level": "bounded-wide-cm",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-cm-pages",
            "--seed-step", "8", "--block-mib", "32", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-XX:+UseSerialGC", "-XX:ReservedCodeCacheSize=32m",
            "-XX:MaxMetaspaceSize=64m", "-XX:CICompilerCount=2", "-Xss512k",
        ],
        "heap": "256m",
        "memory": "320m",
    },
    {
        "id": "hypergrass-bounded-native",
        "kind": "hypergrass",
        "level": "bounded-256m-native",
        "options": [
            "--pangenome", "--bounded", "--single-thread", "--anchor-single-page",
            "--seed-step", "8", "--block-mib", "32", "--anchor-mib", "9",
            "--split-tail-block",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
        ],
        "runtime": "native",
        "binary": "bin/native/hypergrass",
        "heap": "192m",
        "memory": "256m",
    },
    # Codec-ablation profiles. These force one payload model at a time so the
    # four methods can be compared on identical input; the default path still
    # chooses the smallest representation automatically.
    {
        "id": "hypergrass-qcrm",
        "kind": "hypergrass",
        "level": "ablation-qcrm",
        "options": ["--lean"],
        "jvm": ["-Dhypergrass.dnaMode=1"],
        "heap": "1g",
        "memory": "1g",
    },
    {
        "id": "hypergrass-auto",
        "kind": "hypergrass",
        "level": "ablation-auto",
        "options": [
            "--pangenome", "--single-thread", "--anchor-single-page",
            "--block-mib", "64", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-Dhypergrass.noCompact=true",
        ],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "hypergrass-bclzg",
        "kind": "hypergrass",
        "level": "ablation-bclzg",
        "options": [
            "--pangenome", "--single-thread", "--anchor-single-page",
            "--block-mib", "64", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-Dhypergrass.dnaMode=2",
        ],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "hypergrass-achmc",
        "kind": "hypergrass",
        "level": "ablation-achmc",
        "options": [
            "--pangenome", "--single-thread", "--anchor-single-page",
            "--block-mib", "64", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-Dhypergrass.dnaMode=3",
        ],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "hypergrass-anchor-cm",
        "kind": "hypergrass",
        "level": "ablation-anchor-cm",
        "options": [
            "--pangenome", "--single-thread", "--anchor-cm-pages",
            "--block-mib", "64", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-Dhypergrass.noCompact=true",
        ],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "hypergrass-aux-no-lz",
        "kind": "hypergrass",
        "level": "ablation-aux-no-lz",
        "options": [
            "--pangenome", "--single-thread", "--anchor-single-page",
            "--block-mib", "64", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-Dhypergrass.lzAux=false",
            "-Dhypergrass.noCompact=true",
        ],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "hypergrass-aux-no-plane",
        "kind": "hypergrass",
        "level": "ablation-aux-no-plane",
        "options": [
            "--pangenome", "--single-thread", "--anchor-single-page",
            "--block-mib", "64", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-Dhypergrass.planeAux=false",
            "-Dhypergrass.noCompact=true",
        ],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "hypergrass-aux-basic",
        "kind": "hypergrass",
        "level": "ablation-aux-basic",
        "options": [
            "--pangenome", "--single-thread", "--anchor-single-page",
            "--block-mib", "64", "--anchor-mib", "9",
        ],
        "jvm": [
            "-Dhypergrass.matchMin=104", "-Dhypergrass.rangeModel=count",
            "-Dhypergrass.planeAux=false", "-Dhypergrass.lzAux=false",
            "-Dhypergrass.rangeAux=false",
            "-Dhypergrass.noCompact=true",
        ],
        "heap": "1536m",
        "memory": "2g",
    },
    {
        "id": "grass",
        "kind": "grass",
        "level": "default-e2",
        "heap": "1500m",
        "memory": "2g",
    },
)

# The GraSS paper build expects these exact upstream revisions. The checkout is
# a local, gitignored build tree; fetch_vendor.ps1 recreates it.
GRASS_COMMIT = "6d93700f5a2de8841f7cf2c20064413c80df2c4a"
LIBBSC_COMMIT = "1dc676540cb40a41562af6f2db04bd1adce97f90"

VERIFIER = Path(__file__).resolve().parent.parent / "verify" / "verify.py"


def independent_verification(
    original: Path,
    restored: Path,
    archive: Path,
    expected_sha256: str,
    report_path: Path,
) -> dict:
    """Run the standalone verifier as its own process (never imports the codec)."""
    command = [
        sys.executable, "-B", str(VERIFIER),
        "--original", str(original),
        "--restored", str(restored),
        "--archive", str(archive),
        "--expect-sha256", expected_sha256,
        "--report", str(report_path),
    ]
    run = subprocess.run(
        command, cwd=PROJECT, capture_output=True, text=True, timeout=7200,
    )
    if report_path.exists():
        return json.loads(report_path.read_text(encoding="utf-8"))
    return {
        "verdict": "fail",
        "problems": [f"verifier produced no report (exit {run.returncode})"],
        "stderr": run.stderr[-2000:],
    }


def slice_fasta(source: Path, target_bases: int, destination: Path) -> tuple[int, int, int]:
    """Copy sequence until target_bases have been written.

    The last record may be truncated, which keeps the slice close to the
    requested size even when individual genomes are several Mbp.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    records = 0
    bases = 0
    with source.open("rb") as stream, destination.open("wb") as out:
        for line in stream:
            if line.startswith(b">"):
                if bases >= target_bases:
                    break
                out.write(line)
                records += 1
                continue
            sequence = line.strip()
            remaining = target_bases - bases
            if len(sequence) > remaining:
                out.write(sequence[:remaining] + b"\n")
                bases += remaining
                break
            out.write(line)
            bases += len(sequence)
    if bases < target_bases:
        raise SystemExit(f"slice has only {bases:,} bases; target was {target_bases:,}")
    return records, bases, destination.stat().st_size


def sequence_only(fasta: Path, destination: Path) -> dict:
    """Extract one wrapped A/C/G/T sequence, masking rare ambiguity codes."""
    chunks = []
    with fasta.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                continue
            stripped = line.strip()
            if stripped:
                chunks.append(stripped)
    upper = b"".join(chunks).upper()
    counts = {base: upper.count(bytes([base])) for base in b"ACGT"}
    non_acgt = len(upper) - sum(counts.values())
    table = bytes((value if value in b"ACGT" else ord("A")) for value in range(256))
    masked = upper.translate(table)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(masked)
    sequence_digest = hashlib.sha256(masked).hexdigest()
    return {
        "path": str(destination.relative_to(PROJECT)),
        "size_bytes": len(masked),
        "bases": len(masked),
        "sha256": sequence_digest,
        "non_acgt_bases_masked": non_acgt,
        "base_counts": {chr(key): value for key, value in counts.items()},
    }


def sequence_sha256(path: Path) -> str:
    payload = path.read_bytes()
    return hashlib.sha256(payload.translate(None, b"\r\n\t ")).hexdigest()


def parse_input(path: Path) -> tuple[int, int]:
    """Return (records, bases) for a FASTA file without loading it into memory."""
    records = 0
    bases = 0
    with path.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                records += 1
            else:
                bases += len(line.strip())
    return records, bases


def _median(values: list[float]) -> float | None:
    usable = [value for value in values if value is not None]
    return statistics.median(usable) if usable else None


def roundtrip_failure(entry: dict) -> str:
    """Explain why a successful compression did not restore its input."""
    decompress = entry.get("decompress") or {}
    verdict = (entry.get("verification") or {}).get("verdict")
    if decompress.get("status") not in (None, "success"):
        return (
            decompress.get("error")
            or f"restore failed with status {decompress['status']}"
        )
    if verdict and verdict != "pass":
        problems = (entry.get("verification") or {}).get("problems") or []
        return "independent verification failed: " + (
            "; ".join(str(problem) for problem in problems)
            if problems else "no problems recorded"
        )
    if not entry.get("decompress"):
        return "restore was not attempted"
    return "restore did not reproduce the input"


def summarize_profile(profile: dict, observations: list[dict], expected: int) -> dict:
    """Collapse repeated observations into one auditable profile summary.

    Size must be identical across repeats; otherwise the archive is not a
    deterministic function of the input. Wall time uses the median, while RSS
    uses the worst case, which is the value a memory gate must respect.
    """
    summary = {
        "profile": profile["id"],
        "kind": profile["kind"],
        "level": profile["level"],
        "config": {
            key: profile[key]
            for key in ("options", "jvm", "heap", "memory", "runtime", "binary")
            if key in profile
        },
        "planned_observations": expected,
        "status": "failed",
        "roundtrip_ok": False,
        "deterministic": False,
        "artifact_bytes": None,
        "artifact_sha256": None,
        "bits_per_base": None,
        "wall_seconds": None,
        "wall_seconds_min": None,
        "wall_seconds_max": None,
        "cpu_seconds": None,
        "peak_rss_bytes": None,
        "peak_rss_median_bytes": None,
        "cgroup_memory_peak_bytes": None,
        "cgroup_oom_kills": None,
        "decompress": None,
        "error": None,
        "observations": observations,
    }
    if observations:
        summary["input_bytes"] = observations[0].get("input_bytes")
        summary["input_sha256"] = observations[0].get("input_sha256")
        summary["input_bases"] = observations[0].get("input_bases")
    ok = [
        entry for entry in observations
        if entry.get("status") == "success" and entry.get("roundtrip_ok")
    ]
    if not ok:
        summary["error"] = next(
            (entry.get("error") for entry in observations if entry.get("error")),
            "no successful observation",
        )
        summary["attempts"] = len(observations)
        summary["failed_observations"] = sum(
            1 for entry in observations
            if entry.get("status") != "success" or not entry.get("roundtrip_ok")
        )
        attempt_walls = [entry.get("wall_seconds") for entry in observations]
        usable_walls = [value for value in attempt_walls if value is not None]
        summary["attempt_wall_seconds_median"] = _median(attempt_walls)
        summary["attempt_wall_seconds_max"] = (
            max(usable_walls) if usable_walls else None
        )
        attempt_rss = [entry.get("peak_rss_bytes") for entry in observations]
        usable_rss = [value for value in attempt_rss if value is not None]
        summary["attempt_peak_rss_bytes"] = (
            max(usable_rss) if usable_rss else None
        )
        cgroup_peaks = [
            entry["cgroup_memory_peak_bytes"] for entry in observations
            if entry.get("cgroup_memory_peak_bytes") is not None
        ]
        oom_kills = [
            entry["cgroup_oom_kills"] for entry in observations
            if entry.get("cgroup_oom_kills") is not None
        ]
        summary["cgroup_memory_peak_bytes"] = (
            max(cgroup_peaks) if cgroup_peaks else None
        )
        summary["cgroup_oom_kills"] = max(oom_kills) if oom_kills else None
        summary["oom_killed_observations"] = sum(
            1 for entry in observations
            if (entry.get("cgroup_oom_kills") or 0) > 0
        )
        return summary

    sizes = {
        entry["artifact_bytes"] for entry in ok
        if entry.get("artifact_bytes") is not None
    }
    shas = {
        (entry.get("verification") or {}).get("archive_sha256") for entry in ok
    }
    shas.discard(None)
    summary["deterministic"] = len(sizes) == 1 and len(shas) <= 1
    summary["roundtrip_ok"] = len(ok) == len(observations)
    summary["status"] = (
        "success"
        if summary["roundtrip_ok"] and summary["deterministic"]
        else "partial"
    )
    summary["error"] = next(
        (entry.get("error") for entry in observations if entry.get("error")),
        None,
    )
    summary["artifact_bytes"] = next(iter(sizes)) if len(sizes) == 1 else None
    summary["artifact_sha256"] = next(iter(shas)) if len(shas) == 1 else None
    if summary["artifact_bytes"] is not None and summary.get("input_bases"):
        summary["bits_per_base"] = (
            summary["artifact_bytes"] * 8 / summary["input_bases"]
        )

    walls = [entry.get("wall_seconds") for entry in ok]
    cpus = [entry.get("cpu_seconds") for entry in ok]
    rss = [entry.get("peak_rss_bytes") for entry in ok]
    summary["wall_seconds"] = _median(walls)
    usable_walls = [value for value in walls if value is not None]
    summary["wall_seconds_min"] = min(usable_walls) if usable_walls else None
    summary["wall_seconds_max"] = max(usable_walls) if usable_walls else None
    summary["cpu_seconds"] = _median(cpus)
    usable_rss = [value for value in rss if value is not None]
    summary["peak_rss_bytes"] = max(usable_rss) if usable_rss else None
    summary["peak_rss_median_bytes"] = _median(usable_rss)
    cgroup_peaks = [
        entry["cgroup_memory_peak_bytes"] for entry in ok
        if entry.get("cgroup_memory_peak_bytes") is not None
    ]
    oom_kills = [
        entry["cgroup_oom_kills"] for entry in ok
        if entry.get("cgroup_oom_kills") is not None
    ]
    summary["cgroup_memory_peak_bytes"] = (
        max(cgroup_peaks) if cgroup_peaks else None
    )
    summary["cgroup_oom_kills"] = max(oom_kills) if oom_kills else None

    restores = [
        entry["decompress"] for entry in ok if entry.get("decompress")
    ]
    if restores:
        restore_rss = [
            value.get("peak_rss_bytes") for value in restores
            if value.get("peak_rss_bytes") is not None
        ]
        summary["decompress"] = {
            "wall_seconds": _median([
                value.get("wall_seconds") for value in restores
            ]),
            "cpu_seconds": _median([
                value.get("cpu_seconds") for value in restores
            ]),
            "peak_rss_bytes": max(restore_rss) if restore_rss else None,
            "cgroup_memory_peak_bytes": max(
                [
                    value["cgroup_memory_peak_bytes"] for value in restores
                    if value.get("cgroup_memory_peak_bytes") is not None
                ],
                default=None,
            ),
        }
    return summary


def evaluate_gate(
    summaries: list[dict],
    baseline_id: str,
    rss_mib: float,
    size_ratio: float,
    wall_ratio: float,
    candidates: list[str] | None = None,
    size_baseline_bytes: int | None = None,
    wall_baseline_seconds: float | None = None,
) -> dict:
    """Evaluate the published bounded-mode gates against one baseline profile.

    RSS is the worst observed peak (compress and restore independently). Size
    and wall use the aggregated archive and median wall. Any failed repeat makes
    the profile ineligible, so a lucky sample cannot pass the gate. ``candidates``
    restricts the gate to the bounded profiles under test; otherwise every other
    selected profile is treated as a candidate.

    ``size_baseline_bytes`` overrides the archive size used for the size ratio.
    ``wall_baseline_seconds`` pins the wall baseline to a previously published
    median, so a faster or slower session cannot silently change the gate. The
    default baseline for both is the comparison profile measured in this run.
    """
    by_id = {summary["profile"]: summary for summary in summaries}
    baseline = by_id.get(baseline_id)
    gate = {
        "baseline": baseline_id,
        "candidates": candidates,
        "thresholds": {
            "rss_mib": rss_mib,
            "size_ratio": size_ratio,
            "wall_ratio": wall_ratio,
            "size_baseline_bytes": size_baseline_bytes,
            "wall_baseline_seconds": wall_baseline_seconds,
        },
        "checks": [],
        "passed": False,
        "error": None,
    }
    if baseline is None:
        gate["error"] = f"baseline profile not selected: {baseline_id}"
        return gate
    if baseline.get("status") != "success":
        gate["error"] = f"baseline profile {baseline_id} did not succeed"
        return gate
    if candidates is None:
        selected = [
            summary for summary in summaries
            if summary["profile"] != baseline_id
        ]
    else:
        unknown = sorted(set(candidates) - set(by_id))
        if unknown:
            gate["error"] = f"gate profiles not selected: {', '.join(unknown)}"
            return gate
        if not candidates:
            gate["error"] = "no gate candidate profiles"
            return gate
        wanted = set(candidates)
        selected = [
            summary for summary in summaries
            if summary["profile"] in wanted
        ]

    for summary in selected:
        profile = summary["profile"]
        gate["checks"].append({
            "profile": profile,
            "check": "status",
            "passed": summary.get("status") == "success",
            "detail": summary.get("status"),
        })
        compress_rss = summary.get("peak_rss_bytes")
        restore_rss = (summary.get("decompress") or {}).get("peak_rss_bytes")
        for label, value in (
            ("compress_rss", compress_rss),
            ("restore_rss", restore_rss),
        ):
            gate["checks"].append({
                "profile": profile,
                "check": label,
                "value_mib": None if value is None else value / 1048576,
                "limit_mib": rss_mib,
                "passed": value is not None and value / 1048576 <= rss_mib,
            })
        size = summary.get("artifact_bytes")
        base_size = (
            size_baseline_bytes if size_baseline_bytes is not None
            else baseline.get("artifact_bytes")
        )
        size_value = (size / base_size) if size and base_size else None
        gate["checks"].append({
            "profile": profile,
            "check": "size_ratio",
            "value": size_value,
            "limit": size_ratio,
            "passed": size_value is not None and size_value <= size_ratio,
        })
        wall = summary.get("wall_seconds")
        base_wall = (
            wall_baseline_seconds if wall_baseline_seconds is not None
            else baseline.get("wall_seconds")
        )
        wall_value = (wall / base_wall) if wall and base_wall else None
        gate["checks"].append({
            "profile": profile,
            "check": "wall_ratio",
            "value": wall_value,
            "limit": wall_ratio,
            "passed": wall_value is not None and wall_value <= wall_ratio,
        })
    gate["passed"] = bool(gate["checks"]) and all(
        check["passed"] for check in gate["checks"]
    )
    return gate


def print_gate(gate: dict) -> None:
    if gate.get("error"):
        print(f"  gate: FAIL ({gate['error']})", flush=True)
        return
    verdict = "PASS" if gate["passed"] else "FAIL"
    print(f"  gate vs {gate['baseline']}: {verdict}", flush=True)
    for check in gate["checks"]:
        if check["check"] == "status":
            value = check.get("detail")
            limit = "success"
        elif check["check"].endswith("_rss"):
            value = f"{check['value_mib']:.1f} MiB" if check["value_mib"] is not None else "n/a"
            limit = f"{check['limit_mib']:.1f} MiB"
        else:
            value = f"{check['value']:.3f}" if check["value"] is not None else "n/a"
            limit = f"{check['limit']:.3f}"
        print(
            f"    {check['profile']:<38} {check['check']:<14} "
            f"{value:>12} / {limit:<10} {'ok' if check['passed'] else 'FAIL'}",
            flush=True,
        )


def apply_overrides(
    profiles: list[dict],
    heap: str | None = None,
    memory: str | None = None,
    jvm_extra: str | None = None,
    runtime: str | None = None,
    binary: str | None = None,
) -> list[dict]:
    """Copy profiles and apply sweep overrides without mutating PROFILES."""
    selected = [dict(profile) for profile in profiles]
    extra = jvm_extra.split() if jvm_extra else []
    for profile in selected:
        if heap and profile["kind"] == "hypergrass":
            profile["heap"] = heap
        if memory and profile["kind"] == "hypergrass":
            profile["memory"] = memory
        if runtime and profile["kind"] == "hypergrass":
            profile["runtime"] = runtime
        if binary and profile["kind"] == "hypergrass":
            profile["binary"] = binary
        if extra and profile["kind"] == "hypergrass" \
                and profile.get("runtime") != "native":
            profile["jvm"] = [*profile.get("jvm", []), *extra]
    return selected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--source", default=str(SOURCE))
    parser.add_argument("--slice-bases", type=int, default=0)
    parser.add_argument("--profiles", default=",".join(p["id"] for p in PROFILES))
    parser.add_argument("--timeout", type=float, default=7200.0)
    parser.add_argument("--memory", default="2g")
    parser.add_argument(
        "--heap-override", default=None,
        help="override -Xmx for selected hypergrass profiles (sweep probe)",
    )
    parser.add_argument(
        "--memory-override", default=None,
        help="override the container memory cap for selected hypergrass profiles",
    )
    parser.add_argument(
        "--jvm-extra", default=None,
        help="extra JVM flags for hypergrass profiles, e.g. NativeMemoryTracking",
    )
    parser.add_argument(
        "--runtime-override", choices=("jvm", "native"), default=None,
        help="run hypergrass profiles on the JVM or a native-image binary",
    )
    parser.add_argument(
        "--binary-override", default=None,
        help="project-relative native binary for --runtime-override native",
    )
    parser.add_argument("--runs-dir", default=str(RUNS))
    parser.add_argument(
        "--repeats", type=int, default=3,
        help="interleaved observations per profile (default 3)",
    )
    parser.add_argument(
        "--gate-baseline", default=None,
        help="profile id to gate against, e.g. jarvis3-l7",
    )
    parser.add_argument(
        "--gate-profiles", default=None,
        help="comma-separated candidate profiles; default: every bounded-* profile",
    )
    parser.add_argument("--gate-rss-mib", type=float, default=256.0)
    parser.add_argument("--gate-size-ratio", type=float, default=1.02)
    parser.add_argument("--gate-size-baseline-bytes", type=int, default=None)
    parser.add_argument("--gate-wall-baseline-seconds", type=float, default=None)
    parser.add_argument("--gate-wall-ratio", type=float, default=2.0)
    args = parser.parse_args(argv)
    if args.repeats < 1:
        raise SystemExit("--repeats must be at least 1")

    source = Path(args.source).resolve()
    if not source.exists():
        raise SystemExit(f"missing source FASTA: {source}")

    image = pangenome_run.image_id(PROJECT)
    wanted = set(args.profiles.split(","))
    base_profiles = apply_overrides(
        [p for p in PROFILES if p["id"] in wanted],
        heap=args.heap_override,
        memory=args.memory_override,
        jvm_extra=args.jvm_extra,
        runtime=args.runtime_override,
        binary=args.binary_override,
    )
    missing = wanted - {p["id"] for p in base_profiles}
    if missing:
        raise SystemExit(f"unknown profiles: {', '.join(sorted(missing))}")
    if not base_profiles:
        raise SystemExit("no profiles selected")
    if args.gate_baseline and args.gate_baseline not in {
        p["id"] for p in base_profiles
    }:
        raise SystemExit(
            f"gate baseline {args.gate_baseline} is not in the selected profiles"
        )

    # Expand repeats before the execution loop so profiles are interleaved:
    # r0/A, r0/B, r1/A, r1/B ... This spreads host drift across all profiles
    # instead of letting one profile absorb a slow or noisy period.
    profiles = []
    for repeat in range(args.repeats):
        for base in base_profiles:
            expanded = dict(base)
            expanded["base_profile"] = base["id"]
            expanded["repeat"] = repeat
            expanded["id"] = f"{base['id']}-r{repeat}"
            profiles.append(expanded)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = args.run_id or f"specialized-{stamp}"
    runs_root = Path(args.runs_dir).resolve()
    runs_root.mkdir(parents=True, exist_ok=True)
    run_dir = runs_root / run_id
    if run_dir.exists():
        raise SystemExit(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)

    if args.slice_bases:
        fasta_path = run_dir / "input.fasta"
        records, bases, size = slice_fasta(source, args.slice_bases, fasta_path)
    else:
        fasta_path = source
        records, bases = parse_input(source)
        size = source.stat().st_size

    fasta_digest = pangenome_run.sha256_file(fasta_path)
    sequence = sequence_only(fasta_path, run_dir / "input.seq")
    report = {
        "schema_version": 2,
        "kind": "specialized-scale-test",
        "run_id": run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "image_id": image,
        "image_tag": pangenome_run.IMAGE_TAG,
        "provenance": pangenome_run.provenance(PROJECT),
        "harness": {
            "repeats": args.repeats,
            "timeout_seconds": args.timeout,
            "memory_default": args.memory,
            "heap_override": args.heap_override,
            "memory_override": args.memory_override,
            "jvm_extra": args.jvm_extra,
            "runtime_override": args.runtime_override,
            "binary_override": args.binary_override,
        },
        "source": str(source.relative_to(PROJECT)),
        "sliced": bool(args.slice_bases),
        "input": {
            "fasta": {
                "path": str(fasta_path.relative_to(PROJECT)),
                "size_bytes": size,
                "sha256": fasta_digest,
                "records": records,
                "bases": bases,
            },
            "sequence": sequence,
            "normalization": (
                "headers removed, letters uppercased, non-ACGT replaced with A; "
                "both tools receive this sequence and are verified against it"
            ),
        },
        "profiles": [p["id"] for p in base_profiles],
        "execution_order": [p["id"] for p in profiles],
        "results": [],
        "status": "running",
    }
    native_binaries = {}
    for profile in base_profiles:
        if profile.get("runtime") != "native":
            continue
        relative = profile.get("binary")
        if not relative:
            raise SystemExit(
                f"{profile['id']}: --runtime-override native needs "
                "--binary-override"
            )
        binary = (PROJECT / relative).resolve()
        if not binary.is_file():
            raise SystemExit(f"missing native binary: {binary}")
        native_binaries[relative] = {
            "bytes": binary.stat().st_size,
            "sha256": pangenome_run.sha256_file(binary),
        }
    if native_binaries:
        report["native_binaries"] = native_binaries
    if "grass" in wanted:
        report["vendor"] = {
            "grass": {
                "repo": "https://github.com/AnshuOishik/GraSS",
                "commit": GRASS_COMMIT,
                "jvm": "OpenJDK 27",
            },
            "libbsc": {
                "repo": "https://github.com/IlyaGrebnov/libbsc",
                "commit": LIBBSC_COMMIT,
                "bsc_version": "3.4.0",
                "note": "GraSS invokes ./bsc e FinalTar.tar FinalBsc.bsc -e2",
            },
        }
    if not args.slice_bases and source.resolve() == SOURCE.resolve():
        baseline_path = PROJECT / "benchmarks/pangenome/runs/pangenome-scale-2/report.json"
        if baseline_path.exists():
            baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
            for entry in baseline.get("results", []):
                if entry.get("profile") == "xz-9" and entry.get("status") == "success":
                    report["baseline"] = {
                        "profile": "xz-9",
                        "archive_bytes": entry["artifact_bytes"],
                        "source": str(baseline_path.relative_to(PROJECT)),
                    }
    pangenome_run.write_json(run_dir / "report.json", report)
    print(
        f"run {run_id}: {len(profiles)} profiles on {sequence['size_bytes']:,} bytes, "
        f"{sequence['bases']:,} bases, {sequence['non_acgt_bases_masked']} masked"
    )

    failures = 0
    for profile in profiles:
        print(f"  {profile['id']} ...", flush=True)
        entry = {
            "profile": profile["id"],
            "base_profile": profile["base_profile"],
            "repeat": profile["repeat"],
            "kind": profile["kind"],
            "level": profile["level"],
            "input_bytes": sequence["size_bytes"],
            "input_sha256": sequence["sha256"],
            "input_bases": sequence["bases"],
            "status": "error",
            "roundtrip_ok": False,
            "artifact_bytes": None,
            "bits_per_base": None,
            "wall_seconds": None,
            "cpu_seconds": None,
            "peak_rss_bytes": None,
            "error": None,
        }
        artifact: Path | None = None
        restored: Path | None = None
        memory = profile.get("memory", args.memory)
        try:
            if profile["kind"] == "jarvis3":
                op_dir = run_dir / "ops" / profile["id"]
                artifact = op_dir / "archive.jc"
                restored = op_dir / "restored.seq"
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    ["JARVIS3", "-l", profile["level"], "-o",
                     pangenome_run.to_container(PROJECT, artifact),
                     pangenome_run.to_container(PROJECT, run_dir / "input.seq")],
                    run_dir, args.timeout, memory,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        ["JARVIS3", "-d", "-o",
                         pangenome_run.to_container(PROJECT, restored),
                         pangenome_run.to_container(PROJECT, artifact)],
                        run_dir, args.timeout, memory,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            elif profile["kind"] == "geco3":
                op_dir = run_dir / "ops" / profile["id"]
                op_dir.mkdir(parents=True, exist_ok=True)
                local_input = op_dir / "input.seq"
                shutil.copyfile(run_dir / "input.seq", local_input)
                artifact = op_dir / "input.seq.co"
                restored = op_dir / "input.seq.de"
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    ["GeCo3", "-l", profile["level"],
                     pangenome_run.to_container(PROJECT, local_input)],
                    op_dir, args.timeout, args.memory,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        ["GeDe3", pangenome_run.to_container(PROJECT, artifact)],
                        op_dir, args.timeout, memory,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            elif profile["kind"] == "xz":
                op_dir = run_dir / "ops" / profile["id"]
                artifact = op_dir / "archive.xz"
                restored = op_dir / "restored.seq"
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    ["xz", "-9", "-c",
                     pangenome_run.to_container(PROJECT, run_dir / "input.seq")],
                    run_dir, args.timeout, memory,
                    stdout_artifact=artifact,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        ["xz", "-dc", pangenome_run.to_container(PROJECT, artifact)],
                        run_dir, args.timeout, memory,
                        stdout_artifact=restored,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            elif profile["kind"] == "hypergrass":
                op_dir = run_dir / "ops" / profile["id"]
                artifact = op_dir / "archive.hgs"
                restored = op_dir / "restored.seq"
                if profile.get("runtime") == "native":
                    relative = profile.get("binary")
                    if not relative:
                        raise SystemExit(
                            f"{profile['id']}: --runtime-override native needs "
                            "--binary-override"
                        )
                    binary = (PROJECT / relative).resolve()
                    if not binary.is_file():
                        raise SystemExit(f"missing native binary: {binary}")
                    command = [
                        pangenome_run.to_container(PROJECT, binary),
                        f"-Xmx{profile.get('heap', '4g')}",
                        *[
                            flag for flag in profile.get("jvm", [])
                            if flag.startswith("-D")
                        ],
                    ]
                else:
                    command = [
                        "java", f"-Xmx{profile.get('heap', '4g')}",
                        "-XX:ActiveProcessorCount=1", *profile.get("jvm", []),
                        "-cp", "/work/bin", "hypergrass.HyperGraSS",
                    ]
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    command + ["compress", *profile["options"],
                               pangenome_run.to_container(PROJECT, run_dir / "input.seq"),
                               pangenome_run.to_container(PROJECT, artifact)],
                    run_dir, args.timeout, memory,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        command + ["decompress",
                                   pangenome_run.to_container(PROJECT, artifact),
                                   pangenome_run.to_container(PROJECT, restored)],
                        run_dir, args.timeout, memory,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            elif profile["kind"] == "grass":
                # GraSS is a FASTA/text codec, but its raw path accepts a headerless
                # single-line sequence and restores it byte-for-byte (the "0" first
                # line in F1.SSG). That lets it consume input.seq unchanged.
                op_dir = run_dir / "ops" / profile["id"] / "compress"
                dec_dir = run_dir / "ops" / profile["id"] / "decompress"
                op_dir.mkdir(parents=True, exist_ok=True)
                dec_dir.mkdir(parents=True, exist_ok=True)
                bsc = PROJECT / "benchmarks/specialized/vendor/bin/bsc"
                if not bsc.is_file():
                    raise FileNotFoundError(
                        "missing vendor/bin/bsc; run vendor/fetch_vendor.ps1"
                    )
                shutil.copyfile(bsc, op_dir / "bsc")
                artifact = op_dir / "FinalBsc.bsc"
                restored = dec_dir / "Final.SSG"
                java = [
                    "java", f"-Xmx{profile.get('heap', '1500m')}",
                    "-XX:ActiveProcessorCount=1",
                    "-cp", "/work/benchmarks/specialized/vendor/GraSS/build",
                    "SSG.SSG",
                ]
                compress = pangenome_run.run_operation(
                    PROJECT, image, run_dir, f"{profile['id']}/compress",
                    java + ["compress",
                            pangenome_run.to_container(PROJECT, run_dir / "input.seq")],
                    op_dir, args.timeout, memory,
                )
                entry.update(compress)
                if compress["status"] == "success" and artifact.exists():
                    shutil.copyfile(artifact, dec_dir / "FinalBsc.bsc")
                    shutil.copyfile(bsc, dec_dir / "bsc")
                    decompress = pangenome_run.run_operation(
                        PROJECT, image, run_dir, f"{profile['id']}/decompress",
                        java + ["decompress"],
                        dec_dir, args.timeout, memory,
                    )
                    entry["decompress"] = decompress
                    entry["roundtrip_ok"] = (
                        decompress["status"] == "success"
                        and restored.exists()
                        and sequence_sha256(restored) == sequence["sha256"]
                    )
            else:
                raise ValueError(f"unknown kind: {profile['kind']}")
        except Exception as error:
            entry["error"] = f"{type(error).__name__}: {error}"

        if entry["status"] == "success" and entry.get("artifact_bytes") is None:
            if profile["kind"] == "jarvis3":
                artifact = run_dir / "ops" / profile["id"] / "archive.jc"
            elif profile["kind"] == "geco3":
                artifact = run_dir / "ops" / profile["id"] / "input.seq.co"
            elif profile["kind"] == "hypergrass":
                artifact = run_dir / "ops" / profile["id"] / "archive.hgs"
            elif profile["kind"] == "grass":
                artifact = run_dir / "ops" / profile["id"] / "compress" / "FinalBsc.bsc"
            else:
                artifact = run_dir / "ops" / profile["id"] / "archive.xz"
            if artifact.exists():
                entry["artifact_bytes"] = artifact.stat().st_size
        if entry.get("artifact_bytes"):
            entry["bits_per_base"] = entry["artifact_bytes"] * 8 / sequence["bases"]
        if artifact is not None and restored is not None:
            if artifact.is_file() and restored.is_file() and restored.stat().st_size > 0:
                verification = independent_verification(
                    run_dir / "input.seq", restored, artifact, sequence["sha256"],
                    run_dir / "ops" / profile["id"] / "verification.json",
                )
                entry["verification"] = verification
                entry["restored_sha256"] = verification.get("restored_sha256")
                entry["roundtrip_ok"] = bool(entry["roundtrip_ok"]) and (
                    verification.get("verdict") == "pass"
                )
        if (
            entry["status"] == "success"
            and not entry["roundtrip_ok"]
            and not entry.get("error")
        ):
            entry["error"] = roundtrip_failure(entry)
        if entry["status"] != "success" or not entry["roundtrip_ok"]:
            failures += 1
            print(f"    FAILED: status={entry['status']} roundtrip={entry['roundtrip_ok']} "
                  f"error={entry.get('error')}", flush=True)
        else:
            print(
                f"    {entry['artifact_bytes']:,} B, "
                f"{entry['bits_per_base']:.3f} bits/base, "
                f"{entry['wall_seconds']:.1f}s, "
                f"{entry['peak_rss_bytes'] / 1048576:.0f} MiB RSS",
                flush=True,
            )
        report["results"].append(entry)
        pangenome_run.write_json(run_dir / "report.json", report)

    raw_observations = report["results"]
    report["results"] = [
        summarize_profile(
            base,
            [
                entry for entry in raw_observations
                if entry.get("base_profile") == base["id"]
            ],
            args.repeats,
        )
        for base in base_profiles
    ]
    failures = sum(
        1 for summary in report["results"]
        if summary["status"] != "success"
    )
    print("  summary:", flush=True)
    for summary in report["results"]:
        failed = summary["wall_seconds"] is None
        wall = (
            summary.get("attempt_wall_seconds_median")
            if failed else summary["wall_seconds"]
        )
        rss = (
            summary.get("attempt_peak_rss_bytes")
            if failed else summary["peak_rss_bytes"]
        )
        label = "attempt" if failed else "median"
        oom = (
            f"  oom {summary.get('oom_killed_observations', 0)}/"
            f"{summary.get('attempts', 0)}"
            if failed else ""
        )
        print(
            f"    {summary['profile']:<38} {summary['status']:<8} "
            f"{summary['artifact_bytes'] if summary['artifact_bytes'] else '-':>10} B  "
            f"{label} {wall if wall is not None else 0:6.1f}s  "
            f"worst {rss / 1048576 if rss is not None else 0:5.1f} MiB  "
            f"cgroup {summary['cgroup_memory_peak_bytes'] / 1048576 if summary['cgroup_memory_peak_bytes'] else 0:5.1f} MiB"
            f"{oom}",
            flush=True,
        )

    if args.gate_baseline:
        if args.gate_profiles:
            gate_candidates = [
                item.strip()
                for item in args.gate_profiles.split(",")
                if item.strip()
            ]
        else:
            gate_candidates = [
                summary["profile"] for summary in report["results"]
                if str(summary.get("level", "")).startswith("bounded")
            ]
        report["gate"] = evaluate_gate(
            report["results"],
            args.gate_baseline,
            args.gate_rss_mib,
            args.gate_size_ratio,
            args.gate_wall_ratio,
            candidates=gate_candidates,
            size_baseline_bytes=args.gate_size_baseline_bytes,
            wall_baseline_seconds=args.gate_wall_baseline_seconds,
        )
        print_gate(report["gate"])

    report["status"] = "success" if failures == 0 else "partial"
    report["failures"] = failures
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    pangenome_run.write_json(run_dir / "report.json", report)
    print(f"report: {run_dir / 'report.json'}")
    if failures:
        return 2
    if args.gate_baseline and not report["gate"]["passed"]:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
