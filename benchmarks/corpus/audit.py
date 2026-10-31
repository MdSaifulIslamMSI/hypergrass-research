"""Audit a finished corpus report against the files it claims to describe.

This is the automated first review gate. It does not trust the summary or the
recorded hashes: it re-hashes every artifact and restored file on disk, checks
the arithmetic, verifies the expected result matrix, and recomputes each
summary row. A second, human-level gate is the lossless verification suite plus
review of the deviations and conclusions.

Usage:
    python -B benchmarks/corpus/audit.py benchmarks/corpus/runs/<run-id>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

READ_CHUNK = 1024 * 1024
ARTIFACTS = {
    "hypergrass": "archive.hgs",
    "hypergrass-optimize": "archive.hgs",
    "hypergrass-lean": "archive.hgs",
    "jarvis3-l7": "archive.jc",
    "jarvis3-l26": "archive.jc",
    "gzip-9": "archive.gz",
    "xz-9": "archive.xz",
    "geco3-l9": "input.seq.co",
    "geco3-l6": "input.seq.co",
    "geco3-l3": "input.seq.co",
}
RESTORED = {
    "hypergrass": "restored.seq",
    "hypergrass-optimize": "restored.seq",
    "hypergrass-lean": "restored.seq",
    "jarvis3-l7": "restored.seq",
    "jarvis3-l26": "restored.seq",
    "gzip-9": "restored.seq",
    "xz-9": "restored.seq",
    "geco3-l9": "input.seq.de",
    "geco3-l6": "input.seq.de",
    "geco3-l3": "input.seq.de",
}
EXPECTED_REPS = {profile: 3 for profile in ARTIFACTS}
EXPECTED_REPS["geco3-l9"] = 1
ACGT_ONLY = {"jarvis3-l7", "jarvis3-l26", "geco3-l9", "geco3-l6", "geco3-l3"}
# Profiles whose configured effort can exceed the 180 s cap or the 2 GB
# container limit on MB-scale inputs. Their incomplete cells are retained as
# infeasible; no size is inferred from a partial run.
HIGH_EFFORT = {"geco3-l9", "jarvis3-l26"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(READ_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


class Audit:
    def __init__(self):
        self.errors: list[str] = []
        self.checks = 0

    def check(self, condition: bool, message: str) -> None:
        self.checks += 1
        if not condition:
            self.errors.append(message)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir")
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="audit the results that exist in a run that was stopped early; "
             "never treats the partial matrix as complete",
    )
    args = parser.parse_args(argv)
    run_dir = Path(args.run_dir).resolve()
    corpus = run_dir.parent.parent
    report = json.loads((run_dir / "report.json").read_text(encoding="utf-8"))
    protocol = json.loads((corpus / "protocol.json").read_text(encoding="utf-8"))
    data_lock = json.loads((corpus / "data.lock.json").read_text(encoding="utf-8"))
    datasets = {item["id"]: item for item in data_lock["datasets"]}
    audit = Audit()

    audit.check(
        report.get("status") in ("success", "partial")
        or (args.allow_incomplete and report.get("status") in ("running", "stopped")),
        f"report status is {report.get('status')}",
    )
    expected_infeasible = sum(
        1 for entry in report["results"]
        if entry["profile"] in HIGH_EFFORT and entry["status"] != "success"
        and entry["status"] != "not_applicable"
    )
    if not args.allow_incomplete:
        audit.check(
            report.get("failures") == expected_infeasible,
            f"report failures {report.get('failures')} != expected high-effort infeasible {expected_infeasible}",
        )
        audit.check(
            report.get("status") == ("success" if expected_infeasible == 0 else "partial"),
            "report status does not match the expected infeasible count",
        )
    audit.check(len(report.get("deviations", [])) >= 3, "deviations were not recorded")
    audit.check(
        any("geco3-l9" in text and "once" in text for text in report.get("deviations", [])),
        "the geco3-l9 repetition deviation is missing",
    )

    expected = {}
    for dataset in datasets.values():
        for profile in ARTIFACTS:
            expected[(dataset["id"], profile)] = EXPECTED_REPS[profile]

    results = {}
    for entry in report["results"]:
        key = (entry["dataset"], entry["profile"])
        results.setdefault(key, []).append(entry)
        dataset = datasets[entry["dataset"]]
        if entry["status"] == "not_applicable":
            audit.check(entry["profile"] in ACGT_ONLY, f"{key} unexpectedly not applicable")
            audit.check(not dataset["sequence"]["raw_acgt"], f"{key} flagged n/a on ACGT input")
            continue
        audit.check(entry["input_bytes"] == dataset["sequence"]["size_bytes"], f"{key} input size mismatch")
        audit.check(entry["input_sha256"] == dataset["sequence"]["sha256"], f"{key} input hash mismatch")
        if entry["status"] != "success":
            compress = entry.get("compress", {})
            decompress = entry.get("decompress", {})
            allowed = entry["profile"] in HIGH_EFFORT and (
                compress.get("status") == "timeout"
                or (
                    compress.get("status") == "nonzero"
                    and compress.get("returncode") in (137, -9)
                )
                or decompress.get("status") == "timeout"
                or (
                    decompress.get("status") == "nonzero"
                    and decompress.get("returncode") in (137, -9)
                )
            )
            audit.check(allowed, f"{key} unexpected failure: {entry.get('error')}")
            continue
        if entry["profile"] in ACGT_ONLY:
            audit.check(bool(dataset["sequence"]["raw_acgt"]), f"{key} DNA tool ran on non-ACGT input")
        audit.check(entry.get("decompress", {}).get("status") == "success", f"{key} decompression status")

        work = run_dir / "work" / entry["dataset"] / entry["profile"] / f"rep-{entry['repetition']}"
        artifact = work / ARTIFACTS[entry["profile"]]
        restored = work / RESTORED[entry["profile"]]
        audit.check(artifact.exists(), f"{key} artifact missing: {artifact}")
        audit.check(restored.exists(), f"{key} restored file missing: {restored}")
        if artifact.exists():
            audit.check(artifact.stat().st_size == entry["artifact_bytes"], f"{key} artifact size mismatch")
            audit.check(sha256_file(artifact) == entry["artifact_sha256"], f"{key} artifact hash mismatch")
        if restored.exists():
            audit.check(restored.stat().st_size == entry["restored_bytes"], f"{key} restored size mismatch")
            audit.check(sha256_file(restored) == entry["restored_sha256"], f"{key} restored hash mismatch")
            audit.check(sha256_file(restored) == entry["input_sha256"], f"{key} restored bytes differ from input")
        audit.check(entry.get("roundtrip_ok") is True, f"{key} roundtrip_ok is not true")
        if entry["artifact_bytes"]:
            ratio = entry["input_bytes"] / entry["artifact_bytes"]
            bpb = entry["artifact_bytes"] * 8.0 / entry["input_bytes"]
            audit.check(abs(ratio - entry["ratio"]) < 1e-9, f"{key} ratio arithmetic")
            audit.check(abs(bpb - entry["bits_per_base"]) < 1e-9, f"{key} bpb arithmetic")
        for query in entry.get("queries", []):
            audit.check(query["status"] == "success", f"{key} query status")
            audit.check(query["exact_after_trailing_newline"], f"{key} query mismatch")

    for key, reps in results.items():
        if args.allow_incomplete:
            audit.check(
                len(reps) <= expected.get(key),
                f"{key} has {len(reps)} results, more than the {expected.get(key)} planned",
            )
        else:
            audit.check(len(reps) == expected.get(key), f"{key} result count {len(reps)} != {expected.get(key)}")

    # Recompute summary rows from raw results.
    for row in report["summary"]:
        key = (row["dataset"], row["profile"])
        group = results.get(key, [])
        if group and all(entry["status"] == "not_applicable" for entry in group):
            audit.check(row["attempts"] == len(group), f"{key} n/a attempts")
            audit.check(row["successful"] == 0, f"{key} n/a successful")
            audit.check(row["artifact_bytes"] is None, f"{key} n/a artifact bytes")
            audit.check(row["status"] == "not_applicable", f"{key} n/a status")
            continue
        applicable = [entry for entry in group if entry["status"] != "not_applicable"]
        successful = [entry for entry in applicable if entry["status"] == "success"]
        sizes = sorted({entry["artifact_bytes"] for entry in successful if entry.get("artifact_bytes")})
        audit.check(row["attempts"] == len(applicable), f"{key} summary attempts")
        audit.check(row["successful"] == len(successful), f"{key} summary successful")
        audit.check(
            row["artifact_bytes"] == (sizes[0] if len(sizes) == 1 else None),
            f"{key} summary artifact bytes",
        )
        audit.check(not row["nondeterministic_size"], f"{key} artifact size changed between repetitions")
        audit.check(
            row["roundtrip_all_ok"] == (bool(successful) and all(entry["roundtrip_ok"] for entry in successful)),
            f"{key} summary roundtrip flag",
        )
        if row["status"] == "success":
            audit.check(len(successful) == len(applicable), f"{key} success status with failures")
        elif row["status"] == "timeout":
            audit.check(
                any(
                    entry.get("compress", {}).get("status") == "timeout"
                    or entry.get("decompress", {}).get("status") == "timeout"
                    for entry in applicable
                ),
                f"{key} timeout status",
            )
        elif row["status"] == "failed":
            audit.check(
                row["profile"] in HIGH_EFFORT and any(entry["status"] != "success" for entry in applicable),
                f"{key} failed status",
            )

    audit.check(
        report["environment"]["image_id"].startswith("sha256:"),
        "image id is not immutable",
    )
    audit.check(
        report["environment"]["data_protocol_sha256"] == data_lock.get("protocol_sha256"),
        "protocol hash mismatch",
    )
    audit.check(report["protocol"] == protocol, "embedded protocol differs from protocol.json")

    planned = sum(expected.values())
    present = len(report["results"])
    complete_datasets = sum(
        1 for dataset in datasets.values()
        if all(len(results.get((dataset["id"], profile), [])) == reps
               for profile, reps in EXPECTED_REPS.items())
    )
    coverage = (
        f"coverage {present}/{planned} results, "
        f"{complete_datasets}/{len(datasets)} datasets complete"
    )
    if audit.errors:
        print(f"AUDIT FAILED: {len(audit.errors)} of {audit.checks} checks; {coverage}")
        for error in audit.errors[:40]:
            print(f"  - {error}")
        if len(audit.errors) > 40:
            print(f"  ... {len(audit.errors) - 40} more")
        return 1
    print(f"AUDIT PASSED: {audit.checks} checks, {len(results)} dataset/profile cells; {coverage}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
