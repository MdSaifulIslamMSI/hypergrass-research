"""Run the frozen held-out corpus against HyperGraSS and pinned baselines.

Every measured operation executes in one fresh process inside the locked Linux
container. Docker startup is excluded from the reported wall time: the Python
wrapper inside the container starts its clock immediately before spawning the
tool and reads child rusage after reaping. Sizes, SHA-256 hashes, and byte-exact
round trips are verified for every completed operation.

Usage (from the HyperGraSS project root):

    javac --release 21 -d bin src/hypergrass/*.java
    python -B benchmarks/corpus/run.py

The frozen protocol names one GeCo3 profile, ``geco3-l9``. A 1 MiB pilot
exceeded the protocol's 180 s operation cap, so that profile is measured once
per dataset. ``geco3-l3`` and ``geco3-l6`` are supplementary feasible profiles
and are labelled as such in the report. No engine parameters are tuned after
seeing results.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import sys
import time

from process_metrics import run_process


SCHEMA_VERSION = 1
IMAGE_TAG = "hypergrass-corpus:locked"
READ_CHUNK = 1024 * 1024

PROFILES = (
    {
        "id": "hypergrass",
        "kind": "hgs",
        "heap": "1g",
        "options": [],
        "requires_acgt": False,
        "repetitions": None,
        "origin": "frozen",
    },
    {
        "id": "hypergrass-optimize",
        "kind": "hgs",
        "heap": "1g",
        "options": ["--optimize"],
        "requires_acgt": False,
        "repetitions": None,
        "origin": "frozen",
    },
    {
        "id": "hypergrass-lean",
        "kind": "hgs",
        "heap": "8m",
        "options": ["--lean"],
        "requires_acgt": False,
        "repetitions": None,
        "origin": "frozen",
    },
    {"id": "gzip-9", "kind": "gzip", "requires_acgt": False, "repetitions": None, "origin": "frozen"},
    {"id": "xz-9", "kind": "xz", "requires_acgt": False, "repetitions": None, "origin": "frozen"},
    {
        "id": "jarvis3-l7",
        "kind": "jarvis3",
        "level": "7",
        "requires_acgt": True,
        "repetitions": None,
        "origin": "frozen",
    },
    {
        "id": "jarvis3-l26",
        "kind": "jarvis3",
        "level": "26",
        "requires_acgt": True,
        "repetitions": None,
        "origin": "frozen",
    },
    {
        "id": "geco3-l9",
        "kind": "geco3",
        "level": "9",
        "requires_acgt": True,
        "repetitions": 1,
        "origin": "frozen",
    },
    {
        "id": "geco3-l6",
        "kind": "geco3",
        "level": "6",
        "requires_acgt": True,
        "repetitions": None,
        "origin": "supplementary",
    },
    {
        "id": "geco3-l3",
        "kind": "geco3",
        "level": "3",
        "requires_acgt": True,
        "repetitions": None,
        "origin": "supplementary",
    },
)

DEVIATIONS = (
    "Frozen profile geco3-l9 is measured once per dataset, not three times: a "
    "1,048,576-byte pilot exceeded 120 s and the protocol cap is 180 s. Repeating "
    "an operation that cannot finish adds no information.",
    "geco3-l3 and geco3-l6 are supplementary profiles added after the l9 pilot to "
    "provide feasible GeCo3 baselines. They are marked origin=supplementary and "
    "are not part of the frozen protocol.",
    "JARVIS3 and GeCo3 run only when every sequence byte is A/C/G/T, per the "
    "protocol. gzip, xz, and HyperGraSS run on all sequence inputs.",
    "Every tool consumes the exact .seq bytes. FASTA metadata is not part of this "
    "comparison; complete-FASTA results are reported separately by the PanMAN "
    "benchmark.",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(READ_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def to_container(corpus: Path, path: Path) -> str:
    """Translate a host path under the corpus mount to its container path."""
    relative = Path(path).resolve().relative_to(corpus.resolve())
    return "/work/" + relative.as_posix()


def median_or_none(values: list[float | int | None]) -> float | int | None:
    present = [value for value in values if value is not None]
    return statistics.median(present) if present else None


def select_slices(length: int) -> list[tuple[int, int]]:
    if length <= 0:
        return []
    width = min(128, length)
    middle_start = max(0, min(length - width, length // 2 - width // 2))
    return [(0, width), (middle_start, middle_start + width), (length - width, length)]


def profile_repetitions(profile: dict, default: int) -> int:
    return profile["repetitions"] if profile["repetitions"] is not None else default


def docker_prefix(corpus: Path, image_id: str) -> list[str]:
    return [
        "docker",
        "run",
        "--rm",
        "--pull=never",
        "--platform",
        "linux/amd64",
        "--network=none",
        "--cpus",
        "2",
        "--memory",
        "2g",
        "--mount",
        f"type=bind,source={corpus},target=/work",
        "--workdir",
        "/work",
        "--entrypoint",
        "python3",
        image_id,
        "-B",
        "/work/container_run.py",
    ]


def run_operation(
    corpus: Path,
    image_id: str,
    run_dir: Path,
    operation: str,
    argv: list[str],
    cwd: Path,
    timeout: float,
    stdout_artifact: Path | None = None,
) -> dict:
    """Execute one tool invocation in a fresh locked-container process."""
    op_dir = run_dir / "ops" / operation
    op_dir.mkdir(parents=True, exist_ok=True)
    container_dir = "/work/" + op_dir.relative_to(corpus).as_posix()
    spec = {
        "argv": argv,
        "cwd": to_container(corpus, cwd),
        "stdout_log": f"{container_dir}/stdout.log",
        "stderr_log": f"{container_dir}/stderr.log",
        "timeout_seconds": timeout,
        "metrics_path": f"{container_dir}/metrics.json",
    }
    if stdout_artifact is not None:
        spec["stdout_artifact"] = to_container(corpus, stdout_artifact)
    spec_path = op_dir / "spec.json"
    write_json(spec_path, spec)

    docker = docker_prefix(corpus, image_id) + [
        "/work/" + spec_path.relative_to(corpus).as_posix()
    ]
    host = run_process(
        docker,
        cwd=corpus,
        stdout_path=op_dir / "docker.stdout.log",
        stderr_path=op_dir / "docker.stderr.log",
        timeout_seconds=timeout + 90.0,
    )
    metrics_path = op_dir / "metrics.json"
    metrics = load_json(metrics_path) if metrics_path.exists() else None
    result = {
        "operation": operation,
        "argv": argv,
        "docker_status": host["status"],
        "docker_returncode": host["returncode"],
        "metrics": metrics,
        "metrics_present": metrics is not None,
        "logs": {
            "stdout": (op_dir / "stdout.log").relative_to(run_dir).as_posix(),
            "stderr": (op_dir / "stderr.log").relative_to(run_dir).as_posix(),
            "docker_stdout": (op_dir / "docker.stdout.log").relative_to(run_dir).as_posix(),
            "docker_stderr": (op_dir / "docker.stderr.log").relative_to(run_dir).as_posix(),
        },
    }
    if metrics is None:
        result["status"] = "error"
        result["error"] = host.get("error") or "container produced no metrics file"
    else:
        result["status"] = metrics["status"]
        result["wall_seconds"] = metrics["wall_seconds"]
        result["cpu_seconds"] = metrics["cpu_seconds"]
        result["peak_rss_bytes"] = metrics["peak_rss_bytes"]
        result["error"] = metrics.get("error")
    return result


def hgs_java(classes_container_path: str, heap: str) -> list[str]:
    return [
        "java",
        f"-Xmx{heap}",
        "-XX:ActiveProcessorCount=1",
        "-cp",
        classes_container_path,
        "hypergrass.HyperGraSS",
    ]


def compress_and_verify(
    corpus: Path,
    image_id: str,
    run_dir: Path,
    dataset: dict,
    profile: dict,
    repetition: int,
    timeout: float,
    input_path: Path,
    classes: Path,
) -> dict:
    base = f"{dataset['id']}/{profile['id']}/rep-{repetition}"
    work = run_dir / "work" / base
    work.mkdir(parents=True, exist_ok=True)
    seq_path = corpus / dataset["sequence"]["path"]
    input_size = dataset["sequence"]["size_bytes"]
    input_hash = dataset["sequence"]["sha256"]
    kind = profile["kind"]

    artifact: Path
    decompressed: Path
    compress_argv: list[str]
    decompress_argv: list[str]
    compress_artifact_stdout = False
    decompress_artifact_stdout = False
    decompress_cwd = work
    work_c = to_container(corpus, work)
    input_c = to_container(corpus, input_path)
    classes_c = to_container(corpus, classes)

    if kind == "hgs":
        artifact = work / "archive.hgs"
        decompressed = work / "restored.seq"
        compress_argv = hgs_java(classes_c, profile["heap"]) + [
            "compress",
            *profile["options"],
            input_c,
            to_container(corpus, artifact),
        ]
        decompress_argv = hgs_java(classes_c, profile["heap"]) + [
            "decompress",
            to_container(corpus, artifact),
            to_container(corpus, decompressed),
        ]
    elif kind in ("gzip", "xz"):
        tool = "gzip" if kind == "gzip" else "xz"
        artifact = work / ("archive.gz" if kind == "gzip" else "archive.xz")
        decompressed = work / "restored.seq"
        compress_argv = [tool, "-9", "-c", input_c]
        decompress_argv = [tool, "-dc", to_container(corpus, artifact)]
        compress_artifact_stdout = True
        decompress_artifact_stdout = True
    elif kind == "jarvis3":
        artifact = work / "archive.jc"
        decompressed = work / "restored.seq"
        compress_argv = [
            "JARVIS3",
            "-l",
            profile["level"],
            "-o",
            to_container(corpus, artifact),
            input_c,
        ]
        decompress_argv = [
            "JARVIS3",
            "-d",
            "-o",
            to_container(corpus, decompressed),
            to_container(corpus, artifact),
        ]
    elif kind == "geco3":
        local_input = work / "input.seq"
        shutil.copyfile(seq_path, local_input)
        artifact = work / "input.seq.co"
        decompressed = work / "input.seq.de"
        compress_argv = ["GeCo3", "-l", profile["level"], to_container(corpus, local_input)]
        decompress_argv = ["GeDe3", to_container(corpus, artifact)]
    else:
        raise ValueError(f"unknown profile kind: {kind}")

    started_utc = datetime.now(timezone.utc).isoformat()
    compression = run_operation(
        corpus,
        image_id,
        run_dir,
        f"{base}/compress",
        compress_argv,
        work,
        timeout,
        stdout_artifact=artifact if compress_artifact_stdout else None,
    )
    entry = {
        "dataset": dataset["id"],
        "profile": profile["id"],
        "profile_origin": profile["origin"],
        "repetition": repetition,
        "input_contract": "sequence",
        "input_path": dataset["sequence"]["path"],
        "input_bytes": input_size,
        "input_sha256": input_hash,
        "started_utc": started_utc,
        "compress": compression,
        "status": "error",
        "roundtrip_ok": False,
        "artifact_bytes": None,
        "artifact_sha256": None,
        "ratio": None,
        "bits_per_base": None,
        "wall_seconds": compression.get("wall_seconds"),
        "cpu_seconds": compression.get("cpu_seconds"),
        "peak_rss_bytes": compression.get("peak_rss_bytes"),
        "queries": [],
    }

    if compression["status"] != "success" or not artifact.exists():
        entry["error"] = f"compression {compression['status']}"
        return entry

    entry["artifact_bytes"] = artifact.stat().st_size
    entry["artifact_sha256"] = sha256_file(artifact)
    if entry["artifact_bytes"] <= 0:
        entry["error"] = "empty artifact"
        return entry
    entry["ratio"] = input_size / entry["artifact_bytes"]
    entry["bits_per_base"] = entry["artifact_bytes"] * 8.0 / input_size

    decompression = run_operation(
        corpus,
        image_id,
        run_dir,
        f"{base}/decompress",
        decompress_argv,
        decompress_cwd,
        timeout,
        stdout_artifact=decompressed if decompress_artifact_stdout else None,
    )
    entry["decompress"] = decompression
    if decompression["status"] != "success" or not decompressed.exists():
        entry["error"] = f"decompression {decompression['status']}"
        return entry

    restored_hash = sha256_file(decompressed)
    restored_size = decompressed.stat().st_size
    entry["restored_bytes"] = restored_size
    entry["restored_sha256"] = restored_hash
    entry["roundtrip_ok"] = restored_hash == input_hash and restored_size == input_size
    if not entry["roundtrip_ok"]:
        entry["error"] = "round-trip bytes differ"
        return entry

    if kind == "hgs":
        sequence = seq_path.read_bytes()
        for index, (start, end) in enumerate(select_slices(len(sequence))):
            query_file = work / f"query-{index}.seq"
            query = run_operation(
                corpus,
                image_id,
                run_dir,
                f"{base}/query-{index}",
                hgs_java(classes_c, profile["heap"])
                + ["query", to_container(corpus, artifact), "0", str(start), str(end)],
                work,
                timeout,
                stdout_artifact=query_file,
            )
            expected = sequence[start:end]
            observed = query_file.read_bytes() if query_file.exists() else b""
            exact = observed.rstrip(b"\n") == expected
            entry["queries"].append(
                {
                    "start": start,
                    "end": end,
                    "expected_sha256": hashlib.sha256(expected).hexdigest(),
                    "observed_sha256": hashlib.sha256(observed).hexdigest(),
                    "exact_after_trailing_newline": exact,
                    "status": query["status"],
                }
            )
            if query["status"] != "success" or not exact:
                entry["error"] = f"slice query {index} failed"
                return entry

    entry["status"] = "success"
    return entry


def summarise(results: list[dict], datasets: list[dict], profiles: list[dict]) -> list[dict]:
    rows = []
    for dataset in datasets:
        for profile in profiles:
            group = [
                item
                for item in results
                if item["dataset"] == dataset["id"] and item["profile"] == profile["id"]
            ]
            if not group:
                continue
            if all(item["status"] == "not_applicable" for item in group):
                rows.append({
                    "dataset": dataset["id"],
                    "profile": profile["id"],
                    "origin": profile["origin"],
                    "attempts": len(group),
                    "successful": 0,
                    "artifact_bytes": None,
                    "nondeterministic_size": False,
                    "input_bytes": dataset["sequence"]["size_bytes"],
                    "roundtrip_all_ok": False,
                    "median_wall_seconds": None,
                    "median_cpu_seconds": None,
                    "median_peak_rss_bytes": None,
                    "ratio": None,
                    "bits_per_base": None,
                    "status": "not_applicable",
                })
                continue
            applicable = [item for item in group if item["status"] != "not_applicable"]
            successful = [item for item in applicable if item["status"] == "success"]
            sizes = sorted({item["artifact_bytes"] for item in successful if item["artifact_bytes"]})
            row = {
                "dataset": dataset["id"],
                "profile": profile["id"],
                "origin": profile["origin"],
                "attempts": len(applicable),
                "successful": len(successful),
                "artifact_bytes": sizes[0] if len(sizes) == 1 else None,
                "nondeterministic_size": len(sizes) > 1,
                "input_bytes": dataset["sequence"]["size_bytes"],
                "roundtrip_all_ok": bool(successful) and all(item["roundtrip_ok"] for item in successful),
                "median_wall_seconds": median_or_none([item.get("wall_seconds") for item in successful]),
                "median_cpu_seconds": median_or_none([item.get("cpu_seconds") for item in successful]),
                "median_peak_rss_bytes": median_or_none([item.get("peak_rss_bytes") for item in successful]),
                "status": "success" if len(successful) == len(applicable) else (
                    "timeout" if any(item.get("compress", {}).get("status") == "timeout" for item in applicable) else "failed"
                ),
            }
            if row["artifact_bytes"]:
                row["ratio"] = row["input_bytes"] / row["artifact_bytes"]
                row["bits_per_base"] = row["artifact_bytes"] * 8.0 / row["input_bytes"]
            else:
                row["ratio"] = None
                row["bits_per_base"] = None
            rows.append(row)
    return rows


def verify_environment(corpus: Path, classpath: Path) -> tuple[str, dict]:
    data_lock = load_json(corpus / "data.lock.json")
    image_lock = load_json(corpus / "image.lock.json")
    image_id = image_lock.get("image_id")
    if not image_id or not image_id.startswith("sha256:"):
        raise SystemExit("image.lock.json has no immutable image_id")
    inspect = subprocess.run(
        ["docker", "image", "inspect", image_id, "--format", "{{.Os}}/{{.Architecture}}"],
        capture_output=True,
        text=True,
        timeout=60,
    )
    if inspect.returncode != 0:
        raise SystemExit(f"locked image is unavailable: {inspect.stderr.strip()}")
    if inspect.stdout.strip() != "linux/amd64":
        raise SystemExit(f"locked image platform mismatch: {inspect.stdout.strip()}")
    if not classpath.is_dir() or not any(classpath.glob("hypergrass/*.class")):
        raise SystemExit(f"compiled classes not found: {classpath}")
    classes = []
    for path in sorted(classpath.glob("hypergrass/*.class")):
        classes.append(
            {
                "path": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    environment = {
        "image_id": image_id,
        "image_platform": inspect.stdout.strip(),
        "image_versions": image_lock.get("versions"),
        "data_protocol_sha256": data_lock.get("protocol_sha256"),
        "datasets": len(data_lock.get("datasets", [])),
        "hypergrass_classes": classes,
    }
    return image_id, environment


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classpath", default=None)
    parser.add_argument("--repetitions", type=int, default=None)
    parser.add_argument("--timeout", type=float, default=None)
    parser.add_argument("--profiles", default=None, help="comma-separated profile ids")
    parser.add_argument("--datasets", default=None, help="comma-separated dataset ids")
    parser.add_argument("--run-id", default=None)
    args = parser.parse_args(argv)

    corpus = Path(__file__).resolve().parent
    project = corpus.parent.parent
    classpath = Path(args.classpath).resolve() if args.classpath else project / "bin"
    protocol = load_json(corpus / "protocol.json")
    data_lock = load_json(corpus / "data.lock.json")
    image_id, environment = verify_environment(corpus, classpath)

    repetitions = args.repetitions or protocol["measurement"]["repetitions"]
    timeout = args.timeout or protocol["resources"]["timeout_seconds_per_operation"]
    if repetitions < 1 or timeout <= 0:
        raise SystemExit("repetitions and timeout must be positive")

    selected_profiles = [
        profile
        for profile in PROFILES
        if args.profiles is None or profile["id"] in args.profiles.split(",")
    ]
    if not selected_profiles:
        raise SystemExit("no profiles selected")
    by_id = {dataset["id"]: dataset for dataset in data_lock["datasets"]}
    if args.datasets:
        missing = [item for item in args.datasets.split(",") if item not in by_id]
        if missing:
            raise SystemExit(f"unknown datasets: {missing}")
        datasets = [by_id[item] for item in args.datasets.split(",")]
    else:
        datasets = data_lock["datasets"]

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = args.run_id or f"{stamp}-{random.randrange(16**6):06x}"
    run_dir = corpus / "runs" / run_id
    if run_dir.exists():
        raise SystemExit(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)
    classes = run_dir / "classes"
    shutil.copytree(classpath / "hypergrass", classes / "hypergrass")
    environment["staged_classes"] = sorted(
        str(path.relative_to(run_dir)) for path in classes.rglob("*.class")
    )
    write_json(run_dir / "spec.json", {
        "run_id": run_id,
        "protocol_sha256": environment["data_protocol_sha256"],
        "image_id": image_id,
        "repetitions": repetitions,
        "timeout_seconds": timeout,
        "profiles": [profile["id"] for profile in selected_profiles],
        "datasets": [dataset["id"] for dataset in datasets],
        "started_utc": datetime.now(timezone.utc).isoformat(),
    })

    report = {
        "schema_version": SCHEMA_VERSION,
        "status": "running",
        "run_id": run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": protocol,
        "environment": environment,
        "deviations": list(DEVIATIONS),
        "results": [],
        "summary": [],
    }
    write_json(run_dir / "report.json", report)

    def log(message: str) -> None:
        line = f"{datetime.now(timezone.utc).isoformat()} {message}\n"
        with (run_dir / "run.log").open("a", encoding="utf-8") as stream:
            stream.write(line)
        print(line, end="", flush=True)

    log(f"run {run_id}: {len(datasets)} datasets x {len(selected_profiles)} profiles")
    failures = 0
    for dataset in datasets:
        seq_path = corpus / dataset["sequence"]["path"]
        raw_acgt = bool(dataset["sequence"]["raw_acgt"])
        for repetition in range(repetitions):
            rotation = (datasets.index(dataset) + repetition) % len(selected_profiles)
            order = selected_profiles[rotation:] + selected_profiles[:rotation]
            for profile in order:
                if profile["repetitions"] is not None and repetition >= profile["repetitions"]:
                    continue
                if profile["requires_acgt"] and not raw_acgt:
                    report["results"].append({
                        "dataset": dataset["id"],
                        "profile": profile["id"],
                        "profile_origin": profile["origin"],
                        "repetition": repetition,
                        "status": "not_applicable",
                        "reason": "sequence contains non-ACGT bytes; DNA-specialized tool not applicable",
                    })
                    write_json(run_dir / "report.json", {**report, "results": report["results"]})
                    continue
                log(f"{dataset['id']} / {profile['id']} / rep {repetition}")
                try:
                    entry = compress_and_verify(
                        corpus,
                        image_id,
                        run_dir,
                        dataset,
                        profile,
                        repetition,
                        timeout,
                        seq_path,
                        classes,
                    )
                except Exception as error:  # retain the failure, keep the run going
                    entry = {
                        "dataset": dataset["id"],
                        "profile": profile["id"],
                        "profile_origin": profile["origin"],
                        "repetition": repetition,
                        "status": "error",
                        "error": f"{type(error).__name__}: {error}",
                    }
                report["results"].append(entry)
                if entry["status"] != "success":
                    failures += 1
                    log(f"  FAILED: {entry.get('error')}")
                report["summary"] = summarise(report["results"], datasets, selected_profiles)
                write_json(run_dir / "report.json", report)

    report["summary"] = summarise(report["results"], datasets, selected_profiles)
    report["status"] = "success" if failures == 0 else "partial"
    report["failures"] = failures
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(run_dir / "report.json", report)
    log(f"finished: status={report['status']} failures={failures}")
    print(f"report: {run_dir / 'report.json'}")
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
