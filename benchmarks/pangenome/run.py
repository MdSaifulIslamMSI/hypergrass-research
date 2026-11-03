"""Scale test: HyperGraSS and baselines on a >100 MB real pangenome.

This is deliberately separate from the frozen held-out corpus. It answers two
questions the corpus is too small to answer:

1. Does ``--lean`` with a small Java heap round-trip a 100 MB+ multi-FASTA?
2. What are the size, time, and peak-RSS costs at that scale?

Every operation runs in a fresh process inside the same locked Linux image.
Docker startup is excluded from the measured wall time by ``container_run.py``.

Usage:
    python -B benchmarks/pangenome/run.py --run-id pangenome-scale-1
"""

from __future__ import annotations

import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "corpus"))
from process_metrics import run_process  # noqa: E402

IMAGE_TAG = "hypergrass-corpus:jdk27"
SCHEMA_VERSION = 2
READ_CHUNK = 1024 * 1024
LOCK_NAME = ".harness-container.lock"

PROFILES = (
    {"id": "hypergrass", "kind": "hgs", "heap": "1g", "options": []},
    {"id": "hypergrass-pangenome", "kind": "hgs", "heap": "1g", "options": ["--pangenome"]},
    {"id": "hypergrass-lean", "kind": "hgs", "heap": "8m", "options": ["--lean"]},
    {"id": "gzip-9", "kind": "gzip"},
    {"id": "xz-9", "kind": "xz"},
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(READ_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def to_container(project: Path, path: Path) -> str:
    return "/work/" + path.resolve().relative_to(project.resolve()).as_posix()


def _capture(argv: list[str], cwd: Path) -> str | None:
    try:
        completed = subprocess.run(
            argv, cwd=cwd, capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def provenance(project: Path) -> dict:
    """Record the code, tool, and host state that produced a report."""
    status = _capture(["git", "status", "--porcelain"], project)
    return {
        "git_commit": _capture(["git", "rev-parse", "HEAD"], project),
        "git_dirty": bool(status) if status is not None else None,
        "git_status": (status or "")[:4000],
        "docker_version": _capture(
            ["docker", "version", "--format", "{{.Server.Version}}"], project,
        ),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }


@contextlib.contextmanager
def container_lock(project: Path, timeout_seconds: float = 900.0):
    """Hold a cross-process lock for the duration of one container operation.

    The lock is advisory and released by the operating system if the harness
    dies, so it cannot become a permanent stale lock. Its purpose is to enforce
    the one-HyperGraSS-container-at-a-time rule across concurrent harnesses.
    """
    path = project / LOCK_NAME
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch(exist_ok=True)
    handle = path.open("r+b")
    try:
        if path.stat().st_size == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        deadline = time.monotonic() + timeout_seconds
        while True:
            try:
                if os.name == "nt":
                    import msvcrt

                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl

                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError:
                if time.monotonic() >= deadline:
                    raise TimeoutError(
                        f"another container operation holds {path}"
                    ) from None
                time.sleep(0.5)
        try:
            yield
        finally:
            if os.name == "nt":
                import msvcrt

                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    finally:
        handle.close()


def _container_name(run_dir: Path, name: str) -> str:
    raw = f"hgs-{run_dir.name}-{name}-{os.getpid()}".lower()
    clean = "".join(
        char if (char.isalnum() or char in "-_.") else "-" for char in raw
    )
    return clean[:63]


def image_id(project: Path) -> str:
    lock = json.loads((project / "benchmarks/corpus/image.lock.json").read_text(encoding="utf-8"))
    value = lock.get("image_id", "")
    if not value.startswith("sha256:"):
        raise SystemExit("image.lock.json has no immutable image_id")
    return value


def docker_prefix(
    project: Path, image: str, memory: str, name: str | None = None,
) -> list[str]:
    flags = [
        "docker", "run", "--rm", "--pull=never", "--platform", "linux/amd64",
        "--network=none", "--cpus", "2",
        "--memory", memory, "--memory-swap", memory,
    ]
    if name:
        flags += ["--name", name]
    return flags + [
        "--mount", f"type=bind,source={project},target=/work",
        "--workdir", "/work", "--entrypoint", "python3", image,
        "-B", "/work/benchmarks/corpus/container_run.py",
    ]


def run_operation(
    project: Path,
    image: str,
    run_dir: Path,
    name: str,
    argv: list[str],
    cwd: Path,
    timeout: float,
    memory: str = "2g",
    stdout_artifact: Path | None = None,
) -> dict:
    op_dir = run_dir / "ops" / name
    op_dir.mkdir(parents=True, exist_ok=True)
    spec = {
        "argv": argv,
        "cwd": to_container(project, cwd),
        "stdout_log": to_container(project, op_dir / "stdout.log"),
        "stderr_log": to_container(project, op_dir / "stderr.log"),
        "timeout_seconds": timeout,
        "metrics_path": to_container(project, op_dir / "metrics.json"),
    }
    if stdout_artifact is not None:
        spec["stdout_artifact"] = to_container(project, stdout_artifact)
    spec_path = op_dir / "spec.json"
    write_json(spec_path, spec)
    container = _container_name(run_dir, name)
    docker = docker_prefix(
        project, image, memory, container,
    ) + [to_container(project, spec_path)]
    with container_lock(project):
        try:
            host = run_process(
                docker,
                cwd=project,
                stdout_path=op_dir / "docker.stdout.log",
                stderr_path=op_dir / "docker.stderr.log",
                timeout_seconds=timeout + 120.0,
            )
        finally:
            # --rm covers a normal exit; this covers a harness timeout, a killed
            # Docker client, or any path that leaves the container behind.
            try:
                subprocess.run(
                    ["docker", "rm", "-f", container],
                    capture_output=True, text=True, timeout=60,
                )
            except (OSError, subprocess.TimeoutExpired):
                pass
    metrics_path = op_dir / "metrics.json"
    metrics = json.loads(metrics_path.read_text(encoding="utf-8")) if metrics_path.exists() else None
    if metrics is None:
        return {
            "status": host.get("status") or "error",
            "returncode": host.get("returncode"),
            "wall_seconds": host.get("wall_seconds"),
            "cpu_seconds": host.get("cpu_seconds"),
            "peak_rss_bytes": host.get("peak_rss_bytes"),
            "error": host.get("error") or "no metrics",
            "container_name": container,
            "memory_limit": memory,
        }
    return {
        "status": metrics["status"],
        "returncode": metrics["returncode"],
        "wall_seconds": metrics["wall_seconds"],
        "cpu_seconds": metrics["cpu_seconds"],
        "peak_rss_bytes": metrics["peak_rss_bytes"],
        "memory_method": metrics.get("memory_method"),
        "cgroup_memory_peak_bytes": metrics.get("cgroup_memory_peak_bytes"),
        "cgroup_oom_kills": metrics.get("cgroup_oom_kills"),
        "cgroup_memory_method": metrics.get("cgroup_memory_method"),
        "error": metrics.get("error"),
        "stderr_tail": metrics.get("stderr_tail"),
        "stderr_bytes": metrics.get("stderr_bytes"),
        "container_name": container,
        "memory_limit": memory,
    }


def parse_fasta(path: Path) -> tuple[list[str], list[int], list[int]]:
    """Return headers, per-record sequence start offsets, and lengths in bytes.

    The offsets are positions in the sequence stream (headers and line endings
    removed), matching HyperGraSS record coordinates.
    """
    headers: list[str] = []
    starts: list[int] = []
    lengths: list[int] = []
    total = 0
    current = -1
    with path.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                headers.append(line[1:].rstrip(b"\r\n").decode("ascii", "replace"))
                current += 1
                starts.append(total)
                lengths.append(0)
            else:
                clean = line.rstrip(b"\r\n")
                total += len(clean)
                lengths[current] += len(clean)
    return headers, starts, lengths


def sequence_slice(path: Path, start: int, length: int) -> bytes:
    """Read a byte range from the concatenated sequence stream."""
    out = bytearray()
    position = 0
    with path.open("rb") as stream:
        for line in stream:
            if line.startswith(b">"):
                continue
            clean = line.rstrip(b"\r\n")
            end = position + len(clean)
            if end > start and position < start + length:
                left = max(0, start - position)
                right = min(len(clean), start + length - position)
                out.extend(clean[left:right])
            position = end
            if position >= start + length:
                break
    return bytes(out)


def run_profile(
    project: Path,
    image: str,
    run_dir: Path,
    fasta: Path,
    digest: str,
    size: int,
    profile: dict,
    timeout: float,
) -> dict:
    name = profile["id"]
    work = run_dir / "work" / name
    work.mkdir(parents=True, exist_ok=True)
    kind = profile["kind"]
    input_c = to_container(project, fasta)
    memory = "256m" if profile.get("heap") == "8m" else "2g"

    if kind == "hgs":
        artifact = work / "archive.hgs"
        restored = work / "restored.fasta"
        compress = (
            ["java", f"-Xmx{profile['heap']}", "-XX:ActiveProcessorCount=1",
             "-cp", "/work/bin", "hypergrass.HyperGraSS", "compress",
             *profile["options"], input_c, to_container(project, artifact)]
        )
        decompress = (
            ["java", f"-Xmx{profile['heap']}", "-XX:ActiveProcessorCount=1",
             "-cp", "/work/bin", "hypergrass.HyperGraSS", "decompress",
             to_container(project, artifact), to_container(project, restored)]
        )
    elif kind in ("gzip", "xz"):
        tool = "gzip" if kind == "gzip" else "xz"
        artifact = work / ("archive.gz" if kind == "gzip" else "archive.xz")
        restored = work / "restored.fasta"
        compress = [tool, "-9", "-c", input_c]
        decompress = [tool, "-dc", to_container(project, artifact)]
    else:
        raise ValueError(kind)

    started = datetime.now(timezone.utc).isoformat()
    entry = {
        "profile": name,
        "profile_kind": kind,
        "started_utc": started,
        "input_bytes": size,
        "input_sha256": digest,
        "status": "error",
    }
    compression = run_operation(
        project, image, run_dir, f"{name}/compress", compress, work, timeout,
        memory=memory, stdout_artifact=artifact if kind in ("gzip", "xz") else None,
    )
    entry["compress"] = compression
    if compression["status"] != "success" or not artifact.exists():
        entry["error"] = f"compression {compression['status']}"
        return entry
    entry["artifact_bytes"] = artifact.stat().st_size
    entry["artifact_sha256"] = sha256_file(artifact)
    entry["ratio"] = size / entry["artifact_bytes"]
    entry["bits_per_byte"] = entry["artifact_bytes"] * 8.0 / size
    entry["wall_seconds"] = compression["wall_seconds"]
    entry["cpu_seconds"] = compression["cpu_seconds"]
    entry["peak_rss_bytes"] = compression["peak_rss_bytes"]

    decompression = run_operation(
        project, image, run_dir, f"{name}/decompress", decompress, work, timeout,
        memory=memory, stdout_artifact=restored if kind in ("gzip", "xz") else None,
    )
    entry["decompress"] = decompression
    if decompression["status"] != "success" or not restored.exists():
        entry["error"] = f"decompression {decompression['status']}"
        return entry
    entry["restored_bytes"] = restored.stat().st_size
    entry["restored_sha256"] = sha256_file(restored)
    entry["roundtrip_ok"] = entry["restored_sha256"] == digest and entry["restored_bytes"] == size
    if not entry["roundtrip_ok"]:
        entry["error"] = "round-trip bytes differ"
        return entry

    if kind == "hgs":
        headers, starts, lengths = parse_fasta(fasta)
        probes = []
        indexes = sorted({0, len(headers) // 2, len(headers) - 1})
        for index in indexes:
            start = starts[index] + min(1000, max(0, lengths[index] - 512))
            length = min(512, lengths[index] - (start - starts[index]))
            expected = sequence_slice(fasta, start, length)
            query_file = work / f"query-{index}.seq"
            query = run_operation(
                project, image, run_dir, f"{name}/query-{index}",
                ["java", f"-Xmx{profile['heap']}", "-XX:ActiveProcessorCount=1",
                 "-cp", "/work/bin", "hypergrass.HyperGraSS", "query",
                 to_container(project, artifact), str(index), str(start - starts[index]),
                 str(start - starts[index] + length)],
                work, timeout, memory=memory, stdout_artifact=query_file,
            )
            observed = query_file.read_bytes().rstrip(b"\n") if query_file.exists() else b""
            probes.append({
                "record": index,
                "offset": start - starts[index],
                "length": length,
                "status": query["status"],
                "exact": observed == expected,
            })
            if query["status"] != "success" or observed != expected:
                entry["error"] = f"query record {index} failed"
                return entry
        entry["queries"] = probes

    entry["status"] = "success"
    return entry


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--profiles", default=None)
    args = parser.parse_args(argv)

    project = Path(__file__).resolve().parent.parent.parent
    pangenome = Path(__file__).resolve().parent
    manifest = json.loads((pangenome / "manifest.json").read_text(encoding="utf-8"))
    fasta = pangenome / manifest["fasta"]["path"]
    if not fasta.exists():
        raise SystemExit(f"missing pangenome FASTA: {fasta}")
    if fasta.stat().st_size != manifest["fasta"]["size_bytes"]:
        raise SystemExit("pangenome size does not match manifest")
    digest = sha256_file(fasta)
    if digest != manifest["fasta"]["sha256"]:
        raise SystemExit("pangenome SHA-256 does not match manifest")
    if fasta.stat().st_size < 100_000_000:
        raise SystemExit("scale test requires a 100 MB+ pangenome")

    image = image_id(project)
    inspect = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Os}}/{{.Architecture}}"],
        capture_output=True, text=True, timeout=60,
    )
    if inspect.returncode != 0 or inspect.stdout.strip() != "linux/amd64":
        raise SystemExit("locked image is not available or not linux/amd64")
    if not (project / "bin/hypergrass/HyperGraSS.class").exists():
        raise SystemExit("compiled classes not found; run build.ps1 or build.sh")

    profiles = [p for p in PROFILES if args.profiles is None or p["id"] in args.profiles.split(",")]
    if not profiles:
        raise SystemExit("no profiles selected")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = args.run_id or f"pangenome-{stamp}"
    run_dir = pangenome / "runs" / run_id
    if run_dir.exists():
        raise SystemExit(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)

    report = {
        "schema_version": SCHEMA_VERSION,
        "kind": "scale-test",
        "run_id": run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "image_id": image,
        "image_tag": IMAGE_TAG,
        "provenance": provenance(project),
        "manifest": manifest,
        "fasta": str(fasta.relative_to(project)),
        "profiles": [p["id"] for p in profiles],
        "results": [],
        "status": "running",
    }
    write_json(run_dir / "report.json", report)
    print(f"run {run_id}: {len(profiles)} profiles on {manifest['fasta']['size_bytes']:,} bytes")

    failures = 0
    for profile in profiles:
        print(f"  {profile['id']} ...", flush=True)
        try:
            entry = run_profile(
                project, image, run_dir, fasta, digest,
                manifest["fasta"]["size_bytes"], profile, args.timeout,
            )
        except Exception as error:
            entry = {"profile": profile["id"], "status": "error", "error": f"{type(error).__name__}: {error}"}
        report["results"].append(entry)
        if entry["status"] != "success":
            failures += 1
            print(f"    FAILED: {entry.get('error')}", flush=True)
        else:
            print(
                f"    {entry['artifact_bytes']:,} B, {entry['ratio']:.2f}x, "
                f"{entry['wall_seconds']:.1f}s, {entry['peak_rss_bytes'] / 1048576:.0f} MiB RSS",
                flush=True,
            )
        write_json(run_dir / "report.json", report)

    report["status"] = "success" if failures == 0 else "partial"
    report["failures"] = failures
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(run_dir / "report.json", report)
    print(f"report: {run_dir / 'report.json'}")
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
