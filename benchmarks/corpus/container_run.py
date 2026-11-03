"""Run one benchmark operation inside the locked corpus container.

The caller writes a JSON spec and invokes this file with ``python3 -B``. One
wrapper process executes exactly one child, so ``wait4`` rusage describes that
child alone. The wrapper also records the container cgroup's memory high-water
mark and OOM counter, tails the child's stderr, and classifies failures instead
of reporting every nonzero exit as a generic error.

Spec keys (all paths are container paths):
    argv              list of arguments, required, nonempty
    cwd               working directory, required
    stdout_log        file for child stdout (ignored when stdout_artifact set)
    stderr_log        file for child stderr
    stdout_artifact   optional file receiving child stdout verbatim
    timeout_seconds   positive float; the wrapper kills the child on expiry
    metrics_path      file receiving the JSON metrics object
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from process_metrics import (  # noqa: E402
    classify_exit,
    read_cgroup_memory,
    run_process,
    tail_text,
)


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def failure_message(status: str, returncode: int | None) -> str:
    """Explain a classified failure instead of leaving ``error`` null."""
    if status == "oom-killed":
        return f"Container cgroup recorded an OOM kill (exit {returncode})."
    if status == "oom-heap":
        return f"JVM reported OutOfMemoryError (exit {returncode})."
    if status == "killed":
        return (
            f"Child was killed (exit {returncode}) without cgroup OOM evidence."
        )
    return f"Child exited with status {returncode}."


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: container_run.py SPEC.json", file=sys.stderr)
        return 2
    spec = json.loads(Path(argv[1]).read_text(encoding="utf-8"))
    command = spec["argv"]
    cwd = Path(spec["cwd"])
    timeout = float(spec["timeout_seconds"])
    metrics_path = Path(spec["metrics_path"])

    stdout_log = Path(spec["stdout_log"])
    stderr_log = Path(spec["stderr_log"])
    stdout_artifact = spec.get("stdout_artifact")
    if stdout_artifact:
        stdout_path = Path(stdout_artifact)
        stdout_log.parent.mkdir(parents=True, exist_ok=True)
        stdout_log.write_bytes(b"")
    else:
        stdout_path = stdout_log

    started = datetime.now(timezone.utc).isoformat()
    result = run_process(command, cwd, stdout_path, stderr_log, timeout)
    result["argv"] = command
    result["cwd"] = str(cwd)
    result["started_utc"] = started
    result["finished_utc"] = datetime.now(timezone.utc).isoformat()

    stderr_tail, stderr_bytes = tail_text(stderr_log)
    result["stderr_tail"] = stderr_tail
    result["stderr_bytes"] = stderr_bytes
    result.update(read_cgroup_memory())

    if result["status"] == "error" and result.get("returncode") is not None:
        result["status"] = classify_exit(
            result["returncode"],
            stderr_tail,
            result.get("cgroup_oom_kills"),
        )
    if (
        result["status"] != "success"
        and not result.get("error")
        and result.get("returncode") is not None
    ):
        result["error"] = failure_message(result["status"], result["returncode"])

    write_json(metrics_path, result)
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
