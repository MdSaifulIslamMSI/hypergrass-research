"""Run one benchmark operation inside the locked corpus container.

The caller writes a JSON spec and invokes this file with ``python3 -B``. The
child is spawned directly (never through a shell), its stdout/stderr go to
disk, and Linux child rusage is recorded after reaping. One wrapper process
executes exactly one child so ``RUSAGE_CHILDREN`` describes that child alone.

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

import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: container_run.py SPEC.json", file=sys.stderr)
        return 2
    spec = json.loads(Path(argv[1]).read_text(encoding="utf-8"))
    command = spec["argv"]
    cwd = Path(spec["cwd"])
    timeout = float(spec["timeout_seconds"])
    metrics_path = Path(spec["metrics_path"])
    metrics_path.parent.mkdir(parents=True, exist_ok=True)

    stdout_log = Path(spec["stdout_log"])
    stderr_log = Path(spec["stderr_log"])
    stdout_artifact = spec.get("stdout_artifact")
    for path in (stdout_log, stderr_log):
        path.parent.mkdir(parents=True, exist_ok=True)

    result = {
        "argv": command,
        "cwd": str(cwd),
        "status": "error",
        "returncode": None,
        "wall_seconds": 0.0,
        "cpu_seconds": None,
        "peak_rss_bytes": None,
        "error": None,
    }
    started = time.perf_counter()
    try:
        with stderr_log.open("wb") as stderr:
            if stdout_artifact:
                artifact = Path(stdout_artifact)
                artifact.parent.mkdir(parents=True, exist_ok=True)
                stdout_stream = artifact.open("wb")
                stdout_log.write_bytes(b"")
            else:
                stdout_stream = stdout_log.open("wb")
            try:
                completed = subprocess.run(
                    command,
                    cwd=cwd,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout_stream,
                    stderr=stderr,
                    timeout=timeout,
                )
                result["returncode"] = completed.returncode
                result["status"] = "success" if completed.returncode == 0 else "nonzero"
            finally:
                stdout_stream.close()
    except subprocess.TimeoutExpired:
        result["status"] = "timeout"
        result["error"] = f"exceeded {timeout:g}s"
    except OSError as error:
        result["error"] = str(error)
    finally:
        result["wall_seconds"] = time.perf_counter() - started
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        result["cpu_seconds"] = usage.ru_utime + usage.ru_stime
        # Linux reports ru_maxrss in KiB.
        result["peak_rss_bytes"] = int(usage.ru_maxrss) * 1024

    metrics_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
