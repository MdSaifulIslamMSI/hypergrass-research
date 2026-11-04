"""One-off smoke artifact; run inside the locked image with only corpus mounted."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
SOURCE = b"ACGT" * 1024


def digest(content):
    return hashlib.sha256(content).hexdigest()


def run(command, directory, label):
    print("COMMAND " + json.dumps({"argv": command, "cwd": str(directory)}), flush=True)
    try:
        result = subprocess.run(
            command, cwd=directory, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180,
        )
    except subprocess.TimeoutExpired as error:
        output = error.stdout or b""
        (ROOT / f"{label}.log").write_bytes(output)
        print(output.decode("utf-8", errors="replace"), flush=True)
        print("TIMEOUT 180 seconds", flush=True)
        raise
    (ROOT / f"{label}.log").write_bytes(result.stdout)
    print(result.stdout.decode("utf-8", errors="replace"), flush=True)
    print("EXIT " + str(result.returncode), flush=True)
    result.check_returncode()


def verify(name, directory, restored):
    original = directory / "input.seq"
    if original.read_bytes() != SOURCE:
        raise AssertionError(f"{name}: original input was modified")
    content = restored.read_bytes()
    print("VERIFY " + json.dumps({
        "profile": name,
        "input": str(original),
        "output": str(restored),
        "input_size_bytes": len(SOURCE),
        "output_size_bytes": len(content),
        "input_sha256": digest(SOURCE),
        "output_sha256": digest(content),
        "byte_equal": content == SOURCE,
    }), flush=True)
    if content != SOURCE:
        raise AssertionError(f"{name}: decompression was not byte-exact")


def main():
    # Exclusive creation prevents a repeat invocation from overwriting evidence.
    with (ROOT / "acgt.seq").open("xb") as stream:
        stream.write(SOURCE)
    if set(SOURCE) != set(b"ACGT"):
        raise AssertionError("Smoke input must contain only A/C/G/T, with no newline")
    print(f"SOURCE size_bytes={len(SOURCE)} sha256={digest(SOURCE)}", flush=True)
    failures = []
    for name, level in (("jarvis3-l7", 7), ("jarvis3-l26", 26), ("geco3-l9", 9)):
        directory = ROOT / name
        directory.mkdir()
        (directory / "input.seq").write_bytes(SOURCE)
        try:
            if name.startswith("jarvis3"):
                run(["JARVIS3", "-l", str(level), "-o", "input.jc", "input.seq"], directory, name + ".compress")
                if not (directory / "input.jc").is_file():
                    raise AssertionError("JARVIS3 did not create the explicit archive")
                run(["JARVIS3", "-d", "-o", "restored.seq", "input.jc"], directory, name + ".decompress")
                restored = directory / "restored.seq"
            else:
                run(["GeCo3", "-l", str(level), "input.seq"], directory, name + ".compress")
                if not (directory / "input.seq.co").is_file():
                    raise AssertionError("GeCo3 did not create input.seq.co")
                before = set(directory.iterdir())
                run(["GeDe3", "input.seq.co"], directory, name + ".decompress")
                created = sorted(set(directory.iterdir()) - before)
                print("GEDE3_CREATED " + json.dumps([str(path) for path in created]), flush=True)
                if len(created) != 1 or not created[0].is_file():
                    raise AssertionError("Expected exactly one implicit GeDe3 output file")
                restored = created[0]
            verify(name, directory, restored)
            print("PASS " + name, flush=True)
        except (OSError, AssertionError, subprocess.SubprocessError) as error:
            failures.append(name)
            print(f"FAIL {name}: {error}", flush=True)
    print("SUMMARY " + json.dumps({"profiles": 3, "failures": failures}), flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
