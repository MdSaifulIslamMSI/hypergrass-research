"""Build the source-verified corpus image; uses only Python's standard library."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import urlsplit


BASE_IMAGE = "eclipse-temurin@sha256:56f10ddc224d7b900698422ffe82432ed0e178f606bd66859ed1d179825dfbe3"
PLATFORM = "linux/amd64"
IMAGE_TAG = "hypergrass-corpus:locked"
SOURCE_COMMITS = {
    "jarvis3": "4dcb37b05b1a1482c2bf864d4968244817b0ad0a",
    "geco3": "69642a0a913e543599f35aee372e300601e4cbbe",
}
CFLAGS = {"jarvis3": "-Wall -Ofast -DPROGRESS", "geco3": "-Wall -O3"}
EXECUTABLES = ("JARVIS3", "GeCo3", "GeDe3", "gcc", "java", "python3", "gzip", "xz", "make", "bash")
VERSION_COMMANDS = {
    "gcc": ["gcc", "--version"],
    "java": ["java", "-version"],
    "python": ["python3", "--version"],
    "gzip": ["gzip", "--version"],
    "xz": ["xz", "--version"],
}
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
CONTEXT_FILES = ("Dockerfile", ".dockerignore", "build_image.py", "sources.lock.json")


class BuildError(RuntimeError):
    pass


def sha256_file(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise BuildError(f"Cannot read {path}: {error}") from error


def validate_sources(root):
    """Check the exact two pinned archives without fetching or modifying anything."""
    root = Path(root).resolve()
    manifest = read_json(root / "sources.lock.json")
    if not isinstance(manifest, dict) or not isinstance(manifest.get("sources"), list):
        raise BuildError("sources.lock.json must contain a sources list")
    sources = manifest["sources"]
    if len(sources) != len(SOURCE_COMMITS):
        raise BuildError("sources.lock.json must contain exactly jarvis3 and geco3")
    verified = {}
    for source in sources:
        if not isinstance(source, dict):
            raise BuildError("Each source lock entry must be an object")
        source_id = source.get("id")
        if not isinstance(source_id, str) or source_id not in SOURCE_COMMITS or source_id in verified:
            raise BuildError(f"Unknown or duplicate source id: {source_id!r}")
        if source.get("commit") != SOURCE_COMMITS[source_id]:
            raise BuildError(f"Incorrect pinned commit for {source_id}")
        relative = f"vendor/{source_id}.tar.gz"
        if source.get("path") != relative:
            raise BuildError(f"{source_id} path must be {relative}")
        url = source.get("url")
        if not isinstance(url, str) or urlsplit(url).scheme != "https" or not urlsplit(url).netloc:
            raise BuildError(f"{source_id} must have an HTTPS provenance URL")
        size = source.get("size_bytes")
        digest = source.get("sha256")
        if type(size) is not int or size <= 0:
            raise BuildError(f"Invalid size_bytes for {source_id}")
        if not isinstance(digest, str) or not SHA256.fullmatch(digest):
            raise BuildError(f"Invalid sha256 for {source_id}")
        archive = root / relative
        if archive.is_symlink() or not archive.resolve().is_relative_to(root) or not archive.is_file():
            raise BuildError(f"Missing or unsafe source archive: {archive}")
        if archive.stat().st_size != size:
            raise BuildError(f"Source size mismatch: {relative}")
        if sha256_file(archive) != digest:
            raise BuildError(f"Source SHA-256 mismatch: {relative}")
        verified[source_id] = {
            key: source[key] for key in ("id", "commit", "url", "path", "size_bytes", "sha256")
        }
    return [verified[key] for key in sorted(verified)]


def build_inputs(root, sources):
    return {
        "base_image": BASE_IMAGE,
        "platform": PLATFORM,
        "cflags": CFLAGS,
        "sources": sources,
        "context_sha256": {name: sha256_file(root / name) for name in CONTEXT_FILES},
    }


def record_environment(output):
    """Executed inside the image, after compiling and installing the tools."""
    sources = validate_sources(Path("/opt/benchmark"))
    flags = {key: os.environ[f"{key.upper()}_CFLAGS"] for key in CFLAGS}
    if flags != CFLAGS:
        raise BuildError("Build flags differ from the portable benchmark flags")
    versions = {}
    for name, command in VERSION_COMMANDS.items():
        versions[name] = subprocess.run(
            command, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", timeout=30,
        ).stdout.strip()
    executables = {}
    for name in EXECUTABLES:
        executable = shutil.which(name)
        if executable is None:
            raise BuildError(f"Missing installed executable: {name}")
        executables[name] = {"path": executable, "sha256": sha256_file(Path(executable))}
    packages = subprocess.run(
        ["apk", "info", "-v"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, encoding="utf-8", timeout=30,
    ).stdout.splitlines()
    environment = {
        "schema_version": 1,
        "base_image": BASE_IMAGE,
        "platform": PLATFORM,
        "alpine_release": Path("/etc/alpine-release").read_text(encoding="utf-8").strip(),
        "sources": sources,
        "cflags": flags,
        "versions": versions,
        "apk_packages": sorted(packages),
        "executables": executables,
    }
    validate_environment(environment, sources)
    Path(output).write_text(json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_environment(environment, sources):
    if not isinstance(environment, dict):
        raise BuildError("Image environment must be an object")
    expected = {
        "schema_version": 1, "base_image": BASE_IMAGE, "platform": PLATFORM,
        "sources": sources, "cflags": CFLAGS,
    }
    for key, value in expected.items():
        if environment.get(key) != value:
            raise BuildError(f"Image environment has incompatible {key}")
    release = environment.get("alpine_release")
    if not isinstance(release, str) or not release.startswith("3.24."):
        raise BuildError("Expected Alpine 3.24 in the pinned image")
    versions = environment.get("versions")
    if not isinstance(versions, dict) or any(
        not isinstance(versions.get(name), str) or not versions[name].strip()
        for name in VERSION_COMMANDS
    ):
        raise BuildError("Image environment is missing tool versions")
    packages = environment.get("apk_packages")
    if not isinstance(packages, list) or not packages or any(
        not isinstance(package, str) or not package.strip() for package in packages
    ):
        raise BuildError("Image environment is missing APK package versions")
    executables = environment.get("executables")
    if not isinstance(executables, dict):
        raise BuildError("Image environment is missing executable hashes")
    for name in EXECUTABLES:
        item = executables.get(name)
        if not isinstance(item, dict) or not isinstance(item.get("sha256"), str) or not SHA256.fullmatch(item["sha256"]):
            raise BuildError(f"Image environment is missing the SHA-256 of {name}")
        path = item.get("path")
        if not isinstance(path, str) or not path.startswith("/"):
            raise BuildError(f"Image environment is missing the installed path of {name}")
        if name in ("JARVIS3", "GeCo3", "GeDe3") and path != f"/usr/local/bin/{name}":
            raise BuildError(f"{name} was not installed in /usr/local/bin")


def docker(command, log, *, capture=False, timeout=60):
    """Log every invocation; never use a shell or merge diagnostics into JSON."""
    log.write("+ " + shlex.join(command) + "\n")
    log.flush()
    result = subprocess.run(
        command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE if capture else log,
        stderr=log, check=True, text=True, encoding="utf-8", errors="replace", timeout=timeout,
    )
    if capture:
        log.write(result.stdout + "\n")
        log.flush()
        return result.stdout
    return ""


def build_image(root):
    root = Path(root).resolve()
    sources = validate_sources(root)
    inputs = build_inputs(root, sources)
    lock_path = root / "image.lock.json"
    existing = read_json(lock_path) if lock_path.exists() else None
    if lock_path.exists() and (
        not isinstance(existing, dict) or existing.get("schema_version") != 1
        or existing.get("tag") != IMAGE_TAG or existing.get("build_inputs") != inputs
    ):
        raise BuildError("Incompatible image.lock.json; preserve/archive it before explicitly replacing it")

    log_path = root / "build.log"
    print(f"Build output: {log_path}", flush=True)
    # Build by ID, not the locked tag, so a rejected rebuild cannot retag it.
    with log_path.open("a", encoding="utf-8", newline="\n") as log:
        log.write("\n=== corpus image build ===\n")
        with tempfile.TemporaryDirectory(prefix="hypergrass-image-") as temporary:
            iidfile = Path(temporary) / "image.id"
            docker([
                "docker", "build", "--platform", PLATFORM, "--network", "default",
                "--progress", "plain", "--iidfile", str(iidfile), str(root),
            ], log, timeout=1800)
            image_id = iidfile.read_text(encoding="utf-8").strip()
        if not image_id.startswith("sha256:") or not SHA256.fullmatch(image_id[7:]):
            raise BuildError("Docker did not produce a valid immutable image ID")
        inspection = json.loads(docker(["docker", "image", "inspect", image_id], log, capture=True))
        if (
            not isinstance(inspection, list) or len(inspection) != 1
            or not isinstance(inspection[0], dict)
            or inspection[0].get("Id") != image_id
            or inspection[0].get("Os") != "linux" or inspection[0].get("Architecture") != "amd64"
        ):
            raise BuildError("Built image is not the expected linux/amd64 image")
        # This short-lived container only reads image metadata: no mounts or network.
        environment = json.loads(docker([
            "docker", "run", "--rm", "--pull=never", "--network=none", "--platform", PLATFORM,
            "--entrypoint", "cat", image_id, "/opt/benchmark-environment.json",
        ], log, capture=True))
        validate_environment(environment, sources)
        candidate = {
            "schema_version": 1,
            "tag": IMAGE_TAG,
            "image_id": image_id,
            "build_inputs": inputs,
            "environment": environment,
        }
        if existing is not None and existing != candidate:
            raise BuildError("Rebuilt image differs from image.lock.json; the lock and tag were not replaced")
        # Recheck the context before accepting the result if another process edited it.
        if build_inputs(root, validate_sources(root)) != inputs:
            raise BuildError("Build inputs changed during the Docker build; refusing to lock the result")
        docker(["docker", "image", "tag", image_id, IMAGE_TAG], log)
        if existing is None:
            # Exclusive creation also refuses a lock created concurrently.
            with lock_path.open("x", encoding="utf-8", newline="\n") as lock_file:
                json.dump(candidate, lock_file, indent=2, sort_keys=True)
                lock_file.write("\n")
    return candidate


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    internal = parser.add_mutually_exclusive_group()
    internal.add_argument("--verify-sources", type=Path, help=argparse.SUPPRESS)
    internal.add_argument("--record-environment", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    try:
        if args.verify_sources is not None:
            validate_sources(args.verify_sources)
        elif args.record_environment is not None:
            record_environment(args.record_environment)
        else:
            root = Path(__file__).resolve().parent
            lock = build_image(root)
            print(f"Locked {IMAGE_TAG} at {lock['image_id']}\nImage lock/environment: {root / 'image.lock.json'}")
    except (BuildError, OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
