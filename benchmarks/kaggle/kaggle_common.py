"""Shared Kaggle cloud harness for large-scale HGS2/EIMDL validation.

Builds OpenJDK (>=21), compiles HyperGraSS + EIMDL, self-verifies Compression Laws
I-VII and the Fusion Engine, builds pinned JARVIS3 and AGC baselines, prepares
multi-hundred-MB and gigabyte-scale Human (GRCh38 + T2T-CHM13v2.0) and SARS-CoV-2
datasets, and measures wall time, peak RSS, compression ratio, and SHA-256 parity.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request


JARVIS3_COMMIT = "4dcb37b05b1a1482c2bf864d4968244817b0ad0a"
JARVIS3_URL = (
    f"https://codeload.github.com/cobilab/jarvis3/tar.gz/{JARVIS3_COMMIT}"
)
JARVIS3_SHA256 = (
    "03cf90a31333ffb7b0c49cd160aeae974780dd44b07508ed0d322f4d31e70325"
)

CHM13_URL = (
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/"
    "assemblies/analysis_set/chm13v2.0.fa.gz"
)

WORK_ROOT = Path("/tmp/hgs_bench")
BIN_DIR = Path("/tmp/hgs_bin")
RESULTS_DIR = Path("/kaggle/working")


def log(msg: str) -> None:
    ts = time.strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def run_cmd(
    cmd: list[str] | str,
    *,
    cwd: Path | None = None,
    timeout: int = 1800,
    check: bool = True,
    shell: bool = False,
    measure_rss: bool = False,
) -> dict:
    """Runs a command, optionally wrapping with /usr/bin/time -v to capture Max RSS."""
    time_bin = Path("/usr/bin/time")
    actual_cmd = cmd
    rss_log = None
    if measure_rss and time_bin.exists() and not shell and isinstance(cmd, list):
        rss_log = WORK_ROOT / f"time_{os.getpid()}_{time.monotonic_ns()}.txt"
        actual_cmd = [str(time_bin), "-v", "-o", str(rss_log), *cmd]

    t0 = time.monotonic()
    proc = subprocess.run(
        actual_cmd,
        cwd=str(cwd) if cwd else None,
        shell=shell,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    elapsed = time.monotonic() - t0
    peak_rss_kb = None
    if rss_log and rss_log.exists():
        txt = rss_log.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", txt)
        if m:
            peak_rss_kb = int(m.group(1))
        rss_log.unlink(missing_ok=True)

    if check and proc.returncode != 0:
        raise RuntimeError(
            f"Command failed (rc={proc.returncode}): {cmd}\n"
            f"STDOUT:\n{proc.stdout[-2000:]}\nSTDERR:\n{proc.stderr[-2000:]}"
        )
    return {
        "rc": proc.returncode,
        "elapsed": elapsed,
        "peak_rss_kb": peak_rss_kb,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    }


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect_system_info() -> dict:
    info: dict = {
        "cpu_count": os.cpu_count(),
        "python": sys.version.split()[0],
    }
    for key, cmd in [
        ("uname", ["uname", "-a"]),
        ("lscpu", ["lscpu"]),
        ("meminfo", ["head", "-n", "6", "/proc/meminfo"]),
        ("df", ["df", "-h", "/kaggle/working", "/tmp"]),
        ("nvidia_smi", ["nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"]),
    ]:
        try:
            res = run_cmd(cmd, check=False, timeout=30)
            info[key] = res["stdout"].strip() if res["rc"] == 0 else res["stderr"].strip()
        except Exception as exc:
            info[key] = f"unavailable: {exc}"
    return info


def _parse_java_major(version_output: str) -> int:
    m = re.search(r'"(\d+)(?:\.\d+)?', version_output)
    if m:
        return int(m.group(1))
    m2 = re.search(r"javac\s+(\d+)", version_output)
    return int(m2.group(1)) if m2 else 0


def ensure_jdk_and_tools(src_bundle_root: Path) -> dict:
    """Installs apt utilities + OpenJDK (>=21), compiles HGS2, and builds JARVIS3/AGC."""
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    BIN_DIR.mkdir(parents=True, exist_ok=True)

    log("Installing Linux benchmark dependencies (time, zstd, xz-utils, build tools, rustc/cargo)...")
    run_cmd(
        "apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq "
        "time zstd xz-utils build-essential cmake git curl zlib1g-dev pigz cargo rustc",
        shell=True,
        check=False,
        timeout=600,
    )

    java_bin = shutil.which("java")
    javac_bin = shutil.which("javac")
    major = 0
    if javac_bin:
        out = run_cmd([javac_bin, "-version"], check=False)
        major = _parse_java_major(out["stdout"] + "\n" + out["stderr"])

    if major < 21:
        log(f"System javac major={major} < 21; downloading Eclipse Temurin JDK 21 LTS...")
        jdk_tgz = WORK_ROOT / "temurin21.tar.gz"
        jdk_dir = WORK_ROOT / "jdk21"
        existing_java = list(jdk_dir.glob("*/bin/java"))
        if not existing_java:
            url = (
                "https://api.adoptium.net/v3/binary/latest/21/ga/linux/x64/"
                "jdk/hotspot/normal/eclipse?project=jdk"
            )
            urllib.request.urlretrieve(url, jdk_tgz)
            jdk_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(jdk_tgz, "r:gz") as tf:
                tf.extractall(jdk_dir)
            jdk_tgz.unlink(missing_ok=True)
            existing_java = list(jdk_dir.glob("*/bin/java"))
        java_bin = str(existing_java[0])
        javac_bin = str(existing_java[0].parent / "javac")

    java_ver = run_cmd([java_bin, "-version"], check=True)
    log(f"Using Java: {(java_ver['stdout'] + java_ver['stderr']).strip().splitlines()[0]}")

    src_tgz_list = list(src_bundle_root.rglob("hypergrass_eimdl_src.tar.gz"))
    if src_tgz_list:
        src_extract = WORK_ROOT / "src_repo"
        if src_extract.exists():
            shutil.rmtree(src_extract)
        src_extract.mkdir(parents=True, exist_ok=True)
        with tarfile.open(src_tgz_list[0], "r:gz") as tf:
            tf.extractall(src_extract)
        java_dir = src_extract / "src" / "hypergrass"
    else:
        cli_java = next(src_bundle_root.rglob("EimdlCli.java"))
        java_dir = cli_java.parent

    java_files = sorted(str(p) for p in java_dir.glob("*.java"))
    log(f"Compiling {len(java_files)} HyperGraSS/EIMDL Java sources (--release 21 -Xlint:all)...")
    run_cmd([javac_bin, "--release", "21", "-Xlint:all", "-d", str(BIN_DIR), *java_files], timeout=300)

    grass_fixture = next(src_bundle_root.rglob("GCF_001884535.1.fna"), None)
    edge_jvm = [f"-Dgrass.src={grass_fixture.parent}"] if grass_fixture else []

    suites = [
        ("VerifyAllEdgeCases", edge_jvm, 300),
        ("VerifyDonorIndex", [], 300),
        ("VerifyAnchorCodecs", ["-Xmx512m"], 300),
        ("VerifyArchiveAccounting", ["-Xmx512m"], 300),
        ("VerifyDonorCosts", [], 300),
        ("VerifyCompressionLaws", ["-Xmx2g"], 300),
        ("VerifyFusionEngine", ["-Xmx2g"], 300),
    ]
    suite_results: dict = {}
    for cls_name, jvm_flags, suite_timeout in suites:
        log(f"Running cloud self-test: {cls_name}...")
        sres = run_cmd(
            [java_bin, *jvm_flags, "-cp", str(BIN_DIR), f"hypergrass.{cls_name}"],
            cwd=WORK_ROOT,
            timeout=suite_timeout,
        )
        lines = [ln.strip() for ln in sres["stdout"].splitlines() if ln.strip()]
        suite_results[cls_name] = {
            "status": "ok",
            "elapsed_s": round(sres["elapsed"], 3),
            "tail": lines[-1] if lines else "",
        }

    tools: dict = {
        "java": java_bin,
        "javac": javac_bin,
        "java_version": (java_ver["stdout"] + java_ver["stderr"]).strip().splitlines()[0],
        "verification_suites": suite_results,
        "verify_laws_tail": suite_results["VerifyCompressionLaws"]["tail"],
        "verify_fusion_tail": suite_results["VerifyFusionEngine"]["tail"],
        "eimdl_native": None,
        "jarvis3": None,
        "agc": None,
        "zstd": shutil.which("zstd"),
        "xz": shutil.which("xz"),
    }

    # Build Native C11/C++20/Rust engine (eimdl)
    native_dir = java_dir.parent.parent / "native" / "eimdl-rs"
    if native_dir.is_dir():
        try:
            cargo_bin = shutil.which("cargo")
            if not cargo_bin:
                log("Cargo missing via apt; installing minimal rustup toolchain...")
                run_cmd(
                    "curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal",
                    shell=True,
                    timeout=300,
                )
                cargo_cand = Path.home() / ".cargo" / "bin" / "cargo"
                if cargo_cand.exists():
                    cargo_bin = str(cargo_cand)
                    os.environ["PATH"] = f"{cargo_cand.parent}:{os.environ.get('PATH', '')}"
            if cargo_bin:
                run_cmd([cargo_bin, "build", "--release", "--offline"], cwd=native_dir, timeout=300)
                run_cmd([cargo_bin, "test", "--release", "--offline"], cwd=native_dir, timeout=300)
                nb = native_dir / "target" / "release" / "eimdl"
                if nb.exists():
                    tools["eimdl_native"] = str(nb)
                    log(f"Built & verified Native C11/C++20/Rust engine at {nb}")
        except Exception as exc:
            log(f"Warning: Native C/C++/Rust engine build failed: {exc}")

    # Build pinned JARVIS3
    try:
        j_tgz = WORK_ROOT / "jarvis3.tar.gz"
        urllib.request.urlretrieve(JARVIS3_URL, j_tgz)
        got = sha256_file(j_tgz)
        if got != JARVIS3_SHA256:
            raise RuntimeError(f"JARVIS3 SHA-256 mismatch: {got}")
        j_dir = WORK_ROOT / "jarvis3"
        j_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(j_tgz, "r:gz") as tf:
            tf.extractall(j_dir)
        inner_src = next(j_dir.glob("*/src"))
        run_cmd(["make", "-C", str(inner_src), "CC=gcc", "CFLAGS=-Wall -Ofast -DPROGRESS"], timeout=300)
        j_bin = inner_src / "JARVIS3"
        if j_bin.exists():
            tools["jarvis3"] = str(j_bin)
            log(f"Built pinned JARVIS3 at {j_bin}")
    except Exception as exc:
        log(f"Warning: JARVIS3 build failed: {exc}")

    # Build AGC (refresh-bio/agc)
    try:
        agc_dir = WORK_ROOT / "agc_repo"
        if not agc_dir.exists():
            run_cmd(
                ["git", "clone", "--depth", "1", "--recurse-submodules", "https://github.com/refresh-bio/agc", str(agc_dir)],
                timeout=300,
            )
        cpus = max(2, min(16, os.cpu_count() or 4))
        run_cmd(["make", "-C", str(agc_dir), f"-j{cpus}"], timeout=600)
        for cand in [agc_dir / "bin" / "agc", agc_dir / "agc"]:
            if cand.exists():
                tools["agc"] = str(cand)
                log(f"Built AGC at {cand}")
                break
    except Exception as exc:
        log(f"Warning: AGC build failed: {exc}")

    return tools


def discover_kaggle_datasets() -> dict[str, Path | None]:
    """Finds attached source bundle, GRCh38 DNA, and COVID-19 datasets under /kaggle/input."""
    inp = Path("/kaggle/input")
    found: dict[str, Path | None] = {"src": None, "grch38": None, "covid19": None}
    for hit in sorted(inp.rglob("kaggle_common.py")):
        found["src"] = hit.parent
        break
    if found["src"] is None:
        for hit in sorted(inp.rglob("EimdlCli.java")):
            found["src"] = hit.parent
            break
    for hit in sorted(inp.rglob("Homo_sapiens.GRCh38.dna.chromosome.1.fa")):
        found["grch38"] = hit.parent
        break
    for pat in ("group_1.fasta", "group_22.fasta"):
        hits = sorted(inp.rglob(pat))
        if hits:
            found["covid19"] = hits[0].parent
            break
    return found


def fasta_to_flat_seq(fasta_path: Path, out_seq: Path) -> dict:
    """Streams FASTA into canonical uppercase ACGT flat .seq (non-ACGT -> A) identically to run_trio.py."""
    table = bytes((v if v in b"ACGT" else ord("A")) for v in range(256))
    total = 0
    non_acgt = 0
    with fasta_path.open("rb") as fin, out_seq.open("wb") as fout:
        for line in fin:
            if line.startswith(b">"):
                continue
            s = line.strip()
            if not s:
                continue
            u = s.upper()
            acgt_cnt = u.count(b"A") + u.count(b"C") + u.count(b"G") + u.count(b"T")
            non_acgt += len(u) - acgt_cnt
            mapped = u.translate(table)
            fout.write(mapped)
            total += len(mapped)
    return {
        "seq_path": str(out_seq),
        "bases": total,
        "non_acgt_masked": non_acgt,
        "sha256": sha256_file(out_seq),
    }


def download_and_split_chm13(dest_dir: Path) -> dict[str, Path]:
    """Downloads T2T-CHM13v2.0.fa.gz and splits chr1..chr22, chrX, chrY cleanly."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    marker = dest_dir / "_done.json"
    if marker.exists():
        mapping = json.loads(marker.read_text(encoding="utf-8"))
        return {k: Path(v) for k, v in mapping.items()}

    gz_path = WORK_ROOT / "chm13v2.0.fa.gz"
    log("Downloading T2T-CHM13v2.0 assembly from AWS Open Data...")
    run_cmd(["curl", "-sL", "--fail", CHM13_URL, "-o", str(gz_path)], timeout=900)
    log(f"Downloaded chm13v2.0.fa.gz ({gz_path.stat().st_size / (1<<20):.1f} MiB); splitting chromosomes...")

    import gzip

    out_map: dict[str, str] = {}
    cur_fh = None
    cur_name = None
    wanted = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    with gzip.open(gz_path, "rb") as gz:
        for line in gz:
            if line.startswith(b">"):
                if cur_fh is not None:
                    cur_fh.close()
                    cur_fh = None
                hdr = line[1:].strip().split()[0].decode("ascii", errors="ignore")
                if hdr in wanted:
                    p = dest_dir / f"{hdr}.fa"
                    out_map[hdr] = str(p)
                    cur_fh = p.open("wb")
                    cur_fh.write(line)
                    cur_name = hdr
                else:
                    cur_name = None
            elif cur_fh is not None and cur_name is not None:
                cur_fh.write(line)
    if cur_fh is not None:
        cur_fh.close()
    gz_path.unlink(missing_ok=True)
    marker.write_text(json.dumps(out_map, indent=2), encoding="utf-8")
    log(f"Split {len(out_map)} T2T-CHM13v2.0 chromosomes cleanly.")
    return {k: Path(v) for k, v in out_map.items()}


def run_hgs2_compress_verify(
    tools: dict,
    target: Path,
    *,
    mode: str = "collection",
    threads: int = 4,
    xmx: str = "-Xmx24g",
    bounded: bool = False,
    panel: Path | None = None,
    timeout: int = 1200,
) -> dict:
    """Runs HGS2 compress + decompress + exact SHA-256 verification and returns structured metrics."""
    raw_bytes = target.stat().st_size
    in_sha = sha256_file(target)
    arc = WORK_ROOT / f"hgs_{os.getpid()}_{time.monotonic_ns()}.hgs2"
    restored = WORK_ROOT / f"dec_{os.getpid()}_{time.monotonic_ns()}.out"

    cmd = [
        tools["java"],
        xmx,
        "-cp",
        str(BIN_DIR),
        "hypergrass.EimdlCli",
        "compress",
        "--target",
        str(target),
        "--out",
        str(arc),
        "--threads",
        str(threads),
    ]
    if mode == "collection":
        cmd.append("--collection")
    elif mode == "single":
        cmd.append("--single")
    if bounded:
        cmd.append("--bounded")
    if panel is not None:
        cmd.extend(["--panel", str(panel)])

    try:
        enc_res = run_cmd(cmd, timeout=timeout, measure_rss=True)
        telemetry = json.loads(enc_res["stdout"].strip().splitlines()[-1])
        dec_cmd = [
            tools["java"],
            xmx,
            "-cp",
            str(BIN_DIR),
            "hypergrass.EimdlCli",
            "decompress",
            "--archive",
            str(arc),
            "--out",
            str(restored),
        ]
        if panel is not None:
            dec_cmd.extend(["--panel", str(panel)])
        dec_res = run_cmd(dec_cmd, timeout=timeout, measure_rss=True)
        dec_tel = json.loads(dec_res["stdout"].strip().splitlines()[-1])
        out_sha = sha256_file(restored)
        if out_sha != in_sha:
            raise RuntimeError(f"SHA-256 mismatch on decode of {target.name}!")

        arch_bytes = arc.stat().st_size
        return {
            "status": "ok",
            "codec": f"HGS2_{mode}",
            "xmx": xmx,
            "bounded": bounded,
            "threads": threads,
            "raw_bytes": raw_bytes,
            "target_bases": telemetry["target_bases"],
            "archive_bytes": arch_bytes,
            "bpb": (8.0 * arch_bytes) / raw_bytes,
            "ratio": raw_bytes / arch_bytes,
            "enc_wall_s": round(enc_res["elapsed"], 4),
            "enc_inproc_s": telemetry["compress_seconds"],
            "enc_peak_rss_mb": round(enc_res["peak_rss_kb"] / 1024.0, 2) if enc_res["peak_rss_kb"] else None,
            "dec_wall_s": round(dec_res["elapsed"], 4),
            "dec_inproc_s": dec_tel["decompress_seconds"],
            "dec_peak_rss_mb": round(dec_res["peak_rss_kb"] / 1024.0, 2) if dec_res["peak_rss_kb"] else None,
            "sha256_verified": True,
            "telemetry": telemetry,
        }
    except Exception as exc:
        return {
            "status": "error",
            "codec": f"HGS2_{mode}",
            "xmx": xmx,
            "bounded": bounded,
            "threads": threads,
            "raw_bytes": raw_bytes,
            "error": str(exc)[-600:],
        }
    finally:
        arc.unlink(missing_ok=True)
        restored.unlink(missing_ok=True)


def run_hgs2_native_compress_verify(
    tools: dict,
    target: Path,
    *,
    mode: str = "collection",
    threads: int = 4,
    bounded: bool = False,
    panel: Path | None = None,
    timeout: int = 1200,
) -> dict:
    """Runs Native C11/C++20/Rust eimdl compress + decompress + SHA-256 verification."""
    label = f"HGS2_Native_{mode}" + ("_bounded" if bounded else "")
    if not tools.get("eimdl_native"):
        return {"status": "unavailable", "codec": label}
    raw_bytes = target.stat().st_size
    in_sha = sha256_file(target)
    arc = WORK_ROOT / f"nhgs_{os.getpid()}_{time.monotonic_ns()}.hgs2"
    restored = WORK_ROOT / f"ndec_{os.getpid()}_{time.monotonic_ns()}.out"

    cmd = [
        tools["eimdl_native"],
        "compress",
        "--target",
        str(target),
        "--out",
        str(arc),
        "--threads",
        str(threads),
    ]
    if mode == "collection":
        cmd.append("--collection")
    elif mode == "single":
        cmd.append("--single")
    if bounded:
        cmd.append("--bounded")
    if panel is not None:
        cmd.extend(["--panel", str(panel)])

    try:
        enc_res = run_cmd(cmd, timeout=timeout, measure_rss=True)
        telemetry = json.loads(enc_res["stdout"].strip().splitlines()[-1])
        dec_cmd = [
            tools["eimdl_native"],
            "decompress",
            "--archive",
            str(arc),
            "--out",
            str(restored),
            "--threads",
            str(threads),
        ]
        if bounded:
            dec_cmd.append("--bounded")
        if panel is not None:
            dec_cmd.extend(["--panel", str(panel)])
        dec_res = run_cmd(dec_cmd, timeout=timeout, measure_rss=True)
        dec_tel = json.loads(dec_res["stdout"].strip().splitlines()[-1])
        out_sha = sha256_file(restored)
        if out_sha != in_sha:
            raise RuntimeError(f"Native SHA-256 mismatch on decode of {target.name}!")
        arch_bytes = arc.stat().st_size
        return {
            "status": "ok",
            "codec": label,
            "bounded": bounded,
            "threads": threads,
            "raw_bytes": raw_bytes,
            "target_bases": telemetry["target_bases"],
            "archive_bytes": arch_bytes,
            "bpb": (8.0 * arch_bytes) / raw_bytes,
            "ratio": raw_bytes / arch_bytes,
            "enc_wall_s": round(enc_res["elapsed"], 4),
            "enc_inproc_s": telemetry["compress_seconds"],
            "enc_peak_rss_mb": round(enc_res["peak_rss_kb"] / 1024.0, 2) if enc_res["peak_rss_kb"] else None,
            "dec_wall_s": round(dec_res["elapsed"], 4),
            "dec_inproc_s": dec_tel["decompress_seconds"],
            "dec_peak_rss_mb": round(dec_res["peak_rss_kb"] / 1024.0, 2) if dec_res["peak_rss_kb"] else None,
            "sha256_verified": True,
            "telemetry": telemetry,
        }
    except Exception as exc:
        return {
            "status": "error",
            "codec": label,
            "bounded": bounded,
            "threads": threads,
            "raw_bytes": raw_bytes,
            "error": str(exc)[-600:],
        }
    finally:
        arc.unlink(missing_ok=True)
        restored.unlink(missing_ok=True)


def run_zstd_baseline(tools: dict, target: Path, *, level: int = 19, timeout: int = 900) -> dict:
    if not tools.get("zstd"):
        return {"status": "unavailable", "codec": f"zstd-{level}"}
    raw_bytes = target.stat().st_size
    in_sha = sha256_file(target)
    arc = WORK_ROOT / f"z_{time.monotonic_ns()}.zst"
    dec = WORK_ROOT / f"zd_{time.monotonic_ns()}.out"
    try:
        enc = run_cmd([tools["zstd"], f"-{level}", "-T0", "-q", "-f", str(target), "-o", str(arc)], timeout=timeout, measure_rss=True)
        d_res = run_cmd([tools["zstd"], "-d", "-q", "-f", str(arc), "-o", str(dec)], timeout=timeout, measure_rss=True)
        assert sha256_file(dec) == in_sha
        arch = arc.stat().st_size
        return {
            "status": "ok",
            "codec": f"zstd-{level}",
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "enc_peak_rss_mb": round(enc["peak_rss_kb"] / 1024.0, 2) if enc["peak_rss_kb"] else None,
            "dec_wall_s": round(d_res["elapsed"], 4),
            "sha256_verified": True,
        }
    except Exception as exc:
        return {"status": "error", "codec": f"zstd-{level}", "error": str(exc)[-400:]}
    finally:
        arc.unlink(missing_ok=True)
        dec.unlink(missing_ok=True)


def run_xz_baseline(tools: dict, target: Path, *, level: int = 6, timeout: int = 900) -> dict:
    if not tools.get("xz"):
        return {"status": "unavailable", "codec": f"xz-{level}"}
    raw_bytes = target.stat().st_size
    in_sha = sha256_file(target)
    arc = WORK_ROOT / f"x_{time.monotonic_ns()}.xz"
    dec = WORK_ROOT / f"xd_{time.monotonic_ns()}.out"
    try:
        enc = run_cmd(f'"{tools["xz"]}" -{level} -T0 -c "{target}" > "{arc}"', shell=True, timeout=timeout)
        d_res = run_cmd(f'"{tools["xz"]}" -d -c "{arc}" > "{dec}"', shell=True, timeout=timeout)
        assert sha256_file(dec) == in_sha
        arch = arc.stat().st_size
        return {
            "status": "ok",
            "codec": f"xz-{level}",
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "dec_wall_s": round(d_res["elapsed"], 4),
            "sha256_verified": True,
        }
    except Exception as exc:
        return {"status": "error", "codec": f"xz-{level}", "error": str(exc)[-400:]}
    finally:
        arc.unlink(missing_ok=True)
        dec.unlink(missing_ok=True)


def run_jarvis3_baseline(tools: dict, seq_path: Path, *, level: int = 7, timeout: int = 900) -> dict:
    """Runs JARVIS3 on canonical flat ACGT .seq inside an isolated temp directory + SHA-256 verify."""
    if not tools.get("jarvis3"):
        return {"status": "unavailable", "codec": f"JARVIS3-l{level}"}
    raw_bytes = seq_path.stat().st_size
    in_sha = sha256_file(seq_path)
    jdir = WORK_ROOT / f"jwork_{time.monotonic_ns()}"
    jdir.mkdir(parents=True, exist_ok=True)
    local_seq = jdir / "input.seq"
    try:
        shutil.copy2(seq_path, local_seq)
        jc = jdir / f"{local_seq.name}.jc"
        enc = run_cmd(
            [tools["jarvis3"], "-v", "-l", str(level), "-o", jc.name, local_seq.name],
            cwd=jdir,
            timeout=timeout,
            measure_rss=True,
        )
        arch = jc.stat().st_size
        local_seq.unlink(missing_ok=True)
        d_res = run_cmd(
            [tools["jarvis3"], "-d", "-o", local_seq.name, jc.name],
            cwd=jdir,
            timeout=timeout,
            measure_rss=True,
        )
        restored = jdir / local_seq.name
        assert sha256_file(restored) == in_sha
        return {
            "status": "ok",
            "codec": f"JARVIS3-l{level}",
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "enc_peak_rss_mb": round(enc["peak_rss_kb"] / 1024.0, 2) if enc["peak_rss_kb"] else None,
            "dec_wall_s": round(d_res["elapsed"], 4),
            "dec_peak_rss_mb": round(d_res["peak_rss_kb"] / 1024.0, 2) if d_res["peak_rss_kb"] else None,
            "sha256_verified": True,
        }
    except Exception as exc:
        return {"status": "error", "codec": f"JARVIS3-l{level}", "error": str(exc)[-500:]}
    finally:
        shutil.rmtree(jdir, ignore_errors=True)


def _fasta_seq_sha256(path: Path) -> str:
    """Computes SHA-256 over uppercased non-header DNA bases (line-wrap invariant)."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for line in f:
            if line.startswith(b">"):
                continue
            s = line.strip()
            if s:
                h.update(s.upper())
    return h.hexdigest()


def run_agc_collection(
    tools: dict,
    ref_fa: Path,
    other_fas: list[Path],
    *,
    adaptive: bool = False,
    threads: int = 4,
    timeout: int = 900,
) -> dict:
    """Runs AGC create + full getset sequence SHA-256 verification across sample FASTAs."""
    if not tools.get("agc"):
        return {"status": "unavailable", "codec": "AGC-a" if adaptive else "AGC"}
    all_fas = [ref_fa, *other_fas]
    raw_bytes = sum(p.stat().st_size for p in all_fas)
    expected_hashes = {p.stem: _fasta_seq_sha256(p) for p in all_fas}
    awork = WORK_ROOT / f"agc_{time.monotonic_ns()}"
    awork.mkdir(parents=True, exist_ok=True)
    arc = awork / "out.agc"
    lst = awork / "others.list"
    lst.write_text("\n".join(str(p) for p in other_fas) + "\n", encoding="utf-8")
    label = "AGC-a" if adaptive else "AGC"
    cmd = [tools["agc"], "create"]
    if adaptive:
        cmd.append("-a")
    cmd.extend(["-t", str(threads), "-i", str(lst), "-o", str(arc), str(ref_fa)])
    try:
        enc = run_cmd(cmd, timeout=timeout, measure_rss=True)
        arch = arc.stat().st_size
        ls_res = run_cmd([tools["agc"], "listset", str(arc)], timeout=120)
        samples = [x.strip() for x in ls_res["stdout"].splitlines() if x.strip()]
        t_dec0 = time.monotonic()
        verified_samples = 0
        for s in samples:
            safe_s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
            restored = awork / f"{safe_s}.restored.fa"
            with restored.open("wb") as fh:
                subprocess.run(
                    [tools["agc"], "getset", str(arc), s],
                    stdout=fh,
                    check=True,
                    timeout=300,
                )
            got_hash = _fasta_seq_sha256(restored)
            restored.unlink(missing_ok=True)
            key = s.split()[0]
            if key in expected_hashes and got_hash != expected_hashes[key]:
                raise RuntimeError(f"AGC sequence SHA-256 mismatch for {key}")
            verified_samples += 1
        dec_elapsed = time.monotonic() - t_dec0
        return {
            "status": "ok",
            "codec": label,
            "samples": len(all_fas),
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "enc_peak_rss_mb": round(enc["peak_rss_kb"] / 1024.0, 2) if enc["peak_rss_kb"] else None,
            "dec_wall_s": round(dec_elapsed, 4),
            "verified_getset_samples": verified_samples,
            "sha256_verified": verified_samples == len(all_fas),
        }
    except Exception as exc:
        return {"status": "error", "codec": label, "error": str(exc)[-500:]}
    finally:
        shutil.rmtree(awork, ignore_errors=True)
