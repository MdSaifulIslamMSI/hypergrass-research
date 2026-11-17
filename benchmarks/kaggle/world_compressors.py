"""World-leading genomic and general-purpose compressor builders and verified runners.

Designed strictly for Kaggle CPU/GPU/TPU cloud sandboxes so massive end-to-end
shootouts and 7-suite verification run safely within Kaggle free-tier resource limits
(32 GiB RAM, bounded /tmp disk with immediate cleanup, zero local test runs).
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time
from urllib.parse import quote
from urllib.request import Request, urlopen
import urllib.request

import kaggle_common as kc

GECO3_URL = "https://codeload.github.com/cobilab/geco3/tar.gz/69642a0a913e543599f35aee372e300601e4cbbe"
GECO3_SHA256 = "b870ca010a6faede7e4ddb32881629638d895295f6e1ce9e1be23bffa321f90a"
GRASS_REPO = "https://github.com/AnshuOishik/GraSS.git"
GRASS_COMMIT = "6d93700f5a2de8841f7cf2c20064413c80df2c4a"

WORK_ROOT = kc.WORK_ROOT


def _smoke_seq_and_fasta(base_dir: Path) -> tuple[Path, Path, Path, Path]:
    """Creates deterministic 24 KiB FASTA/seq targets + homologous ref for preflight checks."""
    base_dir.mkdir(parents=True, exist_ok=True)
    seq_a = []
    # Deterministic mix of repeat + non-trivial SNPs (80-bp wrapped well-formed FASTA)
    alphabet = "ACGT"
    x = 0x9E3779B9
    for _ in range(24_000):
        x = (1664525 * x + 1013904223) & 0xFFFFFFFF
        seq_a.append(alphabet[(x >> 24) & 3])
    s1 = "".join(seq_a[:12_000])
    # Second record differs slightly (~2% SNPs) for coalescent/pangenome smoke tests
    s2_chars = list(s1)
    for i in range(0, len(s2_chars), 47):
        s2_chars[i] = alphabet[(alphabet.index(s2_chars[i]) + 1) & 3]
    s2 = "".join(s2_chars)

    def write_fa(p: Path, recs: list[tuple[str, str]]) -> None:
        lines: list[str] = []
        for hdr, body in recs:
            lines.append(f">{hdr}")
            for i in range(0, len(body), 80):
                lines.append(body[i : i + 80])
        p.write_text("\n".join(lines) + "\n", encoding="ascii")

    ref_seq = base_dir / "smoke_ref.seq"
    tgt_seq = base_dir / "smoke_tgt.seq"
    ref_fa = base_dir / "smoke_ref.fa"
    multi_fa = base_dir / "smoke_multi.fa"
    ref_seq.write_text(s1, encoding="ascii")
    tgt_seq.write_text(s2, encoding="ascii")
    write_fa(ref_fa, [("ref_chr", s1)])
    write_fa(multi_fa, [("sample_0", s1), ("sample_1", s2)])
    return ref_seq, tgt_seq, ref_fa, multi_fa


def ensure_world_tools(base_tools: dict) -> dict:
    """Installs APT archivers, compiles open-source genomic world leaders, and smoke-verifies each."""
    tools = dict(base_tools)
    kc.log("Installing extra Linux world-leading compressors (brotli, bzip2, gzip, zpaq, cmake, zstd-dev)...")
    try:
        kc.run_cmd(
            "DEBIAN_FRONTEND=noninteractive apt-get update -qq && "
            "DEBIAN_FRONTEND=noninteractive apt-get install -y -qq "
            "cmake g++ gcc make git curl time zlib1g-dev libzstd-dev "
            "zstd xz-utils brotli bzip2 gzip zpaq || true",
            shell=True,
            check=False,
            timeout=300,
        )
    except Exception as exc:
        kc.log(f"Warning: APT install skipped/partial: {exc}")

    for name in ("brotli", "bzip2", "gzip", "zpaq", "zstd", "xz"):
        tools[name] = shutil.which(name) or tools.get(name)

    cpus = max(2, min(8, os.cpu_count() or 4))

    # 1) Pinned GeCo3 & GeDe3 (SHA-256 verified)
    try:
        g_tgz = WORK_ROOT / "geco3.tar.gz"
        if not g_tgz.exists():
            urllib.request.urlretrieve(GECO3_URL, g_tgz)
        got = kc.sha256_file(g_tgz)
        if got != GECO3_SHA256:
            raise RuntimeError(f"GeCo3 SHA-256 mismatch: {got}")
        g_dir = WORK_ROOT / "geco3"
        g_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(g_tgz, "r:gz") as tf:
            tf.extractall(g_dir)
        inner_src = next(g_dir.glob("*/src"))
        kc.run_cmd(["make", "-C", str(inner_src), "CC=gcc", "CFLAGS=-Wall -O3"], timeout=300)
        gc_bin, gd_bin = inner_src / "GeCo3", inner_src / "GeDe3"
        if gc_bin.exists() and gd_bin.exists():
            tools["geco3"] = str(gc_bin)
            tools["gede3"] = str(gd_bin)
            kc.log(f"Built pinned GeCo3 at {gc_bin}")
    except Exception as exc:
        kc.log(f"Warning: GeCo3 build failed: {exc}")

    # 2) MBGC (kowallus/mbgc)
    try:
        m_dir = WORK_ROOT / "mbgc_repo"
        if not m_dir.exists():
            kc.run_cmd(
                ["git", "clone", "--depth", "1", "https://github.com/kowallus/mbgc.git", str(m_dir)],
                timeout=300,
            )
        m_bld = m_dir / "build"
        kc.run_cmd(
            ["cmake", "-S", str(m_dir), "-B", str(m_bld), "-DCMAKE_BUILD_TYPE=Release"],
            timeout=300,
        )
        kc.run_cmd(["cmake", "--build", str(m_bld), f"-j{cpus}"], timeout=600)
        for cand in [m_bld / "mbgc", *m_bld.rglob("mbgc")]:
            if cand.is_file() and cand.stat().st_mode & 0o111:
                tools["mbgc"] = str(cand)
                try:
                    rev = kc.run_cmd(["git", "-C", str(m_dir), "rev-parse", "--short", "HEAD"], check=False)
                    tools["mbgc_commit"] = rev["stdout"].strip()
                except Exception:
                    pass
                kc.log(f"Built MBGC2 ({tools.get('mbgc_commit', 'head')}) at {cand}")
                break
    except Exception as exc:
        kc.log(f"Warning: MBGC build failed: {exc}")

    # 2b) FFC (kowallus/ffc, Bioinformatics 2026)
    try:
        ffc_dir = WORK_ROOT / "ffc_repo"
        if not ffc_dir.exists():
            kc.run_cmd(
                ["git", "clone", "--depth", "1", "https://github.com/kowallus/ffc.git", str(ffc_dir)],
                timeout=300,
            )
        zstd_lib = ffc_dir / "zstd" / "lib"
        if zstd_lib.is_dir():
            kc.run_cmd(
                [
                    "make",
                    "-C",
                    str(zstd_lib),
                    "ZSTD_LEGACY_SUPPORT=0",
                    "ZSTD_LIB_DEPRECATED=0",
                    "ZSTD_LIB_DICTBUILDER=0",
                    "libzstd.a",
                ],
                timeout=300,
            )
        ffc_bld = ffc_dir / "build"
        kc.run_cmd(
            ["cmake", "-S", str(ffc_dir), "-B", str(ffc_bld), "-DCMAKE_BUILD_TYPE=Release"],
            timeout=300,
        )
        kc.run_cmd(["cmake", "--build", str(ffc_bld), f"-j{cpus}"], timeout=600)
        for cand in [ffc_bld / "ffc", *ffc_bld.rglob("ffc")]:
            if cand.is_file() and cand.stat().st_mode & 0o111:
                tools["ffc"] = str(cand)
                kc.log(f"Built FFC at {cand}")
                break
    except Exception as exc:
        kc.log(f"Warning: FFC build failed: {exc}")

    # 3) NAF (KirillKryukov/naf: ennaf & unnaf)
    try:
        n_dir = WORK_ROOT / "naf_repo"
        if not n_dir.exists():
            kc.run_cmd(
                [
                    "git",
                    "clone",
                    "--depth",
                    "1",
                    "--recurse-submodules",
                    "https://github.com/KirillKryukov/naf.git",
                    str(n_dir),
                ],
                timeout=300,
            )
        kc.run_cmd(["make", "-C", str(n_dir), f"-j{cpus}"], timeout=600)
        en = next((p for p in n_dir.rglob("ennaf") if p.is_file() and p.stat().st_mode & 0o111), None)
        un = next((p for p in n_dir.rglob("unnaf") if p.is_file() and p.stat().st_mode & 0o111), None)
        if en and un:
            tools["ennaf"] = str(en)
            tools["unnaf"] = str(un)
            kc.log(f"Built NAF at {en}")
    except Exception as exc:
        kc.log(f"Warning: NAF build failed: {exc}")

    # 4) MFCompress (cobilab/mfcompress)
    try:
        mf_dir = WORK_ROOT / "mfc_repo"
        if not mf_dir.exists():
            kc.run_cmd(
                ["git", "clone", "--depth", "1", "https://github.com/cobilab/mfcompress.git", str(mf_dir)],
                timeout=300,
            )
        try:
            kc.run_cmd(["make", "-C", str(mf_dir), "-f", "Makefile.linux", f"-j{cpus}"], timeout=300)
        except Exception:
            kc.run_cmd(["make", "-C", str(mf_dir), "-f", "Makefile.linux"], timeout=300)
        mc = next((p for p in mf_dir.rglob("MFCompressC") if p.is_file() and p.stat().st_mode & 0o111), None)
        md = next((p for p in mf_dir.rglob("MFCompressD") if p.is_file() and p.stat().st_mode & 0o111), None)
        if mc and md:
            tools["mfcompress_c"] = str(mc)
            tools["mfcompress_d"] = str(md)
            kc.log(f"Built MFCompress at {mc}")
    except Exception as exc:
        kc.log(f"Warning: MFCompress build failed: {exc}")

    # 5) libbsc (IlyaGrebnov/libbsc CPU build)
    try:
        b_dir = WORK_ROOT / "libbsc_repo"
        if not b_dir.exists():
            kc.run_cmd(
                ["git", "clone", "--depth", "1", "https://github.com/IlyaGrebnov/libbsc.git", str(b_dir)],
                timeout=300,
            )
        b_bld = b_dir / "build"
        kc.run_cmd(
            [
                "cmake",
                "-S",
                str(b_dir),
                "-B",
                str(b_bld),
                "-DCMAKE_BUILD_TYPE=Release",
                "-DBSC_ENABLE_CUDA=OFF",
            ],
            timeout=300,
        )
        kc.run_cmd(["cmake", "--build", str(b_bld), f"-j{cpus}"], timeout=300)
        for cand in [b_bld / "bsc", *b_bld.rglob("bsc")]:
            if cand.is_file() and cand.stat().st_mode & 0o111:
                tools["bsc"] = str(cand)
                kc.log(f"Built libbsc at {cand}")
                break
    except Exception as exc:
        kc.log(f"Warning: libbsc build failed: {exc}")

    # 5b) Original GraSS (AnshuOishik/GraSS pinned commit, requires ./bsc in cwd)
    try:
        gr_repo = WORK_ROOT / "grass_orig_repo"
        if not gr_repo.exists():
            kc.run_cmd(["git", "clone", GRASS_REPO, str(gr_repo)], timeout=300)
        kc.run_cmd(["git", "-C", str(gr_repo), "checkout", "--force", GRASS_COMMIT], timeout=120)
        ssg_java = next(gr_repo.rglob("SSG.java"), None)
        if ssg_java is not None:
            gr_cls = WORK_ROOT / "grass_orig_bin"
            gr_cls.mkdir(parents=True, exist_ok=True)
            gr_srcs = sorted(str(p) for p in ssg_java.parent.glob("*.java"))
            kc.run_cmd([tools["javac"], "-d", str(gr_cls), *gr_srcs], timeout=180)
            tools["grass_cp"] = str(gr_cls)
            tools["grass_commit"] = GRASS_COMMIT
            kc.log(f"Built Original GraSS ({GRASS_COMMIT[:7]}) at {gr_cls}")
    except Exception as exc:
        kc.log(f"Warning: Original GraSS build failed: {exc}")

    # 6) Preflight 24 KiB smoke tests on every tool to catch flag bugs immediately!
    smoke_dir = WORK_ROOT / f"smoke_{time.monotonic_ns()}"
    tools["smoke_checks"] = {}
    try:
        ref_seq, tgt_seq, ref_fa, multi_fa = _smoke_seq_and_fasta(smoke_dir)
        tools["smoke_checks"]["HGS2"] = kc.run_hgs2_compress_verify(tools, multi_fa, mode="collection", threads=2)["status"]
        tools["smoke_checks"]["HGS2_Native"] = kc.run_hgs2_native_compress_verify(tools, multi_fa, mode="collection", threads=2)["status"]
        tools["smoke_checks"]["HGS3_FASTA_Native"] = kc.run_hgs3_omics_verify(
            tools, multi_fa, engine="native", fmt="fasta", cds=True, threads=2, cross_check_other=True
        )["status"]
        tools["smoke_checks"]["HGS3_FASTA_Java"] = kc.run_hgs3_omics_verify(
            tools, multi_fa, engine="java", fmt="fasta", cds=True, threads=2, cross_check_other=True
        )["status"]
        fq_src = "".join(tgt_seq.read_text(encoding="ascii").split())
        smoke_fq = smoke_dir / "smoke_reads.fastq"
        with smoke_fq.open("w", encoding="ascii", newline="\n") as fh:
            rlen = 150
            for idx in range(min(160, len(fq_src) // rlen)):
                rd = fq_src[idx * rlen : (idx + 1) * rlen]
                q = ("I" * 120) + ("H" * 15) + ("F" * 10) + ("#" * 5)
                fh.write(f"@SMOKE_READ:{idx} length={rlen}\n{rd}\n+\n{q}\n")
        smoke_hgs3 = smoke_dir / "smoke_fq.hgs3"
        tools["smoke_checks"]["HGS3_FASTQ_Native"] = kc.run_hgs3_omics_verify(
            tools, smoke_fq, engine="native", fmt="fastq", panel=ref_fa, cds=True, threads=2,
            cross_check_other=True, keep_archive_path=smoke_hgs3,
        )["status"]
        tools["smoke_checks"]["HGS3_FASTQ_Java"] = kc.run_hgs3_omics_verify(
            tools, smoke_fq, engine="java", fmt="fastq", panel=ref_fa, cds=True, threads=2,
            cross_check_other=True,
        )["status"]
        smoke_query = fq_src[120:140] if len(fq_src) >= 140 else "ACGTACGTACGTACGT"
        tools["smoke_checks"]["CompressedDomain_Smoke"] = kc.run_compressed_domain_bench(
            tools, smoke_hgs3, panel=ref_fa, query_seq=smoke_query, slice_start=50, slice_end=200, threads=2
        )["status"]
        tools["smoke_checks"]["HGS1"] = run_hgs1_baseline(tools, multi_fa, threads=2)["status"]
        tools["smoke_checks"]["GraSS_seq"] = run_grass_baseline(tools, tgt_seq)["status"]
        tools["smoke_checks"]["GraSS_fa"] = run_grass_baseline(tools, multi_fa)["status"]
        tools["smoke_checks"]["AGC"] = kc.run_agc_collection(tools, ref_fa, [multi_fa], threads=2)["status"]
        tools["smoke_checks"]["JARVIS3"] = kc.run_jarvis3_baseline(tools, tgt_seq, level=7)["status"]
        tools["smoke_checks"]["MBGC"] = run_mbgc_baseline(tools, multi_fa, threads=2)["status"]

        # Probe working FFC thread flag once on tiny smoke data
        tools["_ffc_tflag"] = None
        if tools.get("ffc"):
            for tf in ["-t", "-T", ""]:
                tools["_ffc_tflag"] = tf
                if run_ffc_baseline(tools, multi_fa, threads=2)["status"] == "ok":
                    break
            else:
                tools["_ffc_tflag"] = None
        tools["smoke_checks"]["FFC"] = "ok" if tools.get("ffc") and tools["_ffc_tflag"] is not None else "unavailable"
        tools["smoke_checks"]["NAF"] = run_naf_baseline(tools, multi_fa, level=22, long_win=20)["status"]
        tools["smoke_checks"]["MFCompress"] = run_mfcompress_baseline(tools, multi_fa, threads=2)["status"]

        # Detect working GeCo3 referential CLI flags once on tiny smoke data
        tools["_geco3_ref_flags"] = None
        for rflags in [
            ["-v", "-l", "5"],
            ["-v", "-rm", "20:500:1:35:0.95/3:100:0.95", "-rm", "13:200:1:1:0.95/0:0:0", "-lr", "0.03", "-hs", "64"],
            ["-v"],
        ]:
            tools["_geco3_ref_flags"] = rflags
            st = run_geco3_baseline(tools, tgt_seq, ref_seq=ref_seq)["status"]
            if st == "ok":
                break
        else:
            tools["_geco3_ref_flags"] = None
        tools["smoke_checks"]["GeCo3_ModeS"] = run_geco3_baseline(tools, tgt_seq, level=5)["status"]
        tools["smoke_checks"]["GeCo3_ModeR"] = "ok" if tools["_geco3_ref_flags"] is not None else "error"

        # Detect working libbsc CLI flags once on tiny smoke data
        tools["_bsc_flags"] = None
        for bflags in [["-b64p"], ["-b64"], []]:
            tools["_bsc_flags"] = bflags
            if run_bsc_baseline(tools, multi_fa)["status"] == "ok":
                break
        else:
            tools["_bsc_flags"] = None
        tools["smoke_checks"]["libbsc"] = "ok" if tools["_bsc_flags"] is not None else "error"

        for cname in ("zstd", "xz", "brotli", "bzip2", "gzip", "zpaq"):
            tools["smoke_checks"][cname] = run_named_general_baseline(tools, multi_fa, cname)["status"]
        kc.log(f"Preflight 24 KiB smoke results: {json.dumps(tools['smoke_checks'])}")
    finally:
        shutil.rmtree(smoke_dir, ignore_errors=True)

    return tools


def fetch_ncbi_locked_manifest(manifest_path: Path, out_fa: Path, *, max_records: int | None = None) -> dict:
    """Downloads exact accession versions listed in a locked repo manifest with SHA-256 check."""
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    records = manifest["records"][:max_records] if max_records else manifest["records"]
    accessions = [r["accession"] for r in records]
    ua = "HyperGraSS-Kaggle-Shootout/1.0"
    eutils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def http_get(url: str) -> bytes:
        for attempt in range(4):
            try:
                req = Request(url, headers={"User-Agent": ua})
                with urlopen(req, timeout=120) as resp:
                    return resp.read()
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(3.0 * (attempt + 1))
        raise RuntimeError("unreachable")

    uids: list[str] = []
    for acc in accessions:
        q = quote(f"{acc}[accn]")
        js = json.loads(http_get(f"{eutils}/esearch.fcgi?db=nuccore&retmode=json&retmax=5&term={q}"))
        ids = js["esearchresult"]["idlist"]
        if not ids:
            raise RuntimeError(f"No NCBI UID found for {acc}")
        uids.append(ids[0])
        time.sleep(0.45)

    out_fa.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_fa.with_suffix(".fasta.tmp")
    with tmp.open("wb") as fh:
        for start in range(0, len(uids), 5):
            payload = http_get(
                f"{eutils}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={','.join(uids[start:start+5])}"
            )
            if not payload.startswith(b">"):
                raise RuntimeError(f"Unexpected NCBI efetch payload at batch {start}")
            fh.write(payload if payload.endswith(b"\n") else payload + b"\n")
            time.sleep(0.45)
    tmp.replace(out_fa)

    got_bytes = out_fa.stat().st_size
    got_sha = kc.sha256_file(out_fa)
    if max_records is None and "fasta" in manifest:
        exp = manifest["fasta"]
        if got_bytes != exp["size_bytes"] or got_sha != exp["sha256"]:
            raise RuntimeError(
                f"Manifest mismatch for {manifest_path.name}: got {got_bytes}/{got_sha}, expected {exp['size_bytes']}/{exp['sha256']}"
            )
    return {"bytes": got_bytes, "sha256": got_sha, "records": len(records)}


def split_fasta_records(fasta_path: Path, out_dir: Path, *, max_split: int = 128) -> list[Path]:
    """Splits a multi-record FASTA into per-record .fa files (or grouped bins if > max_split)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    # Count headers cleanly first
    rec_paths: list[Path] = []
    idx = -1
    fh = None
    with fasta_path.open("rb") as fin:
        for line in fin:
            if line.startswith(b">"):
                if idx + 1 < max_split:
                    if fh is not None:
                        fh.close()
                    idx += 1
                    p = out_dir / f"sample_{idx:04d}.fa"
                    rec_paths.append(p)
                    fh = p.open("wb")
            if fh is not None:
                fh.write(line)
    if fh is not None:
        fh.close()
    return rec_paths


def run_grass_baseline(
    tools: dict,
    target: Path,
    *,
    xmx: str = "-Xmx16g",
    timeout: int = 900,
) -> dict:
    """Runs Original GraSS (SSG.SSG + ./bsc -e2) compress/decompress + SHA-256 verification."""
    label = "GraSS-e2"
    if not (tools.get("grass_cp") and tools.get("bsc")):
        return {"status": "unavailable", "codec": label}
    raw_bytes = target.stat().st_size
    in_raw_sha = kc.sha256_file(target)
    in_seq_sha = kc._fasta_seq_sha256(target)
    gdir = WORK_ROOT / f"grass_{time.monotonic_ns()}"
    comp_dir = gdir / "comp"
    dec_dir = gdir / "dec"
    comp_dir.mkdir(parents=True, exist_ok=True)
    dec_dir.mkdir(parents=True, exist_ok=True)
    bsc_src = Path(tools["bsc"])
    try:
        shutil.copy2(bsc_src, comp_dir / "bsc")
        os.chmod(comp_dir / "bsc", 0o755)
        enc = kc.run_cmd(
            [tools["java"], xmx, "-cp", tools["grass_cp"], "SSG.SSG", "compress", str(target)],
            cwd=comp_dir,
            timeout=timeout,
            measure_rss=True,
        )
        arc = comp_dir / "FinalBsc.bsc"
        if not arc.exists() or arc.stat().st_size == 0:
            raise RuntimeError(
                f"GraSS produced no FinalBsc.bsc; out={(enc['stdout'] + ' ' + enc['stderr']).strip()[-250:]}"
            )
        arch = arc.stat().st_size

        shutil.copy2(arc, dec_dir / "FinalBsc.bsc")
        shutil.copy2(bsc_src, dec_dir / "bsc")
        os.chmod(dec_dir / "bsc", 0o755)
        d_res = kc.run_cmd(
            [tools["java"], xmx, "-cp", tools["grass_cp"], "SSG.SSG", "decompress"],
            cwd=dec_dir,
            timeout=timeout,
            measure_rss=True,
        )
        restored = dec_dir / "Final.SSG"
        if not restored.exists() or restored.stat().st_size == 0:
            raise RuntimeError(
                f"GraSS produced no Final.SSG; out={(d_res['stdout'] + ' ' + d_res['stderr']).strip()[-250:]}"
            )
        raw_ok = kc.sha256_file(restored) == in_raw_sha
        seq_ok = kc._fasta_seq_sha256(restored) == in_seq_sha
        if not seq_ok:
            raise RuntimeError(
                f"GraSS sequence SHA-256 mismatch (raw_ok={raw_ok}, restored_bytes={restored.stat().st_size}, input_bytes={raw_bytes})"
            )
        return {
            "status": "ok",
            "codec": label,
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "enc_peak_rss_mb": round(enc["peak_rss_kb"] / 1024.0, 2) if enc["peak_rss_kb"] else None,
            "dec_wall_s": round(d_res["elapsed"], 4),
            "dec_peak_rss_mb": round(d_res["peak_rss_kb"] / 1024.0, 2) if d_res["peak_rss_kb"] else None,
            "raw_sha256_verified": raw_ok,
            "seq_sha256_verified": seq_ok,
            "sha256_verified": seq_ok,
        }
    except Exception as exc:
        return {"status": "error", "codec": label, "error": str(exc)[-500:]}
    finally:
        shutil.rmtree(gdir, ignore_errors=True)


def run_hgs1_baseline(
    tools: dict,
    target: Path,
    *,
    single: bool = False,
    threads: int = 4,
    xmx: str = "-Xmx16g",
    timeout: int = 1200,
) -> dict:
    """Runs Generation-1 HyperGraSS v1 (HGS1) compress + decompress + exact SHA-256 verification."""
    label = "HyperGraSS_v1_single" if single else "HyperGraSS_v1"
    raw_bytes = target.stat().st_size
    in_sha = kc.sha256_file(target)
    arc = WORK_ROOT / f"hgs1_{time.monotonic_ns()}.hgs"
    dec = WORK_ROOT / f"hgs1dec_{time.monotonic_ns()}.out"
    cmd = [
        tools["java"],
        xmx,
        "-cp",
        str(kc.BIN_DIR),
        "hypergrass.HyperGraSS",
        "compress",
        "--threads",
        str(threads),
    ]
    if single:
        cmd.append("--single")
    cmd.extend([str(target), str(arc)])
    try:
        enc = kc.run_cmd(cmd, cwd=WORK_ROOT, timeout=timeout, measure_rss=True)
        d_res = kc.run_cmd(
            [tools["java"], xmx, "-cp", str(kc.BIN_DIR), "hypergrass.HyperGraSS", "decompress", str(arc), str(dec)],
            cwd=WORK_ROOT,
            timeout=timeout,
            measure_rss=True,
        )
        assert kc.sha256_file(dec) == in_sha
        arch = arc.stat().st_size
        return {
            "status": "ok",
            "codec": label,
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
        return {"status": "error", "codec": label, "error": str(exc)[-400:]}
    finally:
        arc.unlink(missing_ok=True)
        dec.unlink(missing_ok=True)


def run_geco3_baseline(
    tools: dict,
    seq_path: Path,
    *,
    level: int = 5,
    ref_seq: Path | None = None,
    timeout: int = 900,
) -> dict:
    """Runs GeCo3 / GeDe3 (Mode S or Mode R with -r) on flat .seq + exact SHA-256 verification."""
    label = "GeCo3-ref" if ref_seq is not None else f"GeCo3-l{level}"
    if not (tools.get("geco3") and tools.get("gede3")):
        return {"status": "unavailable", "codec": label}
    if ref_seq is not None and "_geco3_ref_flags" in tools and tools["_geco3_ref_flags"] is None:
        return {"status": "unavailable", "codec": label}
    raw_bytes = seq_path.stat().st_size
    in_sha = kc.sha256_file(seq_path)
    gdir = WORK_ROOT / f"geco_{time.monotonic_ns()}"
    gdir.mkdir(parents=True, exist_ok=True)
    local_tgt = gdir / "target.seq"
    local_ref = gdir / "ref.seq"
    try:
        shutil.copy2(seq_path, local_tgt)
        if ref_seq is not None:
            shutil.copy2(ref_seq, local_ref)
            flags = tools.get("_geco3_ref_flags") or ["-v", "-l", "5"]
            enc_cmd = [tools["geco3"], *flags, "-r", local_ref.name, local_tgt.name]
            dec_cmd = [tools["gede3"], "-r", local_ref.name, f"{local_tgt.name}.co"]
        else:
            enc_cmd = [tools["geco3"], "-v", "-l", str(level), local_tgt.name]
            dec_cmd = [tools["gede3"], f"{local_tgt.name}.co"]
        enc = kc.run_cmd(enc_cmd, cwd=gdir, timeout=timeout, measure_rss=True)
        co_file = gdir / f"{local_tgt.name}.co"
        if not co_file.exists():
            co_cands = sorted(gdir.glob("*.co"))
            assert co_cands, f"GeCo3 produced no .co output in {gdir}"
            co_file = co_cands[0]
        arch = co_file.stat().st_size
        local_tgt.unlink(missing_ok=True)
        if ref_seq is not None:
            dec_cmd = [tools["gede3"], "-r", local_ref.name, co_file.name]
        else:
            dec_cmd = [tools["gede3"], co_file.name]
        d_res = kc.run_cmd(dec_cmd, cwd=gdir, timeout=timeout, measure_rss=True)
        de_file = gdir / f"{local_tgt.name}.de"
        if not de_file.exists():
            # Fallback if GeDe3 stripped extension differently
            cands = sorted(gdir.glob("*.de"))
            assert cands, f"GeDe3 produced no .de output in {gdir}"
            de_file = cands[0]
        assert kc.sha256_file(de_file) == in_sha
        return {
            "status": "ok",
            "codec": label,
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
        return {"status": "error", "codec": label, "error": str(exc)[-500:]}
    finally:
        shutil.rmtree(gdir, ignore_errors=True)


def run_mbgc_baseline(
    tools: dict,
    fasta_path: Path,
    *,
    mode: int = 3,
    threads: int = 4,
    timeout: int = 900,
) -> dict:
    """Runs MBGC single-file multi-FASTA compress/decompress and verifies raw & sequence hashes."""
    label = f"MBGC-m{mode}"
    if not tools.get("mbgc"):
        return {"status": "unavailable", "codec": label}
    raw_bytes = fasta_path.stat().st_size
    in_raw_sha = kc.sha256_file(fasta_path)
    in_seq_sha = kc._fasta_seq_sha256(fasta_path)
    mdir = WORK_ROOT / f"mbgc_{time.monotonic_ns()}"
    mdir.mkdir(parents=True, exist_ok=True)
    in_copy = mdir / "input.fa"
    arc = mdir / "out.mbgc"
    out_dir = mdir / "dec"
    lossy_flag = False
    try:
        shutil.copy2(fasta_path, in_copy)
        cmd = [tools["mbgc"], "c", f"-m{mode}", "-t", str(threads), "-i", str(in_copy), str(arc)]
        try:
            enc = kc.run_cmd(cmd, cwd=mdir, timeout=timeout, measure_rss=True)
        except Exception:
            lossy_flag = True
            cmd_l = [tools["mbgc"], "c", f"-m{mode}", "-L", "-t", str(threads), "-i", str(in_copy), str(arc)]
            enc = kc.run_cmd(cmd_l, cwd=mdir, timeout=timeout, measure_rss=True)
        arch = arc.stat().st_size
        in_copy.unlink(missing_ok=True)
        d_res = kc.run_cmd([tools["mbgc"], "d", "-f", str(arc), str(out_dir)], cwd=mdir, timeout=timeout, measure_rss=True)
        restored_files = [p for p in out_dir.rglob("*") if p.is_file()]
        assert restored_files, "MBGC produced no decompressed file"
        restored = restored_files[0]
        raw_ok = kc.sha256_file(restored) == in_raw_sha
        seq_ok = kc._fasta_seq_sha256(restored) == in_seq_sha
        assert seq_ok, "MBGC sequence SHA-256 mismatch"
        return {
            "status": "ok",
            "codec": f"{label}{'-L' if lossy_flag else ''}",
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "enc_peak_rss_mb": round(enc["peak_rss_kb"] / 1024.0, 2) if enc["peak_rss_kb"] else None,
            "dec_wall_s": round(d_res["elapsed"], 4),
            "raw_sha256_verified": raw_ok,
            "seq_sha256_verified": seq_ok,
            "sha256_verified": seq_ok,
        }
    except Exception as exc:
        return {"status": "error", "codec": label, "error": str(exc)[-500:]}
    finally:
        shutil.rmtree(mdir, ignore_errors=True)


def run_ffc_baseline(
    tools: dict,
    fasta_path: Path,
    *,
    threads: int = 4,
    timeout: int = 900,
) -> dict:
    """Runs FFC (Bioinformatics 2026) FASTA compress/decompress and verifies raw & sequence SHA-256."""
    label = "FFC"
    if not tools.get("ffc"):
        return {"status": "unavailable", "codec": label}
    raw_bytes = fasta_path.stat().st_size
    in_raw_sha = kc.sha256_file(fasta_path)
    in_seq_sha = kc._fasta_seq_sha256(fasta_path)
    fdir = WORK_ROOT / f"ffc_{time.monotonic_ns()}"
    fdir.mkdir(parents=True, exist_ok=True)
    arc = fdir / "out.ffc"
    dec = fdir / "restored.fa"
    try:
        tflag = tools.get("_ffc_tflag")
        enc_cmd = [tools["ffc"], "-i", str(fasta_path), "-o", str(arc)]
        if tflag:
            enc_cmd.extend([tflag, str(threads)])
        enc = kc.run_cmd(enc_cmd, cwd=fdir, timeout=timeout, measure_rss=True)
        arch = arc.stat().st_size
        dec_cmd = [tools["ffc"], "-d", "-f", "-i", str(arc), "-o", str(dec)]
        if tflag:
            dec_cmd.extend([tflag, str(threads)])
        d_res = kc.run_cmd(dec_cmd, cwd=fdir, timeout=timeout, measure_rss=True)
        raw_ok = kc.sha256_file(dec) == in_raw_sha
        seq_ok = kc._fasta_seq_sha256(dec) == in_seq_sha
        assert raw_ok or seq_ok, "FFC sequence SHA-256 mismatch"
        return {
            "status": "ok",
            "codec": label,
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "enc_peak_rss_mb": round(enc["peak_rss_kb"] / 1024.0, 2) if enc["peak_rss_kb"] else None,
            "dec_wall_s": round(d_res["elapsed"], 4),
            "raw_sha256_verified": raw_ok,
            "seq_sha256_verified": seq_ok,
            "sha256_verified": seq_ok,
        }
    except Exception as exc:
        return {"status": "error", "codec": label, "error": str(exc)[-500:]}
    finally:
        shutil.rmtree(fdir, ignore_errors=True)


def run_naf_baseline(
    tools: dict,
    fasta_path: Path,
    *,
    level: int = 22,
    long_win: int = 31,
    timeout: int = 900,
) -> dict:
    """Runs NAF (ennaf/unnaf) on FASTA and verifies both raw and sequence SHA-256."""
    label = f"NAF-l{level}"
    if not (tools.get("ennaf") and tools.get("unnaf")):
        return {"status": "unavailable", "codec": label}
    raw_bytes = fasta_path.stat().st_size
    in_raw_sha = kc.sha256_file(fasta_path)
    in_seq_sha = kc._fasta_seq_sha256(fasta_path)
    ndir = WORK_ROOT / f"naf_{time.monotonic_ns()}"
    ndir.mkdir(parents=True, exist_ok=True)
    arc = ndir / "out.naf"
    dec = ndir / "dec.fa"
    cmd = [tools["ennaf"], "--level", str(level), "--fasta", "--temp-dir", str(ndir), "-o", str(arc), str(fasta_path)]
    if long_win > 0 and level >= 19:
        cmd.extend(["--long", str(long_win)])
    try:
        enc = kc.run_cmd(cmd, cwd=ndir, timeout=timeout, measure_rss=True)
        arch = arc.stat().st_size
        d_res = kc.run_cmd([tools["unnaf"], "--fasta", "-o", str(dec), str(arc)], cwd=ndir, timeout=timeout, measure_rss=True)
        raw_ok = kc.sha256_file(dec) == in_raw_sha
        seq_ok = kc._fasta_seq_sha256(dec) == in_seq_sha
        assert seq_ok, "NAF sequence SHA-256 mismatch"
        return {
            "status": "ok",
            "codec": label,
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "enc_peak_rss_mb": round(enc["peak_rss_kb"] / 1024.0, 2) if enc["peak_rss_kb"] else None,
            "dec_wall_s": round(d_res["elapsed"], 4),
            "raw_sha256_verified": raw_ok,
            "seq_sha256_verified": seq_ok,
            "sha256_verified": seq_ok,
        }
    except Exception as exc:
        return {"status": "error", "codec": label, "error": str(exc)[-500:]}
    finally:
        shutil.rmtree(ndir, ignore_errors=True)


def run_mfcompress_baseline(
    tools: dict,
    fasta_path: Path,
    *,
    level: int = 3,
    threads: int = 4,
    timeout: int = 900,
) -> dict:
    """Runs MFCompress (MFCompressC / MFCompressD) on FASTA + SHA-256 verification."""
    label = f"MFCompress-{level}"
    if not (tools.get("mfcompress_c") and tools.get("mfcompress_d")):
        return {"status": "unavailable", "codec": label}
    raw_bytes = fasta_path.stat().st_size
    in_raw_sha = kc.sha256_file(fasta_path)
    in_seq_sha = kc._fasta_seq_sha256(fasta_path)
    wdir = WORK_ROOT / f"mfc_{time.monotonic_ns()}"
    wdir.mkdir(parents=True, exist_ok=True)
    arc = wdir / "out.mfc"
    dec = wdir / "dec.fa"
    try:
        enc = kc.run_cmd(
            [tools["mfcompress_c"], f"-{level}", "-t", str(threads), "-o", str(arc), str(fasta_path)],
            cwd=wdir,
            timeout=timeout,
            measure_rss=True,
        )
        arch = arc.stat().st_size
        d_res = kc.run_cmd(
            [tools["mfcompress_d"], "-t", str(threads), "-o", str(dec), str(arc)],
            cwd=wdir,
            timeout=timeout,
            measure_rss=True,
        )
        raw_ok = kc.sha256_file(dec) == in_raw_sha
        seq_ok = kc._fasta_seq_sha256(dec) == in_seq_sha
        assert seq_ok, "MFCompress sequence SHA-256 mismatch"
        return {
            "status": "ok",
            "codec": label,
            "raw_bytes": raw_bytes,
            "archive_bytes": arch,
            "bpb": (8.0 * arch) / raw_bytes,
            "ratio": raw_bytes / arch,
            "enc_wall_s": round(enc["elapsed"], 4),
            "enc_peak_rss_mb": round(enc["peak_rss_kb"] / 1024.0, 2) if enc["peak_rss_kb"] else None,
            "dec_wall_s": round(d_res["elapsed"], 4),
            "raw_sha256_verified": raw_ok,
            "seq_sha256_verified": seq_ok,
            "sha256_verified": seq_ok,
        }
    except Exception as exc:
        return {"status": "error", "codec": label, "error": str(exc)[-500:]}
    finally:
        shutil.rmtree(wdir, ignore_errors=True)


def run_bsc_baseline(tools: dict, target: Path, *, timeout: int = 900) -> dict:
    """Runs libbsc (BWT + QLFC) compress + decompress + exact byte SHA-256 verify."""
    label = "libbsc"
    if not tools.get("bsc") or ("_bsc_flags" in tools and tools["_bsc_flags"] is None):
        return {"status": "unavailable", "codec": label}
    raw_bytes = target.stat().st_size
    in_sha = kc.sha256_file(target)
    arc = WORK_ROOT / f"bsc_{time.monotonic_ns()}.bsc"
    dec = WORK_ROOT / f"bscd_{time.monotonic_ns()}.out"
    flags = tools.get("_bsc_flags") or ["-b64p"]
    try:
        enc = kc.run_cmd(
            [tools["bsc"], "e", str(target), str(arc), *flags],
            cwd=WORK_ROOT,
            timeout=timeout,
            measure_rss=True,
        )
        d_res = kc.run_cmd(
            [tools["bsc"], "d", str(arc), str(dec)],
            cwd=WORK_ROOT,
            timeout=timeout,
            measure_rss=True,
        )
        assert kc.sha256_file(dec) == in_sha
        arch = arc.stat().st_size
        return {
            "status": "ok",
            "codec": label,
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
        return {"status": "error", "codec": label, "error": str(exc)[-400:]}
    finally:
        arc.unlink(missing_ok=True)
        dec.unlink(missing_ok=True)


def _run_zpaq_isolated(tools: dict, target: Path, timeout: int = 600) -> dict:
    """Runs zpaq -m5 inside an isolated cwd so path restoration is deterministic."""
    label = "zpaq-m5"
    if not tools.get("zpaq"):
        return {"status": "unavailable", "codec": label}
    raw_bytes = target.stat().st_size
    in_sha = kc.sha256_file(target)
    zdir = WORK_ROOT / f"zpaq_{time.monotonic_ns()}"
    zdir.mkdir(parents=True, exist_ok=True)
    local_in = zdir / "input.dat"
    arc = zdir / "out.zpaq"
    try:
        shutil.copy2(target, local_in)
        enc = kc.run_cmd(
            [tools["zpaq"], "a", arc.name, local_in.name, "-m5", "-threads", "4"],
            cwd=zdir,
            timeout=timeout,
            measure_rss=True,
        )
        arch = arc.stat().st_size
        local_in.unlink(missing_ok=True)
        d_res = kc.run_cmd([tools["zpaq"], "x", arc.name], cwd=zdir, timeout=timeout, measure_rss=True)
        assert kc.sha256_file(local_in) == in_sha
        return {
            "status": "ok",
            "codec": label,
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
        return {"status": "error", "codec": label, "error": str(exc)[-400:]}
    finally:
        shutil.rmtree(zdir, ignore_errors=True)


def run_named_general_baseline(
    tools: dict,
    target: Path,
    codec_name: str,
    *,
    timeout: int = 900,
) -> dict:
    """Shared verified runner for brotli-11, bzip2-9, gzip-9, zstd-ultra-22, xz-9e, and zpaq-m5."""
    if codec_name == "zpaq":
        return _run_zpaq_isolated(tools, target, timeout=timeout)
    specs = {
        "brotli": ("brotli-11", ".br", lambda a, b: [tools["brotli"], "-q", "11", "-f", "-o", str(b), str(a)], lambda a, b: [tools["brotli"], "-d", "-f", "-o", str(b), str(a)]),
        "bzip2": ("bzip2-9", ".bz2", lambda a, b: f'"{tools["bzip2"]}" -9 -c "{a}" > "{b}"', lambda a, b: f'"{tools["bzip2"]}" -d -c "{a}" > "{b}"'),
        "gzip": ("gzip-9", ".gz", lambda a, b: f'"{tools["gzip"]}" -9 -c "{a}" > "{b}"', lambda a, b: f'"{tools["gzip"]}" -d -c "{a}" > "{b}"'),
        "zstd": ("zstd-ultra-22", ".zst", lambda a, b: [tools["zstd"], "--ultra", "-22", "--long=31", "-T0", "-q", "-f", str(a), "-o", str(b)], lambda a, b: [tools["zstd"], "-d", "--long=31", "-q", "-f", str(a), "-o", str(b)]),
        "xz": ("xz-9e", ".xz", lambda a, b: f'"{tools["xz"]}" -9 -e -T0 -c "{a}" > "{b}"', lambda a, b: f'"{tools["xz"]}" -d -c "{a}" > "{b}"'),
    }
    label, ext, enc_fn, dec_fn = specs[codec_name]
    if not tools.get(codec_name):
        return {"status": "unavailable", "codec": label}
    raw_bytes = target.stat().st_size
    in_sha = kc.sha256_file(target)
    arc = WORK_ROOT / f"gen_{codec_name}_{time.monotonic_ns()}{ext}"
    dec = WORK_ROOT / f"gendec_{codec_name}_{time.monotonic_ns()}.out"
    try:
        ecmd = enc_fn(target, arc)
        dcmd = dec_fn(arc, dec)
        enc = kc.run_cmd(ecmd, shell=isinstance(ecmd, str), timeout=timeout, measure_rss=not isinstance(ecmd, str))
        d_res = kc.run_cmd(dcmd, shell=isinstance(dcmd, str), timeout=timeout, measure_rss=not isinstance(dcmd, str))
        assert kc.sha256_file(dec) == in_sha
        arch = arc.stat().st_size
        return {
            "status": "ok",
            "codec": label,
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
        return {"status": "error", "codec": label, "error": str(exc)[-400:]}
    finally:
        arc.unlink(missing_ok=True)
        dec.unlink(missing_ok=True)


def run_agc_from_multifasta(
    tools: dict,
    fasta_path: Path,
    *,
    adaptive: bool = False,
    threads: int = 4,
    timeout: int = 900,
) -> dict:
    """Runs AGC on a multi-record FASTA by splitting records into ref + others cleanly."""
    if not tools.get("agc"):
        return {"status": "unavailable", "codec": "AGC-a" if adaptive else "AGC"}
    sdir = WORK_ROOT / f"agc_split_{time.monotonic_ns()}"
    try:
        samples = split_fasta_records(fasta_path, sdir, max_split=64)
        if len(samples) < 2:
            # Single record: AGC can compress a single sample against itself or single input
            return kc.run_agc_collection(tools, samples[0], [], adaptive=adaptive, threads=threads, timeout=timeout)
        return kc.run_agc_collection(
            tools, samples[0], samples[1:], adaptive=adaptive, threads=threads, timeout=timeout
        )
    finally:
        shutil.rmtree(sdir, ignore_errors=True)


def run_agc_marginal_mod_r(
    tools: dict,
    ref_fa: Path,
    panel_other_fas: list[Path],
    target_fa: Path,
    *,
    adaptive: bool = True,
    threads: int = 4,
    timeout: int = 600,
) -> dict:
    """Computes AGC Mode R marginal size = |AGC(D_{-0} U {target})| - |AGC(D_{-0})|."""
    label = "AGC-a_marginal_ModeR" if adaptive else "AGC_marginal_ModeR"
    if not tools.get("agc"):
        return {"status": "unavailable", "codec": label}
    wdir = WORK_ROOT / f"agcmr_{time.monotonic_ns()}"
    wdir.mkdir(parents=True, exist_ok=True)
    try:
        def build_arc(others: list[Path], out_name: str) -> tuple[int, float]:
            lst = wdir / f"{out_name}.list"
            arc = wdir / f"{out_name}.agc"
            lst.write_text("\n".join(str(p) for p in others) + "\n", encoding="utf-8")
            cmd = [tools["agc"], "create"]
            if adaptive:
                cmd.append("-a")
            cmd.extend(["-t", str(threads), "-i", str(lst), "-o", str(arc), str(ref_fa)])
            res = kc.run_cmd(cmd, timeout=timeout, measure_rss=True)
            return arc.stat().st_size, res["elapsed"]

        base_bytes, _ = build_arc(panel_other_fas, "base")
        full_bytes, enc_s = build_arc([*panel_other_fas, target_fa], "with_tgt")
        marg = max(0, full_bytes - base_bytes)
        raw_b = target_fa.stat().st_size
        return {
            "status": "ok",
            "codec": label,
            "raw_bytes": raw_b,
            "base_archive_bytes": base_bytes,
            "full_archive_bytes": full_bytes,
            "archive_bytes": marg,
            "bpb": (8.0 * marg) / raw_b,
            "ratio": raw_b / max(1, marg),
            "enc_wall_s": round(enc_s, 4),
        }
    except Exception as exc:
        return {"status": "error", "codec": label, "error": str(exc)[-400:]}
    finally:
        shutil.rmtree(wdir, ignore_errors=True)
