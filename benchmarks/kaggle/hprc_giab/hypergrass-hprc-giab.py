#!/usr/bin/env python3
"""Real Public HPRC Chr21 Multi-Haplotype & NIST GIAB HG002 Benchmark on Kaggle CPU."""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile


def _locate_kaggle_scripts() -> Path:
    inp = Path("/kaggle/input")
    for p in inp.rglob("world_compressors.py"):
        return p.parent
    raise FileNotFoundError("world_compressors.py not found under /kaggle/input")


sys.path.insert(0, str(_locate_kaggle_scripts()))
import kaggle_common as kc  # noqa: E402
import world_compressors as wc  # noqa: E402

WORK = kc.WORK_ROOT / "hprc_giab_real"
OUT_JSON = kc.RESULTS_DIR / "hprc_giab_v3_2_results.json"

PAN_GRCH38 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/GRCh38/assemblies/pansn/GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz"
PAN_CHM13 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/CHM13/assemblies/analysis_set/pansn/chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz"
PAN_HG002_MAT = "https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/Q100/pansn/hg002v1.1.mat_MT.PanSN.fa.gz"
PAN_HG002_PAT = "https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/Q100/pansn/hg002v1.1.pat.PanSN.fa.gz"
HPRC_CHR21_VG = "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.chroms/chr21.vg"

GIAB_BAM = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
GIAB_VCF = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
GIAB_TBI = GIAB_VCF + ".tbi"
GIAB_BED = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
RTG_URL = "https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.13/rtg-tools-3.13-linux-x64.zip"
VG_URL = "https://github.com/vgteam/vg/releases/download/v1.78.0/vg"


def _install_extras(cpus: int) -> None:
    kc.log("Installing samtools, bcftools, tabix, SPRING, vg, rtg-tools...")
    kc.run_cmd(
        "DEBIAN_FRONTEND=noninteractive apt-get update -qq && "
        "DEBIAN_FRONTEND=noninteractive apt-get install -y -qq "
        "samtools bcftools tabix cmake g++ gcc git curl unzip libzstd-dev zlib1g-dev libbz2-dev liblzma-dev || true",
        shell=True,
        check=False,
        timeout=300,
    )
    if not shutil.which("spring"):
        try:
            sp_dir = WORK / "spring_repo"
            if not sp_dir.exists():
                kc.run_cmd(["git", "clone", "--depth", "1", "https://github.com/shubhamchandak94/Spring.git", str(sp_dir)], timeout=300)
            sp_bld = sp_dir / "build"
            kc.run_cmd(["cmake", "-S", str(sp_dir), "-B", str(sp_bld), "-DCMAKE_BUILD_TYPE=Release"], timeout=300)
            kc.run_cmd(["cmake", "--build", str(sp_bld), f"-j{cpus}"], timeout=600)
            cand = next((p for p in sp_bld.rglob("spring") if p.is_file() and (p.stat().st_mode & 0o111)), None)
            if cand:
                shutil.copy2(cand, kc.BIN_DIR / "spring")
                (kc.BIN_DIR / "spring").chmod(0o755)
                kc.log(f"Built SPRING at {kc.BIN_DIR / 'spring'}")
        except Exception as exc:
            kc.log(f"Warning: SPRING build failed: {exc}")

    if not shutil.which("rtg"):
        try:
            rtg_zip = WORK / "rtg-tools.zip"
            kc.download_file(RTG_URL, rtg_zip, timeout=300)
            rtg_ext = WORK / "rtg_ext"
            with zipfile.ZipFile(rtg_zip) as zf:
                zf.extractall(rtg_ext)
            rtg_bin = next((p for p in rtg_ext.rglob("rtg") if p.is_file()), None)
            if rtg_bin:
                rtg_bin.chmod(0o755)
                os.environ["PATH"] = f"{rtg_bin.parent}:{os.environ['PATH']}"
                kc.log(f"Installed rtg-tools at {rtg_bin}")
        except Exception as exc:
            kc.log(f"Warning: rtg-tools install failed: {exc}")

    if not shutil.which("vg"):
        try:
            vg_bin = kc.BIN_DIR / "vg"
            kc.download_file(VG_URL, vg_bin, timeout=300)
            vg_bin.chmod(0o755)
            kc.log(f"Installed vg binary at {vg_bin}")
        except Exception as exc:
            kc.log(f"Warning: vg install failed: {exc}")


def _link_world_bins(tools: dict) -> None:
    kc.BIN_DIR.mkdir(parents=True, exist_ok=True)
    for k, name in [("agc", "agc"), ("mbgc", "mbgc"), ("bsc", "bsc")]:
        src = tools.get(k)
        if src and Path(src).exists():
            dst = kc.BIN_DIR / name
            if not dst.exists():
                dst.symlink_to(Path(src))
    os.environ["PATH"] = f"{kc.BIN_DIR}:{os.environ.get('PATH', '')}"


def _prepare_grch38_chr21(panel_dir: Path) -> Path:
    panel_dir.mkdir(parents=True, exist_ok=True)
    out = panel_dir / "chr21.fa"
    for c21 in Path("/kaggle/input").rglob("*chr21*.fa*"):
        if c21.is_file() and c21.stat().st_size > 20_000_000:
            lines = []
            with c21.open("rt", encoding="utf-8", errors="ignore") as f:
                hdr = f.readline()
                lines.append(">chr21\n")
                for ln in f:
                    lines.append(ln)
            out.write_text("".join(lines), encoding="utf-8")
            kc.log(f"Loaded GRCh38 chr21 from Kaggle input: {c21} ({out.stat().st_size:,} B)")
            return out
    res = kc.run_cmd(["samtools", "faidx", PAN_GRCH38, "GRCh38#0#chr21"], timeout=300)
    txt = res["stdout"].replace(">GRCh38#0#chr21", ">chr21", 1)
    out.write_text(txt, encoding="utf-8")
    kc.log(f"Fetched GRCh38 chr21 via samtools faidx ({out.stat().st_size:,} B)")
    return out


def _prepare_hprc_chr21(grch38_c21: Path, hprc_out: Path) -> dict:
    seqs: list[tuple[str, str]] = []
    # Always pull CHM13v2.0 + HG002 Q100 Mat & Pat real chr21 haplotypes via indexed PanSN range fetch
    for url, contig, label in [
        (PAN_CHM13, "CHM13#0#chr21", "CHM13#0#chr21"),
        (PAN_HG002_MAT, "HG002#2#chr21", "HG002#2#chr21"),
        (PAN_HG002_PAT, "HG002#1#chr21", "HG002#1#chr21"),
    ]:
        try:
            r = kc.run_cmd(["samtools", "faidx", url, contig], timeout=300)
            body = "".join(ln.strip() for ln in r["stdout"].splitlines() if ln and not ln.startswith(">"))
            if len(body) > 20_000_000:
                seqs.append((label, body))
                kc.log(f"Fetched indexed haplotype {label}: {len(body):,} bp")
        except Exception as exc:
            kc.log(f"Warning: {label} fetch failed: {exc}")

    # Include GRCh38 chr21 as first cohort anchor
    ref_body = "".join(ln.strip() for ln in grch38_c21.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith(">"))
    seqs.insert(0, ("GRCh38#0#chr21", ref_body))

    # Extract real HPRC Year 1 haplotypes from chr21.vg if vg available
    vg_haps_added = 0
    if shutil.which("vg"):
        try:
            vg_path = WORK / "chr21.vg"
            kc.download_file(HPRC_CHR21_VG, vg_path, timeout=600)
            plist = kc.run_cmd(["vg", "paths", "-L", "-x", str(vg_path)], timeout=180)["stdout"].splitlines()
            seen = {lbl.split("#")[0] + "#" + (lbl.split("#")[1] if "#" in lbl else "0") for lbl, _ in seqs}
            wanted_prefixes: set[str] = set()
            for p in plist:
                p = p.strip()
                if not p or p.startswith(("GRCh38", "CHM13", "_")):
                    continue
                prefix = "#".join(p.split("#")[:2]) if "#" in p else p
                if prefix in seen:
                    continue
                seen.add(prefix)
                wanted_prefixes.add(prefix)
                if len(wanted_prefixes) >= 12:
                    break
            if wanted_prefixes:
                # Single-pass extraction from chr21.vg for all selected HPRC haplotypes
                proc = subprocess.Popen(["vg", "paths", "-F", "-x", str(vg_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                best: dict[str, tuple[str, str]] = {}
                cur_hdr, cur_buf = None, []
                assert proc.stdout is not None
                for ln in proc.stdout:
                    if ln.startswith(">"):
                        if cur_hdr:
                            pfx = "#".join(cur_hdr.split("#")[:2]) if "#" in cur_hdr else cur_hdr
                            if pfx in wanted_prefixes:
                                seq = "".join(cur_buf)
                                if len(seq) > len(best.get(pfx, ("", ""))[1]):
                                    best[pfx] = (cur_hdr, seq)
                        cur_hdr = ln[1:].strip().split()[0]
                        cur_buf = []
                    elif ln.strip():
                        cur_buf.append(ln.strip())
                if cur_hdr:
                    pfx = "#".join(cur_hdr.split("#")[:2]) if "#" in cur_hdr else cur_hdr
                    if pfx in wanted_prefixes:
                        seq = "".join(cur_buf)
                        if len(seq) > len(best.get(pfx, ("", ""))[1]):
                            best[pfx] = (cur_hdr, seq)
                proc.wait(timeout=300)
                for pfx in sorted(best):
                    hdr, seq = best[pfx]
                    if len(seq) > 20_000_000:
                        seqs.append((hdr, seq))
                        vg_haps_added += 1
                        kc.log(f"Extracted HPRC haplotype {hdr}: {len(seq):,} bp")
            vg_path.unlink(missing_ok=True)
        except Exception as exc:
            kc.log(f"Warning: chr21.vg haplotype extraction skipped/partial: {exc}")

    with hprc_out.open("w", encoding="ascii", newline="\n") as out:
        for hdr, body in seqs:
            out.write(f">{hdr}\n")
            for i in range(0, len(body), 80):
                out.write(body[i : i + 80] + "\n")

    return {"haplotypes": len(seqs), "vg_haps_added": vg_haps_added, "names": [h for h, _ in seqs], "raw_bytes": hprc_out.stat().st_size}


def _prepare_giab_hg002_chr21_region(region: str = "chr21:20000000-23000000") -> tuple[Path, Path, Path]:
    fq_out = WORK / "HG002_chr21_20M_23M.fastq"
    vcf_gz = WORK / "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    tbi_gz = WORK / "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
    bed_raw = WORK / "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

    kc.download_file(GIAB_VCF, vcf_gz, timeout=300)
    kc.download_file(GIAB_TBI, tbi_gz, timeout=120)
    kc.download_file(GIAB_BED, bed_raw, timeout=120)

    # Subset truth VCF and BED cleanly to chr21:20000000-23000000 (3 Mbp dense euchromatic locus)
    sub_vcf = WORK / "HG002_chr21_20M_23M.truth.vcf.gz"
    kc.run_cmd(["bcftools", "view", "-r", region, "-f", "PASS,.", "-O", "z", "-o", str(sub_vcf), str(vcf_gz)], timeout=180)
    kc.run_cmd(["tabix", "-f", "-p", "vcf", str(sub_vcf)], timeout=120)

    sub_bed = WORK / "HG002_chr21_20M_23M.highconf.bed"
    with bed_raw.open("rt", encoding="utf-8") as fin, sub_bed.open("w", encoding="utf-8", newline="\n") as fout:
        for ln in fin:
            if not ln or ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if p[0] == "chr21":
                lo, hi = max(20_000_000, int(p[1])), min(23_000_000, int(p[2]))
                if lo < hi:
                    fout.write(f"chr21\t{lo}\t{hi}\n")

    # Fetch real NIST/GIAB HG002 2x250bp Illumina primary reads over chr21:20000000-23000000 via HTTP byte-range BAI
    bam_sub = WORK / "HG002_chr21_20M_23M.bam"
    kc.run_cmd(["samtools", "view", "-@4", "-F", "0x904", "-b", GIAB_BAM, region, "-o", str(bam_sub)], timeout=600)
    kc.run_cmd(["samtools", "fastq", "-@4", "-N", "-o", str(fq_out), str(bam_sub)], timeout=300)
    bam_sub.unlink(missing_ok=True)
    kc.log(f"Prepared real GIAB HG002 FASTQ ({region}): {fq_out.stat().st_size:,} B")
    return fq_out, sub_vcf, sub_bed


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    cpus = max(2, min(8, os.cpu_count() or 4))

    inp = Path("/kaggle/input")
    src_root = next(inp.rglob("VerifyCompressionLaws.java")).parents[2]
    base_tools = kc.ensure_jdk_and_tools(src_root)
    tools = wc.ensure_world_tools(base_tools)
    _link_world_bins(tools)
    _install_extras(cpus)

    eimdl = tools.get("eimdl_native")
    if not eimdl:
        raise RuntimeError("Native eimdl binary not built")

    bench_script = next(inp.rglob("run_hprc_giab_bench.py"))
    smoke_script = next(inp.rglob("test_hprc_giab_smoke.py"))
    os.environ["EIMDL_BIN"] = eimdl
    kc.run_cmd([sys.executable, str(smoke_script)], cwd=bench_script.parents[2], timeout=300)
    kc.log("Offline HPRC/GIAB smoke suite passed on Kaggle VM.")

    panel_dir = WORK / "panel"
    grch38_c21 = _prepare_grch38_chr21(panel_dir)
    hprc_fa = WORK / "hprc_chr21_real.fa"
    hprc_meta = _prepare_hprc_chr21(grch38_c21, hprc_fa)

    hprc_json = WORK / "hprc_only.json"
    kc.run_cmd(
        [sys.executable, str(bench_script), "--eimdl", eimdl, "--hprc-fasta", str(hprc_fa), "--threads", str(cpus), "--out", str(hprc_json)],
        timeout=1800,
    )
    hprc_res = json.loads(hprc_json.read_text(encoding="utf-8"))
    hprc_fa.unlink(missing_ok=True)

    giab_fq, truth_vcf, truth_bed = _prepare_giab_hg002_chr21_region()
    giab_json = WORK / "giab_ms2.json"
    kc.run_cmd(
        [
            sys.executable,
            str(bench_script),
            "--eimdl",
            eimdl,
            "--panel",
            str(panel_dir),
            "--giab-fastq",
            str(giab_fq),
            "--giab-truth-vcf",
            str(truth_vcf),
            "--giab-truth-bed",
            str(truth_bed),
            "--min-support",
            "2",
            "--min-qual",
            "10.0",
            "--min-af",
            "0.10",
            "--threads",
            str(cpus),
            "--out",
            str(giab_json),
        ],
        timeout=2400,
    )
    giab_res = json.loads(giab_json.read_text(encoding="utf-8"))["giab_hg002"]

    # Also evaluate stricter (--min-support 3 --min-qual 15 --min-af 0.15) sensitivity cleanly in one extra compress+variants pass
    giab_runs = {"ms2_q10_af0.10": giab_res}
    try:
        sys.path.insert(0, str(bench_script.parent))
        import run_hprc_giab_bench as rh  # noqa: E402
        arc3 = WORK / "giab_eval3.hgs3"
        vcf3 = WORK / "giab_eval3.vcf"
        kc.run_cmd([eimdl, "compress-omics", "--target", str(giab_fq), "--panel", str(panel_dir), "--cds", "--out", str(arc3), "--threads", str(cpus)], timeout=900)
        v3_out = kc.run_cmd([
            eimdl, "variants", "--archive", str(arc3), "--panel", str(panel_dir),
            "--min-support", "3", "--min-qual", "15.0", "--min-af", "0.15",
            "--out", str(vcf3), "--json", "--threads", str(cpus),
        ], timeout=600)
        bed = rh.load_bed_intervals(truth_bed)
        pred3 = rh.parse_vcf_variants(vcf3, bed)
        truth3 = rh.parse_vcf_variants(truth_vcf, bed)
        tp3 = len(pred3 & truth3)
        fp3 = len(pred3 - truth3)
        fn3 = len(truth3 - pred3)
        p3 = tp3 / (tp3 + fp3) if (tp3 + fp3) else 0.0
        r3 = tp3 / (tp3 + fn3) if (tp3 + fn3) else 0.0
        f3 = (2.0 * p3 * r3 / (p3 + r3)) if (p3 + r3) else 0.0
        giab_runs["ms3_q15_af0.15"] = {
            "variants_report": json.loads(v3_out["stdout"]),
            "truth_variants": len(truth3),
            "predicted_variants": len(pred3),
            "tp": tp3,
            "fp": fp3,
            "fn": fn3,
            "precision": round(p3, 6),
            "recall": round(r3, 6),
            "f1": round(f3, 6),
        }
        arc3.unlink(missing_ok=True)
        vcf3.unlink(missing_ok=True)
    except Exception as exc:
        giab_runs["ms3_q15_af0.15"] = {"error": str(exc)[-300:]}

    giab_fq.unlink(missing_ok=True)

    out = {
        "schema_version": "hprc_giab_v3_2_real",
        "system": kc.collect_system_info(),
        "tools_detected": hprc_res["tools_detected"],
        "hprc_chr21_metadata": hprc_meta,
        "hprc_chr21": hprc_res["hprc"],
        "giab_hg002_chr21_region": "chr21:20000000-23000000",
        "giab_hg002_chr21_runs": giab_runs,
    }
    OUT_JSON.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    kc.log(f"DONE! Wrote real HPRC + GIAB HG002 results to {OUT_JSON}")


if __name__ == "__main__":
    main()
