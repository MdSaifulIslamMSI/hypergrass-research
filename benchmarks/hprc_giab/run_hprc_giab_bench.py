#!/usr/bin/env python3
"""Reproducible HPRC Multi-Haplotype Pangenome & GIAB HG002 Benchmark Harness for EIMDL v3."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import time


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_timed(cmd: list[str]) -> tuple[float, int | None, str]:
    t0 = time.monotonic()
    time_bin = Path("/usr/bin/time")
    rss_kb = None
    if time_bin.exists():
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            tlog = Path(tf.name)
        proc = subprocess.run([str(time_bin), "-f", "%M", "-o", str(tlog), *cmd], text=True, capture_output=True)
        elapsed = time.monotonic() - t0
        try:
            rss_kb = int(tlog.read_text(encoding="utf-8").strip().splitlines()[-1])
        except Exception:
            rss_kb = None
        tlog.unlink(missing_ok=True)
    else:
        proc = subprocess.run(cmd, text=True, capture_output=True)
        elapsed = time.monotonic() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr[-1500:]}")
    return round(elapsed, 4), rss_kb, proc.stdout


def load_bed_intervals(bed_path: Path | None) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = {}
    if bed_path is None or not bed_path.exists():
        return intervals
    opener = gzip.open if bed_path.suffix == ".gz" else open
    with opener(bed_path, "rt", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                intervals.setdefault(parts[0], []).append((int(parts[1]), int(parts[2])))
    return intervals


def in_bed(chrom: str, pos1: int, bed: dict[str, list[tuple[int, int]]]) -> bool:
    if not bed:
        return True
    pos0 = pos1 - 1
    for lo, hi in bed.get(chrom, ()):
        if lo <= pos0 < hi:
            return True
    return False


def parse_vcf_variants(vcf_path: Path, bed: dict[str, list[tuple[int, int]]]) -> set[tuple[str, int, str, str]]:
    calls: set[tuple[str, int, str, str]] = set()
    opener = gzip.open if vcf_path.suffix == ".gz" else open
    with opener(vcf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom, pos_s, _vid, ref, alts = cols[:5]
            pos = int(pos_s)
            if not in_bed(chrom, pos, bed):
                continue
            for alt in alts.split(","):
                calls.add((chrom, pos, ref.upper(), alt.upper()))
    return calls


def bench_hprc(eimdl: str, hprc_fasta: Path, threads: int, workdir: Path) -> dict:
    raw_bytes = hprc_fasta.stat().st_size
    exp_sha = sha256_file(hprc_fasta)
    out: dict = {"dataset": hprc_fasta.name, "raw_bytes": raw_bytes, "raw_sha256": exp_sha, "runs": []}

    arc = workdir / f"{hprc_fasta.stem}.hgs2"
    dec = workdir / f"{hprc_fasta.stem}.restored.fa"
    enc_s, enc_rss, _ = run_timed([eimdl, "compress", "--collection", "--target", str(hprc_fasta), "--out", str(arc), "--threads", str(threads)])
    dec_s, dec_rss, _ = run_timed([eimdl, "decompress", "--archive", str(arc), "--out", str(dec), "--threads", str(threads)])
    sha_ok = sha256_file(dec) == exp_sha
    dec.unlink(missing_ok=True)
    out["runs"].append({
        "codec": "EIMDL_Native_collection",
        "archive_bytes": arc.stat().st_size,
        "bpb": round((arc.stat().st_size * 8.0) / raw_bytes, 6),
        "ratio": round(raw_bytes / arc.stat().st_size, 4),
        "enc_wall_s": enc_s,
        "dec_wall_s": dec_s,
        "enc_peak_rss_mb": round(enc_rss / 1024.0, 2) if enc_rss else None,
        "dec_peak_rss_mb": round(dec_rss / 1024.0, 2) if dec_rss else None,
        "sha256_verified": sha_ok,
    })

    if shutil.which("agc"):
        agc_arc = workdir / f"{hprc_fasta.stem}.agc"
        enc_a, rss_a, _ = run_timed(["agc", "create", "-t", str(threads), "-o", str(agc_arc), str(hprc_fasta)])
        out["runs"].append({
            "codec": "AGC",
            "archive_bytes": agc_arc.stat().st_size,
            "bpb": round((agc_arc.stat().st_size * 8.0) / raw_bytes, 6),
            "ratio": round(raw_bytes / agc_arc.stat().st_size, 4),
            "enc_wall_s": enc_a,
            "enc_peak_rss_mb": round(rss_a / 1024.0, 2) if rss_a else None,
        })
    return out


def bench_giab(
    eimdl: str,
    giab_fq: Path,
    panel: Path,
    truth_vcf: Path | None,
    truth_bed: Path | None,
    threads: int,
    workdir: Path,
) -> dict:
    raw_bytes = giab_fq.stat().st_size
    exp_sha = sha256_file(giab_fq)
    arc = workdir / f"{giab_fq.stem}.hgs3"
    dec = workdir / f"{giab_fq.stem}.restored.fastq"
    vcf = workdir / f"{giab_fq.stem}.calls.vcf"

    enc_s, enc_rss, enc_out = run_timed([
        eimdl, "compress-omics", "--target", str(giab_fq), "--panel", str(panel), "--cds", "--out", str(arc), "--threads", str(threads)
    ])
    dec_s, dec_rss, _ = run_timed([
        eimdl, "decompress-omics", "--archive", str(arc), "--panel", str(panel), "--out", str(dec), "--threads", str(threads)
    ])
    sha_ok = sha256_file(dec) == exp_sha
    dec.unlink(missing_ok=True)

    ana_s, ana_rss, ana_out = run_timed([eimdl, "analyze", "--archive", str(arc), "--panel", str(panel), "--json", "--threads", str(threads)])
    vcf_s, vcf_rss, vcf_out = run_timed([eimdl, "variants", "--archive", str(arc), "--panel", str(panel), "--out", str(vcf), "--json"])

    res: dict = {
        "dataset": giab_fq.name,
        "raw_bytes": raw_bytes,
        "hgs3_archive_bytes": arc.stat().st_size,
        "bpb": round((arc.stat().st_size * 8.0) / raw_bytes, 6),
        "ratio": round(raw_bytes / arc.stat().st_size, 4),
        "enc_wall_s": enc_s,
        "dec_wall_s": dec_s,
        "enc_peak_rss_mb": round(enc_rss / 1024.0, 2) if enc_rss else None,
        "dec_peak_rss_mb": round(dec_rss / 1024.0, 2) if dec_rss else None,
        "sha256_verified": sha_ok,
        "encode_telemetry": json.loads(enc_out),
        "compressed_domain": {
            "analyze": {"wall_s": ana_s, "peak_rss_mb": round(ana_rss / 1024.0, 2) if ana_rss else None, "report": json.loads(ana_out)},
            "variants": {"wall_s": vcf_s, "peak_rss_mb": round(vcf_rss / 1024.0, 2) if vcf_rss else None, "report": json.loads(vcf_out)},
        },
    }

    if truth_vcf is not None and truth_vcf.exists():
        bed = load_bed_intervals(truth_bed)
        pred_calls = parse_vcf_variants(vcf, bed)
        truth_calls = parse_vcf_variants(truth_vcf, bed)
        tp = len(pred_calls & truth_calls)
        fp = len(pred_calls - truth_calls)
        fn = len(truth_calls - pred_calls)
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = (2.0 * prec * rec / (prec + rec)) if (prec + rec) else 0.0
        res["giab_exact_concordance"] = {
            "truth_variants": len(truth_calls),
            "predicted_variants": len(pred_calls),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": round(prec, 6),
            "recall": round(rec, 6),
            "f1": round(f1, 6),
        }
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description="HPRC Multi-Haplotype & GIAB HG002 EIMDL v3 Benchmark Suite")
    ap.add_argument("--eimdl", default="eimdl", help="Path to native eimdl executable")
    ap.add_argument("--hprc-fasta", type=Path, help="HPRC multi-haplotype chromosome FASTA")
    ap.add_argument("--panel", type=Path, help="Reference panel file or directory (GRCh38/CHM13)")
    ap.add_argument("--giab-fastq", type=Path, help="GIAB HG002 reads FASTQ")
    ap.add_argument("--giab-truth-vcf", type=Path, help="GIAB HG002 benchmark VCF(.gz)")
    ap.add_argument("--giab-truth-bed", type=Path, help="GIAB HG002 high-confidence BED(.gz)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", type=Path, default=Path("hprc_giab_results.json"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    summary = {
        "eimdl_bin": shutil.which(args.eimdl) or str(args.eimdl),
        "tools_detected": {t: bool(shutil.which(t)) for t in ("agc", "mbgc", "spring", "zstd", "libbsc", "bsc")},
        "hprc": None,
        "giab_hg002": None,
    }
    if args.dry_run:
        print(json.dumps(summary, indent=2))
        return

    with tempfile.TemporaryDirectory(prefix="eimdl_hprc_giab_") as td:
        wd = Path(td)
        if args.hprc_fasta:
            summary["hprc"] = bench_hprc(args.eimdl, args.hprc_fasta, args.threads, wd)
        if args.giab_fastq:
            assert args.panel, "--panel is required for --giab-fastq"
            summary["giab_hg002"] = bench_giab(
                args.eimdl, args.giab_fastq, args.panel, args.giab_truth_vcf, args.giab_truth_bed, args.threads, wd
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] Wrote {args.out}")


if __name__ == "__main__":
    main()
