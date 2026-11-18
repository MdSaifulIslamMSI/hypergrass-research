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


def run_timed(cmd: list[str], *, cwd: Path | None = None) -> tuple[float, int | None, str]:
    t0 = time.monotonic()
    time_bin = Path("/usr/bin/time")
    rss_kb = None
    if time_bin.exists():
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            tlog = Path(tf.name)
        proc = subprocess.run(
            [str(time_bin), "-f", "%M", "-o", str(tlog), *cmd],
            cwd=str(cwd) if cwd else None,
            text=True,
            capture_output=True,
        )
        elapsed = time.monotonic() - t0
        try:
            rss_kb = int(tlog.read_text(encoding="utf-8").strip().splitlines()[-1])
        except Exception:
            rss_kb = None
        tlog.unlink(missing_ok=True)
    else:
        proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True, capture_output=True)
        elapsed = time.monotonic() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr[-1500:]}")
    return round(elapsed, 4), rss_kb, proc.stdout


def norm_chrom(chrom: str) -> str:
    c = chrom.strip()
    for ext in (".fasta", ".fna", ".fa", ".seq"):
        if c.lower().endswith(ext):
            c = c[: -len(ext)]
            break
    low = c.lower()
    if low.startswith("chr"):
        c = c[3:]
    if c.upper() in ("M", "MT"):
        return "MT"
    return c


def norm_allele(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    r = ref.upper()
    a = alt.upper()
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r = r[:-1]
        a = a[:-1]
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r = r[1:]
        a = a[1:]
        pos += 1
    return pos, r, a


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
                intervals.setdefault(norm_chrom(parts[0]), []).append((int(parts[1]), int(parts[2])))
    return intervals


def in_bed(chrom: str, pos1: int, bed: dict[str, list[tuple[int, int]]]) -> bool:
    if not bed:
        return True
    pos0 = pos1 - 1
    for lo, hi in bed.get(norm_chrom(chrom), ()):
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
            nc = norm_chrom(chrom)
            raw_pos = int(pos_s)
            if not in_bed(nc, raw_pos, bed):
                continue
            for alt in alts.split(","):
                npos, nref, nalt = norm_allele(raw_pos, ref, alt)
                calls.add((nc, npos, nref, nalt))
    return calls


def _record_result(
    codec: str,
    raw_bytes: int,
    arch_bytes: int,
    enc_s: float,
    dec_s: float | None,
    enc_rss: int | None,
    dec_rss: int | None,
    sha_ok: bool | None,
) -> dict:
    return {
        "status": "ok",
        "codec": codec,
        "archive_bytes": arch_bytes,
        "bpb": round((arch_bytes * 8.0) / raw_bytes, 6),
        "ratio": round(raw_bytes / arch_bytes, 4),
        "enc_wall_s": enc_s,
        "dec_wall_s": dec_s,
        "enc_peak_rss_mb": round(enc_rss / 1024.0, 2) if enc_rss else None,
        "dec_peak_rss_mb": round(dec_rss / 1024.0, 2) if dec_rss else None,
        "sha256_verified": sha_ok,
    }


def _run_stream_compressor(codec: str, target: Path, raw_bytes: int, exp_sha: str, enc_cmd: list[str], dec_cmd: list[str], arc: Path, dec: Path) -> dict:
    try:
        enc_s, enc_rss, _ = run_timed(enc_cmd)
        dec_s, dec_rss, _ = run_timed(dec_cmd)
        sha_ok = sha256_file(dec) == exp_sha
        return _record_result(codec, raw_bytes, arc.stat().st_size, enc_s, dec_s, enc_rss, dec_rss, sha_ok)
    except Exception as exc:
        return {"status": "error", "codec": codec, "error": str(exc)[-300:]}
    finally:
        arc.unlink(missing_ok=True)
        dec.unlink(missing_ok=True)


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
    out["runs"].append(_record_result("EIMDL_Native_collection", raw_bytes, arc.stat().st_size, enc_s, dec_s, enc_rss, dec_rss, sha_ok))

    if shutil.which("agc"):
        agc_arc = workdir / f"{hprc_fasta.stem}.agc"
        enc_a, rss_a, _ = run_timed(["agc", "create", "-t", str(threads), "-o", str(agc_arc), str(hprc_fasta)])
        out["runs"].append(_record_result("AGC", raw_bytes, agc_arc.stat().st_size, enc_a, None, rss_a, None, None))

    if shutil.which("mbgc"):
        mdir = workdir / "mbgc_hprc"
        mdir.mkdir(exist_ok=True)
        in_copy = mdir / "input.fa"
        m_arc = mdir / "out.mbgc"
        m_dec = mdir / "dec"
        shutil.copy2(hprc_fasta, in_copy)
        try:
            enc_m, rss_m, _ = run_timed(["mbgc", "c", "-m3", "-t", str(threads), "-i", str(in_copy), str(m_arc)], cwd=mdir)
            dec_m, dec_rss_m, _ = run_timed(["mbgc", "d", "-f", str(m_arc), str(m_dec)], cwd=mdir)
            restored = next((p for p in m_dec.rglob("*") if p.is_file()), None)
            m_sha = (sha256_file(restored) == exp_sha) if restored else False
            out["runs"].append(_record_result("MBGC-m3", raw_bytes, m_arc.stat().st_size, enc_m, dec_m, rss_m, dec_rss_m, m_sha))
        except Exception as exc:
            out["runs"].append({"status": "error", "codec": "MBGC-m3", "error": str(exc)[-300:]})
        finally:
            shutil.rmtree(mdir, ignore_errors=True)

    if shutil.which("zstd"):
        z_arc = workdir / f"{hprc_fasta.stem}.fa.zst"
        z_dec = workdir / f"{hprc_fasta.stem}.zst.fa"
        out["runs"].append(_run_stream_compressor(
            "zstd-19", hprc_fasta, raw_bytes, exp_sha,
            ["zstd", "-q", "-19", f"-T{threads}", "-f", "-o", str(z_arc), str(hprc_fasta)],
            ["zstd", "-q", "-d", "-f", "-o", str(z_dec), str(z_arc)],
            z_arc, z_dec,
        ))
    return out


def bench_giab(
    eimdl: str,
    giab_fq: Path,
    panel: Path,
    truth_vcf: Path | None,
    truth_bed: Path | None,
    threads: int,
    min_support: int,
    min_qual: float,
    min_af: float,
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
    vcf_s, vcf_rss, vcf_out = run_timed([
        eimdl, "variants", "--archive", str(arc), "--panel", str(panel),
        "--min-support", str(min_support), "--min-qual", str(min_qual), "--min-af", str(min_af),
        "--out", str(vcf), "--json", "--threads", str(threads),
    ])

    baselines: list[dict] = [_record_result("EIMDL_Native_HGS3_CDS", raw_bytes, arc.stat().st_size, enc_s, dec_s, enc_rss, dec_rss, sha_ok)]
    if shutil.which("spring"):
        sp_arc = workdir / f"{giab_fq.stem}.spring"
        sp_dec = workdir / f"{giab_fq.stem}.spring.fq"
        baselines.append(_run_stream_compressor(
            "SPRING", giab_fq, raw_bytes, exp_sha,
            ["spring", "-c", "-i", str(giab_fq), "-o", str(sp_arc), "-t", str(threads)],
            ["spring", "-d", "-i", str(sp_arc), "-o", str(sp_dec), "-t", str(threads)],
            sp_arc, sp_dec,
        ))
    bsc_bin = shutil.which("bsc") or shutil.which("libbsc")
    if bsc_bin:
        b_arc = workdir / f"{giab_fq.stem}.bsc"
        b_dec = workdir / f"{giab_fq.stem}.bsc.fq"
        baselines.append(_run_stream_compressor(
            "libbsc", giab_fq, raw_bytes, exp_sha,
            [bsc_bin, "e", str(giab_fq), str(b_arc), "-b64p"],
            [bsc_bin, "d", str(b_arc), str(b_dec)],
            b_arc, b_dec,
        ))
    if shutil.which("zstd"):
        z_arc = workdir / f"{giab_fq.stem}.fq.zst"
        z_dec = workdir / f"{giab_fq.stem}.zst.fq"
        baselines.append(_run_stream_compressor(
            "zstd-19", giab_fq, raw_bytes, exp_sha,
            ["zstd", "-q", "-19", f"-T{threads}", "-f", "-o", str(z_arc), str(giab_fq)],
            ["zstd", "-q", "-d", "-f", "-o", str(z_dec), str(z_arc)],
            z_arc, z_dec,
        ))
    if shutil.which("xz"):
        x_arc = workdir / f"{giab_fq.stem}.fq.xz"
        x_dec = workdir / f"{giab_fq.stem}.xz.fq"
        baselines.append(_run_stream_compressor(
            "xz-9", giab_fq, raw_bytes, exp_sha,
            ["xz", "-9", f"-T{threads}", "-k", "-f", str(giab_fq)] if False else ["sh", "-c", f"xz -9 -T{threads} -c '{giab_fq}' > '{x_arc}'"],
            ["sh", "-c", f"xz -d -c '{x_arc}' > '{x_dec}'"],
            x_arc, x_dec,
        ))

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
        "baselines": baselines,
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
        rtg_rep = None
        if shutil.which("rtg") and shutil.which("bgzip") and shutil.which("tabix") and panel is not None:
            rtg_dir = workdir / f"rtg_{giab_fq.stem}"
            try:
                rtg_dir.mkdir(parents=True, exist_ok=True)
                ref_fa = rtg_dir / "ref.fa"
                target_panel_fa = panel if panel.is_file() else next(panel.glob("*.fa*"))
                with target_panel_fa.open("r", encoding="utf-8") as fin, ref_fa.open("w", encoding="utf-8") as fout:
                    for ln in fin:
                        if ln.startswith(">"):
                            fout.write(f">{norm_chrom(ln[1:].strip().split()[0])}\n")
                        else:
                            fout.write(ln)

                def _prep_vcf_gz(src: Path, dst_gz: Path) -> None:
                    raw_hdrs, body = [], []
                    opener = gzip.open if src.suffix == ".gz" else open
                    with opener(src, "rt", encoding="utf-8") as fin:
                        for ln in fin:
                            if ln.startswith("##contig=<ID="):
                                rest = ln[len("##contig=<ID="):]
                                cid, tail = rest.split(",", 1)
                                raw_hdrs.append(f"##contig=<ID={norm_chrom(cid)},{tail}")
                            elif ln.startswith("#"):
                                raw_hdrs.append(ln)
                            elif ln.strip():
                                cols = ln.split("\t")
                                cols[0] = norm_chrom(cols[0])
                                body.append((cols[0], int(cols[1]), "\t".join(cols)))
                    body.sort(key=lambda x: (x[0], x[1]))
                    plain = dst_gz.with_suffix("")
                    plain.write_text("".join(raw_hdrs) + "".join(b[2] for b in body), encoding="utf-8")
                    run_timed(["bgzip", "-f", str(plain)])
                    run_timed(["tabix", "-f", "-p", "vcf", str(dst_gz)])

                calls_gz = rtg_dir / "calls.vcf.gz"
                truth_gz = rtg_dir / "truth.vcf.gz"
                _prep_vcf_gz(vcf, calls_gz)
                _prep_vcf_gz(truth_vcf, truth_gz)
                sdf_dir = rtg_dir / "ref.sdf"
                run_timed(["rtg", "format", "-o", str(sdf_dir), str(ref_fa)])
                vcmd = ["rtg", "vcfeval", "-b", str(truth_gz), "-c", str(calls_gz), "-t", str(sdf_dir), "-o", str(rtg_dir / "eval_out")]
                if truth_bed and truth_bed.exists():
                    bed_norm = rtg_dir / "truth.bed"
                    opener = gzip.open if truth_bed.suffix == ".gz" else open
                    with opener(truth_bed, "rt", encoding="utf-8") as fin, bed_norm.open("w", encoding="utf-8") as fout:
                        for ln in fin:
                            if ln.strip() and not ln.startswith("#"):
                                p = ln.rstrip("\n").split("\t")
                                p[0] = norm_chrom(p[0])
                                fout.write("\t".join(p[:3]) + "\n")
                    vcmd.extend(["--evaluation-regions", str(bed_norm)])
                run_timed(vcmd)
                stxt = (rtg_dir / "eval_out" / "summary.txt").read_text(encoding="utf-8")
                rtg_rep = {"summary": stxt.strip()}
                for line in stxt.splitlines():
                    parts = line.split()
                    if parts and parts[0] == "None" and len(parts) >= 8:
                        rtg_rep.update({
                            "tp_baseline": int(parts[1]),
                            "tp_call": int(parts[2]),
                            "fp": int(parts[3]),
                            "fn": int(parts[4]),
                            "precision": float(parts[5]),
                            "recall": float(parts[6]),
                            "f1": float(parts[7]),
                        })
            except Exception as exc:
                rtg_rep = {"error": str(exc)[-300:]}
            finally:
                shutil.rmtree(rtg_dir, ignore_errors=True)

        res["giab_exact_concordance"] = {
            "min_support": min_support,
            "min_qual": min_qual,
            "min_af": min_af,
            "truth_variants": len(truth_calls),
            "predicted_variants": len(pred_calls),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": round(prec, 6),
            "recall": round(rec, 6),
            "f1": round(f1, 6),
            "rtg_vcfeval": rtg_rep,
            "happy_available": bool(shutil.which("hap.py")),
            "rtg_available": bool(shutil.which("rtg")),
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
    ap.add_argument("--min-support", type=int, default=1, help="Minimum allele support AC for VCF concordance")
    ap.add_argument("--min-qual", type=float, default=0.0, help="Minimum Phred QUAL score for compressed-domain VCF calling")
    ap.add_argument("--min-af", type=float, default=0.0, help="Minimum allele fraction AF for compressed-domain VCF calling")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", type=Path, default=Path("hprc_giab_results.json"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    summary = {
        "eimdl_bin": shutil.which(args.eimdl) or str(args.eimdl),
        "tools_detected": {t: bool(shutil.which(t)) for t in ("agc", "mbgc", "spring", "zstd", "xz", "libbsc", "bsc", "hap.py", "rtg")},
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
                args.eimdl,
                args.giab_fastq,
                args.panel,
                args.giab_truth_vcf,
                args.giab_truth_bed,
                args.threads,
                args.min_support,
                args.min_qual,
                args.min_af,
                wd,
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] Wrote {args.out}")


if __name__ == "__main__":
    main()
