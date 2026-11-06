#!/usr/bin/env python3
"""Comprehensive, reproducible, 100% SHA-256-verified benchmark runner for the
HyperGraSS + RC-MDL Fusion Engine (EIMDL / HGS2).

Every reported archive size, bits/base, compression time, decompression time,
multi-core Amdahl speedup, bounded-heap run, and ablation figure is produced by
executing the compiled Java 27 binary (`bin/`), decompressing the resulting
archive in a fresh process, and asserting exact SHA-256 equality with the raw
input bytes. Temporary archives and restored files are removed immediately after
verification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from generate_synthetic_genomes import generate_all

ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "benchmarks" / "fusion" / "results"


def default_java() -> str:
    cand = Path(os.path.expandvars(r"%USERPROFILE%\tools\jdk-27\bin\java.exe"))
    if cand.is_file():
        return str(cand)
    return "java"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_eimdl_roundtrip(
    java_exe: str,
    target: Path,
    tmp_dir: Path,
    tag: str,
    heap: str = "1800m",
    extra_compress_flags: list[str] | None = None,
    panel: Path | None = None,
    decode_heap: str | None = None,
) -> dict:
    extra = extra_compress_flags or []
    dec_heap = decode_heap or heap
    arc = tmp_dir / f"{tag}.hgs2"
    restored = tmp_dir / f"{tag}.restored"

    cmd_c = [
        java_exe,
        f"-Xmx{heap}",
        "-cp",
        str(ROOT / "bin"),
        "hypergrass.EimdlCli",
        "compress",
        "--target",
        str(target),
        "--out",
        str(arc),
    ]
    if panel is not None:
        cmd_c += ["--panel", str(panel)]
    cmd_c += extra

    t0 = time.perf_counter()
    p_c = subprocess.run(cmd_c, capture_output=True, text=True)
    wall_c = time.perf_counter() - t0
    if p_c.returncode != 0:
        raise RuntimeError(
            f"compress failed ({p_c.returncode}): {' '.join(cmd_c)}\n{p_c.stdout}\n{p_c.stderr}"
        )
    json_lines_c = [ln for ln in p_c.stdout.splitlines() if ln.strip().startswith("{")]
    stats_c = json.loads(json_lines_c[-1])

    cmd_d = [
        java_exe,
        f"-Xmx{dec_heap}",
        "-cp",
        str(ROOT / "bin"),
        "hypergrass.EimdlCli",
        "decompress",
        "--archive",
        str(arc),
        "--out",
        str(restored),
    ]
    if panel is not None:
        cmd_d += ["--panel", str(panel)]

    t1 = time.perf_counter()
    p_d = subprocess.run(cmd_d, capture_output=True, text=True)
    wall_d = time.perf_counter() - t1
    if p_d.returncode != 0:
        raise RuntimeError(
            f"decompress failed ({p_d.returncode}): {' '.join(cmd_d)}\n{p_d.stdout}\n{p_d.stderr}"
        )
    json_lines_d = [ln for ln in p_d.stdout.splitlines() if ln.strip().startswith("{")]
    stats_d = json.loads(json_lines_d[-1])

    exp_sha = sha256_file(target)
    act_sha = sha256_file(restored)
    if exp_sha != act_sha:
        raise RuntimeError(f"SHA-256 mismatch for {target}: expected {exp_sha}, got {act_sha}")

    arch_bytes = arc.stat().st_size
    raw_bytes = target.stat().st_size
    dna_bases = max(1, int(stats_c["target_bases"]))

    arc.unlink(missing_ok=True)
    restored.unlink(missing_ok=True)

    return {
        "target": target.name,
        "raw_bytes": raw_bytes,
        "target_bases": int(stats_c["target_bases"]),
        "sha256": exp_sha,
        "archive_bytes": arch_bytes,
        "bpb": arch_bytes * 8.0 / dna_bases,
        "bits_per_raw_byte": arch_bytes * 8.0 / max(1, raw_bytes),
        "ratio": raw_bytes / max(1, arch_bytes),
        "panel_genomes": int(stats_c["panel_genomes"]),
        "panel_match_bases": int(stats_c["panel_match_bases"]),
        "self_match_bases": int(stats_c["self_match_bases"]),
        "tandem_match_bases": int(stats_c["tandem_match_bases"]),
        "novelty_bases": int(stats_c["literal_bases"]),
        "panel_tokens": int(stats_c["panel_tokens"]),
        "self_tokens": int(stats_c["self_tokens"]),
        "tandem_tokens": int(stats_c["tandem_tokens"]),
        "literal_runs": int(stats_c["literal_runs"]),
        "indel_splices": int(stats_c["indel_splices"]),
        "substitutions": int(stats_c["mismatches"]),
        "compress_inproc_seconds": float(stats_c["compress_seconds"]),
        "compress_wall_seconds": round(wall_c, 4),
        "decompress_inproc_seconds": float(stats_d["decompress_seconds"]),
        "decompress_wall_seconds": round(wall_d, 4),
        "verified_sha256": True,
    }


def clusters_at(jaccard: list[list[float]], threshold: float) -> list[int]:
    n = len(jaccard)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if jaccard[i][j] >= threshold:
                a, b = find(i), find(j)
                if a != b:
                    parent[b] = a
    return [find(i) for i in range(n)]


def run_loo_panel(
    java_exe: str,
    species_key: str,
    genomes_dir: Path,
    tmp_dir: Path,
    rcmdl_results_json: Path | None = None,
    agc_results_json: Path | None = None,
) -> dict:
    genomes = sorted(genomes_dir.glob("*.seq"))
    rc_map: dict[str, int] = {}
    if rcmdl_results_json and rcmdl_results_json.is_file():
        for row in json.loads(rcmdl_results_json.read_text(encoding="utf-8")):
            if "loo" in row:
                rc_map[row["target"]] = int(row["loo"]["archive_bytes"])
            elif "archive_bytes" in row:
                rc_map[row["target"]] = int(row["archive_bytes"])

    agc_map: dict[str, int] = {}
    if agc_results_json and agc_results_json.is_file():
        for row in json.loads(agc_results_json.read_text(encoding="utf-8")):
            stem = Path(row["target"]).stem + ".seq"
            agc_map[stem] = int(row["marginal_bytes"])

    records = []
    for idx, g in enumerate(genomes):
        rec = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=g,
            tmp_dir=tmp_dir,
            tag=f"loo_{species_key}_{idx:02d}",
            heap="1600m",
            extra_compress_flags=["--exclude", g.name],
            panel=genomes_dir,
        )
        if g.name in rc_map:
            rec["rcmdl_bytes"] = rc_map[g.name]
            rec["vs_rcmdl_pct"] = round((rec["archive_bytes"] - rc_map[g.name]) * 100.0 / rc_map[g.name], 2)
        if g.name in agc_map:
            rec["agc_bytes"] = agc_map[g.name]
            rec["agc_over_eimdl"] = round(agc_map[g.name] / max(1, rec["archive_bytes"]), 2)
        records.append(rec)
        print(
            f"  [LOO {species_key}] {g.name:8s}  EIMDL={rec['archive_bytes']:8,d} B ({rec['bpb']:.5f} bpb)  "
            f"RC-MDL={rec.get('rcmdl_bytes', -1):8,d}  AGC={rec.get('agc_bytes', -1):8,d}  "
            f"enc={rec['compress_inproc_seconds']:.2f}s",
            flush=True,
        )

    tot_bases = sum(r["target_bases"] for r in records)
    tot_arch = sum(r["archive_bytes"] for r in records)
    bpbs = [r["bpb"] for r in records]
    summary = {
        "species": species_key,
        "genomes": len(records),
        "total_bases": tot_bases,
        "total_archive_bytes": tot_arch,
        "total_bpb": tot_arch * 8.0 / tot_bases,
        "median_bpb": statistics.median(bpbs),
        "total_compress_seconds": round(sum(r["compress_inproc_seconds"] for r in records), 3),
        "total_decompress_seconds": round(sum(r["decompress_inproc_seconds"] for r in records), 3),
        "all_verified_sha256": all(r["verified_sha256"] for r in records),
        "records": records,
    }
    if agc_map:
        summary["agc_wins"] = sum(1 for r in records if r.get("agc_bytes", 10**12) < r["archive_bytes"])
        summary["eimdl_wins_vs_agc"] = sum(
            1 for r in records if r["archive_bytes"] < r.get("agc_bytes", 0)
        )
    out_file = RESULTS_DIR / f"loo_{species_key}.json"
    out_file.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def run_cluster_loo_panel(
    java_exe: str,
    species_key: str,
    genomes_dir: Path,
    sketch_path: Path,
    threshold: float,
    rcmdl_cluster_baseline_bytes: int,
    tmp_dir: Path,
) -> dict:
    genomes = sorted(genomes_dir.glob("*.seq"))
    sketch = json.loads(sketch_path.read_text(encoding="utf-8"))
    roots = clusters_at(sketch["jaccard"], threshold)
    members: dict[int, list[int]] = {}
    for i, r in enumerate(roots):
        members.setdefault(r, []).append(i)

    records = []
    panel_list_file = tmp_dir / f"cluster_panel_{species_key}.txt"
    for idx, target in enumerate(genomes):
        excluded = set(members[roots[idx]])
        panel_names = [genomes[i].name for i in range(len(genomes)) if i not in excluded]
        panel_list_file.write_text("\n".join(panel_names) + "\n", encoding="utf-8")
        rec = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=target,
            tmp_dir=tmp_dir,
            tag=f"cloo_{species_key}_{idx:02d}",
            heap="1600m",
            extra_compress_flags=["--panel-list", str(panel_list_file)],
            panel=genomes_dir,
        )
        rec["excluded_cluster_size"] = len(excluded)
        records.append(rec)
        print(
            f"  [Cluster-LOO {species_key}] {target.name:8s} excl={len(excluded):2d} panel={len(panel_names):2d}  "
            f"{rec['archive_bytes']:9,d} B ({rec['bpb']:.5f} bpb)",
            flush=True,
        )
    panel_list_file.unlink(missing_ok=True)

    tot_bases = sum(r["target_bases"] for r in records)
    tot_arch = sum(r["archive_bytes"] for r in records)
    summary = {
        "species": species_key,
        "jaccard_threshold": threshold,
        "num_clusters": len(members),
        "largest_cluster": max(len(v) for v in members.values()),
        "genomes": len(records),
        "total_bases": tot_bases,
        "total_archive_bytes": tot_arch,
        "total_bpb": tot_arch * 8.0 / tot_bases,
        "rcmdl_cluster_bytes": rcmdl_cluster_baseline_bytes,
        "delta_vs_rcmdl_pct": round(
            (tot_arch - rcmdl_cluster_baseline_bytes) * 100.0 / rcmdl_cluster_baseline_bytes, 2
        ),
        "all_verified_sha256": all(r["verified_sha256"] for r in records),
        "records": records,
    }
    out_file = RESULTS_DIR / f"cluster_loo_{species_key}.json"
    out_file.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def run_mode_s_and_synthetic(java_exe: str, tmp_dir: Path) -> dict:
    syn_dir = ROOT / "benchmarks" / "fusion" / "synthetic_data"
    syn_paths = generate_all(syn_dir)

    real_targets = [
        ("ecoli_00_fa", ROOT / "prototype/rcmdl/work/ecoli/fasta/00.fa"),
        ("ecoli_00_seq", ROOT / "prototype/rcmdl/work/ecoli/genomes/00.seq"),
        ("mtb_00_seq", ROOT / "prototype/rcmdl/work/mtb/genomes/00.seq"),
        ("saureus_00_seq", ROOT / "prototype/rcmdl/work/saureus/genomes/00.seq"),
    ]
    real_results = []
    for label, p in real_targets:
        if not p.is_file():
            continue
        r = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=p,
            tmp_dir=tmp_dir,
            tag=f"modes_{label}",
            heap="1600m",
            extra_compress_flags=["--single"],
        )
        r["label"] = label
        real_results.append(r)
        print(
            f"  [Mode S Real] {label:16s}  {r['archive_bytes']:10,d} B ({r['bpb']:.4f} bpb)  "
            f"enc={r['compress_inproc_seconds']:.2f}s  dec={r['decompress_inproc_seconds']:.2f}s",
            flush=True,
        )

    # Also measure v1 HyperGraSS --single on ecoli_00_fa for direct wall-time & byte comparison on same host!
    fa00 = ROOT / "prototype/rcmdl/work/ecoli/fasta/00.fa"
    v1_single = None
    if fa00.is_file():
        v1_arc = tmp_dir / "v1_00.hgs"
        v1_out = tmp_dir / "v1_00.out"
        t0 = time.perf_counter()
        subprocess.run(
            [java_exe, "-Xmx1600m", "-cp", str(ROOT / "bin"), "hypergrass.HyperGraSS", "compress", "--single", "--no-parallel", str(fa00), str(v1_arc)],
            check=True,
            capture_output=True,
        )
        v1_c_sec = time.perf_counter() - t0
        t1 = time.perf_counter()
        subprocess.run(
            [java_exe, "-Xmx1600m", "-cp", str(ROOT / "bin"), "hypergrass.HyperGraSS", "decompress", str(v1_arc), str(v1_out)],
            check=True,
            capture_output=True,
        )
        v1_d_sec = time.perf_counter() - t1
        assert sha256_file(fa00) == sha256_file(v1_out)
        v1_single = {
            "archive_bytes": v1_arc.stat().st_size,
            "compress_wall_seconds": round(v1_c_sec, 4),
            "decompress_wall_seconds": round(v1_d_sec, 4),
            "verified_sha256": True,
        }
        v1_arc.unlink(missing_ok=True)
        v1_out.unlink(missing_ok=True)
        print(
            f"  [Baseline HyperGraSS --single --no-parallel] {v1_single['archive_bytes']:,d} B  "
            f"enc={v1_single['compress_wall_seconds']:.2f}s",
            flush=True,
        )

    syn_results = []
    for key in [
        "01_euk_centromere",
        "02_viral_rna_acgu",
        "03_chimeric_iupac_crlf",
        "04_xeno_at88",
        "05_xeno_gc86",
        "06_denovo_grammar",
        "07_hypermutator_pangenome",
        "08_crypto_white_noise",
    ]:
        p = syn_paths[key]
        mode_flag = "--collection" if key in ("02_viral_rna_acgu", "07_hypermutator_pangenome") else "--single"
        r = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=p,
            tmp_dir=tmp_dir,
            tag=f"syn_{key}",
            heap="1600m",
            extra_compress_flags=[mode_flag],
        )
        r["corpus_id"] = key
        r["regime"] = "Mode C" if mode_flag == "--collection" else "Mode S"
        # Compare against gzip -9 size on exact same file
        import gzip

        gz_bytes = len(gzip.compress(p.read_bytes(), compresslevel=9))
        r["gzip9_bytes"] = gz_bytes
        syn_results.append(r)
        print(
            f"  [Non-Existing {r['regime']}] {key:26s}  EIMDL={r['archive_bytes']:9,d} B ({r['bpb']:.4f} bpb)  "
            f"gzip-9={gz_bytes:9,d} B",
            flush=True,
        )

    # Also test Mode R on 07_hypermutator_panel (strain_05 vs strains 00..04)
    hm_dir = syn_paths["07_hypermutator_panel_dir"]
    hm_r = run_eimdl_roundtrip(
        java_exe=java_exe,
        target=hm_dir / "strain_05.seq",
        tmp_dir=tmp_dir,
        tag="syn_07_mode_r",
        heap="1600m",
        extra_compress_flags=["--exclude", "strain_05.seq"],
        panel=hm_dir,
    )
    hm_r["corpus_id"] = "07_hypermutator_strain05_loo"
    hm_r["regime"] = "Mode R"
    syn_results.append(hm_r)
    print(
        f"  [Non-Existing Mode R] 07_hypermutator_strain05    EIMDL={hm_r['archive_bytes']:9,d} B ({hm_r['bpb']:.4f} bpb)",
        flush=True,
    )

    out = {
        "real_mode_s": real_results,
        "v1_hypergrass_single_host_control": v1_single,
        "synthetic_and_non_existing": syn_results,
    }
    (RESULTS_DIR / "mode_s_and_synthetic.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    shutil.rmtree(syn_dir, ignore_errors=True)
    return out


def run_mode_c_collections(java_exe: str, tmp_dir: Path) -> dict:
    corpora = [
        ("ecoli_24_flat_seq", ROOT / "prototype/rcmdl/work/ecoli/input.seq"),
        ("ecoli_24_multifasta", ROOT / "prototype/rcmdl/work/hybrid-ecoli-00/pangenome24_blank.fasta"),
        ("mtb_24_multifasta", ROOT / "prototype/rcmdl/work/mtb/panel.fasta"),
        ("saureus_25_multifasta", ROOT / "prototype/rcmdl/work/saureus/panel.fasta"),
    ]
    rows = []
    for label, p in corpora:
        if not p.is_file():
            continue
        r = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=p,
            tmp_dir=tmp_dir,
            tag=f"modec_{label}",
            heap="2200m",
            extra_compress_flags=["--collection"],
        )
        r["label"] = label
        rows.append(r)
        print(
            f"  [Mode C] {label:24s}  {r['archive_bytes']:10,d} B ({r['bpb']:.4f} bpb)  "
            f"enc={r['compress_inproc_seconds']:.2f}s  dec={r['decompress_inproc_seconds']:.2f}s",
            flush=True,
        )
    out = {"collections": rows}
    (RESULTS_DIR / "mode_c_collections.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    return out


def run_amdahl_and_bounded_memory(java_exe: str, tmp_dir: Path) -> dict:
    ecoli_seq = ROOT / "prototype/rcmdl/work/ecoli/input.seq"
    scaling = []
    for t in [1, 2, 4, 8, 12]:
        r = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=ecoli_seq,
            tmp_dir=tmp_dir,
            tag=f"scale_t{t}",
            heap="2200m",
            extra_compress_flags=["--collection", "--threads", str(t)],
        )
        r["threads"] = t
        scaling.append(r)
        print(
            f"  [Amdahl Scaling T={t:2d}] archive={r['archive_bytes']:,d} B  "
            f"enc={r['compress_inproc_seconds']:.2f}s (wall {r['compress_wall_seconds']:.2f}s)  "
            f"dec={r['decompress_inproc_seconds']:.2f}s",
            flush=True,
        )

    t1_sec = scaling[0]["compress_inproc_seconds"]
    for r in scaling:
        r["speedup_vs_1t"] = round(t1_sec / max(1e-6, r["compress_inproc_seconds"]), 2)
        r["speedup_vs_v1_887s"] = round(887.3 / max(1e-6, r["compress_inproc_seconds"]), 2)

    # Strict End-to-End Bounded Memory Runs (-Xmx256m, -Xmx224m, and -Xmx192m on 120.1 MB ecoli/input.seq!)
    bounded_runs = []
    for heap_cap in ["256m", "224m", "192m"]:
        r = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=ecoli_seq,
            tmp_dir=tmp_dir,
            tag=f"bounded_{heap_cap}",
            heap=heap_cap,
            extra_compress_flags=["--collection", "--bounded"],
            decode_heap=heap_cap,
        )
        r["heap_cap"] = f"-Xmx{heap_cap}"
        bounded_runs.append(r)
        print(
            f"  [Bounded Memory -Xmx{heap_cap}] archive={r['archive_bytes']:,d} B ({r['bpb']:.4f} bpb)  "
            f"enc={r['compress_inproc_seconds']:.2f}s  dec={r['decompress_inproc_seconds']:.2f}s  SHA=True",
            flush=True,
        )

    out = {"amdahl_scaling": scaling, "bounded_memory": bounded_runs}
    (RESULTS_DIR / "amdahl_and_memory.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    return out


def run_real_genome_ablations(java_exe: str, tmp_dir: Path) -> dict:
    # Controlled Law I-IV ablations on real E. coli 00.seq LOO + M. tuberculosis 00.seq LOO
    cases = [
        ("full_eimdl", []),
        ("no_kimura_v4_law_I", ["--no-v4"]),
        ("no_codon_wobble_law_II", ["--no-wobble"]),
        ("no_indel_splice_law_III", ["--no-splice"]),
        ("no_v4_no_wobble_no_splice", ["--no-v4", "--no-wobble", "--no-splice"]),
        ("exact_matches_only_no_approx", ["--no-approx", "--no-splice"]),
    ]
    ecoli_dir = ROOT / "prototype/rcmdl/work/ecoli/genomes"
    mtb_dir = ROOT / "prototype/rcmdl/work/mtb/genomes"

    ecoli_rows = []
    for name, flags in cases:
        r = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=ecoli_dir / "00.seq",
            tmp_dir=tmp_dir,
            tag=f"abl_ec_{name}",
            heap="1600m",
            extra_compress_flags=["--exclude", "00.seq"] + flags,
            panel=ecoli_dir,
        )
        r["ablation"] = name
        ecoli_rows.append(r)
        print(f"  [Ablation E.coli 00] {name:28s} -> {r['archive_bytes']:8,d} B ({r['bpb']:.5f} bpb)", flush=True)

    mtb_rows = []
    for name, flags in cases:
        r = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=mtb_dir / "00.seq",
            tmp_dir=tmp_dir,
            tag=f"abl_mtb_{name}",
            heap="1600m",
            extra_compress_flags=["--exclude", "00.seq"] + flags,
            panel=mtb_dir,
        )
        r["ablation"] = name
        mtb_rows.append(r)
        print(f"  [Ablation MTB 00]    {name:28s} -> {r['archive_bytes']:8,d} B ({r['bpb']:.5f} bpb)", flush=True)

    # Diverged single-donor control (E. coli 00.seq vs distant lineage 20.seq)
    single_donor_list = tmp_dir / "single_donor_20.txt"
    single_donor_list.write_text("20.seq\n", encoding="utf-8")
    diverged_rows = []
    for name, flags in cases:
        r = run_eimdl_roundtrip(
            java_exe=java_exe,
            target=ecoli_dir / "00.seq",
            tmp_dir=tmp_dir,
            tag=f"abl_div_{name}",
            heap="1600m",
            extra_compress_flags=["--panel-list", str(single_donor_list)] + flags,
            panel=ecoli_dir,
        )
        r["ablation"] = name
        diverged_rows.append(r)
        print(f"  [Ablation Div 00v20] {name:28s} -> {r['archive_bytes']:8,d} B ({r['bpb']:.5f} bpb)", flush=True)

    out = {
        "ecoli_00_loo_ablations": ecoli_rows,
        "mtb_00_loo_ablations": mtb_rows,
        "ecoli_00_vs_distant_donor20_ablations": diverged_rows,
    }
    (RESULTS_DIR / "ablations_real_genomes.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--suite",
        default="all",
        choices=[
            "all",
            "mode_s",
            "mode_c",
            "scaling_mem",
            "loo_ecoli",
            "loo_mtb",
            "loo_saureus",
            "cluster_loo",
            "ablations",
        ],
    )
    parser.add_argument("--java", default=default_java())
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="hgs2_bench_") as td:
        tmp_dir = Path(td)
        s = args.suite
        if s in ("all", "mode_s"):
            run_mode_s_and_synthetic(args.java, tmp_dir)
        if s in ("all", "ablations"):
            run_real_genome_ablations(args.java, tmp_dir)
        if s in ("all", "loo_mtb"):
            run_loo_panel(
                args.java,
                "mtb",
                ROOT / "prototype/rcmdl/work/mtb/genomes",
                tmp_dir,
                ROOT / "prototype/rcmdl/work/runs/loo24-mtb/results.json",
                None,
            )
        if s in ("all", "loo_saureus"):
            run_loo_panel(
                args.java,
                "saureus",
                ROOT / "prototype/rcmdl/work/saureus/genomes",
                tmp_dir,
                ROOT / "prototype/rcmdl/work/runs/loo25-saureus/results.json",
                None,
            )
        if s in ("all", "loo_ecoli"):
            run_loo_panel(
                args.java,
                "ecoli",
                ROOT / "prototype/rcmdl/work/ecoli/genomes",
                tmp_dir,
                ROOT / "prototype/rcmdl/work/runs/loo24-masked/results.json",
                ROOT / "prototype/rcmdl/work/agc-loo-masked/results.json",
            )
        if s in ("all", "cluster_loo"):
            run_cluster_loo_panel(
                args.java,
                "ecoli",
                ROOT / "prototype/rcmdl/work/ecoli/genomes",
                ROOT / "prototype/rcmdl/work/orders/sketch.json",
                0.98,
                2_331_696,
                tmp_dir,
            )
            run_cluster_loo_panel(
                args.java,
                "mtb",
                ROOT / "prototype/rcmdl/work/mtb/genomes",
                ROOT / "prototype/rcmdl/work/orders-mtb/sketch.json",
                0.99,
                153_182,
                tmp_dir,
            )
            run_cluster_loo_panel(
                args.java,
                "saureus",
                ROOT / "prototype/rcmdl/work/saureus/genomes",
                ROOT / "prototype/rcmdl/work/orders-saureus/sketch.json",
                0.98,
                741_650,
                tmp_dir,
            )
        if s in ("all", "mode_c"):
            run_mode_c_collections(args.java, tmp_dir)
        if s in ("all", "scaling_mem"):
            run_amdahl_and_bounded_memory(args.java, tmp_dir)


if __name__ == "__main__":
    main()
