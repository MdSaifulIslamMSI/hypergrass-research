"""Kaggle CPU-Only End-to-End World-Leading Compressor Shootout for HGS2 / EIMDL.

Executes with zero GPU/TPU quota (Kaggle free CPU VM: 4 vCPUs, 32 GiB RAM):
  0. Compiles all Java sources (--release 25 -Xlint:all) + runs ALL 7 Verify* suites
     + builds Native C23/C++23/Rust 2024 `eimdl` + compiles & smoke-verifies world compressors (incl. FFC, MBGC2, GraSS -e2).
  1. Brutal Adversarial Torture & Law IV/VII Bounded-Heap Stress (100M@-Xmx128m .. 1Gbp@-Xmx512m).
  2. Real NCBI SHA-256-Locked Bacterial Pangenome Shootout (E. coli 24g, S. enterica 25g, K. pneumoniae 23g)
     across Mode S, Mode R (vs GeCo3 -r & AGC marginal), and Mode C (vs AGC, MBGC, NAF, MFCompress, JARVIS3, GeCo3, libbsc, zpaq, xz, zstd, brotli, bzip2, gzip).
  3. SARS-CoV-2 Viral Pangenome Scaling (1.7 MiB -> 60 MiB -> 309 MiB -> ~1.4 GiB).
  4. Human Genomes (GRCh38 Chr21/Chr22/Chr1) & Cross-Human (T2T-CHM13v2.0 | GRCh38 Mode R & Diploid Mode C).

Every run verifies SHA-256 lossless integrity, captures /usr/bin/time -v Peak RSS,
enforces soft wall-clock and /tmp storage budgets, and saves JSON progressively.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import random
import shutil
import sys
import time


def _locate_kaggle_scripts() -> Path:
    inp = Path("/kaggle/input")
    for p in inp.rglob("world_compressors.py"):
        return p.parent
    raise FileNotFoundError("world_compressors.py not found under /kaggle/input")


sys.path.insert(0, str(_locate_kaggle_scripts()))
import kaggle_common as kc  # noqa: E402
import world_compressors as wc  # noqa: E402

OUT_JSON = kc.RESULTS_DIR / "world_shootout_results.json"
SESSION_START = time.monotonic()
SOFT_BUDGET_S = 8.0 * 3600.0  # Safe margin below Kaggle 9-12h commit cutoff


def budget_remaining_s() -> float:
    return SOFT_BUDGET_S - (time.monotonic() - SESSION_START)


def save_state(state: dict) -> None:
    state["elapsed_total_s"] = round(time.monotonic() - SESSION_START, 2)
    OUT_JSON.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def _extract_first_record_fasta(multi_fa: Path, first_fa: Path, rest_fa: Path) -> None:
    """Splits record 0 into first_fa and records 1..N-1 into rest_fa for Mode R panel splits."""
    hdr_count = 0
    with multi_fa.open("rb") as fin, first_fa.open("wb") as f0, rest_fa.open("wb") as frest:
        for line in fin:
            if line.startswith(b">"):
                hdr_count += 1
            if hdr_count <= 1:
                f0.write(line)
            else:
                frest.write(line)


def run_adversarial_and_bounded_suite(tools: dict, state: dict, cov_dir: Path | None, cpus: int) -> None:
    """Runs brutal pathology inputs (random white noise + chimeric IUPAC/CRLF/tandem + 1Gbp@-Xmx512m)."""
    adv = state["adversarial_and_bounded"]

    # 1) Law IV White Noise (10 MB canonical ACGT deterministic Mersenne Twister stream)
    wn_seq = kc.WORK_ROOT / "white_noise_10m.seq"
    rng = random.Random(0x12345678)
    lut4 = bytes(b"ACGT"[b & 3] for b in range(256))
    with wn_seq.open("wb") as f:
        for _ in range(10):
            f.write(rng.randbytes(1_000_000).translate(lut4))
    kc.log("=== Brutal Torture 1A: 10 MB Crypto-Random DNA (Law IV Ceiling <= 2.0012 bpb) ===")
    wn_runs = [
        kc.run_hgs2_compress_verify(tools, wn_seq, mode="single", threads=cpus, xmx="-Xmx4g"),
        kc.run_hgs2_compress_verify(tools, wn_seq, mode="collection", threads=cpus, xmx="-Xmx4g"),
        kc.run_hgs2_native_compress_verify(tools, wn_seq, mode="single", threads=cpus),
        kc.run_hgs2_native_compress_verify(tools, wn_seq, mode="collection", threads=cpus),
        wc.run_hgs1_baseline(tools, wn_seq, single=True, threads=cpus, xmx="-Xmx4g"),
        wc.run_grass_baseline(tools, wn_seq, xmx="-Xmx4g", timeout=300),
        kc.run_jarvis3_baseline(tools, wn_seq, level=7, timeout=300),
        wc.run_geco3_baseline(tools, wn_seq, level=5, timeout=300),
        wc.run_bsc_baseline(tools, wn_seq, timeout=300),
        wc.run_named_general_baseline(tools, wn_seq, "zstd", timeout=300),
        wc.run_named_general_baseline(tools, wn_seq, "xz", timeout=300),
        wc.run_named_general_baseline(tools, wn_seq, "brotli", timeout=300),
        wc.run_named_general_baseline(tools, wn_seq, "gzip", timeout=300),
    ]
    adv["white_noise_10mb"] = wn_runs
    wn_seq.unlink(missing_ok=True)
    save_state(state)

    # 2) Brutal Pathology FASTA (~18 MB: CRLF + ragged widths + lowercase mask + RNA U + IUPAC + 200K Ns + CAG tandem + empty recs)
    path_fa = kc.WORK_ROOT / "brutal_pathology.fa"
    with path_fa.open("wb") as f:
        f.write(b">empty_leading_header\r\n")
        f.write(b">homopolymer_A_2M\r\n" + (b"A" * 80 + b"\r\n") * 25_000)
        f.write(b">trinucleotide_CAG_expansion_3M\n" + (b"CAG" * 20 + b"\n") * 50_000)
        f.write(b">ambiguous_iupac_rna_softmask_ragged\r\n")
        f.write(b"acgtACGTuuuuUUUUttttTTTTRYSWKMBDHVryswkmbdhv\r\n")
        f.write(b"A\n" + b"CCGGTTAA\r\n" + b"N" * 200_000 + b"\n")
        snp_base = bytearray((b"ACGT" * 250) * 1000)
        for rep in range(10):
            rec = bytearray(snp_base)
            for p in range(0, len(rec), 37):
                rec[p] = b"CGTA"[rec[p] & 3]
            f.write(f">chimeric_strain_{rep}\n".encode("ascii"))
            for i in range(0, len(rec), 73):
                f.write(bytes(rec[i : i + 73]) + b"\r\n")
        f.write(b">trailing_empty_no_nl")

    kc.log(f"=== Brutal Torture 1B: Pathology FASTA ({path_fa.stat().st_size / (1<<20):.1f} MiB) ===")
    adv["pathology_fasta"] = [
        kc.run_hgs2_compress_verify(tools, path_fa, mode="collection", threads=cpus, xmx="-Xmx4g"),
        kc.run_hgs2_compress_verify(tools, path_fa, mode="collection", threads=2, xmx="-Xmx128m", bounded=True),
        kc.run_hgs2_native_compress_verify(tools, path_fa, mode="collection", threads=cpus),
        kc.run_hgs2_native_compress_verify(tools, path_fa, mode="collection", threads=2, bounded=True),
        wc.run_hgs1_baseline(tools, path_fa, threads=cpus, xmx="-Xmx4g"),
        wc.run_grass_baseline(tools, path_fa, xmx="-Xmx4g", timeout=300),
        wc.run_mbgc_baseline(tools, path_fa, mode=3, threads=cpus),
        wc.run_ffc_baseline(tools, path_fa, threads=cpus),
        wc.run_naf_baseline(tools, path_fa, level=22, long_win=31),
        wc.run_mfcompress_baseline(tools, path_fa, level=3, threads=cpus),
        wc.run_bsc_baseline(tools, path_fa),
        wc.run_named_general_baseline(tools, path_fa, "zstd"),
        wc.run_named_general_baseline(tools, path_fa, "xz"),
        wc.run_named_general_baseline(tools, path_fa, "brotli"),
        wc.run_named_general_baseline(tools, path_fa, "bzip2"),
        wc.run_named_general_baseline(tools, path_fa, "gzip"),
    ]
    path_fa.unlink(missing_ok=True)
    save_state(state)

    # 3) Option 2 Cloud Bounded-Heap Law VII Stress Matrix (100M -> 1Gbp on real SARS-CoV-2 .seq)
    if cov_dir is not None:
        pool_fa = kc.WORK_ROOT / "cov_pool.fa"
        with pool_fa.open("wb") as fout:
            for sh in ["group_36.fasta", "group_34.fasta", "group_22.fasta"]:
                sp = cov_dir / sh
                if sp.exists():
                    with sp.open("rb") as fin:
                        shutil.copyfileobj(fin, fout, 8 << 20)
        pool_seq = kc.WORK_ROOT / "cov_pool.seq"
        kc.fasta_to_flat_seq(pool_fa, pool_seq)
        pool_fa.unlink(missing_ok=True)
        avail = pool_seq.stat().st_size

        for mb_tgt, xmx in [(100, "-Xmx128m"), (300, "-Xmx128m"), (500, "-Xmx256m"), (1000, "-Xmx512m"), (1000, "-Xmx1g")]:
            cap = mb_tgt * 1_000_000
            if avail < cap:
                continue
            sub_seq = kc.WORK_ROOT / f"heap_{mb_tgt}M.seq"
            with pool_seq.open("rb") as fin, sub_seq.open("wb") as fout:
                rem = cap
                while rem > 0:
                    buf = fin.read(min(rem, 8 << 20))
                    if not buf:
                        break
                    fout.write(buf)
                    rem -= len(buf)
            kc.log(f"=== Bounded-Heap Law VII Stress: {mb_tgt} Mbp @ {xmx} ===")
            hres = kc.run_hgs2_compress_verify(
                tools, sub_seq, mode="collection", threads=min(cpus, 4), xmx=xmx, bounded=True, timeout=1800
            )
            hres["dataset_mbp"] = mb_tgt
            adv["bounded_heap_matrix"].append(hres)
            if xmx in ("-Xmx128m", "-Xmx512m"):
                nres = kc.run_hgs2_native_compress_verify(
                    tools, sub_seq, mode="collection", threads=min(cpus, 4), bounded=True, timeout=1800
                )
                nres["dataset_mbp"] = mb_tgt
                adv["bounded_heap_matrix"].append(nres)
            sub_seq.unlink(missing_ok=True)
            save_state(state)
        pool_seq.unlink(missing_ok=True)


def run_bacterial_world_shootout(tools: dict, state: dict, src_root: Path, cpus: int) -> None:
    """Runs complete Mode S, Mode R, and Mode C shootouts on locked NCBI bacterial pangenomes."""
    manifests = {}
    for mp in sorted(src_root.rglob("manifest.json")):
        if mp.parent.name in ("pangenome", "salmonella", "klebsiella"):
            manifests[mp.parent.name] = mp

    for sp_key in ("pangenome", "salmonella", "klebsiella"):
        mp = manifests.get(sp_key)
        if not mp or budget_remaining_s() < 1800:
            continue
        sp_name = {"pangenome": "E_coli_24g", "salmonella": "S_enterica_25g", "klebsiella": "K_pneumoniae_23g"}[sp_key]
        kc.log(f"=== Fetching Locked NCBI Bacterial Pangenome: {sp_name} ===")
        sp_fa = kc.WORK_ROOT / f"{sp_key}.fasta"
        try:
            fmeta = wc.fetch_ncbi_locked_manifest(mp, sp_fa)
        except Exception as exc:
            kc.log(f"Warning: NCBI fetch failed for {sp_name}: {exc}")
            state["bacterial_shootout"].append({"species": sp_name, "status": "fetch_error", "error": str(exc)[-400:]})
            save_state(state)
            continue

        sp_seq = kc.WORK_ROOT / f"{sp_key}.seq"
        kc.fasta_to_flat_seq(sp_fa, sp_seq)

        rec0_fa = kc.WORK_ROOT / f"{sp_key}_r0.fa"
        rest_fa = kc.WORK_ROOT / f"{sp_key}_panel.fa"
        _extract_first_record_fasta(sp_fa, rec0_fa, rest_fa)
        rec0_seq = kc.WORK_ROOT / f"{sp_key}_r0.seq"
        rest_seq = kc.WORK_ROOT / f"{sp_key}_panel.seq"
        kc.fasta_to_flat_seq(rec0_fa, rec0_seq)
        kc.fasta_to_flat_seq(rest_fa, rest_seq)

        # Pick 1st donor record for GeCo3 -r pairwise neural referential comparison
        donor0_fa = kc.WORK_ROOT / f"{sp_key}_d0.fa"
        _extract_first_record_fasta(rest_fa, donor0_fa, kc.WORK_ROOT / f"{sp_key}_drop.fa")
        (kc.WORK_ROOT / f"{sp_key}_drop.fa").unlink(missing_ok=True)
        donor0_seq = kc.WORK_ROOT / f"{sp_key}_d0.seq"
        kc.fasta_to_flat_seq(donor0_fa, donor0_seq)
        donor0_fa.unlink(missing_ok=True)

        kc.log(f"=== {sp_name}: Mode S (Target 0, {rec0_seq.stat().st_size / 1e6:.2f} MB) — ALL 14+ Codecs ===")
        mode_s_seq = [
            kc.run_hgs2_compress_verify(tools, rec0_seq, mode="single", threads=cpus, xmx="-Xmx8g"),
            kc.run_hgs2_native_compress_verify(tools, rec0_seq, mode="single", threads=cpus),
            wc.run_hgs1_baseline(tools, rec0_seq, single=True, threads=cpus, xmx="-Xmx8g"),
            wc.run_grass_baseline(tools, rec0_seq, xmx="-Xmx8g", timeout=600),
            kc.run_jarvis3_baseline(tools, rec0_seq, level=7, timeout=600),
            kc.run_jarvis3_baseline(tools, rec0_seq, level=31, timeout=600),
            wc.run_geco3_baseline(tools, rec0_seq, level=5, timeout=600),
            wc.run_bsc_baseline(tools, rec0_seq, timeout=300),
            wc.run_named_general_baseline(tools, rec0_seq, "zpaq", timeout=600),
            wc.run_named_general_baseline(tools, rec0_seq, "xz", timeout=300),
            wc.run_named_general_baseline(tools, rec0_seq, "zstd", timeout=300),
            wc.run_named_general_baseline(tools, rec0_seq, "brotli", timeout=300),
            wc.run_named_general_baseline(tools, rec0_seq, "bzip2", timeout=300),
            wc.run_named_general_baseline(tools, rec0_seq, "gzip", timeout=300),
        ]
        mode_s_fa = [
            kc.run_hgs2_compress_verify(tools, rec0_fa, mode="single", threads=cpus, xmx="-Xmx8g"),
            kc.run_hgs2_native_compress_verify(tools, rec0_fa, mode="single", threads=cpus),
            wc.run_grass_baseline(tools, rec0_fa, xmx="-Xmx8g", timeout=600),
            wc.run_mfcompress_baseline(tools, rec0_fa, level=3, threads=cpus),
            wc.run_naf_baseline(tools, rec0_fa, level=22, long_win=31),
            wc.run_mbgc_baseline(tools, rec0_fa, mode=3, threads=cpus),
            wc.run_ffc_baseline(tools, rec0_fa, threads=cpus),
            wc.run_bsc_baseline(tools, rec0_fa, timeout=300),
            wc.run_named_general_baseline(tools, rec0_fa, "zstd", timeout=300),
            wc.run_named_general_baseline(tools, rec0_fa, "xz", timeout=300),
        ]

        # Mode R: HGS2 conditioned on individual N-1 donor .seq files + single-donor GeCo3 -r + AGC marginal!
        pdir = kc.WORK_ROOT / f"panel_{sp_key}"
        pdir_fa = kc.WORK_ROOT / f"panel_fa_{sp_key}"
        split_fas = wc.split_fasta_records(rest_fa, pdir_fa, max_split=64)
        pdir.mkdir(parents=True, exist_ok=True)
        for di, dfa in enumerate(split_fas):
            kc.fasta_to_flat_seq(dfa, pdir / f"donor_{di:03d}.seq")
        p1dir = kc.WORK_ROOT / f"panel1_{sp_key}"
        p1dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(donor0_seq, p1dir / "donor0.seq")
        kc.log(f"=== {sp_name}: Mode R (Target 0 | {len(split_fas)}-Donor Panel vs GeCo3 -r & AGC marginal) ===")
        mode_r_runs = [
            kc.run_hgs2_compress_verify(tools, rec0_seq, mode="ref", threads=cpus, xmx="-Xmx8g", panel=pdir),
            kc.run_hgs2_native_compress_verify(tools, rec0_seq, mode="ref", threads=cpus, panel=pdir),
            kc.run_hgs2_compress_verify(tools, rec0_seq, mode="ref", threads=cpus, xmx="-Xmx8g", panel=p1dir),
            kc.run_hgs2_native_compress_verify(tools, rec0_seq, mode="ref", threads=cpus, panel=p1dir),
            wc.run_geco3_baseline(tools, rec0_seq, ref_seq=donor0_seq, timeout=600),
        ]
        if len(split_fas) >= 2:
            mode_r_runs.append(
                wc.run_agc_marginal_mod_r(tools, split_fas[0], split_fas[1:], rec0_fa, adaptive=False, threads=cpus)
            )
            mode_r_runs.append(
                wc.run_agc_marginal_mod_r(tools, split_fas[0], split_fas[1:], rec0_fa, adaptive=True, threads=cpus)
            )
        shutil.rmtree(pdir, ignore_errors=True)
        shutil.rmtree(pdir_fa, ignore_errors=True)
        shutil.rmtree(p1dir, ignore_errors=True)

        kc.log(f"=== {sp_name}: Mode C Full Collection ({fmeta['bytes'] / (1<<20):.1f} MiB) — World Shootout ===")
        mode_c_fa = [
            kc.run_hgs2_compress_verify(tools, sp_fa, mode="collection", threads=cpus, xmx="-Xmx16g"),
            kc.run_hgs2_compress_verify(tools, sp_fa, mode="collection", threads=min(cpus, 4), xmx="-Xmx256m", bounded=True),
            kc.run_hgs2_native_compress_verify(tools, sp_fa, mode="collection", threads=cpus),
            kc.run_hgs2_native_compress_verify(tools, sp_fa, mode="collection", threads=min(cpus, 4), bounded=True),
            wc.run_hgs1_baseline(tools, sp_fa, threads=cpus, xmx="-Xmx16g"),
            wc.run_grass_baseline(tools, sp_fa, xmx="-Xmx16g", timeout=900),
            wc.run_agc_from_multifasta(tools, sp_fa, adaptive=False, threads=cpus),
            wc.run_agc_from_multifasta(tools, sp_fa, adaptive=True, threads=cpus),
            wc.run_mbgc_baseline(tools, sp_fa, mode=3, threads=cpus),
            wc.run_ffc_baseline(tools, sp_fa, threads=cpus),
            wc.run_naf_baseline(tools, sp_fa, level=22, long_win=31),
            wc.run_mfcompress_baseline(tools, sp_fa, level=3, threads=cpus),
            wc.run_bsc_baseline(tools, sp_fa),
            wc.run_named_general_baseline(tools, sp_fa, "zstd"),
            kc.run_zstd_baseline(tools, sp_fa, level=19),
            wc.run_named_general_baseline(tools, sp_fa, "xz"),
            wc.run_named_general_baseline(tools, sp_fa, "brotli"),
            wc.run_named_general_baseline(tools, sp_fa, "bzip2"),
            wc.run_named_general_baseline(tools, sp_fa, "gzip"),
        ]
        mode_c_seq = [
            kc.run_hgs2_compress_verify(tools, sp_seq, mode="collection", threads=cpus, xmx="-Xmx16g"),
            kc.run_hgs2_native_compress_verify(tools, sp_seq, mode="collection", threads=cpus),
            wc.run_grass_baseline(tools, sp_seq, xmx="-Xmx16g", timeout=900),
            kc.run_jarvis3_baseline(tools, sp_seq, level=7, timeout=900),
            wc.run_geco3_baseline(tools, sp_seq, level=5, timeout=900),
        ]

        state["bacterial_shootout"].append(
            {
                "species": sp_name,
                "fasta_meta": fmeta,
                "mode_s_seq": mode_s_seq,
                "mode_s_fa": mode_s_fa,
                "mode_r": mode_r_runs,
                "mode_c_fa": mode_c_fa,
                "mode_c_seq": mode_c_seq,
            }
        )
        for tmp_p in (sp_fa, sp_seq, rec0_fa, rest_fa, rec0_seq, rest_seq, donor0_seq):
            tmp_p.unlink(missing_ok=True)
        save_state(state)


def run_viral_world_shootout(tools: dict, state: dict, cov_dir: Path | None, cpus: int) -> None:
    """Runs SARS-CoV-2 viral pangenome world shootout from 1.7 MiB up to ~1.4 GiB."""
    if cov_dir is None:
        return
    for shard in ["group_19.fasta", "group_3.fasta", "group_22.fasta"]:
        if budget_remaining_s() < 1500:
            break
        p = cov_dir / shard
        if not p.exists():
            continue
        sz = p.stat().st_size
        kc.log(f"=== Viral World Shootout: {shard} ({sz / (1<<20):.1f} MiB) ===")
        v_fa = [
            kc.run_hgs2_compress_verify(tools, p, mode="collection", threads=cpus, xmx="-Xmx16g"),
            kc.run_hgs2_compress_verify(tools, p, mode="collection", threads=min(cpus, 4), xmx="-Xmx256m", bounded=True),
            kc.run_hgs2_native_compress_verify(tools, p, mode="collection", threads=cpus),
            kc.run_hgs2_native_compress_verify(tools, p, mode="collection", threads=min(cpus, 4), bounded=True),
            wc.run_hgs1_baseline(tools, p, threads=cpus, xmx="-Xmx16g"),
            wc.run_agc_from_multifasta(tools, p, adaptive=True, threads=cpus),
            wc.run_mbgc_baseline(tools, p, mode=3, threads=cpus),
            wc.run_ffc_baseline(tools, p, threads=cpus),
            wc.run_naf_baseline(tools, p, level=22, long_win=31),
            wc.run_mfcompress_baseline(tools, p, level=2 if sz > 150_000_000 else 3, threads=cpus),
            wc.run_bsc_baseline(tools, p),
            wc.run_named_general_baseline(tools, p, "zstd"),
            wc.run_named_general_baseline(tools, p, "xz"),
            wc.run_named_general_baseline(tools, p, "brotli"),
            wc.run_named_general_baseline(tools, p, "bzip2"),
            wc.run_named_general_baseline(tools, p, "gzip"),
        ]
        if sz <= 320_000_000:
            v_fa.append(wc.run_grass_baseline(tools, p, xmx="-Xmx18g", timeout=900))
        if sz <= 5_000_000:
            v_fa.append(wc.run_named_general_baseline(tools, p, "zpaq", timeout=600))

        v_seq_runs: list[dict] = []
        v_seq = kc.WORK_ROOT / f"{p.stem}.seq"
        vmeta = kc.fasta_to_flat_seq(p, v_seq)
        v_seq_runs.append(kc.run_hgs2_compress_verify(tools, v_seq, mode="collection", threads=cpus, xmx="-Xmx16g"))
        v_seq_runs.append(kc.run_hgs2_native_compress_verify(tools, v_seq, mode="collection", threads=cpus))
        if sz <= 165_000_000:
            v_seq_runs.append(wc.run_grass_baseline(tools, v_seq, xmx="-Xmx16g", timeout=900))
            v_seq_runs.append(kc.run_jarvis3_baseline(tools, v_seq, level=7, timeout=900))
            v_seq_runs.append(wc.run_geco3_baseline(tools, v_seq, level=5, timeout=900))
        if sz <= 65_000_000:
            v_seq_runs.append(kc.run_jarvis3_baseline(tools, v_seq, level=31, timeout=900))
        v_seq.unlink(missing_ok=True)

        state["viral_shootout"].append(
            {"shard": shard, "raw_bytes": sz, "flat_bases": vmeta["bases"], "runs_fa": v_fa, "runs_seq": v_seq_runs}
        )
        save_state(state)

    # Massive 1.4 GiB Viral Cohort (~4 groups)
    if budget_remaining_s() >= 2400:
        viral_mega = kc.WORK_ROOT / "covid_mega.fasta"
        with viral_mega.open("wb") as fout:
            for sh in ["group_36.fasta", "group_34.fasta", "group_22.fasta", "group_29.fasta"]:
                sp = cov_dir / sh
                if sp.exists():
                    with sp.open("rb") as fin:
                        shutil.copyfileobj(fin, fout, 8 << 20)
        msz = viral_mega.stat().st_size
        if msz > 200_000_000:
            kc.log(f"=== Viral Mega-Cohort Shootout ({msz / (1<<20):.1f} MiB) ===")
            mruns = [
                kc.run_hgs2_compress_verify(tools, viral_mega, mode="collection", threads=cpus, xmx="-Xmx22g", timeout=1800),
                kc.run_hgs2_compress_verify(tools, viral_mega, mode="collection", threads=min(cpus, 4), xmx="-Xmx512m", bounded=True, timeout=1800),
                kc.run_hgs2_native_compress_verify(tools, viral_mega, mode="collection", threads=cpus, timeout=1800),
                kc.run_hgs2_native_compress_verify(tools, viral_mega, mode="collection", threads=min(cpus, 4), bounded=True, timeout=1800),
                wc.run_mbgc_baseline(tools, viral_mega, mode=2, threads=cpus, timeout=1800),
                wc.run_ffc_baseline(tools, viral_mega, threads=cpus, timeout=1800),
                wc.run_naf_baseline(tools, viral_mega, level=22, long_win=31, timeout=1800),
                kc.run_zstd_baseline(tools, viral_mega, level=19, timeout=900),
                wc.run_bsc_baseline(tools, viral_mega, timeout=1800),
            ]
            state["viral_shootout"].append({"shard": "covid_mega_4groups", "raw_bytes": msz, "runs_fa": mruns})
            save_state(state)
        viral_mega.unlink(missing_ok=True)


def run_human_world_shootout(tools: dict, state: dict, g_dir: Path | None, cpus: int) -> None:
    """Runs Human GRCh38 single chromosomes, CHM13|GRCh38 Mode R, and Diploid Mode C shootouts."""
    if g_dir is None:
        return
    for cname in ["21", "22", "1"]:
        if budget_remaining_s() < 1500:
            break
        fa = g_dir / f"Homo_sapiens.GRCh38.dna.chromosome.{cname}.fa"
        if not fa.exists():
            continue
        fsz = fa.stat().st_size
        kc.log(f"=== Human GRCh38 Chr{cname} ({fsz / (1<<20):.1f} MiB) Shootout ===")
        runs_fa = [
            kc.run_hgs2_compress_verify(tools, fa, mode="collection", threads=cpus, xmx="-Xmx16g"),
            kc.run_hgs2_compress_verify(tools, fa, mode="collection", threads=min(cpus, 4), xmx="-Xmx256m", bounded=True),
            kc.run_hgs2_native_compress_verify(tools, fa, mode="collection", threads=cpus),
            kc.run_hgs2_native_compress_verify(tools, fa, mode="collection", threads=min(cpus, 4), bounded=True),
            wc.run_ffc_baseline(tools, fa, threads=cpus),
            wc.run_naf_baseline(tools, fa, level=22, long_win=31),
            wc.run_bsc_baseline(tools, fa),
            kc.run_zstd_baseline(tools, fa, level=19),
            kc.run_xz_baseline(tools, fa, level=6),
        ]
        if fsz <= 320_000_000:
            runs_fa.append(wc.run_grass_baseline(tools, fa, xmx="-Xmx18g", timeout=900))
        if fsz <= 75_000_000:
            runs_fa.extend([
                wc.run_mfcompress_baseline(tools, fa, level=3, threads=cpus),
                wc.run_named_general_baseline(tools, fa, "zstd"),
                wc.run_named_general_baseline(tools, fa, "brotli"),
                wc.run_named_general_baseline(tools, fa, "bzip2"),
                wc.run_named_general_baseline(tools, fa, "gzip"),
            ])

        runs_seq: list[dict] = []
        seq = kc.WORK_ROOT / f"grch38_chr{cname}.seq"
        kc.fasta_to_flat_seq(fa, seq)
        runs_seq.append(kc.run_hgs2_compress_verify(tools, seq, mode="collection", threads=cpus, xmx="-Xmx16g"))
        runs_seq.append(kc.run_hgs2_compress_verify(tools, seq, mode="collection", threads=min(cpus, 4), xmx="-Xmx256m", bounded=True))
        runs_seq.append(kc.run_hgs2_native_compress_verify(tools, seq, mode="collection", threads=cpus))
        runs_seq.append(kc.run_hgs2_native_compress_verify(tools, seq, mode="collection", threads=min(cpus, 4), bounded=True))
        if fsz <= 75_000_000:
            runs_seq.append(wc.run_grass_baseline(tools, seq, xmx="-Xmx16g", timeout=900))
            runs_seq.append(kc.run_jarvis3_baseline(tools, seq, level=7, timeout=900))
            runs_seq.append(kc.run_jarvis3_baseline(tools, seq, level=31, timeout=1200))
            runs_seq.append(wc.run_geco3_baseline(tools, seq, level=5, timeout=900))
        seq.unlink(missing_ok=True)
        state["human_shootout"]["single_chromosomes"].append({"target": f"GRCh38_chr{cname}", "raw_bytes": fsz, "runs_fa": runs_fa, "runs_seq": runs_seq})
        save_state(state)

    if budget_remaining_s() < 2400:
        return

    chm13_dir = kc.WORK_ROOT / "chm13_chr"
    chm13_map = {}
    try:
        chm13_map = kc.download_and_split_chm13(chm13_dir)
    except Exception as exc:
        kc.log(f"Warning: CHM13 download skipped: {exc}")

    if chm13_map:
        for cname in ["21", "22", "1"]:
            if budget_remaining_s() < 1500:
                break
            g_fa = g_dir / f"Homo_sapiens.GRCh38.dna.chromosome.{cname}.fa"
            c_fa = chm13_map.get(f"chr{cname}")
            if not (g_fa.exists() and c_fa and c_fa.exists()):
                continue
            # Cross-Human Mode R: HGS2 vs GeCo3 -r (on Chr21/Chr22 <= 65M .seq)
            pdir = kc.WORK_ROOT / f"hr_panel_{cname}"
            pdir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(g_fa, pdir / f"grch38_chr{cname}.fa")
            kc.log(f"=== Cross-Human Mode R: CHM13 chr{cname} | GRCh38 chr{cname} ===")
            mr_runs = [
                kc.run_hgs2_compress_verify(tools, c_fa, mode="ref", threads=cpus, xmx="-Xmx16g", panel=pdir),
                kc.run_hgs2_native_compress_verify(tools, c_fa, mode="ref", threads=cpus, panel=pdir),
            ]
            shutil.rmtree(pdir, ignore_errors=True)
            if cname in ("21", "22"):
                g_seq = kc.WORK_ROOT / f"mr_g{cname}.seq"
                c_seq = kc.WORK_ROOT / f"mr_c{cname}.seq"
                kc.fasta_to_flat_seq(g_fa, g_seq)
                kc.fasta_to_flat_seq(c_fa, c_seq)
                mr_runs.append(wc.run_geco3_baseline(tools, c_seq, ref_seq=g_seq, timeout=900))
                g_seq.unlink(missing_ok=True)
                c_seq.unlink(missing_ok=True)
            state["human_shootout"]["cross_human_mode_r"].append({"pair": f"CHM13_chr{cname}_given_GRCh38_chr{cname}", "runs": mr_runs})
            save_state(state)

        for cname in ["21", "22", "1"]:
            if budget_remaining_s() < 1500:
                break
            g_fa = g_dir / f"Homo_sapiens.GRCh38.dna.chromosome.{cname}.fa"
            c_fa = chm13_map.get(f"chr{cname}")
            if not (g_fa.exists() and c_fa and c_fa.exists()):
                continue
            merged = kc.WORK_ROOT / f"diploid_chr{cname}.fa"
            with merged.open("wb") as fout:
                for src in (g_fa, c_fa):
                    with src.open("rb") as fin:
                        shutil.copyfileobj(fin, fout, 8 << 20)
            kc.log(f"=== Human Diploid Mode C: Chr{cname} ({merged.stat().st_size / (1<<20):.1f} MiB) ===")
            d_runs = [
                kc.run_hgs2_compress_verify(tools, merged, mode="collection", threads=cpus, xmx="-Xmx22g", timeout=1800),
                kc.run_hgs2_native_compress_verify(tools, merged, mode="collection", threads=cpus, timeout=1800),
                kc.run_agc_collection(tools, g_fa, [c_fa], adaptive=False, threads=cpus),
                kc.run_agc_collection(tools, g_fa, [c_fa], adaptive=True, threads=cpus),
                wc.run_mbgc_baseline(tools, merged, mode=3, threads=cpus),
                wc.run_ffc_baseline(tools, merged, threads=cpus),
                wc.run_naf_baseline(tools, merged, level=22, long_win=31),
                wc.run_bsc_baseline(tools, merged),
                kc.run_zstd_baseline(tools, merged, level=19),
            ]
            if cname in ("21", "22"):
                d_runs.append(wc.run_grass_baseline(tools, merged, xmx="-Xmx16g", timeout=900))
                d_runs.append(wc.run_mfcompress_baseline(tools, merged, level=3, threads=cpus))
                d_seq = kc.WORK_ROOT / f"diploid_chr{cname}.seq"
                kc.fasta_to_flat_seq(merged, d_seq)
                d_runs.append(kc.run_jarvis3_baseline(tools, d_seq, level=7, timeout=900))
                d_runs.append(wc.run_geco3_baseline(tools, d_seq, level=5, timeout=900))
                d_seq.unlink(missing_ok=True)
            merged.unlink(missing_ok=True)
            state["human_shootout"]["diploid_mode_c"].append({"label": f"Diploid_Chr{cname}", "runs": d_runs})
            save_state(state)

        shutil.rmtree(chm13_dir, ignore_errors=True)


def main() -> None:
    ds = kc.discover_kaggle_datasets()
    assert ds["src"] is not None, "hypergrass-eimdl-src dataset not mounted"
    state: dict = {
        "benchmark": "kaggle_world_leading_shootout",
        "policy": "kaggle_only_zero_local_testing",
        "system": kc.collect_system_info(),
        "datasets_mounted": {k: str(v) for k, v in ds.items()},
        "tools": {},
        "adversarial_and_bounded": {"white_noise_10mb": [], "pathology_fasta": [], "bounded_heap_matrix": []},
        "bacterial_shootout": [],
        "viral_shootout": [],
        "human_shootout": {"single_chromosomes": [], "cross_human_mode_r": [], "diploid_mode_c": []},
    }
    save_state(state)

    base_tools = kc.ensure_jdk_and_tools(ds["src"])
    tools = wc.ensure_world_tools(base_tools)
    state["tools"] = {k: v for k, v in tools.items() if not k.startswith("_")}
    save_state(state)

    cpus = max(1, os.cpu_count() or 4)
    kc.log(f"Kaggle CPU Shootout Ready: {cpus} CPUs | smoke={tools.get('smoke_checks')}")

    cov_dir = next(ds["covid19"].rglob("group_1.fasta")).parent if ds["covid19"] is not None else None
    g_dir = next(ds["grch38"].rglob("Homo_sapiens.GRCh38.dna.chromosome.1.fa")).parent if ds["grch38"] is not None else None

    run_adversarial_and_bounded_suite(tools, state, cov_dir, cpus)
    run_bacterial_world_shootout(tools, state, ds["src"], cpus)
    run_viral_world_shootout(tools, state, cov_dir, cpus)
    run_human_world_shootout(tools, state, g_dir, cpus)

    kc.log(f"WORLD SHOOTOUT COMPLETE in {round(time.monotonic() - SESSION_START, 1)}s -> {OUT_JSON}")


if __name__ == "__main__":
    main()
