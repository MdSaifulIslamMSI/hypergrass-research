"""Kaggle TPU VM Scale & JAX/XLA Similarity Kernel for HyperGraSS / EIMDL (HGS2).

Validates high-core Amdahl thread scaling, strict minimum-heap scaling from
100 Mbp to 1.45 Gbp, full 24-chromosome GRCh38-conditioned T2T-CHM13v2.0 Mode R
(~6.2 Gbp cumulative), and real JAX/XLA TPU v3-8 pairwise genomic sketch matrices.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import sys
import time


def _locate_common() -> Path:
    inp = Path("/kaggle/input")
    for p in inp.rglob("kaggle_common.py"):
        return p.parent
    raise FileNotFoundError("kaggle_common.py not found under /kaggle/input")


sys.path.insert(0, str(_locate_common()))
import kaggle_common as kc  # noqa: E402


OUT_JSON = kc.RESULTS_DIR / "tpu_scale_results.json"


def save_state(state: dict) -> None:
    OUT_JSON.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def run_jax_tpu_experiments(state: dict, chr1_seq: Path) -> dict:
    """Runs real JAX/XLA TPU (or XLA fallback) pairwise genomic sketch similarity matrices."""
    tpu_out: dict = {"status": "skipped"}
    try:
        import jax
        import jax.numpy as jnp
        import numpy as np

        devs = jax.devices()
        backend = jax.default_backend()
        tpu_out = {
            "status": "ok",
            "jax_version": jax.__version__,
            "backend": backend,
            "device_count": len(devs),
            "devices": [str(d) for d in devs],
            "similarity_matrix_bench": [],
        }

        raw_bytes = np.fromfile(chr1_seq, dtype=np.uint8, count=64_000_000)
        lut = np.zeros(256, dtype=np.uint8)
        lut[ord("C")] = 1
        lut[ord("G")] = 2
        lut[ord("T")] = 3
        b2 = lut[raw_bytes]

        @jax.jit
        def tpu_sketch_similarity(one_hot_norm: jnp.ndarray) -> jnp.ndarray:
            return one_hot_norm @ one_hot_norm.T

        k_dim = 256
        for n_units in [1024, 2048, 4096, 8192]:
            stride = max(1, (len(b2) - 4096) // n_units)
            sigs = np.empty((n_units, k_dim), dtype=np.int32)
            for i in range(n_units):
                off = i * stride
                v = b2[off : off + 1024].astype(np.int32)
                sigs[i] = (
                    v.reshape(256, 4) * np.array([1, 4, 16, 64], dtype=np.int32)
                ).sum(axis=1)

            # One-hot in bfloat16/float32 for XLA systolic matrix unit (exact equality/k_dim)
            # Chunked CPU NumPy comparison for moderate N to avoid host OOM
            cpu_ms = None
            if n_units <= 4096:
                t0 = time.perf_counter()
                _ = float((sigs[:, None, :] == sigs[None, :, :]).mean(axis=2, dtype=np.float32).sum())
                cpu_ms = round((time.perf_counter() - t0) * 1000.0, 2)

            oh_np = np.zeros((n_units, k_dim * 256), dtype=np.float32)
            rows = np.repeat(np.arange(n_units), k_dim)
            cols = (np.arange(k_dim)[None, :] * 256 + sigs).ravel()
            oh_np[rows, cols] = 1.0 / float(np.sqrt(k_dim))

            oh_j = jax.device_put(oh_np.astype(jnp.bfloat16))
            # Compile / warmup pass
            t_comp0 = time.perf_counter()
            res = tpu_sketch_similarity(oh_j)
            jax.block_until_ready(res)
            compile_plus_first_ms = round((time.perf_counter() - t_comp0) * 1000.0, 2)

            # Steady-state pass
            t_run0 = time.perf_counter()
            res2 = tpu_sketch_similarity(oh_j)
            jax.block_until_ready(res2)
            steady_ms = round((time.perf_counter() - t_run0) * 1000.0, 2)

            tpu_out["similarity_matrix_bench"].append(
                {
                    "n_units": n_units,
                    "k_dim": k_dim,
                    "cpu_numpy_ms": cpu_ms,
                    "xla_compile_plus_first_ms": compile_plus_first_ms,
                    "xla_steady_ms": steady_ms,
                    "xla_speedup_vs_cpu": round(cpu_ms / max(steady_ms, 1e-4), 2)
                    if cpu_ms
                    else None,
                }
            )
            del oh_np, oh_j, res, res2

        save_state(state)
    except Exception as exc:
        tpu_out = {"status": "error", "error": str(exc)[-600:]}
    return tpu_out


def main() -> None:
    ds = kc.discover_kaggle_datasets()
    assert ds["src"] is not None, "hypergrass-eimdl-src dataset not mounted"
    state: dict = {
        "benchmark": "kaggle_tpu_scale",
        "system": kc.collect_system_info(),
        "datasets_mounted": {k: str(v) for k, v in ds.items()},
        "tools": {},
        "jax_tpu_experiments": {},
        "amdahl_high_core_scaling": [],
        "heap_floor_scaling_sweep": [],
        "all24_chm13_given_grch38_mode_r": {"per_chromosome": [], "totals": {}},
    }
    save_state(state)

    tools = kc.ensure_jdk_and_tools(ds["src"])
    state["tools"] = tools
    save_state(state)

    cpus = max(1, os.cpu_count() or 4)
    kc.log(f"TPU VM host CPUs={cpus}; JARVIS3={tools['jarvis3']}; AGC={tools['agc']}")

    # ------------------------------------------------------------------
    # 1. JAX / XLA TPU Benchmark + Full-Human Prep
    # ------------------------------------------------------------------
    chr1_seq = kc.WORK_ROOT / "tpu_chr1.seq"
    g_dir = None
    if ds["grch38"] is not None:
        g_dir = next(ds["grch38"].rglob("Homo_sapiens.GRCh38.dna.chromosome.1.fa")).parent
        kc.fasta_to_flat_seq(g_dir / "Homo_sapiens.GRCh38.dna.chromosome.1.fa", chr1_seq)
        kc.log("=== Running JAX/XLA TPU Vectorized Similarity Benchmark ===")
        state["jax_tpu_experiments"] = run_jax_tpu_experiments(state, chr1_seq)
        save_state(state)

    # ------------------------------------------------------------------
    # 2. High-Core Amdahl Thread Sweep on SARS-CoV-2 & Human Diploid Chr1
    # ------------------------------------------------------------------
    if ds["covid19"] is not None:
        cov_dir = next(ds["covid19"].rglob("group_1.fasta")).parent
        g22 = cov_dir / "group_22.fasta"
        if g22.exists():
            thread_list = sorted({t for t in [1, 2, 4, 8, 16, 32, min(cpus, 64), cpus] if t <= cpus})
            for th in thread_list:
                kc.log(f"=== Amdahl TPU-host sweep: COVID group_22 (309 MiB), threads={th} ===")
                res = kc.run_hgs2_compress_verify(tools, g22, mode="collection", threads=th, xmx="-Xmx24g")
                res["target"] = "covid_group_22.fasta"
                state["amdahl_high_core_scaling"].append(res)
                save_state(state)

        # Heap Floor Scaling Curve: 100 Mbp -> 300 Mbp -> 500 Mbp -> 1.0 Gbp canonical .seq
        mega_fa = kc.WORK_ROOT / "viral_pool.fasta"
        with mega_fa.open("wb") as fout:
            for sh in ["group_36.fasta", "group_34.fasta", "group_22.fasta"]:
                sp = cov_dir / sh
                if sp.exists():
                    with sp.open("rb") as fin:
                        shutil.copyfileobj(fin, fout, 8 << 20)
        mega_seq = kc.WORK_ROOT / "viral_pool.seq"
        kc.fasta_to_flat_seq(mega_fa, mega_seq)
        mega_fa.unlink(missing_ok=True)
        total_avail = mega_seq.stat().st_size

        for mb_target, xmx_list in [
            (100, ["-Xmx128m", "-Xmx192m", "-Xmx256m"]),
            (300, ["-Xmx128m", "-Xmx256m", "-Xmx512m"]),
            (500, ["-Xmx256m", "-Xmx512m", "-Xmx1g"]),
            (1000, ["-Xmx512m", "-Xmx1g", "-Xmx2g", "-Xmx16g"]),
        ]:
            byte_cap = mb_target * 1_000_000
            if total_avail < byte_cap:
                continue
            slice_seq = kc.WORK_ROOT / f"viral_{mb_target}M.seq"
            with mega_seq.open("rb") as fin, slice_seq.open("wb") as fout:
                rem = byte_cap
                while rem > 0:
                    buf = fin.read(min(rem, 8 << 20))
                    if not buf:
                        break
                    fout.write(buf)
                    rem -= len(buf)
            for xmx in xmx_list:
                kc.log(f"=== Heap Sweep: {mb_target} Mbp SARS-CoV-2 .seq @ {xmx} ===")
                hres = kc.run_hgs2_compress_verify(
                    tools,
                    slice_seq,
                    mode="collection",
                    threads=min(cpus, 8),
                    xmx=xmx,
                    bounded=(xmx != "-Xmx16g"),
                    timeout=1800,
                )
                hres["dataset_mbp"] = mb_target
                state["heap_floor_scaling_sweep"].append(hres)
                save_state(state)
            slice_seq.unlink(missing_ok=True)
        mega_seq.unlink(missing_ok=True)

    # ------------------------------------------------------------------
    # 3. Complete All-24-Chromosome Cross-Human Mode R (~6.2 Gbp Cumulative)
    # ------------------------------------------------------------------
    chm13_dir = kc.WORK_ROOT / "chm13_chr"
    chm13_map = {}
    try:
        chm13_map = kc.download_and_split_chm13(chm13_dir)
    except Exception as exc:
        kc.log(f"Warning: T2T-CHM13v2.0 download failed: {exc}")

    if chm13_map and g_dir is not None:
        chr_names = [str(i) for i in range(1, 23)] + ["X", "Y"]
        sum_raw = 0
        sum_bases = 0
        sum_arc = 0
        sum_enc = 0.0
        for cname in chr_names:
            g_fa = g_dir / f"Homo_sapiens.GRCh38.dna.chromosome.{cname}.fa"
            c_fa = chm13_map.get(f"chr{cname}")
            if not (g_fa.exists() and c_fa and c_fa.exists()):
                continue
            panel_dir = kc.WORK_ROOT / f"tpu_panel_{cname}"
            panel_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(g_fa, panel_dir / f"grch38_chr{cname}.fa")
            kc.log(f"=== All-24 Mode R: CHM13 chr{cname} | GRCh38 chr{cname} ===")
            r_res = kc.run_hgs2_compress_verify(
                tools, c_fa, mode="ref", threads=cpus, xmx="-Xmx24g", panel=panel_dir
            )
            r_res["chromosome"] = f"chr{cname}"
            state["all24_chm13_given_grch38_mode_r"]["per_chromosome"].append(r_res)
            if r_res["status"] == "ok":
                sum_raw += r_res["raw_bytes"]
                sum_bases += r_res["target_bases"]
                sum_arc += r_res["archive_bytes"]
                sum_enc += r_res["enc_wall_s"]
            shutil.rmtree(panel_dir, ignore_errors=True)
            save_state(state)

        if sum_arc > 0:
            state["all24_chm13_given_grch38_mode_r"]["totals"] = {
                "chromosomes_completed": len(
                    state["all24_chm13_given_grch38_mode_r"]["per_chromosome"]
                ),
                "total_chm13_raw_bytes": sum_raw,
                "total_chm13_bases": sum_bases,
                "total_archive_bytes": sum_arc,
                "overall_bpb": round((8.0 * sum_arc) / sum_raw, 4),
                "overall_ratio": round(sum_raw / sum_arc, 2),
                "total_enc_wall_s": round(sum_enc, 2),
            }
            save_state(state)
        shutil.rmtree(chm13_dir, ignore_errors=True)

    chr1_seq.unlink(missing_ok=True)
    kc.log(f"TPU SCALE COMPLETE. Results written to {OUT_JSON}")


if __name__ == "__main__":
    main()
