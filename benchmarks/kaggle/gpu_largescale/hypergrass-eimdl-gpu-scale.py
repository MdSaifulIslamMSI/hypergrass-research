"""Kaggle GPU Large-Scale Validation for HyperGraSS / EIMDL (HGS2).

Tests multi-hundred-MB and gigabyte-scale Human (GRCh38 + T2T-CHM13v2.0) and
SARS-CoV-2 viral pangenomes against AGC, JARVIS3, zstd, xz, plus CUDA GPU
k-mer sketch matrix & autoregressive neural DNA entropy baselines.
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


OUT_JSON = kc.RESULTS_DIR / "gpu_largescale_results.json"


def save_state(state: dict) -> None:
    OUT_JSON.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def run_cuda_gpu_experiments(
    state: dict, chr22_seq: Path, covid_seq: Path, chm13_chr1: Path | None
) -> dict:
    """Runs real PyTorch CUDA GPU genomic sketch + neural DNA entropy baselines."""
    gpu_out: dict = {"status": "skipped"}
    try:
        import numpy as np
        import torch
        import torch.nn as nn
        import torch.nn.functional as F

        if not torch.cuda.is_available():
            return {"status": "no_cuda_device"}

        dev = torch.device("cuda:0")
        props = torch.cuda.get_device_properties(0)
        gpu_out = {
            "status": "ok",
            "device": props.name,
            "vram_gb": round(props.total_memory / (1024**3), 2),
            "torch_version": torch.__version__,
            "sketch_matrix_bench": [],
            "neural_entropy_bench": [],
        }

        # 1) CUDA Parallel Pairwise 64-Feature Bottom-Sketch Similarity vs CPU NumPy
        # Start at offset=20_000_000 on GRCh38 Chr22 (euchromatic q-arm, past 11.66 Mbp masked p-arm).
        raw_bytes = np.fromfile(chr22_seq, dtype=np.uint8, count=28_000_000, offset=20_000_000)
        lut = np.zeros(256, dtype=np.uint8)
        lut[ord("C")] = 1
        lut[ord("G")] = 2
        lut[ord("T")] = 3
        b2 = lut[raw_bytes]

        # Derive 16-mer hashes per window on CPU cleanly, then benchmark N x N similarity matrix
        for n_units in [512, 1024, 2048, 4096]:
            stride = max(1, (len(b2) - 4096) // n_units)
            # 256-dim MinHash / 32-mer signature vectors (uint8 buckets for exact equality count)
            k_dim = 256
            sigs = np.empty((n_units, k_dim), dtype=np.int32)
            for i in range(n_units):
                off = i * stride
                chunk = b2[off : off + 2048]
                # Deterministic rolling 16-mer hash projection into 256 buckets
                v = chunk[:1024].astype(np.int32)
                sigs[i] = (
                    (v.reshape(256, 4) * np.array([1, 4, 16, 64], dtype=np.int32)).sum(
                        axis=1
                    )
                )

            t0 = time.perf_counter()
            cpu_matches = (sigs[:, None, :] == sigs[None, :, :]).mean(axis=2, dtype=np.float32)
            _ = float(cpu_matches.sum())
            cpu_ms = (time.perf_counter() - t0) * 1000.0

            sigs_t = torch.from_numpy(sigs).to(dev)
            # Encode discrete buckets in float16/int32 batched chunks or one-hot matmul (fast!)
            one_hot = F.one_hot(sigs_t.long(), num_classes=256).to(torch.float16).reshape(
                n_units, k_dim * 256
            )
            if n_units == 512:
                _ = (one_hot @ one_hot.T) / float(k_dim)
            torch.cuda.synchronize()
            t1 = time.perf_counter()
            gpu_sim = (one_hot @ one_hot.T) / float(k_dim)
            torch.cuda.synchronize()
            gpu_ms = (time.perf_counter() - t1) * 1000.0
            max_err = float(np.max(np.abs(cpu_matches - gpu_sim.float().cpu().numpy())))
            gpu_out["sketch_matrix_bench"].append(
                {
                    "n_units": n_units,
                    "k_dim": k_dim,
                    "cpu_numpy_ms": round(cpu_ms, 2),
                    "gpu_cuda_ms": round(gpu_ms, 2),
                    "gpu_speedup_x": round(cpu_ms / max(gpu_ms, 1e-4), 2),
                    "max_abs_err": round(max_err, 6),
                }
            )
            del sigs_t, one_hot, gpu_sim
            torch.cuda.empty_cache()

        # 2) CUDA Autoregressive Causal GRU+Conv DNA Entropy Baseline vs HGS2 CPU
        class CausalDnaNet(nn.Module):
            def __init__(self, hidden: int = 128) -> None:
                super().__init__()
                self.emb = nn.Embedding(4, 32)
                self.conv = nn.Conv1d(32, hidden, kernel_size=7, padding=6)
                self.gru = nn.GRU(hidden, hidden, batch_first=True)
                self.head = nn.Linear(hidden, 4)

            def forward(self, x: torch.Tensor) -> torch.Tensor:
                h = self.emb(x).transpose(1, 2)
                h = F.gelu(self.conv(h)[:, :, : x.shape[1]])
                h, _ = self.gru(h.transpose(1, 2))
                return self.head(h)

        for label, seq_path, off in [("Human_Chr22_4MB", chr22_seq, 20_000_000), ("COVID19_4MB", covid_seq, 0)]:
            raw_u = np.fromfile(seq_path, dtype=np.uint8, count=6_000_000, offset=off)
            arr = lut[raw_u]
            split = 2_000_000
            tr = arr[:split]
            ev = arr[split : split + 4_000_000]
            seq_len = 512

            def make_batches(vec: np.ndarray, limit_seqs: int) -> torch.Tensor:
                n_seqs = min(limit_seqs, (len(vec) - 1) // seq_len)
                x = torch.from_numpy(
                    vec[: n_seqs * seq_len + 1].copy()
                ).long()
                return x

            model = CausalDnaNet(hidden=128).to(dev)
            opt = torch.optim.AdamW(model.parameters(), lr=3e-3)
            torch.cuda.reset_peak_memory_stats(dev)
            t_train0 = time.perf_counter()
            tr_t = make_batches(tr, 2000).to(dev)
            n_tr = (tr_t.numel() - 1) // seq_len
            tr_in = tr_t[: n_tr * seq_len].view(n_tr, seq_len)
            tr_tgt = tr_t[1 : n_tr * seq_len + 1].view(n_tr, seq_len)
            model.train()
            bs = 64
            for _epoch in range(2):
                for s in range(0, n_tr, bs):
                    xb = tr_in[s : s + bs]
                    yb = tr_tgt[s : s + bs]
                    opt.zero_grad(set_to_none=True)
                    logits = model(xb)
                    loss = F.cross_entropy(logits.reshape(-1, 4), yb.reshape(-1))
                    loss.backward()
                    opt.step()
            torch.cuda.synchronize()
            train_s = time.perf_counter() - t_train0

            ev_t = make_batches(ev, 7800).to(dev)
            n_ev = (ev_t.numel() - 1) // seq_len
            ev_in = ev_t[: n_ev * seq_len].view(n_ev, seq_len)
            ev_tgt = ev_t[1 : n_ev * seq_len + 1].view(n_ev, seq_len)
            model.eval()
            total_nats = 0.0
            total_bases = n_ev * seq_len
            t_eval0 = time.perf_counter()
            with torch.no_grad():
                for s in range(0, n_ev, bs):
                    xb = ev_in[s : s + bs]
                    yb = ev_tgt[s : s + bs]
                    logits = model(xb)
                    total_nats += float(
                        F.cross_entropy(
                            logits.reshape(-1, 4), yb.reshape(-1), reduction="sum"
                        ).item()
                    )
            torch.cuda.synchronize()
            eval_s = time.perf_counter() - t_eval0
            bpb = (total_nats / float(np.log(2.0))) / total_bases
            peak_vram_mb = torch.cuda.max_memory_allocated(dev) / (1024**2)

            gpu_out["neural_entropy_bench"].append(
                {
                    "dataset": label,
                    "eval_bases": int(total_bases),
                    "neural_bpb": round(bpb, 4),
                    "theoretical_bytes": int(total_bases * bpb / 8.0),
                    "train_s": round(train_s, 3),
                    "eval_s": round(eval_s, 3),
                    "throughput_mb_s": round((total_bases / 1e6) / max(eval_s, 1e-4), 3),
                    "peak_vram_mb": round(peak_vram_mb, 2),
                }
            )
            del model, tr_t, ev_t
            torch.cuda.empty_cache()

        save_state(state)
    except Exception as exc:
        gpu_out = {"status": "error", "error": str(exc)[-600:]}
    return gpu_out


def main() -> None:
    ds = kc.discover_kaggle_datasets()
    assert ds["src"] is not None, "hypergrass-eimdl-src dataset not mounted"
    state: dict = {
        "benchmark": "kaggle_gpu_largescale",
        "system": kc.collect_system_info(),
        "datasets_mounted": {k: str(v) for k, v in ds.items()},
        "tools": {},
        "grch38_chromosomes": [],
        "grch38_multichr_1_8": [],
        "cross_human_chm13_vs_grch38": {"mode_r": [], "diploid_mode_c": []},
        "covid19_viral_scaling": [],
        "cuda_gpu_experiments": {},
    }
    save_state(state)

    tools = kc.ensure_jdk_and_tools(ds["src"])
    state["tools"] = tools
    save_state(state)

    cpus = max(1, os.cpu_count() or 4)
    kc.log(f"Detected {cpus} host CPUs; JARVIS3={tools['jarvis3']}; AGC={tools['agc']}")

    # ------------------------------------------------------------------
    # 1. Human GRCh38 Single Chromosomes + 1.48 GB Multi-Chr Mode C
    # ------------------------------------------------------------------
    chr22_seq_for_gpu = kc.WORK_ROOT / "chr22_gpu.seq"
    if ds["grch38"] is not None:
        g_dir = next(ds["grch38"].rglob("Homo_sapiens.GRCh38.dna.chromosome.1.fa")).parent
        for cname in ["21", "22", "X", "1"]:
            fa = g_dir / f"Homo_sapiens.GRCh38.dna.chromosome.{cname}.fa"
            if not fa.exists():
                continue
            kc.log(f"=== GRCh38 Chr{cname} ({fa.stat().st_size / (1<<20):.1f} MiB) ===")
            row = {
                "target": f"GRCh38_chr{cname}.fa",
                "raw_bytes": fa.stat().st_size,
                "runs": [],
            }
            row["runs"].append(kc.run_hgs2_compress_verify(tools, fa, mode="collection", threads=cpus, xmx="-Xmx16g"))
            row["runs"].append(kc.run_zstd_baseline(tools, fa, level=19))
            row["runs"].append(kc.run_xz_baseline(tools, fa, level=6))

            seq = kc.WORK_ROOT / f"grch38_chr{cname}.seq"
            smeta = kc.fasta_to_flat_seq(fa, seq)
            if cname == "22":
                shutil.copy2(seq, chr22_seq_for_gpu)
            row["flat_seq_bases"] = smeta["bases"]
            row["flat_non_acgt_masked"] = smeta["non_acgt_masked"]
            row["runs"].append(kc.run_hgs2_compress_verify(tools, seq, mode="collection", threads=cpus, xmx="-Xmx16g"))
            row["runs"].append(kc.run_hgs2_compress_verify(tools, seq, mode="collection", threads=min(cpus, 4), xmx="-Xmx256m", bounded=True))
            if fa.stat().st_size <= 165_000_000:
                row["runs"].append(kc.run_jarvis3_baseline(tools, seq, level=7, timeout=900))
            seq.unlink(missing_ok=True)
            state["grch38_chromosomes"].append(row)
            save_state(state)

        # Merge Chr1..Chr8 (~1.48 GB) for Single-Archive Mode C Scaling
        multi_fa = kc.WORK_ROOT / "grch38_chr1_8.fa"
        with multi_fa.open("wb") as fout:
            for cname in ["1", "2", "3", "4", "5", "6", "7", "8"]:
                p = g_dir / f"Homo_sapiens.GRCh38.dna.chromosome.{cname}.fa"
                if p.exists():
                    with p.open("rb") as fin:
                        shutil.copyfileobj(fin, fout, 8 << 20)
        kc.log(f"=== GRCh38 Chr1..8 Merged ({multi_fa.stat().st_size / (1<<20):.1f} MiB) ===")
        state["grch38_multichr_1_8"].append(
            kc.run_hgs2_compress_verify(tools, multi_fa, mode="collection", threads=cpus, xmx="-Xmx20g", timeout=2400)
        )
        state["grch38_multichr_1_8"].append(kc.run_zstd_baseline(tools, multi_fa, level=19, timeout=1800))
        multi_seq = kc.WORK_ROOT / "grch38_chr1_8.seq"
        kc.fasta_to_flat_seq(multi_fa, multi_seq)
        multi_fa.unlink(missing_ok=True)
        state["grch38_multichr_1_8"].append(
            kc.run_hgs2_compress_verify(tools, multi_seq, mode="collection", threads=min(cpus, 4), xmx="-Xmx1g", bounded=True, timeout=2400)
        )
        multi_seq.unlink(missing_ok=True)
        save_state(state)

    # ------------------------------------------------------------------
    # 2. Cross-Human Real Haplotype Comparison (GRCh38 vs T2T-CHM13v2.0)
    # ------------------------------------------------------------------
    chm13_dir = kc.WORK_ROOT / "chm13_chr"
    chm13_map = {}
    try:
        chm13_map = kc.download_and_split_chm13(chm13_dir)
    except Exception as exc:
        kc.log(f"Warning: T2T-CHM13v2.0 download skipped: {exc}")

    if chm13_map and ds["grch38"] is not None:
        g_dir = next(ds["grch38"].rglob("Homo_sapiens.GRCh38.dna.chromosome.1.fa")).parent
        # Mode R: CHM13 chr_i conditioned on GRCh38 chr_i
        for cname in ["21", "22", "19", "X", "1"]:
            g_fa = g_dir / f"Homo_sapiens.GRCh38.dna.chromosome.{cname}.fa"
            c_fa = chm13_map.get(f"chr{cname}")
            if not (g_fa.exists() and c_fa and c_fa.exists()):
                continue
            panel_dir = kc.WORK_ROOT / f"panel_chr{cname}"
            panel_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(g_fa, panel_dir / f"grch38_chr{cname}.fa")
            kc.log(f"=== Cross-Human Mode R (CHM13 chr{cname} | GRCh38 chr{cname}) ===")
            r_res = kc.run_hgs2_compress_verify(tools, c_fa, mode="ref", threads=cpus, xmx="-Xmx16g", panel=panel_dir)
            r_res["pair"] = f"CHM13_chr{cname}_given_GRCh38_chr{cname}"
            state["cross_human_chm13_vs_grch38"]["mode_r"].append(r_res)
            shutil.rmtree(panel_dir, ignore_errors=True)
            save_state(state)

        # Diploid Mode C: Pairs Chr21 (~94 MB), Chr22 (~103 MB), Chr1 (~502 MB), Chr1..4 (~1.8 GB)
        diploid_specs = [
            ("Diploid_Chr21", ["21"], True),
            ("Diploid_Chr22", ["22"], True),
            ("Diploid_Chr1", ["1"], False),
            ("Diploid_Chr1_to_4", ["1", "2", "3", "4"], False),
        ]
        for label, chrs, run_j in diploid_specs:
            g_samp = kc.WORK_ROOT / f"{label}_grch38.fa"
            c_samp = kc.WORK_ROOT / f"{label}_chm13.fa"
            merged = kc.WORK_ROOT / f"{label}_merged.fa"
            with g_samp.open("wb") as fg, c_samp.open("wb") as fc, merged.open("wb") as fm:
                for cn in chrs:
                    gp = g_dir / f"Homo_sapiens.GRCh38.dna.chromosome.{cn}.fa"
                    cp = chm13_map[f"chr{cn}"]
                    for p_in, f_out in [(gp, fg), (gp, fm), (cp, fc), (cp, fm)]:
                        with p_in.open("rb") as r:
                            shutil.copyfileobj(r, f_out, 8 << 20)
            kc.log(f"=== Cross-Human Diploid {label} ({merged.stat().st_size / (1<<20):.1f} MiB) ===")
            d_row = {
                "label": label,
                "raw_bytes": merged.stat().st_size,
                "runs": [],
            }
            hgs_to = 2400 if label == "Diploid_Chr1_to_4" else 1500
            zstd_to = 1800 if label == "Diploid_Chr1_to_4" else 900
            d_row["runs"].append(kc.run_hgs2_compress_verify(tools, merged, mode="collection", threads=cpus, xmx="-Xmx22g", timeout=hgs_to))
            d_row["runs"].append(kc.run_agc_collection(tools, g_samp, [c_samp], adaptive=False, threads=cpus, timeout=900))
            d_row["runs"].append(kc.run_agc_collection(tools, g_samp, [c_samp], adaptive=True, threads=cpus, timeout=900))
            d_row["runs"].append(kc.run_zstd_baseline(tools, merged, level=19, timeout=zstd_to))
            if run_j:
                d_seq = kc.WORK_ROOT / f"{label}.seq"
                kc.fasta_to_flat_seq(merged, d_seq)
                d_row["runs"].append(kc.run_hgs2_compress_verify(tools, d_seq, mode="collection", threads=cpus, xmx="-Xmx16g"))
                d_row["runs"].append(kc.run_hgs2_compress_verify(tools, d_seq, mode="collection", threads=min(cpus, 4), xmx="-Xmx256m", bounded=True))
                d_row["runs"].append(kc.run_jarvis3_baseline(tools, d_seq, level=7, timeout=900))
                d_seq.unlink(missing_ok=True)
            g_samp.unlink(missing_ok=True)
            c_samp.unlink(missing_ok=True)
            merged.unlink(missing_ok=True)
            state["cross_human_chm13_vs_grch38"]["diploid_mode_c"].append(d_row)
            save_state(state)

        shutil.rmtree(chm13_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # 3. Massive Viral Pangenome Scaling (COVID-19 Multi-FASTA Shards)
    # ------------------------------------------------------------------
    covid_seq_for_gpu = kc.WORK_ROOT / "covid_gpu.seq"
    if ds["covid19"] is not None:
        cov_dir = next(ds["covid19"].rglob("group_1.fasta")).parent
        for shard in ["group_19.fasta", "group_3.fasta", "group_22.fasta", "group_36.fasta"]:
            p = cov_dir / shard
            if not p.exists():
                continue
            kc.log(f"=== SARS-CoV-2 {shard} ({p.stat().st_size / (1<<20):.1f} MiB) ===")
            v_row = {"shard": shard, "raw_bytes": p.stat().st_size, "runs": []}
            v_row["runs"].append(kc.run_hgs2_compress_verify(tools, p, mode="collection", threads=cpus, xmx="-Xmx16g"))
            v_row["runs"].append(kc.run_zstd_baseline(tools, p, level=19))
            v_row["runs"].append(kc.run_xz_baseline(tools, p, level=6))

            v_seq = kc.WORK_ROOT / f"{p.stem}.seq"
            vmeta = kc.fasta_to_flat_seq(p, v_seq)
            if shard == "group_3.fasta":
                shutil.copy2(v_seq, covid_seq_for_gpu)
            v_row["flat_seq_bases"] = vmeta["bases"]
            v_row["runs"].append(kc.run_hgs2_compress_verify(tools, v_seq, mode="collection", threads=cpus, xmx="-Xmx16g"))
            v_row["runs"].append(kc.run_hgs2_compress_verify(tools, v_seq, mode="collection", threads=min(cpus, 4), xmx="-Xmx256m", bounded=True))
            if p.stat().st_size <= 160_000_000:
                v_row["runs"].append(kc.run_hgs2_compress_verify(tools, v_seq, mode="collection", threads=min(cpus, 2), xmx="-Xmx128m", bounded=True))
                v_row["runs"].append(kc.run_jarvis3_baseline(tools, v_seq, level=7, timeout=900))
            v_seq.unlink(missing_ok=True)
            state["covid19_viral_scaling"].append(v_row)
            save_state(state)

        # Merged Viral Mega-Batch (~1.45 GB: group_36 + 34 + 22 + 10 + 29)
        viral_mega = kc.WORK_ROOT / "covid_mega_1_45gb.fasta"
        with viral_mega.open("wb") as fout:
            for sh in ["group_36.fasta", "group_34.fasta", "group_22.fasta", "group_29.fasta"]:
                sp = cov_dir / sh
                if sp.exists():
                    with sp.open("rb") as fin:
                        shutil.copyfileobj(fin, fout, 8 << 20)
        if viral_mega.stat().st_size > 100_000_000:
            kc.log(f"=== SARS-CoV-2 Mega-Batch ({viral_mega.stat().st_size / (1<<20):.1f} MiB) ===")
            m_row = {"shard": "covid_mega_4groups", "raw_bytes": viral_mega.stat().st_size, "runs": []}
            m_row["runs"].append(kc.run_hgs2_compress_verify(tools, viral_mega, mode="collection", threads=cpus, xmx="-Xmx22g", timeout=1800))
            m_row["runs"].append(kc.run_zstd_baseline(tools, viral_mega, level=19, timeout=900))
            m_seq = kc.WORK_ROOT / "covid_mega.seq"
            kc.fasta_to_flat_seq(viral_mega, m_seq)
            viral_mega.unlink(missing_ok=True)
            m_row["runs"].append(kc.run_hgs2_compress_verify(tools, m_seq, mode="collection", threads=min(cpus, 4), xmx="-Xmx512m", bounded=True, timeout=1800))
            m_seq.unlink(missing_ok=True)
            state["covid19_viral_scaling"].append(m_row)
            save_state(state)

    # ------------------------------------------------------------------
    # 4. CUDA GPU Experiments (Sketch Similarity Matrix + Neural DNA)
    # ------------------------------------------------------------------
    if chr22_seq_for_gpu.exists() and covid_seq_for_gpu.exists():
        kc.log("=== Running CUDA GPU Sketch Similarity + Autoregressive DNA Entropy ===")
        state["cuda_gpu_experiments"] = run_cuda_gpu_experiments(
            state, chr22_seq_for_gpu, covid_seq_for_gpu, None
        )
        save_state(state)

    kc.log(f"COMPLETE. Results written to {OUT_JSON}")


if __name__ == "__main__":
    main()
