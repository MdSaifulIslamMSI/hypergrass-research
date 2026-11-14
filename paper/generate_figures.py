#!/usr/bin/env python3
"""
Generate all 7 publication-grade multi-panel figures for the HGS2 / EIMDL
end-to-end research paper directly from committed JSON benchmark artifacts
and verified experimental records across Generation 1 (HyperGraSS v1),
Generation 2 (RC-MDL v2), and Generation 3 (HGS2 / EIMDL).
Zero synthetic/interpolated numbers; zero overlapping labels or legends.
"""

import json
import os
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

BASE_DIR = r"C:\Users\mdsai\.gemini\antigravity\scratch\hypergrass"
FIG_DIR = os.path.join(BASE_DIR, "paper", "figures")
ARTIFACT_DIR = r"C:\Users\mdsai\.gemini\antigravity\brain\8fd6ada9-0a0f-403f-95e4-ee038b5c4fa3"

os.makedirs(FIG_DIR, exist_ok=True)

# Publication typography and styling
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif", "Computer Modern Roman"],
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 7.8,
    "figure.titlesize": 12.0,
    "figure.titleweight": "bold",
    "axes.grid": True,
    "grid.alpha": 0.28,
    "grid.linestyle": "--",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

PALETTE = {
    "eimdl": "#0B4F6C",       # Deep Navy Teal
    "eimdl_alt": "#01BAEF",   # Bright Cyan
    "rcmdl": "#20BF55",       # Emerald Green
    "hypergrass1": "#5B8C5A", # Sage Green
    "jarvis_l7": "#E07A5F",   # Terra Cotta
    "jarvis_l31": "#D62828",  # Crimson Red
    "agc_single": "#F4A261",  # Warm Amber
    "agc_oracle": "#8338EC",  # Royal Violet
    "grass": "#6C757D",       # Slate Gray
    "geco": "#B56576",        # Mauve
    "xz": "#9A8C98",          # Muted Gray-Purple
}


def load_json(rel_path):
    with open(os.path.join(BASE_DIR, rel_path), "r", encoding="utf-8") as f:
        return json.load(f)


def save_fig(fig, name):
    pdf_path = os.path.join(FIG_DIR, f"{name}.pdf")
    png_path = os.path.join(FIG_DIR, f"{name}.png")
    fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    if os.path.exists(ARTIFACT_DIR):
        shutil.copy2(png_path, os.path.join(ARTIFACT_DIR, f"{name}.png"))
    print(f"[OK] Saved {pdf_path} and {png_path}")


def make_fig1_regime_overview():
    mode_s_syn = load_json("benchmarks/fusion/results/mode_s_and_synthetic.json")
    mode_c_col = load_json("benchmarks/fusion/results/mode_c_collections.json")
    amdahl_mem = load_json("benchmarks/fusion/results/amdahl_and_memory.json")
    loo_ecoli = load_json("benchmarks/fusion/results/loo_ecoli.json")
    loo_mtb = load_json("benchmarks/fusion/results/loo_mtb.json")
    loo_saur = load_json("benchmarks/fusion/results/loo_saureus.json")
    sal_trio = load_json("benchmarks/independent/runs/salmonella/trio-1/report.json")
    kleb_trio = load_json("benchmarks/independent/runs/klebsiella/trio-1/report.json")
    sal_modes = load_json("benchmarks/independent/runs/salmonella/modes-subset-1/report.json")
    kleb_modes = load_json("benchmarks/independent/runs/klebsiella/modes-subset-1/report.json")
    ablation = load_json("benchmarks/fusion/results/ablations_real_genomes.json")

    fig = plt.figure(figsize=(13.6, 5.2))
    gs = GridSpec(1, 2, width_ratios=[1.35, 1.05], wspace=0.36)

    # Panel A: 3-Regime Compression Rate (bpb, log scale) across all 5 bacterial species
    ax1 = fig.add_subplot(gs[0])
    species_labels = [
        "E. coli\n(24g, 120.1 Mb)",
        "M. tuberculosis\n(24g, 105.7 Mb)",
        "S. aureus\n(25g, 71.0 Mb)",
        "S. enterica\n(25g, 120.5 Mb)",
        "K. pneumoniae\n(23g, 122.7 Mb)",
    ]

    mode_s_map = {r["label"]: r["bpb"] for r in mode_s_syn["real_mode_s"]}
    mode_s_bpb = [
        mode_s_map["ecoli_00_seq"],
        mode_s_map["mtb_00_seq"],
        mode_s_map["saureus_00_seq"],
        sal_modes["summaries"]["S/eimdl"]["total_bpb"],
        kleb_modes["summaries"]["S/eimdl"]["total_bpb"],
    ]

    mode_c_map = {r["label"]: r["bpb"] for r in mode_c_col["collections"]}
    ecoli_t1_bpb = amdahl_mem["amdahl_scaling"][0]["bpb"]
    sal_c_bpb = [r for r in sal_trio["results"] if r["profile"] == "eimdl-flat"][0]["bits_per_base"]
    kleb_c_bpb = [r for r in kleb_trio["results"] if r["profile"] == "eimdl-flat"][0]["bits_per_base"]
    mode_c_bpb = [
        ecoli_t1_bpb,
        mode_c_map["mtb_24_multifasta"],
        mode_c_map["saureus_25_multifasta"],
        sal_c_bpb,
        kleb_c_bpb,
    ]

    mode_r_bpb = [
        loo_ecoli["total_bpb"],
        loo_mtb["total_bpb"],
        loo_saur["total_bpb"],
        sal_modes["summaries"]["R/eimdl"]["total_bpb"],
        kleb_modes["summaries"]["R/eimdl"]["total_bpb"],
    ]

    x = np.arange(len(species_labels))
    w = 0.26
    b1 = ax1.bar(x - w, mode_s_bpb, width=w, color=PALETTE["jarvis_l7"], edgecolor="black", linewidth=0.6,
                 label="Mode S (Reference-Free Single Genome)")
    b2 = ax1.bar(x, mode_c_bpb, width=w, color=PALETTE["eimdl_alt"], edgecolor="black", linewidth=0.6,
                 label="Mode C (Self-Contained Pangenome Archive)")
    b3 = ax1.bar(x + w, mode_r_bpb, width=w, color=PALETTE["eimdl"], edgecolor="black", linewidth=0.6,
                 label="Mode R (Reference-Conditioned LOO)")

    ax1.axhline(2.0, color="#D62828", linestyle=":", linewidth=1.3, label="2.000 bpb Naive 2-Bit Ceiling")
    ax1.set_yscale("log")
    # Generous upper headroom (up to 38) so the legend sits completely above the 2.000 line
    ax1.set_ylim(0.002, 38.0)
    ax1.set_xticks(x)
    ax1.set_xticklabels(species_labels)
    ax1.set_ylabel("Compression Rate (bits / base, log scale)")
    ax1.set_title("(A) Multi-Regime Description Length Across 5 Bacterial Pangenomes", pad=10)

    for bar in b1:
        h = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2.0, h * 1.14, f"{h:.2f}",
                 ha="center", va="bottom", fontsize=7.0)
    for bar in b2:
        h = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2.0, h * 1.14, f"{h:.3f}",
                 ha="center", va="bottom", fontsize=7.0)
    for bar in b3:
        h = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2.0 + 0.03, h * 1.16, f"{h:.4f}",
                 ha="center", va="bottom", fontsize=6.8, rotation=28, fontweight="bold")

    ax1.legend(loc="upper right", ncol=1, frameon=True, facecolor="white", framealpha=0.95, fontsize=7.5)

    # Panel B: Bit-Accounting Progression on E. coli 00.seq (4,856,080 bp)
    ax2 = fig.add_subplot(gs[1])
    ecoli_00_mode_s = [r for r in mode_s_syn["real_mode_s"] if r["label"] == "ecoli_00_seq"][0]["archive_bytes"]
    ab_distant = {r["ablation"]: r["archive_bytes"] for r in ablation["ecoli_00_vs_distant_donor20_ablations"]}
    ab_loo = {r["ablation"]: r["archive_bytes"] for r in ablation["ecoli_00_loo_ablations"]}

    stages = [
        "Raw ASCII\n(00.seq)",
        "2-Bit Packed\nCeiling",
        "Mode S\n(FastCRM)",
        "Mode R (1 Donor)\nExact Only",
        "Mode R (1 Donor)\nFull EIMDL",
        "Mode R (23 Donors)\nExact Only",
        "Mode R (23 Donors)\nFull EIMDL",
    ]
    sizes_kb = [
        4856080 / 1000.0,
        (4856080 / 4.0) / 1000.0,
        ecoli_00_mode_s / 1000.0,
        ab_distant["exact_matches_only_no_approx"] / 1000.0,
        ab_distant["full_eimdl"] / 1000.0,
        ab_loo["exact_matches_only_no_approx"] / 1000.0,
        ab_loo["full_eimdl"] / 1000.0,
    ]
    colors = ["#ADB5BD", "#6C757D", PALETTE["jarvis_l7"], PALETTE["agc_single"],
              PALETTE["eimdl_alt"], PALETTE["rcmdl"], PALETTE["eimdl"]]

    y_pos = np.arange(len(stages))
    bars2 = ax2.barh(y_pos, sizes_kb, color=colors, edgecolor="black", linewidth=0.6, height=0.62)
    ax2.set_xscale("log")
    ax2.set_xlim(20, 32000)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(stages)
    ax2.invert_yaxis()
    ax2.set_xlabel("Archive Size on E. coli 00.seq (KB, log scale)")
    ax2.set_title("(B) E. coli 00.seq (4.86 Mbp) Information Decomposition", pad=10)

    for bar, kb in zip(bars2, sizes_kb):
        bpb_val = (kb * 1000.0 * 8.0) / 4856080.0
        ax2.text(kb * 1.14, bar.get_y() + bar.get_height() / 2.0,
                 f"{kb:,.1f} KB ({bpb_val:.4f} bpb)",
                 va="center", ha="left", fontsize=7.5, fontweight="bold" if kb < 50 else "normal")

    save_fig(fig, "fig1_regime_overview")


def make_fig2_mode_r_loo_73genomes():
    loo_ecoli = load_json("benchmarks/fusion/results/loo_ecoli.json")
    loo_mtb = load_json("benchmarks/fusion/results/loo_mtb.json")
    loo_saur = load_json("benchmarks/fusion/results/loo_saureus.json")
    mtb_agc = load_json("prototype/rcmdl/work/agc-sweep-mtb-full/report.json")
    saur_oracle = load_json("prototype/rcmdl/work/agc-verify-saureus-oracle/results.json")
    cl_ecoli = load_json("benchmarks/fusion/results/cluster_loo_ecoli.json")
    cl_mtb = load_json("benchmarks/fusion/results/cluster_loo_mtb.json")
    cl_saur = load_json("benchmarks/fusion/results/cluster_loo_saureus.json")

    fig = plt.figure(figsize=(14.2, 9.4))
    gs = GridSpec(2, 2, hspace=0.38, wspace=0.25)

    # Panel A: E. coli 24 genomes per-target Mode R LOO
    ax_a = fig.add_subplot(gs[0, 0])
    ecoli_recs = loo_ecoli["records"]
    idx_ec = np.arange(len(ecoli_recs))
    ec_hgs2 = [r["archive_bytes"] / 1000.0 for r in ecoli_recs]
    ec_rcmdl = [r["rcmdl_bytes"] / 1000.0 for r in ecoli_recs]
    ec_agc = [r["agc_bytes"] / 1000.0 for r in ecoli_recs]
    ec_hgs2_kb = sum(ec_hgs2)
    ec_rcmdl_kb = sum(ec_rcmdl)
    ec_agc_kb = sum(ec_agc)
    ec_wins = loo_ecoli.get("eimdl_wins_vs_agc", len(ecoli_recs))
    ec_fac = ec_agc_kb / ec_hgs2_kb

    ax_a.plot(idx_ec, ec_agc, "s--", color=PALETTE["agc_single"], markersize=4, linewidth=1.2,
              label=f"AGC 3.2.5 Best Single (-s 200000): {ec_agc_kb:,.1f} KB")
    ax_a.plot(idx_ec, ec_rcmdl, "^-", color=PALETTE["rcmdl"], markersize=4, linewidth=1.2,
              label=f"RC-MDL v2 (Unaligned): {ec_rcmdl_kb:,.1f} KB")
    ax_a.plot(idx_ec, ec_hgs2, "o-", color=PALETTE["eimdl"], markersize=4.5, linewidth=1.6,
              label=f"EIMDL / HGS2 Mode R: {ec_hgs2_kb:,.1f} KB ({ec_wins}/{len(ecoli_recs)} wins)")
    ax_a.set_yscale("log")
    # Peak is 767 KB; setting ylim up to 11,000 KB puts the upper band completely above all curves
    ax_a.set_ylim(0.5, 11000)
    ax_a.set_xticks(idx_ec[::2])
    ax_a.set_xticklabels([r["target"].replace(".seq", "") for r in ecoli_recs[::2]])
    ax_a.set_xlabel("E. coli Target Genome Index (00..23)")
    ax_a.set_ylabel("Marginal Archive Size (KB, log scale)")
    ax_a.set_title(f"(A) E. coli Leave-One-Out (24 Genomes, 120.11 Mbp): {ec_fac:.2f}x Smaller than AGC")
    ax_a.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.4)

    # Panel B: M. tuberculosis 24 genomes per-target Mode R LOO (with AGC 13-Config Oracle)
    ax_b = fig.add_subplot(gs[0, 1])
    mtb_recs = loo_mtb["records"]
    idx_mtb = np.arange(len(mtb_recs))
    mtb_hgs2 = [r["archive_bytes"] for r in mtb_recs]
    mtb_rcmdl = [r["rcmdl_bytes"] for r in mtb_recs]
    mtb_oracle_map = {r["target"]: r["best_bytes"] for r in mtb_agc["targets_detail"]}
    mtb_agc_oracle = [mtb_oracle_map[r["target"]] for r in mtb_recs]
    mtb_hgs2_sum = loo_mtb["total_archive_bytes"]
    mtb_rcmdl_sum = sum(mtb_rcmdl)
    mtb_oracle_sum = sum(mtb_agc_oracle)
    mtb_o_fac = mtb_oracle_sum / mtb_hgs2_sum

    ax_b.plot(idx_mtb, mtb_agc_oracle, "d-.", color=PALETTE["agc_oracle"], markersize=4, linewidth=1.2,
              label=f"AGC 13-Config Per-Target Oracle: {mtb_oracle_sum:,} B")
    ax_b.plot(idx_mtb, mtb_rcmdl, "^-", color=PALETTE["rcmdl"], markersize=4, linewidth=1.1,
              label=f"RC-MDL v2 (Unaligned): {mtb_rcmdl_sum:,} B")
    ax_b.plot(idx_mtb, mtb_hgs2, "o-", color=PALETTE["eimdl"], markersize=4.5, linewidth=1.6,
              label=f"EIMDL / HGS2 Mode R: {mtb_hgs2_sum:,} B ({mtb_o_fac:.2f}x < Oracle)")
    ax_b.axhline(292170 / 24.0, color=PALETTE["agc_single"], linestyle="--", linewidth=1.1,
                 label="AGC Best Single (-s 100000) Mean: 12,174 B/genome")
    ax_b.set_yscale("log")
    # Peak is 21,796 B; setting ylim up to 180,000 B keeps the top region completely clear
    ax_b.set_ylim(180, 180000)
    ax_b.set_xticks(idx_mtb[::2])
    ax_b.set_xticklabels([r["target"].replace(".seq", "") for r in mtb_recs[::2]])
    ax_b.set_xlabel("M. tuberculosis Target Genome Index (00..23)")
    ax_b.set_ylabel("Marginal Archive Size (Bytes, log scale)")
    ax_b.set_title(f"(B) M. tuberculosis LOO (24 Genomes, 105.75 Mbp): {loo_mtb['total_bpb']:.5f} bpb Total")
    ax_b.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.4)

    # Panel C: S. aureus 25 genomes per-target Mode R LOO
    ax_c = fig.add_subplot(gs[1, 0])
    saur_recs = loo_saur["records"]
    idx_sa = np.arange(len(saur_recs))
    sa_hgs2 = [r["archive_bytes"] / 1000.0 for r in saur_recs]
    sa_rcmdl = [r["rcmdl_bytes"] / 1000.0 for r in saur_recs]
    sa_oracle_map = {r["target"].replace(".fa", ".seq"): r["marginal_bytes"] for r in saur_oracle}
    sa_agc_oracle = [sa_oracle_map[r["target"]] / 1000.0 for r in saur_recs]
    sa_hgs2_kb = loo_saur["total_archive_bytes"] / 1000.0
    sa_rcmdl_kb = sum(sa_rcmdl)
    sa_oracle_kb = sum(sa_agc_oracle)
    sa_o_fac = sa_oracle_kb / sa_hgs2_kb
    sa_o_wins = sum(1 for a, o in zip(sa_hgs2, sa_agc_oracle) if a < o)

    ax_c.axhline((2979078 / 25.0) / 1000.0, color=PALETTE["agc_single"], linestyle="--", linewidth=1.1,
                 label="AGC Best Single (-s 50000) Mean: 119.2 KB/genome (2,979.1 KB total)")
    ax_c.plot(idx_sa, sa_agc_oracle, "d-.", color=PALETTE["agc_oracle"], markersize=4, linewidth=1.2,
              label=f"AGC 13-Config Per-Target Oracle: {sa_oracle_kb:,.1f} KB")
    ax_c.plot(idx_sa, sa_rcmdl, "^-", color=PALETTE["rcmdl"], markersize=4, linewidth=1.1,
              label=f"RC-MDL v2: {sa_rcmdl_kb:,.1f} KB")
    ax_c.plot(idx_sa, sa_hgs2, "o-", color=PALETTE["eimdl"], markersize=4.5, linewidth=1.6,
              label=f"EIMDL / HGS2 Mode R: {sa_hgs2_kb:,.1f} KB ({sa_o_wins}/25 vs Oracle, 25/25 vs Single)")
    ax_c.set_yscale("log")
    # Peak is 274 KB; setting ylim up to 6,000 KB keeps the top band completely clear of all lines
    ax_c.set_ylim(0.08, 6000)
    ax_c.set_xticks(idx_sa[::2])
    ax_c.set_xticklabels([r["target"].replace(".seq", "") for r in saur_recs[::2]])
    ax_c.set_xlabel("S. aureus Target Genome Index (00..24)")
    ax_c.set_ylabel("Marginal Archive Size (KB, log scale)")
    ax_c.set_title(f"(C) S. aureus LOO (25 Genomes, 71.04 Mbp): {sa_o_fac:.2f}x Smaller than AGC Oracle")
    ax_c.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.3)

    # Panel D: Leave-Cluster-Out (Unseen Lineage) vs Standard LOO vs AGC
    ax_d = fig.add_subplot(gs[1, 1])
    sp_names = ["E. coli\n(22 clusters, J=0.98)", "M. tuberculosis\n(13 clusters, J=0.99)", "S. aureus\n(14 clusters, J=0.98)"]
    loo_bpb = [
        loo_ecoli["total_bpb"],
        loo_mtb["total_bpb"],
        loo_saur["total_bpb"],
    ]
    cl_bpb = [
        cl_ecoli["total_bpb"],
        cl_mtb["total_bpb"],
        cl_saur["total_bpb"],
    ]
    rcmdl_cl_bpb = [
        cl_ecoli["rcmdl_cluster_bytes"] * 8.0 / cl_ecoli["total_bases"],
        cl_mtb["rcmdl_cluster_bytes"] * 8.0 / cl_mtb["total_bases"],
        cl_saur["rcmdl_cluster_bytes"] * 8.0 / cl_saur["total_bases"],
    ]
    agc_cl_bpb = [
        9203994 * 8.0 / cl_ecoli["total_bases"],
        1188909 * 8.0 / cl_mtb["total_bases"],
        4132771 * 8.0 / cl_saur["total_bases"],
    ]

    x_d = np.arange(3)
    w_d = 0.20
    b_d1 = ax_d.bar(x_d - 1.5 * w_d, loo_bpb, width=w_d, color=PALETTE["eimdl"], edgecolor="black", linewidth=0.6,
                    label="EIMDL Standard LOO")
    b_d2 = ax_d.bar(x_d - 0.5 * w_d, cl_bpb, width=w_d, color=PALETTE["eimdl_alt"], edgecolor="black", linewidth=0.6,
                    label="EIMDL Leave-Cluster-Out")
    b_d3 = ax_d.bar(x_d + 0.5 * w_d, rcmdl_cl_bpb, width=w_d, color=PALETTE["rcmdl"], edgecolor="black", linewidth=0.6,
                    label="RC-MDL v2 Leave-Cluster-Out")
    b_d4 = ax_d.bar(x_d + 1.5 * w_d, agc_cl_bpb, width=w_d, color=PALETTE["agc_single"], edgecolor="black", linewidth=0.6,
                    label="AGC 3.2.5 (-a) Leave-Cluster-Out")

    for bars in [b_d1, b_d2, b_d3, b_d4]:
        for bar in bars:
            h = bar.get_height()
            ax_d.text(bar.get_x() + bar.get_width() / 2.0, h * 1.12, f"{h:.3f}",
                      ha="center", va="bottom", fontsize=6.6, rotation=30)

    ax_d.set_yscale("log")
    # Highest bar is 0.612 bpb; ylim up to 14.0 gives clean whitespace for the legend at upper left
    ax_d.set_ylim(0.002, 14.0)
    ax_d.set_xticks(x_d)
    ax_d.set_xticklabels(sp_names)
    ax_d.set_ylabel("Compression Rate (bits / base, log scale)")
    ax_d.set_title("(D) Unseen-Lineage Generalization: Leave-Cluster-Out Controls")
    ax_d.legend(loc="upper left", ncol=2, frameon=True, facecolor="white", framealpha=0.95, fontsize=7.2)

    save_fig(fig, "fig2_mode_r_loo_73genomes")


def make_fig3_mode_c_pangenomes():
    mode_c_col = load_json("benchmarks/fusion/results/mode_c_collections.json")
    amdahl_mem = load_json("benchmarks/fusion/results/amdahl_and_memory.json")
    sal_trio = load_json("benchmarks/independent/runs/salmonella/trio-1/report.json")
    kleb_trio = load_json("benchmarks/independent/runs/klebsiella/trio-1/report.json")

    fig = plt.figure(figsize=(13.8, 5.3))
    gs = GridSpec(1, 2, width_ratios=[1.18, 1.18], wspace=0.32)

    # Panel A: Canonical 24-Genome E. coli Pangenome (120,108,223 bp) across all baselines
    ax1 = fig.add_subplot(gs[0])
    ec_t1 = amdahl_mem["amdahl_scaling"][0]["archive_bytes"]
    ec_def = [r for r in mode_c_col["collections"] if r["label"] == "ecoli_24_flat_seq"][0]["archive_bytes"]

    codecs = [
        "EIMDL / HGS2 (1T Exact)",
        "EIMDL / HGS2 (Default 12T)",
        "HyperGraSS v1 (BWT+CRM)",
        "JARVIS3 v3.7 (-l 7)",
        "RC-MDL v2 (Prim Tree)",
        "xz -9 (Flat Concatenated)",
        "xz -9 (Multi-FASTA)",
        "GeCo3 (-l 6)",
        "AGC 3.2.5 (-a -t 1)",
        "GeCo3 (-l 3)",
        "Original GraSS (-e2)",
    ]
    bytes_list = [
        ec_t1,
        ec_def,
        3810409,
        3896361,
        3972975,
        5663640,
        8394900,
        8409358,
        8540341,
        14389043,
        24777592,
    ]
    mb_list = [b / 1e6 for b in bytes_list]
    bar_colors = [
        PALETTE["eimdl"],
        PALETTE["eimdl_alt"],
        PALETTE["hypergrass1"],
        PALETTE["jarvis_l7"],
        PALETTE["rcmdl"],
        PALETTE["xz"],
        PALETTE["xz"],
        PALETTE["geco"],
        PALETTE["agc_single"],
        PALETTE["geco"],
        PALETTE["grass"],
    ]

    y = np.arange(len(codecs))
    bars = ax1.barh(y, mb_list, color=bar_colors, edgecolor="black", linewidth=0.6, height=0.68)
    ax1.set_yticks(y)
    ax1.set_yticklabels(codecs)
    ax1.invert_yaxis()
    ax1.set_xlim(0, 34.5)
    ax1.set_xlabel("Archive Size on 24-Genome E. coli Collection (MB, 120.11 Mbp raw)")
    ax1.set_title("(A) 24-Genome E. coli Pangenome Collection Benchmark", pad=10)

    for bar, b_val in zip(bars, bytes_list):
        bpb_val = (b_val * 8.0) / 120108223.0
        ratio_val = 120108223.0 / b_val
        ax1.text(bar.get_width() + 0.35, bar.get_y() + bar.get_height() / 2.0,
                 f"{b_val / 1e6:.2f} MB ({bpb_val:.3f} bpb, {ratio_val:.1f}x)",
                 va="center", ha="left", fontsize=7.4,
                 fontweight="bold" if b_val <= ec_def else "normal")

    # Panel B: 5-Species Mode C Pangenome Compression (bits/base)
    ax2 = fig.add_subplot(gs[1])
    mc_map = {r["label"]: r for r in mode_c_col["collections"]}
    sal_map = {r["profile"]: r["bits_per_base"] for r in sal_trio["results"]}
    kleb_map = {r["profile"]: r["bits_per_base"] for r in kleb_trio["results"]}

    species_5 = [
        "E. coli\n(24g, 120.1 Mb)",
        "M. tuberculosis\n(24g, 105.7 Mb)",
        "S. aureus\n(25g, 71.0 Mb)",
        "S. enterica\n(25g, 120.5 Mb)",
        "K. pneumoniae\n(23g, 122.7 Mb)",
    ]
    eimdl_5 = [
        amdahl_mem["amdahl_scaling"][0]["bpb"],
        mc_map["mtb_24_multifasta"]["bpb"],
        mc_map["saureus_25_multifasta"]["bpb"],
        sal_map["eimdl-flat"],
        kleb_map["eimdl-flat"],
    ]
    comp_secondary = [
        (3896361 * 8.0) / 120108223.0,  # JARVIS3 l7 on E. coli
        (1148486 * 8.0) / 105745015.0,  # RC-MDL v2 on Mtb
        (1253650 * 8.0) / 71043114.0,   # RC-MDL v2 on S. aureus
        sal_map["jarvis3-l7"],          # JARVIS3 l7 on Salmonella
        kleb_map["jarvis3-l7"],         # JARVIS3 l7 on Klebsiella
    ]
    comp_grass = [
        (24777592 * 8.0) / 120108223.0, # Original GraSS on E. coli
        (3318508 * 8.0) / 105745015.0,  # AGC 3.2.5 on Mtb collection
        (5637705 * 8.0) / 71043114.0,   # AGC 3.2.5 on S. aureus collection
        sal_map["grass-e2"],            # Original GraSS on Salmonella
        kleb_map["grass-e2"],           # Original GraSS on Klebsiella
    ]

    x5 = np.arange(len(species_5))
    w5 = 0.27
    b_e = ax2.bar(x5 - w5, eimdl_5, width=w5, color=PALETTE["eimdl"], edgecolor="black", linewidth=0.6,
                  label="EIMDL / HGS2 Mode C")
    b_s = ax2.bar(x5, comp_secondary, width=w5, color=PALETTE["jarvis_l7"], edgecolor="black", linewidth=0.6,
                  label="Best Runner-Up (JARVIS3 -l 7 / RC-MDL v2)")
    b_g = ax2.bar(x5 + w5, comp_grass, width=w5, color=PALETTE["grass"], edgecolor="black", linewidth=0.6,
                  label="Prior Baseline (Original GraSS -e2 / AGC 3.2.5)")

    ax2.set_yscale("log")
    # Highest bar is 1.650 bpb; setting ylim up to 11.0 keeps the top region completely clear for the legend
    ax2.set_ylim(0.05, 11.0)
    ax2.set_xticks(x5)
    ax2.set_xticklabels(species_5)
    ax2.set_ylabel("Collection Compression Rate (bits / base, log scale)")
    ax2.set_title("(B) Multi-Species Mode C Pangenome Collections (5 Species)", pad=10)

    # Vertical bar-top labels so adjacent bars never touch horizontally
    for bar in b_e:
        h = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2.0, h * 1.10, f"{h:.3f}",
                 ha="center", va="bottom", fontsize=6.8, rotation=90, fontweight="bold", color=PALETTE["eimdl"])
    for bar in b_s:
        h = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2.0, h * 1.10, f"{h:.3f}",
                 ha="center", va="bottom", fontsize=6.8, rotation=90)
    for bar in b_g:
        h = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2.0, h * 1.10, f"{h:.3f}",
                 ha="center", va="bottom", fontsize=6.8, rotation=0)

    ax2.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.4)
    save_fig(fig, "fig3_mode_c_pangenomes")


def make_fig4_amdahl_and_bounded_memory():
    amdahl_mem = load_json("benchmarks/fusion/results/amdahl_and_memory.json")
    mode_c_col = load_json("benchmarks/fusion/results/mode_c_collections.json")

    fig = plt.figure(figsize=(13.6, 5.1))
    gs = GridSpec(1, 2, width_ratios=[1.1, 1.15], wspace=0.30)

    # Panel A: Amdahl Multi-Core Scaling (1, 2, 4, 8, 12 threads)
    ax1 = fig.add_subplot(gs[0])
    scaling = amdahl_mem["amdahl_scaling"]
    threads = [r["threads"] for r in scaling]
    enc_s = [r["compress_inproc_seconds"] for r in scaling]
    dec_s = [r["decompress_inproc_seconds"] for r in scaling]
    tot_s = [e + d for e, d in zip(enc_s, dec_s)]

    ax1.plot(threads, tot_s, "s-", color=PALETTE["eimdl"], linewidth=2.0, markersize=6,
             label=f"Total Round-Trip ({tot_s[0]:.2f}s -> {tot_s[-1]:.2f}s, {tot_s[0]/tot_s[-1]:.2f}x)")
    ax1.plot(threads, enc_s, "o--", color=PALETTE["jarvis_l7"], linewidth=1.6, markersize=5,
             label=f"Encode Time ({enc_s[0]:.2f}s -> {enc_s[-1]:.2f}s, {enc_s[0]/enc_s[-1]:.2f}x)")
    ax1.plot(threads, dec_s, "^-.", color=PALETTE["rcmdl"], linewidth=1.6, markersize=5,
             label=f"Decode Time ({dec_s[0]:.2f}s -> {min(dec_s):.2f}s best / {dec_s[-1]:.2f}s 12T)")

    for t, tot in zip(threads, tot_s):
        ax1.annotate(f"{tot:.1f}s", (t, tot), textcoords="offset points", xytext=(0, 7),
                     ha="center", fontsize=7.5, fontweight="bold")

    ax1.set_xticks(threads)
    ax1.set_xlabel("Worker Thread Count (--threads T)")
    ax1.set_ylabel("In-Process Execution Time (seconds)")
    ax1.set_ylim(0, max(tot_s) * 1.32)
    ax1.set_title("(A) Laws V & VII: Amdahl Multi-Core Scaling on 120.11 Mbp E. coli", pad=10)
    ax1.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.4)

    # Panel B: Bounded-Heap JVM Execution (-Xmx128m .. -Xmx256m, Default)
    ax2 = fig.add_subplot(gs[1])
    mem_recs = amdahl_mem["bounded_memory"]
    default_rec = [r for r in mode_c_col["collections"] if r["label"] == "ecoli_24_flat_seq"][0]

    cap_cfg = {
        "-Xmx256m": ("-Xmx256m\n(s=1, 1/2)", 120.108223 * (0.25 + 0.5)),
        "-Xmx224m": ("-Xmx224m\n(s=1, 1/2)", 120.108223 * (0.25 + 0.5)),
        "-Xmx192m": ("-Xmx192m\n(s=2, 1/4)", 120.108223 * (0.25 + 0.25)),
        "-Xmx160m": ("-Xmx160m\n(s=2, L-VII)", 120.108223 * (0.25 + 0.25)),
        "-Xmx128m": ("-Xmx128m\n(s=2, L-VII)", 120.108223 * (0.25 + 0.25)),
    }
    labels_m = ["Default\n(s=0, full)"] + [cap_cfg.get(r["heap_cap"], (r["heap_cap"], 60.1))[0] for r in mem_recs]
    bank_plus_csr_mb = [120.108223 * (0.25 + 1.0)] + [cap_cfg.get(r["heap_cap"], (r["heap_cap"], 60.1))[1] for r in mem_recs]
    archive_mb = [default_rec["archive_bytes"] / 1e6] + [r["archive_bytes"] / 1e6 for r in mem_recs]

    x_m = np.arange(len(labels_m))
    w_m = 0.36
    bars_mem = ax2.bar(x_m - w_m / 2, bank_plus_csr_mb, width=w_m, color=PALETTE["eimdl_alt"],
                       edgecolor="black", linewidth=0.6, label="Stage-1 DNA Bank + CSR Index (MB)")
    ax2_r = ax2.twinx()
    ax2_r.spines["top"].set_visible(False)
    bars_arc = ax2_r.bar(x_m + w_m / 2, archive_mb, width=w_m, color=PALETTE["eimdl"],
                         edgecolor="black", linewidth=0.6, label="Resulting Mode C Archive Size (MB)")
    ax2_r.axhline(3.810409, color=PALETTE["hypergrass1"], linestyle="--", linewidth=1.2,
                  label="HyperGraSS v1 (3.810 MB @ 1,452 MB RSS)")
    ax2_r.axhline(3.896361, color=PALETTE["jarvis_l7"], linestyle=":", linewidth=1.3,
                  label="JARVIS3 -l 7 (3.896 MB @ 549 MB RSS)")

    ax2.set_xticks(x_m)
    ax2.set_xticklabels(labels_m, fontsize=7.0)
    ax2.set_ylabel("Stage-1 Core Data Structure Footprint (MB)", color=PALETTE["eimdl_alt"])
    ax2_r.set_ylabel("Mode C Compressed Archive Size (MB)", color=PALETTE["eimdl"])
    ax2.set_ylim(0, 225)
    ax2_r.set_ylim(3.10, 4.55)
    ax2.set_title("(B) Laws V & VII: Pure-JVM Bounded Heap (-Xmx128m to -Xmx256m)", pad=10)

    for bar, v in zip(bars_mem, bank_plus_csr_mb):
        ax2.text(bar.get_x() + bar.get_width() / 2.0, v + 3.5, f"{v:.1f}\nMB",
                 ha="center", va="bottom", fontsize=6.6, fontweight="bold")
    for bar, v in zip(bars_arc, archive_mb):
        ax2_r.text(bar.get_x() + bar.get_width() / 2.0 + 0.03, v + 0.02, f"{v:.3f}\nMB",
                   ha="center", va="bottom", fontsize=6.5, fontweight="bold", color=PALETTE["eimdl"])

    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_r.get_legend_handles_labels()
    ax2_r.legend(lines1 + lines2, labels1 + labels2, loc="upper right", frameon=True,
                 facecolor="white", framealpha=0.95, fontsize=6.9)

    save_fig(fig, "fig4_amdahl_and_bounded_memory")


def make_fig5_laws_ablations_and_synthetic():
    ablation = load_json("benchmarks/fusion/results/ablations_real_genomes.json")
    mode_s_syn = load_json("benchmarks/fusion/results/mode_s_and_synthetic.json")
    amdahl_mem = load_json("benchmarks/fusion/results/amdahl_and_memory.json")
    loo_ecoli = load_json("benchmarks/fusion/results/loo_ecoli.json")

    # Spacious 2x2 grid so every panel has plenty of horizontal and vertical breathing room
    fig = plt.figure(figsize=(14.2, 9.6))
    gs = GridSpec(2, 2, hspace=0.40, wspace=0.34)

    # Panel A [Top-Left]: Law-by-Law Ablation Trajectory on E. coli 00.seq
    ax1 = fig.add_subplot(gs[0, 0])
    dist_map = {r["ablation"]: r for r in ablation["ecoli_00_vs_distant_donor20_ablations"]}
    loo_map = {r["ablation"]: r for r in ablation["ecoli_00_loo_ablations"]}

    order_keys = [
        "exact_matches_only_no_approx",
        "no_v4_no_wobble_no_splice",
        "no_kimura_v4_law_I",
        "no_codon_wobble_law_II",
        "no_indel_splice_law_III",
        "full_eimdl",
    ]
    step_labels = [
        "Exact Only\n(No Approx)",
        "No Laws\nI,II,IIIc",
        "No Law I\n(Flat SNP)",
        "No Law II\n(Flat Gap)",
        "No Law IIIc\n(No Splice)",
        "Full EIMDL\n(All Laws ON)",
    ]
    dist_base = dist_map["exact_matches_only_no_approx"]["archive_bytes"]
    loo_base = loo_map["exact_matches_only_no_approx"]["archive_bytes"]
    dist_full = dist_map["full_eimdl"]["archive_bytes"]
    loo_full = loo_map["full_eimdl"]["archive_bytes"]

    dist_pct = [(dist_map[k]["archive_bytes"] / dist_base) * 100.0 for k in order_keys]
    loo_pct = [(loo_map[k]["archive_bytes"] / loo_base) * 100.0 for k in order_keys]

    x_a = np.arange(len(step_labels))
    ax1.plot(x_a, dist_pct, "o-", color=PALETTE["eimdl"], linewidth=2.0, markersize=6,
             label=f"Distant Single Donor (00|20: {dist_base/1e3:.1f} KB -> {dist_full/1e3:.1f} KB, {(dist_full/dist_base - 1)*100:.1f}%)")
    ax1.plot(x_a, loo_pct, "s--", color=PALETTE["rcmdl"], linewidth=1.8, markersize=5.5,
             label=f"Dense 23-Donor LOO (00|23: {loo_base/1e3:.1f} KB -> {loo_full/1e3:.1f} KB, {(loo_full/loo_base - 1)*100:.1f}%)")

    for i, k in enumerate(order_keys):
        db = dist_map[k]["archive_bytes"]
        lb = loo_map[k]["archive_bytes"]
        if i == 0:
            # Separate the two 100% labels cleanly at x=0
            ax1.annotate(f"{db/1000:.1f}K", (i, dist_pct[i]), textcoords="offset points",
                         xytext=(-6, -15), ha="center", fontsize=7.1, fontweight="bold", color=PALETTE["eimdl"])
            ax1.annotate(f"{lb/1000:.2f}K", (i, loo_pct[i]), textcoords="offset points",
                         xytext=(8, 6), ha="left", fontsize=7.1, fontweight="bold", color=PALETTE["rcmdl"])
        else:
            ax1.annotate(f"{db/1000:.1f}K", (i, dist_pct[i]), textcoords="offset points",
                         xytext=(0, -14 if i in [1, 3, 5] else 7), ha="center", fontsize=7.1, color=PALETTE["eimdl"])
            ax1.annotate(f"{lb/1000:.2f}K", (i, loo_pct[i]), textcoords="offset points",
                         xytext=(0, 7 if i in [1, 3, 5] else -14), ha="center", fontsize=7.1, color=PALETTE["rcmdl"])

    ax1.set_xticks(x_a)
    ax1.set_xticklabels(step_labels, rotation=15, ha="center")
    ax1.set_xlim(-0.35, 5.35)
    ax1.set_ylim(77, 111)
    ax1.set_ylabel("Relative Archive Size (% of Exact-Only Baseline)")
    ax1.set_title("(A) Real-Genome Law-by-Law Ablation Trajectory")
    ax1.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.3)

    # Panel B [Top-Right]: Donor Panel Size Scaling Curve (N = 1, 2, 4, 8, 16, 23 donors on E. coli 00.seq)
    ax2 = fig.add_subplot(gs[0, 1])
    n_donors = [1, 2, 4, 8, 16, 23]
    rcmdl_scale_kb = [273.654, 197.390, 185.481, 79.595, 65.671, 57.540]
    hgs2_endpoints_n = [1, 23]
    hgs2_endpoints_kb = [dist_map["full_eimdl"]["archive_bytes"] / 1000.0,
                         loo_map["full_eimdl"]["archive_bytes"] / 1000.0]

    ax2.plot(n_donors, rcmdl_scale_kb, "^--", color=PALETTE["rcmdl"], linewidth=1.7, markersize=6,
             label="RC-MDL v2 Panel Scaling (273.7 KB -> 57.5 KB)")
    ax2.plot(hgs2_endpoints_n, hgs2_endpoints_kb, "o-", color=PALETTE["eimdl"], linewidth=2.1, markersize=6.5,
             label=f"EIMDL / HGS2 Endpoints ({hgs2_endpoints_kb[0]:.1f} KB -> {hgs2_endpoints_kb[1]:.1f} KB)")

    for n, kb in zip(n_donors, rcmdl_scale_kb):
        if n == 8:
            ax2.annotate(f"{kb:.1f}K", (n, kb), textcoords="offset points", xytext=(-14, -14),
                         ha="center", fontsize=7.2, color="#147A36")
        elif n == 1:
            ax2.annotate(f"{kb:.1f}K", (n, kb), textcoords="offset points", xytext=(8, 6),
                         ha="center", fontsize=7.2, color="#147A36")
        else:
            ax2.annotate(f"{kb:.1f}K", (n, kb), textcoords="offset points", xytext=(0, 7),
                         ha="center", fontsize=7.2, color="#147A36")
    ax2.annotate(f"{hgs2_endpoints_kb[0]:.1f} KB\n({dist_map['full_eimdl']['bpb']:.3f} bpb)", (1, hgs2_endpoints_kb[0]),
                 textcoords="offset points", xytext=(4, -28),
                 ha="left", fontsize=7.5, fontweight="bold", color=PALETTE["eimdl"])
    ax2.annotate(f"{hgs2_endpoints_kb[1]:.1f} KB ({loo_map['full_eimdl']['bpb']:.3f} bpb)", (23, hgs2_endpoints_kb[1]),
                 textcoords="offset points", xytext=(-10, -14),
                 ha="right", fontsize=7.5, fontweight="bold", color=PALETTE["eimdl"])

    ax2.set_xscale("log", base=2)
    ax2.set_xticks(n_donors)
    ax2.set_xticklabels([str(n) for n in n_donors])
    ax2.set_xlim(0.8, 30)
    ax2.set_ylim(10, 345)
    ax2.set_xlabel("Donor Panel Size N (genomes)")
    ax2.set_ylabel("E. coli 00.seq Marginal Archive Size (KB)")
    ax2.set_title("(B) Law III: Coalescent Donor Panel Scaling")
    ax2.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.3)

    # Panel C [Bottom-Left]: Gen-1 -> Gen-3 24-Genome E. coli Progression & 100:1 Entropy Floor Proof
    ax3 = fig.add_subplot(gs[1, 0])
    prog_labels = [
        "HyperGraSS v1 Default\n(4 MiB Bank / 2 MiB Blk)",
        "xz -9 Baseline\n(Multi-FASTA)",
        "HyperGraSS v1 --pangenome\n(16 MiB Bank / 8 MiB Blk)",
        "HyperGraSS v1 + LZ77 Gap\n(2.26M Anchor + 1.55M Mod)",
        "RC-MDL v2 Prototype\n(Prim MST Collection)",
        "HGS2 / EIMDL Mode C\n(Self-Contained 1T)",
        "HGS2 / EIMDL Mode R\n(24-Genome LOO Sum)",
    ]
    prog_mb = [
        9215764 / 1e6,
        8394900 / 1e6,
        4885851 / 1e6,
        3810409 / 1e6,
        3972975 / 1e6,
        amdahl_mem["amdahl_scaling"][0]["archive_bytes"] / 1e6,
        loo_ecoli["total_archive_bytes"] / 1e6,
    ]
    prog_colors = [
        "#8D99AE",
        PALETTE["xz"],
        PALETTE["hypergrass1"],
        PALETTE["hypergrass1"],
        PALETTE["rcmdl"],
        PALETTE["eimdl_alt"],
        PALETTE["eimdl"],
    ]
    y_p = np.arange(len(prog_labels))
    bars_p = ax3.barh(y_p, prog_mb, color=prog_colors, edgecolor="black", linewidth=0.6, height=0.62)
    ax3.axvline(1.201082, color="#D62828", linestyle="--", linewidth=1.4,
                label="100:1 Compression Target (1.201 MB = 0.080 bpb)")
    ax3.set_yticks(y_p)
    ax3.set_yticklabels(prog_labels)
    ax3.invert_yaxis()
    ax3.set_xlim(0, 12.2)
    ax3.set_xlabel("Archive Size on 24-Genome E. coli (MB, 120.11 Mbp raw)")
    ax3.set_title("(C) 100:1 Entropy-Floor Proof & Gen 1 -> Gen 3 Progression")

    for bar, mb in zip(bars_p, prog_mb):
        ratio = 120.108223 / mb
        ax3.text(mb + 0.16, bar.get_y() + bar.get_height() / 2.0,
                 f"{mb:.2f} MB ({ratio:.1f}x)",
                 va="center", ha="left", fontsize=7.3,
                 fontweight="bold" if mb < 3.5 else "normal")
    ax3.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.2)

    # Panel D [Bottom-Right]: Law IV Synthetic, Xenobiological & Cryptographic Genomes
    ax4 = fig.add_subplot(gs[1, 1])
    syn_recs = mode_s_syn["synthetic_and_non_existing"]
    syn_map = {r["corpus_id"]: r["bits_per_raw_byte"] for r in syn_recs}
    syn_order = [
        ("07_hypermutator_strain05_loo", "07_strain05 (Mode R)"),
        ("02_viral_rna_acgu", "02_viral_rna_acgu"),
        ("07_hypermutator_pangenome", "07_pangenome (Mode C)"),
        ("06_denovo_grammar", "06_denovo_grammar"),
        ("01_euk_centromere", "01_euk_centromere"),
        ("04_xeno_at88", "04_xeno_at88"),
        ("05_xeno_gc86", "05_xeno_gc86"),
        ("03_chimeric_iupac_crlf", "03_chimeric_iupac"),
        ("08_crypto_white_noise", "08_crypto_noise"),
    ]
    syn_names = [lbl for _, lbl in syn_order]
    syn_bpb = [syn_map[cid] for cid, _ in syn_order]

    y_s = np.arange(len(syn_names))
    c_s = [PALETTE["eimdl"] if v < 1.0 else (PALETTE["eimdl_alt"] if v < 1.9 else "#D62828") for v in syn_bpb]
    bars_s = ax4.barh(y_s, syn_bpb, color=c_s, edgecolor="black", linewidth=0.6, height=0.65)
    ax4.axvline(2.0, color="#D62828", linestyle=":", linewidth=1.4, label="2.0000 bpb Shannon Ceiling")
    ax4.set_yticks(y_s)
    ax4.set_yticklabels(syn_names)
    ax4.invert_yaxis()
    ax4.set_xlim(0, 2.45)
    ax4.set_xlabel("Bits per Raw Input Byte")
    ax4.set_title("(D) Law IV: Non-Existing & Adversarial Genomes")

    for bar, v in zip(bars_s, syn_bpb):
        ax4.text(v + 0.03, bar.get_y() + bar.get_height() / 2.0, f"{v:.4f}",
                 va="center", ha="left", fontsize=7.3, fontweight="bold" if v > 1.99 or v < 0.35 else "normal")

    ax4.legend(loc="center right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.3)
    save_fig(fig, "fig5_laws_ablations_and_synthetic")


def make_fig6_independent_modes_pareto():
    sal_modes = load_json("benchmarks/independent/runs/salmonella/modes-subset-1/report.json")
    kleb_modes = load_json("benchmarks/independent/runs/klebsiella/modes-subset-1/report.json")
    sal_l31 = load_json("benchmarks/independent/runs/salmonella/jarvis-l31-subset-1/report.json")
    kleb_l31 = load_json("benchmarks/independent/runs/klebsiella/jarvis-l31-subset-1/report.json")

    fig = plt.figure(figsize=(13.8, 5.3))
    gs = GridSpec(1, 2, width_ratios=[1.05, 1.25], wspace=0.28)

    # Panel A: 12-Target Aggregate Pareto Frontier (Peak RSS MiB vs Compression Rate bpb)
    ax1 = fig.add_subplot(gs[0])
    total_bases = sal_modes["summaries"]["S/eimdl"]["total_bases"] + kleb_modes["summaries"]["S/eimdl"]["total_bases"]

    def agg_summary(key, src_sal, src_kleb):
        s = src_sal["summaries"][key]
        k = src_kleb["summaries"][key]
        arc = s["total_archive_bytes"] + k["total_archive_bytes"]
        bpb = (arc * 8.0) / total_bases
        rss_mib = max(s["compress_cgroup_peak"]["max_bytes"], k["compress_cgroup_peak"]["max_bytes"]) / (1024.0 * 1024.0)
        time_s = s["compress_wall_seconds_total"] + k["compress_wall_seconds_total"]
        return bpb, rss_mib, time_s, arc

    pts = [
        ("EIMDL Mode R (24-Donor LOO)", *agg_summary("R/eimdl", sal_modes, kleb_modes), PALETTE["eimdl"], "D", (12, 0), "left", "center"),
        ("EIMDL Mode S (--threads 1)", *agg_summary("S/eimdl", sal_modes, kleb_modes), PALETTE["eimdl_alt"], "o", (0, -24), "center", "top"),
        ("JARVIS3 -l 31 (Tuned High-RAM)", *agg_summary("S/jarvis3-l31", sal_l31, kleb_l31), PALETTE["jarvis_l31"], "s", (-10, -22), "center", "top"),
        ("JARVIS3 -l 7 (Default)", *agg_summary("S/jarvis3", sal_modes, kleb_modes), PALETTE["jarvis_l7"], "^", (-10, 8), "right", "bottom"),
        ("Original GraSS (-e2)", *agg_summary("S/grass", sal_modes, kleb_modes), PALETTE["grass"], "v", (10, -8), "left", "top"),
    ]

    for label, bpb, rss, t_s, arc, col, mk, xy_off, ha_align, va_align in pts:
        ax1.scatter(rss, bpb, color=col, marker=mk, s=95, edgecolor="black", linewidth=0.8, zorder=4,
                    label=f"{label}\n({arc/1e6:.2f} MB, {bpb:.4f} bpb, {rss:.0f} MiB, {t_s:.0f}s)")
        ax1.annotate(f"{bpb:.4f} bpb\n({rss:.0f} MiB)", (rss, bpb),
                     textcoords="offset points", xytext=xy_off, ha=ha_align, va=va_align,
                     fontsize=7.4, fontweight="bold" if "EIMDL" in label else "normal")

    ax1.set_xscale("log")
    ax1.set_xlim(65, 2300)
    ax1.set_ylim(-0.12, 2.38)
    ax1.set_xlabel("Peak Container Memory Footprint (cgroup RSS MiB, log scale)")
    ax1.set_ylabel("12-Genome Aggregate Rate (bits / base)")
    ax1.set_title("(A) Independent 12-Genome Rate vs. Memory Pareto Frontier", pad=10)
    # Place legend in the middle-left empty band (y = 0.25..1.25) completely clear of all 5 points
    ax1.legend(loc="center left", bbox_to_anchor=(0.02, 0.43), frameon=True,
               facecolor="white", framealpha=0.95, fontsize=6.9)

    # Panel B: Per-Target Archive Size across all 12 Independent Genomes (6 Salmonella + 6 Klebsiella)
    ax2 = fig.add_subplot(gs[1])

    def target_map(report, mode, tool):
        return {r["target"]: r["archive_bytes"] / 1000.0
                for r in report["results"] if r["mode"] == mode and r["tool"] == tool}

    sal_targets = [r["target"] for r in sal_l31["results"]]
    kleb_targets = [r["target"] for r in kleb_l31["results"]]
    all_targets = sal_targets + kleb_targets
    short_labels = [
        *[f"Sal-{t.split('_')[0]}" for t in sal_targets],
        *[f"Kleb-{t.split('_')[0]}" for t in kleb_targets],
    ]

    m_s_eimdl = {**target_map(sal_modes, "S", "eimdl"), **target_map(kleb_modes, "S", "eimdl")}
    m_r_eimdl = {**target_map(sal_modes, "R", "eimdl"), **target_map(kleb_modes, "R", "eimdl")}
    m_s_jl7 = {**target_map(sal_modes, "S", "jarvis3"), **target_map(kleb_modes, "S", "jarvis3")}
    m_s_jl31 = {**target_map(sal_l31, "S", "jarvis3-l31"), **target_map(kleb_l31, "S", "jarvis3-l31")}
    m_s_grass = {**target_map(sal_modes, "S", "grass"), **target_map(kleb_modes, "S", "grass")}

    x12 = np.arange(len(all_targets))
    ax2.plot(x12, [m_s_grass[t] for t in all_targets], "v:", color=PALETTE["grass"], markersize=4.5, linewidth=1.1,
             label="Original GraSS (-e2) Mode S: 15,797.5 KB")
    ax2.plot(x12, [m_s_jl7[t] for t in all_targets], "^--", color=PALETTE["jarvis_l7"], markersize=4.5, linewidth=1.2,
             label="JARVIS3 (-l 7) Mode S: 15,204.3 KB")
    ax2.plot(x12, [m_s_jl31[t] for t in all_targets], "s-.", color=PALETTE["jarvis_l31"], markersize=4.5, linewidth=1.2,
             label="JARVIS3 (-l 31, 1.16 GiB RAM) Mode S: 14,197.4 KB")
    ax2.plot(x12, [m_s_eimdl[t] for t in all_targets], "o-", color=PALETTE["eimdl_alt"], markersize=5, linewidth=1.6,
             label="EIMDL Mode S (113 MiB RAM): 13,868.4 KB (Wins 12/12 vs l31)")
    ax2.plot(x12, [m_r_eimdl[t] for t in all_targets], "D-", color=PALETTE["eimdl"], markersize=5, linewidth=1.8,
             label="EIMDL Mode R (24-Donor LOO): 433.5 KB (31.99x < Mode S)")

    ax2.axvline(5.5, color="black", linestyle="--", alpha=0.4, linewidth=1.0)
    ax2.set_yscale("log")
    # Mode S curves are at ~1,300 KB; setting ylim up to 95,000 KB puts the legend cleanly ABOVE all curves!
    ax2.set_ylim(0.18, 95000)
    ax2.set_xticks(x12)
    ax2.set_xticklabels(short_labels, rotation=32, ha="right")
    ax2.set_ylabel("Archive Size per Genome (KB, log scale)")
    ax2.set_title("(B) Per-Genome Archive Size Across 12 Independent Enterobacteriaceae Targets", pad=10)
    ax2.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.1)

    save_fig(fig, "fig6_independent_modes_pareto")


def make_fig7_gen1_gen2_foundations():
    """
    Dedicated 4-panel figure covering Generation 1 (HyperGraSS v1) and Generation 2 (RC-MDL v2)
    experimental records:
      (A) Generation 1 Cross-Taxon Natural Collections (6 bacterial species, BROAD.md)
      (B) Generation 1 Frozen Corpus & 100 Mbp Human hg38 Chromosome Slices (BENCHMARKS.md & SPECIALIZED.md)
      (C) Generation 1 End-to-End 13-Slice (2,000 bp) Indexed Query Amortization (END-TO-END.md)
      (D) Generation 2 (RC-MDL v2) Collection Ordering & Anchor Topologies (prototype/rcmdl/RESULTS.md)
    """
    fig = plt.figure(figsize=(14.2, 9.4))
    gs = GridSpec(2, 2, hspace=0.40, wspace=0.30)

    # Panel A: Gen 1 Cross-Taxon Natural Collections (BROAD.md)
    ax1 = fig.add_subplot(gs[0, 0])
    broad_species = [
        "C. jejuni\n(72g, 122.4 Mb)",
        "C. trachomatis\n(66g, 69.8 Mb)",
        "M. tuberculosis\n(28g, 123.7 Mb)",
        "S. pneumoniae\n(57g, 121.0 Mb)",
        "P. aeruginosa\n(19g, 123.2 Mb)",
        "S. aureus\n(43g, 123.0 Mb)",
    ]
    raw_mb_broad = [122.419352, 69.754471, 123.736606, 121.007978, 123.241439, 123.028922]
    hg1_bytes = [1535296, 333507, 1139535, 1754114, 3083542, 1807399]
    j3_bytes = [2017046, 472224, 1418155, 1899711, 3350844, 1895866]
    xz_bytes = [4211348, 862388, 2505964, 5960712, 10847460, 5672772]

    hg1_bpb = [(b * 8.0) / (m * 1e6) for b, m in zip(hg1_bytes, raw_mb_broad)]
    j3_bpb = [(b * 8.0) / (m * 1e6) for b, m in zip(j3_bytes, raw_mb_broad)]
    xz_bpb = [(b * 8.0) / (m * 1e6) for b, m in zip(xz_bytes, raw_mb_broad)]

    x_b = np.arange(len(broad_species))
    w_b = 0.26
    b_hg1 = ax1.bar(x_b - w_b, hg1_bpb, width=w_b, color=PALETTE["hypergrass1"], edgecolor="black", linewidth=0.6,
                    label="HyperGraSS v1 (--pangenome)")
    b_j3 = ax1.bar(x_b, j3_bpb, width=w_b, color=PALETTE["jarvis_l7"], edgecolor="black", linewidth=0.6,
                   label="JARVIS3 (-l 7)")
    b_xz = ax1.bar(x_b + w_b, xz_bpb, width=w_b, color=PALETTE["xz"], edgecolor="black", linewidth=0.6,
                   label="xz -9")

    for bar, b1, b2 in zip(b_hg1, hg1_bytes, j3_bytes):
        gain = (1.0 - b1 / b2) * 100.0
        h = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2.0, h * 1.10, f"-{gain:.1f}%",
                 ha="center", va="bottom", fontsize=6.7, fontweight="bold", color="#2B580C", rotation=90)

    ax1.set_yscale("log")
    ax1.set_ylim(0.02, 2.8)
    ax1.set_xticks(x_b)
    ax1.set_xticklabels(broad_species, fontsize=8.0)
    ax1.set_ylabel("Compression Rate (bits / base, log scale)")
    ax1.set_title("(A) Gen 1: Cross-Taxon Natural Bacterial Collections (6 Species)")
    ax1.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.3)

    # Panel B: Gen 1 Single-Genome QCRM & 100 Mbp Human hg38 Chromosome Slices (SPECIALIZED.md)
    ax2 = fig.add_subplot(gs[0, 1])
    single_datasets = [
        "Yeast 3-chr\n(1.36 Mbp)",
        "B. subtilis 168\n(4.22 Mbp)",
        "E. coli K-12\n(4.64 Mbp)",
        "P. aeruginosa\n(7.32 Mbp)",
        "hg38 chr1\n(100.0 Mbp)",
        "hg38 chr14\n(100.0 Mbp)",
    ]
    hg1_single_bpb = [1.855, 1.888, 1.869, 1.673, 1.650, 1.367]
    hg1_t2_bpb = [1.867, 1.895, 1.880, 1.693, 1.651, 0.0]  # chr14 measured only --single vs GraSS -e2
    grass_single_bpb = [1.963, 1.962, 1.968, 1.848, 1.781, 1.483]

    x_s = np.arange(len(single_datasets))
    w_s = 0.27
    b_s1 = ax2.bar(x_s - w_s, hg1_single_bpb, width=w_s, color=PALETTE["hypergrass1"], edgecolor="black", linewidth=0.6,
                   label="HyperGraSS v1 (--single QCRM)")
    b_s2 = ax2.bar(x_s[:5], hg1_t2_bpb[:5], width=w_s, color=PALETTE["eimdl_alt"], edgecolor="black", linewidth=0.6,
                   label="HyperGraSS v1 (--single --threads 2)")
    b_s3 = ax2.bar(x_s + w_s, grass_single_bpb, width=w_s, color=PALETTE["grass"], edgecolor="black", linewidth=0.6,
                   label="Original GraSS (-e2)")

    for bar, v1, v3 in zip(b_s1, hg1_single_bpb, grass_single_bpb):
        gain = (1.0 - v1 / v3) * 100.0
        ax2.text(bar.get_x() + bar.get_width() / 2.0, v1 + 0.02, f"-{gain:.1f}%",
                 ha="center", va="bottom", fontsize=6.7, fontweight="bold", color="#2B580C", rotation=90)

    ax2.set_ylim(1.15, 2.38)
    ax2.set_xticks(x_s)
    ax2.set_xticklabels(single_datasets, fontsize=8.0)
    ax2.set_ylabel("Compression Rate (bits / base)")
    ax2.set_title("(B) Gen 1: Single-Genome QCRM & 100 Mbp Human hg38 Slices")
    ax2.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.2)

    # Panel C: Gen 1 End-to-End 13-Slice (2,000 bp) Indexed Query Amortization (END-TO-END.md)
    ax3 = fig.add_subplot(gs[1, 0])
    q_counts = np.array([1, 10, 58, 91, 200, 500, 859, 1200])
    # Compress time + Q * per-query latency (from END-TO-END.md: 13-slice batch = 0.73s lean, 1.35s pangenome, 11.12s xz, 34.14s gzip)
    # Wait: in END-TO-END.md, break-even queries are computed per single query (or 13-slice workload).
    # Let's plot the exact 13-slice batch latency & Archive Size + Break-Even Query Count bars from END-TO-END.md!
    wf_labels = [
        "HyperGraSS v1\n(--pangenome)",
        "HyperGraSS v1\n(--lean)",
        "xz -9\n(Full Decompress)",
        "gzip -9\n(Full Decompress)",
    ]
    batch_lat_s = [1.35, 0.73, 11.12, 34.14]
    rss_mib_wf = [74.06, 43.91, 32.89, 3.45]
    arc_mb_wf = [4.909412, 8.952202, 8.394900, 34.236764]

    x_w = np.arange(len(wf_labels))
    w_w = 0.34
    bars_lat = ax3.bar(x_w - w_w / 2, batch_lat_s, width=w_w, color=PALETTE["eimdl"], edgecolor="black",
                       linewidth=0.6, label="13-Slice (2,000 bp) Query Batch Latency (s)")
    bars_arc_w = ax3.bar(x_w + w_w / 2, arc_mb_wf, width=w_w, color=PALETTE["eimdl_alt"], edgecolor="black",
                         linewidth=0.6, label="Compressed Pangenome Archive Size (MB)")

    for bar, v in zip(bars_lat, batch_lat_s):
        ax3.text(bar.get_x() + bar.get_width() / 2.0, v + 0.8, f"{v:.2f}s",
                 ha="center", va="bottom", fontsize=7.3, fontweight="bold", color=PALETTE["eimdl"])
    for bar, v in zip(bars_arc_w, arc_mb_wf):
        ax3.text(bar.get_x() + bar.get_width() / 2.0, v + 0.8, f"{v:.2f}M",
                 ha="center", va="bottom", fontsize=7.3, color="#007799")

    ax3.set_xticks(x_w)
    ax3.set_xticklabels(wf_labels)
    ax3.set_ylim(0, 43)
    ax3.set_ylabel("Seconds (Batch Query) / Megabytes (Archive)")
    ax3.set_title("(C) Gen 1: Indexed Random-Access 13-Slice Query Workload")
    ax3.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.3)

    # Panel D: Gen 2 (RC-MDL) Collection Ordering & Anchor Topologies (prototype/rcmdl/RESULTS.md)
    ax4 = fig.add_subplot(gs[1, 1])
    topo_labels = [
        "HGS2 Mode C (1T)",
        "HyperGraSS v1 (BWT+CRM)",
        "RC-MDL Prim MST",
        "RC-MDL Cluster-Prim",
        "RC-MDL Greedy Marginal",
        "RC-MDL File Order",
        "Fixed Anchors (N=20)",
        "Fixed Anchors (N=16)",
        "Fixed Anchors (N=12)",
        "Fixed Anchors (N=8)",
        "Fixed Anchors (N=4)",
        "AGC 3.2.5 (-a -s 200k)",
    ]
    topo_mb = [
        load_json("benchmarks/fusion/results/amdahl_and_memory.json")["amdahl_scaling"][0]["archive_bytes"] / 1e6,
        3810409 / 1e6,
        3972975 / 1e6,
        3982813 / 1e6,
        3990249 / 1e6,
        4056419 / 1e6,
        4115368 / 1e6,
        4921615 / 1e6,
        5494235 / 1e6,
        5818621 / 1e6,
        6621279 / 1e6,
        8540341 / 1e6,
    ]
    topo_cols = [
        PALETTE["eimdl"],
        PALETTE["hypergrass1"],
        PALETTE["rcmdl"],
        PALETTE["rcmdl"],
        "#52B788",
        "#74C69D",
        "#95D5B2",
        "#95D5B2",
        "#95D5B2",
        "#95D5B2",
        "#95D5B2",
        PALETTE["agc_single"],
    ]
    y_t = np.arange(len(topo_labels))
    bars_t = ax4.barh(y_t, topo_mb, color=topo_cols, edgecolor="black", linewidth=0.6, height=0.66)
    ax4.set_yticks(y_t)
    ax4.set_yticklabels(topo_labels, fontsize=7.8)
    ax4.invert_yaxis()
    ax4.set_xlim(0, 10.8)
    ax4.set_xlabel("24-Genome E. coli Archive Size (MB)")
    ax4.set_title("(D) Gen 2 (RC-MDL): Tree Topologies vs Fixed Anchor Sets")

    for bar, mb in zip(bars_t, topo_mb):
        bpb = (mb * 1e6 * 8.0) / 120108223.0
        ax4.text(mb + 0.12, bar.get_y() + bar.get_height() / 2.0,
                 f"{mb:.2f} MB ({bpb:.3f} bpb)",
                 va="center", ha="left", fontsize=7.0,
                 fontweight="bold" if mb < 3.98 else "normal")

    save_fig(fig, "fig7_gen1_gen2_foundations")


def make_fig8_kaggle_world_shootout_and_human():
    ws = load_json("benchmarks/kaggle/results/world_shootout_v2_results.json")
    tpu = load_json("benchmarks/kaggle/results/tpu_v3_scale_results.json")
    gpu = load_json("benchmarks/kaggle/results/gpu_v3_largescale_results.json")

    fig = plt.figure(figsize=(14.2, 10.6))
    gs = GridSpec(2, 2, hspace=0.40, wspace=0.30)

    # Panel A: 15-Codec Bacterial Mode C Shootout (Archive MB, log scale)
    ax1 = fig.add_subplot(gs[0, 0])
    bact_list = ws["bacterial_shootout"]
    sp_names = ["E. coli 24g\n(121.8 MB)", "S. enterica 25g\n(122.3 MB)", "K. pneumoniae 23g\n(124.5 MB)"]

    def get_run(runs, codec, xmx=None):
        for r in runs:
            if r.get("codec") == codec and (xmx is None or r.get("xmx") == xmx):
                return r
        return None

    codecs_a = [
        ("HGS2 (-Xmx16g)", lambda b: get_run(b["mode_c_fa"], "HGS2_collection", "-Xmx16g"), PALETTE["eimdl"]),
        ("HGS2 (-Xmx256m)", lambda b: get_run(b["mode_c_fa"], "HGS2_collection", "-Xmx256m"), PALETTE["eimdl_alt"]),
        ("MBGC -m3 -L*", lambda b: get_run(b["mode_c_fa"], "MBGC-m3-L") or get_run(b["mode_c_fa"], "MBGC-m3"), PALETTE["rcmdl"]),
        ("JARVIS3 -l 7 (.seq)", lambda b: get_run(b["mode_c_seq"], "JARVIS3-l7"), PALETTE["jarvis_l7"]),
        ("AGC 3.2.5", lambda b: get_run(b["mode_c_fa"], "AGC"), PALETTE["agc_single"]),
        ("xz -9e", lambda b: get_run(b["mode_c_fa"], "xz-9e"), PALETTE["xz"]),
        ("NAF -l 22*", lambda b: get_run(b["mode_c_fa"], "NAF-l22"), PALETTE["geco"]),
        ("Orig. GraSS -e2 (.seq)", lambda b: get_run(b["mode_c_seq"], "GraSS-e2"), PALETTE["grass"]),
    ]

    x_sp = np.arange(len(bact_list))
    n_c = len(codecs_a)
    bw = 0.098
    for idx, (lbl, fn, col) in enumerate(codecs_a):
        vals = [fn(b)["archive_bytes"] / 1e6 for b in bact_list]
        offsets = x_sp + (idx - (n_c - 1) / 2.0) * bw
        bars = ax1.bar(offsets, vals, width=bw * 0.92, color=col, edgecolor="black", linewidth=0.45, label=lbl)
        if idx in (0, 3, 4, 7):
            for bar, v in zip(bars, vals):
                ax1.text(bar.get_x() + bar.get_width() / 2.0, v * 1.07,
                         f"{v:.2f}", ha="center", va="bottom", fontsize=6.2,
                         fontweight="bold" if idx == 0 else "normal", rotation=90)

    ax1.set_yscale("log")
    ax1.set_ylim(1.5, 118.0)
    ax1.set_xticks(x_sp)
    ax1.set_xticklabels(sp_names, fontsize=8.4)
    ax1.set_ylabel("Archive Size (MB, log scale)")
    ax1.set_title("(A) Bacterial Mode C World Shootout (*seq-only FASTA lossless: NAF/MBGC)")
    ax1.legend(loc="upper left", ncol=2, frameon=True, facecolor="white", framealpha=0.95, fontsize=6.7)

    # Panel B: Human Diploid Mode C (Chr21, Chr22, Chr1, Chr1..4) & All-24 CHM13 Mode R
    ax2 = fig.add_subplot(gs[0, 1])
    dip_groups = gpu["cross_human_chm13_vs_grch38"]["diploid_mode_c"]
    dip_labels = ["Diploid Chr21\n(93.1 MB)", "Diploid Chr22\n(103.6 MB)", "Diploid Chr1\n(504.6 MB)", "Diploid Chr1..4\n(1,791.2 MB)"]
    hgs_dip = [get_run(g["runs"], "HGS2_collection", "-Xmx22g")["archive_bytes"] / 1e6 for g in dip_groups]
    agc_dip = [get_run(g["runs"], "AGC")["archive_bytes"] / 1e6 for g in dip_groups]
    zst_dip = [get_run(g["runs"], "zstd-19")["archive_bytes"] / 1e6 for g in dip_groups]

    x_d = np.arange(len(dip_labels))
    w_d = 0.25
    b_h = ax2.bar(x_d - w_d, hgs_dip, width=w_d, color=PALETTE["eimdl"], edgecolor="black", linewidth=0.5, label="EIMDL / HGS2 (-Xmx22g)")
    b_a = ax2.bar(x_d, agc_dip, width=w_d, color=PALETTE["agc_single"], edgecolor="black", linewidth=0.5, label="AGC 3.2.5 (Multi-FASTA)")
    b_z = ax2.bar(x_d + w_d, zst_dip, width=w_d, color=PALETTE["xz"], edgecolor="black", linewidth=0.5, label="zstd -19")
    for bar, v in zip(b_h, hgs_dip):
        ax2.text(bar.get_x() + bar.get_width() / 2.0, v * 1.08, f"{v:.1f}M", ha="center", va="bottom", fontsize=6.3, fontweight="bold", color=PALETTE["eimdl"])
    for bar, v in zip(b_a, agc_dip):
        ax2.text(bar.get_x() + bar.get_width() / 2.0, v * 1.08, f"{v:.1f}M", ha="center", va="bottom", fontsize=6.1)
    for bar, v in zip(b_z, zst_dip):
        ax2.text(bar.get_x() + bar.get_width() / 2.0, v * 1.08, f"{v:.1f}M", ha="center", va="bottom", fontsize=6.1)

    chm_tot = tpu["all24_chm13_given_grch38_mode_r"]["totals"]
    ax2.text(
        0.03, 0.96,
        f"All-24 T2T-CHM13v2.0 | GRCh38 Mode R (tpu_v3):\n"
        f"{chm_tot['total_chm13_raw_bytes']/1e9:.3f} GB -> {chm_tot['total_archive_bytes']/1e6:.2f} MB "
        f"({chm_tot['overall_ratio']:.2f}x, {chm_tot['overall_bpb']:.4f} bpb, 24/24 SHA-256)",
        transform=ax2.transAxes, va="top", ha="left", fontsize=7.2,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#F0F8FF", edgecolor=PALETTE["eimdl"], linewidth=0.9),
    )
    ax2.set_yscale("log")
    ax2.set_ylim(5.0, 1450.0)
    ax2.set_xticks(x_d)
    ax2.set_xticklabels(dip_labels, fontsize=8.0)
    ax2.set_ylabel("Archive Size (MB, log scale)")
    ax2.set_title("(B) Human Diploid Mode C & All-24 CHM13v2.0 Mode R")
    ax2.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.2)

    # Panel C: Bounded-Heap Scaling Matrix (100M to 1.41 Gbp)
    ax3 = fig.add_subplot(gs[1, 0])
    sweep = [r for r in tpu["heap_floor_scaling_sweep"] if r.get("status") == "ok"]
    covid_mega_512 = get_run(gpu["covid19_viral_scaling"][-1]["runs"], "HGS2_collection", "-Xmx512m")
    labels_c = [
        f"{int(round(r['raw_bytes']/1e6))}M@{r['xmx'].replace('-Xmx','')}" for r in sweep
    ] + ["1409M@512m\n(COVID Mega)"]
    bpb_c = [r["bpb"] * 1000.0 for r in sweep] + [covid_mega_512["bpb"] * 1000.0]  # millibits/base
    rss_c = [r["enc_peak_rss_mb"] for r in sweep] + [covid_mega_512["enc_peak_rss_mb"]]

    x_c = np.arange(len(labels_c))
    ax3_r = ax3.twinx()
    bars_r = ax3.bar(x_c - 0.18, bpb_c, width=0.36, color=PALETTE["eimdl"], edgecolor="black", linewidth=0.5, label="Compression Rate (mbpb)")
    line_rss = ax3_r.plot(x_c + 0.18, rss_c, color=PALETTE["jarvis_l31"], marker="o", linewidth=1.7, markersize=4.5, label="Peak RSS (MB)")
    ax3.set_xticks(x_c)
    ax3.set_xticklabels(labels_c, rotation=42, ha="right", fontsize=6.8)
    ax3.set_ylim(0, 15.5)
    ax3_r.set_yscale("log")
    ax3_r.set_ylim(80, 7000)
    ax3.set_ylabel("Rate (millibits / base)", color=PALETTE["eimdl"])
    ax3_r.set_ylabel("Encoder Peak RSS (MB, log)", color=PALETTE["jarvis_l31"])
    ax3.set_title("(C) Law VII Bounded-Heap Scaling Sweep (13/13 Matrix + 1.41 Gbp@-Xmx512m)")
    h1, l1 = ax3.get_legend_handles_labels()
    h2, l2 = ax3_r.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, loc="upper left", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.0)

    # Panel D: Tesla T4 CUDA Acceleration & Adversarial Security Checks
    ax4 = fig.add_subplot(gs[1, 1])
    cuda_rows = gpu["cuda_gpu_experiments"]["sketch_matrix_bench"]
    dims = [f"{r['n_units']}x{r['k_dim']}" for r in cuda_rows]
    cpu_ms = [r["cpu_numpy_ms"] for r in cuda_rows]
    gpu_ms = [r["gpu_cuda_ms"] for r in cuda_rows]
    sp_x = [r["gpu_speedup_x"] for r in cuda_rows]

    x_g = np.arange(len(dims))
    w_g = 0.34
    b_cpu = ax4.bar(x_g - w_g / 2, cpu_ms, width=w_g, color=PALETTE["jarvis_l7"], edgecolor="black", linewidth=0.5, label="CPU NumPy (ms)")
    b_gpu = ax4.bar(x_g + w_g / 2, gpu_ms, width=w_g, color=PALETTE["eimdl"], edgecolor="black", linewidth=0.5, label="Tesla T4 CUDA (ms)")
    for bar, sp in zip(b_gpu, sp_x):
        ax4.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height() * 1.18, f"{sp:.1f}x", ha="center", va="bottom", fontsize=7.4, fontweight="bold", color=PALETTE["eimdl"])
    ax4.set_yscale("log")
    ax4.set_ylim(0.8, 160000)
    ax4.set_xticks(x_g)
    ax4.set_xticklabels(dims, fontsize=8.4)
    ax4.set_ylabel("All-Pairs Sketch Matrix Latency (ms, log)")
    ax4.set_title("(D) 2x Tesla T4 CUDA Speedup & Adversarial Torture Integrity")
    ax4.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.1)

    wn_hgs = get_run(ws["adversarial_and_bounded"]["white_noise_10mb"], "HGS2_single")
    path_hgs = get_run(ws["adversarial_and_bounded"]["pathology_fasta"], "HGS2_collection", "-Xmx4g")
    ax4.text(
        0.98, 0.96,
        f"Adversarial Torture (world_shootout_v2):\n"
        f"- 10M White Noise: HGS2 {wn_hgs['archive_bytes']:,} B ({wn_hgs['bpb']:.4f} bpb)\n"
        f"- 15.6M Pathology FASTA: HGS2 {path_hgs['archive_bytes']:,} B (SHA-256 PASS)\n"
        f"  GraSS-e2 & NAF-l22: sequence SHA-256 MISMATCH",
        transform=ax4.transAxes, va="top", ha="right", fontsize=6.6,
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#FFF9F0", edgecolor=PALETTE["jarvis_l7"], linewidth=0.9),
    )

    save_fig(fig, "fig8_kaggle_world_shootout_and_human")


if __name__ == "__main__":
    make_fig1_regime_overview()
    make_fig2_mode_r_loo_73genomes()
    make_fig3_mode_c_pangenomes()
    make_fig4_amdahl_and_bounded_memory()
    make_fig5_laws_ablations_and_synthetic()
    make_fig6_independent_modes_pareto()
    make_fig7_gen1_gen2_foundations()
    make_fig8_kaggle_world_shootout_and_human()
    print("[ALL 8 FIGURES GENERATED SUCCESSFULLY]")
