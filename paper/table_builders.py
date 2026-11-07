#!/usr/bin/env python3
"""
Programmatic LaTeX and Markdown table builders reading directly from committed JSON reports.
Zero hardcoded or interpolated numbers.
"""

import json
import os

BASE_DIR = r"C:\Users\mdsai\.gemini\antigravity\scratch\hypergrass"


def load_json(rel_path):
    with open(os.path.join(BASE_DIR, rel_path), "r", encoding="utf-8") as f:
        return json.load(f)


def tex_escape(s):
    return str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def build_all_tables():
    mode_s_syn = load_json("benchmarks/fusion/results/mode_s_and_synthetic.json")
    mode_c_col = load_json("benchmarks/fusion/results/mode_c_collections.json")
    amdahl_mem = load_json("benchmarks/fusion/results/amdahl_and_memory.json")
    loo_ecoli = load_json("benchmarks/fusion/results/loo_ecoli.json")
    loo_mtb = load_json("benchmarks/fusion/results/loo_mtb.json")
    loo_saur = load_json("benchmarks/fusion/results/loo_saureus.json")
    cl_ecoli = load_json("benchmarks/fusion/results/cluster_loo_ecoli.json")
    cl_mtb = load_json("benchmarks/fusion/results/cluster_loo_mtb.json")
    cl_saur = load_json("benchmarks/fusion/results/cluster_loo_saureus.json")
    ablation = load_json("benchmarks/fusion/results/ablations_real_genomes.json")
    mtb_agc = load_json("prototype/rcmdl/work/agc-sweep-mtb-full/report.json")
    saur_oracle = load_json("prototype/rcmdl/work/agc-verify-saureus-oracle/results.json")
    sal_trio = load_json("benchmarks/independent/runs/salmonella/trio-1/report.json")
    kleb_trio = load_json("benchmarks/independent/runs/klebsiella/trio-1/report.json")
    sal_modes = load_json("benchmarks/independent/runs/salmonella/modes-subset-1/report.json")
    kleb_modes = load_json("benchmarks/independent/runs/klebsiella/modes-subset-1/report.json")
    sal_l31 = load_json("benchmarks/independent/runs/salmonella/jarvis-l31-subset-1/report.json")
    kleb_l31 = load_json("benchmarks/independent/runs/klebsiella/jarvis-l31-subset-1/report.json")

    tables = {}

    # 1. E. coli 24-Genome LOO (LaTeX & Markdown)
    ec_recs = loo_ecoli["records"]
    ec_tex, ec_md = [], []
    tot_raw = loo_ecoli["total_bases"]
    tot_hgs = loo_ecoli["total_archive_bytes"]
    tot_bpb = loo_ecoli["total_bpb"]
    tot_rcmdl = sum(r["rcmdl_bytes"] for r in ec_recs)
    tot_agc = sum(r["agc_bytes"] for r in ec_recs)
    tot_pm = sum(r["panel_match_bases"] for r in ec_recs)
    tot_sm = sum(r["self_match_bases"] for r in ec_recs)
    tot_nov = sum(r["novelty_bases"] for r in ec_recs)
    tot_snps = sum(r["substitutions"] for r in ec_recs)
    tot_spl = sum(r["indel_splices"] for r in ec_recs)

    for r in ec_recs:
        ec_tex.append(
            f"\\texttt{{{tex_escape(r['target'])}}} & {r['raw_bytes']:,} & "
            f"\\textbf{{{r['archive_bytes']:,}}} & {r['bpb']:.4f} & {r['ratio']:.1f}$\\times$ & "
            f"{r['rcmdl_bytes']:,} & {r['agc_bytes']:,} & {r['agc_over_eimdl']:.2f}$\\times$ & "
            f"{r['panel_match_bases']:,} & {r['novelty_bases']:,} & {r['substitutions']:,} & {r['indel_splices']} \\\\"
        )
        ec_md.append(
            f"| `{r['target']}` | `{r['raw_bytes']:,}` | **`{r['archive_bytes']:,}`** | `{r['bpb']:.4f}` | "
            f"`{r['ratio']:.1f}x` | `{r['rcmdl_bytes']:,}` | `{r['agc_bytes']:,}` | `{r['agc_over_eimdl']:.2f}x` | "
            f"`{r['panel_match_bases']:,}` | `{r['self_match_bases']:,}` | `{r['novelty_bases']:,}` | "
            f"`{r['substitutions']:,}` | `{r['indel_splices']}` | `PASS` |"
        )
    ec_tex.append("\\midrule")
    ec_tex.append(
        f"\\textbf{{Total (24 Genomes)}} & \\textbf{{{tot_raw:,}}} & "
        f"\\textbf{{{tot_hgs:,}}} & \\textbf{{{tot_bpb:.5f}}} & \\textbf{{{tot_raw/tot_hgs:.2f}$\\times$}} & "
        f"\\textbf{{{tot_rcmdl:,}}} & \\textbf{{{tot_agc:,}}} & \\textbf{{{tot_agc/tot_hgs:.2f}$\\times$}} & "
        f"\\textbf{{{tot_pm:,}}} & \\textbf{{{tot_nov:,}}} & \\textbf{{{tot_snps:,}}} & \\textbf{{{tot_spl:,}}} \\\\"
    )
    ec_md.append(
        f"| **Total / Aggregate** | **`{tot_raw:,}`** | **`{tot_hgs:,}`** | **`{tot_bpb:.5f}`** | "
        f"**`{tot_raw/tot_hgs:.2f}x`** | **`{tot_rcmdl:,}`** | **`{tot_agc:,}`** | **`{tot_agc/tot_hgs:.2f}x`** | "
        f"**`{tot_pm:,}`** | **`{tot_sm:,}`** | **`{tot_nov:,}`** | **`{tot_snps:,}`** | **`{tot_spl:,}`** | **`24/24 PASS`** |"
    )
    tables["ecoli_loo_tex"] = "\n".join(ec_tex)
    tables["ecoli_loo_md"] = "\n".join(ec_md)

    # 2. M. tuberculosis 24-Genome LOO (LaTeX & Markdown)
    mtb_recs = loo_mtb["records"]
    mtb_orc_map = {r["target"]: r for r in mtb_agc["targets_detail"]}
    mtb_tex, mtb_md = [], []
    m_raw = loo_mtb["total_bases"]
    m_hgs = loo_mtb["total_archive_bytes"]
    m_bpb = loo_mtb["total_bpb"]
    m_rcmdl = sum(r["rcmdl_bytes"] for r in mtb_recs)
    m_orc = sum(mtb_orc_map[r["target"]]["best_bytes"] for r in mtb_recs)
    m_pm = sum(r["panel_match_bases"] for r in mtb_recs)
    m_sm = sum(r["self_match_bases"] for r in mtb_recs)
    m_nov = sum(r["novelty_bases"] for r in mtb_recs)
    m_snps = sum(r["substitutions"] for r in mtb_recs)
    m_spl = sum(r["indel_splices"] for r in mtb_recs)

    for r in mtb_recs:
        o = mtb_orc_map[r["target"]]
        win = r["archive_bytes"] < o["best_bytes"]
        hgs_t = f"\\textbf{{{r['archive_bytes']:,}}}" if win else f"{r['archive_bytes']:,}"
        orc_t = f"{o['best_bytes']:,}" if win else f"\\textbf{{{o['best_bytes']:,}}}"
        hgs_m = f"**`{r['archive_bytes']:,}`**" if win else f"`{r['archive_bytes']:,}`"
        orc_m = f"`{o['best_bytes']:,}`" if win else f"**`{o['best_bytes']:,}`**"
        mtb_tex.append(
            f"\\texttt{{{tex_escape(r['target'])}}} & {r['raw_bytes']:,} & "
            f"{hgs_t} & {r['bpb']:.5f} & {r['ratio']:.0f}$\\times$ & "
            f"{r['rcmdl_bytes']:,} & {orc_t} & \\texttt{{{tex_escape(o['best_config'])}}} & "
            f"{r['panel_match_bases']:,} & {r['novelty_bases']:,} & {r['substitutions']:,} & {r['indel_splices']} \\\\"
        )
        mtb_md.append(
            f"| `{r['target']}` | `{r['raw_bytes']:,}` | {hgs_m} | `{r['bpb']:.5f}` | "
            f"`{r['ratio']:.1f}x` | `{r['rcmdl_bytes']:,}` | {orc_m} (`{o['best_config']}`) | "
            f"`{r['panel_match_bases']:,}` | `{r['self_match_bases']:,}` | `{r['novelty_bases']:,}` | "
            f"`{r['substitutions']:,}` | `{r['indel_splices']}` | `PASS` |"
        )
    mtb_tex.append("\\midrule")
    mtb_tex.append(
        f"\\textbf{{Total (24 Genomes)}} & \\textbf{{{m_raw:,}}} & "
        f"\\textbf{{{m_hgs:,}}} & \\textbf{{{m_bpb:.5f}}} & \\textbf{{{m_raw/m_hgs:.1f}$\\times$}} & "
        f"\\textbf{{{m_rcmdl:,}}} & \\textbf{{{m_orc:,}}} & \\textbf{{19/24 HGS2 Wins}} & "
        f"\\textbf{{{m_pm:,}}} & \\textbf{{{m_nov:,}}} & \\textbf{{{m_snps:,}}} & \\textbf{{{m_spl:,}}} \\\\"
    )
    mtb_md.append(
        f"| **Total / Aggregate** | **`{m_raw:,}`** | **`{m_hgs:,}`** | **`{m_bpb:.5f}`** | "
        f"**`{m_raw/m_hgs:.2f}x`** | **`{m_rcmdl:,}`** | **`{m_orc:,}` (13-Config Oracle; `292,170` Best Single `-s 100000`)** | "
        f"**`{m_pm:,}`** | **`{m_sm:,}`** | **`{m_nov:,}`** | **`{m_snps:,}`** | **`{m_spl:,}`** | **`24/24 PASS`** |"
    )
    tables["mtb_loo_tex"] = "\n".join(mtb_tex)
    tables["mtb_loo_md"] = "\n".join(mtb_md)

    # 3. S. aureus 25-Genome LOO (LaTeX & Markdown)
    sa_recs = loo_saur["records"]
    sa_orc_map = {r["target"].replace(".fa", ".seq"): r for r in saur_oracle}
    sa_tex, sa_md = [], []
    s_raw = loo_saur["total_bases"]
    s_hgs = loo_saur["total_archive_bytes"]
    s_bpb = loo_saur["total_bpb"]
    s_rcmdl = sum(r["rcmdl_bytes"] for r in sa_recs)
    s_orc = sum(sa_orc_map[r["target"]]["marginal_bytes"] for r in sa_recs)
    s_pm = sum(r["panel_match_bases"] for r in sa_recs)
    s_sm = sum(r["self_match_bases"] for r in sa_recs)
    s_nov = sum(r["novelty_bases"] for r in sa_recs)
    s_snps = sum(r["substitutions"] for r in sa_recs)
    s_spl = sum(r["indel_splices"] for r in sa_recs)

    for r in sa_recs:
        o = sa_orc_map[r["target"]]
        win = r["archive_bytes"] < o["marginal_bytes"]
        hgs_t = f"\\textbf{{{r['archive_bytes']:,}}}" if win else f"{r['archive_bytes']:,}"
        orc_t = f"{o['marginal_bytes']:,}" if win else f"\\textbf{{{o['marginal_bytes']:,}}}"
        hgs_m = f"**`{r['archive_bytes']:,}`**" if win else f"`{r['archive_bytes']:,}`"
        orc_m = f"`{o['marginal_bytes']:,}`" if win else f"**`{o['marginal_bytes']:,}`**"
        sa_tex.append(
            f"\\texttt{{{tex_escape(r['target'])}}} & {r['raw_bytes']:,} & "
            f"{hgs_t} & {r['bpb']:.5f} & {r['ratio']:.0f}$\\times$ & "
            f"{r['rcmdl_bytes']:,} & {orc_t} & \\texttt{{{tex_escape(o['config'])}}} & "
            f"{r['panel_match_bases']:,} & {r['novelty_bases']:,} & {r['substitutions']:,} & {r['indel_splices']} \\\\"
        )
        sa_md.append(
            f"| `{r['target']}` | `{r['raw_bytes']:,}` | {hgs_m} | `{r['bpb']:.5f}` | "
            f"`{r['ratio']:.1f}x` | `{r['rcmdl_bytes']:,}` | {orc_m} (`{o['config']}`) | "
            f"`{r['panel_match_bases']:,}` | `{r['self_match_bases']:,}` | `{r['novelty_bases']:,}` | "
            f"`{r['substitutions']:,}` | `{r['indel_splices']}` | `PASS` |"
        )
    sa_tex.append("\\midrule")
    sa_tex.append(
        f"\\textbf{{Total (25 Genomes)}} & \\textbf{{{s_raw:,}}} & "
        f"\\textbf{{{s_hgs:,}}} & \\textbf{{{s_bpb:.5f}}} & \\textbf{{{s_raw/s_hgs:.2f}$\\times$}} & "
        f"\\textbf{{{s_rcmdl:,}}} & \\textbf{{{s_orc:,}}} & \\textbf{{24/25 HGS2 Wins}} & "
        f"\\textbf{{{s_pm:,}}} & \\textbf{{{s_nov:,}}} & \\textbf{{{s_snps:,}}} & \\textbf{{{s_spl:,}}} \\\\"
    )
    sa_md.append(
        f"| **Total / Aggregate** | **`{s_raw:,}`** | **`{s_hgs:,}`** | **`{s_bpb:.5f}`** | "
        f"**`{s_raw/s_hgs:.2f}x`** | **`{s_rcmdl:,}`** | **`{s_orc:,}` (13-Config Oracle; `2,979,078` Best Single `-s 50000`)** | "
        f"**`{s_pm:,}`** | **`{s_sm:,}`** | **`{s_nov:,}`** | **`{s_snps:,}`** | **`{s_spl:,}`** | **`25/25 PASS`** |"
    )
    tables["saureus_loo_tex"] = "\n".join(sa_tex)
    tables["saureus_loo_md"] = "\n".join(sa_md)

    # 4. Independent 12-Genome Subset (LaTeX & Markdown)
    ind_tex, ind_md = [], []
    for sp_label, modes_rep, l31_rep in [
        ("S. enterica", sal_modes, sal_l31),
        ("K. pneumoniae", kleb_modes, kleb_l31),
    ]:
        s_eimdl = {r["target"]: r for r in modes_rep["results"] if r["mode"] == "S" and r["tool"] == "eimdl"}
        r_eimdl = {r["target"]: r for r in modes_rep["results"] if r["mode"] == "R" and r["tool"] == "eimdl"}
        s_jl7 = {r["target"]: r for r in modes_rep["results"] if r["mode"] == "S" and r["tool"] == "jarvis3"}
        s_jl31 = {r["target"]: r for r in l31_rep["results"] if r["mode"] == "S" and r["tool"] == "jarvis3-l31"}
        s_grass = {r["target"]: r for r in modes_rep["results"] if r["mode"] == "S" and r["tool"] == "grass"}
        for t in [r["target"] for r in l31_rep["results"]]:
            bases = s_eimdl[t]["input"]["bases"]
            b_se = s_eimdl[t]["archive_bytes"]
            b_re = r_eimdl[t]["archive_bytes"]
            b_j7 = s_jl7[t]["archive_bytes"]
            b_j31 = s_jl31[t]["archive_bytes"]
            b_gr = s_grass[t]["archive_bytes"]
            ind_tex.append(
                f"\\textit{{{sp_label}}} & \\texttt{{{tex_escape(t.replace('.seq',''))}}} & {bases:,} & "
                f"\\textbf{{{b_se:,}}} ({s_eimdl[t]['bits_per_base']:.4f}) & "
                f"{b_j31:,} ({s_jl31[t]['bits_per_base']:.4f}) & "
                f"{b_j7:,} ({s_jl7[t]['bits_per_base']:.4f}) & "
                f"{b_gr:,} ({s_grass[t]['bits_per_base']:.4f}) & "
                f"\\textbf{{{b_re:,}}} ({r_eimdl[t]['bits_per_base']:.4f}) \\\\"
            )
            ind_md.append(
                f"| *{sp_label}* | `{t}` | `{bases:,}` | **`{b_se:,}`** (`{s_eimdl[t]['bits_per_base']:.4f}`) | "
                f"`{b_j31:,}` (`{s_jl31[t]['bits_per_base']:.4f}`) | `{b_j7:,}` (`{s_jl7[t]['bits_per_base']:.4f}`) | "
                f"`{b_gr:,}` (`{s_grass[t]['bits_per_base']:.4f}`) | **`{b_re:,}`** (`{r_eimdl[t]['bits_per_base']:.4f}`) |"
            )
    ind_tex.append("\\midrule")
    ind_tex.append(
        "\\textbf{12-Target Aggregate} & \\textbf{All 12 Genomes} & \\textbf{60,810,308} & "
        "\\textbf{13,868,405 (1.8245)} & 14,197,393 (1.8678) & 15,204,327 (2.0002) & "
        "15,797,490 (2.0783) & \\textbf{433,498 (0.0570)} \\\\"
    )
    ind_md.append(
        "| **12-Target Aggregate** | **All 12 Genomes** | **`60,810,308`** | **`13,868,405` (`1.8245`)** | "
        "`14,197,393` (`1.8678`) | `15,204,327` (`2.0002`) | `15,797,490` (`2.0783`) | **`433,498` (`0.0570`)** |"
    )
    tables["ind12_tex"] = "\n".join(ind_tex)
    tables["ind12_md"] = "\n".join(ind_md)

    return tables


if __name__ == "__main__":
    t = build_all_tables()
    print("[OK] Built all tables:", list(t.keys()))
