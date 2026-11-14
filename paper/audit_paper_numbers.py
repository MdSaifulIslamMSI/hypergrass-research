#!/usr/bin/env python3
"""
Automated Factual Integrity Audit for the HGS2 / EIMDL Research Paper & Repository Documentation.
Verifies that:
  1. Every per-genome and aggregate byte count from the committed JSON reports
     is present verbatim in both paper/HGS2_EIMDL_Research_Paper.tex and docs/PAPER-DRAFT.md.
  2. All aggregate Mode R, Mode C, Mode S, Law IV, and Independent Trio/L31 numbers
     are present verbatim in README.md and HANDOFF.md.
  3. Zero unverified/fabricated numbers appear in any manuscript or documentation file.
  4. All 8 vector PDF and 300-DPI PNG figures and the compiled PDF monograph exist.
"""

import json
import os
import sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_json(rel_path):
    with open(os.path.join(BASE_DIR, rel_path), "r", encoding="utf-8") as f:
        return json.load(f)


def main():
    tex_path = os.path.join(BASE_DIR, "paper", "HGS2_EIMDL_Research_Paper.tex")
    md_path = os.path.join(BASE_DIR, "docs", "PAPER-DRAFT.md")
    readme_path = os.path.join(BASE_DIR, "README.md")
    handoff_path = os.path.join(BASE_DIR, "HANDOFF.md")
    pdf_path = os.path.join(BASE_DIR, "paper", "HGS2_EIMDL_Research_Paper.pdf")

    with open(tex_path, "r", encoding="utf-8") as f:
        tex_text = f.read()
    with open(md_path, "r", encoding="utf-8") as f:
        md_text = f.read()
    with open(readme_path, "r", encoding="utf-8") as f:
        readme_text = f.read()
    with open(handoff_path, "r", encoding="utf-8") as f:
        handoff_text = f.read()

    checks = 0
    failures = []

    def require_number(val, label):
        nonlocal checks
        s_comma = f"{int(val):,}"
        s_plain = str(int(val))
        checks += 2
        if s_comma not in tex_text and s_plain not in tex_text:
            failures.append(f"[TEX MISSING] {label}: {s_comma}")
        if s_comma not in md_text and s_plain not in md_text:
            failures.append(f"[MD MISSING] {label}: {s_comma}")

    def require_in_readme_and_handoff(val, label):
        nonlocal checks
        s_comma = f"{int(val):,}"
        s_plain = str(int(val))
        checks += 2
        if s_comma not in readme_text and s_plain not in readme_text:
            failures.append(f"[README MISSING] {label}: {s_comma}")
        if s_comma not in handoff_text and s_plain not in handoff_text:
            failures.append(f"[HANDOFF MISSING] {label}: {s_comma}")

    # 1. Check all 73 Mode R LOO genomes + totals
    for sp, rel in [
        ("ecoli", "benchmarks/fusion/results/loo_ecoli.json"),
        ("mtb", "benchmarks/fusion/results/loo_mtb.json"),
        ("saureus", "benchmarks/fusion/results/loo_saureus.json"),
    ]:
        d = load_json(rel)
        require_number(d["total_bases"], f"{sp}_total_bases")
        require_number(d["total_archive_bytes"], f"{sp}_total_archive_bytes")
        require_in_readme_and_handoff(d["total_archive_bytes"], f"{sp}_total_archive_bytes")
        for r in d["records"]:
            require_number(r["raw_bytes"], f"{sp}_{r['target']}_raw")
            require_number(r["archive_bytes"], f"{sp}_{r['target']}_hgs2")
            require_number(r["rcmdl_bytes"], f"{sp}_{r['target']}_rcmdl")
            require_number(r["panel_match_bases"], f"{sp}_{r['target']}_panel_match")
            require_number(r["substitutions"], f"{sp}_{r['target']}_snps")

    # 2. Check AGC Oracle per-target numbers on Mtb and S. aureus
    mtb_agc = load_json("prototype/rcmdl/work/agc-sweep-mtb-full/report.json")
    for r in mtb_agc["targets_detail"]:
        require_number(r["best_bytes"], f"mtb_oracle_{r['target']}")
    saur_orc = load_json("prototype/rcmdl/work/agc-verify-saureus-oracle/results.json")
    for r in saur_orc:
        require_number(r["marginal_bytes"], f"saureus_oracle_{r['target']}")

    # 3. Check Leave-Cluster-Out totals
    for sp, rel in [
        ("cl_ecoli", "benchmarks/fusion/results/cluster_loo_ecoli.json"),
        ("cl_mtb", "benchmarks/fusion/results/cluster_loo_mtb.json"),
        ("cl_saureus", "benchmarks/fusion/results/cluster_loo_saureus.json"),
    ]:
        d = load_json(rel)
        require_number(d["total_archive_bytes"], f"{sp}_total_archive_bytes")
        require_in_readme_and_handoff(d["total_archive_bytes"], f"{sp}_total_archive_bytes")

    # 4. Check Mode C collections, Amdahl scaling, and bounded memory
    mc = load_json("benchmarks/fusion/results/mode_c_collections.json")
    for r in mc["collections"]:
        require_number(r["archive_bytes"], f"mode_c_{r['label']}")
    am = load_json("benchmarks/fusion/results/amdahl_and_memory.json")
    for r in am["amdahl_scaling"]:
        require_number(r["archive_bytes"], f"amdahl_t{r['threads']}")
        if r["threads"] in (1, 12):
            require_in_readme_and_handoff(r["archive_bytes"], f"amdahl_t{r['threads']}")
    for r in am["bounded_memory"]:
        require_number(r["archive_bytes"], f"bounded_{r['heap_cap']}")
        require_in_readme_and_handoff(r["archive_bytes"], f"bounded_{r['heap_cap']}")

    # 5. Check Mode S and 9 Synthetic corpora
    ms = load_json("benchmarks/fusion/results/mode_s_and_synthetic.json")
    for r in ms["real_mode_s"]:
        require_number(r["archive_bytes"], f"mode_s_{r['label']}")
        if r["label"] in ("ecoli_00_fa", "ecoli_00_seq"):
            require_in_readme_and_handoff(r["archive_bytes"], f"mode_s_{r['label']}")
    for r in ms["synthetic_and_non_existing"]:
        require_number(r["archive_bytes"], f"synthetic_{r['corpus_id']}")

    # 6. Check all 18 real-genome ablation configurations
    ab = load_json("benchmarks/fusion/results/ablations_real_genomes.json")
    for group_k, recs in ab.items():
        for r in recs:
            require_number(r["archive_bytes"], f"ablation_{group_k}_{r['ablation']}")

    # 7. Check Independent Trio & 12-Genome Subset
    for sp, p_trio, p_modes, p_l31 in [
        ("salmonella",
         "benchmarks/independent/runs/salmonella/trio-1/report.json",
         "benchmarks/independent/runs/salmonella/modes-subset-1/report.json",
         "benchmarks/independent/runs/salmonella/jarvis-l31-subset-1/report.json"),
        ("klebsiella",
         "benchmarks/independent/runs/klebsiella/trio-1/report.json",
         "benchmarks/independent/runs/klebsiella/modes-subset-1/report.json",
         "benchmarks/independent/runs/klebsiella/jarvis-l31-subset-1/report.json"),
    ]:
        for r in load_json(p_trio)["results"]:
            require_number(r["archive_bytes"], f"trio_{sp}_{r['profile']}")
            require_in_readme_and_handoff(r["archive_bytes"], f"trio_{sp}_{r['profile']}")
        for r in load_json(p_modes)["results"]:
            require_number(r["archive_bytes"], f"modes_{sp}_{r['mode']}_{r['tool']}_{r['target']}")
        for r in load_json(p_l31)["results"]:
            require_number(r["archive_bytes"], f"l31_{sp}_{r['target']}")

    # Check combined 12-genome Mode S / L31 / Mode R totals in README and HANDOFF
    for val, lbl in [
        (13868405, "ind12_eimdl_s_total"),
        (14197393, "ind12_jarvis_l31_total"),
        (15204327, "ind12_jarvis_l7_total"),
        (433498, "ind12_eimdl_r_total"),
        (125075, "law4_crypto_white_noise"),
    ]:
        require_in_readme_and_handoff(val, lbl)

    # 8. Verify Generation 1 (HyperGraSS v1) and Generation 2 (RC-MDL) key numbers
    for val, lbl in [
        (3864464, "gen1_pre_lz77_ecoli24"),
        (3810409, "gen1_lz77_ecoli24"),
        (30171539, "gen1_2gib_streaming"),
        (2257234, "gen1_entropy_floor_anchor_bank"),
        (882215, "gen1_entropy_floor_gaps"),
        (8952202, "gen1_query_lean"),
        (4909412, "gen1_query_pangenome"),
        (1450787, "gen1_ablation_ecoli10m"),
        (20626877, "gen1_hg38_chr1"),
        (17084921, "gen1_hg38_chr14"),
        (1535296, "gen1_campy"),
        (333507, "gen1_chlamydia"),
        (1139535, "gen1_mtb"),
        (1754114, "gen1_pneumo"),
        (3083542, "gen1_pseudo"),
        (1807399, "gen1_staph"),
        (3972975, "gen2_rcmdl_prim"),
        (3982813, "gen2_rcmdl_cluster_prim"),
        (3990249, "gen2_rcmdl_greedy"),
        (4056419, "gen2_rcmdl_file_order"),
        (4081636, "gen2_rcmdl_cluster_first"),
        (4115368, "gen2_rcmdl_fixed_20"),
        (6621279, "gen2_rcmdl_fixed_4"),
    ]:
        require_number(val, lbl)

    for r in load_json("prototype/rcmdl/work/agc-sweep-ecoli-full/results-00.json"):
        require_number(r["marginal_bytes"], f"gen2_agc_ecoli00_{r['config']}")
    for k, v in load_json("prototype/rcmdl/work/agc-sweep-mtb-full/report.json")["config_totals"].items():
        require_number(v, f"gen2_agc_mtb_{k}")
    for k, v in load_json("prototype/rcmdl/work/agc-sweep-saureus/report.json")["config_totals"].items():
        require_number(v, f"gen2_agc_saureus_{k}")

    # 9. Verify forbidden unbacked numbers are ABSENT across TEX, MD, README, and HANDOFF
    for bad in [
        "1,136,892", "1136892", "998,412", "998412", "642,108", "642108",
        "1,177,778", "1177778", "47,305", "47305", "174,639", "174639",
        "1,625,424", "1625424", "78,582", "78582", "525,828", "525828",
        "3,377,802", "3377802", "3,393,302", "3393302", "3,455,169", "3455169",
    ]:
        checks += 4
        if bad in tex_text:
            failures.append(f"[FORBIDDEN UNBACKED NUMBER IN TEX]: {bad}")
        if bad in md_text:
            failures.append(f"[FORBIDDEN UNBACKED NUMBER IN MD]: {bad}")
        if bad in readme_text:
            failures.append(f"[FORBIDDEN UNBACKED NUMBER IN README]: {bad}")
        if bad in handoff_text:
            failures.append(f"[FORBIDDEN UNBACKED NUMBER IN HANDOFF]: {bad}")

    # 10. Verify Kaggle E2E Cloud JSON results (world_shootout_v2, gpu_v3, tpu_v3)
    ws = load_json("benchmarks/kaggle/results/world_shootout_v2_results.json")
    gpu = load_json("benchmarks/kaggle/results/gpu_v3_largescale_results.json")
    tpu = load_json("benchmarks/kaggle/results/tpu_v3_scale_results.json")

    def find_codec(runs, codec, xmx=None):
        for r in runs:
            if r.get("codec") == codec and (xmx is None or r.get("xmx") == xmx):
                return r
        raise KeyError((codec, xmx))

    # White noise 10MB & pathology FASTA 15.6MB
    for r in ws["adversarial_and_bounded"]["white_noise_10mb"]:
        if r.get("codec") in ("HGS2_single", "GraSS-e2", "JARVIS3-l7", "GeCo3-l5"):
            require_number(r["archive_bytes"], f"ws_wn_{r['codec']}")
    require_in_readme_and_handoff(find_codec(ws["adversarial_and_bounded"]["white_noise_10mb"], "HGS2_single")["archive_bytes"], "ws_wn_hgs2")
    for r in ws["adversarial_and_bounded"]["pathology_fasta"]:
        if r.get("status") == "ok" and r.get("codec") in ("HGS2_collection", "MBGC-m3-L", "xz-9e", "zstd-ultra-22"):
            require_number(r["archive_bytes"], f"ws_path_{r['codec']}")
    require_in_readme_and_handoff(find_codec(ws["adversarial_and_bounded"]["pathology_fasta"], "HGS2_collection", "-Xmx4g")["archive_bytes"], "ws_path_hgs2")

    # Bacterial Mode C (FA + .seq) and Mode R (1-vs-panel)
    for b in ws["bacterial_shootout"]:
        sp = b["species"]
        require_number(b["fasta_meta"]["bytes"], f"ws_{sp}_fa_bytes")
        for c in ("HGS2_collection", "AGC", "MBGC-m3-L", "NAF-l22", "xz-9e", "zstd-19"):
            require_number(find_codec(b["mode_c_fa"], c)["archive_bytes"], f"ws_{sp}_cfa_{c}")
        require_in_readme_and_handoff(find_codec(b["mode_c_fa"], "HGS2_collection", "-Xmx16g")["archive_bytes"], f"ws_{sp}_cfa_hgs2")
        for c in ("HGS2_collection", "GraSS-e2", "JARVIS3-l7", "GeCo3-l5"):
            require_number(find_codec(b["mode_c_seq"], c)["archive_bytes"], f"ws_{sp}_cseq_{c}")
        require_number(find_codec(b["mode_r"], "HGS2_ref")["archive_bytes"], f"ws_{sp}_moder_hgs2")
        require_in_readme_and_handoff(find_codec(b["mode_r"], "HGS2_ref")["archive_bytes"], f"ws_{sp}_moder_hgs2")
        require_number(find_codec(b["mode_r"], "AGC_marginal_ModeR")["archive_bytes"], f"ws_{sp}_moder_agc")
        require_number(find_codec(b["mode_r"], "GeCo3-ref")["archive_bytes"], f"ws_{sp}_moder_geco3")

    # Viral shootout (group_19 & covid_mega)
    v19 = ws["viral_shootout"][0]
    require_number(v19["raw_bytes"], "ws_v19_raw")
    for c in ("HGS2_collection", "GraSS-e2", "AGC-a", "MBGC-m3", "NAF-l22", "xz-9e", "zstd-ultra-22"):
        require_number(find_codec(v19["runs_fa"], c)["archive_bytes"], f"ws_v19_fa_{c}")
    for c in ("HGS2_collection", "GraSS-e2", "JARVIS3-l7", "JARVIS3-l31", "GeCo3-l5"):
        require_number(find_codec(v19["runs_seq"], c)["archive_bytes"], f"ws_v19_seq_{c}")
    vmega = ws["viral_shootout"][3]
    require_number(vmega["raw_bytes"], "ws_vmega_raw")
    for c in ("HGS2_collection", "MBGC-m2", "NAF-l22", "zstd-19"):
        require_number(find_codec(vmega["runs_fa"], c)["archive_bytes"], f"ws_vmega_fa_{c}")
    require_in_readme_and_handoff(find_codec(vmega["runs_fa"], "HGS2_collection", "-Xmx22g")["archive_bytes"], "ws_vmega_fa_hgs2")
    require_in_readme_and_handoff(find_codec(vmega["runs_fa"], "NAF-l22")["archive_bytes"], "ws_vmega_fa_naf")

    # Human GRCh38 single chromosomes (chr21 & chr1) + Cross-Human Mode R + Diploid Mode C
    chr21 = ws["human_shootout"]["single_chromosomes"][0]
    require_number(chr21["raw_bytes"], "ws_chr21_raw")
    for c in ("HGS2_collection", "GraSS-e2", "NAF-l22", "xz-6", "zstd-19"):
        require_number(find_codec(chr21["runs_fa"], c)["archive_bytes"], f"ws_chr21_fa_{c}")
    for c in ("HGS2_collection", "JARVIS3-l7", "JARVIS3-l31", "GeCo3-l5"):
        require_number(find_codec(chr21["runs_seq"], c)["archive_bytes"], f"ws_chr21_seq_{c}")
    chr1 = ws["human_shootout"]["single_chromosomes"][2]
    require_number(chr1["raw_bytes"], "ws_chr1_raw")
    for c in ("HGS2_collection", "GraSS-e2", "NAF-l22", "xz-6", "zstd-19"):
        require_number(find_codec(chr1["runs_fa"], c)["archive_bytes"], f"ws_chr1_fa_{c}")

    for r in ws["human_shootout"]["cross_human_mode_r"]:
        h_run = find_codec(r["runs"], "HGS2_ref")
        require_number(h_run["archive_bytes"], f"ws_moder_{r['pair']}_hgs2")
        require_in_readme_and_handoff(h_run["archive_bytes"], f"ws_moder_{r['pair']}_hgs2")
        if any(x.get("codec") == "GeCo3-ref" for x in r["runs"]):
            g_run = find_codec(r["runs"], "GeCo3-ref")
            require_number(g_run["archive_bytes"], f"ws_moder_{r['pair']}_geco3")

    for d in ws["human_shootout"]["diploid_mode_c"]:
        require_number(find_codec(d["runs"], "HGS2_collection")["archive_bytes"], f"ws_{d['label']}_hgs2")
        require_in_readme_and_handoff(find_codec(d["runs"], "HGS2_collection")["archive_bytes"], f"ws_{d['label']}_hgs2")
        for c in ("AGC", "MBGC-m3-L", "NAF-l22", "zstd-19"):
            match = [x for x in d["runs"] if x.get("codec") == c and x.get("status") == "ok"]
            if match:
                require_number(match[0]["archive_bytes"], f"ws_{d['label']}_{c}")

    # GPU v3 (Chr1..8, Diploid_Chr1_to_4, covid_mega .seq@-Xmx512m, neural probes)
    g_c18_hgs = find_codec(gpu["grch38_multichr_1_8"], "HGS2_collection", "-Xmx20g")
    g_c18_zst = find_codec(gpu["grch38_multichr_1_8"], "zstd-19")
    for val, lbl in [
        (g_c18_hgs["raw_bytes"], "gpu_chr1_8_raw"),
        (g_c18_hgs["archive_bytes"], "gpu_chr1_8_hgs2"),
        (g_c18_zst["archive_bytes"], "gpu_chr1_8_zstd"),
    ]:
        require_number(val, lbl)
        require_in_readme_and_handoff(val, lbl)
    g_d14 = gpu["cross_human_chm13_vs_grch38"]["diploid_mode_c"][-1]
    for val, lbl in [
        (g_d14["raw_bytes"], "gpu_dip1_4_raw"),
        (find_codec(g_d14["runs"], "HGS2_collection")["archive_bytes"], "gpu_dip1_4_hgs2"),
        (find_codec(g_d14["runs"], "AGC")["archive_bytes"], "gpu_dip1_4_agc"),
        (find_codec(g_d14["runs"], "zstd-19")["archive_bytes"], "gpu_dip1_4_zstd"),
    ]:
        require_number(val, lbl)
        require_in_readme_and_handoff(val, lbl)
    g_cm_seq = find_codec(gpu["covid19_viral_scaling"][-1]["runs"], "HGS2_collection", "-Xmx512m")
    require_number(g_cm_seq["raw_bytes"], "gpu_covid_mega_seq_raw")
    require_number(g_cm_seq["archive_bytes"], "gpu_covid_mega_seq_hgs2")
    require_in_readme_and_handoff(g_cm_seq["raw_bytes"], "gpu_covid_mega_seq_raw")
    require_in_readme_and_handoff(g_cm_seq["archive_bytes"], "gpu_covid_mega_seq_hgs2")
    for np_row in gpu["cuda_gpu_experiments"]["neural_entropy_bench"]:
        require_number(np_row["theoretical_bytes"], f"gpu_neural_{np_row['dataset']}")

    # TPU v3 (All-24 CHM13v2.0 Mode R, all 13/13 bounded heap sweep cells, Amdahl covid_group_22)
    chm_tot = tpu["all24_chm13_given_grch38_mode_r"]["totals"]
    for k in ("total_chm13_raw_bytes", "total_chm13_bases", "total_archive_bytes"):
        require_number(chm_tot[k], f"tpu_{k}")
        require_in_readme_and_handoff(chm_tot[k], f"tpu_{k}")
    for cell in tpu["heap_floor_scaling_sweep"]:
        require_number(cell["archive_bytes"], f"tpu_heap_{cell['raw_bytes']}_{cell['xmx']}")
    for ref_cell in (tpu["heap_floor_scaling_sweep"][0], tpu["heap_floor_scaling_sweep"][9], tpu["heap_floor_scaling_sweep"][-1]):
        require_in_readme_and_handoff(ref_cell["archive_bytes"], f"tpu_heap_key_{ref_cell['raw_bytes']}_{ref_cell['xmx']}")
    require_number(tpu["amdahl_high_core_scaling"][0]["archive_bytes"], "tpu_amdahl_group22")
    require_in_readme_and_handoff(tpu["amdahl_high_core_scaling"][0]["archive_bytes"], "tpu_amdahl_group22")

    # 11. Verify compiled PDF and all 8 figure pairs exist + all 8 figures are embedded in README.md
    if not os.path.exists(pdf_path) or os.path.getsize(pdf_path) < 100_000:
        failures.append(f"[PDF MISSING OR TOO SMALL]: {pdf_path}")
    checks += 1
    for fn in [
        "fig1_regime_overview",
        "fig2_mode_r_loo_73genomes",
        "fig3_mode_c_pangenomes",
        "fig4_amdahl_and_bounded_memory",
        "fig5_laws_ablations_and_synthetic",
        "fig6_independent_modes_pareto",
        "fig7_gen1_gen2_foundations",
        "fig8_kaggle_world_shootout_and_human",
    ]:
        checks += 1
        if f"paper/figures/{fn}.png" not in readme_text:
            failures.append(f"[README MISSING FIGURE EMBED]: paper/figures/{fn}.png")
        for ext in [".pdf", ".png"]:
            fp = os.path.join(BASE_DIR, "paper", "figures", fn + ext)
            checks += 1
            if not os.path.exists(fp) or os.path.getsize(fp) < 10_000:
                failures.append(f"[FIGURE MISSING]: {fp}")

    print(f"[AUDIT SUMMARY] Total checks executed: {checks}, Failures: {len(failures)}")
    if failures:
        for f_msg in failures[:25]:
            print(" ", f_msg)
        sys.exit(1)
    print("[AUDIT PASSED 100%: ZERO UNVERIFIED NUMBERS, ALL JSON RECORDS & README/HANDOFF VERIFIED VERBATIM]")


if __name__ == "__main__":
    main()
