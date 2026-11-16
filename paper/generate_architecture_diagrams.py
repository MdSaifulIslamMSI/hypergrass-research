#!/usr/bin/env python3
"""
Generate 4 massive publication-grade architectural blueprints (300 DPI PNG + vector PDF)
for the HGS2 / EIMDL genomic compression research monograph and repository documentation:

  1. arch1_eimdl_system_pipeline     : Unified HGS2 / EIMDL 3-Regime Runtime & Disjoint Stage-1/2 Lifecycle
  2. arch2_mosaic_evolutionary_codec : Laws I-III & V: Two-Pass CSR Syncmer Index, Two-Gear Matcher, Indel Splice, V4 + Wobble
  3. arch3_fast_novelty_coder        : Laws IIb, IV & VI: 6F-ORF-HMM, 19/22-Expert FastCRM, Division-Free Mixing, Raw 2-Bit Bypass
  4. arch4_container_and_lineage     : Byte-Level .hgs2 v2.0 / .hgs v1.0 Binary Wire Specification & 3-Generation Evolution

Every symbol, constant, threshold, and metric matches the committed Java 21 source code and JSON reports verbatim.
"""

import os
import shutil
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon

BASE_DIR = r"C:\Users\mdsai\.gemini\antigravity\scratch\hypergrass"
FIG_DIR = os.path.join(BASE_DIR, "paper", "figures")
ARTIFACT_DIR = r"C:\Users\mdsai\.gemini\antigravity\brain\8fd6ada9-0a0f-403f-95e4-ee038b5c4fa3"

os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Cohesive high-contrast engineering palette (colorblind-friendly)
C_NAVY = "#0B4F6C"
C_NAVY_BG = "#F1F7FA"
C_TEAL = "#0F766E"
C_TEAL_BG = "#EFFAF7"
C_BLUE = "#1D4ED8"
C_BLUE_BG = "#EFF5FF"
C_VIOLET = "#6D28D9"
C_VIOLET_BG = "#F6F2FE"
C_AMBER = "#B45309"
C_AMBER_BG = "#FFF8EC"
C_CORAL = "#B91C1C"
C_CORAL_BG = "#FEF2F2"
C_SLATE = "#334155"
C_SLATE_BG = "#F4F6F9"
C_GREEN = "#15803D"
C_GREEN_BG = "#F0FDF4"
C_DARK = "#0F172A"
C_MUTED = "#475569"


def new_canvas(figsize=(15.4, 10.2)):
    fig = plt.figure(figsize=figsize, facecolor="white")
    ax = fig.add_axes([0.006, 0.006, 0.988, 0.988])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    return fig, ax


def save_arch(fig, name):
    pdf_path = os.path.join(FIG_DIR, f"{name}.pdf")
    png_path = os.path.join(FIG_DIR, f"{name}.png")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig)
    if os.path.exists(ARTIFACT_DIR):
        shutil.copy2(png_path, os.path.join(ARTIFACT_DIR, f"{name}.png"))
    print(f"[OK] Saved {pdf_path} and {png_path}")


def draw_title(ax, title, subtitle):
    bar = FancyBboxPatch(
        (0.6, 94.2), 98.8, 5.2,
        boxstyle="round,pad=0.0,rounding_size=0.6",
        facecolor=C_DARK, edgecolor="none", zorder=1
    )
    ax.add_patch(bar)
    ax.text(50.0, 97.55, title, ha="center", va="center",
            fontsize=13.4, fontweight="bold", color="white", zorder=3)
    ax.text(50.0, 95.3, subtitle, ha="center", va="center",
            fontsize=9.1, color="#CBD5E1", zorder=3)


def draw_zone(ax, x, y, w, h, title, border_color, bg_color, linestyle="-",
              badge=None, two_line_col=False):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.0,rounding_size=0.7",
        facecolor=bg_color, edgecolor=border_color, linewidth=1.6,
        linestyle=linestyle, zorder=1
    )
    ax.add_patch(box)
    if two_line_col:
        # Clean centered column header + subtitle pill (never overflows narrow columns)
        hdr = FancyBboxPatch(
            (x + 0.5, y + h - 3.6), w - 1.0, 3.1,
            boxstyle="round,pad=0.0,rounding_size=0.45",
            facecolor=border_color, edgecolor="none", zorder=2
        )
        ax.add_patch(hdr)
        ax.text(x + w / 2.0, y + h - 1.45, title, ha="center", va="center",
                fontsize=9.5, fontweight="bold", color="white", zorder=4)
        if badge:
            ax.text(x + w / 2.0, y + h - 2.9, badge, ha="center", va="center",
                    fontsize=7.7, color="#E2E8F0", family="monospace", zorder=4)
    else:
        ax.text(x + 1.1, y + h - 1.4, title, ha="left", va="center",
                fontsize=10.1, fontweight="bold", color=border_color, zorder=3)
        if badge:
            bw = min(w * 0.34, max(11.0, len(badge) * 0.44 + 2.0))
            bp = FancyBboxPatch(
                (x + w - bw - 0.8, y + h - 2.3), bw, 1.85,
                boxstyle="round,pad=0.0,rounding_size=0.4",
                facecolor=border_color, edgecolor="none", zorder=3
            )
            ax.add_patch(bp)
            ax.text(x + w - bw / 2 - 0.8, y + h - 1.38, badge,
                    ha="center", va="center", fontsize=7.9, fontweight="bold",
                    color="white", zorder=4)


def draw_card(ax, x, y, w, h, header, body_lines, color, bg="white",
              header_h=2.35, title_size=9.1, body_size=8.1, mono_footer=None, center_body=False):
    outer = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.0,rounding_size=0.55",
        facecolor=bg, edgecolor=color, linewidth=1.35, zorder=2
    )
    ax.add_patch(outer)
    hdr = FancyBboxPatch(
        (x, y + h - header_h), w, header_h,
        boxstyle="round,pad=0.0,rounding_size=0.55",
        facecolor=color, edgecolor=color, linewidth=1.0, zorder=3
    )
    ax.add_patch(hdr)
    ax.add_patch(plt.Rectangle((x, y + h - header_h), w, min(0.65, header_h * 0.35),
                               facecolor=color, edgecolor="none", zorder=3))
    ax.text(x + w / 2.0, y + h - header_h / 2.0, header,
            ha="center", va="center", fontsize=title_size, fontweight="bold",
            color="white", zorder=4)

    body_top = y + h - header_h - 0.4
    body_bot = y + (2.25 if mono_footer else 0.35)
    n = len(body_lines)
    if n > 0:
        avail = max(1.0, body_top - body_bot)
        step = avail / n
        for idx, line in enumerate(body_lines):
            cy = body_top - (idx + 0.5) * step
            if center_body:
                ax.text(x + w / 2.0, cy, line, ha="center", va="center",
                        fontsize=body_size, color=C_DARK, zorder=4)
            else:
                ax.text(x + 0.72, cy, line, ha="left", va="center",
                        fontsize=body_size, color=C_DARK, zorder=4)

    if mono_footer:
        fp = FancyBboxPatch(
            (x + 0.4, y + 0.32), w - 0.8, 1.7,
            boxstyle="round,pad=0.0,rounding_size=0.35",
            facecolor="#F1F5F9", edgecolor="#CBD5E1", linewidth=0.8, zorder=3
        )
        ax.add_patch(fp)
        ax.text(x + w / 2.0, y + 1.17, mono_footer,
                ha="center", va="center", fontsize=7.5, family="monospace",
                fontweight="bold", color=color, zorder=4)


def draw_diamond(ax, cx, cy, hw, hh, label, color, bg="#FFFBEB", fontsize=8.0):
    pts = [(cx, cy + hh), (cx + hw, cy), (cx, cy - hh), (cx - hw, cy)]
    poly = Polygon(pts, closed=True, facecolor=bg, edgecolor=color, linewidth=1.45, zorder=4)
    ax.add_patch(poly)
    ax.text(cx, cy, label, ha="center", va="center", fontsize=fontsize,
            fontweight="bold", color=C_DARK, zorder=5)


def draw_arrow(ax, x1, y1, x2, y2, color=C_SLATE, label=None, rad=0.0,
               lw=1.55, ls="-", lpos=0.5, lx_off=0.0, ly_off=0.0, label_bg="white"):
    arr = FancyArrowPatch(
        (x1, y1), (x2, y2),
        connectionstyle=f"arc3,rad={rad}",
        arrowstyle="-|>,head_length=5.5,head_width=3.2",
        color=color, linewidth=lw, linestyle=ls, zorder=6
    )
    ax.add_patch(arr)
    if label:
        mx = x1 + (x2 - x1) * lpos + lx_off
        my = y1 + (y2 - y1) * lpos + ly_off
        ax.text(
            mx, my, label, ha="center", va="center", fontsize=7.5, fontweight="bold",
            color=color, zorder=8,
            bbox=dict(boxstyle="round,pad=0.18", facecolor=label_bg, edgecolor=color, linewidth=0.75)
        )


def draw_poly_arrow(ax, pts, color=C_SLATE, label=None, lpt_idx=0, lx_off=0.0, ly_off=0.0, lw=1.55, ls="-"):
    """Draw multi-segment orthogonal polyline ending in a crisp arrowhead without clipping text."""
    for i in range(len(pts) - 2):
        ax.plot([pts[i][0], pts[i + 1][0]], [pts[i][1], pts[i + 1][1]],
                color=color, linewidth=lw, linestyle=ls, zorder=6, solid_capstyle="round")
    p_a, p_b = pts[-2], pts[-1]
    arr = FancyArrowPatch(
        p_a, p_b,
        arrowstyle="-|>,head_length=5.5,head_width=3.2",
        color=color, linewidth=lw, linestyle=ls, zorder=6
    )
    ax.add_patch(arr)
    if label:
        pa, pb = pts[lpt_idx], pts[lpt_idx + 1]
        mx = 0.5 * (pa[0] + pb[0]) + lx_off
        my = 0.5 * (pa[1] + pb[1]) + ly_off
        ax.text(
            mx, my, label, ha="center", va="center", fontsize=7.5, fontweight="bold",
            color=color, zorder=8,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=color, linewidth=0.75)
        )


def draw_footer_strip(ax, items, y=0.6, h=5.2):
    strip = FancyBboxPatch(
        (0.6, y), 98.8, h,
        boxstyle="round,pad=0.0,rounding_size=0.55",
        facecolor="#0F172A", edgecolor=C_NAVY, linewidth=1.2, zorder=2
    )
    ax.add_patch(strip)
    w_each = 98.8 / len(items)
    for i, (kpi, desc) in enumerate(items):
        cx = 0.6 + (i + 0.5) * w_each
        if i > 0:
            ax.plot([0.6 + i * w_each, 0.6 + i * w_each], [y + 0.7, y + h - 0.7],
                    color="#334155", linewidth=1.0, zorder=3)
        ax.text(cx, y + h * 0.66, kpi, ha="center", va="center",
                fontsize=9.4, fontweight="bold", color="#38BDF8", zorder=4)
        ax.text(cx, y + h * 0.28, desc, ha="center", va="center",
                fontsize=7.5, color="#E2E8F0", zorder=4)


# ==============================================================================
# DIAGRAM 1: END-TO-END UNIFIED HGS2 / EIMDL ENGINE ARCHITECTURE
# ==============================================================================
def make_arch1():
    fig, ax = new_canvas((15.4, 10.3))
    draw_title(
        ax,
        "HGS2 / EIMDL Unified System Architecture & Strictly Disjoint Two-Stage Lifecycle",
        "End-to-End Bijective Ingestion (Phi), Three-Regime Router (Mode S / R / C), Zero-Materialization CSR Stage 1, Explicit Memory Barrier & Parallel Stage 2"
    )

    # Band 0: Input & 6-Channel Bijective Factorization (y=75.8..93.2)
    draw_zone(ax, 0.8, 75.8, 98.4, 17.4,
              "0. Streaming Input Ingestion & Six-Channel Exact Bijection (Phi / Phi^-1)",
              C_SLATE, C_SLATE_BG, badge="CRC-32 + 100% SHA-256")
    draw_card(ax, 2.0, 76.9, 17.0, 11.4, "Raw Sequence Input", [
        "• Flat ACGT / ACGU (.seq)",
        "• FASTA / Multi-FASTA",
        "• LF or CRLF + soft-mask",
        "• IUPAC + chimeric T/U",
        "• Streaming IEEE CRC-32",
    ], C_SLATE, mono_footer="InputStream / Byte[]")

    draw_diamond(ax, 27.2, 82.6, 5.2, 4.1, "Flat ACGT/U\n<= 2^31-65?", C_TEAL, "#ECFDF5", fontsize=7.8)

    draw_card(ax, 35.6, 76.9, 21.8, 11.4, "tryStreamFlatFile (64 KiB)", [
        "• Direct 64 KiB I/O chunk scan",
        "• Zero full-file byte[] allocation",
        "• FLAG_RAW_FLAT_ACGT (0x04)",
        "• FLAG_RAW_RNA_U (0x08)",
        "• Metadata overhead: 0 Bytes",
    ], C_TEAL, mono_footer="PageBuilder -> 256 KiB Pages")

    draw_card(ax, 60.0, 76.9, 37.8, 11.4, "FastaStreamFactorizer (6-Channel Bijection)", [
        "• Ch 1: Front-coded + Deflate FASTA headers & line deltas",
        "• Ch 2: RLE sequence line widths + LF/CRLF/NONE markers",
        "• Ch 3: Alternating uppercase / lowercase soft-mask intervals",
        "• Ch 4: Scaffold 'N'-run gap/len + Ch 5: IUPAC / T+U exceptions",
        "• Ch 6: 2-Bit Bank (A=00, C=01, G=10, T/U=11; 256 KiB < 512 KiB G1GC)",
    ], C_NAVY, mono_footer="Packed2BitSequence (0.25 B/bp)")

    # Clean dedicated upper routing lane at y=89.6 over tryStreamFlatFile (top=88.3)
    draw_arrow(ax, 19.0, 82.6, 22.0, 82.6, C_SLATE)
    draw_arrow(ax, 32.4, 82.6, 35.6, 82.6, C_TEAL, "Yes", lpos=0.45, ly_off=1.35)
    draw_poly_arrow(ax, [(27.2, 86.7), (27.2, 89.6), (78.9, 89.6), (78.9, 88.3)],
                    C_NAVY, label="No (FASTA / Mixed)", lpt_idx=1, ly_off=0.0)

    # Band 1: Three-Regime Conditioning & Causal Ordering (y=57.2..74.0)
    draw_zone(ax, 0.8, 57.2, 98.4, 16.8,
              "1. Strict Three-Regime Information Router (mode in {0, 1, 2})",
              C_BLUE, C_BLUE_BG, badge="EimdlCli / FusionEngine")

    draw_card(ax, 2.0, 58.3, 29.0, 11.5, "Mode S: Single Genome (mode=0)", [
        "• Reference R = empty (0 B pin)",
        "• Target self-indexed as donor 0",
        "• Bidirectional (+/-) self-homology",
        "• Tandem periods P in [1..64]",
        "• Feeds Mode-S 19-expert FastCRM",
    ], C_TEAL, mono_footer="Single-Genome Profile (<= 6M)")

    draw_card(ax, 33.2, 58.3, 31.8, 11.5, "Mode R: Reference Panel (mode=1)", [
        "• Loads external --panel donors D[0..M-1]",
        "• Header stores Deflate donor name table",
        "• Mandatory 32-B SHA-256 Merkle manifest",
        "• Decoder aborts on 1-bp panel mismatch",
        "• Targets conditioned on all M donors",
    ], C_BLUE, mono_footer="verifyAndLoadPanelForDecode()")

    draw_card(ax, 67.2, 58.3, 30.6, 11.5, "Mode C: Pangenome Collection (mode=2)", [
        "• Bottom-2048 canonical open syncmers",
        "• Pairwise Jaccard -> Medoid root r = 0",
        "• Prim MST orders ranks r = 0 .. U-1",
        "• Strict causal cap maxAncestorRank[r] = r",
        "• Similarity list passed to Stage-1 DAG",
    ], C_VIOLET, mono_footer="buildPhylogeneticPrimOrder()")

    # Inter-band arrows stop neatly at Band 1 border (avoiding zone titles completely)
    draw_arrow(ax, 46.5, 76.9, 46.5, 74.0, C_TEAL, lw=1.6)
    draw_arrow(ax, 78.9, 76.9, 78.9, 74.0, C_NAVY, lw=1.6)

    # Band 2: STAGE 1 (y=35.4..55.4)
    draw_zone(ax, 0.8, 35.4, 98.4, 20.0,
              "STAGE 1: Zero-Materialization CSR Syncmer Index & Parallel Mosaic Emission (Laws I–III, V)",
              C_NAVY, C_NAVY_BG, badge="Residency: M_Stage1")

    draw_card(ax, 2.0, 36.6, 34.4, 15.4, "EimdlSyncmerIndex (2-Pass CSR)", [
        "• Canonical open syncmers: k=16, s=10, w=6, t=3",
        "• Pass 1: count into csrHeads[<=1,048,576] (4 MiB)",
        "• Prefix sum -> exact in-place bucket offsets",
        "• Pass 2: fill csrPositions[M] (4B) + csrTags[M] (2B)",
        "• Subsample shift s in {0, 1, 2} from JVM maxMem:",
        "   -Xmx16g -> s=0 (1.58 B/bp) | -Xmx256m -> s=1 (0.93)",
        "   -Xmx128m / ultraTight (<=205 MiB) -> s=2 (0.59 B/bp)",
    ], C_NAVY, mono_footer="int[] csrHeads + csrPos + short[] tags")

    draw_card(ax, 39.2, 36.6, 36.8, 15.4, "EimdlMosaicCodec (Parallel Rank DAG Workers)", [
        "• Workers throttled: extremeTight<=145M -> 1T; ultraTight -> <=2T",
        "• Target syncmers narrowed long[] -> int[] (-38% worker RAM)",
        "• Priority probe: MTF[0..1] -> MST parent -> CSR bucket (48 probes)",
        "• Two-Gear matcher: Clonal Gear 1 vs Diverged Gear 2 + Horizon (384 bp)",
        "• 2D Micro-Indel Splice [1..256] x [1..64] (>=12 bp anchor)",
        "• Law I (V4 Kimura XOR delta) + Law IIa (gap+1 mod 3 wobble)",
        "• BitRangeEnc (scale 4096, tau=20 Taper -> p_min = 11/4096)",
    ], C_TEAL, mono_footer="ForkJoinPool / Intra-Unit Chunks")

    draw_card(ax, 78.2, 36.6, 19.6, 15.4, "Stage-1 Outputs", [
        "• rankMosaics[0..U-1]",
        "  (order-1 range streams)",
        "",
        "• unitNovSlices[0..U-1]",
        "  (unmatched literal +",
        "   insertion bases)",
    ], C_VIOLET, mono_footer="EncodedMosaic[]")

    draw_arrow(ax, 16.5, 57.2, 16.5, 55.4, C_TEAL, lw=1.4)
    draw_arrow(ax, 49.1, 57.2, 49.1, 55.4, C_BLUE, lw=1.4)
    draw_arrow(ax, 82.5, 57.2, 82.5, 55.4, C_VIOLET, lw=1.4)
    draw_arrow(ax, 36.4, 44.3, 39.2, 44.3, C_NAVY, label="CSR", ly_off=1.35)
    draw_arrow(ax, 76.0, 44.3, 78.2, 44.3, C_TEAL)

    # Barrier Strip: Law VII Disjoint Deallocation (y=29.6..34.0)
    barrier = FancyBboxPatch(
        (0.8, 29.6), 98.4, 4.3,
        boxstyle="round,pad=0.0,rounding_size=0.5",
        facecolor=C_AMBER_BG, edgecolor=C_AMBER, linewidth=2.0, linestyle="--", zorder=3
    )
    ax.add_patch(barrier)
    ax.text(50.0, 32.5,
            "LAW VII DISJOINT LIFECYCLE BARRIER:  sharedDagIdx = null;  mstParentRank = null;  rankDonorSimilarityOrder = null;",
            ha="center", va="center", fontsize=9.1, fontweight="bold", family="monospace", color=C_AMBER, zorder=4)
    ax.text(50.0, 30.65,
            "Stage-1 CSR Index & 2-Bit Bank Freed Before Stage 2 Starts  ->  M_peak = max(M_Stage1, M_Stage2) << M_Stage1 + M_Stage2  |  rankMosaics[r] = null streamed incrementally",
            ha="center", va="center", fontsize=8.2, fontweight="bold", color=C_DARK, zorder=4)

    draw_arrow(ax, 88.0, 36.6, 88.0, 33.9, C_AMBER, lw=1.7)

    # Band 3: STAGE 2 (y=6.8..28.2)
    draw_zone(ax, 0.8, 6.8, 98.4, 21.4,
              "STAGE 2: Unified Pan-Genome Novelty Superstream (FastNoveltyCoder) & .hgs2 Emission (Laws IIb, IV, VI)",
              C_VIOLET, C_VIOLET_BG, badge="Residency: M_Stage2 (~1.1-6.8 MB/worker)")

    draw_card(ax, 2.0, 8.1, 22.8, 16.4, "Novelty Superstream", [
        "• Gathers unitNovSlices[0..U-1]",
        "• Frees slices once concatenated",
        "• Auto-selects capacity profile:",
        "  <= 6 Mbp or bounded:",
        "    19 experts (87 KiB ord0..7)",
        "  > 6 Mbp PanSuper:",
        "    22 experts (+ord8/9 + sparse12)",
    ], C_VIOLET, mono_footer="totalNoveltyBases + Profile")

    draw_card(ax, 27.0, 8.1, 30.0, 16.4, "Warm-Primed FastCrmModel Stripes", [
        "• Seg 0 populates immutable WarmSnapshot",
        "• 6F-ORF-HMM tracks {+0,+1,+2,-0,-1,-2} frames",
        "• Hierarchical Pitman-Yor suffix backoff chain",
        "• Dual combiner: LN/EXP Bayes + weights[8][E]",
        "• Division-free RECIP[65600] normalization",
        "• 1,536-bin SSE -> 15-bit carryless RangeCoder",
    ], C_NAVY, mono_footer="4 Mbp Stripes | ForkJoinPool")

    draw_diamond(ax, 65.5, 16.3, 6.1, 5.1, "Law IV Gate\n|C_j| >= ceil(n_j/4)?", C_CORAL, C_CORAL_BG, fontsize=7.8)

    draw_card(ax, 74.0, 17.0, 23.8, 7.5, "SEG_RAW_2BIT (Mode 0) | Safety", [
        "• Incompressible / white-noise guard",
        "• Raw 2-bit packed: ceil(segLen/4) B",
        "• Bounds crypto noise <= 2.0012 bpb",
    ], C_CORAL, body_size=8.0)

    draw_card(ax, 74.0, 8.1, 23.8, 7.7, ".hgs2 Binary Container (v2.0)", [
        "• Writes 0x48475332 + Mode streams",
        "• Verified lossless Phi^-1 (64 KiB buf)",
        "• Bit-exact SHA-256 FASTA/RNA match",
    ], C_GREEN, body_size=8.0)

    # Clean non-overlapping drop at right gutter (x=76..88) into Stage 2 border
    draw_arrow(ax, 88.0, 29.6, 88.0, 28.2, C_VIOLET, lw=1.7)
    draw_arrow(ax, 24.8, 16.3, 27.0, 16.3, C_VIOLET)
    draw_arrow(ax, 57.0, 16.3, 59.4, 16.3, C_NAVY)
    draw_arrow(ax, 69.8, 19.3, 74.0, 20.8, C_CORAL, "Yes")
    draw_arrow(ax, 69.8, 13.3, 74.0, 12.0, C_GREEN, "No (CRM)")
    draw_arrow(ax, 85.9, 17.0, 85.9, 15.8, C_SLATE)

    draw_footer_strip(ax, [
        ("1,174,113 B (102.30x)", "Mode R E. coli 24g LOO (7.14x < AGC)"),
        ("47,339 B (2,233.78x)", "Mode R MTB 24g LOO (2.04x < AGC Oracle)"),
        ("3,371,489 / 3,388,639 B", "Mode C E. coli 1T / 12T (19.65s, 45.2x < v1)"),
        ("3,452,703 B @ -Xmx128m", "Pure JVM 1.06 B/bp on 120.11 Mbp"),
        ("<= 2.0012 bpb | 100% SHA-256", "Crypto White Noise & 1,705/1,705 Audited"),
    ])

    save_arch(fig, "arch1_eimdl_system_pipeline")


# ==============================================================================
# DIAGRAM 2: MOSAIC EVOLUTIONARY ENGINE MICRO-ARCHITECTURE (LAWS I-III, V, VIa)
# ==============================================================================
def make_arch2():
    fig, ax = new_canvas((15.4, 10.2))
    draw_title(
        ax,
        "Evolutionary-Information Mosaic Engine Micro-Architecture (Laws I–III, V & VIa)",
        "Two-Pass CSR Open-Syncmer Indexing, Priority Donor Probing, Horizon-Bounded Two-Gear Coalescent Matcher, 2D Micro-Indel Splicing & V4 / Wobble Coding"
    )

    # Wide gutters (x=31.0..34.4 and x=66.4..69.8) so arrows & pills never touch box borders
    draw_zone(ax, 0.8, 7.2, 30.2, 85.6,
              "A. Open-Syncmer CSR & Probing", C_NAVY, C_NAVY_BG,
              badge="Laws IIIa, V  |  0.50-1.25 B/bp", two_line_col=True)

    draw_card(ax, 1.8, 71.2, 28.2, 16.5, "Canonical Open Syncmers (k=16, s=10)", [
        "• Window w = k - s_mer = 6 overlapping 10-mers",
        "• Canonical 16-mer hash: min(h_fwd, h_rc)",
        "• Selected iff min 10-mer falls at center t = 3",
        "• Exact base density = 1 / (6 * 2^s_sub)",
        "• Target syncmers narrowed: long[] -> int[] (-38%)",
    ], C_NAVY, mono_footer="EimdlSyncmerIndex.java")

    draw_card(ax, 1.8, 46.8, 28.2, 21.6, "Two-Pass Zero-Materialization CSR Index", [
        "• Pass 1: Hash into csrHeads[2^B] (<= 1M -> 4 MiB)",
        "• In-place prefix sum converts counts -> offsets",
        "• Pass 2: Write positions & tags directly into CSR:",
        "   - csrPositions[M] : 4 Bytes/entry (32-bit pos)",
        "   - csrTags[M]      : 2 Bytes/entry (16-bit verify)",
        "• Bank + CSR = N * (0.25 + 2^-s) + 4 MiB",
        "• Zero intermediate ArrayList/Entry object GC",
    ], C_TEAL, mono_footer="s=0: 1.25 | s=1: 0.75 | s=2: 0.50 B/bp")

    draw_card(ax, 1.8, 8.6, 28.2, 35.4, "Strict Candidate Probe Priority Chain", [
        "Priority 1 — Sticky Collinear Donors (mtf[0..1]):",
        "  • Tracks lastRefEnd[ds] + (pos - lastTgtEnd[ds])",
        "  • Keeps clonal frame locked past minor gaps/SNPs",
        "Priority 2 — Causal MST Parent & Siblings:",
        "  • Explores Prim MST parent before general buckets",
        "Priority 3 — Distinct-Donor CSR Bucket Scan:",
        "  • Up to MAX_BUCKET_PROBES = 48 unique donors",
        "  • Prevents multi-copy IS repeats flooding bucket",
        "Priority 4 — Law IV Local Self & Tandem Tiers:",
        "  • Bidirectional (+/-) self-homology ring",
        "  • Overlapping tandem periods P in [1 .. 64]",
    ], C_BLUE, mono_footer="probePanelStart() / probePanelBucket()")

    draw_arrow(ax, 15.9, 71.2, 15.9, 68.4, C_NAVY)
    draw_arrow(ax, 15.9, 46.8, 15.9, 44.0, C_TEAL)

    # Column 2 (x=34.4..66.4): Two-Gear Coalescent Matcher & Micro-Indel Splice
    draw_zone(ax, 34.4, 7.2, 32.0, 85.6,
              "B. Two-Gear Coalescent & Splice Engine", C_TEAL, C_TEAL_BG,
              badge="Laws IIIb-c  |  EimdlMosaicCodec", two_line_col=True)

    draw_card(ax, 35.6, 69.2, 29.6, 18.5, "Gear 1: Strict Clonal Haplotype Copying", [
        "• clonalGear = true  (vertical clonal lineage)",
        "• Scoring: baseGain = +0.008 b/bp, snpPenalty = -11.2 b",
        "• xDrop = -16.0 b, exact >= 128 bp tail guard",
        "• Rejects consecMm >= 2 or windowMm > 0 (512 bp decay)",
        "• 64-bit SWAR word XOR (countExactMatch 32 bases/op)",
    ], C_TEAL, mono_footer="Fast Exit if exact span >= 2,048 bp")

    draw_diamond(ax, 50.4, 59.8, 12.2, 5.5,
                 "Burst >1 SNP or short (<96 bp)?\nTrigger findClonalReentryHorizon()",
                 C_AMBER, C_AMBER_BG, fontsize=7.9)

    draw_card(ax, 35.6, 32.2, 29.6, 18.8, "Gear 2: Horizon-Bounded Diverged Synteny", [
        "• Probes 384 bp (4 x 96 bp) ahead for >= 192 bp clonal match",
        "• Clamps extension to [pos .. gear2Limit) strictly",
        "• Scoring: baseGain = +1.865, snpPenalty = -7.8, xDrop = -22.0",
        "• Allows consecMm < 4 & windowMm < 28 (snp*7 <= k)",
        "• Halts recombination bulldozing -> exact clonal reentry",
    ], C_AMBER, mono_footer="E. coli 00.seq: 29,623 -> 5,706 SNPs")

    draw_card(ax, 35.6, 8.6, 29.6, 20.6, "Law IIIc: 2D Affine Micro-Indel Splicing", [
        "• Scans downstream 2D grid upon match termination:",
        "   Delta_ref in [1 .. 256]  x  Delta_tgt in [1 .. 64]",
        "• Requires MIN_SPLICE_ANCHOR >= 12 bp continuation",
        "• Encodes splice offset in ~10-18 b vs 60-80 b new token",
        "• hasSplice conditioned on 8-state (kind, spanBin) prior",
        "• Chains up to MAX_SPLICE_CHAIN = 1,024 sub-blocks",
    ], C_VIOLET, mono_footer="attemptMicroIndelSplice()")

    # Orthogonal route in wide x=32.7 gutter
    draw_poly_arrow(ax, [(30.0, 26.0), (32.7, 26.0), (32.7, 78.5), (35.6, 78.5)],
                    C_BLUE, label="Cands", lpt_idx=1, ly_off=0.0)
    draw_arrow(ax, 50.4, 69.2, 50.4, 65.3, C_TEAL)
    draw_arrow(ax, 50.4, 54.3, 50.4, 51.0, C_AMBER, "Clamp [pos..horizon)")
    draw_arrow(ax, 50.4, 32.2, 50.4, 29.2, C_VIOLET, "At Match End")

    # Column 3 (x=69.8..99.2): Evolutionary Entropy Coders
    draw_zone(ax, 69.8, 7.2, 29.4, 85.6,
              "C. Evolutionary Entropy & Range Coding", C_VIOLET, C_VIOLET_BG,
              badge="Laws I, IIa, VIa  |  BitRangeEnc", two_line_col=True)

    draw_card(ax, 71.0, 59.8, 27.0, 27.9, "Law I: Klein V4 Kimura-CpG SNPs", [
        "• Packed 2-bit: A=00, C=01, G=10, T=11",
        "• Complement b_bar = b XOR 3 (11_2)",
        "• Diff delta = tgt XOR ref in {1, 2, 3}:",
        "   delta = 2 (10): Transition (A<->G, C<->T)",
        "   delta = 1 (01): Keto-Amino Transversion",
        "   delta = 3 (11): WC-Complement Transversion",
        "• Strand-Invariant Proof:",
        "   (tgt^3)^(ref^3) = tgt^ref = delta",
        "• Context: CpG (prevRef, ref) + mod3Zero",
    ], C_NAVY, mono_footer="1.10-1.32 b/SNP (vs 1.585 uniform)")

    draw_card(ax, 71.0, 33.4, 27.0, 23.2, "Law IIa: Bijective Wobble Gap Split", [
        "• Forward map phi(gap) = (q, phi):",
        "   q   = floor(gap / 3)",
        "   phi = (gap + 1) mod 3  in {0, 1, 2}",
        "• Inverse: gap = 3*q + ((phi + 2) mod 3)",
        "• Zero overhead if non-coding (uniform phi)",
        "• Captures 3rd-codon wobble (gap=2 mod 3)",
    ], C_TEAL, mono_footer="phi = 0 Resonance under dN/dS << 1")

    draw_card(ax, 71.0, 8.6, 27.0, 21.6, "Law VIa: Two-Adic Taper Range Coder", [
        "• Binary adaptive probability scale M = 4096",
        "• Classical truncation stalls at 31/4096",
        "• Law VIa symmetric tau = 20 rounding:",
        "   b=1: p += (4096 - p + 20) >> 5",
        "   b=0: p -= (p + 20) >> 5",
        "• Fixed-point p_min = 31 - 20 = 11/4096",
    ], C_VIOLET, mono_footer="2.83x Lower Floor: 0.00388 b/flag")

    # Orthogonal routes in wide gutter x=66.4..69.8 (non-overlapping vertical lanes)
    draw_poly_arrow(ax, [(65.2, 46.5), (68.5, 46.5), (68.5, 74.0), (71.0, 74.0)],
                    C_NAVY, label="SNPs", lpt_idx=1, ly_off=0.0)
    draw_poly_arrow(ax, [(65.2, 38.5), (67.3, 38.5), (67.3, 45.0), (71.0, 45.0)],
                    C_TEAL, label="Gaps", lpt_idx=0, lx_off=0.0, ly_off=-1.6)
    draw_arrow(ax, 84.5, 59.8, 84.5, 56.6, C_NAVY)
    draw_arrow(ax, 84.5, 33.4, 84.5, 30.2, C_TEAL)

    draw_footer_strip(ax, [
        ("-6,177 B (-3.13%)", "Law I Ablation: Distant E. coli 00|20"),
        ("-1,861 B Wobble", "Law IIa Codon Split on E. coli 00|20"),
        ("-16,225 B (-23,917 SNPs)", "Law IIIb Two-Gear on E. coli 00.seq LOO"),
        ("261 -> 18 SNPs (256 B)", "MTB 03.seq LOO vs Gen 2 Single-Gear"),
        ("-32.3% (6,818 vs 10,067 B)", "Law IIIc Micro-Indel Splice (Hypermutator)"),
    ])

    save_arch(fig, "arch2_mosaic_evolutionary_codec")


# ==============================================================================
# DIAGRAM 3: FAST NOVELTY CODER MICRO-ARCHITECTURE (LAWS IIb, IV, VI)
# ==============================================================================
def make_arch3():
    fig, ax = new_canvas((15.4, 10.2))
    draw_title(
        ax,
        "FastNoveltyCoder & FastCrmModel Micro-Architecture (Laws IIb, IV & VI)",
        "Causal 6-Frame ORF-HMM, Adaptive 19/22-Expert Pitman-Yor Hierarchy, Division-Free Log-Bayes & Logistic Mixing, SSE & Law IV Raw 2-Bit Safety Bypass"
    )

    # Column 1 (x=0.8..24.2)
    draw_zone(ax, 0.8, 7.2, 23.4, 85.6,
              "A. Causal Features & Rings", C_NAVY, C_NAVY_BG,
              badge="O(1) / base  |  FastNoveltyCoder", two_line_col=True)

    draw_card(ax, 1.7, 71.8, 21.6, 16.0, "64-Bit Rolling Histories", [
        "• fwdHist: 32-base shift window",
        "• rcHist: reverse-complement mirror",
        "• Auto Capacity Switch:",
        "  <=6M or bounded: 19 experts",
        "  >6M panSuper: 22 experts (+ord8/9,sp12)",
    ], C_NAVY, mono_footer="L2/L3-Resident vs PanSuper")

    draw_card(ax, 1.7, 39.2, 21.6, 29.6, "Law IIb: 6-Frame Causal ORF-HMM", [
        "• Tracks 6 reading frames simultaneously:",
        "  Fwd: +0, +1, +2  |  Rev: -0, -1, -2",
        "• Stop priors at phase c = 2:",
        "  + stops: TAA, TAG, TGA (ctx2=12,14)",
        "  - stops: TTA, CTA, TCA (ctx2=7,13,15)",
        "• Updates integer log-posteriors frameLog[6]",
        "• Crossing gene boundaries resets frame,",
        "  routes bestFrame -> activeRegime / codon",
    ], C_TEAL, mono_footer="orfPlus2/3/5 + orfMinus2/3/5")

    draw_card(ax, 1.7, 8.6, 21.6, 27.6, "Tag-Verified Repeat Seed Rings", [
        "• Exact 32-bit tag + 32-bit pos packing:",
        "  repFineFwd/Rev (K_FINE=13, 2-way)",
        "  repLongFwd/Rev (K_LONG=18, 2-way)",
        "• Zero false positives (no scan loops)",
        "• Spawns 6 MicroRepeatHead[0..5] cursors",
        "• Post-SNP instantaneous confidence restore",
    ], C_BLUE, mono_footer="[tag32 | pos32] -> 6 Repeat Heads")

    draw_arrow(ax, 12.5, 71.8, 12.5, 68.8, C_NAVY)
    draw_arrow(ax, 12.5, 39.2, 12.5, 36.2, C_TEAL)

    # Column 2 (x=26.6..52.2)
    draw_zone(ax, 26.6, 7.2, 25.6, 85.6,
              "B. 22-Expert Pitman-Yor Bank", C_VIOLET, C_VIOLET_BG,
              badge="Law VIb  |  Suffix Summary Backoff", two_line_col=True)

    draw_card(ax, 27.5, 54.6, 23.8, 33.2, "Dense & Sparse Markov Suffix Chain", [
        "• ord0 [4] .. ord7 [16,384*4] (87 KiB dense)",
        "• ord8 [65,536 int] + ord9 [262,144 int] (panSuper)",
        "• sparse12 [2^18 long] open-addressed 12-mer",
        "• Law VIb Hierarchical Pitman-Yor Backoff:",
        "   P_k(a|c_k) = (n_k + alpha_k * P_{k-1}) / (N_k + alpha_k)",
        "• Pre-evaluated lower-order summary prior P_{k-1}",
        "• Prevents low-count dilution vs flat 1/4 prior",
    ], C_VIOLET, mono_footer="Experts 0..7 + 13..15 (ord0..sparse12)")

    draw_card(ax, 27.5, 31.8, 23.8, 19.8, "Codon-Phase & Spaced-Seed Experts", [
        "• codon1, codon3, codon5 (Experts 8, 9, 10)",
        "  marginalized over active 6F-ORF-HMM frameLog",
        "• wobble8 [65,536*4] (Expert 11: mod-3 stride)",
        "• gapped8 [65,536*4] (Expert 12: spaced 9-mer)",
    ], C_TEAL, mono_footer="Experts 8..12")

    draw_card(ax, 27.5, 8.6, 23.8, 20.2, "6 Multi-Scale Kimura-XOR Repeat Heads", [
        "• Forward (+) & reverse-complement (-) cursors",
        "• Predicts next 2-bit base & tracks misses",
        "• On mismatch, switches smoothly without reset",
        "• When misses == 0, recovers high repeat gain",
    ], C_BLUE, mono_footer="Experts 13..18 / 16..21 (heads[0..5])")

    draw_arrow(ax, 23.3, 79.5, 27.5, 76.0, C_NAVY)
    draw_arrow(ax, 23.3, 53.5, 27.5, 42.0, C_TEAL)
    draw_arrow(ax, 23.3, 22.0, 27.5, 19.0, C_BLUE)

    # Column 3 (x=54.8..75.4): keep right-side vertical channel at x=73.6..74.8 for Log-Bayes -> SSE!
    draw_zone(ax, 54.8, 7.2, 20.6, 85.6,
              "C. Mixers, SSE & RangeCoder", C_NAVY, C_NAVY_BG,
              badge="Law VIc  |  Zero Hardware idiv", two_line_col=True)

    draw_card(ax, 55.7, 67.2, 18.8, 20.6, "Log-Bayes Probability Pool", [
        "• LN_LUT[16385] (Q12 ln)",
        "• EXP_NEG_LUT[65537]",
        "• Tracks bayesScores[E]",
        "• Emits bayesP[0..3]",
    ], C_NAVY, mono_footer="metaBayes Weighting")

    # Narrower Logistic Mixer so Log-Bayes routes down right gutter (x=74.3) cleanly outside Mixer!
    draw_card(ax, 55.7, 43.2, 16.4, 21.0, "8-Regime Logistic Mixer", [
        "• weights[8][E] selected by",
        "  (repReg, hiReg, orfActive)",
        "• Log-odds gradient update",
        "• Emits neuralP[0..3]",
    ], C_VIOLET, body_size=7.9, mono_footer="metaNeural Weighting")

    draw_card(ax, 55.7, 25.0, 18.8, 15.2, "Reciprocal Norm + SSE", [
        "• (num * RECIP[d]) >>> 30",
        "• Fuses bayesP + neuralP",
        "• sseCounts[1536*4] refine",
    ], C_TEAL, mono_footer="RECIP[65600] Lookup")

    draw_card(ax, 55.7, 8.6, 18.8, 13.4, "Carryless RangeCoder", [
        "• TOT = 32,768 (15-bit prec.)",
        "• TOP_VALUE = 1 << 24",
    ], C_SLATE, mono_footer="Exact Int Arithmetic")

    draw_arrow(ax, 51.3, 72.0, 55.7, 77.0, C_VIOLET)
    draw_arrow(ax, 51.3, 42.0, 55.7, 55.0, C_TEAL)
    draw_arrow(ax, 51.3, 19.0, 55.7, 49.0, C_BLUE, rad=-0.06)
    # Clean orthogonal bypass from Log-Bayes down x=73.5 right gutter into Reciprocal Norm + SSE
    draw_poly_arrow(ax, [(73.5, 67.2), (73.5, 40.2)], C_NAVY, lw=1.4)
    draw_arrow(ax, 63.9, 43.2, 63.9, 40.2, C_VIOLET)
    draw_arrow(ax, 65.1, 25.0, 65.1, 22.0, C_TEAL)

    # Column 4 (x=78.4..99.2): clean gutter x=75.4..78.4 AND right bypass lane x=95.6..98.6
    draw_zone(ax, 78.4, 7.2, 20.8, 85.6,
              "D. Stripes & Law IV Gate", C_CORAL, C_CORAL_BG,
              badge="Safety <= 2.0012 bpb", two_line_col=True)

    draw_card(ax, 79.3, 67.2, 16.2, 20.6, "WarmSnapshot Stripes", [
        "• Seg 0 warms model first",
        "• Clones immutable snapshot",
        "• Frees warmModel = null",
        "• ForkJoinPool encodes Seg",
        "  1..S-1 (4 Mbp stripes)",
    ], C_NAVY, body_size=7.9, mono_footer="TURBO_SEG = 4 Mbp")

    draw_diamond(ax, 87.4, 52.8, 7.4, 5.8,
                 "Law IV Gate:\n|C_j| >= ceil(n_j/4)?", C_CORAL, "#FFF1F2", fontsize=7.9)

    draw_card(ax, 79.3, 27.6, 15.6, 13.8, "SEG_RAW_2BIT (0)", [
        "• Random / crypto guard",
        "• 4 bases/byte raw pack",
        "• Cap: <= 2.0012 bpb",
    ], C_CORAL, body_size=7.8)

    draw_card(ax, 79.3, 8.6, 15.6, 14.4, "SEG_CRM (1)", [
        "• Compressed sub-block",
        "• Segmented byte framing",
        "  segLens[] + compLens[]",
    ], C_GREEN, body_size=7.8)

    # Orthogonal route inside dedicated gutter x=76.9 between Col C and Col D
    draw_poly_arrow(ax, [(74.5, 15.3), (76.9, 15.3), (76.9, 77.5), (79.3, 77.5)],
                    C_SLATE, label="Segs", lpt_idx=1, ly_off=0.0)
    draw_arrow(ax, 87.4, 67.2, 87.4, 58.6, C_NAVY)
    draw_arrow(ax, 87.4, 47.0, 87.4, 41.4, C_CORAL, "Yes (Bypass)")
    # Orthogonal right bypass channel outside SEG_RAW_2BIT (x=96.9)
    draw_poly_arrow(ax, [(94.8, 52.8), (97.2, 52.8), (97.2, 15.8), (94.9, 15.8)],
                    C_GREEN, label="No", lpt_idx=1, lx_off=0.0, ly_off=2.0)

    draw_footer_strip(ax, [
        ("1,122,777 B (6.75 s)", "Mode S E. coli 00.fa (vs 1,136,069 B v1)"),
        ("12 / 12 Genomes Won", "vs JARVIS3 -l 31 (-2.32%, 1.26x faster)"),
        ("108 MiB vs 1,159 MiB", "10.7x Less Compress RAM vs JARVIS3 -l 31"),
        ("87 KiB Dense Core", "ord0..7 Fits L2 Cache in Bounded Mode"),
        ("125,075 B (2.0012 bpb)", "500,000 bp SHA-256 Crypto Noise Gate"),
    ])

    save_arch(fig, "arch3_fast_novelty_coder")


# ==============================================================================
# DIAGRAM 4: .HGS2 / .HGS BINARY WIRE FORMAT & THREE-GENERATION LINEAGE
# ==============================================================================
def make_arch4():
    fig, ax = new_canvas((15.4, 10.3))
    draw_title(
        ax,
        ".hgs2 v2.0 / .hgs v1.0 Binary Wire Specification & Three-Generation Architecture Evolution",
        "Byte-Exact Bitstream Anatomy (0x48475332 & 0x48475301/HGF1), Cryptographic Panel Pinning, and Architectural Lineage from Baseline GraSS to EIMDL"
    )

    # Top Section (y=46.2..93.2): Binary Container Formats
    draw_zone(ax, 0.8, 46.2, 98.4, 47.0,
              "Part I — Physical .hgs2 v2.0 Wire Container Specification (docs/HGS-FORMAT.md)",
              C_SLATE, C_SLATE_BG, badge="Magic: 0x48475332 (HGS2)")

    # Row 1 of wire layout: compact proportional field strip (y=77.0..89.2)
    wire_blocks = [
        (2.0, 12.8, "0x00..0x03 (4 B)", [
            'magic = "HGS2"',
            "48 47 53 32",
            "Always verified",
        ], C_DARK),
        (15.8, 17.6, "0x04 (1 B) headerFlags", [
            "bits 0..1: mode (0,1,2)",
            "b2=FLAT | b3=RNA_U",
            "bits 4..7 must equal 0",
        ], C_NAVY),
        (34.4, 16.2, "0x05..0x10 (12 B)", [
            "rawByteLen (8B LE)",
            "rawCrc32 (4B LE)",
            "Full input integrity",
        ], C_TEAL),
        (51.6, 12.8, "uvarint x 2", [
            "totalRecords",
            "numUnits (U)",
            "FASTA unit cap",
        ], C_SLATE),
        (65.4, 16.4, "[mode==1] Mode R Pin", [
            "panelCount + names",
            "Deflate donor table",
            "32B SHA-256 Merkle Pin",
        ], C_BLUE),
        (82.8, 15.2, "[mode==2] Mode C DAG", [
            "Prim MST permutation",
            "rankToUnit[0..U-1]",
            "maxAncestorRank[r]<=r",
        ], C_VIOLET),
    ]
    for bx, bw, hdr_lbl, lines, col in wire_blocks:
        draw_card(ax, bx, 77.0, bw, 12.0, hdr_lbl, lines, col,
                  header_h=2.4, title_size=8.2, body_size=7.9, center_body=True)

    for i in range(len(wire_blocks) - 1):
        x_end = wire_blocks[i][0] + wire_blocks[i][1]
        x_next = wire_blocks[i + 1][0]
        draw_arrow(ax, x_end, 83.0, x_next, 83.0, C_SLATE, lw=1.2)

    # Row 2 of wire layout: Framing, Mosaic Bitstream, Novelty Superstream (y=54.6..73.6)
    draw_card(ax, 2.0, 54.6, 27.6, 19.0, "Unit Framing & 6-Ch Metadata", [
        "• If FLAG_RAW_FLAT_ACGT (0x04):",
        "   unitDnaLens[0..U-1] (0 B meta!)",
        "• Else 6-Channel Meta (per unit u):",
        "   firstRecIdx, recordCount, metaLen,",
        "   {Headers, RLE Widths/EOL, Case,",
        "    N-runs, IUPAC/T+U exceptions}",
    ], C_TEAL, mono_footer="Exact Phi^-1 Reconstruction")

    draw_card(ax, 31.4, 54.6, 33.6, 19.0, "rankMosaics[0..U-1] (EimdlMosaicCodec)", [
        "• configByte (1B): flags (Approx,Splice,V4,Wobble,Self,Tandem,Panel)",
        "• targetBases (uvarint) + noveltyBases (uvarint)",
        "• Order-1 BitRangeEnc (scale 4096, tau=20 Two-Adic Taper):",
        "   PhraseKind (LIT/PANEL/SELF/TANDEM) -> MTF Donor -> Strand (+/-)",
        "   hasSplice (8-state spanBin) -> (q, phi mod 3) -> V4 Kimura (CpG)",
    ], C_BLUE, mono_footer="Self-Contained per Rank r = 0..U-1")

    draw_card(ax, 66.8, 54.6, 31.2, 19.0, "Novelty Superstream (FastNoveltyCoder)", [
        "• totalNoveltyBases + payloadLen (uvarint)",
        "• totalSegs + modelFlags (bounded|panSuper|turbo)",
        "• segLens[0..S-1] + compLens[0..S-1]",
        "• Per Segment s: segMode (0=RAW_2BIT | 1=CRM) + body",
        "• Seg 0 warms WarmSnapshot; Seg 1..S-1 parallel",
    ], C_VIOLET, mono_footer="Strict EOF Guard: 0 Trailing Bytes")

    draw_arrow(ax, 58.0, 77.0, 15.8, 73.6, C_TEAL, lw=1.2)
    draw_arrow(ax, 29.6, 64.1, 31.4, 64.1, C_BLUE, lw=1.5)
    draw_arrow(ax, 65.0, 64.1, 66.8, 64.1, C_VIOLET, lw=1.5)

    # Legacy .hgs v1.0 comparison bar (y=47.5..53.0)
    leg = FancyBboxPatch(
        (2.0, 47.4), 96.0, 5.4,
        boxstyle="round,pad=0.0,rounding_size=0.45",
        facecolor="#F8FAFC", edgecolor="#64748B", linewidth=1.1, linestyle="--", zorder=3
    )
    ax.add_patch(leg)
    ax.text(50.0, 51.0,
            "Legacy .hgs v1.0 Container Comparison (Magic: 0x48475301 = \"HGS\\x01\", hypergrass.Main):",
            ha="center", va="center", fontsize=8.8, fontweight="bold", color=C_SLATE, zorder=4)
    ax.text(50.0, 48.7,
            "(a) Compact Header (byte 5 >= 32: auxFormat 32|33, 24-bit dnaLen)   |   (b) Multi-Block Indexed: 65,536-bp Anchor Bank + Modes 0..3 (Empty/QCRM/BC-LZG/AC-HMC) + 28-B HGF1 Footer Directory (O(1) Seek)",
            ha="center", va="center", fontsize=7.85, color=C_MUTED, zorder=4)

    # Bottom Section (y=7.0..44.6): Three-Generation Evolution
    draw_zone(ax, 0.8, 7.0, 98.4, 37.6,
              "Part II — Three-Generation Architectural Evolution & Resolved Structural Floors",
              C_NAVY, C_NAVY_BG, badge="E. coli 24g: 24.78 MB -> 3.37 MB")

    draw_card(ax, 2.0, 8.3, 22.0, 32.2, "G0: Baseline GraSS (2025)", [
        "Roy, Maity & Mukhopadhyay",
        "Briefings in Functional Genomics",
        "",
        "• Java ASCII string transforms:",
        "   {0,1,2} pairs -> {P..Y}",
        "• Monolithic tar + ./bsc -e2",
        "• Fixed vec_size caps; exits 0",
        "• 24g Mode C: 24,777,592 B",
        "  (1.6503 bpb); GRCh38 7.49 GB",
        "• 6/6 Bacterial FASTA crashes",
    ], C_CORAL, body_size=7.9, mono_footer="1.6503 bpb | Non-Pangenomic")

    draw_card(ax, 26.2, 8.3, 22.6, 32.2, "G1: HyperGraSS v1 (.hgs)", [
        "Generation 1 Multi-Block Engine",
        "hypergrass.Main / SyncmerEngine",
        "",
        "• 5-ch FASTA + 16 MiB Anchor Bank",
        "• AC-HMC + BC-LZG + QCRM",
        "• HGF1 footer + 2.08 GiB 64-bit",
        "• GraalVM native 256 MiB cgroup",
        "• E. coli Mode C: 3,810,409 B",
        "• Bottleneck: 62% serial QCRM",
        "  (887.45 s), 7.1 B/bp hash table",
    ], C_GREEN, body_size=7.9, mono_footer="0.2538 bpb | Pure JVM OOM 256M")

    draw_card(ax, 51.0, 8.3, 22.6, 32.2, "G2: RC-MDL v2 Prototype", [
        "Generation 2 Conditioning Study",
        "prototype/rcmdl/Rcmdl.java",
        "",
        "• Unaligned Mode R panel LOO",
        "• Syncmer Jaccard Prim MST",
        "• E. coli LOO: 1,636,367 B",
        "• MTB LOO: 92,506 B (-25% vs v1)",
        "• S. aureus LOO: 243,807 B",
        "• Gaps: single-gear bulldozing,",
        "  no splice (twins 13/14 lost), 1.8 GB",
    ], C_BLUE, body_size=7.9, mono_footer="Proved Mode R; Lost Mode S/C")

    draw_card(ax, 75.8, 8.3, 22.2, 32.2, "G3: Unified HGS2 / EIMDL", [
        "Generation 3 (.hgs2, Laws I–VII)",
        "HyperGraSSFusionEngine",
        "",
        "• 2-pass CSR (0.50-1.25 B/bp)",
        "• Two-Gear + Horizon + 2D Splice",
        "• V4 Kimura-CpG + Wobble + 6F-ORF",
        "• Two-Adic tau=20 + Pitman-Yor",
        "• Mode R: 1,174,113 B (E. coli)",
        "  47,339 B (MTB, 2.04x < AGC Orc)",
        "• Mode C: 3,371,489 B / 19.65s 12T",
    ], C_NAVY, body_size=7.9, mono_footer="3,452,703 B @ -Xmx128m Pure JVM")

    draw_arrow(ax, 24.0, 24.4, 26.2, 24.4, C_SLATE, lw=1.6)
    draw_arrow(ax, 48.8, 24.4, 51.0, 24.4, C_GREEN, lw=1.6)
    draw_arrow(ax, 73.6, 24.4, 75.8, 24.4, C_NAVY, lw=1.8)

    draw_footer_strip(ax, [
        ("6-Ch Phi -> 0 B on .seq", "Exact FASTA / RNA / IUPAC Identity"),
        ("32-B SHA-256 Panel Pin", "Cryptographic Mode R Integrity"),
        ("Causal DAG: anc[r] <= r", "Cycle-Free Mode C Reconstruction"),
        ("7.35x Smaller Than GraSS", "24.78 MB -> 3.37 MB on E. coli 24g"),
        ("45.2x Faster Than Gen 1", "887.45 s -> 19.65 s (12T Amdahl)"),
    ])

    save_arch(fig, "arch4_container_and_lineage")


def main():
    make_arch1()
    make_arch2()
    make_arch3()
    make_arch4()
    print("[DONE] Generated all 4 massive architectural diagrams (PDF + 300 DPI PNG).")


if __name__ == "__main__":
    main()
