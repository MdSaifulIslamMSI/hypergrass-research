# Research Monograph Build & Automated Numerical Audit (`paper/`)

This directory contains the reproducible publication pipeline for the **33-page single-column LaTeX research monograph** on the four-generation **`HyperGraSS` / `RC-MDL` / `HGS2` / `HGS3` (`EIMDL`)** genomic compression architecture (Laws I–VII, OpenAI `2026` ETH transfer + §1.1.3 OVC scope boundary, Next-Gen EIMDL v3 `.hgs3`/Native/compressed-domain genomics, and Kaggle E2E cloud verification).

## Contents

| File | Description |
|---|---|
| [`EIMDL_Main_Manuscript.pdf`](EIMDL_Main_Manuscript.pdf) | Compiled 7-page flagship bioRxiv / journal manuscript (built via `build_main_manuscript.py`, `0` `Overfull \hbox` warnings) |
| [`EIMDL_Lean4_Paper.pdf`](EIMDL_Lean4_Paper.pdf) | Compiled 4-page standalone Lean 4 formal-methods companion pearl (built via `build_lean4_paper.py`, `0` `Overfull \hbox` warnings) |
| [`HGS2_EIMDL_Research_Paper.pdf`](HGS2_EIMDL_Research_Paper.pdf) | Compiled 33-page single-column PDF monograph (`0` `Overfull \hbox` warnings, `24` tables, `14` figures) |
| [`HGS2_EIMDL_Research_Paper.tex`](HGS2_EIMDL_Research_Paper.tex) | Generated LaTeX source file (built deterministically from committed JSON evidence) |
| [`references.bib`](references.bib) | BibTeX bibliography covering genomic compression, syncmers, Kimura/codon models, MDL, `roy2025grass`, and `openai2026eth` |
| [`generate_figures.py`](generate_figures.py) | Generates all 9 publication-grade 300 DPI data figures (`figures/fig1_regime_overview.png` .. `fig9_nextgen_eimdl_v3_paradigm.png`) directly from committed JSON results |
| [`generate_architecture_diagrams.py`](generate_architecture_diagrams.py) | Generates all 5 massive 300 DPI + vector PDF system/wire architectural blueprints (`arch1_eimdl_system_pipeline` .. `arch5_eimdl_v3_paradigm`) |
| [`table_builders.py`](table_builders.py) | Builds the per-genome Mode R LOO, Cluster-LOO, Mode C, Mode S, and independent benchmark LaTeX tables directly from JSON |
| [`build_main_manuscript.py`](build_main_manuscript.py) | Generates & compiles `EIMDL_Main_Manuscript.tex` / `.pdf` |
| [`build_lean4_paper.py`](build_lean4_paper.py) | Generates & compiles `EIMDL_Lean4_Paper.tex` / `.pdf` |
| [`build_latex_paper.py`](build_latex_paper.py) | Generates `HGS2_EIMDL_Research_Paper.tex` |
| [`build_markdown_paper.py`](build_markdown_paper.py) | Generates the full Markdown monograph [`../docs/PAPER-DRAFT.md`](../docs/PAPER-DRAFT.md) |
| [`audit_paper_numbers.py`](audit_paper_numbers.py) | Executes automated assertions verifying every byte count, base count, bit-rate, PDF existence, Lean theorem (`0 sorry`), and all 14 figure/blueprint embeds across `.tex`, `docs/PAPER-DRAFT.md`, `README.md`, and `HANDOFF.md` |
| [`figures/`](figures/) | `fig1_regime_overview` .. `fig9_nextgen_eimdl_v3_paradigm` and `arch1_eimdl_system_pipeline` .. `arch5_eimdl_v3_paradigm` (PNG + PDF vector outputs) |

## Rebuilding Everything End-to-End

From the repository root (`C:\Users\mdsai\.gemini\antigravity\scratch\hypergrass`):

```powershell
# 1. Regenerate all 9 data figures + 5 architectural blueprints (PNG + PDF)
python paper/generate_figures.py
python paper/generate_architecture_diagrams.py

# 2. Rebuild the LaTeX monograph (paper/HGS2_EIMDL_Research_Paper.tex) and Markdown monograph (docs/PAPER-DRAFT.md)
python paper/build_latex_paper.py
python paper/build_markdown_paper.py

# 3. Compile the 33-page PDF via pdflatex + bibtex
Push-Location paper
pdflatex -interaction=nonstopmode HGS2_EIMDL_Research_Paper.tex
bibtex HGS2_EIMDL_Research_Paper
pdflatex -interaction=nonstopmode HGS2_EIMDL_Research_Paper.tex
pdflatex -interaction=nonstopmode HGS2_EIMDL_Research_Paper.tex
Pop-Location

# 4. Run the automated numerical audit (verifies 100% of numbers against JSON benchmark reports)
python paper/audit_paper_numbers.py
```
