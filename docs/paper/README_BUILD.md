# Building the manuscript

The paper is a single, self-contained LaTeX source: **`manuscript.tex`**.
It has no external `.bib` (the bibliography is inline via `thebibliography`) and
pulls its figures from `figures/*.png`, so it compiles anywhere with a standard
TeX install — nothing else in the repository is needed to build the PDF.

## Option A — Overleaf (recommended for editing/sharing)

Overleaf compiles committed **LaTeX (`.tex`) files**; it does **not** generate a
manuscript from the project's Python/JS code. The workflow is:

1. In Overleaf: **New Project -> Import from GitHub** and select this repository
   (or upload a zip of the `docs/paper/` folder).
2. Overleaf auto-detects the main document (the file containing
   `\documentclass`, i.e. `docs/paper/manuscript.tex`). If it doesn't, set it as
   the main file in *Menu -> Main document*.
3. Compile. The figures resolve because `manuscript.tex` sets
   `\graphicspath{{figures/}}` and the PNGs live in `docs/paper/figures/`.

Everything the paper needs is already in the repo, so the import is one step and
then it just compiles. (If you import the whole repo rather than just
`docs/paper/`, point Overleaf's main document at `docs/paper/manuscript.tex`.)

## Option B — local

```
# one-off dependencies (Debian/Ubuntu):
sudo apt-get install -y --no-install-recommends \
  texlive-latex-base texlive-latex-recommended texlive-latex-extra \
  texlive-fonts-recommended texlive-science lmodern latexmk

# then, from the repo root:
make paper
# or directly:
cd docs/paper && latexmk -pdf manuscript.tex
```

The compiled `docs/paper/manuscript.pdf` (14 pages) is committed for convenience;
`make paper` regenerates it. Intermediate files (`.aux`, `.log`, `.out`, `.toc`,
`.fls`, `.fdb_latexmk`) are git-ignored.

## Figures

The figures are produced from `results/*.json` by
`scripts/make_paper_figures.py` (matplotlib, deterministic — no randomness). To
regenerate them:

```
.venv/bin/python scripts/make_paper_figures.py     # writes docs/paper/figures/*.png
```

The convergence figure is produced separately by
`scripts/convergence_analyze.py roney`.

## What is where

| section of the paper            | grounded in                                    |
|---------------------------------|------------------------------------------------|
| §3.1 predictive null            | `results/gm1_expanded_metrics.json`, `gm4_*`   |
| §4.2 validity radius            | `results/gm3_metrics.json`, `docs/SFI_THEORY.md` |
| §4.3 rho scaling law            | `results/rho_scaling.json`                     |
| §4.4 localization transfer      | `results/gm4_*_metrics.json`                   |
| §4.5 falsification protocol     | `results/falsification_protocol.json`, `power_analysis.json` |
| §4.6 uncertainty quantification | `results/e6_metrics.json`                      |
| §4.8 resolution convergence     | `results/convergence_summary_roney.json`       |
| Appendix A (Coq theorems)       | `formal/sfi_edge_identity_Q.v`                 |

Every number in the manuscript traces to one of these files; the narrative
mirror is `PREPRINT.md`.
