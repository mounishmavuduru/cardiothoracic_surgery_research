# Reproducibility guide

Everything here is deterministic given the frozen config; the only randomness is explicitly seeded.
The pre-registration (`docs/PRE_REGISTRATION.md`) was frozen before labels; the dated audit trail is
`notebooks/lab_notebook.md`.

## 1. Environment

```
python -m venv .venv && . .venv/bin/activate
pip install -e .           # numpy, scipy, scikit-learn, pandas, meshio; matplotlib for figures
```

## 2. Data (public, no PHI)

- **Roney LA virtual cohort** — Zenodo 5801337 → `data/roney/Mesh_*.vtk` (loader `asb.substrate.roney`).
- **UW/Boyle LGE-MRI cohort** — Dryad `10.5061/dryad.kkwh70sg0`. Behind an Anubis proof-of-work wall;
  `scripts/dryad_fetch.py` solves it and downloads to `data/uw_boyle/` (loader `asb.substrate.uw_boyle`).
  The 0.5 mm consolidated version is sufficient (the pipeline coarsens to ~2000-node graphs); the native
  6.88 GB version matters only for the deferred openCARP pass.

## 3. Pipeline → results

| experiment | command | output |
|---|---|---|
| GM1 predictive test (real cohorts) + expanded/UW | `python -m asb.experiments.gm1` ; `python scripts/run_uw_cohort.py` | `results/gm1_*_metrics.json` |
| GM3 validity radius | `python -m asb.experiments.gm3` | `results/gm3_metrics.json` |
| GM4 transfer (FHN) + 100k scale-up | `python scripts/run_gm4_100k.py` | `results/gm4_100k_metrics.json`, `outputs/scaled100k/` |
| GM4 third medium (Kuramoto) | `python scripts/run_gm4_kuramoto.py` | `results/gm4_kuramoto_metrics.json` |
| E6 uncertainty quantification | `python -m asb.experiments.e6_uq` | `results/e6_metrics.json` |
| Power analysis | `python scripts/power_analysis.py` | `results/power_analysis.json` |
| ρ scaling law (novelty 2) | `python scripts/rho_scaling.py` | `results/rho_scaling.json` |
| Falsification protocol (novelty 3) | `python scripts/falsification_protocol.py` | `results/falsification_protocol.json` |
| Source–sink isthmus (novelty 1) | `python scripts/run_monodomain_isthmus.py` | `results/monodomain_isthmus_metrics.json` |
| Figures | `python scripts/make_paper_figures.py` | `docs/paper/figures/*.png` |

## 4. Long-run robustness

The 100k scale-up is sharded (2500 networks/shard, atomic writes to `outputs/scaled100k/`) and resumed
by `scripts/supervise_100k.sh`, which re-launches the builder on any death and skips completed shards.
This survived a workflow-cleanup kill (at 55k) and a full container restart (at 85k) with zero loss.
Re-running any script is idempotent where a cache/output exists.

## 5. What to read, in order

1. `docs/PRE_REGISTRATION.md` — the frozen, falsifiable hypothesis and endpoints.
2. `docs/SFI_THEORY.md` — the mathematics, validity radius (§9), ρ scaling law (§10), protocol (§11).
3. `docs/paper/PREPRINT.md` — the writeup; every number traces to `results/*.json`.
4. `notebooks/lab_notebook.md` — the dated record, including the corrected GM4 audit and the honest
   negative for novelty experiment 1.
