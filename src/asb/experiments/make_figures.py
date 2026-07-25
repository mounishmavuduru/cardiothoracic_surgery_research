r"""Publication figures for AtrialSpectralBench, rendered from cached result JSON.

One-command figure generation for the paper/poster. Each function guards on the
availability of its inputs (``results/*_metrics.json``) and is a pure renderer.
The two centerpieces are :func:`fig_validity_radius` (GM3 — the honest core) and
:func:`fig_cross_medium` (GM4 — the transfer). Uses the non-interactive Agg backend.
"""
from __future__ import annotations

import json
import os
from typing import Dict, Optional

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

__all__ = ["make_all_figures"]

_RESULTS = "results"
_FIGDIR = "results/figures"


def _load(name: str) -> Optional[dict]:
    path = os.path.join(_RESULTS, name)
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


# --------------------------------------------------------------------------- #
# Fig 5 (centerpiece) — the validity radius (GM3).
# --------------------------------------------------------------------------- #
def fig_validity_radius(out: str) -> Optional[str]:
    m = _load("gm3_metrics.json")
    if not m:
        return None
    sw = m["validity_radius_sweep"]["sweep"]
    rho = np.array([r["rho"] for r in sw])
    first = np.array([r["rel_err_first"] for r in sw])
    sub = np.array([r["rel_err_subspace"] for r in sw])
    rho_star = m["validity_radius_sweep"].get("rho_star_10pct")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(rho, 100 * first, "o-", color="#c0392b", label="first-order SFI error")
    ax.plot(rho, 100 * np.clip(sub, 0, 5), "s--", color="#2980b9", label="subspace SFI error")
    ax.axhline(10, color="gray", ls=":", lw=1)
    if rho_star:
        ax.axvline(rho_star, color="#27ae60", lw=1.5)
        ax.text(rho_star * 1.1, 60, f"ρ* ≈ {rho_star:.1f}\n(10% error)", color="#27ae60")
    # Real-cohort operating points.
    for label, rv, col in [("atria (ρ≈2422)", 2422, "#8e44ad"),
                           ("neural net (ρ≈1852)", 1852, "#d35400")]:
        ax.axvline(rv, color=col, ls="-", lw=1, alpha=0.6)
        ax.text(rv, 3, label, rotation=90, va="bottom", ha="right", color=col, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel(r"perturbation ratio  ρ = ‖ΔL‖ / (λ₃ − λ₂)")
    ax.set_ylabel("relative error of the Δλ₂ prediction (%)")
    ax.set_title("GM3 — validity radius of the linear Spectral Fragility Index")
    ax.legend(loc="center left")
    ax.set_ylim(-3, 100)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Fig 6 (centerpiece) — cross-medium keep/delete transfer (GM2 + GM4).
# --------------------------------------------------------------------------- #
def fig_cross_medium(out: str) -> Optional[str]:
    g2 = _load("gm2_metrics.json")
    g4 = _load("gm4_metrics.json")
    if not (g2 and g4):
        return None
    order = ["grad_phi2", "perron", "combined", "wdegree", "fibrosis", "fibrosis_grad"]
    labels = ["|∇φ₂|", "Perron", "|∇φ₂|∩Perron", "w-degree", "substrate", "∇substrate"]
    atr = [g2["variants"].get(k, {}).get("mean_origin_score_rank", np.nan) for k in order]
    neu = [g4["localize"]["fields"].get(k, {}).get("mean_origin_rank", np.nan) for k in order]

    x = np.arange(len(order)); w = 0.38
    fig, ax = plt.subplots(figsize=(8, 5))
    b1 = ax.bar(x - w / 2, atr, w, label="atrium (cardiac reentry)", color="#8e44ad")
    b2 = ax.bar(x + w / 2, neu, w, label="neural net (FHN seizure focus)", color="#d35400")
    ax.axhline(0.5, color="gray", ls="--", lw=1, label="chance (spatial null)")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("mean origin score-rank  (lower = localizes better)")
    ax.set_title("GM4 — the Fiedler-gradient localizer transfers across media\n"
                 "(|∇φ₂| localizes the origin in both; centrality/degree do not)")
    ax.legend()
    # KEEP annotations.
    for bars, vals in [(b1, atr), (b2, neu)]:
        for rect, v in zip(bars, vals):
            if np.isfinite(v) and v < 0.4:
                ax.text(rect.get_x() + rect.get_width() / 2, v + 0.01, "KEEP",
                        ha="center", fontsize=7, color="#27ae60", fontweight="bold")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Fig 4 — GM1 ΔAUC (the honest null).
# --------------------------------------------------------------------------- #
def fig_gm1_auc(out: str) -> Optional[str]:
    m = _load("gm1_metrics.json")
    if not m:
        return None
    comps = m["comparisons"]
    names, base, plus = [], [], []
    for cname, comp in comps.items():
        for clf, r in comp.items():
            names.append(f"{cname.replace('_vs_+SFI','')}\n{clf}")
            base.append(r["grouped_auc_base"]); plus.append(r["grouped_auc_sfi"])
    x = np.arange(len(names)); w = 0.38
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - w / 2, base, w, label="competitors", color="#7f8c8d")
    ax.bar(x + w / 2, plus, w, label="competitors + SFI", color="#2980b9")
    ax.axhline(0.5, color="gray", ls=":", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("grouped (patient-held-out) AUC")
    ax.set_ylim(0.5, 1.0)
    ax.set_title("GM1 — adding SFI does not beat the competitor set (pre-registered null)")
    ax.legend()
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Fig 7 — E6 sensitivity (Morris μ* + Δw robustness).
# --------------------------------------------------------------------------- #
def fig_e6_sensitivity(out: str) -> Optional[str]:
    m = _load("e6_metrics.json")
    if not m:
        return None
    ee = m["morris_label_screen"]["elementary_effects"]
    names = list(ee.keys())
    mu = [ee[n]["mu_star"] for n in names]
    sig = [ee[n]["sigma"] for n in names]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5))
    a1.barh(names, mu, color="#16a085")
    a1.set_xlabel("Morris μ* (label sensitivity)")
    a1.set_title("Ground-truth label: which physics knob matters")
    dw = m["delta_w_sensitivity"]
    ks = list(dw.keys())
    lr = [dw[k]["lr"]["grouped_delta_auc"] for k in ks]
    gbt = [dw[k]["gbt"]["grouped_delta_auc"] for k in ks]
    x = np.arange(len(ks)); w = 0.38
    a2.bar(x - w / 2, lr, w, label="LR", color="#7f8c8d")
    a2.bar(x + w / 2, gbt, w, label="GBT", color="#2980b9")
    a2.axhline(0.05, color="#c0392b", ls="--", lw=1, label="endpoint 0.05")
    a2.axhline(0.0, color="gray", lw=0.8)
    a2.set_xticks(x); a2.set_xticklabels([f"Δw={k}" for k in ks])
    a2.set_ylabel("GM1 grouped ΔAUC"); a2.set_title("Conclusion robust to Δw"); a2.legend()
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


def make_all_figures(figdir: str = _FIGDIR) -> Dict[str, Optional[str]]:
    """Render every available figure; skip those whose inputs are missing."""
    os.makedirs(figdir, exist_ok=True)
    figs = {
        "fig5_validity_radius": fig_validity_radius(os.path.join(figdir, "fig5_validity_radius.png")),
        "fig6_cross_medium": fig_cross_medium(os.path.join(figdir, "fig6_cross_medium.png")),
        "fig4_gm1_auc": fig_gm1_auc(os.path.join(figdir, "fig4_gm1_auc.png")),
        "fig7_e6_sensitivity": fig_e6_sensitivity(os.path.join(figdir, "fig7_e6_sensitivity.png")),
    }
    return figs
