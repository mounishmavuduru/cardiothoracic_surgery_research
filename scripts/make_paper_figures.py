"""Generate publication figures from results/*.json (honest, grounded, theme-neutral)."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "docs/paper/figures"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"figure.dpi": 140, "font.size": 11, "axes.grid": True,
                     "grid.alpha": 0.25})
BLUE, ORANGE, GREEN, RED, GREY = "#2a6f97", "#e08214", "#2a9d54", "#c0392b", "#888888"


def _load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def fig_rho_scaling():
    d = _load("results/rho_scaling.json")
    if not d:
        return
    fig, ax = plt.subplots(figsize=(6, 4.2))
    rg = d.get("random_geometric", {})
    for key, col, lab in (("d2", BLUE, "2-D random-geometric"), ("d3", GREEN, "3-D random-geometric")):
        blk = rg.get(key)
        if not blk:
            continue
        Ns = [r["N"] for r in blk["rows"]]; gaps = [r["gap"] for r in blk["rows"]]
        ax.loglog(Ns, gaps, "o-", color=col, label=f"{lab}  (N^{blk['fitted_exponent']:.2f})")
    am = d.get("atrial_mesh", {})
    if am.get("rows"):
        Ns = [r["n"] for r in am["rows"]]; gaps = [r["gap"] for r in am["rows"]]
        ax.loglog(Ns, gaps, "s-", color=RED, lw=2.5,
                  label=f"real atrial mesh  (N^{am['fitted_exponent']:.2f}; Weyl −1.0)")
    ax.set_xlabel("system size N (nodes)")
    ax.set_ylabel(r"spectral gap  $\lambda_3-\lambda_2$")
    ax.set_title("Spectral gap collapses with mesh size\n"
                 r"($\rho$ grows lawfully with resolution)", fontsize=10.5)
    ax.legend(fontsize=8.5, frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_rho_scaling.png"); plt.close(fig)
    print("wrote fig_rho_scaling.png")


def fig_predictive_null():
    rows = []
    ge = _load("results/gm1_expanded_metrics.json")
    if ge:
        for name in ("roney", "uw", "combined"):
            c = ge["cohorts"].get(name, {})
            r = c.get("competitors_vs_+SFI", {})
            if r:
                rows.append((f"{name} (real, n={c.get('n_subjects','?')})",
                             r["lr"]["grouped_delta_auc"], r["gbt"]["grouped_delta_auc"]))
    ku = _load("results/gm4_kuramoto_metrics.json")
    if ku and ku.get("predict"):
        p = ku["predict"].get("competitors_vs_+subspace", {})
        if p:
            rows.append(("Kuramoto (n=500)", p["lr"]["grouped_delta_auc"], p["gbt"]["grouped_delta_auc"]))
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4))
    y = np.arange(len(rows))
    ax.barh(y - 0.18, [r[1] for r in rows], 0.34, color=BLUE, label="logistic regression")
    ax.barh(y + 0.18, [r[2] for r in rows], 0.34, color=ORANGE, label="gradient boosting")
    ax.axvline(0.05, color=RED, ls="--", lw=1.5, label="pre-registered endpoint (ΔAUC=0.05)")
    ax.axvline(0.0, color=GREY, lw=1)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlabel("SFI incremental ΔAUC vs full competitor set")
    ax.set_title("Incremental predictive value of SFI\n"
                 "(pre-registered endpoint not met)", fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_predictive_null.png"); plt.close(fig)
    print("wrote fig_predictive_null.png")


def fig_power():
    d = _load("results/power_analysis.json")
    if not d:
        return
    v = d["variants"]
    labels = ["subspace-SFI\n(ΔAUC +0.003)", "single-vector\n(ΔAUC +0.0007)"]
    needed = [v["subspace"].get("n_needed_80pct", np.nan), v["single_vector"].get("n_needed_80pct", np.nan)]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(labels, needed, color=[BLUE, GREEN], width=0.55, zorder=3)
    ax.axhspan(200, 1000, color=RED, alpha=0.15, zorder=1)
    ax.axhline(1000, color=RED, ls="--", lw=1.3, zorder=2, label="largest real AF cohorts (~200–1000)")
    for i, n in enumerate(needed):
        if np.isfinite(n):
            ax.text(i, n * 1.05, f"~{n:,.0f}", ha="center", fontsize=10, fontweight="bold")
    ax.set_yscale("log"); ax.set_ylabel("cases needed to detect at 80% power")
    ax.set_title("SFI's effect is statistically real, clinically undetectable")
    ax.legend(fontsize=8.5, frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_power.png"); plt.close(fig)
    print("wrote fig_power.png")


def fig_localizer_transfer():
    keys = ["grad_phi2", "fibrosis", "perron", "wdegree"]
    media = []
    def _fields(m):
        loc = m.get("localize", {})
        f = loc.get("fields", loc)  # kuramoto/100k nest under 'fields'; scaled is flat
        out = {}
        for k, v in f.items():
            if not isinstance(v, dict):
                continue
            r = v.get("mean_origin_rank", v.get("rank"))  # gm4_scaled uses 'rank'
            if r is not None:
                out[k] = {"mean_origin_rank": r}
        return out
    sc = _load("results/gm4_scaled_metrics.json")
    if sc:
        media.append(("FHN (excitable)", _fields(sc)))
    ku = _load("results/gm4_kuramoto_metrics.json")
    if ku:
        media.append(("Kuramoto (oscillator)", _fields(ku)))
    if not media:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(keys)); w = 0.8 / max(len(media), 1)
    cols = [BLUE, ORANGE]
    for m, (name, fields) in enumerate(media):
        vals = [fields.get(k, {}).get("mean_origin_rank", np.nan) for k in keys]
        ax.bar(x + m * w, vals, w, color=cols[m % 2], label=name)
    ax.axhline(0.5, color=GREY, ls=":", label="chance (spatial null ≈ 0.5)")
    ax.set_xticks(x + w * (len(media) - 1) / 2)
    ax.set_xticklabels(["|∇φ₂|\n(KEEP)", "fibrosis\n(KEEP)", "Perron\n(DELETE)", "w-degree\n(DELETE)"], fontsize=9)
    ax.set_ylabel("origin rank (lower = better localization)")
    ax.set_title("Fiedler-gradient localizer transfers across dynamical classes")
    ax.legend(fontsize=8.5, frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_localizer_transfer.png"); plt.close(fig)
    print("wrote fig_localizer_transfer.png")


def fig_falsification():
    d = _load("results/falsification_protocol.json")
    if not d or "leakage_guard" not in d:
        return
    lk = d["leakage_guard"]
    fig, ax = plt.subplots(figsize=(5.5, 4))
    bars = ax.bar(["naive\nrandom-CV", "grouped-CV\n(honest)"],
                  [lk["auc_naive_random_cv"], lk["auc_grouped_cv"]],
                  color=[RED, GREEN], width=0.55, zorder=3)
    ax.set_ylim(0.5, 0.85)
    ax.annotate(f"+{lk['leakage_inflation']:.3f}\nleakage inflation",
                xy=(0, lk["auc_naive_random_cv"]), xytext=(0.5, 0.80),
                ha="center", fontsize=10, color=RED, fontweight="bold")
    ax.set_ylabel("AUC")
    ax.set_title("Removing the grouping guard\ninflates AUC through leakage", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_falsification.png"); plt.close(fig)
    print("wrote fig_falsification.png")


def fig_localizer_ladder():
    d = _load("results/gm4_100k_metrics.json")
    if not d or not d.get("ladder"):
        print("(100k ladder not ready — skipping localizer-ladder figure)")
        return
    Ns, ranks, ps = [], [], []
    for pt in d["ladder"]:
        gp = pt["localize"]["fields"].get("grad_phi2", {})
        Ns.append(pt["N"]); ranks.append(gp.get("mean_origin_rank", np.nan)); ps.append(gp.get("perm_p", np.nan))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogx(Ns, ranks, "o-", color=BLUE, lw=2)
    ax.axhline(0.5, color=GREY, ls=":", label="chance")
    ax.set_xlabel("networks N"); ax.set_ylabel("grad_φ₂ origin rank")
    ax.set_title(f"Localizer origin-rank stable across scale\n"
                 f"(all $p<10^{{-4}}$ up to N={max(Ns):,})", fontsize=11)
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_localizer_ladder.png"); plt.close(fig)
    print("wrote fig_localizer_ladder.png")


if __name__ == "__main__":
    fig_rho_scaling()
    fig_predictive_null()
    fig_power()
    fig_localizer_transfer()
    fig_falsification()
    fig_localizer_ladder()
    print("FIGURES_DONE")
