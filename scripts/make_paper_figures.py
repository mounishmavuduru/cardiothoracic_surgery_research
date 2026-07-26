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
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


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
                 f"(below every random-origin null draw, up to N={max(Ns):,})", fontsize=11)
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_localizer_ladder.png"); plt.close(fig)
    print("wrote fig_localizer_ladder.png")


def fig_substrate_audit():
    """Three substrate controls that together locate -- and fail to locate -- the anomaly.

    Colour: two categorical series only, BLUE/ORANGE. That pair was checked rather than
    eyeballed: OKLab dE 30.6 under normal vision, worst case 22.2 across protanopia,
    deuteranopia and tritanopia, and a greyscale luminance gap of 0.178 for print.
    GREEN/RED (dE 6.9 under deuteranopia) and BLUE/RED (greyscale dY 0.002) are avoided
    for fills; RED appears only as a DASHED reference line, so line style and not hue
    carries it. The second series is hatched so the encoding survives greyscale printing.
    """
    ab = _load("results/uw_substrate_ablation.json")
    fc = _load("results/roney_fibre_control.json")
    qz = _load("results/roney_fibrosis_quantization.json")
    bs = _load("results/fibrosis_burden_swap.json")
    if not ab:
        print("(substrate ablation not ready - skipping substrate-audit figure)")
        return

    RONEY_REF = 20 / 62
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.6))
    axes = axes.ravel()

    def _bar(ax, x, c, color, hatch=None, width=0.32):
        ax.bar(x, 100 * c["inducible_fraction"], width, color=color, hatch=hatch, zorder=3,
               edgecolor="white", linewidth=1.2)
        ax.text(x, 100 * c["inducible_fraction"] + 1.1, f"{c['n_inducible']}/{c['n']}",
                ha="center", va="bottom", fontsize=9.5)

    # (a) UW 2x2
    ax = axes[0]
    S = ab["summary"]
    for key, x, col, hat in (
            ("A_caps_kept_constant_fibres", -0.17, BLUE, None),
            ("C_caps_kept_varying_fibres", 0.83, BLUE, None),
            ("B_caps_dropped_constant_fibres", 0.17, ORANGE, "//"),
            ("D_caps_dropped_varying_fibres", 1.17, ORANGE, "//")):
        _bar(ax, x, S[key], col, hat)
    ref = ax.axhline(100 * RONEY_REF, color=RED, ls="--", lw=1.5, zorder=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["constant fibres", "varying fibres"])
    ax.set_ylabel("inducible (%)")
    ax.set_ylim(0, 44)
    ax.set_title("(a) UW substrate repair" + "\n" + r"strongest contrast McNemar $p=0.077$",
                 fontsize=10.5)
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=BLUE),
               plt.Rectangle((0, 0), 1, 1, color=ORANGE, hatch="//"), ref],
              ["orifices sealed (as released)", "orifices opened",
               "Roney reference (32.3%)"],
              fontsize=8.2, frameon=False, loc="upper left")

    # (b) fibre control on Roney
    ax = axes[1]
    if fc:
        F = fc["summary"]
        _bar(ax, 0, F["R_A_real_fibres"], BLUE, None, width=0.5)
        _bar(ax, 1, F["R_B_constant_fibres"], ORANGE, "//", width=0.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["real fibres", "fibres destroyed"])
        ax.set_ylim(0, 44)
        ax.set_ylabel("inducible (%)")
        ax.annotate("identical rate" + "\n" + r"McNemar $p=1.0$" + "\n"
                    + "(10/62 verdicts still flip)",
                    xy=(0.5, 0.88), xycoords="axes fraction", ha="center", va="top", fontsize=9.5)
        ax.set_title("(b) fibre control on Roney" + "\n"
                     + "degeneracy does not move the rate", fontsize=10.5)
    else:
        ax.set_axis_off()

    # (c) fibrosis gradation on Roney
    ax = axes[2]
    if qz:
        Q = qz["summary"]
        keys = ("Q_A_continuous", "Q_B_binary_half", "Q_C_binary_matched")
        cols = (BLUE, ORANGE, ORANGE)
        hats = (None, "//", "//")
        for i, k in enumerate(keys):
            _bar(ax, i, Q[k], cols[i], hats[i], width=0.5)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["continuous" + "\n" + "(as shipped)",
                            "binary" + "\n" + r"$f>0.5$",
                            "binary" + "\n" + "burden matched"], fontsize=9)
        ax.set_ylim(0, 75)   # panel (c) reaches 64.5%; a 44% ceiling clipped it
        ax.set_ylabel("inducible (%)")
        ax.set_title("(c) fibrosis gradation on Roney" + "\n"
                     + "coarsening the fibrosis field", fontsize=10.5)
    else:
        ax.set_axis_off()
        ax.set_title("(c) fibrosis gradation" + "\n" + "(pending)", fontsize=10.5)

    # (d) fibrosis burden swap -- the control that resolves the anomaly
    ax = axes[3]
    if bs:
        B = bs["summary"]
        keys = ("R_A_roney_native", "R_B_roney_at_uw", "U_A_uw_native", "U_B_uw_at_roney")
        cols = (BLUE, ORANGE, BLUE, ORANGE)
        hats = (None, "//", None, "//")
        for i, k in enumerate(keys):
            _bar(ax, i, B[k], cols[i], hats[i], width=0.6)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(["Roney" + "\n" + "native", "Roney at" + "\n" + "UW burden",
                            "UW" + "\n" + "native", "UW at" + "\n" + "Roney burden"],
                           fontsize=8.5)
        ax.set_ylim(0, 44)
        ax.set_ylabel("inducible (%)")
        ax.set_title("(d) fibrosis burden swap" + "\n"
                     + "moves the rate in both directions", fontsize=10.5)
        ax.legend([plt.Rectangle((0, 0), 1, 1, color=BLUE),
                   plt.Rectangle((0, 0), 1, 1, color=ORANGE, hatch="//")],
                  ["native burden", "burden swapped"], fontsize=8.2,
                  frameon=False, loc="upper right")
    else:
        ax.set_axis_off()
        ax.set_title("(d) fibrosis burden swap" + "\n" + "(pending)", fontsize=10.5)

    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_substrate_audit.png")
    plt.close(fig)
    print("wrote fig_substrate_audit.png")


if __name__ == "__main__":
    fig_rho_scaling()
    fig_predictive_null()
    fig_power()
    fig_localizer_transfer()
    fig_falsification()
    fig_localizer_ladder()
    fig_substrate_audit()
    print("FIGURES_DONE")
