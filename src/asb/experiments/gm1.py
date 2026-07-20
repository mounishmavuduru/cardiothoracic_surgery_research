r"""GM1 (E2) — PREDICT: does the SFI beat the real competitors head-to-head?

Runs the pre-registered primary endpoint on the real Roney cohort + monodomain-MS
labels: adding the per-region Spectral Fragility Index to the published competitor
feature set must raise grouped (patient-held-out) inducibility-classification AUC by
**ΔAUC >= 0.05 with DeLong p < 0.05**. Grouped AND naive AUC are both reported;
ΔAUC gets a >=10 000-resample bootstrap CI. A clean null (SFI == re-encoded fibrosis)
is an accepted, reportable outcome and is never tuned away.

Competitor baseline (pre-registration 7.5): fibrosis burden / spatial entropy /
patch size, deterministic min-cut, percolation threshold, and lambda2-alone. SFI is a
strict add-on/ablation on the *same* labels. Two classifiers are used: logistic
regression (with nested GroupKFold C-selection) and gradient-boosted trees.
"""
from __future__ import annotations

import json
import os
from typing import Dict, List, Sequence

import numpy as np

from asb.evaluation import delong_test
from asb.experiments.realcohort import cohort_records, records_to_frame

__all__ = ["COMPETITOR_COLS", "FIBHET_COLS", "run_gm1"]

FIBHET_COLS = (
    "fibhet_fibrosis_burden",
    "fibhet_fibrosis_spatial_entropy",
    "fibhet_fibrosis_patch_size",
)
CONN_COLS = (
    "conn_min_cut_value",
    "conn_percolation_threshold",
    "conn_lambda2_alone",
)
#: The full pre-registered competitor set the SFI must beat.
COMPETITOR_COLS = FIBHET_COLS + CONN_COLS


def _sfi_cols(columns: Sequence[str]) -> List[str]:
    return sorted(c for c in columns if c.startswith("sfi_"))


# --------------------------------------------------------------------------- #
# Classifiers.
# --------------------------------------------------------------------------- #
def _lr_pipeline(C: float, seed: int):
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    return Pipeline([
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(max_iter=2000, solver="liblinear",
                                   C=C, random_state=seed)),
    ])


def _gbt(seed: int):
    from sklearn.ensemble import GradientBoostingClassifier

    # Fixed, deliberately shallow (small N): documented as NOT tuned, to avoid
    # overfitting the ~dozens of subjects. Nesting is applied to LR only.
    return GradientBoostingClassifier(
        n_estimators=120, max_depth=2, learning_rate=0.05, subsample=0.8,
        random_state=seed,
    )


def _oof_probs_lr(X, y, groups, cols, splits, seed) -> np.ndarray:
    """Nested-GroupKFold OOF probabilities for logistic regression."""
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import roc_auc_score

    Xf = X.loc[:, list(cols)].to_numpy(dtype=float)
    y = np.asarray(y, int)
    groups = np.asarray(groups)
    oof = np.full(len(y), np.nan)
    C_grid = [0.03, 0.1, 0.3, 1.0, 3.0]
    for tr, te in splits:
        gtr = groups[tr]
        best_C, best_auc = 1.0, -1.0
        n_inner = min(3, len(np.unique(gtr)))
        if n_inner >= 2 and len(np.unique(y[tr])) == 2:
            inner = list(GroupKFold(n_splits=n_inner).split(Xf[tr], y[tr], gtr))
            for C in C_grid:
                probs = np.full(len(tr), np.nan)
                for itr, ite in inner:
                    if len(np.unique(y[tr][itr])) < 2:
                        probs[ite] = float(y[tr][itr].mean()) if len(itr) else 0.5
                        continue
                    clf = _lr_pipeline(C, seed)
                    clf.fit(Xf[tr][itr], y[tr][itr])
                    probs[ite] = clf.predict_proba(Xf[tr][ite])[:, 1]
                if len(np.unique(y[tr])) == 2:
                    auc = roc_auc_score(y[tr], np.nan_to_num(probs, nan=0.5))
                    if auc > best_auc:
                        best_auc, best_C = auc, C
        if len(np.unique(y[tr])) < 2:
            oof[te] = float(y[tr].mean()) if len(tr) else 0.5
            continue
        clf = _lr_pipeline(best_C, seed)
        clf.fit(Xf[tr], y[tr])
        oof[te] = clf.predict_proba(Xf[te])[:, 1]
    return np.where(np.isnan(oof), float(y.mean()), oof)


def _oof_probs_gbt(X, y, groups, cols, splits, seed) -> np.ndarray:
    Xf = X.loc[:, list(cols)].to_numpy(dtype=float)
    y = np.asarray(y, int)
    oof = np.full(len(y), np.nan)
    for tr, te in splits:
        if len(np.unique(y[tr])) < 2:
            oof[te] = float(y[tr].mean()) if len(tr) else 0.5
            continue
        clf = _gbt(seed)
        clf.fit(Xf[tr], y[tr])
        oof[te] = clf.predict_proba(Xf[te])[:, 1]
    return np.where(np.isnan(oof), float(y.mean()), oof)


_OOF = {"lr": _oof_probs_lr, "gbt": _oof_probs_gbt}


def _bootstrap_delta_auc(y, prob_base, prob_sfi, n_boot, seed) -> Dict[str, float]:
    """Stratified bootstrap CI for ΔAUC = AUC(base+SFI) − AUC(base)."""
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(seed)
    y = np.asarray(y, int)
    pos = np.flatnonzero(y == 1)
    neg = np.flatnonzero(y == 0)
    deltas = np.empty(n_boot)
    count = 0
    for _ in range(n_boot):
        idx = np.concatenate([
            pos[rng.integers(0, len(pos), len(pos))],
            neg[rng.integers(0, len(neg), len(neg))],
        ])
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            continue
        a = roc_auc_score(yb, prob_sfi[idx])
        b = roc_auc_score(yb, prob_base[idx])
        deltas[count] = a - b
        count += 1
    deltas = deltas[:count]
    return {
        "delta_mean": float(deltas.mean()) if count else float("nan"),
        "ci_low": float(np.percentile(deltas, 2.5)) if count else float("nan"),
        "ci_high": float(np.percentile(deltas, 97.5)) if count else float("nan"),
        "p_delta_le_0": float(np.mean(deltas <= 0.0)) if count else float("nan"),
        "n_boot_valid": int(count),
    }


def _evaluate_pair(X, y, groups, base_cols, sfi_cols, cfg) -> Dict[str, object]:
    """Grouped+naive AUC, DeLong and bootstrap for one baseline vs baseline+SFI."""
    from sklearn.model_selection import GroupKFold, KFold
    from sklearn.metrics import roc_auc_score

    y = np.asarray(y, int)
    groups = np.asarray(groups)
    n_groups = len(np.unique(groups))
    k = max(2, min(cfg["n_splits"], n_groups))
    grouped = list(GroupKFold(n_splits=k).split(X, y, groups))
    naive = list(KFold(n_splits=k, shuffle=True, random_state=cfg["seed"]).split(X, y))

    plus = list(base_cols) + list(sfi_cols)
    out: Dict[str, object] = {}
    for clf_name, oof_fn in _OOF.items():
        gb = oof_fn(X, y, groups, base_cols, grouped, cfg["seed"])
        gs = oof_fn(X, y, groups, plus, grouped, cfg["seed"])
        nb = oof_fn(X, y, groups, base_cols, naive, cfg["seed"])
        ns = oof_fn(X, y, groups, plus, naive, cfg["seed"])
        grouped_auc_base = float(roc_auc_score(y, gb))
        grouped_auc_sfi = float(roc_auc_score(y, gs))
        naive_auc_base = float(roc_auc_score(y, nb))
        naive_auc_sfi = float(roc_auc_score(y, ns))
        dl = delong_test(y, gs, gb)  # paired on grouped OOF; delta = SFI - base
        boot = _bootstrap_delta_auc(y, gb, gs, cfg["n_boot"], cfg["seed"])
        out[clf_name] = {
            "grouped_auc_base": grouped_auc_base,
            "grouped_auc_sfi": grouped_auc_sfi,
            "grouped_delta_auc": grouped_auc_sfi - grouped_auc_base,
            "naive_auc_base": naive_auc_base,
            "naive_auc_sfi": naive_auc_sfi,
            "naive_delta_auc": naive_auc_sfi - naive_auc_base,
            "delong_p": dl["p"], "delong_z": dl["z"],
            "bootstrap": boot,
            "endpoint_met": bool((grouped_auc_sfi - grouped_auc_base) >= 0.05
                                 and dl["p"] < 0.05),
        }
    return out


def run_gm1(
    roney_dir: str = "data/roney",
    *, limit: int | None = None, n_jobs: int = 4, n_boot: int = 10000,
    n_splits: int = 5, seed: int = 0, outputs_dir: str = "outputs",
) -> dict:
    """Run the GM1 head-to-head and write ``gm1_metrics.json`` + a report.

    Returns the metrics dict. The primary comparison is the full competitor set
    (``COMPETITOR_COLS``) vs competitor+SFI; a secondary compares the
    fibrosis-heterogeneity subset vs +SFI.
    """
    records = cohort_records(roney_dir, limit=limit, n_jobs=n_jobs,
                             cache_dir=outputs_dir)
    X, y, groups = records_to_frame(records)
    sfi_cols = _sfi_cols(X.columns)

    n_ind = int(y.sum())
    cfg = {"n_splits": n_splits, "seed": seed, "n_boot": n_boot}

    comparisons = {}
    # Require >= 2 in each class so the DeLong variance (which divides by n-1 per
    # class) and the stratified bootstrap are both well-defined.
    if min(n_ind, len(y) - n_ind) >= 2:
        comparisons["competitors_vs_+SFI"] = _evaluate_pair(
            X, y, groups, COMPETITOR_COLS, sfi_cols, cfg)
        comparisons["fibrosis_vs_+SFI"] = _evaluate_pair(
            X, y, groups, FIBHET_COLS, sfi_cols, cfg)

    # Spectral-gap distribution across the cohort (pre-registration 4).
    gaps = [r.features.get("spec_gap", float("nan")) for r in records]

    metrics = {
        "cohort": {
            "n_subjects": int(len(y)),
            "n_inducible": n_ind,
            "inducible_fraction": float(n_ind / len(y)) if len(y) else 0.0,
            "n_groups": int(len(set(groups))),
            "n_features_total": int(X.shape[1]),
            "n_sfi_features": len(sfi_cols),
            "both_classes": bool(0 < n_ind < len(y)),
        },
        "sfi_columns": sfi_cols,
        "competitor_columns": list(COMPETITOR_COLS),
        "comparisons": comparisons,
        "spectral_gap": {
            "min": float(np.nanmin(gaps)), "median": float(np.nanmedian(gaps)),
            "max": float(np.nanmax(gaps)),
        },
        "config": {
            "n_boot": n_boot, "n_splits": n_splits, "seed": seed,
            "label_source": "monodomain_ms",
        },
        "label_disclaimer": (
            "Labels are monodomain Mitchell-Schaeffer simulator verdicts on real "
            "Roney LA anatomy. NOT clinical POAF; NOT openCARP (which is deferred)."
        ),
    }

    os.makedirs(outputs_dir, exist_ok=True)
    with open(os.path.join(outputs_dir, "gm1_metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2)
    _write_report(os.path.join(outputs_dir, "gm1_report.md"), metrics)
    return metrics


def _fmt(x) -> str:
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "n/a" if not np.isfinite(xf) else f"{xf:.4f}"


def _write_report(path: str, m: dict) -> None:
    c = m["cohort"]
    lines = ["# GM1 — SFI vs competitors on real Roney anatomy (monodomain-MS labels)",
             "",
             "> Labels are a monodomain Mitchell-Schaeffer **simulator verdict**, "
             "not clinical POAF and not openCARP (deferred).", "",
             f"- Subjects: **{c['n_subjects']}** (patient-held-out GroupKFold), "
             f"inducible **{c['n_inducible']}** ({_fmt(c['inducible_fraction'])})",
             f"- SFI features added: **{c['n_sfi_features']}**; "
             f"competitors: {len(m['competitor_columns'])}", ""]
    for name, comp in m["comparisons"].items():
        lines += [f"## {name}", "",
                  "| clf | grouped base | grouped +SFI | grouped ΔAUC | DeLong p "
                  "| naive base | naive +SFI | boot ΔAUC [95% CI] | endpoint |",
                  "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for clf, r in comp.items():
            b = r["bootstrap"]
            lines.append(
                f"| {clf} | {_fmt(r['grouped_auc_base'])} | {_fmt(r['grouped_auc_sfi'])} "
                f"| {_fmt(r['grouped_delta_auc'])} | {_fmt(r['delong_p'])} "
                f"| {_fmt(r['naive_auc_base'])} | {_fmt(r['naive_auc_sfi'])} "
                f"| {_fmt(b['delta_mean'])} [{_fmt(b['ci_low'])}, {_fmt(b['ci_high'])}] "
                f"| {'MET' if r['endpoint_met'] else 'not met'} |")
        lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
