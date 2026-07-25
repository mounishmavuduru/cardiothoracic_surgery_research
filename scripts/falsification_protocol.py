"""Novelty experiment 3: the falsification protocol (naive vs rigorous).

The biomarker literature is full of spectral markers that 'work' under naive
methodology. This assembles the contrast: for each guard in our protocol, REMOVING
it manufactures a false positive our protocol catches.

  (a) No grouping  -> shape-family LEAKAGE inflates accuracy (controlled demo here).
  (b) No effect-size gate -> the P-VALUE TRAP: at large N a +0.003 AUC 'improvement'
      is p<1e-4 and looks real (assembled from results/power_analysis.json).

Together: a naive analyst would have reported SFI as a working biomarker; the
pre-registered, grouped, effect-size-gated protocol correctly rejects it.
"""
import json
import os

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from asb.config import Config
from asb.experiments.realcohort import FROZEN_SFI
from asb.features import subject_features
from asb.transfer.network import NetworkConfig, make_excitable_network


def _lr():
    return Pipeline([("s", StandardScaler()),
                     ("c", LogisticRegression(max_iter=2000, solver="liblinear"))])


def _oof(X, y, splits):
    oof = np.full(len(y), np.nan)
    for tr, te in splits:
        if len(np.unique(y[tr])) < 2:
            oof[te] = y[tr].mean() if len(tr) else 0.5
            continue
        m = _lr().fit(X[tr], y[tr]); oof[te] = m.predict_proba(X[te])[:, 1]
    return np.nan_to_num(oof, nan=0.5)


def leakage_demo(n_base=12, n_rep=6, n_nodes=300, seed=0):
    """Replicated shape-family cohort: naive random-CV vs grouped-CV AUC."""
    rng = np.random.default_rng(seed)
    cfg = NetworkConfig(n_nodes=n_nodes, radius=0.11, n_region_grid=3)
    rows, groups, labels = [], [], []
    for b in range(n_base):
        base_burden = float(rng.uniform(0.08, 0.45))
        # label is a noisy function of the BASE (shared by all replicas) -> the leak
        y_base = int(base_burden + 0.15 * rng.standard_normal() > 0.28)
        for r in range(n_rep):
            G = make_excitable_network(seed * 1000 + b * 50 + r, cfg,
                                       lesion_burden=base_burden * float(1 + 0.05 * rng.standard_normal()))
            feats = subject_features(G, Config(seed=0, sfi=FROZEN_SFI), np.random.default_rng(b * 50 + r))
            rows.append(feats); groups.append(b); labels.append(y_base)
    X = pd.DataFrame(rows).fillna(0.0)
    cols = [c for c in X.columns if c.startswith(("fibhet_", "conn_", "spec_"))]
    Xc = X[cols].to_numpy(float)
    y = np.array(labels, int); groups = np.array(groups)
    naive = list(KFold(n_splits=6, shuffle=True, random_state=0).split(Xc))
    grouped = list(GroupKFold(n_splits=min(6, n_base)).split(Xc, y, groups))
    auc_naive = roc_auc_score(y, _oof(Xc, y, naive)) if len(np.unique(y)) == 2 else float("nan")
    auc_group = roc_auc_score(y, _oof(Xc, y, grouped)) if len(np.unique(y)) == 2 else float("nan")
    return {"n_base": n_base, "n_rep": n_rep, "n_rows": len(y), "positives": int(y.sum()),
            "auc_naive_random_cv": float(auc_naive), "auc_grouped_cv": float(auc_group),
            "leakage_inflation": float(auc_naive - auc_group)}


def main():
    out = {}
    print("[leakage] building replicated shape-family cohort ...", flush=True)
    out["leakage_guard"] = leakage_demo()
    lk = out["leakage_guard"]
    print(f"[leakage] naive random-CV AUC={lk['auc_naive_random_cv']:.3f} vs "
          f"grouped-CV AUC={lk['auc_grouped_cv']:.3f}  (inflation "
          f"{lk['leakage_inflation']:+.3f})", flush=True)

    # p-value trap, assembled from the power analysis
    if os.path.exists("results/power_analysis.json"):
        pa = json.load(open("results/power_analysis.json"))
        ss = pa["variants"]["subspace"]
        out["effect_size_guard"] = {
            "naive_conclusion": f"subspace-SFI improves prediction, p={ss['p_at_N']:.1e} at N={pa['N']} -> 'it works'",
            "rigorous_conclusion": f"effect size dAUC={ss['delta_auc']:+.4f}; needs ~{ss['n_needed_80pct']:.0f} cases for 80% power; clinical cohorts ~200-1000 -> undetectable, useless",
            "delta_auc": ss["delta_auc"], "p_at_N": ss["p_at_N"],
            "n_needed_80pct": ss["n_needed_80pct"]}
        print(f"[p-value trap] naive: p={ss['p_at_N']:.1e} 'works'  |  rigorous: dAUC="
              f"{ss['delta_auc']:+.4f}, needs ~{ss['n_needed_80pct']:.0f} cases -> useless", flush=True)

    out["summary"] = ("A naive analyst (random CV + p-value threshold, no effect-size gate) would "
                      "report SFI as a working biomarker. The pre-registered protocol "
                      "(shape-family GroupKFold + DeLong + effect-size gate + spatial nulls) "
                      "correctly rejects it. Each guard, removed, manufactures a specific false positive.")
    json.dump(out, open("results/falsification_protocol.json", "w", encoding="utf-8"), indent=2)
    print("FALSIFICATION_DONE", flush=True)


if __name__ == "__main__":
    main()
