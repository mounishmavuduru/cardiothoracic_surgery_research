"""Power analysis: how many cases are needed to DETECT the SFI predictive edge?

The large-N synthetic cohort pins the SFI predictive contribution to a stable,
tiny effect size (subspace-SFI ΔAUC ~= +0.003; single-vector ~= 0). This script
turns "tiny" into a hard number: using the converged effect and its DeLong
standard error, it computes the minimum sample size to detect that ΔAUC at 80%
and 90% power (two-sided, alpha=0.05), then contrasts it with realistic clinical
atrial-fibrillation cohort sizes. Reuses the project's exact nested-GroupKFold LR
so the effect size matches the reported GM1 numbers.
"""
import json

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from asb.evaluation import delong_test
from asb.experiments.gm1 import COMPETITOR_COLS, _oof_probs_lr
from asb.experiments.gm4_scale import _load

# z for two-sided alpha=0.05 at target power: z_{1-a/2} + z_{power}
Z80 = 1.959964 + 0.841621
Z90 = 1.959964 + 1.281552
# realistic reference cohort sizes (real AF-ablation LGE-MRI studies)
CLINICAL_REF = {"typical_single_center": 200, "large_multicenter": 1000}

recs = _load("outputs/scaled100k")
y = np.array([int(r["inducible"]) for r in recs], int)
groups = [f"snet_{r['seed']}" for r in recs]
rows = [r["features"] for r in recs]
X = pd.DataFrame(rows).reindex(sorted({k for row in rows for k in row}), axis=1).fillna(0.0)
seven = ("total", "region_max", "region_mean", "region_std", "top1", "top2", "top3")
variants = {"single_vector": [f"sfi_{s}" for s in seven if f"sfi_{s}" in X.columns],
            "subspace": sorted(c for c in X.columns if c.startswith("sfisub_"))}

n = len(y)
npos = int(y.sum())
print(f"N={n} positives={npos} ({npos / n:.1%})", flush=True)
splits = list(GroupKFold(n_splits=5).split(X, y, groups))
base = _oof_probs_lr(X, y, groups, COMPETITOR_COLS, splits, 0)

out = {"N": n, "positives": npos, "clinical_reference_sizes": CLINICAL_REF, "variants": {}}
for name, cols in variants.items():
    sfi = _oof_probs_lr(X, y, groups, list(COMPETITOR_COLS) + cols, splits, 0)
    d = delong_test(y, sfi, base)
    delta, z = d["delta"], d["z"]
    rec = {"delta_auc": delta, "z_at_N": z, "p_at_N": d["p"],
           "auc_base": d["auc_b"], "auc_sfi": d["auc_a"]}
    if abs(z) > 1e-9:
        rec["n_needed_80pct"] = float(n * (Z80 / abs(z)) ** 2)
        rec["n_needed_90pct"] = float(n * (Z90 / abs(z)) ** 2)
        rec["detectable_at_200"] = bool(rec["n_needed_80pct"] <= 200)
        rec["detectable_at_1000"] = bool(rec["n_needed_80pct"] <= 1000)
    else:
        rec["n_needed_80pct"] = float("inf")
        rec["n_needed_90pct"] = float("inf")
    out["variants"][name] = rec
    print(f"[{name}] dAUC={delta:+.5f} z@N={z:.2f} p={d['p']:.2e} "
          f"n_needed(80%)={rec['n_needed_80pct']:.0f} n_needed(90%)={rec['n_needed_90pct']:.0f}",
          flush=True)

json.dump(out, open("results/power_analysis.json", "w"), indent=2)
print("POWER_DONE", flush=True)
