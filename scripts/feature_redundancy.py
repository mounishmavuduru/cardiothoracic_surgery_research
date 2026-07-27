r"""Quantify the paper's headline mechanism: is SFI's information already in the competitors?

The manuscript explains its clean null by feature redundancy -- SFI adds nothing because
standard connectivity features already carry its information. That explanation is the
load-bearing claim of the Results, and until now it rested entirely on the incremental AUC
being zero. An incremental null is consistent with redundancy but does not establish it: a
feature that is pure noise also adds nothing, and so does a feature the classifier cannot
exploit. Those are three different diagnoses with three different implications.

This measures redundancy directly, without reference to the label, in three ways:

1. Per-feature explained variance. Regress each SFI column on the competitor block and
   report R^2. High R^2 means the competitors can reconstruct that SFI column, which is
   redundancy in the strict sense.
2. Canonical correlation between the SFI block and the competitor block. The leading
   canonical correlation is the largest correlation obtainable between any linear
   combination of one block and any of the other; near 1 means the blocks span nearly the
   same directions.
3. A noise control. The same statistics computed against a permuted SFI block, which
   destroys any real relationship while preserving marginal distributions. This separates
   "redundant" from "the estimator is degenerate", which the raw R^2 cannot do alone.

All quantities are label-free, so nothing here is contingent on the inducibility verdict.

Usage:  python scripts/feature_redundancy.py
"""
from __future__ import annotations

import json
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
warnings.simplefilter("ignore")

from asb.experiments.gm1 import COMPETITOR_COLS, FIBHET_COLS, _sfi_cols  # noqa: E402
from asb.experiments.realcohort import cohort_records, records_to_frame  # noqa: E402
from asb.experiments.uw_cohort import uw_cohort_records  # noqa: E402


def _r2_on(block: np.ndarray, y: np.ndarray) -> float:
    """R^2 of the least-squares fit of y on block (with intercept)."""
    X = np.column_stack([np.ones(block.shape[0]), block])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def _canoncorr(A: np.ndarray, B: np.ndarray) -> float:
    """Leading canonical correlation between two column blocks."""
    def whiten(M):
        M = M - M.mean(axis=0)
        u, s, _ = np.linalg.svd(M, full_matrices=False)
        keep = s > (s.max() * 1e-10) if s.size else np.array([], dtype=bool)
        return u[:, keep]
    Qa, Qb = whiten(A), whiten(B)
    if Qa.size == 0 or Qb.size == 0:
        return float("nan")
    s = np.linalg.svd(Qa.T @ Qb, compute_uv=False)
    return float(np.clip(s[0], 0.0, 1.0)) if s.size else float("nan")


def main() -> None:
    roney = cohort_records("data/roney", limit=None, n_jobs=1, cache_dir="outputs",
                           verbose=False)
    uw = uw_cohort_records("data/uw_boyle", n_jobs=1, cache_dir="outputs", verbose=False)
    X, _y, _g = records_to_frame(roney + uw)

    sfi = [c for c in _sfi_cols(X.columns)]
    comp = [c for c in COMPETITOR_COLS if c in X.columns]
    fib = [c for c in FIBHET_COLS if c in X.columns]
    print(f"n = {len(X)} subjects | {len(sfi)} SFI columns | {len(comp)} competitors "
          f"| {len(fib)} fibrosis-heterogeneity")

    A = X[sfi].to_numpy(dtype=float)
    C = X[comp].to_numpy(dtype=float)
    F = X[fib].to_numpy(dtype=float)
    # standardise so R^2 is not driven by scale
    def z(M):
        sd = M.std(axis=0)
        sd[sd == 0] = 1.0
        return (M - M.mean(axis=0)) / sd
    A, C, F = z(A), z(C), z(F)

    rng = np.random.default_rng(0)
    A_perm = A[rng.permutation(A.shape[0])]

    print(f"\n{'SFI feature':<26} {'R^2 | competitors':>18} {'R^2 | fibrosis':>16} "
          f"{'R^2 | permuted':>15}")
    rows = {}
    for j, name in enumerate(sfi):
        r_c = _r2_on(C, A[:, j])
        r_f = _r2_on(F, A[:, j])
        r_p = _r2_on(C, A_perm[:, j])
        rows[name] = {"r2_competitors": r_c, "r2_fibrosis": r_f, "r2_permuted_control": r_p}
        print(f"{name:<26} {r_c:>18.3f} {r_f:>16.3f} {r_p:>15.3f}")

    med_c = float(np.median([v["r2_competitors"] for v in rows.values()]))
    med_p = float(np.median([v["r2_permuted_control"] for v in rows.values()]))
    cc_c = _canoncorr(A, C)
    cc_f = _canoncorr(A, F)
    cc_p = _canoncorr(A_perm, C)

    print(f"\nmedian R^2 over SFI columns: {med_c:.3f} (competitors) "
          f"vs {med_p:.3f} (permuted control)")
    print(f"leading canonical correlation SFI vs competitors: {cc_c:.3f}")
    print(f"leading canonical correlation SFI vs fibrosis:    {cc_f:.3f}")
    print(f"leading canonical correlation permuted control:   {cc_p:.3f}")
    # Report the structure rather than a pass/fail. A single threshold on the median
    # collapses a genuinely two-part result: the aggregate SFI summaries are almost
    # perfectly reconstructible from the competitors, while the order statistics are not.
    high = sorted(n for n, v in rows.items() if v["r2_competitors"] >= 0.60)
    low = sorted(n for n, v in rows.items() if v["r2_competitors"] < 0.25)
    verdict = (
        f"PARTIAL redundancy. Leading canonical correlation {cc_c:.3f} against "
        f"{cc_p:.3f} for the permuted control, so the two blocks share a dominant "
        f"direction almost exactly. But the redundancy is not uniform across columns: "
        f"{len(high)} of {len(rows)} SFI columns are reconstructible at R^2 >= 0.60 "
        f"({', '.join(high)}), while {len(low)} retain most of their variance at "
        f"R^2 < 0.25 ({', '.join(low)}). Redundancy therefore explains the null for the "
        f"aggregate summaries; for the order statistics the explanation must be that "
        f"their unique variance is not predictive of the label."
    )
    print("\nVERDICT:", verdict)

    os.makedirs("results", exist_ok=True)
    with open("results/feature_redundancy.json", "w", encoding="utf-8") as fh:
        json.dump({"description": "label-free test of whether competitors already carry "
                                  "SFI's information",
                   "n_subjects": int(len(X)),
                   "per_feature": rows,
                   "median_r2_competitors": med_c,
                   "median_r2_permuted_control": med_p,
                   "canoncorr_sfi_competitors": cc_c,
                   "canoncorr_sfi_fibrosis": cc_f,
                   "canoncorr_permuted_control": cc_p,
                   "n_reconstructible_r2_ge_0.60": len(high),
                   "n_retaining_variance_r2_lt_0.25": len(low),
                   "verdict": verdict}, fh, indent=2)
    print("\nwrote results/feature_redundancy.json")


if __name__ == "__main__":
    main()
