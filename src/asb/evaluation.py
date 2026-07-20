"""Leakage-controlled evaluation and statistics for AtrialSpectralBench.

This module holds the honest-broker statistics used to decide whether the
Spectral Fragility Index (SFI) adds real predictive value over the fibrosis /
connectivity baselines. Its guiding concern is *leakage*: because the synthetic
cohort is built from a small number of base anatomies (shape families) that are
then statistically resampled into near-duplicate variants, a naive random
cross-validation split will place near-copies of the same atrium on both sides
of the split and report an optimistically inflated AUC. :func:`nested_group_kfold_auc`
exposes that gap by reporting a group-aware AUC (folds split on ``shape_family``)
alongside the naive AUC.

The remaining routines are the inferential gates:

* :func:`delong_test` -- paired DeLong test comparing two ROC curves on the same
  labels (fast midrank algorithm), used to test SFI vs baseline.
* :func:`roc_auc_ci` -- bootstrap confidence interval on a single AUC.
* :func:`colocalization_vs_null` -- does the predicted SFI hotspot land near the
  true reentry origin more than a UAC spatial-null (rotational / shift) would?
* :func:`calibration_curve_data` -- reliability-diagram data for a probabilistic
  classifier.

All stochastic routines take an explicit ``seed`` and build a local
``np.random.Generator``; there is no global random state and no wall-clock
seeding.
"""
from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from asb.config import EvalConfig

__all__ = [
    "nested_group_kfold_auc",
    "delong_test",
    "roc_auc_ci",
    "colocalization_vs_null",
    "calibration_curve_data",
]


# --------------------------------------------------------------------------- #
# Fast DeLong machinery (Sun & Xu 2014; DeLong et al. 1988)
# --------------------------------------------------------------------------- #
def _compute_midrank(x: np.ndarray) -> np.ndarray:
    """Midranks of ``x`` (ties share the average rank).

    Parameters
    ----------
    x : (n,) ndarray
        Values to rank.

    Returns
    -------
    (n,) ndarray
        Midrank of each element, in the original order.
    """
    x = np.asarray(x, dtype=float)
    order = np.argsort(x)
    xs = x[order]
    n = len(x)
    tr = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n and xs[j] == xs[i]:
            j += 1
        tr[i:j] = 0.5 * (i + j - 1) + 1.0  # 1-based average rank over the tie block
        i = j
    out = np.empty(n, dtype=float)
    out[order] = tr
    return out


def _fast_delong(predictions_sorted: np.ndarray, m: int):
    """Fast DeLong AUC + covariance for a set of predictors on shared labels.

    Parameters
    ----------
    predictions_sorted : (k, n) ndarray
        ``k`` predictors' scores, with the first ``m`` columns the positive
        (label == 1) samples and the remaining columns the negatives.
    m : int
        Number of positive samples.

    Returns
    -------
    aucs : (k,) ndarray
        AUC of each predictor.
    cov : (k, k) ndarray
        Covariance matrix of the AUC estimates.
    """
    n_total = predictions_sorted.shape[1]
    n = n_total - m
    k = predictions_sorted.shape[0]

    tx = np.empty([k, m], dtype=float)
    ty = np.empty([k, n], dtype=float)
    tz = np.empty([k, n_total], dtype=float)
    for r in range(k):
        tx[r, :] = _compute_midrank(predictions_sorted[r, :m])
        ty[r, :] = _compute_midrank(predictions_sorted[r, m:])
        tz[r, :] = _compute_midrank(predictions_sorted[r, :])

    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    delongcov = np.atleast_2d(delongcov)
    return aucs, delongcov


def _auc_and_variance(y_true: np.ndarray, scores: np.ndarray):
    """AUC and its DeLong variance for a single predictor."""
    order = (-y_true).argsort(kind="mergesort")
    label_1_count = int(y_true.sum())
    sorted_scores = scores[order][np.newaxis, :]
    aucs, cov = _fast_delong(sorted_scores, label_1_count)
    return float(aucs[0]), float(cov[0, 0])


def delong_test(y_true, prob_a, prob_b) -> Dict[str, float]:
    """Paired DeLong test comparing two ROC curves on the same labels.

    Tests the null hypothesis that the two predictors ``prob_a`` and ``prob_b``
    have equal AUC, accounting for the correlation induced by scoring the *same*
    samples. Uses the fast midrank algorithm of Sun & Xu (2014).

    Parameters
    ----------
    y_true : (n,) array_like
        Binary ground-truth labels (0 / 1). Must contain both classes.
    prob_a, prob_b : (n,) array_like
        Predicted scores / probabilities of the positive class from predictors
        A and B, aligned with ``y_true``.

    Returns
    -------
    dict
        ``{"auc_a", "auc_b", "delta", "z", "p"}`` where ``delta = auc_a - auc_b``,
        ``z`` is the DeLong test statistic and ``p`` is the two-sided p-value in
        ``[0, 1]``. When the two score vectors are identical (or their variances
        cancel exactly) ``delta`` is ~0, ``z`` is 0 and ``p`` is 1.
    """
    y_true = np.asarray(y_true).astype(float).ravel()
    prob_a = np.asarray(prob_a, dtype=float).ravel()
    prob_b = np.asarray(prob_b, dtype=float).ravel()

    classes = np.unique(y_true)
    if not np.array_equal(classes, np.array([0.0, 1.0])):
        raise ValueError("y_true must be binary with both classes present (0 and 1).")

    order = (-y_true).argsort(kind="mergesort")
    label_1_count = int(y_true.sum())
    predictions_sorted = np.vstack((prob_a, prob_b))[:, order]
    aucs, cov = _fast_delong(predictions_sorted, label_1_count)

    auc_a, auc_b = float(aucs[0]), float(aucs[1])
    delta = auc_a - auc_b
    var = cov[0, 0] + cov[1, 1] - 2.0 * cov[0, 1]

    if var <= 0.0:
        # Identical (or perfectly correlated equal-variance) predictors.
        z = 0.0
        p = 1.0
    else:
        z = delta / np.sqrt(var)
        # Two-sided normal p-value.
        from scipy.stats import norm

        p = float(2.0 * norm.sf(abs(z)))
        p = min(max(p, 0.0), 1.0)

    return {"auc_a": auc_a, "auc_b": auc_b, "delta": float(delta),
            "z": float(z), "p": float(p)}


# --------------------------------------------------------------------------- #
# Grouped vs naive cross-validated AUC
# --------------------------------------------------------------------------- #
def _make_classifier(seed: int) -> Pipeline:
    """Simple, standardized logistic-regression classifier (deterministic)."""
    return Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    max_iter=1000, solver="liblinear", random_state=seed
                ),
            ),
        ]
    )


def _oof_auc(X, y, splits, seed) -> float:
    """Out-of-fold AUC given a list of (train_idx, test_idx) splits."""
    from sklearn.metrics import roc_auc_score

    y = np.asarray(y).astype(int).ravel()
    oof = np.full(len(y), np.nan, dtype=float)
    for train_idx, test_idx in splits:
        if len(np.unique(y[train_idx])) < 2:
            # Degenerate training fold; fall back to prior probability.
            oof[test_idx] = float(y[train_idx].mean()) if len(train_idx) else 0.5
            continue
        clf = _make_classifier(seed)
        clf.fit(X[train_idx], y[train_idx])
        oof[test_idx] = clf.predict_proba(X[test_idx])[:, 1]
    # Guard against any unfilled entries.
    if np.isnan(oof).any():
        oof = np.where(np.isnan(oof), float(y.mean()), oof)
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, oof))


def nested_group_kfold_auc(
    X, y, groups, feature_sets: Dict[str, Sequence], cfg: EvalConfig, seed: int
) -> Dict[str, Dict[str, float]]:
    """Grouped vs naive cross-validated AUC for each named feature set.

    For every feature set the design matrix is cross-validated twice: once with
    :class:`sklearn.model_selection.GroupKFold` (folds split on ``groups`` so no
    shape family straddles the train/test boundary) and once with an ordinary
    shuffled :class:`~sklearn.model_selection.KFold` that ignores groups. The gap
    between the two AUCs is the leakage the group-aware protocol removes: because
    resampled near-duplicate atria share a ``shape_family``, the naive AUC is
    (weakly) optimistically biased.

    Parameters
    ----------
    X : pandas.DataFrame or (n, d) array_like
        Design matrix. If a DataFrame, ``feature_sets`` values are column names;
        otherwise they are integer column indices.
    y : (n,) array_like
        Binary inducibility labels (from ``mock_ep`` in development -- never
        clinical POAF).
    groups : (n,) array_like
        Group id per subject (typically ``shape_family``) for GroupKFold.
    feature_sets : dict
        Maps a feature-set name to the list of columns it uses, e.g.
        ``{"fibrosis": [...], "fibrosis+SFI": [...]}``.
    cfg : EvalConfig
        Provides ``n_splits``.
    seed : int
        Seed for fold shuffling and the classifier.

    Returns
    -------
    dict
        ``{name: {"grouped_auc": float, "naive_auc": float}}`` for each feature
        set. AUCs are out-of-fold ROC-AUC over the whole cohort.
    """
    y = np.asarray(y).astype(int).ravel()
    groups = np.asarray(groups).ravel()
    n = len(y)

    # Resolve column selection into a numpy matrix per feature set.
    def _select(cols) -> np.ndarray:
        if hasattr(X, "loc") and hasattr(X, "columns"):
            return np.asarray(X.loc[:, list(cols)].to_numpy(), dtype=float)
        return np.asarray(X, dtype=float)[:, list(cols)]

    n_groups = len(np.unique(groups))
    n_splits = int(min(cfg.n_splits, n_groups, n))
    n_splits = max(n_splits, 2)

    group_kf = GroupKFold(n_splits=min(n_splits, n_groups))
    naive_kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)

    results: Dict[str, Dict[str, float]] = {}
    for name, cols in feature_sets.items():
        Xf = _select(cols)

        grouped_splits = list(group_kf.split(Xf, y, groups))
        naive_splits = list(naive_kf.split(Xf, y))

        grouped_auc = _oof_auc(Xf, y, grouped_splits, seed)
        naive_auc = _oof_auc(Xf, y, naive_splits, seed)

        results[name] = {"grouped_auc": grouped_auc, "naive_auc": naive_auc}

    return results


# --------------------------------------------------------------------------- #
# Bootstrap CI on a single AUC
# --------------------------------------------------------------------------- #
def roc_auc_ci(y_true, prob, n_boot: int, seed: int) -> Dict[str, float]:
    """Bootstrap confidence interval for a single ROC-AUC.

    Parameters
    ----------
    y_true : (n,) array_like
        Binary labels (0 / 1).
    prob : (n,) array_like
        Predicted positive-class scores.
    n_boot : int
        Number of stratified bootstrap resamples.
    seed : int
        Seed for the bootstrap generator.

    Returns
    -------
    dict
        ``{"auc", "ci_low", "ci_high", "se"}`` -- the point AUC, the 2.5th and
        97.5th percentiles of the bootstrap distribution, and its standard error.
    """
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true).astype(int).ravel()
    prob = np.asarray(prob, dtype=float).ravel()

    auc = float(roc_auc_score(y_true, prob))
    pos = np.flatnonzero(y_true == 1)
    neg = np.flatnonzero(y_true == 0)

    boots = np.empty(n_boot, dtype=float)
    count = 0
    for _ in range(n_boot):
        # Stratified resample keeps both classes present in every replicate.
        pi = rng.integers(0, len(pos), size=len(pos))
        ni = rng.integers(0, len(neg), size=len(neg))
        idx = np.concatenate([pos[pi], neg[ni]])
        yb = y_true[idx]
        pb = prob[idx]
        if len(np.unique(yb)) < 2:
            continue
        boots[count] = roc_auc_score(yb, pb)
        count += 1

    if count == 0:
        return {"auc": auc, "ci_low": auc, "ci_high": auc, "se": 0.0}

    boots = boots[:count]
    ci_low = float(np.percentile(boots, 2.5))
    ci_high = float(np.percentile(boots, 97.5))
    se = float(np.std(boots, ddof=1)) if count > 1 else 0.0
    return {"auc": auc, "ci_low": ci_low, "ci_high": ci_high, "se": se}


# --------------------------------------------------------------------------- #
# Colocalization of predicted hotspot vs true reentry origin
# --------------------------------------------------------------------------- #
def colocalization_vs_null(
    hotspot_scores, origin_node: int, uac, n_null: int, seed: int
) -> Dict[str, float]:
    """Test whether the predicted hotspot colocalizes with the reentry origin.

    The predicted hotspot is the node with the maximal ``hotspot_scores`` value.
    The observed statistic is its UAC (universal atrial coordinate) distance to
    the true ``origin_node``. To decide whether this is closer than chance, we
    build a spatial null by applying random rotations / shifts to the UAC field
    (treating alpha, beta as periodic on ``[0, 1)``) and, for each, recomputing
    the distance from the (transformed) hotspot location to the fixed origin
    coordinate. A small observed distance relative to the null lower tail means
    the hotspot genuinely tracks the origin rather than reflecting UAC geometry.

    Parameters
    ----------
    hotspot_scores : (n,) array_like
        Per-node SFI hotspot score; larger = more likely reentry initiation.
    origin_node : int
        Node index of the true reentry origin.
    uac : (n, 2) array_like
        Universal atrial coordinates (alpha, beta) in ``[0, 1]``.
    n_null : int
        Number of rotational / shift spatial-null resamples.
    seed : int
        Seed for the null generator.

    Returns
    -------
    dict
        ``{"observed", "null_p5", "p", "passed"}`` -- the observed hotspot-to-
        origin UAC distance, the 5th percentile of the null distance
        distribution, the fraction of null distances ``<=`` observed (a
        percentile / p-value in ``[0, 1]``), and a boolean pass flag that is True
        when ``observed <= null_p5`` (significant colocalization at ~5%).
    """
    rng = np.random.default_rng(seed)
    scores = np.asarray(hotspot_scores, dtype=float).ravel()
    uac = np.asarray(uac, dtype=float)
    origin = np.asarray(uac[int(origin_node)], dtype=float)

    def _torus_dist(a: np.ndarray, b: np.ndarray) -> float:
        d = np.abs(a - b)
        d = np.minimum(d, 1.0 - d)  # periodic (torus) distance on [0, 1)
        return float(np.sqrt(np.sum(d ** 2)))

    hotspot_idx = int(np.argmax(scores))
    observed = _torus_dist(uac[hotspot_idx], origin)

    null_dists = np.empty(n_null, dtype=float)
    for t in range(n_null):
        shift = rng.uniform(0.0, 1.0, size=2)
        # Optional swap of the two UAC axes for a coarse "rotation".
        if rng.random() < 0.5:
            shifted = (uac[:, ::-1] + shift) % 1.0
        else:
            shifted = (uac + shift) % 1.0
        idx = int(np.argmax(scores))  # scores are fixed; geometry is what moves
        null_dists[t] = _torus_dist(shifted[idx], origin)

    null_p5 = float(np.percentile(null_dists, 5.0))
    p = float(np.mean(null_dists <= observed))
    passed = bool(observed <= null_p5)

    return {
        "observed": float(observed),
        "null_p5": null_p5,
        "p": p,
        "passed": passed,
    }


# --------------------------------------------------------------------------- #
# Calibration / reliability-diagram data
# --------------------------------------------------------------------------- #
def calibration_curve_data(y_true, prob, n_bins: int = 10) -> Dict[str, np.ndarray]:
    """Reliability-diagram data for a probabilistic classifier.

    Bins predicted probabilities into ``n_bins`` equal-width bins on ``[0, 1]``
    and, for each non-empty bin, reports the mean predicted probability, the
    observed positive fraction, and the bin count. Also returns the Brier score
    and expected calibration error (ECE).

    Parameters
    ----------
    y_true : (n,) array_like
        Binary labels (0 / 1).
    prob : (n,) array_like
        Predicted positive-class probabilities in ``[0, 1]``.
    n_bins : int, optional
        Number of equal-width probability bins (default 10).

    Returns
    -------
    dict
        ``{"prob_pred", "prob_true", "counts", "bin_edges", "brier", "ece"}``.
        ``prob_pred`` / ``prob_true`` / ``counts`` cover only non-empty bins.
    """
    y_true = np.asarray(y_true).astype(float).ravel()
    prob = np.asarray(prob, dtype=float).ravel()

    edges = np.linspace(0.0, 1.0, n_bins + 1)
    # Bin index in [0, n_bins-1]; clip so prob == 1.0 lands in the last bin.
    binids = np.clip(np.digitize(prob, edges[1:-1], right=False), 0, n_bins - 1)

    prob_pred: List[float] = []
    prob_true: List[float] = []
    counts: List[int] = []
    ece = 0.0
    n = len(y_true)
    for b in range(n_bins):
        mask = binids == b
        c = int(mask.sum())
        if c == 0:
            continue
        mp = float(prob[mask].mean())
        mt = float(y_true[mask].mean())
        prob_pred.append(mp)
        prob_true.append(mt)
        counts.append(c)
        ece += (c / n) * abs(mt - mp)

    brier = float(np.mean((prob - y_true) ** 2)) if n else float("nan")

    return {
        "prob_pred": np.asarray(prob_pred, dtype=float),
        "prob_true": np.asarray(prob_true, dtype=float),
        "counts": np.asarray(counts, dtype=int),
        "bin_edges": edges,
        "brier": brier,
        "ece": float(ece),
    }
