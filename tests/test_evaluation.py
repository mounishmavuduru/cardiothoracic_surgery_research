"""Tests for asb.evaluation -- leakage-controlled evaluation and statistics.

Self-contained: depends only on numpy, pandas, sklearn, asb.config and the
module under test. Synthetic X/y are constructed locally.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from asb.config import EvalConfig
from asb.evaluation import (
    calibration_curve_data,
    colocalization_vs_null,
    delong_test,
    nested_group_kfold_auc,
    roc_auc_ci,
)


# --------------------------------------------------------------------------- #
# DeLong
# --------------------------------------------------------------------------- #
def test_delong_dominance_gives_positive_delta_and_valid_p():
    """When A strictly dominates B, delta > 0 and p is a valid probability."""
    rng = np.random.default_rng(0)
    n = 200
    y = np.concatenate([np.ones(n // 2), np.zeros(n // 2)]).astype(int)

    # A: perfect separation. B: noisy (weak) separation.
    prob_a = np.where(y == 1, rng.uniform(0.6, 1.0, n), rng.uniform(0.0, 0.4, n))
    prob_b = y * 0.15 + rng.uniform(0.0, 1.0, n)

    res = delong_test(y, prob_a, prob_b)
    assert res["auc_a"] > res["auc_b"]
    assert res["delta"] > 0.0
    assert 0.0 <= res["p"] <= 1.0
    assert np.isfinite(res["z"])


def test_delong_identical_inputs_delta_zero_p_one():
    """Identical predictors give delta ~ 0 and p = 1."""
    rng = np.random.default_rng(1)
    y = np.concatenate([np.ones(50), np.zeros(50)]).astype(int)
    prob = rng.uniform(0, 1, 100)

    res = delong_test(y, prob, prob.copy())
    assert res["delta"] == pytest.approx(0.0, abs=1e-9)
    assert res["auc_a"] == pytest.approx(res["auc_b"], abs=1e-12)
    assert res["p"] == pytest.approx(1.0)
    assert 0.0 <= res["p"] <= 1.0


def test_delong_requires_both_classes():
    with pytest.raises(ValueError):
        delong_test(np.ones(10, dtype=int), np.linspace(0, 1, 10), np.linspace(0, 1, 10))


# --------------------------------------------------------------------------- #
# Grouped vs naive AUC on a deliberately leaky dataset
# --------------------------------------------------------------------------- #
def _leaky_dataset(seed=0):
    """Exact-duplicate rows across a random cross-subject label => leakage.

    Each subject has a random 'fingerprint' feature vector and a random label.
    Duplicating each subject's row (tiny noise) makes the label memorizable when
    copies leak into the training fold (naive KFold) but not when a subject is
    held out as a group (GroupKFold).
    """
    rng = np.random.default_rng(seed)
    n_subjects = 24
    d = 15
    reps = 4

    fingerprints = rng.normal(size=(n_subjects, d))
    labels = np.array([0, 1] * (n_subjects // 2))
    rng.shuffle(labels)

    rows, ys, groups = [], [], []
    for s in range(n_subjects):
        for _ in range(reps):
            rows.append(fingerprints[s] + rng.normal(scale=1e-3, size=d))
            ys.append(labels[s])
            groups.append(f"family_{s}")

    X = pd.DataFrame(np.asarray(rows), columns=[f"f{i}" for i in range(d)])
    y = np.asarray(ys, dtype=int)
    return X, y, groups


def test_grouped_auc_le_naive_auc_on_leaky_data():
    X, y, groups = _leaky_dataset(seed=3)
    cfg = EvalConfig(n_splits=5)
    feature_sets = {"all": list(X.columns)}

    res = nested_group_kfold_auc(X, y, groups, feature_sets, cfg, seed=7)
    grouped = res["all"]["grouped_auc"]
    naive = res["all"]["naive_auc"]

    assert 0.0 <= grouped <= 1.0
    assert 0.0 <= naive <= 1.0
    # Leakage inflates the naive estimate; grouped must not exceed it.
    assert grouped <= naive + 1e-9
    # And the gap should be substantial for this constructed leak.
    assert naive > grouped


def test_nested_group_kfold_multiple_feature_sets():
    X, y, groups = _leaky_dataset(seed=1)
    cfg = EvalConfig(n_splits=4)
    feature_sets = {"half": list(X.columns[:7]), "all": list(X.columns)}
    res = nested_group_kfold_auc(X, y, groups, feature_sets, cfg, seed=0)
    assert set(res.keys()) == {"half", "all"}
    for v in res.values():
        assert {"grouped_auc", "naive_auc"} <= set(v.keys())


# --------------------------------------------------------------------------- #
# Bootstrap CI
# --------------------------------------------------------------------------- #
def test_roc_auc_ci_bounds_and_determinism():
    rng = np.random.default_rng(0)
    y = np.concatenate([np.ones(60), np.zeros(60)]).astype(int)
    prob = np.where(y == 1, rng.uniform(0.4, 1.0, 120), rng.uniform(0.0, 0.6, 120))

    r1 = roc_auc_ci(y, prob, n_boot=200, seed=42)
    r2 = roc_auc_ci(y, prob, n_boot=200, seed=42)
    assert r1 == r2  # deterministic in seed
    assert 0.0 <= r1["ci_low"] <= r1["auc"] <= r1["ci_high"] <= 1.0
    assert r1["se"] >= 0.0


# --------------------------------------------------------------------------- #
# Colocalization vs null
# --------------------------------------------------------------------------- #
def test_colocalization_returns_percentile_in_unit_interval():
    rng = np.random.default_rng(0)
    n = 300
    uac = rng.uniform(0, 1, size=(n, 2))
    origin_node = 17
    # Hotspot scores peak near the true origin -> should colocalize.
    d = np.linalg.norm(uac - uac[origin_node], axis=1)
    scores = np.exp(-(d ** 2) / 0.01)

    res = colocalization_vs_null(scores, origin_node, uac, n_null=200, seed=5)
    assert 0.0 <= res["p"] <= 1.0
    assert res["observed"] >= 0.0
    assert res["null_p5"] >= 0.0
    assert isinstance(res["passed"], bool)


def test_colocalization_deterministic():
    rng = np.random.default_rng(2)
    uac = rng.uniform(0, 1, size=(120, 2))
    scores = rng.uniform(0, 1, 120)
    r1 = colocalization_vs_null(scores, 3, uac, n_null=100, seed=9)
    r2 = colocalization_vs_null(scores, 3, uac, n_null=100, seed=9)
    assert r1 == r2


# --------------------------------------------------------------------------- #
# Calibration
# --------------------------------------------------------------------------- #
def test_calibration_curve_data_shapes_and_ranges():
    rng = np.random.default_rng(0)
    prob = rng.uniform(0, 1, 500)
    y = (rng.uniform(0, 1, 500) < prob).astype(int)  # well-calibrated by design

    res = calibration_curve_data(y, prob, n_bins=10)
    assert len(res["prob_pred"]) == len(res["prob_true"]) == len(res["counts"])
    assert res["counts"].sum() == 500
    assert np.all(res["prob_pred"] >= 0.0) and np.all(res["prob_pred"] <= 1.0)
    assert np.all(res["prob_true"] >= 0.0) and np.all(res["prob_true"] <= 1.0)
    assert 0.0 <= res["brier"] <= 1.0
    assert 0.0 <= res["ece"] <= 1.0
    # A well-calibrated predictor has small ECE.
    assert res["ece"] < 0.15
