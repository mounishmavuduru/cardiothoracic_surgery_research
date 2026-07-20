"""Tests for asb.features -- design-matrix assembly.

Builds a tiny synthetic cohort locally and checks the namespaced feature vectors,
the stacked design matrix shape/labels/groups, and determinism.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from asb.config import Config
from asb.features import (
    build_design_matrix,
    sfi_feature_vector,
    spectral_feature_vector,
    subject_features,
)
from asb.labels.mock_ep import induce
from asb.substrate.mesh import mesh_to_graph
from asb.substrate.synthetic import make_cohort


def _tiny_cfg() -> Config:
    cfg = Config()
    cfg.cohort.n_base = 2
    cfg.cohort.n_variants = 2
    cfg.cohort.n_nodes = 42
    return cfg


def _tiny_cohort(cfg: Config):
    meshes = make_cohort(cfg.cohort)
    graphs = [mesh_to_graph(m) for m in meshes]
    labels = [
        induce(G, cfg.labels, np.random.default_rng(cfg.seed + 1 + i))
        for i, G in enumerate(graphs)
    ]
    return graphs, labels


def test_spectral_feature_vector_has_expected_keys_and_finite():
    cfg = _tiny_cfg()
    graphs, _ = _tiny_cohort(cfg)
    feats = spectral_feature_vector(graphs[0])
    for key in ("gap", "spectral_entropy", "cheeger_estimate", "lambda2",
                "perron_ipr", "fiedler_ipr", "fiedler_abs_mean"):
        assert key in feats
    assert all(np.isfinite(v) for v in feats.values())


def test_sfi_feature_vector_keys_nonneg_and_deterministic():
    cfg = _tiny_cfg()
    graphs, _ = _tiny_cohort(cfg)
    G = graphs[0]

    v1 = sfi_feature_vector(G, cfg.sfi, np.random.default_rng(0))
    v2 = sfi_feature_vector(G, cfg.sfi, np.random.default_rng(123))

    for key in ("total", "region_max", "region_mean", "top1", "top2", "top3",
                "edge_frag_max"):
        assert key in v1
    # Analytic expectation is deterministic in G/cfg regardless of rng.
    assert v1 == v2
    # SFI is a non-negative expected drop in connectivity.
    assert v1["total"] >= 0.0
    assert v1["region_max"] >= v1["region_mean"] >= 0.0
    # top-k are sorted descending.
    assert v1["top1"] >= v1["top2"] >= v1["top3"]


def test_subject_features_namespacing():
    cfg = _tiny_cfg()
    graphs, _ = _tiny_cohort(cfg)
    feats = subject_features(graphs[0], cfg, np.random.default_rng(0))
    prefixes = {"fibhet_", "conn_", "spec_", "sfi_"}
    for key in feats:
        assert any(key.startswith(p) for p in prefixes), key
    # Every feature block is represented.
    for p in prefixes:
        assert any(k.startswith(p) for k in feats), p


def test_build_design_matrix_shapes_and_groups():
    cfg = _tiny_cfg()
    graphs, labels = _tiny_cohort(cfg)
    X, y, groups = build_design_matrix(graphs, labels, cfg, np.random.default_rng(0))

    assert isinstance(X, pd.DataFrame)
    assert X.shape[0] == len(graphs)
    assert len(y) == len(graphs)
    assert len(groups) == len(graphs)
    # Labels are binary ints.
    assert set(np.unique(y)).issubset({0, 1})
    # Groups are the shape families.
    assert groups == [G.shape_family for G in graphs]
    assert X.shape[1] > 0
    assert np.isfinite(X.to_numpy()).all()
    # Columns come sorted for a stable order.
    assert list(X.columns) == sorted(X.columns)


def test_build_design_matrix_length_mismatch_raises():
    cfg = _tiny_cfg()
    graphs, labels = _tiny_cohort(cfg)
    with pytest.raises(ValueError):
        build_design_matrix(graphs, labels[:-1], cfg, np.random.default_rng(0))


def test_build_design_matrix_deterministic():
    cfg = _tiny_cfg()
    graphs, labels = _tiny_cohort(cfg)
    X1, y1, g1 = build_design_matrix(graphs, labels, cfg, np.random.default_rng(7))
    X2, y2, g2 = build_design_matrix(graphs, labels, cfg, np.random.default_rng(7))
    pd.testing.assert_frame_equal(X1, X2)
    assert np.array_equal(y1, y2)
    assert g1 == g2
