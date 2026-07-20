"""Tests for :mod:`asb.baselines`.

Inputs are hand-built tiny :class:`asb.types.AtrialGraph` instances so the suite
depends only on numpy/scipy/networkx, ``asb.types`` and ``asb.baselines`` (plus
``asb.spectral`` at runtime for ``lambda2_alone``). A local ``lambda2`` check is
computed independently via a dense eigensolver.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from asb.types import AtrialGraph
from asb.baselines import (
    fibrosis_burden,
    fibrosis_spatial_entropy,
    fibrosis_patch_size,
    min_cut_value,
    percolation_threshold,
    lambda2_alone,
    baseline_feature_vector,
)


# --------------------------------------------------------------------------- #
# Tiny graph builders
# --------------------------------------------------------------------------- #
def _line_graph(n: int, fibrosis=None, weights=None) -> AtrialGraph:
    """Path graph P_n with unit weights unless overridden."""
    edges = np.array([[i, i + 1] for i in range(n - 1)], dtype=np.int64)
    if weights is None:
        weights = np.ones(n - 1)
    if fibrosis is None:
        fibrosis = np.zeros(n)
    coords = np.zeros((n, 3))
    coords[:, 0] = np.arange(n)
    return AtrialGraph(
        coords=coords,
        edges=edges,
        weights=np.asarray(weights, dtype=float),
        fibrosis=np.asarray(fibrosis, dtype=float),
        uac=np.zeros((n, 2)),
        region=np.zeros(n, dtype=np.int64),
        shape_family="line",
    )


def _clique_edges(nodes):
    return [(int(a), int(b)) for i, a in enumerate(nodes) for b in nodes[i + 1 :]]


def _barbell_graph(k: int = 4, bridge_weight: float = 0.01) -> AtrialGraph:
    """Two k-cliques joined by a single thin bridge edge."""
    left = list(range(k))
    right = list(range(k, 2 * k))
    edge_list = _clique_edges(left) + _clique_edges(right)
    weights = [1.0] * len(edge_list)
    edge_list.append((k - 1, k))  # bridge
    weights.append(bridge_weight)
    edges = np.array([sorted(e) for e in edge_list], dtype=np.int64)
    n = 2 * k
    coords = np.zeros((n, 3))
    coords[:k, 0] = 0.0
    coords[k:, 0] = 5.0
    return AtrialGraph(
        coords=coords,
        edges=edges,
        weights=np.asarray(weights, dtype=float),
        fibrosis=np.zeros(n),
        uac=np.zeros((n, 2)),
        region=np.zeros(n, dtype=np.int64),
        shape_family="barbell",
    )


def _local_lambda2(G: AtrialGraph) -> float:
    """Independent dense computation of the Fiedler value."""
    W = G.adjacency().toarray()
    d = W.sum(axis=1)
    L = np.diag(d) - W
    vals = np.linalg.eigvalsh(0.5 * (L + L.T))
    return float(np.sort(vals)[1])


# --------------------------------------------------------------------------- #
# fibrosis_burden
# --------------------------------------------------------------------------- #
def test_fibrosis_burden_clean_is_zero():
    G = _line_graph(6, fibrosis=np.zeros(6))
    assert fibrosis_burden(G) == 0.0


def test_fibrosis_burden_all_fibrotic_is_one():
    G = _line_graph(6, fibrosis=np.ones(6))
    assert abs(fibrosis_burden(G) - 1.0) < 1e-12


def test_fibrosis_burden_mixed_in_range():
    G = _line_graph(6, fibrosis=np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]))
    b = fibrosis_burden(G)
    assert 0.0 <= b <= 1.0
    assert abs(b - 0.5) < 1e-12


# --------------------------------------------------------------------------- #
# fibrosis_spatial_entropy
# --------------------------------------------------------------------------- #
def test_spatial_entropy_uniform_field_is_low():
    G = _line_graph(8, fibrosis=np.full(8, 0.3))
    e = fibrosis_spatial_entropy(G)
    assert 0.0 <= e <= 1.0
    assert e < 1e-9  # all mass in one bin


def test_spatial_entropy_spread_field_is_higher():
    G = _line_graph(16, fibrosis=np.linspace(0.0, 1.0, 16))
    e = fibrosis_spatial_entropy(G, n_bins=16)
    assert 0.0 <= e <= 1.0
    assert e > 0.5


# --------------------------------------------------------------------------- #
# fibrosis_patch_size
# --------------------------------------------------------------------------- #
def test_patch_size_no_fibrosis_is_zero():
    G = _line_graph(6, fibrosis=np.zeros(6))
    assert fibrosis_patch_size(G) == 0.0


def test_patch_size_single_patch():
    # Nodes 1,2,3 fibrotic and connected on the path -> one patch of size 3.
    fib = np.array([0.0, 1.0, 1.0, 1.0, 0.0, 0.0])
    G = _line_graph(6, fibrosis=fib)
    assert fibrosis_patch_size(G) == 3.0


def test_patch_size_two_patches_mean():
    # Two separated fibrotic patches: {0,1} and {4,5} -> mean size 2.
    fib = np.array([1.0, 1.0, 0.0, 0.0, 1.0, 1.0])
    G = _line_graph(6, fibrosis=fib)
    assert fibrosis_patch_size(G) == 2.0


# --------------------------------------------------------------------------- #
# min_cut_value
# --------------------------------------------------------------------------- #
def test_min_cut_barbell_is_small():
    G = _barbell_graph(k=4, bridge_weight=0.01)
    cut = min_cut_value(G)
    assert np.isfinite(cut)
    assert abs(cut - 0.01) < 1e-9  # cutting the thin bridge separates the cliques


def test_min_cut_line_is_min_edge_weight():
    G = _line_graph(5, weights=np.array([1.0, 1.0, 0.5, 1.0]))
    cut = min_cut_value(G)
    assert abs(cut - 0.5) < 1e-9


def test_min_cut_trivial_graph_zero():
    G = _line_graph(2, weights=np.array([2.0]))
    assert min_cut_value(G) == 2.0


# --------------------------------------------------------------------------- #
# percolation_threshold
# --------------------------------------------------------------------------- #
def test_percolation_in_unit_interval():
    G = _barbell_graph(k=5)
    rng = np.random.default_rng(0)
    p = percolation_threshold(G, n_steps=20, rng=rng)
    assert np.isfinite(p)
    assert 0.0 <= p <= 1.0


def test_percolation_deterministic_for_fixed_rng():
    G = _barbell_graph(k=5)
    p1 = percolation_threshold(G, n_steps=20, rng=np.random.default_rng(42))
    p2 = percolation_threshold(G, n_steps=20, rng=np.random.default_rng(42))
    assert p1 == p2


# --------------------------------------------------------------------------- #
# lambda2_alone
# --------------------------------------------------------------------------- #
def test_lambda2_matches_local_dense():
    G = _line_graph(7)
    lam = lambda2_alone(G)
    assert np.isfinite(lam)
    assert abs(lam - _local_lambda2(G)) < 1e-8


def test_lambda2_positive_for_connected_graph():
    G = _barbell_graph(k=4)
    assert lambda2_alone(G) > 0.0


# --------------------------------------------------------------------------- #
# baseline_feature_vector
# --------------------------------------------------------------------------- #
def test_feature_vector_keys_and_finite():
    G = _barbell_graph(k=4)
    fv = baseline_feature_vector(G)
    expected = {
        "fibrosis_burden",
        "fibrosis_spatial_entropy",
        "fibrosis_patch_size",
        "min_cut_value",
        "percolation_threshold",
        "lambda2_alone",
    }
    assert set(fv.keys()) == expected
    for v in fv.values():
        assert np.isfinite(v)


def test_feature_vector_deterministic():
    G = _barbell_graph(k=4)
    fv1 = baseline_feature_vector(G)
    fv2 = baseline_feature_vector(G)
    assert fv1 == fv2
