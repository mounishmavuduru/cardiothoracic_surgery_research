"""Analytic gate tests for :mod:`asb.spectral`.

These tests depend only on numpy/scipy, ``asb.types``, and ``asb.spectral`` -- the
path/ring Laplacians are constructed locally so the suite does not import the
concurrently-written ``asb.graph`` module.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from asb.types import AtrialGraph
from asb.spectral import (
    smallest_eigpairs,
    fiedler,
    spectral_gap,
    perron,
    inverse_participation_ratio,
    spectral_features,
    cheeger_estimate,
)


# --------------------------------------------------------------------------- #
# Local Laplacian constructors (no dependency on asb.graph)
# --------------------------------------------------------------------------- #
def path_laplacian(n: int) -> sp.csr_matrix:
    """Combinatorial Laplacian of the path graph P_n."""
    rows = np.arange(n - 1)
    i = np.concatenate([rows, rows + 1])
    j = np.concatenate([rows + 1, rows])
    W = sp.csr_matrix((np.ones(2 * (n - 1)), (i, j)), shape=(n, n))
    d = np.asarray(W.sum(axis=1)).ravel()
    return (sp.diags(d) - W).tocsr()


def ring_laplacian(n: int) -> sp.csr_matrix:
    """Combinatorial Laplacian of the cycle graph C_n."""
    a = np.arange(n)
    b = (a + 1) % n
    i = np.concatenate([a, b])
    j = np.concatenate([b, a])
    W = sp.csr_matrix((np.ones(2 * n), (i, j)), shape=(n, n))
    d = np.asarray(W.sum(axis=1)).ravel()
    return (sp.diags(d) - W).tocsr()


def path_adjacency(n: int) -> sp.csr_matrix:
    """Adjacency of the path graph P_n."""
    rows = np.arange(n - 1)
    i = np.concatenate([rows, rows + 1])
    j = np.concatenate([rows + 1, rows])
    return sp.csr_matrix((np.ones(2 * (n - 1)), (i, j)), shape=(n, n))


# --------------------------------------------------------------------------- #
# Analytic-gate tests
# --------------------------------------------------------------------------- #
def test_path_laplacian_eigenvalues_closed_form():
    N = 40
    L = path_laplacian(N)
    k = 8
    vals, vecs = smallest_eigpairs(L, k=k)

    analytic = np.sort(2.0 - 2.0 * np.cos(np.arange(N) * np.pi / N))[:k]
    assert vals.shape == (k,)
    assert vecs.shape == (N, k)
    np.testing.assert_allclose(vals, analytic, atol=1e-9)
    # ascending
    assert np.all(np.diff(vals) >= -1e-12)


def test_ring_laplacian_eigenvalues_closed_form():
    N = 40
    L = ring_laplacian(N)
    k = 8
    vals, _ = smallest_eigpairs(L, k=k)

    analytic = np.sort(2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(N) / N))[:k]
    np.testing.assert_allclose(vals, analytic, atol=1e-9)


def test_eigvectors_orthonormal():
    N = 40
    L = path_laplacian(N)
    _, V = smallest_eigpairs(L, k=6)
    gram = V.T @ V
    np.testing.assert_allclose(gram, np.eye(6), atol=1e-9)


def test_fiedler_matches_smallest_eigpairs():
    N = 40
    L = path_laplacian(N)
    lam2, phi2 = fiedler(L)
    analytic = 2.0 - 2.0 * np.cos(np.pi / N)
    assert abs(lam2 - analytic) < 1e-9
    assert phi2.shape == (N,)
    # Fiedler vector of a path is monotone (up to sign): sort by index is monotone.
    signed = phi2 * np.sign(phi2[np.argmax(np.abs(phi2))])
    assert np.all(np.diff(signed) > -1e-6) or np.all(np.diff(signed) < 1e-6)


def test_spectral_gap_path():
    N = 40
    L = path_laplacian(N)
    lam2 = 2.0 - 2.0 * np.cos(np.pi / N)
    lam3 = 2.0 - 2.0 * np.cos(2.0 * np.pi / N)
    assert abs(spectral_gap(L) - (lam3 - lam2)) < 1e-9


def test_first_eigenvalue_zero_and_constant_vector():
    N = 40
    L = path_laplacian(N)
    vals, vecs = smallest_eigpairs(L, k=3)
    assert abs(vals[0]) < 1e-9
    # first eigenvector is constant
    c = vecs[:, 0]
    assert np.allclose(c, c[0], atol=1e-8)


def test_perron_nonnegative_path():
    N = 40
    W = path_adjacency(N)
    rho, v = perron(W)
    # spectral radius of P_N adjacency is 2 cos(pi/(N+1))
    assert rho > 0
    assert v.shape == (N,)
    assert np.all(v >= 0.0)
    # dominant eigenvalue is the largest in magnitude
    assert abs(rho - 2.0 * np.cos(np.pi / (N + 1))) < 1e-8


def test_perron_nonnegative_ring():
    N = 40
    a = np.arange(N)
    b = (a + 1) % N
    i = np.concatenate([a, b])
    j = np.concatenate([b, a])
    W = sp.csr_matrix((np.ones(2 * N), (i, j)), shape=(N, N))
    rho, v = perron(W)
    assert np.all(v >= 0.0)
    assert abs(rho - 2.0) < 1e-8  # ring adjacency spectral radius = 2


def test_inverse_participation_ratio_bounds():
    # delocalized (uniform) vector -> IPR = 1/n
    n = 25
    u = np.ones(n)
    assert abs(inverse_participation_ratio(u) - 1.0 / n) < 1e-12
    # localized (single spike) -> IPR = 1
    e = np.zeros(n)
    e[3] = 5.0
    assert abs(inverse_participation_ratio(e) - 1.0) < 1e-12
    # zero vector -> 0.0
    assert inverse_participation_ratio(np.zeros(n)) == 0.0


def test_cheeger_estimate_between_bounds():
    N = 40
    L = ring_laplacian(N)
    lam2, _ = fiedler(L)
    lower = lam2 / 2.0
    upper = np.sqrt(2.0 * lam2)
    h = cheeger_estimate(L)
    assert lower <= h <= upper


# --------------------------------------------------------------------------- #
# spectral_features on an AtrialGraph
# --------------------------------------------------------------------------- #
def _path_graph(n: int) -> AtrialGraph:
    coords = np.zeros((n, 3))
    coords[:, 0] = np.arange(n)
    edges = np.column_stack([np.arange(n - 1), np.arange(1, n)]).astype(np.int64)
    weights = np.ones(n - 1)
    return AtrialGraph(
        coords=coords,
        edges=edges,
        weights=weights,
        fibrosis=np.zeros(n),
        uac=np.zeros((n, 2)),
        region=np.zeros(n, dtype=np.int64),
        shape_family="test_path",
    )


def test_spectral_features_keys_and_values():
    G = _path_graph(40)
    feats = spectral_features(G)
    expected = {
        "gap",
        "spectral_entropy",
        "cheeger_estimate",
        "n_near_zero",
        "perron_ipr",
        "heat_kernel_trace",
        "lambda2",
    }
    assert expected <= set(feats)
    # connected graph -> exactly one near-zero eigenvalue
    assert feats["n_near_zero"] == 1
    assert feats["gap"] >= 0.0
    assert feats["spectral_entropy"] >= 0.0
    assert 0.0 < feats["perron_ipr"] <= 1.0
    assert feats["heat_kernel_trace"] > 0.0
    lam2_analytic = 2.0 - 2.0 * np.cos(np.pi / 40)
    assert abs(feats["lambda2"] - lam2_analytic) < 1e-9


def test_spectral_features_deterministic():
    G = _path_graph(30)
    f1 = spectral_features(G)
    f2 = spectral_features(G)
    for kk in f1:
        assert f1[kk] == f2[kk]


def test_smallest_eigpairs_sparse_large_path():
    # exercise the eigsh shift-invert branch (n > 256)
    N = 400
    L = path_laplacian(N)
    vals, _ = smallest_eigpairs(L, k=5)
    analytic = np.sort(2.0 - 2.0 * np.cos(np.arange(N) * np.pi / N))[:5]
    np.testing.assert_allclose(vals, analytic, atol=1e-7)
