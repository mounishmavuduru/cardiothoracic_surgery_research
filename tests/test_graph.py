"""Tests for asb.graph — Laplacians and weighted-graph construction.

Self-contained: depends only on asb.types, asb.graph, numpy and scipy. Builds a
tiny AtrialGraph / AtrialMesh by hand; imports no sibling asb.* compute modules.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from asb.types import AtrialGraph, AtrialMesh
from asb.graph import (
    laplacian,
    laplacian_from_adjacency,
    edge_weights_from_fibres,
    cut_edges,
)


def _toy_graph(seed: int = 0) -> AtrialGraph:
    """Small connected weighted graph on 6 nodes with positive weights."""
    rng = np.random.default_rng(seed)
    edges = np.array(
        [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [0, 5], [1, 4], [0, 3]],
        dtype=np.int64,
    )
    weights = rng.uniform(0.2, 1.5, size=len(edges))
    n = 6
    coords = rng.uniform(size=(n, 3))
    return AtrialGraph(
        coords=coords,
        edges=edges,
        weights=weights,
        fibrosis=rng.uniform(0, 1, size=n),
        uac=rng.uniform(0, 1, size=(n, 2)),
        region=np.zeros(n, dtype=np.int64),
        shape_family="toy",
    )


def test_combinatorial_symmetric_psd_zero_rowsum():
    G = _toy_graph()
    L = laplacian(G, "combinatorial")
    A = L.toarray()

    # Symmetric.
    assert np.allclose(A, A.T, atol=1e-12)

    # Exact zero row sums (and column sums).
    assert np.allclose(A.sum(axis=1), 0.0, atol=1e-10)
    assert np.allclose(A.sum(axis=0), 0.0, atol=1e-10)

    # Positive semi-definite: all eigenvalues >= 0.
    evals = np.linalg.eigvalsh(A)
    assert evals.min() >= -1e-9
    # Connected graph -> exactly one zero eigenvalue.
    assert np.isclose(evals[0], 0.0, atol=1e-9)
    assert evals[1] > 1e-6


def test_sym_normalized_spectrum_in_0_2():
    G = _toy_graph(seed=3)
    L = laplacian(G, "sym")
    A = L.toarray()
    assert np.allclose(A, A.T, atol=1e-10)
    evals = np.linalg.eigvalsh(A)
    assert evals.min() >= -1e-9
    assert evals.max() <= 2.0 + 1e-9


def test_rw_laplacian_eigenvalues_in_0_2():
    G = _toy_graph(seed=7)
    L = laplacian(G, "rw")
    # Random-walk Laplacian is similar to the symmetric one -> real spectrum in [0, 2].
    evals = np.linalg.eigvals(L.toarray())
    assert np.allclose(evals.imag, 0.0, atol=1e-8)
    real = np.sort(evals.real)
    assert real.min() >= -1e-8
    assert real.max() <= 2.0 + 1e-8


def test_laplacian_from_adjacency_matches_graph():
    G = _toy_graph(seed=1)
    W = G.adjacency()
    for kind in ("combinatorial", "sym", "rw"):
        a = laplacian(G, kind).toarray()
        b = laplacian_from_adjacency(W, kind).toarray()
        assert np.allclose(a, b, atol=1e-12)


def test_isolated_node_handled():
    # Node 2 has no edges.
    edges = np.array([[0, 1], [0, 3]], dtype=np.int64)
    weights = np.array([1.0, 0.5])
    n = 4
    G = AtrialGraph(
        coords=np.zeros((n, 3)),
        edges=edges,
        weights=weights,
        fibrosis=np.zeros(n),
        uac=np.zeros((n, 2)),
        region=np.zeros(n, dtype=np.int64),
        shape_family="iso",
    )
    Lc = laplacian(G, "combinatorial").toarray()
    assert np.allclose(Lc[2], 0.0)  # zero row for isolated node
    Ls = laplacian(G, "sym").toarray()
    assert np.isclose(Ls[2, 2], 1.0)  # identity row for isolated node
    assert np.all(np.isfinite(Ls))


def test_invalid_kind_raises():
    G = _toy_graph()
    try:
        laplacian(G, "nonsense")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for invalid kind")


def _toy_mesh() -> AtrialMesh:
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
    )
    faces = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64)
    # Fibres all point along +x.
    fibres = np.tile(np.array([1.0, 0.0, 0.0]), (4, 1))
    n = 4
    return AtrialMesh(
        points=points,
        faces=faces,
        fibres=fibres,
        uac=np.zeros((n, 2)),
        fibrosis=np.zeros(n),
        region=np.zeros(n, dtype=np.int64),
        shape_family="m",
    )


def test_edge_weights_anisotropy_and_fibrosis():
    mesh = _toy_mesh()
    # Edge 0-1 runs along +x (fibre aligned); edge 0-2 runs along +y (cross-fibre).
    edges = np.array([[0, 1], [0, 2]], dtype=np.int64)
    w = edge_weights_from_fibres(mesh, edges, along=1.0, cross=0.3, fibrosis_floor=0.05)
    assert w.shape == (2,)
    assert np.isclose(w[0], 1.0)  # perfectly along fibre
    assert np.isclose(w[1], 0.3)  # perfectly cross fibre
    assert w[0] > w[1]

    # With full fibrosis on both endpoints of edge 0-1, factor floors.
    mesh.fibrosis[:] = 1.0
    w2 = edge_weights_from_fibres(mesh, edges, along=1.0, cross=0.3, fibrosis_floor=0.05)
    assert np.isclose(w2[0], 1.0 * 0.05)
    assert np.all(w2 > 0)


def test_cut_edges_zeroes_weights():
    G = _toy_graph(seed=2)
    orig = np.array(G.weights, copy=True)
    ids = [1, 4]
    H = cut_edges(G, ids)
    # Chosen weights zeroed.
    assert np.allclose(H.weights[ids], 0.0)
    # Others untouched.
    keep = [k for k in range(G.n_edges) if k not in ids]
    assert np.allclose(H.weights[keep], orig[keep])
    # Original graph unmodified (copy semantics).
    assert np.allclose(G.weights, orig)
    # Adjacency reflects the cut.
    A = H.adjacency().toarray()
    i, j = G.edges[1]
    assert np.isclose(A[i, j], 0.0)


def test_cut_edges_empty():
    G = _toy_graph()
    H = cut_edges(G, [])
    assert np.allclose(H.weights, G.weights)


def test_smallest_eig_via_scipy():
    # Sanity: eigsh on the combinatorial Laplacian returns a zero smallest eigenvalue.
    G = _toy_graph(seed=5)
    L = laplacian(G, "combinatorial")
    vals = spla.eigsh(L, k=2, which="SM", return_eigenvectors=False)
    vals = np.sort(vals)
    assert np.isclose(vals[0], 0.0, atol=1e-7)
    assert vals[1] > 0
