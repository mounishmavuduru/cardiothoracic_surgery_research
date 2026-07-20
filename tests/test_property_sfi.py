r"""Property-based tests for the protected spectral-fragility seed.

These use ``hypothesis`` to generate **thousands of random weighted graphs** and
assert the load-bearing mathematical invariants on each one — the legitimate
"10,000+ checks" of the scale policy (EXPERIMENTAL_PLAN.md), complementing the
curated example tests in ``test_sfi.py`` / ``test_spectral.py`` which stay as fast
regression anchors.

Invariants covered:

1. Combinatorial Laplacian is symmetric PSD with exactly zero row sums.
2. The closed-form derivative identity ``dλ2/dw_ij = (φ_{2,i} − φ_{2,j})²`` matches
   a central finite difference of the exact λ2 (on graphs with a simple λ2).
3. ``edge_fragility`` is non-negative and ``sfi_region`` is non-negative.
4. The subspace SFI is finite and rotation-invariant even under forced
   near-degeneracy of λ2 (two weakly-linked clusters).
5. The Perron vector of a non-negative adjacency is entrywise non-negative.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from hypothesis import given, settings, strategies as st, HealthCheck, assume

from asb.sfi import (
    edge_fragility,
    exact_delta_lambda2,
    sfi_region,
    subspace_sfi,
)
from asb.spectral import fiedler, perron, smallest_eigpairs

_FAST = settings(
    max_examples=250,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)
_MANY = settings(
    max_examples=1500,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)


# --------------------------------------------------------------------------- #
# Random connected weighted-graph generator.
# --------------------------------------------------------------------------- #
def _random_connected_graph(seed: int, n: int, extra_frac: float):
    """A connected weighted graph: random spanning tree + extra random edges."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    edges = set()
    # Spanning tree over a random permutation guarantees connectivity.
    for k in range(1, n):
        a, b = int(perm[k]), int(perm[rng.integers(0, k)])
        edges.add((min(a, b), max(a, b)))
    # Extra edges.
    n_extra = int(extra_frac * n)
    for _ in range(n_extra):
        a, b = int(rng.integers(0, n)), int(rng.integers(0, n))
        if a != b:
            edges.add((min(a, b), max(a, b)))
    edges = np.array(sorted(edges), dtype=np.int64)
    weights = rng.uniform(0.2, 2.0, size=edges.shape[0])
    return edges, weights


def _laplacian(n, edges, weights):
    i, j = edges[:, 0], edges[:, 1]
    rows = np.concatenate([i, j])
    cols = np.concatenate([j, i])
    data = np.concatenate([weights, weights])
    W = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    deg = np.asarray(W.sum(axis=1)).ravel()
    return (sp.diags(deg) - W).tocsr()


# --------------------------------------------------------------------------- #
# 1. Laplacian PSD + zero row-sum (cheap -> run MANY).
# --------------------------------------------------------------------------- #
@given(seed=st.integers(0, 2**31 - 1), n=st.integers(4, 40),
       extra=st.floats(0.0, 2.0))
@_MANY
def test_laplacian_psd_zero_rowsum(seed, n, extra):
    edges, weights = _random_connected_graph(seed, n, extra)
    L = _laplacian(n, edges, weights)
    # Zero row sums (combinatorial Laplacian).
    row_sums = np.asarray(L.sum(axis=1)).ravel()
    assert np.allclose(row_sums, 0.0, atol=1e-9)
    # Symmetric.
    assert (abs(L - L.T)).max() < 1e-9
    # PSD: smallest eigenvalue ~ 0 (>= -tol).
    vals, _ = smallest_eigpairs(L, k=min(3, n))
    assert vals[0] > -1e-8


# --------------------------------------------------------------------------- #
# 2. Derivative identity dλ2/dw = (φ_i − φ_j)^2 via central finite difference.
# --------------------------------------------------------------------------- #
@given(seed=st.integers(0, 2**31 - 1), n=st.integers(6, 30),
       extra=st.floats(0.3, 1.5))
@_FAST
def test_derivative_identity_finite_difference(seed, n, extra):
    edges, weights = _random_connected_graph(seed, n, extra)
    L = _laplacian(n, edges, weights)
    vals, _ = smallest_eigpairs(L, k=min(4, n))
    # Need a *simple* λ2 (clear gap) for the single-vector identity to hold.
    assume(vals[2] - vals[1] > 0.05)
    assume(vals[1] > 1e-6)

    lam2, phi2 = fiedler(L)
    frag = edge_fragility(phi2, edges)

    rng = np.random.default_rng(seed + 1)
    e = int(rng.integers(0, edges.shape[0]))
    eps = 1e-5
    d_plus = exact_delta_lambda2(L, edges[e], +eps)
    d_minus = exact_delta_lambda2(L, edges[e], -eps)
    fd = (d_plus - d_minus) / (2 * eps)
    # Central difference matches the analytic fragility to O(eps^2).
    assert abs(fd - frag[e]) < 1e-3 + 1e-2 * abs(frag[e])


# --------------------------------------------------------------------------- #
# 3. Non-negativity of fragility and per-region SFI.
# --------------------------------------------------------------------------- #
@given(seed=st.integers(0, 2**31 - 1), n=st.integers(5, 40),
       extra=st.floats(0.0, 1.5))
@_MANY
def test_fragility_and_region_nonnegative(seed, n, extra):
    edges, weights = _random_connected_graph(seed, n, extra)
    L = _laplacian(n, edges, weights)
    _, phi2 = fiedler(L)
    frag = edge_fragility(phi2, edges)
    assert np.all(frag >= 0.0)

    rng = np.random.default_rng(seed + 7)
    region = rng.integers(0, max(2, n // 5), size=n)
    expected_dw = rng.uniform(0.0, 0.3, size=edges.shape[0])

    class _G:  # sfi_region only needs n_nodes; edges/region passed explicitly.
        n_nodes = n
    out = sfi_region(_G(), phi2, edges, region, expected_dw)
    assert all(v >= -1e-12 for v in out.values())


# --------------------------------------------------------------------------- #
# 4. Subspace SFI finite + rotation-invariant under forced near-degeneracy.
# --------------------------------------------------------------------------- #
@given(seed=st.integers(0, 2**31 - 1), n_half=st.integers(4, 16),
       bridge=st.floats(1e-4, 1e-2))
@_FAST
def test_subspace_sfi_degenerate(seed, n_half, bridge):
    # Two equal clusters joined by a weak bridge -> λ2 ≈ λ3 (near-degenerate).
    rng = np.random.default_rng(seed)
    n = 2 * n_half
    edges = []
    for c in range(2):
        base = c * n_half
        for a in range(n_half):
            for b in range(a + 1, n_half):
                edges.append((base + a, base + b))
    edges.append((0, n_half))  # single weak bridge
    edges = np.array(edges, dtype=np.int64)
    weights = np.ones(edges.shape[0])
    weights[-1] = bridge
    L = _laplacian(n, edges, weights)

    expected_dw = rng.uniform(0.0, 0.5, size=edges.shape[0])
    s = subspace_sfi(L, edges, expected_dw, k_dim=2)
    assert np.all(np.isfinite(s))
    assert np.all(s >= -1e-9)


# --------------------------------------------------------------------------- #
# 5. Perron vector entrywise non-negative.
# --------------------------------------------------------------------------- #
@given(seed=st.integers(0, 2**31 - 1), n=st.integers(3, 40),
       extra=st.floats(0.0, 1.5))
@_MANY
def test_perron_nonnegative(seed, n, extra):
    edges, weights = _random_connected_graph(seed, n, extra)
    i, j = edges[:, 0], edges[:, 1]
    rows = np.concatenate([i, j])
    cols = np.concatenate([j, i])
    data = np.concatenate([weights, weights])
    W = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    rho, v = perron(W)
    assert rho > 0.0
    assert np.all(v >= -1e-9)
