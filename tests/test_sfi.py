"""Tests for asb.sfi — the protected novel seed (Spectral Fragility Index).

Self-contained: depends only on asb.types, asb.config, asb.sfi, numpy and scipy.
Eigenpairs used as references are computed locally with numpy/scipy so the tests
never depend on the concurrently-written asb.spectral. Small graphs are built by
hand; every stochastic step uses an explicit seeded Generator.

Hard analytic gates (BUILD_SPEC §2.3):
  1. derivative identity  d(lambda2)/d w_ij == (phi2_i - phi2_j)^2   via finite diff;
  2. analytic sfi_region expectation == sfi_monte_carlo mean within MC standard error;
  3. near-degenerate lambda2 -> validity_radius_ok is False while subspace_sfi stays
     finite and basis-independent where the single-vector formula is unstable.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from asb.config import SFIConfig
from asb.types import AtrialGraph
from asb.sfi import (
    edge_fragility,
    perturbation_field,
    sfi_region,
    sfi_monte_carlo,
    validity_radius_ok,
    subspace_sfi,
    exact_delta_lambda2,
    hotspot_map,
)


# --------------------------------------------------------------------------- #
# Local helpers (no asb.spectral dependency).
# --------------------------------------------------------------------------- #
def _laplacian_dense(n, edges, weights):
    """Combinatorial Laplacian L = D - W as a dense array."""
    W = np.zeros((n, n))
    for (i, j), w in zip(edges, weights):
        W[i, j] += w
        W[j, i] += w
    return np.diag(W.sum(axis=1)) - W


def _fiedler_local(L):
    """(lambda2, phi2) via a dense symmetric eigensolver."""
    A = L.toarray() if sp.issparse(L) else np.asarray(L, dtype=float)
    A = 0.5 * (A + A.T)
    vals, vecs = np.linalg.eigh(A)
    return float(vals[1]), vecs[:, 1]


def _low_eigvals(L, k=4):
    A = L.toarray() if sp.issparse(L) else np.asarray(L, dtype=float)
    A = 0.5 * (A + A.T)
    vals = np.linalg.eigvalsh(A)
    return vals[:k]


def _random_connected_graph(seed=0, n=40, extra=60):
    """Random connected weighted AtrialGraph with a clear lambda2-lambda3 gap.

    A spanning path guarantees connectivity; extra random chords enrich the
    spectrum. Reseeds internally until lambda3 - lambda2 is comfortably positive.
    """
    base = np.random.default_rng(seed)
    for _ in range(200):
        rng = np.random.default_rng(base.integers(1 << 30))
        # Spanning path for guaranteed connectivity.
        path = np.column_stack([np.arange(n - 1), np.arange(1, n)])
        # Extra unique chords.
        chords = set()
        while len(chords) < extra:
            a, b = rng.integers(0, n, size=2)
            if a != b:
                chords.add((min(int(a), int(b)), max(int(a), int(b))))
        chords = np.array(sorted(chords - {(int(p[0]), int(p[1])) for p in path}))
        edges = np.vstack([path, chords]).astype(np.int64)
        # canonicalize i<j and drop dups
        edges = np.unique(np.sort(edges, axis=1), axis=0)
        weights = rng.uniform(0.3, 1.7, size=len(edges))
        coords = rng.uniform(size=(n, 3))
        fibrosis = rng.uniform(0.0, 1.0, size=n)
        region = rng.integers(0, 3, size=n).astype(np.int64)
        G = AtrialGraph(
            coords=coords,
            edges=edges,
            weights=weights,
            fibrosis=fibrosis,
            uac=rng.uniform(size=(n, 2)),
            region=region,
            shape_family="rand",
        )
        L = _laplacian_dense(n, edges, weights)
        v = _low_eigvals(L, k=4)
        # Connected (lambda1 ~ 0 < lambda2) with a clear lambda2-lambda3 gap.
        if v[0] < 1e-9 and v[1] > 1e-6 and (v[2] - v[1]) > 0.05:
            return G, L
    raise RuntimeError("failed to build a graph with a clear spectral gap")


def _ring_graph(n=8):
    """Unweighted cycle C_n: lambda2 == lambda3 exactly (degenerate Fiedler pair)."""
    edges = np.array([[k, (k + 1) % n] for k in range(n)], dtype=np.int64)
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    weights = np.ones(len(edges))
    coords = np.zeros((n, 3))
    G = AtrialGraph(
        coords=coords,
        edges=edges,
        weights=weights,
        fibrosis=np.linspace(0.1, 0.5, n),
        uac=np.zeros((n, 2)),
        region=np.zeros(n, dtype=np.int64),
        shape_family="ring",
    )
    L = _laplacian_dense(n, edges, weights)
    return G, L


# --------------------------------------------------------------------------- #
# 1. Derivative identity: FD d(lambda2)/dw == edge_fragility.
# --------------------------------------------------------------------------- #
def test_edge_fragility_shape_and_values():
    phi2 = np.array([0.0, 1.0, -2.0, 0.5])
    edges = np.array([[0, 1], [1, 2], [0, 3]], dtype=np.int64)
    frag = edge_fragility(phi2, edges)
    assert frag.shape == (3,)
    assert np.allclose(frag, [(0 - 1) ** 2, (1 + 2) ** 2, (0 - 0.5) ** 2])
    assert np.all(frag >= 0)


def test_derivative_identity_finite_difference():
    G, L = _random_connected_graph(seed=1, n=36, extra=50)
    lam2, phi2 = _fiedler_local(L)
    frag = edge_fragility(phi2, G.edges)

    eps = 1e-6
    # A spread of edges: some fragile, some not.
    order = np.argsort(frag)
    probe = np.unique(np.concatenate([order[:3], order[-3:], order[len(order) // 2 : len(order) // 2 + 2]]))

    for e in probe:
        i, j = int(G.edges[e, 0]), int(G.edges[e, 1])
        # exact_delta_lambda2 recomputes lambda2 exactly for the rank-1 edit.
        d_plus = exact_delta_lambda2(L, (i, j), eps)
        d_minus = exact_delta_lambda2(L, (i, j), -eps)
        central = (d_plus - d_minus) / (2 * eps)
        assert abs(central - frag[e]) < 1e-4, (e, central, frag[e])


def test_exact_delta_lambda2_matches_full_rebuild():
    G, L = _random_connected_graph(seed=5, n=30, extra=40)
    lam2_base, _ = _fiedler_local(L)
    e = 7
    i, j = int(G.edges[e, 0]), int(G.edges[e, 1])
    dw = -0.5 * G.weights[e]  # partial uncoupling of one edge

    w2 = np.array(G.weights, copy=True)
    w2[e] += dw
    L2 = _laplacian_dense(G.n_nodes, G.edges, w2)
    lam2_ref, _ = _fiedler_local(L2)

    got = exact_delta_lambda2(L, (i, j), dw)
    assert abs(got - (lam2_ref - lam2_base)) < 1e-10


# --------------------------------------------------------------------------- #
# 2. Analytic sfi_region == sfi_monte_carlo mean (small-Delta w regime).
# --------------------------------------------------------------------------- #
def test_analytic_sfi_matches_monte_carlo():
    G, L = _random_connected_graph(seed=3, n=32, extra=44)
    lam2, phi2 = _fiedler_local(L)

    # Small-Delta w regime so higher-order (O(dw^2)) bias << Monte-Carlo se.
    cfg = SFIConfig(
        n_monte_carlo=4000,
        delta_w_mean_frac=0.004,
        delta_w_cov=0.4,
        fibrosis_coupling=1.0,
    )

    field = perturbation_field(G, cfg, np.random.default_rng(0))
    expected_dw = field["expected_dw"]
    assert np.all(expected_dw >= 0)
    assert field["cov_diag"].shape == expected_dw.shape

    analytic = sfi_region(G, phi2, G.edges, G.region, expected_dw)
    mc = sfi_monte_carlo(G, L, cfg, np.random.default_rng(42))

    # Every analytic region key must appear in the MC output.
    for r, sfi_val in analytic.items():
        assert r in mc
        m = mc[r]
        se = max(m["se"], 1e-12)
        # Within ~5 standard errors of the MC mean (seeded, deterministic).
        assert abs(sfi_val - m["mean"]) <= 5.0 * se + 1e-9, (r, sfi_val, m)

    # Additivity: sum of region SFIs == analytic total == first-order '_all' mean.
    total_analytic = sum(analytic.values())
    all_mean = mc["_all"]["mean"]
    all_se = max(mc["_all"]["se"], 1e-12)
    assert abs(total_analytic - all_mean) <= 6.0 * all_se + 1e-9

    # SFI is a non-negative drop for a pure uncoupling field.
    assert all(v >= -1e-12 for v in analytic.values())


def test_sfi_region_is_first_order_derivative_weighting():
    # SFI(R) must equal sum of expected_dw * (phi_i - phi_j)^2 over R exactly.
    G, L = _random_connected_graph(seed=9, n=24, extra=30)
    _, phi2 = _fiedler_local(L)
    expected_dw = np.linspace(0.01, 0.05, G.n_edges)
    out = sfi_region(G, phi2, G.edges, G.region, expected_dw)

    frag = edge_fragility(phi2, G.edges)
    contrib = expected_dw * frag
    ri = G.region[G.edges[:, 0]]
    rj = G.region[G.edges[:, 1]]
    er = np.where(ri == rj, ri, -1)
    for r, val in out.items():
        assert np.isclose(val, contrib[er == r].sum())


# --------------------------------------------------------------------------- #
# 3. Near-degenerate lambda2: validity switch + subspace SFI well-defined.
# --------------------------------------------------------------------------- #
def test_validity_radius_flips_and_subspace_is_well_defined():
    G, L = _ring_graph(n=8)
    vals = _low_eigvals(L, k=4)
    lam2, lam3 = float(vals[1]), float(vals[2])
    gap = lam3 - lam2
    assert gap < 1e-8  # exactly degenerate Fiedler pair on a ring

    # Any real perturbation exceeds safety * (near-zero gap) -> first-order invalid.
    dL_norm = 0.05
    assert validity_radius_ok(dL_norm, lam2, lam3, safety=0.25) is False
    # A perturbation smaller than a genuine gap is accepted (sanity of the guard).
    assert validity_radius_ok(1e-3, 0.0, 1.0, safety=0.25) is True

    # Single-vector fragility is unstable: it depends on the arbitrary rotation
    # within the degenerate {phi2, phi3} subspace; subspace_sfi does not.
    A = 0.5 * (L + L.T)
    w, V = np.linalg.eigh(A)
    phi2, phi3 = V[:, 1], V[:, 2]

    expected_dw = perturbation_field(G, SFIConfig(), np.random.default_rng(0))["expected_dw"]

    sub = subspace_sfi(L, G.edges, expected_dw, k_dim=2)
    assert sub.shape == (G.n_edges,)
    assert np.all(np.isfinite(sub))
    assert np.all(sub >= -1e-12)

    # Rotate the degenerate pair by an arbitrary angle: a legitimate alternative
    # eigenbasis. Single-vector fragility changes; subspace fragility must not.
    theta = 0.7
    c, s = np.cos(theta), np.sin(theta)
    phi2_rot = c * phi2 + s * phi3
    phi3_rot = -s * phi2 + c * phi3

    frag_a = edge_fragility(phi2, G.edges)
    frag_b = edge_fragility(phi2_rot, G.edges)
    assert not np.allclose(frag_a, frag_b, atol=1e-6)  # single vector is basis-dependent

    # Subspace fragility (sum over the cluster) is rotation-invariant.
    def _cluster_frag(u, v):
        du = u[G.edges[:, 0]] - u[G.edges[:, 1]]
        dv = v[G.edges[:, 0]] - v[G.edges[:, 1]]
        return du * du + dv * dv

    inv_a = _cluster_frag(phi2, phi3)
    inv_b = _cluster_frag(phi2_rot, phi3_rot)
    assert np.allclose(inv_a, inv_b, atol=1e-10)
    # subspace_sfi encodes exactly expected_dw * (rotation-invariant cluster frag).
    assert np.allclose(sub, expected_dw * inv_a, atol=1e-8)


def test_subspace_reduces_to_single_vector_when_gap_is_clear():
    # With a clear gap the k_dim=1 subspace SFI == single-vector SFI field.
    G, L = _random_connected_graph(seed=11, n=28, extra=36)
    _, phi2 = _fiedler_local(L)
    expected_dw = np.full(G.n_edges, 0.02)

    single = expected_dw * edge_fragility(phi2, G.edges)
    sub1 = subspace_sfi(L, G.edges, expected_dw, k_dim=1)
    assert np.allclose(single, sub1, atol=1e-8)


# --------------------------------------------------------------------------- #
# Hotspot map + determinism.
# --------------------------------------------------------------------------- #
def test_hotspot_map_shape_and_nonneg():
    G, L = _random_connected_graph(seed=2, n=26, extra=34)
    _, phi2 = _fiedler_local(L)
    # Perron vector of the adjacency (nonnegative dominant eigenvector).
    A = G.adjacency().toarray()
    w, V = np.linalg.eigh(A)
    v = V[:, -1]
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
    v = np.clip(v, 0, None)

    score = hotspot_map(G, phi2, v)
    assert score.shape == (G.n_nodes,)
    assert np.all(score >= 0)
    assert np.all(np.isfinite(score))
    # A node with zero Perron weight or zero gradient contributes zero.
    assert np.isclose(hotspot_map(G, np.zeros(G.n_nodes), v).sum(), 0.0)


def test_perturbation_field_deterministic_moments():
    G, _ = _random_connected_graph(seed=4, n=20, extra=24)
    cfg = SFIConfig(delta_w_mean_frac=0.1, delta_w_cov=0.5, fibrosis_coupling=1.0)
    f1 = perturbation_field(G, cfg, np.random.default_rng(1))
    f2 = perturbation_field(G, cfg, np.random.default_rng(999))
    # Moments are deterministic (independent of rng).
    assert np.allclose(f1["expected_dw"], f2["expected_dw"])
    assert np.allclose(f1["cov_diag"], f2["cov_diag"])
    # dw grows with local fibrosis: higher endpoint fibrosis -> larger frac drop.
    order = np.argsort(f1["fibrosis_edge"])
    assert f1["frac_drop"][order[0]] <= f1["frac_drop"][order[-1]] + 1e-12
    assert np.all(f1["expected_dw"] >= 0)


def test_monte_carlo_reproducible():
    G, L = _random_connected_graph(seed=6, n=22, extra=26)
    cfg = SFIConfig(n_monte_carlo=64, delta_w_mean_frac=0.01, delta_w_cov=0.4)
    a = sfi_monte_carlo(G, L, cfg, np.random.default_rng(7))
    b = sfi_monte_carlo(G, L, cfg, np.random.default_rng(7))
    for k in a:
        assert np.isclose(a[k]["mean"], b[k]["mean"])
        assert np.isclose(a[k]["std"], b[k]["std"])
