r"""Spectral Fragility Index (SFI) — the protected novel seed of AtrialSpectralBench.

This module implements the closed-form sensitivity of atrial algebraic
connectivity :math:`\lambda_2` (the Fiedler value of the weighted graph
Laplacian) to edge uncoupling, and aggregates it — under a stochastic,
fibrosis-weighted, *diffuse* perioperative-uncoupling field — into a per-region
**Spectral Fragility Index**.

The load-bearing identity
-------------------------
Because the combinatorial Laplacian factors as a sum of rank-one edge terms

.. math::
    L = \sum_{(i,j)\in E} w_{ij}\,(e_i - e_j)(e_i - e_j)^\top ,

the first-order sensitivity of any *simple* eigenvalue :math:`\lambda` with unit
eigenvector :math:`\varphi` is

.. math::
    \frac{\partial \lambda}{\partial w_{ij}}
        = \varphi^\top \frac{\partial L}{\partial w_{ij}} \varphi
        = (\varphi_i - \varphi_j)^2 .

Applied to the Fiedler pair :math:`(\lambda_2, \varphi_2)` this gives the
per-edge **fragility** :math:`(\varphi_{2,i} - \varphi_{2,j})^2`
(:func:`edge_fragility`). Reducing a conductance :math:`w_{ij}` by a small
amount :math:`\Delta w_{ij} \ge 0` therefore *drops* :math:`\lambda_2` by
:math:`\Delta w_{ij}\,(\varphi_{2,i}-\varphi_{2,j})^2` to first order, and the
per-region SFI is the expectation of that drop over the stochastic field
(:func:`sfi_region`, with the Monte-Carlo ground truth in
:func:`sfi_monte_carlo`).

Sign convention
---------------
Throughout this module a perturbation ``Delta w >= 0`` is the *magnitude* of a
conductance **reduction** (uncoupling): the edge weight becomes
``w_ij - Delta w_ij``. Consequently the reported SFI is the expected **drop**
in :math:`\lambda_2`, a non-negative fragility score, and
:func:`sfi_monte_carlo` returns the mean/std of the connectivity *drop*
``lambda2_base - lambda2_perturbed``.

Validity and the near-degenerate regime
---------------------------------------
The single-vector formula is a first-order expansion, valid only while the
perturbation is small relative to the spectral gap
:math:`\lambda_3 - \lambda_2` (:func:`validity_radius_ok`). When
:math:`\lambda_2` is near-degenerate with :math:`\lambda_3` the individual
Fiedler vector is ill-defined (an arbitrary rotation within the degenerate
subspace), but the *invariant subspace* it spans is stable; :func:`subspace_sfi`
uses the projector/trace (Davis–Kahan) generalization, which is basis
independent. For genuinely discrete edits (e.g. a Maze cut) the linear formula
is abandoned entirely and :math:`\lambda_2` is recomputed exactly
(:func:`exact_delta_lambda2`).

All functions here are pure and deterministic; every stochastic routine takes an
explicit ``rng``. This module has no dependency on sibling ``asb.*`` compute
modules — eigenpairs are computed locally with numpy/scipy.
"""
from __future__ import annotations

from typing import Dict, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

from asb.config import SFIConfig
from asb.types import AtrialGraph

__all__ = [
    "edge_fragility",
    "perturbation_field",
    "sfi_region",
    "sfi_monte_carlo",
    "validity_radius_ok",
    "subspace_sfi",
    "exact_delta_lambda2",
    "hotspot_map",
]

# Above this node count a sparse shift-invert solve is used for the low spectrum;
# below it a dense symmetric eigensolver (exact, fast, and robust on the small
# graphs used in the analytic gate tests).
_DENSE_MAX_N = 400
# Shift just below zero so the singular combinatorial Laplacian can be factorized
# by the shift-invert solver without hitting the exact null mode.
_SHIFT_SIGMA = -1e-8


# --------------------------------------------------------------------------- #
# Local linear-algebra helpers (no dependency on asb.spectral).
# --------------------------------------------------------------------------- #
def _dense_sym(L) -> np.ndarray:
    """Return a symmetrized dense array of ``L``."""
    A = L.toarray() if sp.issparse(L) else np.asarray(L, dtype=float)
    return 0.5 * (A + A.T)


def _smallest_eigpairs(L, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return the ``k`` smallest eigenpairs (ascending) of a symmetric ``L``.

    Small matrices use a dense symmetric solver; large sparse ones use
    shift-invert ``eigsh``. Eigenvectors are unit-norm columns.
    """
    n = int(L.shape[0])
    k = int(min(max(k, 1), n))
    use_dense = (not sp.issparse(L)) or n <= _DENSE_MAX_N or k >= n - 1
    if not use_dense:
        try:
            vals, vecs = eigsh(L.astype(float), k=k, sigma=_SHIFT_SIGMA, which="LM")
            order = np.argsort(vals)
            return vals[order], vecs[:, order]
        except Exception:
            use_dense = True
    A = _dense_sym(L)
    vals, vecs = np.linalg.eigh(A)
    return vals[:k].copy(), vecs[:, :k].copy()


def _lambda2(L) -> float:
    """Second-smallest eigenvalue (algebraic connectivity) of ``L``."""
    vals, _ = _smallest_eigpairs(L, k=min(2, int(L.shape[0])))
    if vals.shape[0] < 2:
        raise ValueError("lambda2 requires n >= 2")
    return float(vals[1])


def _weighted_laplacian(n: int, edges: np.ndarray, weights: np.ndarray) -> sp.csr_matrix:
    """Combinatorial Laplacian ``L = B^T diag(weights) B`` for the given edges.

    ``B`` is the signed incidence matrix, so this equals
    ``sum_e weights[e] (e_i - e_j)(e_i - e_j)^T``. Used to build both the base
    Laplacian and the additive perturbation ``Delta L``.
    """
    edges = np.asarray(edges, dtype=np.int64)
    weights = np.asarray(weights, dtype=float)
    i, j = edges[:, 0], edges[:, 1]
    rows = np.concatenate([i, j])
    cols = np.concatenate([j, i])
    data = np.concatenate([weights, weights])
    W = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    deg = np.asarray(W.sum(axis=1)).ravel()
    return (sp.diags(deg) - W).tocsr()


def _edge_regions(edges: np.ndarray, region_ids: np.ndarray) -> np.ndarray:
    """Assign each edge to a region label.

    An edge is labelled with region ``r`` iff *both* endpoints carry region
    ``r``; edges that straddle two regions are labelled ``-1`` (an inter-region
    boundary bucket). This deterministic rule is shared by :func:`sfi_region`
    and :func:`sfi_monte_carlo` so their region keys line up exactly.
    """
    region_ids = np.asarray(region_ids)
    ri = region_ids[edges[:, 0]]
    rj = region_ids[edges[:, 1]]
    return np.where(ri == rj, ri.astype(np.int64), np.int64(-1))


# --------------------------------------------------------------------------- #
# Core closed-form sensitivity.
# --------------------------------------------------------------------------- #
def edge_fragility(phi2: np.ndarray, edges: np.ndarray) -> np.ndarray:
    r"""Per-edge Fiedler fragility ``(phi2_i - phi2_j)**2 == d lambda2 / d w_ij``.

    Parameters
    ----------
    phi2 : (n,) array_like
        Fiedler vector (second-smallest Laplacian eigenvector). Any scaling is
        allowed but the closed-form derivative identity holds only for the
        unit-norm eigenvector.
    edges : (m, 2) array_like of int
        Undirected endpoint index pairs.

    Returns
    -------
    frag : (m,) ndarray
        ``(phi2_i - phi2_j)**2`` for each edge, the exact first-order derivative
        of :math:`\lambda_2` with respect to that edge's conductance.
    """
    phi2 = np.asarray(phi2, dtype=float).ravel()
    edges = np.asarray(edges, dtype=np.int64)
    d = phi2[edges[:, 0]] - phi2[edges[:, 1]]
    return d * d


# --------------------------------------------------------------------------- #
# Stochastic, fibrosis-weighted, diffuse uncoupling field.
# --------------------------------------------------------------------------- #
def perturbation_field(
    G: AtrialGraph, cfg: SFIConfig, rng: np.random.Generator
) -> Dict[str, np.ndarray]:
    r"""Per-edge mean and (diagonal) covariance of the diffuse uncoupling field.

    Models diffuse perioperative stress (inflammatory edema, transient ischemia,
    stretch) as a distribution of *small, non-negative* fractional conductance
    reductions that **grow with local fibrosis**. For edge :math:`(i,j)` with
    endpoint-mean fibrosis :math:`\bar f_{ij}`, the expected fractional drop is

    .. math::
        \text{frac}_{ij}
            = \texttt{delta\_w\_mean\_frac}\,\bigl(1 + \texttt{fibrosis\_coupling}\,
              \bar f_{ij}\bigr),

    and the expected weight reduction is
    :math:`\mathbb{E}[\Delta w_{ij}] = \text{frac}_{ij}\, w_{ij} \ge 0`. Edge
    perturbations are modelled as independent, so the covariance is diagonal with
    per-edge coefficient of variation ``delta_w_cov`` (variance
    :math:`(\texttt{delta\_w\_cov}\cdot\mathbb{E}[\Delta w_{ij}])^2`).

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph; supplies ``weights``, ``edges`` and per-node
        ``fibrosis``.
    cfg : SFIConfig
        Field parameters (``delta_w_mean_frac``, ``delta_w_cov``,
        ``fibrosis_coupling``).
    rng : numpy.random.Generator
        Present for interface uniformity. The returned *expectation* and
        *covariance* are deterministic functions of ``G`` and ``cfg`` (only the
        sampling in :func:`sfi_monte_carlo` consumes randomness), so ``rng`` is
        accepted but not required to reproduce these moments.

    Returns
    -------
    field : dict
        ``expected_dw`` : (m,) non-negative expected reduction magnitudes;
        ``cov_diag``    : (m,) per-edge variances (diagonal covariance);
        ``frac_drop``   : (m,) expected fractional drops;
        ``fibrosis_edge``: (m,) endpoint-mean fibrosis.
    """
    del rng  # moments are deterministic; sampling happens in sfi_monte_carlo.
    edges = np.asarray(G.edges, dtype=np.int64)
    w = np.asarray(G.weights, dtype=float)
    fib = np.asarray(G.fibrosis, dtype=float)

    fib_edge = 0.5 * (fib[edges[:, 0]] + fib[edges[:, 1]])
    frac = float(cfg.delta_w_mean_frac) * (1.0 + float(cfg.fibrosis_coupling) * fib_edge)
    expected_dw = np.clip(frac, 0.0, None) * w

    cov = float(cfg.delta_w_cov)
    cov_diag = (cov * expected_dw) ** 2

    return {
        "expected_dw": expected_dw,
        "cov_diag": cov_diag,
        "frac_drop": frac,
        "fibrosis_edge": fib_edge,
    }


def _sample_dw(
    expected_dw: np.ndarray,
    cov_diag: np.ndarray,
    w: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw one non-negative reduction field with the given mean/variance.

    Uses a Gamma marginal per edge (non-negative, matches mean and variance),
    then clips so a reduction never exceeds the available conductance. Because
    the first-order drop is *linear* in ``Delta w``, only the mean governs the
    expectation matched by :func:`sfi_region`; the distributional choice affects
    only the Monte-Carlo variance.
    """
    mean = np.asarray(expected_dw, dtype=float)
    var = np.asarray(cov_diag, dtype=float)
    out = np.array(mean, dtype=float, copy=True)

    pos = (mean > 0) & (var > 0)
    if np.any(pos):
        shape = mean[pos] ** 2 / var[pos]
        scale = var[pos] / mean[pos]
        out[pos] = rng.gamma(shape, scale)
    # Cannot uncouple more conductance than exists.
    return np.minimum(out, 0.999 * np.asarray(w, dtype=float))


# --------------------------------------------------------------------------- #
# Analytic per-region SFI (first-order expectation).
# --------------------------------------------------------------------------- #
def sfi_region(
    G: AtrialGraph,
    phi2: np.ndarray,
    edges: np.ndarray,
    region_ids: np.ndarray,
    expected_dw: np.ndarray,
) -> Dict[int, float]:
    r"""First-order analytic per-region Spectral Fragility Index.

    .. math::
        \mathrm{SFI}(R)
            = \sum_{(i,j)\in R} \mathbb{E}[\Delta w_{ij}]\,
              (\varphi_{2,i} - \varphi_{2,j})^2
            \;=\; \mathbb{E}\bigl[\text{drop in }\lambda_2\text{ from region }R\bigr],

    the expected drop in algebraic connectivity caused by the region's share of
    the diffuse uncoupling field, to first order.

    Parameters
    ----------
    G : AtrialGraph
        Weighted graph (supplies ``n_nodes``; ``edges``/``region_ids`` are passed
        explicitly so this stays a pure function of its arguments).
    phi2 : (n,) array_like
        Unit-norm Fiedler vector.
    edges : (m, 2) array_like of int
        Endpoint index pairs, aligned with ``expected_dw``.
    region_ids : (n,) array_like of int
        Per-node region labels. An edge belongs to region ``r`` iff both
        endpoints are in ``r`` (straddling edges fall in the ``-1`` boundary
        bucket); see :func:`_edge_regions`.
    expected_dw : (m,) array_like
        Non-negative expected reduction magnitudes (e.g.
        ``perturbation_field(...)['expected_dw']``).

    Returns
    -------
    sfi : dict[int, float]
        Region label -> ``SFI(R)`` (non-negative). Includes every region label
        present among the edges, including the ``-1`` boundary bucket if any edge
        straddles two regions.
    """
    edges = np.asarray(edges, dtype=np.int64)
    expected_dw = np.asarray(expected_dw, dtype=float)
    frag = edge_fragility(phi2, edges)
    contrib = expected_dw * frag

    er = _edge_regions(edges, region_ids)
    out: Dict[int, float] = {}
    for r in np.unique(er):
        out[int(r)] = float(contrib[er == r].sum())
    return out


# --------------------------------------------------------------------------- #
# Monte-Carlo ground truth for the analytic expectation.
# --------------------------------------------------------------------------- #
def sfi_monte_carlo(
    G: AtrialGraph, L, cfg: SFIConfig, rng: np.random.Generator
) -> Dict[Union[int, str], Dict[str, float]]:
    r"""Monte-Carlo per-region drop in :math:`\lambda_2` (ground-truth check).

    For each region, repeatedly samples the region's share of the diffuse
    uncoupling field, applies the reduction to that region's edges *exactly*
    (``L' = L - \sum_{e\in R} \Delta w_e\, s_e s_e^\top``), recomputes
    :math:`\lambda_2` exactly, and records the connectivity **drop**
    ``lambda2_base - lambda2'``. The per-region mean converges (in the small-Δw
    regime) to the first-order analytic :func:`sfi_region` value; higher-order
    (resolvent) corrections are :math:`O(\Delta w^2)` and vanish as the field
    shrinks.

    Perturbing each region in isolation makes the Monte-Carlo mean directly
    comparable to the additive analytic per-region expectation. A ``'_all'``
    entry additionally reports the drop when every edge is perturbed together.

    Parameters
    ----------
    G : AtrialGraph
        Weighted graph (supplies ``edges``, ``weights``, ``region``, ``fibrosis``).
    L : (n, n) sparse or dense
        Base combinatorial Laplacian of ``G`` (its :math:`\lambda_2` is the
        reference connectivity).
    cfg : SFIConfig
        Field + Monte-Carlo parameters (``n_monte_carlo`` draws).
    rng : numpy.random.Generator
        Explicit RNG; the only source of randomness.

    Returns
    -------
    mc : dict
        Keys are region labels (int) plus the string ``'_all'``. Each value is a
        dict ``{"mean", "std", "se", "n"}`` where ``mean``/``std`` are of the
        :math:`\lambda_2` drop across draws and ``se = std / sqrt(n)`` is the
        standard error of the mean.
    """
    edges = np.asarray(G.edges, dtype=np.int64)
    w = np.asarray(G.weights, dtype=float)
    n = int(G.n_nodes)
    n_mc = int(cfg.n_monte_carlo)

    field = perturbation_field(G, cfg, rng)
    expected_dw = field["expected_dw"]
    cov_diag = field["cov_diag"]

    lam2_base = _lambda2(L)
    er = _edge_regions(edges, G.region)

    # Region label -> boolean edge mask; plus the all-edges case.
    groups: Dict[Union[int, str], np.ndarray] = {
        int(r): (er == r) for r in np.unique(er)
    }
    groups["_all"] = np.ones(edges.shape[0], dtype=bool)

    out: Dict[Union[int, str], Dict[str, float]] = {}
    for key, mask in groups.items():
        if not np.any(mask):
            out[key] = {"mean": 0.0, "std": 0.0, "se": 0.0, "n": 0}
            continue
        e_sub = edges[mask]
        dw_mean = expected_dw[mask]
        dw_var = cov_diag[mask]
        w_sub = w[mask]

        drops = np.empty(n_mc, dtype=float)
        for t in range(n_mc):
            dw = _sample_dw(dw_mean, dw_var, w_sub, rng)
            dL = _weighted_laplacian(n, e_sub, dw)
            lam2_new = _lambda2(L - dL)
            drops[t] = lam2_base - lam2_new

        mean = float(drops.mean())
        std = float(drops.std(ddof=1)) if n_mc > 1 else 0.0
        se = std / np.sqrt(n_mc) if n_mc > 0 else 0.0
        out[key] = {"mean": mean, "std": std, "se": float(se), "n": n_mc}
    return out


# --------------------------------------------------------------------------- #
# Validity radius + subspace (near-degenerate) generalization.
# --------------------------------------------------------------------------- #
def validity_radius_ok(
    dL_norm: float, lambda2: float, lambda3: float, safety: float
) -> bool:
    r"""Whether the first-order single-vector SFI is trustworthy.

    Returns ``True`` iff ``dL_norm <= safety * (lambda3 - lambda2)`` — i.e. the
    perturbation is small relative to the spectral gap, the regime in which the
    single-Fiedler-vector expansion is valid. When the :math:`\lambda_2`
    –:math:`\lambda_3` gap is small (near-degenerate) the right-hand side
    collapses and this returns ``False``, signalling the caller to fall back to
    :func:`subspace_sfi`.

    Parameters
    ----------
    dL_norm : float
        Norm of the perturbation :math:`\|\Delta L\|` (operator or Frobenius;
        the caller chooses, consistently).
    lambda2, lambda3 : float
        Second- and third-smallest Laplacian eigenvalues.
    safety : float
        Safety factor in ``(0, 1]`` (e.g. ``SFIConfig.validity_safety``).

    Returns
    -------
    ok : bool
    """
    gap = float(lambda3) - float(lambda2)
    return bool(float(dL_norm) <= float(safety) * gap)


def subspace_sfi(
    L, edges: np.ndarray, expected_dw: np.ndarray, k_dim: int = 2
) -> np.ndarray:
    r"""Projector / Davis–Kahan subspace fragility (near-degenerate :math:`\lambda_2`).

    When :math:`\lambda_2` is (near-)degenerate the individual Fiedler vector is
    an arbitrary rotation within the low invariant subspace, so
    :func:`edge_fragility` is unstable. The trace of the perturbation restricted
    to the ``k_dim``-dimensional invariant subspace
    :math:`P = \sum_{c} \varphi_c \varphi_c^\top` is basis independent, giving a
    well-defined per-edge score

    .. math::
        s_{ij}^\top P\, s_{ij}
            = \sum_{c}(\varphi_{c,i} - \varphi_{c,j})^2 ,
        \qquad s_{ij} = e_i - e_j,

    where the sum runs over the ``k_dim`` non-trivial low modes
    :math:`\{\varphi_2, \dots, \varphi_{k\_dim+1}\}` (the null mode is skipped).
    The returned field is ``expected_dw`` times this subspace fragility — the
    first-order contribution of each edge to the drop in the *sum* of the
    clustered eigenvalues, and the correct generalization of the single-vector
    SFI when the gap closes.

    Parameters
    ----------
    L : (n, n) sparse or dense
        Base combinatorial Laplacian (assumed connected: exactly one null mode).
    edges : (m, 2) array_like of int
        Endpoint index pairs, aligned with ``expected_dw``.
    expected_dw : (m,) array_like
        Non-negative expected reduction magnitudes.
    k_dim : int, optional
        Dimension of the invariant subspace above the null mode, by default 2
        (so the cluster is :math:`\{\varphi_2, \varphi_3\}`).

    Returns
    -------
    sfi : (m,) ndarray
        Per-edge subspace SFI contributions (non-negative, always finite even
        when :math:`\lambda_2 = \lambda_3`).
    """
    edges = np.asarray(edges, dtype=np.int64)
    expected_dw = np.asarray(expected_dw, dtype=float)
    k_dim = int(max(k_dim, 1))

    vals, vecs = _smallest_eigpairs(L, k=k_dim + 1)
    # Skip the trivial null (constant) mode; keep the next k_dim modes.
    cluster = vecs[:, 1 : 1 + k_dim]

    i, j = edges[:, 0], edges[:, 1]
    diffs = cluster[i, :] - cluster[j, :]  # (m, k_dim)
    subspace_frag = np.sum(diffs * diffs, axis=1)  # rotation-invariant
    return expected_dw * subspace_frag


# --------------------------------------------------------------------------- #
# Exact recompute (honest boundary for discrete / large edits).
# --------------------------------------------------------------------------- #
def exact_delta_lambda2(
    L, edge_id: Union[Sequence[int], np.ndarray], dw: float
) -> float:
    r"""Exact change in :math:`\lambda_2` from editing one edge by ``dw``.

    Applies the rank-one update
    :math:`L' = L + dw\,(e_i - e_j)(e_i - e_j)^\top` (increasing the conductance
    of edge :math:`(i,j)` by ``dw``; pass a negative ``dw`` for a reduction /
    cut), recomputes :math:`\lambda_2` **exactly**, and returns
    :math:`\lambda_2(L') - \lambda_2(L)`. This is the honest boundary used when
    the linear formula is invalid (Maze cuts, large edits) and the reference for
    the finite-difference derivative-identity test.

    Parameters
    ----------
    L : (n, n) sparse or dense
        Base combinatorial Laplacian.
    edge_id : sequence of two int
        The node pair ``(i, j)`` identifying the edge to edit.
    dw : float
        Signed conductance change added to that edge.

    Returns
    -------
    dlambda2 : float
        Exact :math:`\lambda_2(L + dw\, s s^\top) - \lambda_2(L)`.
    """
    i, j = int(edge_id[0]), int(edge_id[1])
    n = int(L.shape[0])
    dL = _weighted_laplacian(n, np.array([[i, j]], dtype=np.int64), np.array([float(dw)]))

    A = _dense_sym(L) if not sp.issparse(L) else L
    lam2_base = _lambda2(A)
    lam2_new = _lambda2(A + dL if sp.issparse(A) else A + dL.toarray())
    return float(lam2_new - lam2_base)


# --------------------------------------------------------------------------- #
# Reentry-initiation hotspot map.
# --------------------------------------------------------------------------- #
def hotspot_map(
    G: AtrialGraph, phi2: np.ndarray, perron_v: np.ndarray
) -> np.ndarray:
    r"""Per-node reentry-initiation hotspot score.

    Combines the local Fiedler gradient with the Perron (dominant-activation)
    localization:

    .. math::
        \text{score}_i
            = \Bigl(\tfrac{1}{|\mathcal{N}(i)|}
              \sum_{j\in\mathcal{N}(i)} |\varphi_{2,i} - \varphi_{2,j}|\Bigr)
              \cdot |v_i| ,

    where the first factor is the edge-averaged :math:`|\nabla \varphi_2|` at
    node :math:`i` and :math:`v` is the Perron eigenvector of the adjacency.
    Hotspots (large score) are where a sharp Fiedler transition co-localizes with
    strong activation localization — the hypothesized reentry-initiation sites.

    Parameters
    ----------
    G : AtrialGraph
        Weighted graph (supplies ``edges`` and ``n_nodes``).
    phi2 : (n,) array_like
        Fiedler vector.
    perron_v : (n,) array_like
        Perron eigenvector of the adjacency (entrywise non-negative).

    Returns
    -------
    score : (n,) ndarray
        Non-negative per-node hotspot score (0 for isolated nodes).
    """
    phi2 = np.asarray(phi2, dtype=float).ravel()
    perron_v = np.abs(np.asarray(perron_v, dtype=float).ravel())
    edges = np.asarray(G.edges, dtype=np.int64)
    n = int(G.n_nodes)

    i, j = edges[:, 0], edges[:, 1]
    absdiff = np.abs(phi2[i] - phi2[j])

    grad_sum = np.zeros(n, dtype=float)
    count = np.zeros(n, dtype=float)
    np.add.at(grad_sum, i, absdiff)
    np.add.at(grad_sum, j, absdiff)
    np.add.at(count, i, 1.0)
    np.add.at(count, j, 1.0)

    grad_mag = np.zeros(n, dtype=float)
    nz = count > 0
    grad_mag[nz] = grad_sum[nz] / count[nz]

    return grad_mag * perron_v
