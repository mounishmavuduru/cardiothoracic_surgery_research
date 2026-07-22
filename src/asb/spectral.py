"""Eigenpairs and spectral features for atrial conduction graphs.

This module computes the low end of the graph-Laplacian spectrum (the physically
meaningful slow/global modes), the Fiedler value/vector, the Perron eigenpair of a
nonnegative adjacency matrix, and a bundle of scalar spectral descriptors used as
features downstream.

All routines are pure (no I/O, no global random state). Sparse solves use
``scipy.sparse.linalg.eigsh`` in shift-invert mode for the smallest eigenvalues
(``sigma`` set just below zero as a small regularizer so the singular combinatorial
Laplacian can be factorized) and ``which='LM'`` for the Perron eigenpair. Small
matrices fall back to a dense symmetric eigensolver, which is exact for the analytic
gate tests (path / ring Laplacians).

Notes
-----
Contract (verified in ``tests/test_spectral.py``):

* On the path graph :math:`P_N` the combinatorial-Laplacian eigenvalues are
  :math:`\\lambda_k = 2 - 2\\cos(k\\pi/N)`, :math:`k = 0, \\dots, N-1`.
* On the ring graph :math:`C_N` they are :math:`2 - 2\\cos(2\\pi k / N)`.
* The Perron vector of a nonnegative matrix is entrywise nonnegative.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

from asb.types import AtrialGraph

__all__ = [
    "smallest_eigpairs",
    "fiedler",
    "spectral_gap",
    "perron",
    "inverse_participation_ratio",
    "spectral_features",
    "cheeger_estimate",
]

# Regularizer used as the shift-invert target for the smallest eigenvalues. The
# combinatorial Laplacian is singular (lambda_0 = 0), so factorizing ``L - sigma I``
# at sigma exactly 0 is ill-posed; a tiny negative sigma keeps L - sigma I positive
# definite while still selecting the eigenvalues nearest zero.
_SHIFT_SIGMA = -1e-8
# Threshold below which an eigenvalue is treated as numerically zero.
_ZERO_TOL = 1e-9


def _as_dense_symmetric(L) -> np.ndarray:
    """Return a symmetrized dense array of ``L``."""
    A = L.toarray() if sp.issparse(L) else np.asarray(L, dtype=float)
    return 0.5 * (A + A.T)


def smallest_eigpairs(L, k: int = 6) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the ``k`` smallest eigenpairs of a symmetric matrix ``L``.

    Parameters
    ----------
    L : (n, n) array_like or scipy.sparse matrix
        Symmetric positive-semidefinite matrix (typically a graph Laplacian).
    k : int, optional
        Number of smallest eigenpairs to return, by default 6. Clipped to ``n``.

    Returns
    -------
    vals : (k,) ndarray
        Eigenvalues in ascending order.
    vecs : (n, k) ndarray
        Corresponding orthonormal eigenvectors as columns, aligned with ``vals``.

    Notes
    -----
    Small matrices (or requests for nearly the full spectrum) use a dense
    symmetric eigensolver; larger ones use ``eigsh`` in shift-invert mode with the
    shift just below zero. On the path graph :math:`P_N` this reproduces
    :math:`\\lambda_k = 2 - 2\\cos(k\\pi/N)`.
    """
    n = int(L.shape[0])
    if k < 1:
        raise ValueError("k must be >= 1")
    k = min(k, n)

    use_dense = (not sp.issparse(L)) or n <= 256 or k >= n - 1
    if not use_dense:
        try:
            # Fixed deterministic start vector so eigenvectors are reproducible
            # across processes even when the spectrum is near-degenerate.
            v0 = np.random.default_rng(0).standard_normal(n)
            vals, vecs = eigsh(L.astype(float), k=k, sigma=_SHIFT_SIGMA, which="LM", v0=v0)
            order = np.argsort(vals)
            return vals[order], vecs[:, order]
        except Exception:
            use_dense = True  # noqa: F841 (fall through to dense path)

    A = _as_dense_symmetric(L)
    vals, vecs = np.linalg.eigh(A)
    return vals[:k].copy(), vecs[:, :k].copy()


def fiedler(L) -> Tuple[float, np.ndarray]:
    """Return the Fiedler value and vector (2nd-smallest eigenpair) of ``L``.

    Parameters
    ----------
    L : (n, n) array_like or scipy.sparse matrix
        Symmetric graph Laplacian.

    Returns
    -------
    lambda2 : float
        Second-smallest eigenvalue (algebraic connectivity).
    phi2 : (n,) ndarray
        Corresponding eigenvector.
    """
    vals, vecs = smallest_eigpairs(L, k=min(2, int(L.shape[0])))
    if vals.shape[0] < 2:
        raise ValueError("fiedler requires n >= 2")
    return float(vals[1]), vecs[:, 1]


def spectral_gap(L) -> float:
    """Return the spectral gap ``lambda3 - lambda2`` of ``L``.

    Parameters
    ----------
    L : (n, n) array_like or scipy.sparse matrix
        Symmetric graph Laplacian.

    Returns
    -------
    gap : float
        Difference between the third- and second-smallest eigenvalues. Zero if
        ``n < 3``.
    """
    n = int(L.shape[0])
    if n < 3:
        return 0.0
    vals, _ = smallest_eigpairs(L, k=3)
    return float(vals[2] - vals[1])


def perron(W) -> Tuple[float, np.ndarray]:
    """Return the Perron eigenpair (dominant eigenvalue and vector) of ``W``.

    Parameters
    ----------
    W : (n, n) array_like or scipy.sparse matrix
        Symmetric nonnegative matrix (e.g. a weighted adjacency).

    Returns
    -------
    rho : float
        Spectral radius (largest eigenvalue, the Perron root).
    v : (n,) ndarray
        Corresponding eigenvector, sign-normalized to be entrywise nonnegative
        (Perron-Frobenius guarantees a nonnegative dominant eigenvector for a
        nonnegative matrix).
    """
    n = int(W.shape[0])
    if n == 1:
        val = float(_as_dense_symmetric(W)[0, 0])
        return val, np.array([1.0])

    if sp.issparse(W) and n > 256:
        try:
            # Positive start vector (Perron vector is entrywise positive) for a
            # deterministic, well-conditioned dominant-eigenpair solve.
            v0 = np.ones(n, dtype=float)
            vals, vecs = eigsh(W.astype(float), k=1, which="LM", v0=v0)
            rho, v = float(vals[0]), vecs[:, 0]
        except Exception:
            A = _as_dense_symmetric(W)
            vals, vecs = np.linalg.eigh(A)
            rho, v = float(vals[-1]), vecs[:, -1]
    else:
        A = _as_dense_symmetric(W)
        vals, vecs = np.linalg.eigh(A)
        rho, v = float(vals[-1]), vecs[:, -1]

    # Sign-normalize: pick the sign of the entry with largest magnitude so the
    # dominant (Perron) eigenvector comes out nonnegative.
    pivot = v[np.argmax(np.abs(v))]
    if pivot < 0:
        v = -v
    # Clamp tiny negative numerical noise to zero so the vector is exactly >= 0.
    v = np.where(v < 0, np.maximum(v, 0.0), v)
    return rho, v


def inverse_participation_ratio(v) -> float:
    """Inverse participation ratio ``sum(v**4) / sum(v**2)**2`` of a vector.

    Parameters
    ----------
    v : (n,) array_like
        Real vector (e.g. an eigenvector).

    Returns
    -------
    ipr : float
        IPR in ``(0, 1]``. Values near ``1/n`` indicate a delocalized vector;
        values near ``1`` indicate localization on a single node. Returns 0.0 for
        the all-zero vector.
    """
    v = np.asarray(v, dtype=float)
    s2 = float(np.sum(v ** 2))
    if s2 <= 0.0:
        return 0.0
    return float(np.sum(v ** 4) / (s2 ** 2))


def cheeger_estimate(L) -> float:
    """Estimate the Cheeger constant from the Fiedler value via its bounds.

    Cheeger's inequality bounds the isoperimetric (Cheeger) constant ``h`` by
    :math:`\\lambda_2 / 2 \\le h \\le \\sqrt{2 \\lambda_2}`. This returns the
    midpoint of those two bounds as a cheap scalar estimate.

    Parameters
    ----------
    L : (n, n) array_like or scipy.sparse matrix
        Symmetric graph Laplacian.

    Returns
    -------
    h_est : float
        Midpoint of the lower and upper Cheeger bounds.
    """
    lam2, _ = fiedler(L)
    lam2 = max(lam2, 0.0)
    lower = lam2 / 2.0
    upper = np.sqrt(2.0 * lam2)
    return float(0.5 * (lower + upper))


def _combinatorial_laplacian(G: AtrialGraph) -> sp.csr_matrix:
    """Build the combinatorial Laplacian ``L = D - W`` from a graph's adjacency."""
    W = G.adjacency()
    d = np.asarray(W.sum(axis=1)).ravel()
    L = sp.diags(d) - W
    return L.tocsr()


def spectral_features(G: AtrialGraph, *, k: int = 6, t: float = 1.0) -> dict:
    """Compute scalar spectral descriptors of an atrial conduction graph.

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph. The combinatorial Laplacian ``L = D - W`` is
        built directly from ``G.adjacency()`` (this module does not depend on
        ``asb.graph``).
    k : int, optional
        Number of low eigenvalues to use for entropy / heat-kernel summaries,
        by default 6.
    t : float, optional
        Diffusion time for the heat-kernel trace, by default 1.0.

    Returns
    -------
    features : dict
        Dictionary with keys:

        ``gap``
            Spectral gap ``lambda3 - lambda2``.
        ``spectral_entropy``
            Shannon entropy of the (normalized) low eigenvalue distribution.
        ``cheeger_estimate``
            Midpoint Cheeger-constant estimate from ``lambda2``.
        ``n_near_zero``
            Count of eigenvalues below a small tolerance (connected-component
            multiplicity of ``lambda = 0``).
        ``perron_ipr``
            Inverse participation ratio of the adjacency Perron vector (a
            localization measure).
        ``heat_kernel_trace``
            Partial heat-kernel trace ``sum_i exp(-t * lambda_i)`` over the
            computed low eigenvalues.
        ``lambda2``
            The Fiedler value itself, exposed for convenience.
    """
    n = G.n_nodes
    L = _combinatorial_laplacian(G)
    kk = min(max(k, 3), n)
    vals, _ = smallest_eigpairs(L, k=kk)

    lam2 = float(vals[1]) if kk >= 2 else 0.0
    lam3 = float(vals[2]) if kk >= 3 else lam2
    gap = float(lam3 - lam2)

    # Spectral entropy of the normalized low-eigenvalue distribution (clip tiny
    # negatives from floating point, drop the trivially-zero modes).
    lv = np.clip(vals, 0.0, None)
    total = float(lv.sum())
    if total > 0.0:
        p = lv / total
        nz = p[p > 0.0]
        spectral_entropy = float(-np.sum(nz * np.log(nz)))
    else:
        spectral_entropy = 0.0

    n_near_zero = int(np.sum(vals < _ZERO_TOL))
    heat_kernel_trace = float(np.sum(np.exp(-t * lv)))

    _, perron_v = perron(G.adjacency())
    perron_ipr = inverse_participation_ratio(perron_v)

    return {
        "gap": gap,
        "spectral_entropy": spectral_entropy,
        "cheeger_estimate": cheeger_estimate(L),
        "n_near_zero": n_near_zero,
        "perron_ipr": perron_ipr,
        "heat_kernel_trace": heat_kernel_trace,
        "lambda2": lam2,
    }
