"""Graph Laplacians and weighted-graph construction for AtrialSpectralBench.

This module turns an :class:`asb.types.AtrialGraph` (or a raw weighted adjacency)
into the graph Laplacians that the spectral-fragility machinery operates on, and
provides small helpers for building conduction edge weights from a fibre field
and for applying discrete surgical cuts (e.g. Maze lesions).

All functions are pure and deterministic; no I/O, no global RNG state.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Iterable

import numpy as np
import scipy.sparse as sp

from asb.types import AtrialGraph, AtrialMesh

__all__ = [
    "laplacian",
    "laplacian_from_adjacency",
    "edge_weights_from_fibres",
    "cut_edges",
]

_VALID_KINDS = ("combinatorial", "sym", "rw")


def laplacian_from_adjacency(W, kind: str = "combinatorial") -> sp.csr_matrix:
    """Build a graph Laplacian from a weighted adjacency matrix.

    Parameters
    ----------
    W : array_like or scipy.sparse matrix, shape (n, n)
        Symmetric, nonnegative weighted adjacency. The diagonal is ignored
        (self-loops do not contribute to the Laplacian).
    kind : {"combinatorial", "sym", "rw"}, optional
        Laplacian variant:

        - ``"combinatorial"`` : ``L = D - W``.
        - ``"sym"``           : ``L = I - D^{-1/2} W D^{-1/2}`` (symmetric-normalized).
        - ``"rw"``            : ``L = I - D^{-1} W`` (random-walk normalized).

        Isolated nodes (zero degree) contribute an identity row for the
        normalized kinds and a zero row for the combinatorial kind.

    Returns
    -------
    scipy.sparse.csr_matrix, shape (n, n)
        The requested Laplacian.

    Notes
    -----
    For ``kind="combinatorial"`` the result is symmetric, positive
    semi-definite, and has exactly zero row sums. For ``kind="sym"`` the
    spectrum lies in ``[0, 2]``.
    """
    if kind not in _VALID_KINDS:
        raise ValueError(f"kind must be one of {_VALID_KINDS}, got {kind!r}")

    W = sp.csr_matrix(W).astype(float)
    if W.shape[0] != W.shape[1]:
        raise ValueError(f"adjacency must be square, got shape {W.shape}")

    # Drop the diagonal so self-loops never enter the degree / Laplacian.
    W = W - sp.diags(W.diagonal())
    W.eliminate_zeros()

    n = W.shape[0]
    deg = np.asarray(W.sum(axis=1)).ravel()

    if kind == "combinatorial":
        L = sp.diags(deg) - W
        return sp.csr_matrix(L)

    identity = sp.identity(n, format="csr", dtype=float)

    if kind == "rw":
        inv_deg = np.zeros(n)
        nz = deg > 0
        inv_deg[nz] = 1.0 / deg[nz]
        Dinv = sp.diags(inv_deg)
        L = identity - Dinv @ W
        return sp.csr_matrix(L)

    # kind == "sym"
    inv_sqrt = np.zeros(n)
    nz = deg > 0
    inv_sqrt[nz] = 1.0 / np.sqrt(deg[nz])
    Dinv_sqrt = sp.diags(inv_sqrt)
    L = identity - Dinv_sqrt @ W @ Dinv_sqrt
    # Symmetrize to remove tiny floating-point asymmetry from the products.
    L = 0.5 * (L + L.T)
    return sp.csr_matrix(L)


def laplacian(G: AtrialGraph, kind: str = "combinatorial") -> sp.csr_matrix:
    """Graph Laplacian of an :class:`~asb.types.AtrialGraph`.

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph.
    kind : {"combinatorial", "sym", "rw"}, optional
        Laplacian variant; see :func:`laplacian_from_adjacency`.

    Returns
    -------
    scipy.sparse.csr_matrix, shape (n_nodes, n_nodes)
        The requested Laplacian, built from ``G.adjacency()``.
    """
    return laplacian_from_adjacency(G.adjacency(), kind)


def edge_weights_from_fibres(
    mesh: AtrialMesh,
    edges: np.ndarray,
    *,
    along: float = 1.0,
    cross: float = 0.3,
    fibrosis_floor: float = 0.05,
) -> np.ndarray:
    """Conduction weights for mesh edges from the local fibre field.

    Conductance is anisotropic: high when the edge runs along the local fibre
    direction, low when it runs across it. It is then attenuated by local
    fibrosis via a ``(1 - fibrosis)`` factor with a floor so that heavily
    fibrotic tissue keeps a small residual conductance rather than zero.

    Parameters
    ----------
    mesh : AtrialMesh
        Surface carrying ``points`` (n, 3), ``fibres`` (n, 3) and
        ``fibrosis`` (n,).
    edges : array_like, shape (m, 2)
        Undirected endpoint index pairs.
    along : float, optional
        Conductance for a perfectly fibre-aligned edge.
    cross : float, optional
        Conductance for a perfectly cross-fibre edge.
    fibrosis_floor : float, optional
        Lower bound on the ``(1 - fibrosis)`` attenuation factor, in ``[0, 1]``.

    Returns
    -------
    np.ndarray, shape (m,)
        Positive edge conductances.

    Notes
    -----
    Let ``a = |cos theta|`` be the absolute alignment between the unit edge
    vector and the mean unit fibre direction of its endpoints. The anisotropic
    base conductance is ``cross + (along - cross) * a**2``. The fibrosis factor
    is ``clip(1 - mean_fibrosis, fibrosis_floor, 1)``.
    """
    edges = np.asarray(edges, dtype=np.int64)
    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError(f"edges must have shape (m, 2), got {edges.shape}")

    pts = np.asarray(mesh.points, dtype=float)
    fib = np.asarray(mesh.fibres, dtype=float)
    fibrosis = np.asarray(mesh.fibrosis, dtype=float)

    i, j = edges[:, 0], edges[:, 1]

    # Unit edge directions.
    edge_vec = pts[j] - pts[i]
    edge_len = np.linalg.norm(edge_vec, axis=1)
    safe_len = np.where(edge_len > 0, edge_len, 1.0)
    edge_hat = edge_vec / safe_len[:, None]

    # Mean fibre direction of the two endpoints, unit-normalized.
    fib_mean = 0.5 * (fib[i] + fib[j])
    fib_norm = np.linalg.norm(fib_mean, axis=1)
    safe_fnorm = np.where(fib_norm > 0, fib_norm, 1.0)
    fib_hat = fib_mean / safe_fnorm[:, None]

    align = np.abs(np.sum(edge_hat * fib_hat, axis=1))
    # Degenerate directions -> isotropic midpoint conductance.
    align = np.where((edge_len > 0) & (fib_norm > 0), align, 0.0)

    base = cross + (along - cross) * align**2

    fib_edge = 0.5 * (fibrosis[i] + fibrosis[j])
    factor = np.clip(1.0 - fib_edge, fibrosis_floor, 1.0)

    return base * factor


def cut_edges(G: AtrialGraph, edge_ids: Iterable[int]) -> AtrialGraph:
    """Return a copy of ``G`` with the given edges surgically cut.

    Models a discrete surgical lesion (e.g. a Maze incision): the selected
    edges keep their topology but have their conductance set to zero, severing
    conduction across them.

    Parameters
    ----------
    G : AtrialGraph
        Source graph (left unmodified).
    edge_ids : iterable of int
        Indices into ``G.edges`` / ``G.weights`` to cut.

    Returns
    -------
    AtrialGraph
        A copy with the chosen edge weights set to zero. All other fields are
        deep-copied so the original is untouched.
    """
    edge_ids = np.asarray(list(edge_ids), dtype=np.int64)
    weights = np.array(G.weights, dtype=float, copy=True)
    if edge_ids.size:
        if edge_ids.min() < 0 or edge_ids.max() >= weights.shape[0]:
            raise IndexError("edge_ids out of range for G.weights")
        weights[edge_ids] = 0.0

    return replace(
        G,
        coords=np.array(G.coords, copy=True),
        edges=np.array(G.edges, copy=True),
        weights=weights,
        fibrosis=np.array(G.fibrosis, copy=True),
        uac=np.array(G.uac, copy=True),
        region=np.array(G.region, copy=True),
        meta=dict(G.meta),
    )
