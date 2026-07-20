"""Competitor / baseline features for AtrialSpectralBench.

These are the "strawman" descriptors that the Spectral Fragility Index (SFI) must
beat: simple fibrosis-burden statistics and classical graph-connectivity measures
(global min-cut, percolation threshold, the Fiedler value on its own). They are
assembled into a single named feature vector by :func:`baseline_feature_vector` and
consumed by :mod:`asb.features` / :mod:`asb.evaluation`.

All routines are pure and side-effect free. The only stochastic routine,
:func:`percolation_threshold`, takes an explicit ``rng`` so results are fully
reproducible; there is no global random state and no wall-clock seeding.

Notes
-----
``lambda2_alone`` is the algebraic connectivity (second-smallest combinatorial
Laplacian eigenvalue). It is computed here via :mod:`asb.spectral` from a locally
built ``L = D - W`` so this module needs no other sibling module at import time.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import networkx as nx
import scipy.sparse as sp

from asb.types import AtrialGraph

__all__ = [
    "fibrosis_burden",
    "fibrosis_spatial_entropy",
    "fibrosis_patch_size",
    "min_cut_value",
    "percolation_threshold",
    "lambda2_alone",
    "baseline_feature_vector",
]

# Nodes with fibrosis at or above this fraction are treated as "fibrotic" for the
# discrete patch-size / connectivity descriptors.
_FIBROTIC_THRESHOLD = 0.5


def fibrosis_burden(G: AtrialGraph) -> float:
    """Mean per-node fibrosis fraction.

    Parameters
    ----------
    G : AtrialGraph
        Atrial conduction graph.

    Returns
    -------
    float
        Average fibrosis over all nodes, in ``[0, 1]``. Approximately ``1`` for a
        fully fibrotic graph and ``0`` for a clean one. Returns ``0.0`` for an
        empty graph.
    """
    fib = np.asarray(G.fibrosis, dtype=float)
    if fib.size == 0:
        return 0.0
    return float(np.mean(fib))


def fibrosis_spatial_entropy(G: AtrialGraph, n_bins: int = 16) -> float:
    """Shannon entropy of the fibrosis-value distribution (normalized to ``[0, 1]``).

    The per-node fibrosis fractions are binned into ``n_bins`` equal bins over
    ``[0, 1]`` and the Shannon entropy of the resulting empirical distribution is
    computed and divided by ``log(n_bins)``. A spatially uniform / single-valued
    fibrosis field has entropy near ``0``; a field spread evenly across the range
    approaches ``1``.

    Parameters
    ----------
    G : AtrialGraph
        Atrial conduction graph.
    n_bins : int, optional
        Number of histogram bins over ``[0, 1]``, by default 16.

    Returns
    -------
    float
        Normalized entropy in ``[0, 1]``. Returns ``0.0`` for an empty graph or
        ``n_bins <= 1``.
    """
    fib = np.asarray(G.fibrosis, dtype=float)
    if fib.size == 0 or n_bins <= 1:
        return 0.0
    counts, _ = np.histogram(np.clip(fib, 0.0, 1.0), bins=n_bins, range=(0.0, 1.0))
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts.astype(float) / total
    nz = p > 0
    entropy = -np.sum(p[nz] * np.log(p[nz]))
    return float(entropy / np.log(n_bins))


def _fibrotic_subgraph_components(G: AtrialGraph) -> list[int]:
    """Return sizes of connected components of the fibrotic-node subgraph.

    A node is fibrotic if its fibrosis fraction is ``>= _FIBROTIC_THRESHOLD``.
    Edges are kept only when both endpoints are fibrotic.
    """
    fib = np.asarray(G.fibrosis, dtype=float)
    fibrotic = fib >= _FIBROTIC_THRESHOLD
    fib_nodes = np.nonzero(fibrotic)[0]
    if fib_nodes.size == 0:
        return []

    H = nx.Graph()
    H.add_nodes_from(int(u) for u in fib_nodes)
    edges = np.asarray(G.edges, dtype=np.int64)
    if edges.size:
        for a, b in edges:
            if fibrotic[a] and fibrotic[b]:
                H.add_edge(int(a), int(b))
    return [len(c) for c in nx.connected_components(H)]


def fibrosis_patch_size(G: AtrialGraph) -> float:
    """Mean size of connected fibrotic patches.

    Fibrotic nodes (fibrosis ``>= 0.5``) are connected through the graph edges;
    the mean number of nodes per connected component is returned.

    Parameters
    ----------
    G : AtrialGraph
        Atrial conduction graph.

    Returns
    -------
    float
        Mean connected-fibrotic-patch size (number of nodes). ``0.0`` when there
        are no fibrotic nodes.
    """
    sizes = _fibrotic_subgraph_components(G)
    if not sizes:
        return 0.0
    return float(np.mean(sizes))


def min_cut_value(G: AtrialGraph) -> float:
    """Deterministic global minimum cut of the weighted graph (Stoer-Wagner).

    Uses :func:`networkx.stoer_wagner` on the weighted conduction graph. If the
    graph is empty, has a single node, or is already disconnected, the minimum cut
    is ``0.0`` (no edge weight need be removed to separate it).

    Parameters
    ----------
    G : AtrialGraph
        Atrial conduction graph.

    Returns
    -------
    float
        Weight of the global minimum edge cut. Small for barbell / two-clique
        graphs joined by a thin bridge.
    """
    n = G.n_nodes
    if n < 2 or G.n_edges == 0:
        return 0.0

    H = nx.Graph()
    H.add_nodes_from(range(n))
    edges = np.asarray(G.edges, dtype=np.int64)
    weights = np.asarray(G.weights, dtype=float)
    for (a, b), w in zip(edges, weights):
        # Accumulate parallel edges into a single positive weight.
        if H.has_edge(int(a), int(b)):
            H[int(a)][int(b)]["weight"] += float(w)
        else:
            H.add_edge(int(a), int(b), weight=float(w))

    if not nx.is_connected(H):
        return 0.0

    cut_value, _ = nx.stoer_wagner(H, weight="weight")
    return float(cut_value)


def _largest_component_fraction(n: int, edges: np.ndarray) -> float:
    """Fraction of nodes in the largest connected component."""
    if n == 0:
        return 0.0
    H = nx.Graph()
    H.add_nodes_from(range(n))
    if edges.size:
        H.add_edges_from((int(a), int(b)) for a, b in edges)
    if H.number_of_nodes() == 0:
        return 0.0
    largest = max((len(c) for c in nx.connected_components(H)), default=0)
    return largest / n


def percolation_threshold(
    G: AtrialGraph, n_steps: int = 25, rng: Optional[np.random.Generator] = None
) -> float:
    """Edge-removal fraction at which the giant component breaks.

    Edges are removed in a single random order (drawn from ``rng``); the removed
    fraction is swept over ``n_steps`` levels in ``[0, 1]``. The threshold is the
    smallest removed fraction at which the largest connected component drops below
    half of the nodes.

    Parameters
    ----------
    G : AtrialGraph
        Atrial conduction graph.
    n_steps : int, optional
        Number of removal levels to sweep, by default 25.
    rng : numpy.random.Generator, optional
        Random generator controlling the edge-removal order. If ``None`` a fresh
        default generator is used; pass an explicit ``rng`` for reproducibility.

    Returns
    -------
    float
        Removed-edge fraction in ``[0, 1]`` at which the giant component (>= half
        the nodes) first breaks. Returns ``0.0`` if the graph is already
        fragmented and ``1.0`` if it never fragments within the sweep.
    """
    if rng is None:
        rng = np.random.default_rng()
    n = G.n_nodes
    m = G.n_edges
    if n < 2 or m == 0:
        return 0.0

    edges = np.asarray(G.edges, dtype=np.int64)
    perm = rng.permutation(m)
    edges = edges[perm]

    if n_steps < 1:
        n_steps = 1
    fractions = np.linspace(0.0, 1.0, n_steps + 1)
    for f in fractions:
        n_remove = int(round(f * m))
        kept = edges[n_remove:]
        if _largest_component_fraction(n, kept) < 0.5:
            return float(f)
    return 1.0


def lambda2_alone(G: AtrialGraph) -> float:
    """Algebraic connectivity: the Fiedler value of the combinatorial Laplacian.

    The combinatorial Laplacian ``L = D - W`` is built locally from the graph's
    symmetric adjacency and its second-smallest eigenvalue is returned via
    :func:`asb.spectral.fiedler`.

    Parameters
    ----------
    G : AtrialGraph
        Atrial conduction graph.

    Returns
    -------
    float
        Second-smallest combinatorial-Laplacian eigenvalue (``lambda_2``). ``0.0``
        for graphs with fewer than two nodes. Note this is ``0`` (to numerical
        tolerance) for a disconnected graph.
    """
    from asb.spectral import fiedler

    n = G.n_nodes
    if n < 2:
        return 0.0
    W = G.adjacency()
    deg = np.asarray(W.sum(axis=1)).ravel()
    L = sp.csr_matrix(sp.diags(deg) - W)
    lam2, _ = fiedler(L)
    return float(lam2)


def baseline_feature_vector(G: AtrialGraph) -> Dict[str, float]:
    """Assemble all baseline features into a single named dictionary.

    Parameters
    ----------
    G : AtrialGraph
        Atrial conduction graph.

    Returns
    -------
    dict of str -> float
        Named baseline features:

        - ``"fibrosis_burden"``
        - ``"fibrosis_spatial_entropy"``
        - ``"fibrosis_patch_size"``
        - ``"min_cut_value"``
        - ``"percolation_threshold"``
        - ``"lambda2_alone"``

    Notes
    -----
    ``percolation_threshold`` here uses a fixed internal generator
    (``default_rng(0)``) so the assembled vector is deterministic; call
    :func:`percolation_threshold` directly with your own ``rng`` when a different
    stream is desired.
    """
    return {
        "fibrosis_burden": fibrosis_burden(G),
        "fibrosis_spatial_entropy": fibrosis_spatial_entropy(G),
        "fibrosis_patch_size": fibrosis_patch_size(G),
        "min_cut_value": min_cut_value(G),
        "percolation_threshold": percolation_threshold(G, rng=np.random.default_rng(0)),
        "lambda2_alone": lambda2_alone(G),
    }
