"""Lightweight excitable-media EP inducibility surrogate (DEVELOPMENT STAND-IN).

.. warning::

    Everything in this module is a **development stand-in for openCARP**, used
    only to exercise the downstream AtrialSpectralBench pipeline until real
    electrophysiology labels exist. The ``inducible`` verdict it produces is
    **NOT clinical post-operative atrial fibrillation (POAF)** and must never be
    interpreted, reported, or published as such. The label ``source`` is always
    ``'mock_ep'`` so consumers can tell it apart from real ``'opencarp'``
    ground truth.

Physical model
--------------
Rather than integrate a full monodomain reaction-diffusion system, we use a
cheap **eikonal + wavelength / source-sink** surrogate on the conduction graph:

1. **Eikonal activation.** Conduction velocity on an edge scales with the square
   root of its (normalized) conductance, ``CV_ij = sqrt(w_ij / max_w)`` — the
   standard cable-theory ``CV ~ sqrt(D)`` relation. The traversal time of an
   edge is its geodesic length divided by ``CV``. A wavefront paced from a few
   sites gives per-node activation times via a multi-source Dijkstra solve.

2. **Wavelength.** The reentry wavelength is ``lambda = CV * ERP`` where the
   effective refractory period ``ERP`` is *shortened* by local fibrosis
   (fibrotic / remodelled tissue has a shorter ERP), so fibrosis shrinks the
   wavelength through both slowed conduction and shortened refractoriness.

3. **Source-sink mismatch.** At each node the wavefront must charge its
   downstream ("sink") neighbours from its upstream ("source") neighbours.
   Where a poorly-coupled isthmus expands into well-coupled tissue the sink
   conductance exceeds the source conductance and unidirectional block — the
   seed of reentry — becomes likely.

The per-node **reentry vulnerability** multiplies the local source-sink
mismatch by the inverse local wavelength (short wavelength => small sustaining
circuit fits). An atrium is flagged *inducible* when its peak vulnerability
exceeds, by a protocol-dependent factor, the peak vulnerability of the same
anatomy with a *healthy* (uniform-conductance, no ERP-shortening) substrate.
The reentry origin is the node attaining that peak.

Because the criterion is built from **wavelength and source-sink dispersion**
— never from the Fiedler value ``lambda_2`` — it is *not* circular with the
Spectral Fragility Index it is meant to help validate. It does, by
construction, tend to flag low-conductance / high-fibrosis, spatially
heterogeneous atria as more inducible so the downstream ML task has signal.

All functions are deterministic: the only stochasticity is pacing-site
selection, which is driven entirely by the caller-supplied ``rng``.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra

from asb.config import LabelConfig
from asb.types import AtrialGraph, InducibilityLabel

__all__ = ["induce", "paced_activation"]

# --- Surrogate physiology constants (dimensionless; documented stand-ins) ---
#: Baseline effective refractory period (arbitrary time units; only ratios matter).
_ERP0 = 1.0
#: Fraction by which unit fibrosis shortens the local ERP (AF-like remodelling).
_ERP_FIBROSIS_COEFF = 0.6
#: Numerical floor to keep conduction velocities / lengths strictly positive.
_EPS = 1e-12
#: Peak-vulnerability ratio above the healthy baseline that trips inducibility,
#: keyed by pacing protocol (burst pacing is the more aggressive inducer).
_REENTRY_THRESHOLD = {"S1S2": 1.30, "burst": 1.15}
_DEFAULT_THRESHOLD = 1.30


def _edge_geometry(G: AtrialGraph) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-edge endpoint indices, geodesic length, and normalized weight.

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph.

    Returns
    -------
    edges : np.ndarray, shape (m, 2) int
        Endpoint index pairs.
    length : np.ndarray, shape (m,) float
        Euclidean length of each edge (floored at ``_EPS``).
    w_norm : np.ndarray, shape (m,) float
        Edge conductances scaled by the maximum, in ``(0, 1]``.
    """
    edges = np.asarray(G.edges, dtype=np.int64)
    coords = np.asarray(G.coords, dtype=float)
    weights = np.asarray(G.weights, dtype=float)

    i, j = edges[:, 0], edges[:, 1]
    length = np.linalg.norm(coords[i] - coords[j], axis=1)
    length = np.maximum(length, _EPS)

    wmax = float(weights.max()) if weights.size else 1.0
    wmax = wmax if wmax > 0 else 1.0
    w_norm = np.clip(weights / wmax, _EPS, 1.0)
    return edges, length, w_norm


def _cost_matrix(
    edges: np.ndarray, cost: np.ndarray, n: int
) -> sp.csr_matrix:
    """Build a symmetric sparse edge-cost matrix for a Dijkstra solve.

    Parameters
    ----------
    edges : np.ndarray, shape (m, 2) int
        Endpoint index pairs.
    cost : np.ndarray, shape (m,) float
        Nonnegative traversal cost (time) per edge.
    n : int
        Number of nodes.

    Returns
    -------
    scipy.sparse.csr_matrix, shape (n, n)
        Symmetric cost matrix suitable for :func:`scipy.sparse.csgraph.dijkstra`.
    """
    i, j = edges[:, 0], edges[:, 1]
    rows = np.concatenate([i, j])
    cols = np.concatenate([j, i])
    data = np.concatenate([cost, cost])
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n))


def paced_activation(G: AtrialGraph, site: int, cfg: LabelConfig) -> np.ndarray:
    """Per-node eikonal activation times for a wavefront paced from ``site``.

    DEVELOPMENT STAND-IN for openCARP — see the module docstring. This is a
    cheap eikonal surrogate, not a monodomain solve, and its output is not
    clinical.

    Conduction velocity on an edge is ``sqrt(w_ij / max_w)`` and its traversal
    time is ``length / CV``; activation times are the graph-geodesic (least
    conduction-time) distances from ``site``.

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph.
    site : int
        Index of the paced (stimulus) node.
    cfg : LabelConfig
        Label configuration (unused by the eikonal helper but kept for a
        uniform pacing-site API and future protocol-specific timing).

    Returns
    -------
    np.ndarray, shape (n_nodes,) float
        Activation time of every node; ``np.inf`` for nodes unreachable from
        ``site``. The paced node has activation time ``0``.
    """
    n = G.n_nodes
    site = int(site)
    if not 0 <= site < n:
        raise IndexError(f"site {site} out of range for {n} nodes")

    edges, length, w_norm = _edge_geometry(G)
    cv = np.sqrt(w_norm)
    cost = length / np.maximum(cv, _EPS)
    C = _cost_matrix(edges, cost, n)
    times = dijkstra(C, directed=False, indices=site)
    return np.asarray(times, dtype=float).ravel()


def _pacing_sites(n: int, n_sites: int, rng: np.random.Generator) -> np.ndarray:
    """Choose distinct pacing-site node indices deterministically from ``rng``."""
    n_sites = max(1, min(int(n_sites), n))
    return rng.choice(n, size=n_sites, replace=False)


def _vulnerability_field(
    G: AtrialGraph,
    activation: np.ndarray,
    w_norm: np.ndarray,
    *,
    use_fibrosis_erp: bool,
) -> np.ndarray:
    """Per-node reentry-vulnerability score from wavelength and source-sink.

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph.
    activation : np.ndarray, shape (n,)
        Node activation times from the paced wavefront (``inf`` if unreached).
    w_norm : np.ndarray, shape (m,)
        Normalized edge conductances in ``(0, 1]`` (edge order matches
        ``G.edges``).
    use_fibrosis_erp : bool
        If ``True``, shorten the local ERP by local fibrosis (the fibrotic
        substrate); if ``False``, use the uniform baseline ERP (healthy null).

    Returns
    -------
    np.ndarray, shape (n,) float
        Nonnegative vulnerability score per node. Nodes the wavefront did not
        pass through (no upstream or no downstream neighbour) score ``0``.

    Notes
    -----
    ``vulnerability(v) = (sink / source) * 1 / (CV_local(v) * ERP(v))``. The
    first factor is the source-sink mismatch (block at conduction expansions);
    the second is the inverse local wavelength (short wavelength sustains small
    reentrant circuits).
    """
    n = G.n_nodes
    edges = np.asarray(G.edges, dtype=np.int64)
    i, j = edges[:, 0], edges[:, 1]
    fib = np.clip(np.asarray(G.fibrosis, dtype=float), 0.0, 1.0)

    # Local conduction velocity per node: mean CV over incident edges.
    cv_edge = np.sqrt(w_norm)
    cv_sum = np.zeros(n)
    deg = np.zeros(n)
    np.add.at(cv_sum, i, cv_edge)
    np.add.at(cv_sum, j, cv_edge)
    np.add.at(deg, i, 1.0)
    np.add.at(deg, j, 1.0)
    cv_local = np.where(deg > 0, cv_sum / np.maximum(deg, 1.0), _EPS)
    cv_local = np.maximum(cv_local, _EPS)

    # Effective refractory period -> local wavelength.
    if use_fibrosis_erp:
        erp = _ERP0 * (1.0 - _ERP_FIBROSIS_COEFF * fib)
    else:
        erp = np.full(n, _ERP0)
    erp = np.maximum(erp, _EPS)
    inv_wavelength = 1.0 / (cv_local * erp)

    # Source (upstream) and sink (downstream) conductance per node, using the
    # activation-time ordering to orient each edge along the wavefront.
    source = np.zeros(n)
    sink = np.zeros(n)
    ti, tj = activation[i], activation[j]
    finite = np.isfinite(ti) & np.isfinite(tj)

    # Edge oriented i -> j when i activates first.
    i_first = finite & (ti < tj)
    j_first = finite & (tj < ti)

    # i -> j: j receives from i (source of j), i drives j (sink of i).
    np.add.at(source, j[i_first], w_norm[i_first])
    np.add.at(sink, i[i_first], w_norm[i_first])
    # j -> i.
    np.add.at(source, i[j_first], w_norm[j_first])
    np.add.at(sink, j[j_first], w_norm[j_first])

    active = np.isfinite(activation)
    passed = active & (source > 0) & (sink > 0)

    mismatch = np.zeros(n)
    mismatch[passed] = sink[passed] / (source[passed] + _EPS)

    vuln = mismatch * inv_wavelength
    vuln[~passed] = 0.0
    return vuln


def induce(
    G: AtrialGraph, cfg: LabelConfig, rng: np.random.Generator
) -> InducibilityLabel:
    """Mock-EP inducibility verdict for one conduction graph.

    DEVELOPMENT STAND-IN for openCARP. The returned ``inducible`` flag is **not
    clinical POAF** and must never be treated as real ground truth; it exists
    only to give the downstream AtrialSpectralBench pipeline a learnable signal
    until real electrophysiology labels are available. The returned label has
    ``source='mock_ep'``.

    An eikonal wavefront is paced from a few ``rng``-chosen sites; the atrium is
    flagged inducible when its peak wavelength / source-sink reentry
    vulnerability exceeds, by a protocol-dependent factor, that of the same
    anatomy on a healthy (uniform-conductance, unshortened-ERP) substrate. The
    reentry origin is the node attaining that peak. The criterion is built from
    conduction wavelength and source-sink dispersion, deliberately **not** from
    the Fiedler value ``lambda_2``, so it is not circular with the Spectral
    Fragility Index.

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph of one atrium.
    cfg : LabelConfig
        Pacing configuration; uses ``protocol`` (``'S1S2'`` or ``'burst'``) and
        ``n_pacing_sites``.
    rng : numpy.random.Generator
        Source of randomness for pacing-site selection. All other computation
        is deterministic, so a fixed ``rng`` gives a fixed label.

    Returns
    -------
    InducibilityLabel
        ``inducible`` (bool), ``reentry_origin`` (int node index if inducible
        else ``None``), ``protocol`` (echoing ``cfg.protocol``), and
        ``source='mock_ep'``.
    """
    n = G.n_nodes
    protocol = cfg.protocol
    threshold = _REENTRY_THRESHOLD.get(protocol, _DEFAULT_THRESHOLD)

    meta: dict = {
        "note": "mock_ep development stand-in for openCARP; NOT clinical POAF",
        "protocol": protocol,
    }

    if n == 0 or G.n_edges == 0:
        return InducibilityLabel(
            inducible=False,
            reentry_origin=None,
            protocol=protocol,
            source="mock_ep",
            meta={**meta, "reason": "empty or edgeless graph"},
        )

    edges, length, w_norm = _edge_geometry(G)
    cv = np.sqrt(w_norm)
    cost = length / np.maximum(cv, _EPS)
    C = _cost_matrix(edges, cost, n)

    sites = _pacing_sites(n, cfg.n_pacing_sites, rng)
    # Combined wavefront: earliest activation over all pacing sites.
    activation = dijkstra(C, directed=False, indices=sites, min_only=True)
    activation = np.asarray(activation, dtype=float).ravel()

    # Fibrotic-substrate vulnerability vs. a healthy (uniform) null on the same
    # anatomy and same wavefront ordering.
    vuln = _vulnerability_field(G, activation, w_norm, use_fibrosis_erp=True)
    healthy = _vulnerability_field(
        G, activation, np.ones_like(w_norm), use_fibrosis_erp=False
    )

    peak = float(vuln.max()) if vuln.size else 0.0
    healthy_peak = float(healthy.max()) if healthy.size else 0.0
    reentry_index = peak / (healthy_peak + _EPS)

    inducible = bool(peak > 0.0 and reentry_index >= threshold)
    origin = int(np.argmax(vuln)) if inducible else None

    meta.update(
        {
            "reentry_index": reentry_index,
            "peak_vulnerability": peak,
            "healthy_peak_vulnerability": healthy_peak,
            "threshold": threshold,
            "n_pacing_sites": int(sites.size),
        }
    )
    return InducibilityLabel(
        inducible=inducible,
        reentry_origin=origin,
        protocol=protocol,
        source="mock_ep",
        meta=meta,
    )
