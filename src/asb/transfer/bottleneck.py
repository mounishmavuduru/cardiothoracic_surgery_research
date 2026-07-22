r"""Source--sink / connectivity-bottleneck substrate generator (novelty experiment 1).

Builds excitable networks whose vulnerable site is a **connectivity bottleneck** — a
narrow isthmus of well-coupled tissue joining two dense regions — that is deliberately
*decoupled* from the fibrosis field. This isolates the clinically real regime where
reentry originates at a source--sink mismatch (a protected isthmus of surviving tissue
between scar), which LGE-MRI fibrosis density is structurally blind to.

The isthmus is **not** a low-conductance lesion: its edges carry normal weight. It is
sparse in the *graph-cut* sense (few parallel paths), so the Fiedler vector `φ₂` separates
the two communities across it and `|∇φ₂|` peaks there — while the fibrosis field, placed a
controllable distance away, points elsewhere. The instability *origin is decided by the
simulator*, not assumed; `meta` records the isthmus nodes and the fibrosis peak so the
dissociation (origin-vs-bottleneck vs origin-vs-fibrosis) can be measured honestly.

The key knob is ``lesion_offset ∈ [0,1]``: 0 puts the fibrosis patch on the isthmus
(confounded, like baseline GM4), 1 puts it deep in a community interior (fully dissociated).
Sweeping it turns "fibrosis also localizes" into a controlled test.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from asb.types import AtrialGraph

__all__ = ["BottleneckConfig", "make_bottleneck_network"]


@dataclass
class BottleneckConfig:
    """Two-community + narrow-isthmus excitable substrate parameters."""

    n_per_blob: int = 220
    blob_std: float = 0.085        # community spatial spread
    radius: float = 0.085          # within-community connection radius
    n_bridge: int = 5              # isthmus chain length (narrow => strong bottleneck)
    lesion_burden: float = 0.30    # mean fibrosis intensity
    lesion_scale: float = 0.10     # fibrosis patch spatial scale
    lesion_offset: float = 1.0     # 0 = fibrosis ON isthmus; 1 = deep in a community
    lesion_floor: float = 0.05     # residual coupling inside full lesion
    n_region_grid: int = 3


_C0 = np.array([0.30, 0.5])
_C1 = np.array([0.70, 0.5])
_ISTHMUS_MID = np.array([0.5, 0.5])


def _largest_component_mask(n: int, edges: np.ndarray) -> np.ndarray:
    if edges.shape[0] == 0:
        m = np.zeros(n, bool); m[: min(1, n)] = True; return m
    i, j = edges[:, 0], edges[:, 1]
    A = sp.csr_matrix((np.ones(2 * len(i)),
                       (np.concatenate([i, j]), np.concatenate([j, i]))), shape=(n, n))
    ncomp, lab = connected_components(A, directed=False)
    if ncomp == 1:
        return np.ones(n, bool)
    return lab == int(np.argmax(np.bincount(lab, minlength=ncomp)))


def _blob_edges(pts: np.ndarray, idx: np.ndarray, radius: float) -> List[np.ndarray]:
    """Radius-graph edges within one community (returned in global indices)."""
    if idx.size < 2:
        return []
    tree = cKDTree(pts[idx])
    pairs = tree.query_pairs(r=radius, output_type="ndarray")
    if pairs.shape[0] == 0:
        _, nn = tree.query(pts[idx], k=min(6, idx.size))
        pairs = np.array([[i, j] for i, row in enumerate(nn[:, 1:]) for j in row], np.int64)
    return [idx[pairs]]


def make_bottleneck_network(seed: int, cfg: BottleneckConfig) -> AtrialGraph:
    """One two-community network joined by a narrow isthmus, fibrosis offset by knob."""
    rng = np.random.default_rng(seed)
    n0 = cfg.n_per_blob
    p0 = _C0 + cfg.blob_std * rng.standard_normal((n0, 2))
    p1 = _C1 + cfg.blob_std * rng.standard_normal((n0, 2))
    # isthmus chain along the midline between the communities
    bx = np.linspace(0.40, 0.60, cfg.n_bridge)
    pb = np.column_stack([bx, np.full(cfg.n_bridge, 0.5)])
    pts = np.clip(np.vstack([p0, p1, pb]), 0.02, 0.98)

    i0 = np.arange(0, n0)
    i1 = np.arange(n0, 2 * n0)
    ib = np.arange(2 * n0, 2 * n0 + cfg.n_bridge)

    edge_blocks = _blob_edges(pts, i0, cfg.radius) + _blob_edges(pts, i1, cfg.radius)
    # isthmus: chain the bridge nodes, attach ends to the nearest node in each community
    chain = np.column_stack([ib[:-1], ib[1:]])
    end0 = np.array([[ib[0], i0[np.argmin(np.sum((pts[i0] - pts[ib[0]]) ** 2, axis=1))]]])
    end1 = np.array([[ib[-1], i1[np.argmin(np.sum((pts[i1] - pts[ib[-1]]) ** 2, axis=1))]]])
    edge_blocks += [chain, end0, end1]
    edges = np.unique(np.sort(np.vstack(edge_blocks), axis=1), axis=0).astype(np.int64)

    n_all = pts.shape[0]
    keep = _largest_component_mask(n_all, edges)
    remap = -np.ones(n_all, np.int64)
    remap[keep] = np.arange(int(keep.sum()))
    emask = keep[edges[:, 0]] & keep[edges[:, 1]]
    edges = remap[edges[emask]]
    pts = pts[keep]
    n = pts.shape[0]
    isthmus_nodes = remap[ib[keep[ib]]]

    # fibrosis patch: centered between the isthmus (offset 0) and community-0 core (offset 1)
    center = _ISTHMUS_MID + cfg.lesion_offset * (_C0 - _ISTHMUS_MID)
    d2 = np.sum((pts - center) ** 2, axis=1)
    field = np.exp(-d2 / (2 * cfg.lesion_scale ** 2))
    m = field.mean()
    if m > 0:
        field *= cfg.lesion_burden / m
    fibrosis = np.clip(field, 0.0, 1.0)

    # coupling: base, reduced by endpoint-mean fibrosis (isthmus edges keep ~full weight,
    # since the isthmus is deliberately low-fibrosis when offset > 0)
    les_edge = 0.5 * (fibrosis[edges[:, 0]] + fibrosis[edges[:, 1]])
    weights = np.clip(1.0 - les_edge, cfg.lesion_floor, 1.0)

    coords = np.column_stack([pts, np.zeros(n)])
    uac = np.clip(pts, 0.0, 1.0)
    g = cfg.n_region_grid
    rb = np.clip((uac[:, 1] * g).astype(np.int64), 0, g - 1)
    cb = np.clip((uac[:, 0] * g).astype(np.int64), 0, g - 1)
    region = (rb * g + cb).astype(np.int64)

    return AtrialGraph(
        coords=coords, edges=edges, weights=weights, fibrosis=fibrosis,
        uac=uac, region=region, shape_family=f"btl_{seed}",
        meta={"medium": "bottleneck_network", "seed": int(seed),
              "lesion_offset": float(cfg.lesion_offset),
              "isthmus_nodes": [int(x) for x in isthmus_nodes if x >= 0],
              "fibrosis_peak": int(np.argmax(fibrosis)),
              "isthmus_xy": [0.5, 0.5], "n_nodes": int(n)},
    )
