r"""Second excitable medium — a 2-D neural-network substrate cohort.

Generates spatially-embedded random-geometric excitable networks (a stand-in for a
sheet of coupled neural populations / an epilepsy-like network), each carrying a
heterogeneous **low-coupling "lesion"** field — the structural analogue of cardiac
fibrosis. Every network is returned as an :class:`asb.types.AtrialGraph` so the
*entire* spectral / SFI / baseline stack runs on it **unchanged**: that literal code
reuse is the whole point of GM4 (the fragility calculus is graph-universal, not
cardiac-specific).

Field mapping (heart → neural network):

- ``coords``   : node (x, y, 0) position in a unit square (mm-scale arbitrary units);
- ``uac``      : the (x, y) position in [0,1]² — the 2-D coordinate used for the
  localization spatial-null (rotational/shift), exactly as UAC was for the atrium;
- ``fibrosis`` : the **lesion** field in [0,1] — local coupling *reduction* (a
  low-conductance patch that promotes propagation block / reentry);
- ``weights``  : inter-node coupling conductances, reduced inside the lesion;
- ``region``   : coarse spatial grid label (for the per-region SFI aggregation);
- ``shape_family`` : a per-network group id (each network its own group — the
  leakage-free grouping, as each patient was).

Deterministic in the seed; no I/O.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from asb.types import AtrialGraph

__all__ = ["NetworkConfig", "make_excitable_network", "make_network_cohort"]


@dataclass
class NetworkConfig:
    """Random-geometric excitable-network generation parameters."""

    n_nodes: int = 900
    radius: float = 0.075       # connection radius in the unit square
    along: float = 1.0          # base coupling conductance (healthy)
    lesion_floor: float = 0.05  # residual coupling inside a full lesion
    lesion_burden: float = 0.2  # target mean lesion intensity (fibrosis analogue)
    n_lesion_patches: int = 4
    patch_scale: float = 0.12   # lesion patch spatial scale (fraction of the square)
    n_region_grid: int = 3      # region = n_region_grid x n_region_grid spatial bins


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


def _paint_lesion(pts: np.ndarray, cfg: NetworkConfig, rng: np.random.Generator) -> np.ndarray:
    """Patchy low-coupling lesion field in [0,1], mean ~= lesion_burden."""
    n = pts.shape[0]
    field = np.zeros(n)
    if cfg.n_lesion_patches > 0:
        centers = rng.choice(n, size=min(cfg.n_lesion_patches, n), replace=False)
        sigma = max(cfg.patch_scale, 1e-6)
        for c in centers:
            d2 = np.sum((pts - pts[c]) ** 2, axis=1)
            field += rng.uniform(0.5, 1.0) * np.exp(-d2 / (2 * sigma ** 2))
    fmax = field.max()
    if fmax > 0:
        field /= fmax
        mean = field.mean()
        if mean > 0:
            field *= cfg.lesion_burden / mean
    return np.clip(field, 0.0, 1.0)


def make_excitable_network(
    seed: int, cfg: NetworkConfig, *, lesion_burden: float | None = None
) -> AtrialGraph:
    """One 2-D random-geometric excitable network with a heterogeneous lesion field.

    Parameters
    ----------
    seed : int
        Determines node placement + lesion patches.
    cfg : NetworkConfig
        Generation parameters.
    lesion_burden : float, optional
        Override the mean lesion intensity for this network (used to sweep a
        stable→unstable substrate gradient across the cohort).

    Returns
    -------
    AtrialGraph
        A weighted excitable-network graph with the neural-network field mapping.
    """
    rng = np.random.default_rng(seed)
    if lesion_burden is not None:
        cfg = NetworkConfig(**{**cfg.__dict__, "lesion_burden": float(lesion_burden)})

    pts2 = rng.uniform(0.0, 1.0, size=(cfg.n_nodes, 2))
    tree = cKDTree(pts2)
    pairs = tree.query_pairs(r=cfg.radius, output_type="ndarray")
    if pairs.shape[0] == 0:
        # Fall back to a k-NN graph if the radius was too small.
        _, idx = tree.query(pts2, k=6)
        pairs = np.array([[i, j] for i, row in enumerate(idx[:, 1:]) for j in row], np.int64)
    edges = np.unique(np.sort(pairs, axis=1), axis=0).astype(np.int64)

    # Keep the largest connected component so the medium is a single sheet.
    keep = _largest_component_mask(cfg.n_nodes, edges)
    remap = -np.ones(cfg.n_nodes, np.int64)
    remap[keep] = np.arange(int(keep.sum()))
    emask = keep[edges[:, 0]] & keep[edges[:, 1]]
    edges = remap[edges[emask]]
    pts2 = pts2[keep]
    n = pts2.shape[0]

    lesion = _paint_lesion(pts2, cfg, rng)

    # Coupling conductance per edge: base, reduced by endpoint-mean lesion (floor).
    les_edge = 0.5 * (lesion[edges[:, 0]] + lesion[edges[:, 1]])
    weights = cfg.along * np.clip(1.0 - les_edge, cfg.lesion_floor, 1.0)

    coords = np.column_stack([pts2, np.zeros(n)])
    uac = np.clip(pts2, 0.0, 1.0)
    g = cfg.n_region_grid
    rb = np.clip((uac[:, 1] * g).astype(np.int64), 0, g - 1)
    cb = np.clip((uac[:, 0] * g).astype(np.int64), 0, g - 1)
    region = (rb * g + cb).astype(np.int64)

    return AtrialGraph(
        coords=coords, edges=edges, weights=weights,
        fibrosis=lesion, uac=uac, region=region,
        shape_family=f"net_{seed}",
        meta={"medium": "fhn_network", "seed": int(seed),
              "lesion_burden": float(cfg.lesion_burden), "n_nodes": int(n)},
    )


def make_network_cohort(
    n_networks: int, cfg: NetworkConfig, *, seed: int = 0,
    burden_lo: float = 0.05, burden_hi: float = 0.45,
) -> List[AtrialGraph]:
    """A cohort of excitable networks spanning a mild→severe lesion-burden gradient.

    The lesion burden is swept across ``[burden_lo, burden_hi]`` so the FHN
    instability labeller produces both stable and unstable networks (both classes),
    mirroring the atrial mild→severe fibrosis gradient.
    """
    rng = np.random.default_rng(seed + 20259)
    out: List[AtrialGraph] = []
    for k in range(n_networks):
        burden = float(rng.uniform(burden_lo, burden_hi))
        out.append(make_excitable_network(seed + k, cfg, lesion_burden=burden))
    return out
