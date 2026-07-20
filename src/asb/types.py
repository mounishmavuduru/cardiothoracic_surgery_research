"""Canonical shared data types for AtrialSpectralBench.

FROZEN CONTRACT: every module imports these types; do not change field names or
shapes without updating BUILD_SPEC.md and every consumer.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import scipy.sparse as sp


@dataclass
class AtrialGraph:
    """Weighted conduction graph of one atrium.

    coords       : (n, 3) float   node xyz positions.
    edges        : (m, 2) int     undirected endpoints, stored with i < j.
    weights      : (m,)   float    positive edge conductance w_ij.
    fibrosis     : (n,)   float    per-node fibrosis fraction in [0, 1].
    uac          : (n, 2) float    universal atrial coords (alpha, beta) in [0, 1].
    region       : (n,)   int      region label per node.
    shape_family : str             id of the base anatomy this graph derives from
                                   (used to group train/test folds so PCA-resampled
                                   near-duplicates never leak across the split).
    """

    coords: np.ndarray
    edges: np.ndarray
    weights: np.ndarray
    fibrosis: np.ndarray
    uac: np.ndarray
    region: np.ndarray
    shape_family: str
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.coords = np.asarray(self.coords, dtype=float)
        self.edges = np.asarray(self.edges, dtype=np.int64)
        self.weights = np.asarray(self.weights, dtype=float)
        self.fibrosis = np.asarray(self.fibrosis, dtype=float)
        self.uac = np.asarray(self.uac, dtype=float)
        self.region = np.asarray(self.region, dtype=np.int64)

    @property
    def n_nodes(self) -> int:
        return int(self.coords.shape[0])

    @property
    def n_edges(self) -> int:
        return int(self.edges.shape[0])

    def adjacency(self) -> sp.csr_matrix:
        """Symmetric weighted adjacency W (csr, shape n x n)."""
        n = self.n_nodes
        i, j = self.edges[:, 0], self.edges[:, 1]
        rows = np.concatenate([i, j])
        cols = np.concatenate([j, i])
        data = np.concatenate([self.weights, self.weights])
        return sp.csr_matrix((data, (rows, cols)), shape=(n, n))

    def degree(self) -> np.ndarray:
        """Weighted degree vector d_i = sum_j w_ij, shape (n,)."""
        return np.asarray(self.adjacency().sum(axis=1)).ravel()


@dataclass
class AtrialMesh:
    """Triangulated atrial surface with physiological fields.

    points       : (n, 3) float   vertex positions.
    faces        : (f, 3) int     triangle vertex indices.
    fibres       : (n, 3) float   unit fibre direction per vertex.
    uac          : (n, 2) float   universal atrial coords (alpha, beta).
    fibrosis     : (n,)   float    per-vertex fibrosis fraction in [0, 1].
    region       : (n,)   int      region label per vertex.
    shape_family : str             base-anatomy id (for grouped CV).
    """

    points: np.ndarray
    faces: np.ndarray
    fibres: np.ndarray
    uac: np.ndarray
    fibrosis: np.ndarray
    region: np.ndarray
    shape_family: str
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.points = np.asarray(self.points, dtype=float)
        self.faces = np.asarray(self.faces, dtype=np.int64)
        self.fibres = np.asarray(self.fibres, dtype=float)
        self.uac = np.asarray(self.uac, dtype=float)
        self.fibrosis = np.asarray(self.fibrosis, dtype=float)
        self.region = np.asarray(self.region, dtype=np.int64)

    @property
    def n_points(self) -> int:
        return int(self.points.shape[0])


@dataclass
class InducibilityLabel:
    """Ground-truth EP verdict for one (graph, protocol).

    inducible      : bool    whether a self-sustaining reentry was induced.
    reentry_origin : int|None node index of earliest reentry site (None if not inducible).
    protocol       : str     pacing protocol name, e.g. 'S1S2' or 'burst'.
    source         : str     'mock_ep' (in-repo development stand-in) or 'opencarp'
                             (real ground truth). NEVER treat 'mock_ep' as real POAF.
    """

    inducible: bool
    reentry_origin: Optional[int]
    protocol: str
    source: str
    meta: dict = field(default_factory=dict)
