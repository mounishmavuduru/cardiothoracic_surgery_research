"""Mesh <-> weighted-graph conversion for AtrialSpectralBench.

A triangulated :class:`asb.types.AtrialMesh` is turned into a weighted
conduction :class:`asb.types.AtrialGraph` whose edges are the unique triangle
edges and whose weights come from the fibre-anisotropic, fibrosis-attenuated
conductance rule in :func:`asb.graph.edge_weights_from_fibres`.

All functions are pure and deterministic; no I/O, no global RNG state.
"""
from __future__ import annotations

import numpy as np

from asb.types import AtrialGraph, AtrialMesh

__all__ = ["mesh_edges", "geodesic_edge_length", "mesh_to_graph"]


def mesh_edges(faces: np.ndarray) -> np.ndarray:
    """Unique undirected edges of a triangle mesh, stored with ``i < j``.

    Parameters
    ----------
    faces : array_like, shape (f, 3)
        Triangle vertex indices.

    Returns
    -------
    np.ndarray, shape (m, 2) int
        Unique undirected endpoint pairs, each row sorted so column 0 < column 1,
        and the rows lexicographically sorted for determinism.
    """
    faces = np.asarray(faces, dtype=np.int64)
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError(f"faces must have shape (f, 3), got {faces.shape}")

    # The three undirected edges of every triangle.
    e = np.concatenate(
        [faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0
    )
    # Canonical orientation i < j.
    e = np.sort(e, axis=1)
    # Unique rows, lexicographically ordered.
    e = np.unique(e, axis=0)
    return e.astype(np.int64)


def geodesic_edge_length(points: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Edge lengths (Euclidean surface-geodesic proxy) for mesh edges.

    On a fine triangulation the straight-line length of a triangle edge is a
    first-order approximation of the surface geodesic distance between its
    endpoints, so it is used directly here.

    Parameters
    ----------
    points : array_like, shape (n, 3)
        Vertex positions.
    edges : array_like, shape (m, 2)
        Undirected endpoint index pairs.

    Returns
    -------
    np.ndarray, shape (m,)
        Non-negative edge lengths.
    """
    points = np.asarray(points, dtype=float)
    edges = np.asarray(edges, dtype=np.int64)
    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError(f"edges must have shape (m, 2), got {edges.shape}")
    diff = points[edges[:, 1]] - points[edges[:, 0]]
    return np.linalg.norm(diff, axis=1)


def mesh_to_graph(
    mesh: AtrialMesh, *, along: float = 1.0, cross: float = 0.3
) -> AtrialGraph:
    """Build a weighted conduction :class:`AtrialGraph` from a surface mesh.

    Edges are the unique triangle edges. Weights are the fibre-anisotropic,
    fibrosis-attenuated conductances from
    :func:`asb.graph.edge_weights_from_fibres`. The per-vertex fields
    (``fibrosis``, ``uac``, ``region``) and ``shape_family`` are carried through
    unchanged.

    Parameters
    ----------
    mesh : AtrialMesh
        Triangulated atrial surface with ``points``, ``faces``, ``fibres``,
        ``uac``, ``fibrosis``, ``region`` and ``shape_family``.
    along, cross : float, optional
        Along-fibre and cross-fibre base conductances; see
        :func:`asb.graph.edge_weights_from_fibres`.

    Returns
    -------
    AtrialGraph
        Graph with ``coords = mesh.points`` and edges sorted ``i < j``.
    """
    # Imported lazily so importing this module does not require asb.graph.
    from asb.graph import edge_weights_from_fibres

    edges = mesh_edges(mesh.faces)
    weights = edge_weights_from_fibres(mesh, edges, along=along, cross=cross)

    return AtrialGraph(
        coords=np.asarray(mesh.points, dtype=float),
        edges=edges,
        weights=np.asarray(weights, dtype=float),
        fibrosis=np.asarray(mesh.fibrosis, dtype=float),
        uac=np.asarray(mesh.uac, dtype=float),
        region=np.asarray(mesh.region, dtype=np.int64),
        shape_family=mesh.shape_family,
        meta=dict(mesh.meta),
    )
