r"""Real Roney LA virtual-cohort loader + field-preserving coarsening.

Loads the public Roney *et al.* left-atrial virtual-cohort surface meshes
(Zenodo record **5801337**, "Predicting atrial fibrillation recurrence by
combining population data and virtual cohorts of patient-specific left atrial
models") into :class:`asb.types.AtrialMesh` objects with **real** physiological
fields, and coarsens them to a resolution the CPU EP/spectral pipeline can chew
on while carrying every field through.

Each archived ``Mesh_<id>.vtk`` is an ASCII legacy-VTK ``POLYDATA`` surface with:

- ``POINTS``                    : vertex xyz in **microns**;
- point ``SCALARS UAC1``/``UAC2``: universal atrial coordinates :math:`(\alpha,\beta)\in[0,1]`;
- point ``SCALARS IIR``          : LGE image-intensity ratio (fibrosis proxy);
- cell  ``VECTORS fiber_endo``/``fiber_epi``: per-triangle fibre directions.

Fibrosis is derived from the IIR field by the frozen, literature-anchored ramp in
:func:`iir_to_fibrosis` (see ``docs/PRE_REGISTRATION.md``). Fibres are mapped from
per-cell to per-vertex (endocardial layer). Regions are a coarse UAC grid so the
per-region SFI aggregation is well defined.

This module *does* read files, so it lives in the substrate/loader layer (never
imported by the pure compute modules). It is deterministic in its inputs; the only
randomness is the tie-break seed in :func:`coarsen_mesh`.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

from asb.types import AtrialMesh

__all__ = [
    "read_vtk_polydata_full",
    "iir_to_fibrosis",
    "uac_regions",
    "load_roney_mesh",
    "coarsen_mesh",
    # Frozen fibrosis-calibration constants (pre-registered).
    "IIR_HEALTHY",
    "IIR_DENSE",
    "MICRON_TO_MM",
]

#: IIR at/below which tissue is treated as healthy (fibrosis = 0). Khurram et al.
#: 2014 (Heart Rhythm) put normal myocardium at IIR ~= 1.0.
IIR_HEALTHY: float = 1.0
#: IIR at/above which tissue is treated as dense fibrosis (fibrosis = 1). The
#: IIR = 1.32 threshold is the widely used dense-scar cut (Khurram 2014; used in
#: the Roney LA-fibrosis workflow).
IIR_DENSE: float = 1.32
#: Roney meshes store coordinates in microns; convert to millimetres.
MICRON_TO_MM: float = 1.0e-3


# --------------------------------------------------------------------------- #
# VTK parsing
# --------------------------------------------------------------------------- #
def read_vtk_polydata_full(
    path: str,
) -> Dict[str, object]:
    """Parse points, triangle faces, point scalars and cell vectors from a VTK file.

    Supports the ASCII legacy ``POLYDATA`` subset used by the Roney meshes:
    ``POINTS`` + ``POLYGONS`` (triangles) + ``POINT_DATA`` ``SCALARS`` +
    ``CELL_DATA`` ``VECTORS``. Non-triangular polygons are fan-triangulated.

    Parameters
    ----------
    path : str
        Path to an ASCII legacy-VTK PolyData file.

    Returns
    -------
    dict
        ``points`` (n, 3) float, ``faces`` (f, 3) int, ``point_scalars``
        {name: (n,) float}, ``cell_vectors`` {name: (f, 3) float}.
    """
    with open(path) as fh:
        lines = fh.read().splitlines()

    points = np.zeros((0, 3), dtype=float)
    faces_list: list[list[int]] = []
    point_scalars: Dict[str, np.ndarray] = {}
    cell_vectors: Dict[str, np.ndarray] = {}

    n_points = 0
    n_cells = 0
    i = 0
    n_lines = len(lines)
    while i < n_lines:
        parts = lines[i].split()
        if not parts:
            i += 1
            continue
        key = parts[0].upper()

        if key == "POINTS":
            n_points = int(parts[1])
            vals: list[float] = []
            i += 1
            while len(vals) < 3 * n_points and i < n_lines:
                vals.extend(float(x) for x in lines[i].split())
                i += 1
            points = np.asarray(vals[: 3 * n_points], dtype=float).reshape(n_points, 3)
            continue

        if key == "POLYGONS":
            n_cells = int(parts[1])
            i += 1
            read = 0
            while read < n_cells and i < n_lines:
                cell = [int(x) for x in lines[i].split()]
                i += 1
                if not cell:
                    continue
                count, idx = cell[0], cell[1:]
                for k in range(1, count - 1):
                    faces_list.append([idx[0], idx[k], idx[k + 1]])
                read += 1
            continue

        if key == "SCALARS":
            name = parts[1]
            i += 1
            # Optional LOOKUP_TABLE line.
            if i < n_lines and lines[i].split() and lines[i].split()[0].upper() == "LOOKUP_TABLE":
                i += 1
            vals = []
            while len(vals) < n_points and i < n_lines:
                row = lines[i].split()
                if row and row[0].upper() in _SECTION_KEYWORDS:
                    break
                vals.extend(float(x) for x in row)
                i += 1
            point_scalars[name] = np.asarray(vals[:n_points], dtype=float)
            continue

        if key == "VECTORS":
            name = parts[1]
            i += 1
            vals = []
            while len(vals) < 3 * n_cells and i < n_lines:
                row = lines[i].split()
                if row and row[0].upper() in _SECTION_KEYWORDS:
                    break
                vals.extend(float(x) for x in row)
                i += 1
            cell_vectors[name] = np.asarray(vals[: 3 * n_cells], dtype=float).reshape(-1, 3)
            continue

        i += 1

    return {
        "points": points,
        "faces": np.asarray(faces_list, dtype=np.int64),
        "point_scalars": point_scalars,
        "cell_vectors": cell_vectors,
    }


_SECTION_KEYWORDS = {
    "POINTS", "POLYGONS", "VERTICES", "LINES", "TRIANGLE_STRIPS",
    "POINT_DATA", "CELL_DATA", "SCALARS", "VECTORS", "NORMALS", "TENSORS",
    "FIELD", "LOOKUP_TABLE", "COLOR_SCALARS",
}


# --------------------------------------------------------------------------- #
# Field derivation
# --------------------------------------------------------------------------- #
def iir_to_fibrosis(
    iir: np.ndarray, *, healthy: float = IIR_HEALTHY, dense: float = IIR_DENSE
) -> np.ndarray:
    r"""Map the LGE image-intensity ratio to a continuous fibrosis fraction.

    A linear ramp clipped to :math:`[0, 1]`:
    :math:`f = \mathrm{clip}((\mathrm{IIR} - \text{healthy}) / (\text{dense} -
    \text{healthy}), 0, 1)`. Tissue at ``IIR <= healthy`` is healthy (0) and at
    ``IIR >= dense`` is dense fibrosis (1); values between ramp linearly. The
    thresholds are the pre-registered, literature-anchored constants
    :data:`IIR_HEALTHY` (1.0) and :data:`IIR_DENSE` (1.32).

    Parameters
    ----------
    iir : (n,) array_like
        Per-vertex image-intensity ratio.
    healthy, dense : float, optional
        Ramp endpoints (defaults are the frozen constants).

    Returns
    -------
    (n,) ndarray
        Fibrosis fraction in ``[0, 1]``.
    """
    iir = np.asarray(iir, dtype=float)
    denom = float(dense) - float(healthy)
    if denom <= 0:
        raise ValueError("dense must exceed healthy")
    return np.clip((iir - float(healthy)) / denom, 0.0, 1.0)


def uac_regions(uac: np.ndarray, n_alpha: int = 3, n_beta: int = 3) -> np.ndarray:
    """Coarse integer region labels from a UAC grid (``n_beta`` x ``n_alpha``)."""
    uac = np.asarray(uac, dtype=float)
    alpha = np.clip(uac[:, 0], 0.0, 1.0 - 1e-12)
    beta = np.clip(uac[:, 1], 0.0, 1.0 - 1e-12)
    ab = np.clip((alpha * n_alpha).astype(np.int64), 0, n_alpha - 1)
    bb = np.clip((beta * n_beta).astype(np.int64), 0, n_beta - 1)
    return (bb * n_alpha + ab).astype(np.int64)


def _cell_vectors_to_vertex(
    n_points: int, faces: np.ndarray, cell_vec: np.ndarray
) -> np.ndarray:
    """Average per-cell vectors onto vertices and unit-normalize (fallback +x)."""
    acc = np.zeros((n_points, 3), dtype=float)
    cnt = np.zeros(n_points, dtype=float)
    for t in range(faces.shape[0]):
        v = cell_vec[t] if t < cell_vec.shape[0] else np.array([1.0, 0.0, 0.0])
        for k in range(3):
            acc[faces[t, k]] += v
            cnt[faces[t, k]] += 1.0
    nz = cnt > 0
    acc[nz] /= cnt[nz, None]
    norm = np.linalg.norm(acc, axis=1)
    good = norm > 1e-12
    out = np.zeros((n_points, 3), dtype=float)
    out[good] = acc[good] / norm[good, None]
    out[~good] = np.array([1.0, 0.0, 0.0])
    return out


def load_roney_mesh(
    path: str, *, shape_family: Optional[str] = None
) -> AtrialMesh:
    """Load one Roney LA VTK mesh into an :class:`AtrialMesh` with real fields.

    Coordinates are converted microns -> mm, fibrosis is derived from IIR by the
    frozen :func:`iir_to_fibrosis` ramp, fibres are the endocardial per-cell field
    mapped to vertices, UAC is taken verbatim, and regions are a coarse UAC grid.

    Parameters
    ----------
    path : str
        Path to a ``Mesh_<id>.vtk`` file.
    shape_family : str, optional
        Group id (defaults to the file stem, i.e. one group per patient — the
        natural leakage-free grouping for the real cohort).

    Returns
    -------
    AtrialMesh
    """
    parsed = read_vtk_polydata_full(path)
    points = np.asarray(parsed["points"], dtype=float) * MICRON_TO_MM
    faces = np.asarray(parsed["faces"], dtype=np.int64)
    ps: Dict[str, np.ndarray] = parsed["point_scalars"]  # type: ignore[assignment]
    cv: Dict[str, np.ndarray] = parsed["cell_vectors"]  # type: ignore[assignment]

    n = points.shape[0]
    uac1 = ps.get("UAC1", np.zeros(n))
    uac2 = ps.get("UAC2", np.zeros(n))
    uac = np.clip(np.stack([uac1, uac2], axis=1), 0.0, 1.0)

    iir = ps.get("IIR", np.ones(n))
    fibrosis = iir_to_fibrosis(iir)

    fib_cell = cv.get("fiber_endo", cv.get("fiber_epi", None))
    if fib_cell is None:
        fibres = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    else:
        fibres = _cell_vectors_to_vertex(n, faces, np.asarray(fib_cell, dtype=float))

    region = uac_regions(uac)
    fam = shape_family if shape_family is not None else Path(path).stem

    return AtrialMesh(
        points=points,
        faces=faces,
        fibres=fibres,
        uac=uac,
        fibrosis=fibrosis,
        region=region,
        shape_family=fam,
        meta={"source": "roney_5801337", "path": os.path.abspath(path),
              "iir_mean": float(np.mean(iir)), "n_full": int(n)},
    )


# --------------------------------------------------------------------------- #
# Field-preserving coarsening (grid vertex clustering)
# --------------------------------------------------------------------------- #
def _largest_component(n: int, edges: np.ndarray) -> np.ndarray:
    """Boolean mask of the vertices in the largest connected component."""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    if edges.shape[0] == 0:
        keep = np.zeros(n, dtype=bool)
        if n:
            keep[0] = True
        return keep
    i, j = edges[:, 0], edges[:, 1]
    A = sp.csr_matrix(
        (np.ones(len(i) * 2), (np.concatenate([i, j]), np.concatenate([j, i]))),
        shape=(n, n),
    )
    ncomp, labels = connected_components(A, directed=False)
    if ncomp == 1:
        return np.ones(n, dtype=bool)
    sizes = np.bincount(labels, minlength=ncomp)
    biggest = int(np.argmax(sizes))
    return labels == biggest


def coarsen_mesh(
    mesh: AtrialMesh, target_nodes: int, *, seed: int = 0
) -> AtrialMesh:
    """Coarsen a surface mesh to ~``target_nodes`` vertices, carrying every field.

    Uses uniform-grid **vertex clustering**: the bounding box is divided into a
    cubic grid whose cell size is chosen so the number of occupied cells is close
    to ``target_nodes``. Each occupied cell becomes one coarse vertex at the
    centroid of its members; per-vertex fields (``uac``, ``fibrosis``, ``region``)
    are averaged (region by majority vote) and fibres are the normalized mean.
    Faces are remapped to cluster indices with degenerate (collapsed) triangles
    dropped, and only the largest connected component is kept so the result is a
    single valid closed-ish surface.

    Parameters
    ----------
    mesh : AtrialMesh
        Fine input surface.
    target_nodes : int
        Desired coarse vertex count (approximate).
    seed : int
        Unused except to keep a deterministic, explicit signature (clustering is
        deterministic in the geometry).

    Returns
    -------
    AtrialMesh
        Coarsened surface with the same ``shape_family`` and derived fields.
    """
    del seed
    points = np.asarray(mesh.points, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    n = points.shape[0]
    if n <= target_nodes or n == 0:
        return mesh

    lo = points.min(axis=0)
    hi = points.max(axis=0)
    extent = np.maximum(hi - lo, 1e-9)
    # Choose grid resolution R so R^3 * fill ~ target; start from cube-root and
    # adjust by the surface fill fraction empirically via a short search.
    def occupied_count(res: int) -> Tuple[int, np.ndarray]:
        cell = np.floor((points - lo) / extent * res).astype(np.int64)
        cell = np.clip(cell, 0, res - 1)
        key = (cell[:, 0].astype(np.int64) * res + cell[:, 1]) * res + cell[:, 2]
        uniq, inv = np.unique(key, return_inverse=True)
        return uniq.shape[0], inv

    # Binary-search the grid resolution for ~target_nodes occupied cells.
    lo_r, hi_r = 2, 512
    inv = None
    best_inv = None
    best_res = lo_r
    for _ in range(24):
        mid = (lo_r + hi_r) // 2
        cnt, inv = occupied_count(mid)
        best_inv, best_res = inv, mid
        if cnt < target_nodes:
            lo_r = mid + 1
        elif cnt > int(target_nodes * 1.15):
            hi_r = mid - 1
        else:
            break
        if lo_r > hi_r:
            break
    inv = best_inv
    n_clusters = int(inv.max()) + 1

    # Aggregate members into clusters.
    coords = np.zeros((n_clusters, 3), dtype=float)
    fib = np.zeros(n_clusters, dtype=float)
    uac = np.zeros((n_clusters, 2), dtype=float)
    fibre = np.zeros((n_clusters, 3), dtype=float)
    cnt = np.zeros(n_clusters, dtype=float)
    np.add.at(coords, inv, points)
    np.add.at(fib, inv, np.asarray(mesh.fibrosis, dtype=float))
    np.add.at(uac, inv, np.asarray(mesh.uac, dtype=float))
    np.add.at(fibre, inv, np.asarray(mesh.fibres, dtype=float))
    np.add.at(cnt, inv, 1.0)
    nz = cnt > 0
    coords[nz] /= cnt[nz, None]
    fib[nz] /= cnt[nz]
    uac[nz] /= cnt[nz, None]
    fnorm = np.linalg.norm(fibre, axis=1)
    good = fnorm > 1e-12
    fibre[good] /= fnorm[good, None]
    fibre[~good] = np.array([1.0, 0.0, 0.0])
    region = uac_regions(uac)

    # Remap faces; drop degenerate.
    new_faces = inv[faces]
    keep = (
        (new_faces[:, 0] != new_faces[:, 1])
        & (new_faces[:, 1] != new_faces[:, 2])
        & (new_faces[:, 0] != new_faces[:, 2])
    )
    new_faces = new_faces[keep]
    # Deduplicate faces.
    new_faces = np.unique(np.sort(new_faces, axis=1), axis=0)

    # Keep the largest connected component so the surface is a single piece.
    from asb.substrate.mesh import mesh_edges

    edges = mesh_edges(new_faces) if new_faces.shape[0] else np.zeros((0, 2), np.int64)
    comp = _largest_component(n_clusters, edges)
    remap = -np.ones(n_clusters, dtype=np.int64)
    remap[comp] = np.arange(int(comp.sum()))
    face_in = comp[new_faces[:, 0]] & comp[new_faces[:, 1]] & comp[new_faces[:, 2]] \
        if new_faces.shape[0] else np.zeros(0, dtype=bool)
    final_faces = remap[new_faces[face_in]] if new_faces.shape[0] else np.zeros((0, 3), np.int64)

    return AtrialMesh(
        points=coords[comp],
        faces=final_faces,
        fibres=fibre[comp],
        uac=uac[comp],
        fibrosis=fib[comp],
        region=region[comp],
        shape_family=mesh.shape_family,
        meta={**dict(mesh.meta), "coarsened_to": int(comp.sum()),
              "coarsen_grid_res": int(best_res), "n_full": int(n)},
    )
