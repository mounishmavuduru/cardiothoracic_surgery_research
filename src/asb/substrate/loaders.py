"""Real open-data loaders for AtrialSpectralBench (interface only).

These loaders parse public open cardiac geometries into
:class:`asb.types.AtrialMesh` objects. **They are never executed in the build
environment** (network access to Zenodo is blocked by policy and no data is
present). When the requested data root is missing or empty, each loader raises
:class:`DataUnavailableError` naming the exact Zenodo record and the download
steps needed to obtain it.

Datasets
--------
- Rodero et al. 1000 synthetic four-chamber meshes  -> Zenodo 4506930
- Roney et al. LA virtual cohort (fibrosis/UAC/fibres) -> Zenodo 5801337
- Roney et al. Human Atrial Fibre Atlas (CARP)        -> Zenodo 3764917
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from asb.types import AtrialMesh

__all__ = [
    "DataUnavailableError",
    "load_rodero_cohort",
    "load_roney_meshes",
    "load_fibre_atlas",
]


class DataUnavailableError(RuntimeError):
    """Raised when a real dataset root is missing/empty.

    The message names the exact Zenodo record and the download steps so the
    user can obtain the data (these loaders never download anything themselves).
    """


def _download_hint(record: str, title: str) -> str:
    """Build a uniform, actionable download instruction string."""
    return (
        f"Dataset '{title}' is not available at the given root.\n"
        f"It is NOT downloaded automatically (network is blocked at build time).\n"
        f"To obtain it:\n"
        f"  1. Visit https://zenodo.org/record/{record}\n"
        f"  2. Download the archive and extract it to a local directory ROOT.\n"
        f"  3. Re-run passing that ROOT as the loader's `root` argument.\n"
        f"Zenodo record: {record}."
    )


def _has_data(root: str) -> bool:
    """True iff ``root`` exists, is a directory, and contains at least one file."""
    if not root:
        return False
    p = Path(root)
    if not p.is_dir():
        return False
    for _dirpath, _dirnames, filenames in os.walk(p):
        if filenames:
            return True
    return False


def _list_by_suffix(root: str, suffixes: tuple[str, ...]) -> list[Path]:
    """All files under ``root`` whose suffix (lowercased) is in ``suffixes``."""
    out: list[Path] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if Path(fn).suffix.lower() in suffixes:
                out.append(Path(dirpath) / fn)
    return sorted(out)


# --------------------------------------------------------------------------- #
# Minimal format parsers (used only when real files are present)
# --------------------------------------------------------------------------- #
def _read_vtk_polydata(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse points and triangle faces from an ASCII legacy VTK PolyData file.

    Only the subset needed for atrial surfaces (POINTS + POLYGONS of triangles)
    is supported; richer attributes are ignored.
    """
    tokens: list[str] = []
    with open(path) as fh:
        text = fh.read()
    lines = text.splitlines()

    points: list[list[float]] = []
    faces: list[list[int]] = []
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if not parts:
            i += 1
            continue
        key = parts[0].upper()
        if key == "POINTS":
            n = int(parts[1])
            vals: list[float] = []
            i += 1
            while len(vals) < 3 * n:
                vals.extend(float(x) for x in lines[i].split())
                i += 1
            arr = np.asarray(vals[: 3 * n], dtype=float).reshape(n, 3)
            points = arr.tolist()
            continue
        if key in ("POLYGONS", "TRIANGLE_STRIPS"):
            n_cells = int(parts[1])
            i += 1
            read = 0
            while read < n_cells:
                cell = [int(x) for x in lines[i].split()]
                i += 1
                if not cell:
                    continue
                count, idx = cell[0], cell[1:]
                # Fan-triangulate polygons with more than three vertices.
                for k in range(1, count - 1):
                    faces.append([idx[0], idx[k], idx[k + 1]])
                read += 1
            continue
        i += 1

    _ = tokens
    return np.asarray(points, dtype=float), np.asarray(faces, dtype=np.int64)


def _read_carp(
    pts_path: Path, elem_path: Path, lon_path: Path | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Parse CARP ``.pts`` / ``.elem`` / ``.lon`` files into arrays.

    Returns
    -------
    points : (n, 3) float
    faces  : (f, 3) int   (triangular surface elements, tag 'Tr')
    fibres : (f, 3) float or None   per-element fibre direction from ``.lon``.
    """
    with open(pts_path) as fh:
        n = int(fh.readline().split()[0])
        pts = np.array(
            [[float(x) for x in fh.readline().split()[:3]] for _ in range(n)],
            dtype=float,
        )

    faces: list[list[int]] = []
    with open(elem_path) as fh:
        n_el = int(fh.readline().split()[0])
        for _ in range(n_el):
            parts = fh.readline().split()
            if not parts:
                continue
            if parts[0] == "Tr":
                faces.append([int(parts[1]), int(parts[2]), int(parts[3])])
    faces_arr = np.asarray(faces, dtype=np.int64)

    fibres: np.ndarray | None = None
    if lon_path is not None and lon_path.exists():
        rows: list[list[float]] = []
        with open(lon_path) as fh:
            first = fh.readline().split()
            # Some .lon files start with a header count; skip if non-numeric-vec.
            if len(first) >= 3:
                try:
                    rows.append([float(x) for x in first[:3]])
                except ValueError:
                    pass
            for line in fh:
                sp = line.split()
                if len(sp) >= 3:
                    rows.append([float(x) for x in sp[:3]])
        if rows:
            fibres = np.asarray(rows, dtype=float)

    return pts, faces_arr, fibres


def _mesh_from_surface(
    points: np.ndarray,
    faces: np.ndarray,
    shape_family: str,
    *,
    vertex_fibres: np.ndarray | None = None,
) -> AtrialMesh:
    """Assemble an AtrialMesh from bare surface geometry with default fields."""
    n = points.shape[0]
    if vertex_fibres is not None and vertex_fibres.shape[0] == n:
        fibres = np.asarray(vertex_fibres, dtype=float)
    else:
        fibres = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    return AtrialMesh(
        points=points,
        faces=faces,
        fibres=fibres,
        uac=np.zeros((n, 2), dtype=float),
        fibrosis=np.zeros(n, dtype=float),
        region=np.zeros(n, dtype=np.int64),
        shape_family=shape_family,
    )


# --------------------------------------------------------------------------- #
# Public loaders
# --------------------------------------------------------------------------- #
def load_rodero_cohort(root: str) -> list[AtrialMesh]:
    """Load the Rodero et al. 1000 synthetic four-chamber meshes.

    Parameters
    ----------
    root : str
        Directory containing the extracted Zenodo 4506930 archive (VTK meshes).

    Returns
    -------
    list[AtrialMesh]

    Raises
    ------
    DataUnavailableError
        If ``root`` is missing or empty. The message names Zenodo record 4506930
        and the download steps.
    """
    if not _has_data(root):
        raise DataUnavailableError(
            _download_hint("4506930", "Rodero et al. 1000 four-chamber meshes")
        )
    meshes: list[AtrialMesh] = []
    for path in _list_by_suffix(root, (".vtk", ".vtp")):
        points, faces = _read_vtk_polydata(path)
        if points.size and faces.size:
            meshes.append(_mesh_from_surface(points, faces, shape_family=path.stem))
    if not meshes:
        raise DataUnavailableError(
            _download_hint("4506930", "Rodero et al. 1000 four-chamber meshes")
        )
    return meshes


def load_roney_meshes(root: str) -> list[AtrialMesh]:
    """Load the Roney et al. LA virtual cohort (fibrosis + UAC + fibres).

    Parameters
    ----------
    root : str
        Directory containing the extracted Zenodo 5801337 archive.

    Returns
    -------
    list[AtrialMesh]

    Raises
    ------
    DataUnavailableError
        If ``root`` is missing or empty. The message names Zenodo record 5801337
        and the download steps.
    """
    if not _has_data(root):
        raise DataUnavailableError(
            _download_hint("5801337", "Roney et al. LA virtual cohort")
        )
    meshes: list[AtrialMesh] = []
    for path in _list_by_suffix(root, (".vtk", ".vtp")):
        points, faces = _read_vtk_polydata(path)
        if points.size and faces.size:
            meshes.append(_mesh_from_surface(points, faces, shape_family=path.stem))
    if not meshes:
        raise DataUnavailableError(
            _download_hint("5801337", "Roney et al. LA virtual cohort")
        )
    return meshes


def load_fibre_atlas(root: str) -> list[AtrialMesh]:
    """Load the Roney et al. Human Atrial Fibre Atlas (CARP ``.pts/.elem/.lon``).

    Parameters
    ----------
    root : str
        Directory containing the extracted Zenodo 3764917 archive.

    Returns
    -------
    list[AtrialMesh]

    Raises
    ------
    DataUnavailableError
        If ``root`` is missing or empty. The message names Zenodo record 3764917
        and the download steps.
    """
    if not _has_data(root):
        raise DataUnavailableError(
            _download_hint("3764917", "Roney et al. Human Atrial Fibre Atlas")
        )
    meshes: list[AtrialMesh] = []
    for pts_path in _list_by_suffix(root, (".pts",)):
        elem_path = pts_path.with_suffix(".elem")
        lon_path = pts_path.with_suffix(".lon")
        if not elem_path.exists():
            continue
        points, faces, _fibres = _read_carp(
            pts_path, elem_path, lon_path if lon_path.exists() else None
        )
        if points.size and faces.size:
            meshes.append(_mesh_from_surface(points, faces, shape_family=pts_path.stem))
    if not meshes:
        raise DataUnavailableError(
            _download_hint("3764917", "Roney et al. Human Atrial Fibre Atlas")
        )
    return meshes
