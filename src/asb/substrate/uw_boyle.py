r"""UW/Boyle LGE-MRI left-atrial cohort loader (Dryad ``10.5061/dryad.kkwh70sg0``).

Loads the University of Washington / Boyle-lab atrial meshes — 82 distinct AF
patients segmented from late-gadolinium-enhanced MRI, two models each
(pre-ablation fibrotic substrate, post-ablation with procedural scar) — into
:class:`asb.types.AtrialMesh` objects, so they extend the real cohort alongside
the Roney meshes through the *same* graph/spectral/EP pipeline.

Publication: Bifulco, Magoon, Chahine, Kim, Macheret, Akoum, Boyle (2025),
"Predicting arrhythmia recurrence post-ablation in atrial fibrillation using
explainable machine learning" (meshes CC0 via Dryad).

Format differences from the Roney loader (``asb.substrate.roney``)
-----------------------------------------------------------------
- Files are **binary legacy VTK** ``UNSTRUCTURED_GRID`` (triangle cells,
  ``celltype 5``), big-endian, not ASCII ``POLYDATA``. :func:`read_vtk_unstructured_bin`
  parses them directly (``meshio`` rejects the int ``elemTag`` block).
- Fibrosis is **categorical** via the per-cell ``elemTag`` field, not a
  continuous IIR. Tag semantics, established empirically over all 82 patients
  (variance across patients + invariance pre/post ablation, see
  ``notebooks/lab_notebook.md``):

  =======  ==========================  =====================================
  elemTag  tissue                      evidence
  =======  ==========================  =====================================
  111      healthy myocardium          dominant; drops post-ablation
  115      dense (LGE) fibrosis         high inter-patient variance; drops post-abl
  164      remodelled / patchy tissue  distributed (least compact); ablation-invariant
  199      ablation scar               post-ablation only
  =======  ==========================  =====================================

  The frozen mapping :data:`TAG_FIBROSIS` sends healthy→0, dense fibrosis→1,
  remodelled→0.5, scar→1. :func:`load_uw_mesh` accepts ``tag_fibrosis`` overrides
  so the SFI-vs-fibrosis conclusion can be sensitivity-tested against the tag
  interpretation (e.g. 164→0 or 164→1).
- **No UAC and no per-cell UAC** are shipped, so ``uac`` is a *surrogate*: the
  top-2 principal axes of the surface coordinates, min-max normalised to
  ``[0, 1]``. This is a smooth 2-D chart adequate for region binning and the
  spatial null; it is NOT anatomically registered UAC. Flagged in ``meta``.
- Fibres **are** shipped (per-cell unit vectors) and are mapped to vertices.

Coordinates are microns (LA span ~10 cm) → mm via :data:`MICRON_TO_MM`.
"""
from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import Dict, List, Mapping, Optional

import numpy as np

from asb.types import AtrialMesh

__all__ = [
    "read_vtk_unstructured_bin",
    "TAG_FIBROSIS",
    "tags_to_fibrosis",
    "load_uw_mesh",
    "uw_preablation_paths",
    "MICRON_TO_MM",
]

#: UW meshes store coordinates in microns; convert to millimetres.
MICRON_TO_MM: float = 1.0e-3

#: Frozen per-cell ``elemTag`` -> fibrosis-fraction map (see module docstring).
TAG_FIBROSIS: Dict[int, float] = {111: 0.0, 115: 1.0, 164: 0.5, 199: 1.0}


# --------------------------------------------------------------------------- #
# Binary VTK UNSTRUCTURED_GRID parsing
# --------------------------------------------------------------------------- #
def read_vtk_unstructured_bin(path: str) -> Dict[str, object]:
    """Parse points, triangle faces, per-cell scalars and vectors from a binary VTK.

    Supports the big-endian legacy ``UNSTRUCTURED_GRID`` subset used by the UW
    meshes: ``POINTS`` + ``CELLS`` (triangles) + ``CELL_TYPES`` + ``CELL_DATA``
    (``SCALARS`` int/float + ``VECTORS`` float).

    Returns
    -------
    dict
        ``points`` (n, 3) float, ``faces`` (f, 3) int, ``cell_scalars``
        {name: (f,) float}, ``cell_vectors`` {name: (f, 3) float}.
    """
    raw = open(path, "rb").read()

    def line_at(off: int):
        end = raw.find(b"\n", off)
        return raw[off:end].decode("latin-1"), end + 1

    pos = 0
    for _ in range(4):  # 4 header lines: version, title, BINARY, DATASET ...
        _, pos = line_at(pos)

    out: Dict[str, object] = {"cell_scalars": {}, "cell_vectors": {}}
    ncells = 0
    while pos < len(raw):
        ln, nxt = line_at(pos)
        parts = ln.split()
        if not parts:
            pos = nxt
            continue
        key = parts[0].upper()
        if key == "POINTS":
            n = int(parts[1])
            cnt = n * 3
            arr = np.frombuffer(raw, dtype=">f4", count=cnt, offset=nxt).astype(np.float64)
            out["points"] = arr.reshape(n, 3)
            pos = nxt + cnt * 4
        elif key == "CELLS":
            ncells = int(parts[1])
            size = int(parts[2])
            arr = np.frombuffer(raw, dtype=">i4", count=size, offset=nxt).astype(np.int64)
            # UW cells are all triangles: [3, i, j, k] repeated.
            out["faces"] = arr.reshape(ncells, size // ncells)[:, 1:4]
            pos = nxt + size * 4
        elif key == "CELL_TYPES":
            pos = nxt + int(parts[1]) * 4
        elif key == "SCALARS":
            name = parts[1]
            dtype = parts[2] if len(parts) > 2 else "float"
            _, nxt2 = line_at(nxt)  # skip LOOKUP_TABLE line
            npdt = {"int": ">i4", "float": ">f4", "double": ">f8"}.get(dtype, ">f4")
            arr = np.frombuffer(raw, dtype=npdt, count=ncells, offset=nxt2).astype(np.float64)
            out["cell_scalars"][name] = arr  # type: ignore[index]
            pos = nxt2 + ncells * np.dtype(npdt).itemsize
        elif key == "VECTORS":
            name = parts[1]
            cnt = ncells * 3
            arr = np.frombuffer(raw, dtype=">f4", count=cnt, offset=nxt).astype(np.float64)
            out["cell_vectors"][name] = arr.reshape(ncells, 3)  # type: ignore[index]
            pos = nxt + cnt * 4
        else:  # POINT_DATA / CELL_DATA headers and anything else
            pos = nxt
    return out


# --------------------------------------------------------------------------- #
# Field derivation
# --------------------------------------------------------------------------- #
def tags_to_fibrosis(
    tags: np.ndarray, faces: np.ndarray, n_points: int,
    *, mapping: Mapping[int, float] = TAG_FIBROSIS, default: float = 0.0,
) -> np.ndarray:
    """Per-cell ``elemTag`` -> per-vertex fibrosis fraction (area-agnostic mean)."""
    cell_fib = np.array([mapping.get(int(t), default) for t in tags], dtype=float)
    acc = np.zeros(n_points, dtype=float)
    cnt = np.zeros(n_points, dtype=float)
    for t in range(faces.shape[0]):
        val = cell_fib[t]
        for k in range(3):
            acc[faces[t, k]] += val
            cnt[faces[t, k]] += 1.0
    nz = cnt > 0
    acc[nz] /= cnt[nz]
    return np.clip(acc, 0.0, 1.0)


def _cell_vectors_to_vertex(n_points: int, faces: np.ndarray, cell_vec: np.ndarray) -> np.ndarray:
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


def _pca_uac_surrogate(points: np.ndarray) -> np.ndarray:
    """Top-2 principal-axis chart of the surface, min-max normalised to [0,1]^2.

    A surrogate for universal atrial coordinates: smooth and monotone across the
    surface, adequate for coarse region binning and the spatial null, but NOT an
    anatomically registered UAC. Documented as such in ``meta``.
    """
    c = points - points.mean(axis=0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    proj = c @ vt[:2].T
    lo = proj.min(axis=0)
    rng = np.maximum(proj.max(axis=0) - lo, 1e-9)
    return np.clip((proj - lo) / rng, 0.0, 1.0)


def _uac_regions(uac: np.ndarray, n_alpha: int = 3, n_beta: int = 3) -> np.ndarray:
    """Coarse integer region labels from a UAC grid (matches roney.uac_regions)."""
    alpha = np.clip(uac[:, 0], 0.0, 1.0 - 1e-12)
    beta = np.clip(uac[:, 1], 0.0, 1.0 - 1e-12)
    ab = np.clip((alpha * n_alpha).astype(np.int64), 0, n_alpha - 1)
    bb = np.clip((beta * n_beta).astype(np.int64), 0, n_beta - 1)
    return (bb * n_alpha + ab).astype(np.int64)


# --------------------------------------------------------------------------- #
# Loader
# --------------------------------------------------------------------------- #
def load_uw_mesh(
    path: str, *, shape_family: Optional[str] = None,
    tag_fibrosis: Mapping[int, float] = TAG_FIBROSIS,
) -> AtrialMesh:
    """Load one UW/Boyle VTK mesh into an :class:`AtrialMesh` with real fields.

    Coordinates microns→mm, fibrosis from the per-cell ``elemTag`` via
    ``tag_fibrosis`` mapped to vertices, fibres per-cell→per-vertex, UAC a PCA
    surrogate, regions a 3×3 grid over the surrogate UAC. ``shape_family``
    defaults to ``uw_<stem>`` (one group per patient/model).
    """
    parsed = read_vtk_unstructured_bin(path)
    points = np.asarray(parsed["points"], dtype=float) * MICRON_TO_MM
    faces = np.asarray(parsed["faces"], dtype=np.int64)
    n = points.shape[0]
    cs: Dict[str, np.ndarray] = parsed["cell_scalars"]  # type: ignore[assignment]
    cv: Dict[str, np.ndarray] = parsed["cell_vectors"]  # type: ignore[assignment]

    tags = cs.get("elemTag", np.full(faces.shape[0], 111.0))
    fibrosis = tags_to_fibrosis(tags, faces, n, mapping=tag_fibrosis)

    fib_cell = cv.get("fiber", None)
    if fib_cell is None:
        fibres = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    else:
        fibres = _cell_vectors_to_vertex(n, faces, np.asarray(fib_cell, dtype=float))

    uac = _pca_uac_surrogate(points)
    region = _uac_regions(uac)
    fam = shape_family if shape_family is not None else f"uw_{Path(path).stem}"

    u, c = np.unique(np.asarray(tags, int), return_counts=True)
    return AtrialMesh(
        points=points, faces=faces, fibres=fibres, uac=uac,
        fibrosis=fibrosis, region=region, shape_family=fam,
        meta={"source": "uw_boyle_dryad_kkwh70sg0", "path": os.path.abspath(path),
              "uac_is_surrogate": True, "fibrosis_from": "elemTag",
              "tag_fractions": {int(a): float(b / c.sum()) for a, b in zip(u, c)},
              "fibrosis_mean": float(np.mean(fibrosis)), "n_full": int(n)},
    )


def uw_preablation_paths(root: str = "data/uw_boyle") -> List[str]:
    """Sorted list of all ``*_PreAbl_*.vtk`` mesh paths (one substrate per patient)."""
    paths = glob.glob(os.path.join(root, "**", "*_PreAbl_*.vtk"), recursive=True)
    return sorted(paths)
