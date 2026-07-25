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
  continuous IIR. The Dryad README does not document the tag values; the semantics
  below were established by direct measurement over the released meshes
  (2026-07-24, ``notebooks/lab_notebook.md``):

  =======  =========================  ======================================
  elemTag  tissue                     evidence
  =======  =========================  ======================================
  111      healthy myocardium         dominant; 59.0 % → 47.3 % post-ablation
  115      dense (LGE) fibrosis       high inter-patient variance; 20.9 % → 14.9 %
  164      **not myocardium** — the   4–6 large components (5 in 75/82), each a disc
           caps over the pulmonary    (χ=1), 11–42 mm; borders tag 115 on ~0.08 %
           veins and mitral valve     of edges vs 14–20 % by chance; 20.0 % → 19.9 %
  199      ablation scar              post-ablation only, 17.9 %
  =======  =========================  ======================================

  An earlier revision of this file read tag 164 as "remodelled / patchy tissue,
  distributed (least compact)" and :data:`TAG_FIBROSIS` still maps it to 0.5.
  That reading is wrong: the region is maximally compact -- four to six large
  components (five in 75 of the 82 meshes, four in five, six in two, matching normal
  pulmonary-vein variation) -- and it borders fibrosis on only ~0.08% of its incident
  edges against a 14-20% chance rate. :data:`DROP_TAGS` therefore removes tag 164
  before any field is derived, which reopens the venous and mitral orifices (Euler
  characteristic falls by a median of 7, range 2-12).
  ``TAG_FIBROSIS[164]`` is retained only so the discredited mapping can still be
  reproduced via ``drop_tags=()`` for the record.
  :func:`load_uw_mesh` accepts ``tag_fibrosis`` overrides so the SFI-vs-fibrosis
  conclusion can be sensitivity-tested against the tag interpretation.
- **No UAC and no per-cell UAC** are shipped, so ``uac`` is a *surrogate*: the
  top-2 principal axes of the surface coordinates, min-max normalised to
  ``[0, 1]``. This is a smooth 2-D chart adequate for region binning and the
  spatial null; it is NOT anatomically registered UAC. Flagged in ``meta``.
- A per-cell ``fiber`` array **is** shipped, but in all 164 released meshes it is a
  constant ``(1, 0, 0)``: measured mean directional spread is exactly ``0.0``. The
  anisotropy the solver then applies is a fixed coordinate bias, not anatomy. See
  :data:`DEGENERATE_FIBRE_SPREAD`; the condition is flagged in ``meta`` and warned about.

Coordinates are microns (LA span ~10 cm) → mm via :data:`MICRON_TO_MM`.
"""
from __future__ import annotations

import glob
import os
import warnings
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

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

#: Mean per-vertex deviation from the mean fibre direction, below which the field is
#: treated as degenerate (no anatomical anisotropy). The released UW meshes score ~0.
DEGENERATE_FIBRE_SPREAD: float = 1.0e-6

#: Element tags that are NOT myocardium and must be removed before simulation.
#:
#: Tag 164 was identified geometrically on 2026-07-24 (the Dryad README does not
#: document the tag values). Across every mesh examined it:
#:   * resolves into FOUR TO SIX large connected components (>=1% of the region) in all
#:     82 meshes -- five in 75, four in five, six in two. That is the usual variation in
#:     pulmonary-vein anatomy: a left common trunk gives four, a right middle vein six;
#:   * each component is topologically a disc (Euler characteristic chi = 1), i.e. a
#:     cap sealing an opening, 11-42 mm across, sitting 0.5-0.8 of the atrial radius
#:     from the centroid;
#:   * borders fibrotic tag-115 elements on only ~0.08% of its incident edges (median 27
#:     edges) against a chance expectation of 14-20% -- a >100-fold depletion, which rules
#:     out the "fibrosis border zone" reading that ``TAG_FIBROSIS[164]=0.5`` implies;
#:   * occupies 20.0 % of elements pre-ablation and 19.9 % post -- untouched by the
#:     procedure, as anatomy is and as tissue is not;
#:   * removing it opens the surface: Euler characteristic falls by a median of 7
#:     (range 2-12). ID040 and ID049 are topologically pathological and go the other way.
#:
#: An earlier revision of this comment claimed "exactly five components", "EXACTLY ZERO
#: adjacency" and "chi = -3". Those were generalised from one mesh and from a rounded mean;
#: the census figures above replace them. See ``tests/test_uw_boyle.py``.
#:
#: Leaving these caps in as half-conducting tissue does more than distort fibrosis: it
#: seals the atrium's openings, so activation can cross the mitral valve and the vein
#: ostia instead of travelling around them. Reentry anchored on those orifices -- a
#: principal atrial-fibrillation mechanism -- then cannot form at all.
DROP_TAGS: Tuple[int, ...] = (164,)


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
    drop_tags: Sequence[int] = DROP_TAGS,
) -> AtrialMesh:
    """Load one UW/Boyle VTK mesh into an :class:`AtrialMesh` with real fields.

    Coordinates microns→mm, fibrosis from the per-cell ``elemTag`` via
    ``tag_fibrosis`` mapped to vertices, fibres per-cell→per-vertex, UAC a PCA
    surrogate, regions a 3×3 grid over the surrogate UAC. ``shape_family``
    defaults to ``uw_<stem>`` (one group per patient/model).

    ``drop_tags`` removes element classes that are not myocardium before anything
    else is derived; see :data:`DROP_TAGS`. Pass ``drop_tags=()`` to reproduce the
    pre-2026-07-24 behaviour in which the openings were left in as conducting tissue.
    """
    parsed = read_vtk_unstructured_bin(path)
    points = np.asarray(parsed["points"], dtype=float) * MICRON_TO_MM
    faces = np.asarray(parsed["faces"], dtype=np.int64)
    n = points.shape[0]
    cs: Dict[str, np.ndarray] = parsed["cell_scalars"]  # type: ignore[assignment]
    cv: Dict[str, np.ndarray] = parsed["cell_vectors"]  # type: ignore[assignment]

    tags = np.asarray(cs.get("elemTag", np.full(faces.shape[0], 111.0)))
    fib_cell_raw = cv.get("fiber", None)

    # Drop non-myocardial elements (the caps over the PV and mitral-valve openings)
    # BEFORE deriving any field, then compact the vertex indexing.
    n_dropped = 0
    if len(drop_tags):
        keep = ~np.isin(tags.astype(int), np.asarray(drop_tags, dtype=int))
        n_dropped = int((~keep).sum())
        if n_dropped:
            faces, tags = faces[keep], tags[keep]
            if fib_cell_raw is not None:
                fib_cell_raw = np.asarray(fib_cell_raw)[keep]
            used = np.unique(faces)
            remap = np.full(n, -1, dtype=np.int64)
            remap[used] = np.arange(used.size)
            faces = remap[faces]
            points = points[used]
            n = used.size

    fibrosis = tags_to_fibrosis(tags, faces, n, mapping=tag_fibrosis)

    # Fibre field. Verified 2026-07-24 across all 164 released meshes: the ``VECTORS
    # fiber`` array IS present, but it is a constant (1, 0, 0) for every element of every
    # mesh -- i.e. the public Dryad deposit carries a degenerate fibre field, apparently
    # lost when the meshes were downsampled for release. This is upstream data, not a
    # local fallback. It matters because ``edge_weights_from_fibres`` makes conduction
    # anisotropic (along 1.0 / cross 0.3) relative to the local fibre direction: with one
    # global direction the anisotropy becomes a fixed coordinate bias with no anatomical
    # meaning, and all fibre heterogeneity -- a primary substrate for unidirectional block
    # and hence reentry initiation -- disappears. Measured contribution: see
    # results/uw_substrate_ablation.json. Repairing this AND the sealed orifices together
    # lifts inducibility 6/82 -> 14/82, but that is not significant (McNemar p=0.077) and
    # closes only 39% of the gap to Roney's 32.3%. The direct control settles the fibre
    # half: results/roney_fibre_control.json swaps the real Roney fibres for a constant
    # field and inducibility is 20/62 either way (McNemar p=1.0), so the fibre field does
    # not drive the aggregate rate at all -- though per-subject dynamics do move (sustained
    # duration differs on 50/62 subjects, 10/62 verdicts flip). The degenerate fibres are a
    # fidelity defect, not an established cause of the low UW rate, which remains open.
    fib_cell = fib_cell_raw
    if fib_cell is None:
        fibres = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
        fibre_spread = 0.0
    else:
        fibres = _cell_vectors_to_vertex(n, faces, np.asarray(fib_cell, dtype=float))
        fibre_spread = float(np.mean(np.linalg.norm(fibres - fibres.mean(axis=0), axis=1)))
    fibres_are_degenerate = bool(fibre_spread < DEGENERATE_FIBRE_SPREAD)
    if fibres_are_degenerate:
        warnings.warn(
            f"{Path(path).name}: fibre field is degenerate (directional spread "
            f"{fibre_spread:.2e}); every element shares one direction. Anisotropy is a "
            "fixed coordinate bias, not anatomy -- do not report inducibility from this "
            "substrate as a measurement of the patient.",
            RuntimeWarning, stacklevel=2,
        )

    uac = _pca_uac_surrogate(points)
    region = _uac_regions(uac)
    fam = shape_family if shape_family is not None else f"uw_{Path(path).stem}"

    u, c = np.unique(np.asarray(tags, int), return_counts=True)
    return AtrialMesh(
        points=points, faces=faces, fibres=fibres, uac=uac,
        fibrosis=fibrosis, region=region, shape_family=fam,
        meta={"source": "uw_boyle_dryad_kkwh70sg0", "path": os.path.abspath(path),
              "uac_is_surrogate": True, "fibrosis_from": "elemTag",
              "fibres_are_degenerate": fibres_are_degenerate,
              "dropped_tags": tuple(int(t) for t in drop_tags),
              "n_cells_dropped": int(n_dropped),
              "fibre_spread": fibre_spread,
              "tag_fractions": {int(a): float(b / c.sum()) for a, b in zip(u, c)},
              "fibrosis_mean": float(np.mean(fibrosis)), "n_full": int(n)},
    )


def uw_preablation_paths(root: str = "data/uw_boyle") -> List[str]:
    """Sorted list of all ``*_PreAbl_*.vtk`` mesh paths (one substrate per patient)."""
    paths = glob.glob(os.path.join(root, "**", "*_PreAbl_*.vtk"), recursive=True)
    return sorted(paths)
