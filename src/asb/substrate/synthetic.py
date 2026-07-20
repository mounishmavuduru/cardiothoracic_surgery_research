"""Synthetic atrial anatomy generation (no patient data).

Produces closed, left-atrium-like triangulated surfaces with a smooth
rule-based fibre field, UAC-like ``(alpha, beta)`` coordinates in ``[0, 1]``,
region labels, and patchy fibrosis. Everything is fully deterministic in the
supplied ``seed``: the topology comes from a subdivided icosahedron (an
*icosphere*, guaranteeing a valid closed surface in which every vertex is used
by at least one face), which is then radially deformed into an ellipsoidal body
plus a few pulmonary-vein stub bumps.

These meshes are a development stand-in for real open cardiac geometries
(see :mod:`asb.substrate.loaders`); they are not patient data.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from asb.config import CohortConfig
from asb.types import AtrialMesh

__all__ = [
    "make_base_atrium",
    "shape_variant",
    "paint_fibrosis",
    "make_cohort",
]


# --------------------------------------------------------------------------- #
# Icosphere construction
# --------------------------------------------------------------------------- #
def _icosahedron() -> tuple[np.ndarray, np.ndarray]:
    """Unit-radius icosahedron: 12 vertices, 20 triangular faces."""
    t = (1.0 + np.sqrt(5.0)) / 2.0
    verts = np.array(
        [
            [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
            [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
            [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
        ],
        dtype=float,
    )
    verts /= np.linalg.norm(verts, axis=1, keepdims=True)
    faces = np.array(
        [
            [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
            [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
            [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
            [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
        ],
        dtype=np.int64,
    )
    return verts, faces


def _subdivide(
    verts: np.ndarray, faces: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """One loop of triangle midpoint subdivision, re-projected to unit sphere."""
    verts = list(verts)
    cache: dict[tuple[int, int], int] = {}

    def midpoint(a: int, b: int) -> int:
        key = (a, b) if a < b else (b, a)
        idx = cache.get(key)
        if idx is not None:
            return idx
        m = 0.5 * (verts[a] + verts[b])
        m = m / np.linalg.norm(m)
        idx = len(verts)
        verts.append(m)
        cache[key] = idx
        return idx

    new_faces = []
    for f0, f1, f2 in faces:
        a = midpoint(int(f0), int(f1))
        b = midpoint(int(f1), int(f2))
        c = midpoint(int(f2), int(f0))
        new_faces.extend([[f0, a, c], [f1, b, a], [f2, c, b], [a, b, c]])

    return np.asarray(verts, dtype=float), np.asarray(new_faces, dtype=np.int64)


def _icosphere(n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Icosphere with vertex count closest to ``n_nodes``.

    An icosphere at subdivision level ``L`` has ``10 * 4**L + 2`` vertices.
    """
    counts = [10 * 4 ** L + 2 for L in range(0, 8)]
    level = int(np.argmin([abs(c - n_nodes) for c in counts]))
    verts, faces = _icosahedron()
    for _ in range(level):
        verts, faces = _subdivide(verts, faces)
    return verts, faces


# --------------------------------------------------------------------------- #
# Geometry helpers
# --------------------------------------------------------------------------- #
def _vertex_normals(points: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Area-weighted unit vertex normals."""
    p0 = points[faces[:, 0]]
    p1 = points[faces[:, 1]]
    p2 = points[faces[:, 2]]
    face_n = np.cross(p1 - p0, p2 - p0)  # area-weighted face normals

    normals = np.zeros_like(points)
    for k in range(3):
        np.add.at(normals, faces[:, k], face_n)

    norm = np.linalg.norm(normals, axis=1, keepdims=True)
    norm = np.where(norm > 0, norm, 1.0)
    return normals / norm


def _fibre_field(points: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Smooth rule-based tangential fibre field (unit vectors per vertex).

    Fibres follow the tangential projection of a global preferred direction
    (a proxy for the circumferential myofibre arrangement); near the poles the
    projection is degenerate and a fallback axis is used.
    """
    normals = _vertex_normals(points, faces)
    d = np.array([0.0, 0.0, 1.0])
    fib = np.broadcast_to(d, points.shape).copy()
    fib = fib - (np.sum(fib * normals, axis=1, keepdims=True)) * normals

    mag = np.linalg.norm(fib, axis=1)
    degenerate = mag < 1e-6
    if np.any(degenerate):
        d2 = np.array([1.0, 0.0, 0.0])
        alt = np.broadcast_to(d2, points.shape).copy()
        alt = alt - (np.sum(alt * normals, axis=1, keepdims=True)) * normals
        fib[degenerate] = alt[degenerate]
        mag = np.linalg.norm(fib, axis=1)

    mag = np.where(mag > 0, mag, 1.0)
    return fib / mag[:, None]


def _uac(dirs: np.ndarray) -> np.ndarray:
    """UAC-like (alpha, beta) coords in [0, 1] from unit direction vectors."""
    x, y, z = dirs[:, 0], dirs[:, 1], dirs[:, 2]
    azimuth = np.arctan2(y, x)  # [-pi, pi]
    polar = np.arccos(np.clip(z, -1.0, 1.0))  # [0, pi]
    alpha = (azimuth + np.pi) / (2.0 * np.pi)
    beta = polar / np.pi
    return np.clip(np.stack([alpha, beta], axis=1), 0.0, 1.0)


def _regions(uac: np.ndarray) -> np.ndarray:
    """Integer region labels from a coarse UAC grid (3 beta bands x 2 alpha)."""
    alpha, beta = uac[:, 0], uac[:, 1]
    beta_band = np.clip((beta * 3).astype(np.int64), 0, 2)
    alpha_band = np.clip((alpha * 2).astype(np.int64), 0, 1)
    return (beta_band * 2 + alpha_band).astype(np.int64)


def _pv_centers(rng: np.random.Generator, n_pv: int) -> np.ndarray:
    """Unit direction vectors for pulmonary-vein stubs (upper hemisphere)."""
    centers = []
    for _ in range(n_pv):
        azimuth = rng.uniform(0.0, 2.0 * np.pi)
        polar = rng.uniform(0.15 * np.pi, 0.45 * np.pi)  # upper region
        centers.append(
            [
                np.sin(polar) * np.cos(azimuth),
                np.sin(polar) * np.sin(azimuth),
                np.cos(polar),
            ]
        )
    return np.asarray(centers, dtype=float)


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def make_base_atrium(
    seed: int, n_nodes: int = 1200, family: str | None = None
) -> AtrialMesh:
    """Generate a closed LA-like triangulated surface, deterministic in ``seed``.

    The mesh is an icosphere radially deformed into an ellipsoidal atrial body
    with a few pulmonary-vein stub bumps, carrying a smooth rule-based fibre
    field, UAC-like ``(alpha, beta)`` coordinates in ``[0, 1]``, and region
    labels. Fibrosis is initialised to zero (paint it with
    :func:`paint_fibrosis`).

    Parameters
    ----------
    seed : int
        Seed controlling anatomy perturbations (fully deterministic).
    n_nodes : int, optional
        Approximate vertex count; the nearest icosphere resolution is used.
    family : str or None, optional
        Shape-family id stored on the mesh (for grouped CV). Defaults to
        ``f"family_seed{seed}"``.

    Returns
    -------
    AtrialMesh
        Valid closed surface (every vertex used by >= 1 face).
    """
    rng = np.random.default_rng(seed)
    unit, faces = _icosphere(n_nodes)

    # Ellipsoidal body: per-axis radii perturbed by seed.
    axes = np.array([1.0, 0.85, 0.7]) * (1.0 + rng.uniform(-0.1, 0.1, size=3))
    points = unit * axes[None, :]

    # Pulmonary-vein stub bumps: radial outward Gaussian pushes.
    n_pv = 4
    centers = _pv_centers(rng, n_pv)
    for c in centers:
        cos_ang = unit @ c
        bump = np.exp(-((1.0 - cos_ang) / 0.05))  # tight cap around the center
        amp = rng.uniform(0.25, 0.4)
        points = points + amp * bump[:, None] * unit

    # Fields derived from the (unit) direction so they stay smooth and stable.
    fibres = _fibre_field(points, faces)
    uac = _uac(unit)
    region = _regions(uac)
    fibrosis = np.zeros(points.shape[0], dtype=float)

    fam = family if family is not None else f"family_seed{int(seed)}"
    return AtrialMesh(
        points=points,
        faces=faces,
        fibres=fibres,
        uac=uac,
        fibrosis=fibrosis,
        region=region,
        shape_family=fam,
        meta={"unit_dirs": unit, "seed": int(seed)},
    )


def shape_variant(base: AtrialMesh, seed: int) -> AtrialMesh:
    """Perturb the shape modes of ``base``, keeping the same ``shape_family``.

    Low-frequency radial displacements (a proxy for statistical-shape modes) are
    applied along the vertex normals. Topology is unchanged, so the surface stays
    valid; ``fibres`` are recomputed for the deformed geometry while ``uac``,
    ``region`` and ``fibrosis`` are carried through.

    Parameters
    ----------
    base : AtrialMesh
        Parent anatomy.
    seed : int
        Seed for the mode perturbation.

    Returns
    -------
    AtrialMesh
        Variant with ``shape_family == base.shape_family``.
    """
    rng = np.random.default_rng(seed)
    points = np.array(base.points, dtype=float, copy=True)
    normals = _vertex_normals(points, base.faces)

    # A handful of low-frequency modes as functions of direction angles.
    center = points.mean(axis=0)
    dirs = points - center
    dirs = dirs / np.maximum(np.linalg.norm(dirs, axis=1, keepdims=True), 1e-12)
    azimuth = np.arctan2(dirs[:, 1], dirs[:, 0])
    polar = np.arccos(np.clip(dirs[:, 2], -1.0, 1.0))

    disp = np.zeros(points.shape[0], dtype=float)
    for _ in range(4):
        k1 = rng.integers(1, 4)
        k2 = rng.integers(1, 4)
        phase = rng.uniform(0.0, 2.0 * np.pi)
        amp = rng.uniform(-0.08, 0.08)
        disp += amp * np.sin(k1 * azimuth + phase) * np.cos(k2 * polar)

    scale = float(np.mean(np.linalg.norm(points - center, axis=1)))
    points = points + (disp * scale)[:, None] * normals

    fibres = _fibre_field(points, base.faces)
    return replace(
        base,
        points=points,
        fibres=fibres,
        uac=np.array(base.uac, copy=True),
        fibrosis=np.array(base.fibrosis, copy=True),
        region=np.array(base.region, copy=True),
        meta={**dict(base.meta), "variant_seed": int(seed)},
    )


def paint_fibrosis(
    mesh: AtrialMesh,
    seed: int,
    burden: float = 0.2,
    n_patches: int = 5,
    patch_scale: float = 0.15,
) -> AtrialMesh:
    """Return a copy of ``mesh`` with a patchy fibrosis field in ``[0, 1]``.

    A few Gaussian patches centred on random vertices are summed into a smooth
    field, then rescaled so the mean fibrosis matches ``burden`` and clipped to
    ``[0, 1]``. Deterministic in ``seed``.

    Parameters
    ----------
    mesh : AtrialMesh
        Source mesh (left unmodified).
    seed : int
        Seed for patch placement.
    burden : float, optional
        Target mean fibrosis fraction in ``[0, 1]``.
    n_patches : int, optional
        Number of Gaussian fibrotic patches.
    patch_scale : float, optional
        Spatial scale of each patch, as a fraction of the mesh diagonal.

    Returns
    -------
    AtrialMesh
        Copy with ``.fibrosis`` set in ``[0, 1]``.
    """
    rng = np.random.default_rng(seed)
    points = np.asarray(mesh.points, dtype=float)
    n = points.shape[0]

    diag = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0)))
    sigma = max(patch_scale * diag, 1e-9)

    field = np.zeros(n, dtype=float)
    if n_patches > 0 and n > 0:
        centers = rng.choice(n, size=min(n_patches, n), replace=False)
        for c in centers:
            d2 = np.sum((points - points[c]) ** 2, axis=1)
            field += rng.uniform(0.5, 1.0) * np.exp(-d2 / (2.0 * sigma ** 2))

    fmax = field.max()
    if fmax > 0:
        field = field / fmax  # into [0, 1]
        mean = field.mean()
        if mean > 0:
            field = field * (float(burden) / mean)
    fibrosis = np.clip(field, 0.0, 1.0)

    return replace(mesh, fibrosis=fibrosis, meta={**dict(mesh.meta),
                                                  "fibrosis_seed": int(seed)})


def make_cohort(cfg: CohortConfig) -> list[AtrialMesh]:
    """Generate a synthetic cohort: ``n_base`` families x ``n_variants`` each.

    For each base family a base atrium is generated; variant 0 is the base and
    the remaining variants are shape-mode perturbations of it. Every mesh is
    fibrosis-painted and tagged with its base ``shape_family``.

    Parameters
    ----------
    cfg : CohortConfig
        Cohort knobs (``n_base``, ``n_variants``, ``n_nodes``, ``seed``).

    Returns
    -------
    list[AtrialMesh]
        ``n_base * n_variants`` meshes.
    """
    cohort: list[AtrialMesh] = []
    for b in range(cfg.n_base):
        family = f"family_{b}"
        base = make_base_atrium(
            seed=cfg.seed + b, n_nodes=cfg.n_nodes, family=family
        )
        for v in range(cfg.n_variants):
            # Distinct, deterministic per-(base, variant) seeds.
            vseed = (cfg.seed + b) * 10_000 + v
            mesh = base if v == 0 else shape_variant(base, seed=vseed)
            mesh = paint_fibrosis(mesh, seed=vseed + 7)
            cohort.append(mesh)
    return cohort
