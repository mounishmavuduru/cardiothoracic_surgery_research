"""openCARP mesh export, pacing protocols, and solver runner (INTERFACE ONLY).

This module is the clean seam between AtrialSpectralBench and the real openCARP
cardiac-electrophysiology solver (https://opencarp.org). It can:

- write an :class:`asb.types.AtrialMesh` to the openCARP text mesh triple
  (``.pts`` / ``.elem`` / ``.lon``) — pure file I/O, always available;
- build a pacing protocol dictionary from a :class:`asb.config.LabelConfig`;
- *run* openCARP — but the binary is a heavy external dependency and is **not
  installed in the build/CI environment** (the network is blocked at build
  time), so :func:`run_opencarp` raises :class:`OpenCARPUnavailableError` with
  an install pointer rather than executing anything here.

When openCARP is present the runner would return a real ``'opencarp'``-sourced
:class:`~asb.types.InducibilityLabel` (true ground truth), replacing the
``mock_ep`` development stand-in.
"""
from __future__ import annotations

import os
import shutil
from typing import Optional

import numpy as np

from asb.config import LabelConfig
from asb.types import AtrialMesh, InducibilityLabel

__all__ = [
    "OpenCARPUnavailableError",
    "write_carp_mesh",
    "make_pacing_protocol",
    "run_opencarp",
    "reproduce_niederer_benchmark",
]

_INSTALL_POINTER = (
    "openCARP is not installed / not on PATH. Install it from "
    "https://opencarp.org (see the 'Download' page for binary packages and "
    "https://git.opencarp.org/openCARP/openCARP for source), then pass the "
    "solver path via `opencarp_bin=` or put `openCARP` on your PATH. openCARP "
    "is intentionally NOT fetched or run in this build environment (network is "
    "blocked at build time)."
)

# Candidate executable names for the openCARP monodomain/bidomain solver.
_CARP_BINARIES = ("openCARP", "carp.pt", "carp.petsc", "carp")


class OpenCARPUnavailableError(RuntimeError):
    """Raised when the openCARP solver binary cannot be located or run.

    The message always includes an install pointer so a downstream user knows
    exactly how to provide the real solver.
    """


def write_carp_mesh(mesh: AtrialMesh, out_dir: str) -> dict:
    """Write an :class:`~asb.types.AtrialMesh` as an openCARP text mesh triple.

    Emits three sibling files in ``out_dir``:

    - ``mesh.pts`` : header line ``n_points`` then one ``x y z`` per vertex.
    - ``mesh.elem``: header line ``n_elems`` then one ``Tr v0 v1 v2 region``
      per triangle (``Tr`` = triangle element tag; the trailing field is the
      element region label, taken from the first vertex of each triangle).
    - ``mesh.lon`` : header line ``1`` (one fibre direction per element) then
      one ``fx fy fz`` per triangle — the unit-normalized mean of its vertex
      fibre vectors.

    This is pure file I/O; it neither runs nor requires openCARP.

    Parameters
    ----------
    mesh : AtrialMesh
        Surface to export.
    out_dir : str
        Output directory (created if absent).

    Returns
    -------
    dict
        ``{"pts", "elem", "lon"}`` absolute file paths plus ``"n_points"`` and
        ``"n_elems"`` counts.
    """
    os.makedirs(out_dir, exist_ok=True)

    points = np.asarray(mesh.points, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    fibres = np.asarray(mesh.fibres, dtype=float)
    region = np.asarray(mesh.region, dtype=np.int64)

    n_points = int(points.shape[0])
    n_elems = int(faces.shape[0])

    pts_path = os.path.abspath(os.path.join(out_dir, "mesh.pts"))
    elem_path = os.path.abspath(os.path.join(out_dir, "mesh.elem"))
    lon_path = os.path.abspath(os.path.join(out_dir, "mesh.lon"))

    # --- .pts ---
    with open(pts_path, "w") as fh:
        fh.write(f"{n_points}\n")
        for x, y, z in points:
            fh.write(f"{x:.9g} {y:.9g} {z:.9g}\n")

    # --- .elem ---
    with open(elem_path, "w") as fh:
        fh.write(f"{n_elems}\n")
        for tri in faces:
            v0, v1, v2 = int(tri[0]), int(tri[1]), int(tri[2])
            tag = int(region[v0]) if n_points else 0
            fh.write(f"Tr {v0} {v1} {v2} {tag}\n")

    # --- .lon ---
    with open(lon_path, "w") as fh:
        fh.write("1\n")
        for tri in faces:
            f = fibres[tri].mean(axis=0)
            norm = float(np.linalg.norm(f))
            if norm > 0:
                f = f / norm
            else:
                f = np.array([1.0, 0.0, 0.0])
            fh.write(f"{f[0]:.9g} {f[1]:.9g} {f[2]:.9g}\n")

    return {
        "pts": pts_path,
        "elem": elem_path,
        "lon": lon_path,
        "n_points": n_points,
        "n_elems": n_elems,
    }


def make_pacing_protocol(cfg: LabelConfig) -> dict:
    """Build an openCARP-style pacing protocol dictionary from a config.

    Deterministic; a fixed ``cfg`` gives a fixed protocol. Timing values are in
    milliseconds and follow standard clinical inducibility protocols (an S1
    drive train followed by a premature S2 extrastimulus, or a rapid burst).

    Parameters
    ----------
    cfg : LabelConfig
        Uses ``protocol`` (``'S1S2'`` or ``'burst'``) and ``n_pacing_sites``.

    Returns
    -------
    dict
        Protocol specification: ``name``, ``n_sites``, ``s1_cycle_length_ms``,
        ``n_s1``, and either ``s2_coupling_ms`` (for ``'S1S2'``) or
        ``burst_cycle_length_ms`` / ``n_burst`` (for ``'burst'``).
    """
    protocol = cfg.protocol
    base = {
        "name": protocol,
        "n_sites": int(cfg.n_pacing_sites),
        "s1_cycle_length_ms": 600.0,
        "n_s1": 6,
        "stim_strength_uA_per_cm2": 100.0,
        "stim_duration_ms": 2.0,
    }
    if protocol == "burst":
        base.update({"burst_cycle_length_ms": 130.0, "n_burst": 10})
    else:  # 'S1S2' (default) and any other name fall back to an S1-S2 train.
        base.update({"s2_coupling_ms": 250.0, "s2_step_ms": -10.0})
    return base


def _find_opencarp(opencarp_bin: Optional[str]) -> Optional[str]:
    """Resolve an openCARP executable path, or ``None`` if not found."""
    if opencarp_bin:
        if os.path.isfile(opencarp_bin) and os.access(opencarp_bin, os.X_OK):
            return opencarp_bin
        resolved = shutil.which(opencarp_bin)
        return resolved
    for name in _CARP_BINARIES:
        resolved = shutil.which(name)
        if resolved:
            return resolved
    return None


def run_opencarp(
    mesh_dir: str,
    protocol: dict,
    *,
    opencarp_bin: Optional[str] = None,
) -> InducibilityLabel:
    """Run openCARP on an exported mesh and return a real inducibility label.

    In this build/CI environment the solver is **not installed** and the
    network is blocked, so this function always raises
    :class:`OpenCARPUnavailableError`. When openCARP is present it would launch
    the monodomain solve under ``protocol``, detect self-sustaining reentry, and
    return an :class:`~asb.types.InducibilityLabel` with ``source='opencarp'``
    (true ground truth).

    Parameters
    ----------
    mesh_dir : str
        Directory containing the ``.pts`` / ``.elem`` / ``.lon`` triple written
        by :func:`write_carp_mesh`.
    protocol : dict
        Pacing protocol from :func:`make_pacing_protocol`.
    opencarp_bin : str, optional
        Explicit path to (or name of) the openCARP executable. If ``None`` the
        standard binary names are searched on ``PATH``.

    Returns
    -------
    InducibilityLabel
        Real ``'opencarp'``-sourced label (only when the solver is available).

    Raises
    ------
    OpenCARPUnavailableError
        If the openCARP binary cannot be located (always, in this environment).
    FileNotFoundError
        If ``mesh_dir`` does not exist.
    """
    if not os.path.isdir(mesh_dir):
        raise FileNotFoundError(f"mesh_dir does not exist: {mesh_dir!r}")

    binary = _find_opencarp(opencarp_bin)
    if binary is None:
        raise OpenCARPUnavailableError(_INSTALL_POINTER)

    # Real invocation is intentionally not implemented in this environment.
    raise OpenCARPUnavailableError(
        "openCARP binary was located but the solver invocation is not wired up "
        "in this build environment. " + _INSTALL_POINTER
    )


def reproduce_niederer_benchmark(
    out_dir: Optional[str] = None,
    *,
    opencarp_bin: Optional[str] = None,
) -> dict:
    """Deferred code-verification target: the Niederer et al. (2011) benchmark.

    Niederer et al., *Phil. Trans. R. Soc. A* 369:4331-4351 (2011), defines a
    standard cardiac-tissue activation benchmark (a 20x7x3 mm slab paced from a
    corner) used to verify that a monodomain solver reproduces the reference
    activation times across codes and resolutions. Running it here would
    establish that our openCARP wiring is quantitatively correct before any
    ``'opencarp'`` labels are trusted.

    This is a **documented deferred target**: it requires the openCARP solver,
    which is not installed in this environment, so it raises
    :class:`OpenCARPUnavailableError`. The signature and docstring pin down the
    intended verification for when the solver becomes available.

    Parameters
    ----------
    out_dir : str, optional
        Where benchmark outputs would be written.
    opencarp_bin : str, optional
        Explicit openCARP executable path/name; searched on ``PATH`` if ``None``.

    Returns
    -------
    dict
        Benchmark result summary (only when openCARP is available).

    Raises
    ------
    OpenCARPUnavailableError
        Always, in this environment (solver absent; deferred run).
    """
    binary = _find_opencarp(opencarp_bin)
    if binary is None:
        raise OpenCARPUnavailableError(
            "Niederer (2011) benchmark reproduction is deferred: it requires "
            "the openCARP solver. " + _INSTALL_POINTER
        )
    raise OpenCARPUnavailableError(
        "Niederer benchmark reproduction is a deferred code-verification "
        "target and is not implemented in this build environment. "
        + _INSTALL_POINTER
    )
