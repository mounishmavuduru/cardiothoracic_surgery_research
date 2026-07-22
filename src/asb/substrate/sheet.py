r"""2-D triangulated sheet with a scar-bounded healthy isthmus (novelty experiment 1, real model).

The abstract graph models (``asb.transfer.bottleneck``) could not dissociate connectivity from
substrate because in a pure diffusive graph the Fiedler bottleneck *is* the low-conduction lesion.
A continuous reaction--diffusion sheet breaks that tie: a **healthy** (low-fibrosis) narrow isthmus
between two blocks of scar is a genuine *source--sink* bottleneck — a wavefront funnelling through it
faces a curvature/load mismatch it cannot in an abstract graph — so reentry can nucleate at
surviving tissue that fibrosis imaging is blind to.

Layout (a rectangular monolayer, fibres along x):
  - left third and right third: healthy open tissue;
  - centre third: SCAR (blocking, high fibrosis) everywhere EXCEPT a narrow horizontal channel of
    width ``isthmus_w`` at mid-height — the healthy isthmus.
A planar wave from the left must funnel through the central isthmus into the right tissue.
Fibrosis peaks in the scar blocks (off the mid-line); the isthmus is low-fibrosis. The graph Fiedler
vector separates left from right across the isthmus, so ``|∇φ₂|`` peaks *at the isthmus* — spatially
dissociated from the fibrosis peak. Whether reentry actually originates at the isthmus is decided by
the monodomain simulator, not assumed; ``meta`` records the isthmus band and the fibrosis peak.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from asb.types import AtrialMesh

__all__ = ["SheetConfig", "make_isthmus_sheet"]


@dataclass
class SheetConfig:
    """Rectangular sheet + central scar-bounded isthmus."""

    lx: float = 60.0        # mm
    ly: float = 60.0        # mm
    spacing: float = 1.0    # mm between grid nodes (~0.2 mm would be ideal; coarser for CPU)
    isthmus_w: float = 6.0  # mm, healthy-channel height at mid-Ly
    scar_frac: float = 0.85 # fibrosis level assigned to scar (blocking)
    scar_soft_mm: float = 2.0  # smoothing width of the scar boundary


def make_isthmus_sheet(seed: int, cfg: SheetConfig) -> AtrialMesh:
    """Triangulated monolayer with a healthy isthmus between two central scar blocks."""
    del seed  # geometry is deterministic; kept for interface symmetry
    nx = int(round(cfg.lx / cfg.spacing)) + 1
    ny = int(round(cfg.ly / cfg.spacing)) + 1
    xs = np.linspace(0.0, cfg.lx, nx)
    ys = np.linspace(0.0, cfg.ly, ny)
    gx, gy = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    n = pts.shape[0]

    def idx(i, j):  # i over y (row), j over x (col)
        return i * nx + j

    faces = []
    for i in range(ny - 1):
        for j in range(nx - 1):
            a, b, c, d = idx(i, j), idx(i, j + 1), idx(i + 1, j), idx(i + 1, j + 1)
            faces.append([a, b, d])
            faces.append([a, d, c])
    faces = np.asarray(faces, np.int64)

    x, y = pts[:, 0], pts[:, 1]
    in_center = (x >= cfg.lx / 3.0) & (x <= 2.0 * cfg.lx / 3.0)
    dist_from_mid = np.abs(y - cfg.ly / 2.0)
    in_channel = dist_from_mid <= (cfg.isthmus_w / 2.0)
    # scar = central third, outside the healthy channel; smooth the channel edge
    raw_scar = in_center & (~in_channel)
    # soft fibrosis: ramp from 0 (channel) to scar_frac (deep scar) over scar_soft_mm
    edge = np.clip((dist_from_mid - cfg.isthmus_w / 2.0) / max(cfg.scar_soft_mm, 1e-6), 0.0, 1.0)
    fibrosis = np.where(in_center, cfg.scar_frac * edge, 0.0)
    fibrosis = np.clip(fibrosis, 0.0, 1.0)

    fibres = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))  # along-x fibres
    uac = np.column_stack([np.clip(x / cfg.lx, 0, 1), np.clip(y / cfg.ly, 0, 1)])
    g = 3
    rb = np.clip((uac[:, 1] * g).astype(np.int64), 0, g - 1)
    cb = np.clip((uac[:, 0] * g).astype(np.int64), 0, g - 1)
    region = (rb * g + cb).astype(np.int64)

    isthmus_mask = in_center & in_channel
    return AtrialMesh(
        points=pts, faces=faces, fibres=fibres, uac=uac,
        fibrosis=fibrosis, region=region, shape_family=f"sheet_iw{cfg.isthmus_w:.0f}",
        meta={"medium": "monodomain_sheet", "isthmus_w": cfg.isthmus_w,
              "isthmus_nodes": np.flatnonzero(isthmus_mask).tolist(),
              "isthmus_xy": [cfg.lx / 2.0, cfg.ly / 2.0],
              "fibrosis_peak": int(np.argmax(fibrosis)),
              "nx": nx, "ny": ny, "n_nodes": int(n)},
    )
