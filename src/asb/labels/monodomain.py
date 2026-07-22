r"""Monodomain Mitchell--Schaeffer reentry-inducibility labeller.

This is a **genuine nonlinear excitable-media simulator**: a monodomain
reaction--diffusion solve with the two-variable Mitchell & Schaeffer (2003)
membrane model on a triangulated atrial surface, with an anisotropic
cotangent-Laplacian diffusion operator whose local conductance is reduced by
fibrosis and whose refractoriness is shortened by fibrosis. It paces an S1--S2 (or
burst) protocol and detects **self-sustaining reentry** to emit a binary
inducibility label + reentry-origin coordinate.

Scope and honesty
------------------
The verdict is a *simulator's* inducibility verdict, exactly like openCARP's, and
is **never** clinical post-operative AF. The label ``source`` is
``'monodomain_ms'`` so it is distinguishable from both the ``'mock_ep'`` eikonal
stand-in and real ``'opencarp'`` labels. This solver is the CPU-feasible
realization of plan §8.4 ("monodomain + phenomenological Mitchell--Schaeffer on
coarsened meshes"); the openCARP runner in :mod:`asb.labels.opencarp` remains the
wired, deferred target for the full Claude-Science sweep.

Why it is not circular with the SFI
-----------------------------------
Reentry here is governed by conduction **wavelength** (CV x ERP), **source--sink**
mismatch and **unidirectional block** in a nonlinear PDE — never by the Fiedler
value :math:`\lambda_2`. So a predictive link from the linear-spectral SFI to this
nonlinear label is a genuine finding, not a tautology.

Model
-----
Dimensionless transmembrane potential :math:`V\in[0,1]`, gating variable
:math:`h\in[0,1]`, time in ms:

.. math::
    \frac{dV}{dt} &= \nabla\!\cdot(D\nabla V) + \frac{h V^2 (1-V)}{\tau_{in}}
                     - \frac{V}{\tau_{out}} + I_{stim}, \\
    \frac{dh}{dt} &= \begin{cases}(1-h)/\tau_{open} & V < V_{gate}\\
                                   -h/\tau_{close}  & V \ge V_{gate}\end{cases}

Diffusion is discretised with the (conductance-weighted, fibrosis-attenuated)
cotangent Laplacian and barycentric lumped mass, integrated by explicit forward
Euler with reaction/diffusion operator splitting. All randomness (pacing-site
choice) comes from the caller's ``rng``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import scipy.sparse as sp

from asb.substrate.mesh import mesh_edges
from asb.types import AtrialMesh, InducibilityLabel

__all__ = [
    "MonodomainConfig",
    "cotangent_operator",
    "measure_planar_cv",
    "single_cell_apd",
    "simulate_monodomain",
    "induce_monodomain",
]


@dataclass
class MonodomainConfig:
    """Monodomain Mitchell--Schaeffer parameters and pacing protocol.

    Defaults are the standard MS(2003) set tuned toward human-atrial APD/ERP, with
    a diffusivity ``d0`` calibrated (see :func:`measure_planar_cv`) so healthy
    cross-fibre CV is order ~0.4 m/s and along-fibre order ~0.7 m/s. Fibrosis both
    reduces local conductance (already in the graph weights) and shortens
    ``tau_close`` (shorter ERP -> shorter wavelength).
    """

    # --- Mitchell-Schaeffer membrane (ms) ---
    tau_in: float = 0.3
    tau_out: float = 6.0
    tau_open: float = 120.0
    tau_close: float = 130.0
    v_gate: float = 0.13
    # Fibrosis shortens tau_close by up to this fraction (ERP shortening).
    fibrosis_erp_shortening: float = 0.5

    # --- Diffusion ---
    # Global diffusivity scale (mm^2/ms). Multiplies the normalized anisotropic
    # cotangent operator; calibrated to physiological CV.
    d0: float = 0.10
    # Reference conductance that maps graph weight -> relative diffusivity 1.
    w_ref: float = 0.3

    # --- Numerics ---
    dt: float = 0.025          # ms
    # --- Pacing protocol ---
    protocol: str = "S1S2"     # 'S1S2' or 'burst'
    n_pacing_sites: int = 3
    n_s1: int = 3
    s1_cycle_length: float = 300.0   # ms
    s2_coupling: float = 200.0       # ms (S1S2)
    burst_cycle_length: float = 130.0
    n_burst: int = 8
    stim_amp: float = 0.6            # added to V at stim nodes
    stim_duration: float = 2.0       # ms
    stim_radius_mm: float = 4.0      # nodes within this radius of a site are paced
    # Observation window after the last stimulus (ms). Must exceed reentry_min_ms.
    observe_after: float = 1000.0
    # Quiescence: max(V) below this for `quiescent_ms` continuous ms => terminated.
    quiescent_v: float = 0.10
    quiescent_ms: float = 40.0
    # Inducible iff self-sustained supra-threshold activity persists at least this
    # long AFTER the last stimulus (>= ~2-3 reentrant rotations). A normal paced
    # response transits + repolarizes and quiesces well inside this window.
    reentry_min_ms: float = 600.0


# --------------------------------------------------------------------------- #
# Mesh geometry: anisotropic cotangent Laplacian + lumped mass.
# --------------------------------------------------------------------------- #
def cotangent_operator(
    mesh: AtrialMesh, cfg: MonodomainConfig,
    *, along: Optional[float] = None, cross: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""Build the conductance-weighted cotangent Laplacian and lumped mass.

    Returns the unique edges, per-edge transmissibility ``T_ij`` (the cotangent
    weight times the relative fibre/fibrosis conductance), the per-vertex lumped
    (barycentric) mass ``M_i`` in mm^2, and per-edge geodesic length in mm.

    The diffusion operator applied to a field ``V`` is
    ``(D V)_i = (d0 / M_i) * sum_j T_ij (V_j - V_i)``.

    Parameters
    ----------
    mesh : AtrialMesh
        Triangulated surface (mm) with ``fibres`` and ``fibrosis``.
    cfg : MonodomainConfig
        Supplies ``w_ref`` (conductance normalization).

    Returns
    -------
    edges : (m, 2) int
    T : (m,) float
        Non-negative per-edge transmissibility.
    mass : (n,) float
        Per-vertex lumped area (mm^2), strictly positive.
    length : (m,) float
        Edge length (mm).
    """
    from asb.graph import edge_weights_from_fibres

    points = np.asarray(mesh.points, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    n = points.shape[0]

    edges = mesh_edges(faces)
    edge_index = {(int(a), int(b)): k for k, (a, b) in enumerate(edges)}
    m = edges.shape[0]

    cot = np.zeros(m, dtype=float)
    mass = np.zeros(n, dtype=float)

    # Per-triangle: accumulate cotangent of the angle opposite each edge, and
    # one-third of the triangle area into each vertex mass (barycentric lumping).
    p0 = points[faces[:, 0]]
    p1 = points[faces[:, 1]]
    p2 = points[faces[:, 2]]
    tri_normal = np.cross(p1 - p0, p2 - p0)
    area = 0.5 * np.linalg.norm(tri_normal, axis=1)
    area = np.maximum(area, 1e-12)
    for k in range(3):
        np.add.at(mass, faces[:, k], area / 3.0)

    # For each of the three edges of every triangle, the opposite vertex's angle.
    # Edge (a,b) opposite vertex c: cot(C) = dot(a-c, b-c)/ (2*area).
    tri = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]  # (a, b, opposite c)
    for a_i, b_i, c_i in tri:
        a = faces[:, a_i]
        b = faces[:, b_i]
        c = faces[:, c_i]
        va = points[a] - points[c]
        vb = points[b] - points[c]
        dot = np.sum(va * vb, axis=1)
        cotC = dot / (2.0 * area)
        for t in range(faces.shape[0]):
            u, v = int(a[t]), int(b[t])
            key = (u, v) if u < v else (v, u)
            idx = edge_index.get(key)
            if idx is not None:
                cot[idx] += 0.5 * cotC[t]

    # Clamp tiny/negative cotangent weights (obtuse triangles) for a stable
    # M-matrix operator.
    cot = np.maximum(cot, 1e-4)

    # Relative fibre/fibrosis conductance per edge (anisotropy + fibrosis floor).
    # ``along``/``cross`` default to the standard 1.0/0.3 anisotropy; they are exposed
    # as overrides for the E6 fibre-field uncertainty sweep.
    w = edge_weights_from_fibres(
        mesh, edges,
        along=1.0 if along is None else float(along),
        cross=0.3 if cross is None else float(cross),
    )
    c_rel = w / float(cfg.w_ref)

    T = cot * c_rel
    length = np.linalg.norm(points[edges[:, 1]] - points[edges[:, 0]], axis=1)
    mass = np.maximum(mass, 1e-9)
    return edges, T, mass, np.maximum(length, 1e-9)


def _diffusion_matrix(
    n: int, edges: np.ndarray, T: np.ndarray, mass: np.ndarray, d0: float
) -> sp.csr_matrix:
    """Sparse operator A with (A V)_i = (d0/M_i) sum_j T_ij (V_j - V_i)."""
    i, j = edges[:, 0], edges[:, 1]
    # Off-diagonal +T, diagonal -sum T, then scale rows by d0/M.
    rows = np.concatenate([i, j, i, j])
    cols = np.concatenate([j, i, i, j])
    data = np.concatenate([T, T, -T, -T])
    A = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    scale = d0 / mass
    return sp.diags(scale) @ A


# --------------------------------------------------------------------------- #
# Core reaction-diffusion time-stepper.
# --------------------------------------------------------------------------- #
def _site_mask(mesh: AtrialMesh, sites: np.ndarray, radius: float) -> np.ndarray:
    """Boolean node mask: within ``radius`` mm of any of ``sites``."""
    points = np.asarray(mesh.points, dtype=float)
    n = points.shape[0]
    mask = np.zeros(n, dtype=bool)
    for s in np.atleast_1d(sites).astype(int):
        mask |= np.linalg.norm(points - points[s], axis=1) <= radius
    if not mask.any():
        mask[int(np.atleast_1d(sites).astype(int)[0])] = True
    return mask


def build_stimuli(
    mesh: AtrialMesh,
    cfg: MonodomainConfig,
    stim_sites: np.ndarray,
    *,
    s2_coupling: Optional[float] = None,
) -> Tuple[list, float]:
    r"""Build the stimulus schedule ``[(onset_ms, node_mask), ...]`` and total ms.

    Three protocols:

    - ``'burst'``: ``n_burst`` beats at ``burst_cycle_length`` from the site mask.
    - ``'S1S2'``: ``n_s1`` S1 beats then a premature S2, all from the site mask.
    - ``'crossfield'``: an S1 conditioning wave from a **line** at low UAC-alpha,
      then a spatially offset S2 covering the low-UAC-beta half-field at the
      coupling interval — the standard cross-field spiral-wave inducer. The S2
      blocks in still-refractory tissue and propagates elsewhere, seeding a rotor
      wherever the wavelength is short enough (fibrotic, heterogeneous substrate).
    """
    ci = cfg.s2_coupling if s2_coupling is None else float(s2_coupling)
    sites = np.atleast_1d(stim_sites).astype(int)
    r = cfg.stim_radius_mm

    if cfg.protocol == "burst":
        mask = _site_mask(mesh, sites, r)
        stimuli = [(k * cfg.burst_cycle_length, mask) for k in range(cfg.n_burst)]
    elif cfg.protocol == "crossfield":
        uac = np.asarray(mesh.uac, dtype=float)
        # S1 conditioning line at low alpha; wave sweeps in +alpha.
        s1_mask = uac[:, 0] <= 0.12
        if not s1_mask.any():
            s1_mask = _site_mask(mesh, sites, r)
        # S2 half-field at low beta, delivered after the S1 wave.
        s2_mask = uac[:, 1] <= 0.5
        if not s2_mask.any():
            s2_mask = ~s1_mask
        stimuli = [(0.0, s1_mask), (ci, s2_mask)]
    else:  # 'S1S2'
        mask = _site_mask(mesh, sites, r)
        stimuli = [(k * cfg.s1_cycle_length, mask) for k in range(cfg.n_s1)]
        stimuli.append(((cfg.n_s1 - 1) * cfg.s1_cycle_length + ci, mask))

    last = max(on for on, _ in stimuli)
    total = float(last + cfg.stim_duration + cfg.observe_after)
    return stimuli, total


def _react_step(
    V: np.ndarray, h: np.ndarray, cfg: MonodomainConfig,
    tau_close: np.ndarray, dt: float, stim: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """One explicit reaction (+stim) update of (V, h)."""
    j_in = h * V * V * (1.0 - V) / cfg.tau_in
    j_out = -V / cfg.tau_out
    V = V + dt * (j_in + j_out + stim)
    # Gating: open when V<v_gate else close (per-node tau_close from fibrosis).
    below = V < cfg.v_gate
    dh = np.where(below, (1.0 - h) / cfg.tau_open, -h / tau_close)
    h = h + dt * dh
    V = np.clip(V, 0.0, 1.5)
    h = np.clip(h, 0.0, 1.0)
    return V, h


def simulate_monodomain(
    mesh: AtrialMesh,
    cfg: MonodomainConfig,
    stim_sites: np.ndarray,
    *,
    record_activation: bool = True,
    s2_coupling: Optional[float] = None,
    along: Optional[float] = None,
    cross: Optional[float] = None,
) -> Dict[str, object]:
    r"""Run one monodomain Mitchell--Schaeffer pacing simulation.

    Parameters
    ----------
    mesh : AtrialMesh
        Triangulated atrial surface (mm) with fibres + fibrosis.
    cfg : MonodomainConfig
        Model + protocol parameters.
    stim_sites : (s,) int
        Vertex indices of the pacing site centres. All vertices within
        ``cfg.stim_radius_mm`` of a site receive the stimulus current.
    record_activation : bool
        If True, record the last-activation time per node and a max(V) trace.

    Returns
    -------
    dict
        ``inducible`` (bool), ``reentry_origin`` (int|None), ``activation`` (n,),
        ``maxV_trace`` (T,), ``t`` (T,), ``terminated_at`` (float|None),
        ``n_active_end`` (int), plus diagnostic scalars.
    """
    from scipy.sparse.linalg import factorized

    points = np.asarray(mesh.points, dtype=float)
    n = points.shape[0]
    edges, T, mass, _length = cotangent_operator(mesh, cfg, along=along, cross=cross)
    A = _diffusion_matrix(n, edges, T, mass, cfg.d0)
    # Implicit (backward-Euler) diffusion: unconditionally stable, so dt is limited
    # only by reaction accuracy, not the mesh CFL. Factorize (I - dt*A) once.
    identity = sp.identity(n, format="csc", dtype=float)
    B = (identity - cfg.dt * A).tocsc()
    diffuse = factorized(B)

    # Per-node refractoriness: fibrosis shortens tau_close.
    fib = np.clip(np.asarray(mesh.fibrosis, dtype=float), 0.0, 1.0)
    tau_close = cfg.tau_close * (1.0 - cfg.fibrosis_erp_shortening * fib)
    tau_close = np.maximum(tau_close, 1e-3)

    # Stimulus schedule: list of (onset_ms, node_mask).
    stimuli, total = build_stimuli(mesh, cfg, stim_sites, s2_coupling=s2_coupling)
    onsets = np.array([on for on, _ in stimuli], dtype=float)
    masks = [mk for _, mk in stimuli]
    dt = cfg.dt
    n_steps = int(np.ceil(total / dt))

    V = np.zeros(n, dtype=float)
    h = np.ones(n, dtype=float)
    activation = np.full(n, -1.0, dtype=float)
    prevV = V.copy()

    maxV_trace = np.empty(n_steps, dtype=float)
    t_trace = np.empty(n_steps, dtype=float)

    stim_amp = cfg.stim_amp
    stim_dur = cfg.stim_duration
    last_onset = float(onsets[-1])
    post_stim_t = last_onset + stim_dur
    quiescent_run = 0.0
    terminated_at: Optional[float] = None
    last_active_t = post_stim_t  # last time supra-threshold activity was seen post-stim

    for step in range(n_steps):
        t = step * dt
        # Reaction (+stim) explicit half, then implicit diffusion (Godunov split).
        stim = np.zeros(n, dtype=float)
        active = (t >= onsets) & (t < onsets + stim_dur)
        if active.any():
            for a in np.flatnonzero(active):
                stim[masks[a]] = stim_amp
        V, h = _react_step(V, h, cfg, tau_close, dt, stim)
        V = diffuse(V)

        # Activation crossing (upstroke through v_gate).
        crossed = (prevV < cfg.v_gate) & (V >= cfg.v_gate)
        activation[crossed] = t
        prevV = V

        mv = float(V.max())
        maxV_trace[step] = mv
        t_trace[step] = t

        # After the last stimulus, track self-sustained activity and quiescence.
        if t > post_stim_t:
            if mv < cfg.quiescent_v:
                quiescent_run += dt
                if quiescent_run >= cfg.quiescent_ms and terminated_at is None:
                    terminated_at = t
                    break
            else:
                quiescent_run = 0.0
                last_active_t = t

    # Duration of self-sustained supra-threshold activity after pacing ended.
    sustained_ms = float(max(0.0, last_active_t - post_stim_t))
    n_active_end = int((V >= cfg.v_gate).sum())
    # Inducible iff a reentrant response persisted >= reentry_min_ms past the last
    # stimulus (a normal paced response transits + repolarizes well inside it).
    inducible = bool(sustained_ms >= cfg.reentry_min_ms)

    # Reentry origin: among nodes activated in the post-stim observation window, the
    # earliest one (the leading edge of the sustained circuit).
    origin: Optional[int] = None
    if inducible:
        post = activation > post_stim_t
        if post.any():
            cand = np.where(post, activation, np.inf)
            origin = int(np.argmin(cand))
        else:
            origin = int(np.argmax(V))

    out: Dict[str, object] = {
        "inducible": inducible,
        "reentry_origin": origin,
        "sustained_ms": sustained_ms,
        "terminated_at": terminated_at,
        "n_active_end": n_active_end,
        "last_onset": last_onset,
        "total": total,
    }
    if record_activation:
        used = step + 1
        out["activation"] = activation
        out["maxV_trace"] = maxV_trace[:used]
        out["t"] = t_trace[:used]
    return out


def induce_monodomain(
    mesh: AtrialMesh, cfg: MonodomainConfig, rng: np.random.Generator,
    *, burst_cls: Optional[Tuple[float, ...]] = None,
    along: Optional[float] = None, cross: Optional[float] = None,
) -> InducibilityLabel:
    """Monodomain-MS inducibility verdict for one atrium (``source='monodomain_ms'``).

    Runs a small, fixed **burst-pacing induction battery** and flags the atrium
    inducible if *any* attempt provokes self-sustained reentry (activity persisting
    ``>= reentry_min_ms`` after pacing ends). Burst pacing is the standard aggressive
    inducer and, unlike the half-field cross-field protocol (which anchors geometric
    reentry independent of substrate), it discriminates on refractoriness/wavelength:
    pacing faster than healthy ERP but slower than short fibrotic ERP breaks waves
    only where refractoriness is short and heterogeneous. The battery scans a couple
    of burst cycle lengths from ``n_pacing_sites`` random sites; the reentry origin is
    taken from the first inducing attempt. **Simulator verdict, never clinical POAF.**

    Parameters
    ----------
    mesh : AtrialMesh
        Triangulated atrial surface (mm).
    cfg : MonodomainConfig
        Model parameters (``protocol`` is overridden to ``'burst'`` per attempt).
    rng : numpy.random.Generator
        Source of randomness for burst pacing-site selection.
    burst_cls : tuple of float, optional
        Burst cycle lengths (ms) to scan. Defaults to ``(130, 160)``.

    Returns
    -------
    InducibilityLabel
        ``source='monodomain_ms'``.
    """
    from dataclasses import replace as _replace

    n = mesh.n_points
    if burst_cls is None:
        burst_cls = (130.0, 160.0)
    n_sites = max(1, min(int(cfg.n_pacing_sites), n))
    sites = rng.choice(n, size=n_sites, replace=False)

    attempts: list = []
    for cl in burst_cls:
        bcfg = _replace(cfg, protocol="burst", burst_cycle_length=float(cl))
        for s in sites:
            attempts.append(("burst", bcfg, np.array([int(s)]), None))

    any_ind = False
    origin: Optional[int] = None
    best_sustained = 0.0
    detail: list = []
    for kind, acfg, site, ci in attempts:
        res = simulate_monodomain(
            mesh, acfg, site, record_activation=False, s2_coupling=ci,
            along=along, cross=cross,
        )
        sm = float(res["sustained_ms"])
        best_sustained = max(best_sustained, sm)
        detail.append({"kind": kind, "ci": ci, "site": int(site[0]),
                       "inducible": bool(res["inducible"]), "sustained_ms": sm})
        if res["inducible"] and not any_ind:
            any_ind = True
            origin = res["reentry_origin"]

    meta: Dict[str, object] = {
        "note": "monodomain Mitchell-Schaeffer simulator verdict; NOT clinical POAF",
        "battery": detail,
        "best_sustained_ms": best_sustained,
        "sites": sites.tolist(),
    }
    return InducibilityLabel(
        inducible=bool(any_ind),
        reentry_origin=origin,
        protocol="battery",
        source="monodomain_ms",
        meta=meta,
    )


# --------------------------------------------------------------------------- #
# Verification helpers (planar CV, single-cell APD/ERP) — used by E0 gates.
# --------------------------------------------------------------------------- #
def single_cell_apd(
    cfg: MonodomainConfig, *, tau_close: Optional[float] = None,
    apd_frac: float = 0.9, cycle_length: float = 1000.0, n_beats: int = 2,
) -> Dict[str, float]:
    """Single-cell APD_90 (ms) for the MS model (diffusion off).

    Paces an isolated cell and measures the action-potential duration at
    ``apd_frac`` repolarization on the last beat.
    """
    tc = cfg.tau_close if tau_close is None else tau_close
    dt = cfg.dt
    V = 0.0
    h = 1.0
    total = cycle_length * n_beats
    n_steps = int(total / dt)
    apd_start = None
    apd = float("nan")
    thresh = cfg.v_gate
    last_beat_start = (n_beats - 1) * cycle_length
    for step in range(n_steps):
        t = step * dt
        stim = cfg.stim_amp if (t % cycle_length) < cfg.stim_duration else 0.0
        j_in = h * V * V * (1.0 - V) / cfg.tau_in
        j_out = -V / cfg.tau_out
        Vn = V + dt * (j_in + j_out + stim)
        below = Vn < cfg.v_gate
        dh = (1.0 - h) / cfg.tau_open if below else -h / tc
        h = min(max(h + dt * dh, 0.0), 1.0)
        # APD on the last beat: from upstroke to (1-apd_frac) repolarization.
        if t >= last_beat_start:
            vmax_ref = 1.0
            if V < thresh <= Vn and apd_start is None:
                apd_start = t
            if apd_start is not None and Vn <= (1.0 - apd_frac) * vmax_ref and np.isnan(apd):
                # crude: repolarization below (1-frac) of peak amplitude
                if Vn < V:  # repolarizing
                    apd = t - apd_start
        V = min(max(Vn, 0.0), 1.5)
    return {"apd90": float(apd), "tau_close": float(tc)}


def measure_planar_cv(
    cfg: MonodomainConfig, *, length_mm: float = 40.0, dx_mm: float = 0.5,
    width: int = 3,
) -> Dict[str, float]:
    """Planar-wave conduction velocity (m/s) on a regular strip (calibration).

    Builds a thin regular triangulated strip of isotropic healthy tissue, paces
    one end, and measures CV from activation times between two probe positions.
    Returns CV in m/s (== mm/ms).
    """
    nx = int(round(length_mm / dx_mm)) + 1
    ny = int(width)
    xs = np.arange(nx) * dx_mm
    ys = np.arange(ny) * dx_mm
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), np.zeros(nx * ny)], axis=1)

    def vid(ix, iy):
        return ix * ny + iy

    faces = []
    for ix in range(nx - 1):
        for iy in range(ny - 1):
            a, b, c, d = vid(ix, iy), vid(ix + 1, iy), vid(ix + 1, iy + 1), vid(ix, iy + 1)
            faces.append([a, b, c])
            faces.append([a, c, d])
    faces = np.asarray(faces, dtype=np.int64)
    n = pts.shape[0]
    strip = AtrialMesh(
        points=pts, faces=faces,
        fibres=np.tile(np.array([1.0, 0.0, 0.0]), (n, 1)),
        uac=np.zeros((n, 2)), fibrosis=np.zeros(n), region=np.zeros(n, dtype=np.int64),
        shape_family="strip",
    )
    # One S1 beat from the x=0 end; short observation.
    scfg = MonodomainConfig(
        **{**cfg.__dict__, "protocol": "S1S2", "n_s1": 1, "n_pacing_sites": 1,
           "observe_after": 120.0, "s1_cycle_length": 400.0, "stim_radius_mm": dx_mm * 1.5}
    )
    site = int(vid(0, ny // 2))
    res = simulate_monodomain(strip, scfg, np.array([site]), record_activation=True)
    act = np.asarray(res["activation"], dtype=float)
    # Probe at x=10mm and x=30mm on the mid-row.
    x1, x2 = 10.0, 30.0
    i1 = int(vid(int(round(x1 / dx_mm)), ny // 2))
    i2 = int(vid(int(round(x2 / dx_mm)), ny // 2))
    t1, t2 = act[i1], act[i2]
    cv = float("nan")
    if t1 >= 0 and t2 > t1:
        cv = (x2 - x1) / (t2 - t1)  # mm/ms == m/s
    return {"cv_m_per_s": cv, "t1": float(t1), "t2": float(t2)}
