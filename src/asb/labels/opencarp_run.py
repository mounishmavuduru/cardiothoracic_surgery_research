r"""Label a mesh with openCARP -- the pre-registered ground-truth solver.

`asb.labels.opencarp` is the interface stub written when openCARP could not be installed.
This module is the working runner. It exports a mesh, builds a parameter file, invokes the
solver and reads an inducibility verdict back out.

Three details are easy to get wrong and are handled explicitly:

**CALIBRATION STATUS: see `results/opencarp_calibration.json`.**
The plumbing is verified end to end -- a 2177-node atrial mesh labels in ~14 s with real
propagation (up to 1028 of 2177 nodes active) -- but a verdict of "inducible" is not yet
trustworthy. openCARP's default MitchellSchaeffer gives APD90 = 246 ms while the burst
cycle length inherited from the monodomain protocol is 150 ms. Pacing faster than the
action-potential duration produces persistent depolarisation and rate-dependent block,
which :func:`detect_reentry` cannot distinguish from reentry: the first smoke test
returned ``sustained_ms = 1000.0``, i.e. exactly the full observation window, which is
what a tissue that never repolarises looks like. Before any label is used, the
pre-registration's calibration gate must be re-passed on this solver -- inducible fraction
in the 10-40 % band and a positive fibrosis/reentry rank correlation -- with the membrane
parameters or the cycle length adjusted so the two are mutually consistent.

*Units.* openCARP `.pts` files are in MICROMETRES. :class:`~asb.types.AtrialMesh` carries
millimetres (`asb.substrate.roney` and `asb.substrate.uw_boyle` both divide by 1000 on
load). Writing millimetres into a `.pts` shrinks the atrium a thousandfold, which does not
error -- it silently produces a mesh on which nothing propagates. :data:`MM_TO_MICRON`
undoes it.

*Fibrosis discretisation.* openCARP assigns conductivity through discrete `gregion` element
tags, so a continuous fibrosis field must be binned. `results/roney_fibrosis_quantization.json`
shows that binarising fibrosis roughly DOUBLES inducibility on this pipeline (20/62 -> 40/62
at matched burden), so a coarse binning would import exactly that artefact. The default is
therefore deliberately fine, and :func:`bin_sensitivity` exists so a study can demonstrate
insensitivity to the bin count before trusting any label.

*Line endings.* Every file the solver consumes is written with an explicit LF newline. On
Windows, Python text mode emits CRLF; openCARP then reads the ionic-model name with a
trailing carriage return and reports ``Illegal IM specified: MitchellSchaeffer`` -- an
error that names the model and therefore points away from the real cause.

*Membrane model.* openCARP ships `MitchellSchaeffer`, the same model as
`asb.labels.monodomain`, so switching solvers is not also a change of physics. Note it is a
NORMALISED model (Vm ~ 0..1): a physiological 250 uA/cm^2 stimulus makes the parabolic solve
diverge with NaN. See :data:`DEFAULT_STIM_STRENGTH`.
"""
from __future__ import annotations

import os
import re
import shutil
import struct
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from asb.types import AtrialMesh, InducibilityLabel

__all__ = [
    "MM_TO_MICRON",
    "DEFAULT_STIM_STRENGTH",
    "OpenCARPConfig",
    "find_opencarp",
    "export_for_opencarp",
    "build_par",
    "read_igb",
    "detect_reentry",
    "label_with_opencarp",
]

#: openCARP point files are in micrometres; AtrialMesh carries millimetres.
MM_TO_MICRON: float = 1.0e3

#: MitchellSchaeffer is normalised (Vm ~ 0..1). 250 uA/cm^2 diverges; 60 propagates.
DEFAULT_STIM_STRENGTH: float = 60.0

_CARP_NAMES = ("openCARP", "carp.pt", "carp.petsc", "carp")


@dataclass(frozen=True)
class OpenCARPConfig:
    """Frozen openCARP protocol. Mirrors ``FROZEN_MONO`` where the two overlap."""

    imp: str = "MitchellSchaeffer"
    #: Longitudinal / transverse intracellular conductivity of HEALTHY tissue (S/m).
    g_il: float = 0.174
    g_it: float = 0.019
    #: Fibrosis attenuates conductivity by (1 - f), floored so scar still conducts a little.
    fibrosis_floor: float = 0.05
    #: Number of conductivity bins the continuous fibrosis field is discretised into.
    #: Deliberately fine -- see the module docstring on the quantization artefact.
    n_bins: int = 20
    #: Burst pacing, matching the monodomain protocol.
    n_pacing_sites: int = 2
    n_burst: int = 6
    burst_cl_ms: float = 150.0
    stim_strength: float = DEFAULT_STIM_STRENGTH
    stim_duration_ms: float = 2.0
    stim_radius_um: float = 4000.0
    #: Observation window after the last stimulus, and the sustained-activity threshold.
    observe_after_ms: float = 1000.0
    reentry_min_ms: float = 650.0
    #: Solver timestep in MICROseconds, and output sampling in milliseconds.
    dt_us: int = 25
    out_dt_ms: float = 5.0
    #: MitchellSchaeffer membrane parameters, matched to ``asb.labels.monodomain``'s
    #: frozen values so a switch of solver is not also a change of physics. openCARP's
    #: own defaults are tau_close 150 / tau_out 5, which give APD90 = 277 ms.
    #: Calibrated on this build: tau_close 110 -> 209 ms, 55 -> 113 ms, so
    #: tau_close = 115 x (1 - 0.5 f) reproduces the in-house 218 ms healthy / 121 ms
    #: fibrotic within about 2 %.
    ms_tau_close: float = 115.0
    ms_tau_out: float = 6.0
    ms_tau_in: float = 0.3
    ms_tau_open: float = 120.0
    ms_v_gate: float = 0.13
    #: Fibrosis shortens the action potential as well as slowing conduction. Omitting
    #: this was the calibration failure: with ERP fixed, more fibrosis produced only more
    #: block, giving 90 % "inducible" and a NEGATIVE fibrosis-verdict correlation.
    fibrosis_erp_shortening: float = 0.5
    #: Reentry-vs-block discrimination (see :func:`detect_reentry`). The tissue must
    #: repolarise below this depolarised fraction at some point in the window, and at
    #: least this share of nodes must cross threshold upward twice or more.
    max_depol_fraction: float = 0.95
    min_reactivating_nodes: float = 0.05


def find_opencarp(explicit: Optional[str] = None) -> Optional[str]:
    """Resolve an openCARP binary, preferring an extracted AppImage under ``$HOME``."""
    if explicit and os.path.isfile(explicit):
        return explicit
    env = os.environ.get("OPENCARP_BIN")
    if env and os.path.isfile(env):
        return env
    for name in _CARP_NAMES:
        found = shutil.which(name)
        if found:
            return found
    home = os.path.expanduser("~")
    import glob
    for pat in ("opencarp/*/squashfs-root/usr/bin/openCARP",
                "opencarp/*/*/squashfs-root/usr/bin/openCARP"):
        hits = sorted(glob.glob(os.path.join(home, pat)))
        if hits:
            return hits[0]
    # Windows host, solver inside WSL. Returned with a "wsl:" marker so the runner
    # knows to bridge; a bare Linux path is not executable by Windows subprocess.
    if shutil.which("wsl"):
        probe = subprocess.run(
            ["wsl", "-e", "bash", "-lc",
             "ls -d ~/opencarp/*/squashfs-root/usr/bin/openCARP 2>/dev/null | head -1"],
            capture_output=True, text=True)
        cand = (probe.stdout or "").strip().splitlines()
        if cand and cand[0]:
            return "wsl:" + cand[0]
    return None


def _to_wsl_path(win_path: str) -> str:
    r"""``C:\Users\x`` -> ``/mnt/c/Users/x``."""
    p = os.path.abspath(win_path).replace("\\", "/")
    if len(p) > 1 and p[1] == ":":
        return "/mnt/" + p[0].lower() + p[2:]
    return p


def _fibrosis_bins(fib: np.ndarray, n_bins: int) -> Tuple[np.ndarray, np.ndarray]:
    """Map fibrosis in [0,1] to ``1..n_bins`` and return (bin index, bin centre)."""
    f = np.clip(np.asarray(fib, dtype=float), 0.0, 1.0)
    idx = np.clip((f * n_bins).astype(np.int64), 0, n_bins - 1)
    centres = (np.arange(n_bins) + 0.5) / n_bins
    return idx + 1, centres


def export_for_opencarp(mesh: AtrialMesh, out_dir: str, cfg: OpenCARPConfig) -> Dict:
    """Write ``mesh.pts`` / ``.elem`` / ``.lon`` in openCARP units and tag by fibrosis bin."""
    os.makedirs(out_dir, exist_ok=True)
    pts = np.asarray(mesh.points, dtype=float) * MM_TO_MICRON     # mm -> um
    faces = np.asarray(mesh.faces, dtype=np.int64)
    fib = np.asarray(mesh.fibrosis, dtype=float)
    fibres = np.asarray(mesh.fibres, dtype=float)

    vbin, centres = _fibrosis_bins(fib, cfg.n_bins)
    # An element's tag is the bin of its mean vertex fibrosis.
    elem_f = fib[faces].mean(axis=1)
    ebin, _ = _fibrosis_bins(elem_f, cfg.n_bins)

    base = os.path.join(out_dir, "mesh")
    with open(base + ".pts", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"{pts.shape[0]}\n")
        fh.writelines(f"{x:.6f} {y:.6f} {z:.6f}\n" for x, y, z in pts)
    with open(base + ".elem", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"{faces.shape[0]}\n")
        fh.writelines(f"Tr {a} {b} {c} {t}\n"
                      for (a, b, c), t in zip(faces, ebin))
    with open(base + ".lon", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("1\n")
        for tri in faces:
            v = fibres[tri].mean(axis=0)
            n = float(np.linalg.norm(v))
            v = v / n if n > 0 else np.array([1.0, 0.0, 0.0])
            fh.write(f"{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")

    return {"basename": base, "n_points": int(pts.shape[0]),
            "n_elems": int(faces.shape[0]),
            "tags_present": sorted(set(int(t) for t in ebin)),
            "bin_centres": centres.tolist(), "vertex_bins": vbin}


def _stim_block(i: int, centre_um: np.ndarray, cfg: OpenCARPConfig, start_ms: float) -> str:
    r = cfg.stim_radius_um
    x, y, z = centre_um
    return "\n".join([
        f"stimulus[{i}].name = \"burst{i}\"",
        f"stimulus[{i}].stimtype = 0",
        f"stimulus[{i}].strength = {cfg.stim_strength}",
        f"stimulus[{i}].duration = {cfg.stim_duration_ms}",
        f"stimulus[{i}].start = {start_ms}",
        f"stimulus[{i}].npls = {cfg.n_burst}",
        f"stimulus[{i}].bcl = {cfg.burst_cl_ms}",
        f"stimulus[{i}].x0 = {x - r:.1f}", f"stimulus[{i}].xd = {2 * r:.1f}",
        f"stimulus[{i}].y0 = {y - r:.1f}", f"stimulus[{i}].yd = {2 * r:.1f}",
        f"stimulus[{i}].z0 = {z - r:.1f}", f"stimulus[{i}].zd = {2 * r:.1f}",
    ])


def build_par(export: Dict, cfg: OpenCARPConfig, stim_centres_um: Sequence[np.ndarray]) -> str:
    """Assemble the parameter file: one gregion per occupied fibrosis bin."""
    tags = export["tags_present"]
    centres = export["bin_centres"]
    lines: List[str] = [
        "num_phys_regions = 1",
        'phys_region[0].name = "intra"',
        "phys_region[0].ptype = 0",
        f"phys_region[0].num_IDs = {len(tags)}",
    ]
    lines += [f"phys_region[0].ID[{k}] = {t}" for k, t in enumerate(tags)]

    lines += ["", f"num_imp_regions = {len(tags)}"]
    for k, t in enumerate(tags):
        f = centres[t - 1]
        # Fibrosis shortens the action potential, exactly as fibrosis_erp_shortening does
        # in the monodomain solver. Without this the two solvers are not comparable.
        tau_close = cfg.ms_tau_close * (1.0 - cfg.fibrosis_erp_shortening * f)
        params = (f"tau_close={tau_close:.4f},tau_out={cfg.ms_tau_out},"
                  f"tau_in={cfg.ms_tau_in},tau_open={cfg.ms_tau_open},"
                  f"V_gate={cfg.ms_v_gate}")
        lines += [f"imp_region[{k}].im = {cfg.imp}",
                  f'imp_region[{k}].im_param = "{params}"',
                  f"imp_region[{k}].num_IDs = 1",
                  f"imp_region[{k}].ID[0] = {t}"]

    lines += ["", f"num_gregions = {len(tags)}"]
    for k, t in enumerate(tags):
        f = centres[t - 1]
        scale = max(1.0 - f, cfg.fibrosis_floor)
        lines += [f"gregion[{k}].num_IDs = 1",
                  f"gregion[{k}].ID[0] = {t}",
                  f"gregion[{k}].g_il = {cfg.g_il * scale:.6f}",
                  f"gregion[{k}].g_it = {cfg.g_it * scale:.6f}",
                  f"gregion[{k}].g_in = {cfg.g_it * scale:.6f}"]

    lines += ["", f"num_stim = {len(stim_centres_um)}"]
    for i, c in enumerate(stim_centres_um):
        lines.append(_stim_block(i, c, cfg, start_ms=1.0))

    last_stim = 1.0 + cfg.n_burst * cfg.burst_cl_ms
    lines += ["", "bidomain = 0",
              f"tend = {last_stim + cfg.observe_after_ms}",
              f"dt = {cfg.dt_us}",
              f"spacedt = {cfg.out_dt_ms}",
              f"timedt = {cfg.out_dt_ms}"]
    # openCARP's parameter reader rejects blank lines ("Unrecognized keyword"), so the
    # spacers used above for readability must be stripped before writing.
    return "\n".join(l for l in lines if l.strip()) + "\n"


def read_igb(path: str) -> np.ndarray:
    """Read an openCARP ``.igb`` into ``(n_frames, n_nodes)`` float32."""
    with open(path, "rb") as fh:
        raw = fh.read()
    head = raw[:1024].decode("latin-1")
    x = int(re.search(r"x:(\d+)", head).group(1))
    t = int(re.search(r"t:(\d+)", head).group(1))
    n = x * t
    vals = struct.unpack("<%df" % n, raw[1024:1024 + n * 4])
    return np.asarray(vals, dtype=np.float32).reshape(t, x)


def detect_reentry(vm: np.ndarray, cfg: OpenCARPConfig, thresh: float = 0.5) -> Dict:
    """Verdict on the post-pacing window, separating reentry from rate-dependent block.

    Sustained supra-threshold activity alone is not sufficient. Because the burst cycle
    length (150 ms) is shorter than the action-potential duration (openCARP's default
    MitchellSchaeffer gives APD90 = 246 ms), tissue can simply fail to repolarise, which
    a duration-only criterion scores as inducible. Three conditions are therefore required
    together:

    1. activity persists for at least ``reentry_min_ms`` after the last stimulus;
    2. the tissue REPOLARISES -- the depolarised fraction drops below
       ``max_depol_fraction`` at some point, so a wavefront is moving rather than the
       sheet sitting depolarised;
    3. tissue REACTIVATES -- a non-trivial share of nodes cross threshold upward at least
       twice, which a single decaying wave cannot produce.

    Condition 1 alone reproduces the monodomain criterion; 2 and 3 are what the monodomain
    labeller did not need and this solver does.
    """
    n_frames = vm.shape[0]
    dt = cfg.out_dt_ms
    last_stim_ms = 1.0 + cfg.n_burst * cfg.burst_cl_ms
    first = int(np.ceil(last_stim_ms / dt))
    if first >= n_frames:
        return {"inducible": False, "sustained_ms": 0.0,
                "reason": "observation window shorter than the pacing train"}

    post = vm[first:]
    above = post > thresh
    active = above.sum(axis=1)

    best = run = 0
    for a in active:
        run = run + 1 if a > 0 else 0
        best = max(best, run)
    sustained = float(best * dt)

    depol_frac = above.mean(axis=1)
    min_depol = float(depol_frac.min()) if depol_frac.size else 0.0
    repolarises = min_depol < cfg.max_depol_fraction

    ups = ((~above[:-1]) & above[1:]).sum(axis=0) if above.shape[0] > 1 else np.zeros(1)
    n_reactivating = int((ups >= 2).sum())
    reactivates = n_reactivating >= cfg.min_reactivating_nodes * vm.shape[1]

    inducible = bool(sustained >= cfg.reentry_min_ms and repolarises and reactivates)
    return {"inducible": inducible,
            "sustained_ms": sustained,
            "max_active_nodes": int(active.max()) if active.size else 0,
            "min_depolarised_fraction": min_depol,
            "repolarises": bool(repolarises),
            "n_nodes_reactivating": n_reactivating,
            "reactivates": bool(reactivates)}


def label_with_opencarp(
    mesh: AtrialMesh, cfg: OpenCARPConfig, rng: np.random.Generator,
    *, opencarp_bin: Optional[str] = None, workdir: Optional[str] = None,
    keep: bool = False, timeout_s: int = 3600,
) -> InducibilityLabel:
    """Export, run and score one mesh. Raises if the solver is unavailable or fails."""
    binary = find_opencarp(opencarp_bin)
    if binary is None:
        from asb.labels.opencarp import OpenCARPUnavailableError, _INSTALL_POINTER
        raise OpenCARPUnavailableError(_INSTALL_POINTER)

    tmp = workdir or tempfile.mkdtemp(prefix="asb_carp_")
    export = export_for_opencarp(mesh, tmp, cfg)
    pts_um = np.asarray(mesh.points, dtype=float) * MM_TO_MICRON
    sites = rng.choice(pts_um.shape[0], size=min(cfg.n_pacing_sites, pts_um.shape[0]),
                       replace=False)
    par = build_par(export, cfg, [pts_um[s] for s in sites])
    par_path = os.path.join(tmp, "run.par")
    with open(par_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(par)

    sim = os.path.join(tmp, "out")
    if binary.startswith("wsl:"):
        lin = binary[4:]
        lib = os.path.dirname(os.path.dirname(lin)) + "/lib"
        cmd = (f'export LD_LIBRARY_PATH="{lib}:$LD_LIBRARY_PATH"; '
               f'"{lin}" +F "{_to_wsl_path(par_path)}" '
               f'-meshname "{_to_wsl_path(export["basename"])}" '
               f'-simID "{_to_wsl_path(sim)}"')
        proc = subprocess.run(["wsl", "-e", "bash", "-lc", cmd],
                              capture_output=True, text=True, timeout=timeout_s)
    else:
        env = dict(os.environ)
        lib = os.path.normpath(os.path.join(os.path.dirname(binary), "..", "lib"))
        env["LD_LIBRARY_PATH"] = lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
        proc = subprocess.run([binary, "+F", par_path, "-meshname", export["basename"],
                               "-simID", sim],
                              capture_output=True, text=True, timeout=timeout_s, env=env)
    igb = os.path.join(sim, "vm.igb")
    if not os.path.isfile(igb):
        raise RuntimeError(f"openCARP produced no vm.igb (rc={proc.returncode})\n"
                           f"{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}")
    verdict = detect_reentry(read_igb(igb), cfg)
    if not keep and workdir is None:
        shutil.rmtree(tmp, ignore_errors=True)
    return InducibilityLabel(
        inducible=verdict["inducible"],
        reentry_origin=None,
        protocol="burst",
        source="opencarp",
        meta={**verdict, "solver": binary, "imp": cfg.imp, "n_bins": cfg.n_bins,
              "pacing_sites": [int(s) for s in sites]},
    )
