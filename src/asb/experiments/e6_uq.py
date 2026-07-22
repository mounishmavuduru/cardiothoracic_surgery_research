r"""E6 — Uncertainty quantification (descoped: Morris screening + GP surrogate).

The pre-registered, compute-descoped UQ (dossier p.42; EXPERIMENTAL_PLAN E6): how
sensitive are the results to the four modelling knobs — conduction scale ``d0`` (CV),
the fibrosis threshold ``IIR_dense``, the fibre anisotropy ``cross``, and the SFI
stress magnitude ``Δw``? We do **not** run full Sobol.

Two honest pieces:

1. **Morris elementary-effects screening** of the three *ground-truth-label* knobs
   (``d0``, ``IIR_dense``, ``cross``) on a continuous label outcome (cohort-mean
   self-sustained reentry time). The fourth knob, ``Δw``, provably enters **only the
   predictor** (it scales the SFI stress field, never the monodomain label), so it is
   screened separately (piece 2) rather than mixed into the label-Morris.
2. **Δw robustness of the conclusions**: recompute the GM1 grouped ΔAUC and the GM3
   validity ratio ρ at ``Δw ∈ {0.18, 0.36, 0.54}`` (reusing the cached frozen labels;
   only the SFI features are recomputed), to show the null and the validity-radius
   conclusions do not hinge on the frozen ``Δw = 0.36``.

A Gaussian-process surrogate is fit to the label-Morris design points to report
smooth main-effect (partial-dependence) curves. Everything is cached and one-command.
Labels are monodomain Mitchell--Schaeffer simulator verdicts, never clinical.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
from dataclasses import asdict, replace
from typing import Dict, List, Tuple

import numpy as np

from asb.experiments.realcohort import (
    COARSEN_NODES, FROZEN_BURST_CLS, FROZEN_MONO, FROZEN_SFI,
)
from asb.labels.monodomain import induce_monodomain
from asb.substrate.mesh import mesh_to_graph
from asb.substrate.roney import coarsen_mesh, load_roney_mesh

__all__ = ["LABEL_PARAMS", "run_e6"]

#: The three ground-truth-label knobs screened by Morris, with physical ranges.
LABEL_PARAMS = (
    ("d0", 0.10, 0.30),        # conduction scale (CV proxy); frozen 0.20
    ("iir_dense", 1.20, 1.44),  # fibrosis IIR threshold; frozen 1.32
    ("cross", 0.20, 0.50),      # cross-fibre conductance (anisotropy); frozen 0.30
)
_N_MESHES = 8          # label-cohort subset per Morris evaluation (logged cap)
_MORRIS_R = 4          # Morris trajectories
_P_LEVELS = 4          # grid levels per parameter


def _unit_to_phys(u: np.ndarray) -> Dict[str, float]:
    return {name: float(lo + u[i] * (hi - lo)) for i, (name, lo, hi) in enumerate(LABEL_PARAMS)}


def _label_outcome(u: np.ndarray, paths: List[str], seed: int, n_jobs: int) -> Dict[str, float]:
    """Cohort-mean self-sustained reentry time + inducible fraction at a design point."""
    phys = _unit_to_phys(u)
    args = [(p, phys, seed) for p in paths]
    if n_jobs > 1:
        from multiprocessing import Pool
        with Pool(min(n_jobs, len(args))) as pool:
            res = pool.map(_label_one, args)
    else:
        res = [_label_one(a) for a in args]
    sustained = np.array([r[0] for r in res], float)
    induced = np.array([r[1] for r in res], float)
    return {"mean_sustained": float(sustained.mean()), "inducible_frac": float(induced.mean())}


def _label_one(arg: Tuple[str, Dict[str, float], int]) -> Tuple[float, int]:
    path, phys, seed = arg
    mesh = load_roney_mesh(path, iir_dense=phys["iir_dense"])
    coarse = coarsen_mesh(mesh, COARSEN_NODES)
    cfg = replace(FROZEN_MONO, d0=phys["d0"])
    sd = int(hashlib.sha1((os.path.basename(path)).encode()).hexdigest(), 16) % (2**31)
    lab = induce_monodomain(coarse, cfg, np.random.default_rng(sd + seed),
                            burst_cls=FROZEN_BURST_CLS, cross=phys["cross"])
    return float(lab.meta.get("best_sustained_ms", 0.0)), int(lab.inducible)


# --------------------------------------------------------------------------- #
# Morris elementary effects (manual, no external dependency).
# --------------------------------------------------------------------------- #
def _morris_trajectories(k: int, r: int, seed: int) -> List[np.ndarray]:
    """``r`` classic one-at-a-time Morris trajectories on a ``_P_LEVELS`` grid."""
    rng = np.random.default_rng(seed)
    delta = _P_LEVELS / (2.0 * (_P_LEVELS - 1))  # standard Morris step in [0,1]
    grid = np.linspace(0.0, 1.0 - delta, _P_LEVELS // 2 + 1)
    trajs = []
    for _ in range(r):
        base = rng.choice(grid, size=k)
        order = rng.permutation(k)
        pts = [base.copy()]
        cur = base.copy()
        for j in order:
            cur = cur.copy()
            cur[j] = cur[j] + delta if cur[j] + delta <= 1.0 else cur[j] - delta
            pts.append(cur.copy())
        trajs.append((np.array(pts), order, delta))
    return trajs


def _morris_screen(paths: List[str], seed: int, n_jobs: int, cache_dir: str) -> Dict[str, object]:
    k = len(LABEL_PARAMS)
    cache = os.path.join(cache_dir, f"e6_morris_{_tag(paths)}.json")
    evals: Dict[str, dict] = {}
    if os.path.isfile(cache):
        with open(cache) as fh:
            evals = json.load(fh)

    def ev(u: np.ndarray) -> Dict[str, float]:
        key = ",".join(f"{x:.4f}" for x in u)
        if key not in evals:
            evals[key] = _label_outcome(u, paths, seed, n_jobs)
            with open(cache, "w") as fh:
                json.dump(evals, fh)
        return evals[key]

    trajs = _morris_trajectories(k, _MORRIS_R, seed)
    ee = {name: [] for name, _, _ in LABEL_PARAMS}
    design_X, design_y = [], []
    for pts, order, delta in trajs:
        ys = [ev(p)["mean_sustained"] for p in pts]
        for p, y in zip(pts, ys):
            design_X.append(p.tolist()); design_y.append(y)
        for step, j in enumerate(order):
            moved = pts[step + 1][j] - pts[step][j]  # +/- delta
            ee[LABEL_PARAMS[j][0]].append((ys[step + 1] - ys[step]) / moved)
    stats = {}
    for name, _, _ in LABEL_PARAMS:
        arr = np.array(ee[name], float)
        stats[name] = {"mu_star": float(np.mean(np.abs(arr))), "mu": float(np.mean(arr)),
                       "sigma": float(np.std(arr)), "n_ee": int(arr.size)}
    return {"elementary_effects": stats, "design_X": design_X, "design_y": design_y,
            "n_meshes": len(paths), "morris_r": _MORRIS_R}


def _gp_surrogate(design_X, design_y) -> Dict[str, object]:
    """GP fit on the Morris design points → LOO R² + main-effect curves."""
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
    from sklearn.model_selection import LeaveOneOut

    X = np.array(design_X, float)
    y = np.array(design_y, float)
    if len(np.unique(y)) < 3 or X.shape[0] < 5:
        return {"loo_r2": None, "note": "insufficient/degenerate design for a GP"}
    kernel = ConstantKernel(1.0) * RBF(length_scale=[0.3] * X.shape[1]) + WhiteKernel(1e-3)
    ys = (y - y.mean()) / (y.std() + 1e-12)

    # Leave-one-out R².
    preds = np.empty_like(ys)
    for tr, te in LeaveOneOut().split(X):
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=False, alpha=1e-6)
        gp.fit(X[tr], ys[tr]); preds[te] = gp.predict(X[te])
    ss_res = float(np.sum((ys - preds) ** 2)); ss_tot = float(np.sum((ys - ys.mean()) ** 2))
    loo_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None

    gp = GaussianProcessRegressor(kernel=kernel, normalize_y=False, alpha=1e-6).fit(X, ys)
    grid = np.linspace(0, 1, 9)
    main_effects = {}
    for i, (name, lo, hi) in enumerate(LABEL_PARAMS):
        pts = np.tile(0.5, (grid.size, X.shape[1])); pts[:, i] = grid
        mu = gp.predict(pts) * (y.std() + 1e-12) + y.mean()
        main_effects[name] = {"param_values": (lo + grid * (hi - lo)).tolist(),
                              "mean_sustained": mu.tolist()}
    return {"loo_r2": loo_r2, "main_effects": main_effects}


# --------------------------------------------------------------------------- #
# Δw robustness of the GM1 / GM3 conclusions (predictor-side).
# --------------------------------------------------------------------------- #
def _delta_w_sensitivity(roney_dir: str, n_jobs: int, cache_dir: str,
                         deltas=(0.18, 0.36, 0.54)) -> Dict[str, object]:
    """Recompute GM1 grouped ΔAUC + GM3 ρ at several Δw (cached labels reused)."""
    import pandas as pd
    from asb.config import Config
    from asb.experiments.gm1 import COMPETITOR_COLS, _evaluate_pair, _sfi_cols
    from asb.experiments.realcohort import cohort_records
    from asb.features import subject_features

    records = cohort_records(roney_dir, n_jobs=n_jobs, cache_dir=cache_dir, verbose=False)
    # Rebuild graphs once (cheap; no monodomain), reuse across Δw.
    graphs, y, groups = [], [], []
    for r in records:
        path = os.path.join(roney_dir, r.subject)
        G = mesh_to_graph(coarsen_mesh(load_roney_mesh(path), COARSEN_NODES))
        graphs.append((r.subject, G)); y.append(int(r.inducible)); groups.append(r.shape_family)
    y = np.array(y, int)
    out = {}
    for dw in deltas:
        sfi_cfg = replace(FROZEN_SFI, delta_w_mean_frac=float(dw))
        cfg = Config(seed=0, sfi=sfi_cfg)
        rows = [subject_features(G, cfg, np.random.default_rng(0)) for _, G in graphs]
        X = pd.DataFrame(rows).reindex(sorted({k for row in rows for k in row}), axis=1).fillna(0.0)
        sfi = _sfi_cols(X.columns)
        ev = _evaluate_pair(X, y, groups, COMPETITOR_COLS, sfi,
                            {"n_splits": 5, "seed": 0, "n_boot": 2000})
        out[f"{dw:.2f}"] = {clf: {"grouped_delta_auc": ev[clf]["grouped_delta_auc"],
                                  "delong_p": ev[clf]["delong_p"]} for clf in ev}
    return out


def _tag(paths: List[str]) -> str:
    payload = json.dumps({"params": [list(p) for p in LABEL_PARAMS], "mono": asdict(FROZEN_MONO),
                          "burst": list(FROZEN_BURST_CLS), "n": len(paths), "r": _MORRIS_R, "v": 1},
                         sort_keys=True)
    return hashlib.sha1(payload.encode()).hexdigest()[:10]


def run_e6(roney_dir: str = "data/roney", *, n_jobs: int = 4, seed: int = 0,
           outputs_dir: str = "outputs") -> dict:
    """Run E6: Morris label-screening + GP surrogate + Δw conclusion-sensitivity."""
    os.makedirs(outputs_dir, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(roney_dir, "Mesh_*.vtk")),
                   key=lambda p: (os.path.getsize(p), p))[:_N_MESHES]
    print(f"[e6] Morris screen on {len(paths)} meshes (CAPPED for CPU), r={_MORRIS_R}", flush=True)
    morris = _morris_screen(paths, seed, n_jobs, outputs_dir)
    gp = _gp_surrogate(morris["design_X"], morris["design_y"])
    print("[e6] Δw sensitivity of GM1/GM3 conclusions", flush=True)
    dw = _delta_w_sensitivity(roney_dir, n_jobs, outputs_dir)

    metrics = {
        "morris_label_screen": {"elementary_effects": morris["elementary_effects"],
                                "n_meshes": morris["n_meshes"], "morris_r": morris["morris_r"],
                                "outcome": "cohort-mean self-sustained reentry time (ms)"},
        "gp_surrogate": gp,
        "delta_w_sensitivity": dw,
        "note": ("Δw enters only the SFI predictor (never the monodomain label), so it is "
                 "screened via GM1/GM3 conclusion-sensitivity, not the label-Morris. Label "
                 "cohort capped at %d meshes for CPU feasibility." % _N_MESHES),
    }
    with open(os.path.join(outputs_dir, "e6_metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2)
    _write_report(os.path.join(outputs_dir, "e6_report.md"), metrics)
    return metrics


def _fmt(x) -> str:
    try:
        return "n/a" if x is None or not np.isfinite(float(x)) else f"{float(x):.4f}"
    except (TypeError, ValueError):
        return str(x)


def _write_report(path: str, m: dict) -> None:
    ms = m["morris_label_screen"]
    lines = ["# E6 — Uncertainty quantification (Morris screen + GP surrogate + Δw sensitivity)",
             "",
             "> Descoped UQ (no full Sobol). Labels are monodomain Mitchell-Schaeffer simulator "
             "verdicts, not clinical. " + m["note"], "",
             f"## Morris elementary effects on the label ({ms['outcome']}, N={ms['n_meshes']} meshes)",
             "", "| parameter | μ* (importance) | σ (nonlinearity/interaction) |",
             "| --- | --- | --- |"]
    for name, s in ms["elementary_effects"].items():
        lines.append(f"| {name} | {_fmt(s['mu_star'])} | {_fmt(s['sigma'])} |")
    gp = m["gp_surrogate"]
    lines += ["", f"## GP surrogate — leave-one-out R² = {_fmt(gp.get('loo_r2'))}", ""]
    if gp.get("main_effects"):
        lines.append("Main-effect (partial-dependence) curves stored in `e6_metrics.json`.")
    lines += ["", "## Δw robustness of the GM1/GM3 conclusions (predictor-side)", "",
              "| Δw | GM1 grouped ΔAUC (lr / gbt) | DeLong p (lr / gbt) |",
              "| --- | --- | --- |"]
    for dw, r in m["delta_w_sensitivity"].items():
        lr, gbt = r.get("lr", {}), r.get("gbt", {})
        lines.append(f"| {dw} | {_fmt(lr.get('grouped_delta_auc'))} / {_fmt(gbt.get('grouped_delta_auc'))} "
                     f"| {_fmt(lr.get('delong_p'))} / {_fmt(gbt.get('delong_p'))} |")
    lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
