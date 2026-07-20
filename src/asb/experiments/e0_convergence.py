r"""E0 — calculation verification: mesh & Δt convergence of λ₂ and the label.

ASME V&V-40 calculation verification for the two quantities the science rests on:
the Fiedler value :math:`\lambda_2` (spectral engine) and the monodomain-MS
inducibility verdict (ground-truth labeller). Both must be stable under mesh
refinement and time-step reduction, else nothing downstream is interpretable.

This is the CPU-feasible half of E0. The openCARP Niederer N-version benchmark
(PMID 21969679) is the *other* half and is deferred to Claude Science (openCARP is
not installable here); it is wired in :func:`asb.labels.opencarp.reproduce_niederer_benchmark`.
"""
from __future__ import annotations

import json
import os
from typing import Dict, List

import numpy as np
import scipy.sparse as sp

from asb.labels.monodomain import induce_monodomain
from asb.spectral import fiedler, spectral_gap
from asb.substrate.mesh import mesh_to_graph
from asb.substrate.roney import coarsen_mesh, load_roney_mesh
from asb.experiments.realcohort import FROZEN_BURST_CLS, FROZEN_MONO

__all__ = ["lambda2_convergence", "label_dt_convergence", "run_e0"]


def _lambda2_gap(mesh) -> tuple[float, float]:
    G = mesh_to_graph(mesh)
    W = G.adjacency()
    L = sp.csr_matrix(sp.diags(np.asarray(W.sum(axis=1)).ravel()) - W)
    lam2, _ = fiedler(L)
    return float(lam2), float(spectral_gap(L))


def lambda2_convergence(
    path: str, resolutions: List[int]
) -> List[Dict[str, float]]:
    """λ₂ and spectral gap vs coarsening resolution for one mesh.

    Reports the relative change between successive (finer) resolutions; a small
    relative change is mesh-convergence of the spectral quantity.
    """
    mesh = load_roney_mesh(path)
    rows: List[Dict[str, float]] = []
    prev = None
    for r in sorted(resolutions):
        coarse = coarsen_mesh(mesh, r)
        lam2, gap = _lambda2_gap(coarse)
        rel = float("nan") if prev is None else abs(lam2 - prev) / max(abs(prev), 1e-12)
        rows.append({"target": r, "n_nodes": int(coarse.n_points),
                     "lambda2": lam2, "gap": gap, "rel_change_lambda2": rel})
        prev = lam2
    return rows


def label_dt_convergence(
    path: str, dts: List[float], *, coarsen_nodes: int = 2000
) -> List[Dict[str, object]]:
    """Inducibility label + sustained-reentry ms vs time step for one mesh.

    A stable label (and a bounded change in sustained_ms) under Δt reduction is
    the temporal-convergence evidence.
    """
    from dataclasses import replace

    mesh = coarsen_mesh(load_roney_mesh(path), coarsen_nodes)
    rows: List[Dict[str, object]] = []
    for dt in dts:
        cfg = replace(FROZEN_MONO, dt=float(dt))
        lab = induce_monodomain(mesh, cfg, np.random.default_rng(0),
                                burst_cls=FROZEN_BURST_CLS)
        rows.append({"dt": float(dt), "inducible": bool(lab.inducible),
                     "sustained_ms": float(lab.meta["best_sustained_ms"])})
    return rows


def run_e0(
    paths: List[str],
    *, resolutions: List[int] | None = None, dts: List[float] | None = None,
    outputs_dir: str = "outputs",
) -> dict:
    """Run the E0 convergence study on a few meshes and write ``e0_convergence.json``."""
    resolutions = resolutions or [1500, 2000, 3000]
    dts = dts or [0.1, 0.05, 0.025]

    lam2 = {os.path.basename(p): lambda2_convergence(p, resolutions) for p in paths}
    # Δt convergence is expensive (a label per dt) -> just the first mesh.
    dt_conv = label_dt_convergence(paths[0], dts) if paths else []

    # Convergence verdicts.
    lam2_max_rel = max(
        (row["rel_change_lambda2"] for rows in lam2.values() for row in rows
         if np.isfinite(row["rel_change_lambda2"])), default=float("nan"))
    labels_stable = len({r["inducible"] for r in dt_conv}) <= 1 if dt_conv else None

    result = {
        "lambda2_convergence": lam2,
        "label_dt_convergence": dt_conv,
        "lambda2_max_rel_change": float(lam2_max_rel),
        "label_stable_across_dt": labels_stable,
        "note": ("openCARP Niederer benchmark is the deferred other half of E0 "
                 "(solver not installable here); see labels/opencarp.py."),
    }
    os.makedirs(outputs_dir, exist_ok=True)
    with open(os.path.join(outputs_dir, "e0_convergence.json"), "w") as fh:
        json.dump(result, fh, indent=2)
    return result
