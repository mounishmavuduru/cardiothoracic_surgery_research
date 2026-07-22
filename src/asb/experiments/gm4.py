r"""GM4 (E5) — TRANSFER: does the spectral-fragility calculus generalize off the heart?

Ports the *entire* fragility calculus to a **second excitable medium** — a 2-D
FitzHugh--Nagumo neural network with an epileptic-focus (hyperexcitable-lesion)
instability ground truth (:mod:`asb.transfer`) — and asks whether the **pattern of
atrial results replicates**. Crucially the spectral engine, the SFI, the competitor
baselines and the localizer fields are reused **verbatim** (an ``AtrialGraph`` is just
a weighted graph), so GM4 is a literal test of medium-independence.

Three transfer tests mirror GM1/GM2/GM3:

- **Part A (predict, GM1-analog):** does SFI (single-vector / subspace / exact) beat the
  competitor set at classifying FHN instability? Nested GroupKFold + DeLong + bootstrap.
- **Part B (localize, GM2-analog):** does ``|∇φ2|`` (and the ``|∇φ2|∩Perron`` hotspot)
  localize the instability origin above a rotational spatial null? Tie-robust permutation
  keep/delete, per field.
- **Part C (validity, GM3-analog):** the ρ = ‖ΔL‖/(λ3−λ2) distribution across the
  network cohort, against the medium-agnostic validity boundary ρ*≈3.

Honest scope: that :math:`\lambda_2` governs excitable-network stability (Pecora--Carroll
MSF 1998; Fiedler value rises at seizure onset, Bomela 2020) and the
:math:`\partial\lambda_2/\partial w=(\varphi_i-\varphi_j)^2` identity (Ghosh--Boyd 2006)
are prior art. GM4 tests only the **spatial-localization** use and its **cross-medium
transfer**. FHN labels are simulator verdicts, never clinical.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict
from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sp

from asb.config import Config
from asb.experiments.gm1 import COMPETITOR_COLS, FIBHET_COLS, _evaluate_pair
from asb.experiments.gm2 import geodesic_from, localizer_fields, _torus_shift_null_nodes
from asb.experiments.gm3 import gm3_extra_features, validity_radius_sweep
from asb.experiments.realcohort import FROZEN_SFI
from asb.features import subject_features
from asb.sfi import perturbation_field
from asb.spectral import fiedler, spectral_gap
from asb.transfer.fhn import FHNConfig, induce_fhn
from asb.transfer.network import NetworkConfig, make_excitable_network

__all__ = ["FROZEN_NETWORK", "FROZEN_FHN", "run_gm4"]

# Frozen second-medium configuration.
FROZEN_NETWORK = NetworkConfig(n_nodes=900, radius=0.065, n_region_grid=3)
FROZEN_FHN = FHNConfig()  # calibrated: ~0.47 unstable, Spearman(lesion,unstable)=0.68
_N_NETWORKS = 80
_BURDEN_LO, _BURDEN_HI = 0.05, 0.45
_LOCALIZER_KEYS = ("grad_phi2", "perron", "combined", "wdegree", "fibrosis", "fibrosis_grad")


def _laplacian(G) -> sp.csr_matrix:
    W = G.adjacency()
    return sp.csr_matrix(sp.diags(np.asarray(W.sum(1)).ravel()) - W)


def _net_seed(k: int) -> int:
    return k  # networks are deterministic in their integer seed


def _process_network(k: int) -> Dict[str, object]:
    """Build network k, FHN-label it, extract all features (single+subspace+exact SFI)."""
    rng = np.random.default_rng(20259 + k)
    burden = float(rng.uniform(_BURDEN_LO, _BURDEN_HI))
    G = make_excitable_network(_net_seed(k), FROZEN_NETWORK, lesion_burden=burden)
    label = induce_fhn(G, FROZEN_FHN)
    cfg = Config(seed=0, sfi=FROZEN_SFI)
    feats = subject_features(G, cfg, np.random.default_rng(k))
    feats.update(gm3_extra_features(G, k_dim=2))
    return {
        "subject": f"net_{k}", "shape_family": f"net_{k}", "seed": int(k),
        "burden": burden, "n_nodes": int(G.n_nodes),
        "features": {kk: float(v) for kk, v in feats.items()},
        "inducible": bool(label.inducible),
        "reentry_origin": (None if label.reentry_origin is None else int(label.reentry_origin)),
    }


def _tag() -> str:
    payload = json.dumps({"net": asdict(FROZEN_NETWORK), "fhn": asdict(FROZEN_FHN),
                          "sfi": asdict(FROZEN_SFI), "n": _N_NETWORKS, "v": 2},
                         sort_keys=True)
    return hashlib.sha1(payload.encode()).hexdigest()[:10]


def _build_cohort(n_jobs: int, cache_dir: str, verbose: bool) -> List[dict]:
    os.makedirs(cache_dir, exist_ok=True)
    cache = os.path.join(cache_dir, f"gm4_cohort_{_tag()}.json")
    if os.path.isfile(cache):
        with open(cache) as fh:
            return json.load(fh)
    if verbose:
        print(f"[gm4] building + FHN-labelling {_N_NETWORKS} networks", flush=True)
    ks = list(range(_N_NETWORKS))
    if n_jobs > 1:
        from multiprocessing import Pool
        with Pool(n_jobs) as pool:
            recs = pool.map(_process_network, ks)
    else:
        recs = [_process_network(k) for k in ks]
    with open(cache, "w") as fh:
        json.dump(recs, fh)
    return recs


# --------------------------------------------------------------------------- #
# Part B — localization keep/delete with tie-robust permutation null.
# --------------------------------------------------------------------------- #
def _localize(recs: List[dict], n_null: int, seed: int) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    subj = [r for r in recs if r["inducible"] and r["reentry_origin"] is not None]
    per_field: Dict[str, Dict[str, list]] = {k: {"obs_norm": [], "origin_rank": [],
                                                 "perm_rank": [], "pass": []}
                                             for k in _LOCALIZER_KEYS}
    for r in subj:
        G = make_excitable_network(r["seed"], FROZEN_NETWORK, lesion_burden=r["burden"])
        origin = int(r["reentry_origin"])
        fields = localizer_fields(G)
        geo = geodesic_from(G, origin)
        diam = float(geo[np.isfinite(geo)].max()) if np.isfinite(geo).any() else 1.0
        n = G.n_nodes
        for name, field in fields.items():
            hot = int(np.argmax(field))
            obs = float(geo[hot])
            null_nodes = _torus_shift_null_nodes(G.uac, hot, n_null, rng)
            nd = geo[null_nodes]; nd = nd[np.isfinite(nd)]
            per_field[name]["obs_norm"].append(obs / diam if diam > 0 else np.nan)
            per_field[name]["origin_rank"].append(float(np.mean(field >= field[origin])))
            rand = rng.integers(0, n, size=min(300, n))
            per_field[name]["perm_rank"].append(
                [float(np.mean(field >= field[int(v)])) for v in rand])
            per_field[name]["pass"].append(
                bool(obs <= (np.percentile(nd, 5) if nd.size else np.inf)))

    out: Dict[str, object] = {"n_localized": len(subj), "origin_def": "sustained_core",
                              "fields": {}}
    for name in _LOCALIZER_KEYS:
        ranks = np.array(per_field[name]["origin_rank"], float)
        if ranks.size == 0:
            continue
        obs_mean = float(np.mean(ranks))
        perm_stacks = [np.asarray(x, float) for x in per_field[name]["perm_rank"]]
        pmin = min((a.size for a in perm_stacks), default=0)
        if pmin:
            draws = min(pmin, 5000)
            picks = np.stack([a[rng.integers(0, a.size, size=draws)] for a in perm_stacks])
            perm_means = picks.mean(axis=0)
            perm_p = float(np.mean(perm_means <= obs_mean))
            perm_null = float(np.mean(perm_means))
        else:
            perm_p, perm_null = float("nan"), float("nan")
        out["fields"][name] = {
            "mean_origin_rank": obs_mean, "perm_null_mean": perm_null,
            "perm_p": perm_p, "endpoint_pass_frac": float(np.mean(per_field[name]["pass"])),
            "verdict": "KEEP" if (np.isfinite(perm_p) and perm_p < 0.05) else "DELETE",
        }
    return out


# --------------------------------------------------------------------------- #
# Part B robustness — origin-definition sensitivity.
# --------------------------------------------------------------------------- #
def _perm_verdict(origin_ranks: np.ndarray, perm_stacks: list, rng) -> tuple:
    """Cohort mean origin-rank + tie-robust permutation p (random-origin null)."""
    obs = float(np.mean(origin_ranks))
    pmin = min((a.size for a in perm_stacks), default=0)
    if not pmin:
        return obs, float("nan")
    draws = min(pmin, 5000)
    picks = np.stack([a[rng.integers(0, a.size, size=draws)] for a in perm_stacks])
    return obs, float(np.mean(picks.mean(axis=0) <= obs))


def _origin_sensitivity(recs: List[dict], seed: int) -> Dict[str, object]:
    """How the keep/delete verdict depends on the (ambiguous) origin definition.

    Re-simulates each unstable network, extracts three origin definitions
    (``sustained_core`` / ``first_activation`` / ``earliest_last``), and reports the
    permutation-null verdict of the three key localizers under each — the honest
    transparency check the GM4 audit demanded.
    """
    from asb.transfer.fhn import origin_variants, simulate_fhn

    rng = np.random.default_rng(seed)
    keys = ("grad_phi2", "perron", "wdegree")
    defs = ("sustained_core", "first_activation", "earliest_last")
    acc = {d: {k: {"rank": [], "perm": []} for k in keys} for d in defs}
    subj = [r for r in recs if r["inducible"]]
    for r in subj:
        G = make_excitable_network(r["seed"], FROZEN_NETWORK, lesion_burden=r["burden"])
        res = simulate_fhn(G, FROZEN_FHN)
        variants = origin_variants(res)
        fields = localizer_fields(G)
        n = G.n_nodes
        for d in defs:
            o = variants.get(d)
            if o is None:
                continue
            for k in keys:
                fld = fields[k]
                acc[d][k]["rank"].append(float(np.mean(fld >= fld[int(o)])))
                rand = rng.integers(0, n, size=min(300, n))
                acc[d][k]["perm"].append([float(np.mean(fld >= fld[int(v)])) for v in rand])
    table: Dict[str, object] = {}
    for d in defs:
        table[d] = {}
        for k in keys:
            ranks = np.array(acc[d][k]["rank"], float)
            if ranks.size == 0:
                continue
            stacks = [np.asarray(x, float) for x in acc[d][k]["perm"]]
            mean_rank, perm_p = _perm_verdict(ranks, stacks, rng)
            table[d][k] = {"mean_rank": mean_rank, "perm_p": perm_p,
                           "verdict": "KEEP" if (np.isfinite(perm_p) and perm_p < 0.05) else "DELETE"}
    return table


# --------------------------------------------------------------------------- #
# Part C — validity radius on the network cohort.
# --------------------------------------------------------------------------- #
def _validity_scan(recs: List[dict]) -> Dict[str, object]:
    ratios = []
    for r in recs[:24]:
        G = make_excitable_network(r["seed"], FROZEN_NETWORK, lesion_burden=r["burden"])
        L = _laplacian(G)
        lam2, _ = fiedler(L)
        gap = spectral_gap(L)
        field = perturbation_field(G, FROZEN_SFI, np.random.default_rng(0))
        dw = field["expected_dw"]; e = np.asarray(G.edges)
        n = G.n_nodes; i, j = e[:, 0], e[:, 1]
        rows = np.concatenate([i, j]); cols = np.concatenate([j, i]); data = np.concatenate([dw, dw])
        Wd = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
        dL = sp.diags(np.asarray(Wd.sum(1)).ravel()) - Wd
        from scipy.sparse.linalg import eigsh
        dLnorm = float(eigsh(dL.tocsr(), k=1, which="LM", return_eigenvectors=False)[0])
        ratios.append(dLnorm / max(gap, 1e-12))
    ratios = np.array(ratios, float)
    # First-order SFI is trustworthy when ρ = ‖ΔL‖/gap <= validity_safety.
    return {"n_scanned": len(ratios), "rho_median": float(np.median(ratios)),
            "rho_min": float(ratios.min()), "rho_max": float(ratios.max()),
            "frac_in_radius": float(np.mean(ratios <= FROZEN_SFI.validity_safety))}


def run_gm4(
    *, n_jobs: int = 4, n_boot: int = 10000, n_null: int = 2000, n_splits: int = 5,
    seed: int = 0, outputs_dir: str = "outputs",
) -> dict:
    """Run the full GM4 transfer battery and write ``gm4_metrics.json`` + report."""
    import pandas as pd

    recs = _build_cohort(n_jobs, outputs_dir, verbose=True)
    y = np.array([int(r["inducible"]) for r in recs], int)
    groups = [r["shape_family"] for r in recs]
    rows = [r["features"] for r in recs]
    X = pd.DataFrame(rows).reindex(sorted({k for row in rows for k in row}), axis=1).fillna(0.0)
    n_ind = int(y.sum())

    seven = ("total", "region_max", "region_mean", "region_std", "top1", "top2", "top3")
    variants = {
        "single_vector": [f"sfi_{s}" for s in seven if f"sfi_{s}" in X.columns],
        "subspace": sorted(c for c in X.columns if c.startswith("sfisub_")),
        "exact": sorted(c for c in X.columns if c.startswith("sfiex_")),
    }
    cfg = {"n_splits": n_splits, "seed": seed, "n_boot": n_boot}

    predict: Dict[str, object] = {}
    if min(n_ind, len(y) - n_ind) >= 2:
        for vname, cols in variants.items():
            predict[f"competitors_vs_+{vname}"] = _evaluate_pair(
                X, y, groups, COMPETITOR_COLS, cols, cfg)
            predict[f"lesion_vs_+{vname}"] = _evaluate_pair(
                X, y, groups, FIBHET_COLS, cols, cfg)

    localize = _localize(recs, n_null, seed)
    origin_sensitivity = _origin_sensitivity(recs, seed)
    validity = _validity_scan(recs)
    sweep = validity_radius_sweep()

    metrics = {
        "medium": "fhn_network (2-D FitzHugh-Nagumo excitable network)",
        "cohort": {"n_networks": int(len(y)), "n_unstable": n_ind,
                   "unstable_fraction": float(n_ind / len(y)) if len(y) else 0.0},
        "predict": predict,
        "localize": localize,
        "origin_sensitivity": origin_sensitivity,
        "validity_cohort": validity,
        "validity_boundary_synthetic": {"rho_star_10pct": sweep["rho_star_10pct"]},
        "config": {"n_boot": n_boot, "n_null": n_null, "label_source": "fhn_network"},
        "label_disclaimer": ("FitzHugh-Nagumo excitable-network instability, a nonlinear "
                             "simulator verdict on an abstract neural medium. Not clinical."),
    }
    os.makedirs(outputs_dir, exist_ok=True)
    with open(os.path.join(outputs_dir, "gm4_metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2)
    _write_report(os.path.join(outputs_dir, "gm4_report.md"), metrics)
    return metrics


def _fmt(x) -> str:
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "n/a" if not np.isfinite(xf) else f"{xf:.4f}"


def _write_report(path: str, m: dict) -> None:
    c = m["cohort"]
    lines = [
        "# GM4 — TRANSFER: spectral-fragility calculus on a second excitable medium",
        "",
        "> Second medium: a 2-D FitzHugh-Nagumo neural network with an epileptic-focus "
        "(hyperexcitable-lesion) instability. FHN labels are simulator verdicts, not clinical. "
        "The spectral/SFI/baseline code is reused **verbatim** from the cardiac pipeline.",
        "",
        f"- Networks **{c['n_networks']}**, unstable **{c['n_unstable']}** "
        f"({_fmt(c['unstable_fraction'])}), each network its own group.",
        "",
        "## Part A — predict FHN instability (GM1-analog, grouped)",
        "",
        "| add-on | clf | grouped base | grouped +SFI | ΔAUC | DeLong p | endpoint |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for name, comp in m["predict"].items():
        for clf, r in comp.items():
            lines.append(f"| {name.replace('_vs_+',' vs +')} | {clf} | "
                         f"{_fmt(r['grouped_auc_base'])} | {_fmt(r['grouped_auc_sfi'])} | "
                         f"{_fmt(r['grouped_delta_auc'])} | {_fmt(r['delong_p'])} | "
                         f"{'MET' if r['endpoint_met'] else 'not met'} |")
    loc = m["localize"]
    lines += ["", f"## Part B — localize the instability origin (GM2-analog, N={loc['n_localized']}, "
              f"origin = {loc.get('origin_def','sustained_core')})", "",
              "| localizer | mean origin-rank | perm-null | perm p | verdict |",
              "| --- | --- | --- | --- | --- |"]
    for name, r in loc["fields"].items():
        lines.append(f"| {name} | {_fmt(r['mean_origin_rank'])} | {_fmt(r['perm_null_mean'])} "
                     f"| {_fmt(r['perm_p'])} | **{r['verdict']}** |")
    os_ = m.get("origin_sensitivity", {})
    if os_:
        lines += ["", "### Part B robustness — origin-definition sensitivity (perm p; verdict)",
                  "", "| origin definition | grad_phi2 | perron | wdegree |",
                  "| --- | --- | --- | --- |"]
        for d, row in os_.items():
            cells = []
            for k in ("grad_phi2", "perron", "wdegree"):
                rr = row.get(k)
                cells.append(f"{_fmt(rr['perm_p'])} {rr['verdict']}" if rr else "n/a")
            lines.append(f"| {d} | {cells[0]} | {cells[1]} | {cells[2]} |")
    v = m["validity_cohort"]
    lines += ["", "## Part C — validity radius on the network cohort (GM3-analog)", "",
              f"- ρ = ‖ΔL‖/(λ3−λ2): median **{_fmt(v['rho_median'])}** "
              f"(min {_fmt(v['rho_min'])}, max {_fmt(v['rho_max'])}) over {v['n_scanned']} networks.",
              f"- Synthetic validity boundary ρ* (first-order 10% error) ≈ "
              f"{m['validity_boundary_synthetic']['rho_star_10pct']}.", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
