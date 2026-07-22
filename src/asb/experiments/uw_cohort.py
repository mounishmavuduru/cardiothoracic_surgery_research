r"""UW/Boyle real-cohort integration + expanded GM1.

Runs the 82 UW/Boyle pre-ablation LA meshes (``asb.substrate.uw_boyle``) through
the *frozen* Roney pipeline — same coarsening, same spectral/SFI features, same
monodomain Mitchell-Schaeffer inducibility labeller, same leakage-free per-patient
grouping — so they form a second real cohort. :func:`run_gm1_expanded` then reads
the SFI-vs-competitors predictive test on (a) Roney only, (b) UW only, and (c) the
combined ~182-patient cohort, a direct external-generalization check of the GM1
null on anatomy the pipeline was never tuned on.

Nothing here changes the frozen constants; it reuses ``realcohort``/``gm1`` verbatim.
"""
from __future__ import annotations

import glob
import json
import os
from typing import Dict, List, Optional

import numpy as np

from asb.config import Config
from asb.experiments.gm1 import (COMPETITOR_COLS, FIBHET_COLS, _evaluate_pair,
                                 _sfi_cols)
from asb.experiments.realcohort import (COARSEN_NODES, FROZEN_BURST_CLS,
                                        FROZEN_MONO, FROZEN_SFI, SubjectRecord,
                                        _config_tag, _subject_seed,
                                        cohort_records, records_to_frame)
from asb.features import subject_features
from asb.labels.monodomain import induce_monodomain
from asb.substrate.mesh import mesh_to_graph
from asb.substrate.roney import coarsen_mesh
from asb.substrate.uw_boyle import load_uw_mesh, uw_preablation_paths

__all__ = ["uw_cohort_records", "combined_records", "run_gm1_expanded"]


def _process_one_uw(path: str) -> Dict[str, object]:
    """Worker: load+coarsen one UW mesh, compute features, run the induction battery."""
    name = "uw_" + os.path.basename(path)
    mesh = load_uw_mesh(path)
    coarse = coarsen_mesh(mesh, COARSEN_NODES)
    graph = mesh_to_graph(coarse)

    cfg = Config(seed=0, sfi=FROZEN_SFI)
    feats = subject_features(graph, cfg, np.random.default_rng(_subject_seed(name)))
    label = induce_monodomain(
        coarse, FROZEN_MONO, np.random.default_rng(_subject_seed(name)),
        burst_cls=FROZEN_BURST_CLS,
    )
    return {
        "subject": name, "shape_family": name, "n_nodes": int(coarse.n_points),
        "features": {k: float(v) for k, v in feats.items()},
        "inducible": bool(label.inducible),
        "sustained_ms": float(label.meta.get("best_sustained_ms", 0.0)),
        "reentry_origin": (None if label.reentry_origin is None
                           else int(label.reentry_origin)),
        "fib_burden": float(feats.get("fibhet_fibrosis_burden", 0.0)),
    }


def uw_cohort_records(
    root: str = "data/uw_boyle", *, limit: Optional[int] = None, n_jobs: int = 2,
    cache_dir: str = "outputs", verbose: bool = True,
) -> List[SubjectRecord]:
    """Build (or load cached) per-patient records for the UW pre-ablation cohort."""
    paths = uw_preablation_paths(root)
    n_available = len(paths)
    if limit is not None:
        paths = paths[:limit]
    if verbose:
        capped = "" if limit is None or limit >= n_available else \
            f"  (CAPPED from {n_available} -> {len(paths)})"
        print(f"[uw-cohort] {len(paths)} pre-ablation subjects{capped}", flush=True)

    os.makedirs(cache_dir, exist_ok=True)
    cache_path = os.path.join(cache_dir, f"uw_cohort_records_{_config_tag()}.json")
    cached: Dict[str, dict] = {}
    if os.path.isfile(cache_path):
        with open(cache_path) as fh:
            cached = {r["subject"]: r for r in json.load(fh)}

    todo = [p for p in paths if "uw_" + os.path.basename(p) not in cached]
    if verbose:
        print(f"[uw-cohort] {len(cached)} cached, {len(todo)} to label", flush=True)
    if todo:
        if n_jobs > 1:
            from multiprocessing import Pool
            with Pool(min(n_jobs, len(todo))) as pool:
                results = pool.map(_process_one_uw, todo)
        else:
            results = [_process_one_uw(p) for p in todo]
        for r in results:
            cached[r["subject"]] = r
        with open(cache_path, "w") as fh:
            json.dump(list(cached.values()), fh)

    records = []
    for p in paths:
        r = cached["uw_" + os.path.basename(p)]
        records.append(SubjectRecord(
            subject=r["subject"], shape_family=r["shape_family"],
            n_nodes=r["n_nodes"], features=r["features"], inducible=r["inducible"],
            sustained_ms=r["sustained_ms"], reentry_origin=r["reentry_origin"],
            fib_burden=r["fib_burden"]))
    if verbose:
        n_ind = sum(r.inducible for r in records)
        print(f"[uw-cohort] inducible {n_ind}/{len(records)} "
              f"({n_ind / max(1, len(records)):.2f})", flush=True)
    return records


def combined_records(
    roney_dir: str = "data/roney", uw_root: str = "data/uw_boyle",
    *, n_jobs: int = 2, cache_dir: str = "outputs", verbose: bool = True,
) -> Dict[str, List[SubjectRecord]]:
    """Return ``{'roney': [...], 'uw': [...], 'combined': [...]}`` record sets."""
    roney = cohort_records(roney_dir, n_jobs=min(n_jobs, 4), cache_dir=cache_dir,
                           verbose=verbose)
    uw = uw_cohort_records(uw_root, n_jobs=n_jobs, cache_dir=cache_dir, verbose=verbose)
    return {"roney": roney, "uw": uw, "combined": roney + uw}


def _gm1_on(records: List[SubjectRecord], cfg: dict) -> dict:
    """SFI predictive test (competitors_vs_+SFI, fibrosis_vs_+SFI) on a record set."""
    X, y, groups = records_to_frame(records)
    sfi_cols = _sfi_cols(X.columns)
    n_ind = int(y.sum())
    out = {"n_subjects": int(len(y)), "n_inducible": n_ind,
           "inducible_fraction": float(n_ind / len(y)) if len(y) else 0.0,
           "n_groups": int(len(set(groups)))}
    if min(n_ind, len(y) - n_ind) >= 2:
        out["competitors_vs_+SFI"] = _evaluate_pair(X, y, groups, COMPETITOR_COLS,
                                                    sfi_cols, cfg)
        out["fibrosis_vs_+SFI"] = _evaluate_pair(X, y, groups, FIBHET_COLS,
                                                 sfi_cols, cfg)
    return out


def run_gm1_expanded(
    roney_dir: str = "data/roney", uw_root: str = "data/uw_boyle",
    *, n_jobs: int = 2, n_boot: int = 10000, n_splits: int = 5, seed: int = 0,
    outputs_dir: str = "outputs", metrics_path: str = "results/gm1_expanded_metrics.json",
) -> dict:
    """Expanded GM1: SFI predictive test on Roney-only, UW-only, and combined cohorts."""
    sets = combined_records(roney_dir, uw_root, n_jobs=n_jobs, cache_dir=outputs_dir)
    cfg = {"n_splits": n_splits, "seed": seed, "n_boot": n_boot}
    metrics = {"label_source": "monodomain_ms",
               "cohorts": {name: _gm1_on(recs, cfg) for name, recs in sets.items()},
               "config": {"n_boot": n_boot, "n_splits": n_splits, "seed": seed},
               "disclaimer": ("Monodomain Mitchell-Schaeffer verdicts on real LA "
                              "anatomy (Roney + UW/Boyle). NOT clinical POAF; NOT "
                              "openCARP. UW UAC is a PCA surrogate.")}
    os.makedirs(os.path.dirname(metrics_path), exist_ok=True)
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    return metrics
