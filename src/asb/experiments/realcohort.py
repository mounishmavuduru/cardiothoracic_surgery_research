r"""Real Roney cohort: load -> coarsen -> monodomain-MS label -> features (cached).

This is the config-driven bridge from the public Roney LA meshes (Zenodo 5801337)
to a per-subject record carrying the frozen competitor + spectral + SFI features and
the monodomain Mitchell--Schaeffer inducibility label. Labelling is the compute
bottleneck (a nonlinear PDE per subject), so records are cached to disk keyed by a
hash of the frozen configuration; re-runs are instant.

All frozen numbers come from ``docs/PRE_REGISTRATION.md`` section 7 (git-tagged
before any SFI-vs-label analysis). ``mock_ep`` is not used here; labels are
``source='monodomain_ms'`` simulator verdicts, never clinical POAF.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from asb.config import Config, SFIConfig
from asb.features import subject_features
from asb.labels.monodomain import MonodomainConfig, induce_monodomain
from asb.substrate.mesh import mesh_to_graph
from asb.substrate.roney import coarsen_mesh, load_roney_mesh

__all__ = [
    "FROZEN_SFI",
    "FROZEN_MONO",
    "FROZEN_BURST_CLS",
    "COARSEN_NODES",
    "cohort_records",
    "records_to_frame",
]

# --------------------------------------------------------------------------- #
# FROZEN configuration (pre-registration section 7). Do not tune to an endpoint.
# --------------------------------------------------------------------------- #
#: SFI feature field: literature Delta-w = 0.36 (f_CV 0.20 -> f_D = 1-0.8^2).
FROZEN_SFI = SFIConfig(
    n_monte_carlo=256,
    delta_w_mean_frac=0.36,
    delta_w_cov=0.5,
    fibrosis_coupling=1.0,
    validity_safety=0.25,
)
#: Monodomain Mitchell-Schaeffer ground-truth protocol (calibration gate PASSED:
#: inducible fraction 0.33, Spearman(fibrosis, reentry) 0.75).
FROZEN_MONO = MonodomainConfig(
    tau_in=0.3, tau_out=6.0, tau_open=120.0, tau_close=110.0, v_gate=0.13,
    fibrosis_erp_shortening=0.5,
    d0=0.20, w_ref=0.3, dt=0.05,
    protocol="burst", n_pacing_sites=2, n_burst=6,
    stim_amp=0.15, stim_duration=2.0, stim_radius_mm=4.0,
    observe_after=1000.0, reentry_min_ms=650.0,
)
FROZEN_BURST_CLS: Tuple[float, ...] = (150.0,)
#: Vertex-clustering coarsening target.
COARSEN_NODES = 2000


@dataclass
class SubjectRecord:
    """One real subject: features + monodomain-MS label + provenance."""

    subject: str
    shape_family: str
    n_nodes: int
    features: Dict[str, float]
    inducible: bool
    sustained_ms: float
    reentry_origin: Optional[int]
    fib_burden: float


def _config_tag() -> str:
    """Short hash of the frozen config so caches invalidate if any knob changes."""
    payload = json.dumps(
        {
            "sfi": asdict(FROZEN_SFI),
            "mono": asdict(FROZEN_MONO),
            "burst_cls": list(FROZEN_BURST_CLS),
            "coarsen": COARSEN_NODES,
            "v": 1,
        },
        sort_keys=True,
    )
    return hashlib.sha1(payload.encode()).hexdigest()[:10]


def _subject_seed(name: str) -> int:
    """Deterministic per-subject seed from its filename (not wall-clock)."""
    return int(hashlib.sha1(name.encode()).hexdigest(), 16) % (2**31)


def _process_one(path: str) -> Dict[str, object]:
    """Worker: load+coarsen a mesh, compute features, run the induction battery.

    Top-level (picklable) so it can run under ``multiprocessing.Pool``. Returns a
    plain dict (JSON-serializable) for caching.
    """
    name = os.path.basename(path)
    mesh = load_roney_mesh(path)
    coarse = coarsen_mesh(mesh, COARSEN_NODES)
    graph = mesh_to_graph(coarse)

    cfg = Config(seed=0, sfi=FROZEN_SFI)
    feats = subject_features(graph, cfg, np.random.default_rng(_subject_seed(name)))

    label = induce_monodomain(
        coarse, FROZEN_MONO, np.random.default_rng(_subject_seed(name)),
        burst_cls=FROZEN_BURST_CLS,
    )
    return {
        "subject": name,
        "shape_family": name,  # each patient is its own group (leakage-free)
        "n_nodes": int(coarse.n_points),
        "features": {k: float(v) for k, v in feats.items()},
        "inducible": bool(label.inducible),
        "sustained_ms": float(label.meta.get("best_sustained_ms", 0.0)),
        "reentry_origin": (None if label.reentry_origin is None
                           else int(label.reentry_origin)),
        "fib_burden": float(feats.get("fibhet_fibrosis_burden", 0.0)),
    }


def cohort_records(
    roney_dir: str,
    *,
    limit: Optional[int] = None,
    n_jobs: int = 4,
    cache_dir: str = "outputs",
    verbose: bool = True,
) -> List[SubjectRecord]:
    """Build (or load cached) per-subject records for the real Roney cohort.

    Parameters
    ----------
    roney_dir : str
        Directory of ``Mesh_*.vtk`` files (Zenodo 5801337).
    limit : int, optional
        Cap the number of subjects (smallest files first for reproducibility).
        ``None`` uses all present. The cap is logged, never silent.
    n_jobs : int
        Parallel workers for labelling.
    cache_dir : str
        Where to read/write the record cache.
    verbose : bool
        Print progress + cap logging.

    Returns
    -------
    list[SubjectRecord]
    """
    paths = sorted(glob.glob(os.path.join(roney_dir, "Mesh_*.vtk")),
                   key=lambda p: (os.path.getsize(p), p))
    n_available = len(paths)
    if limit is not None:
        paths = paths[:limit]
    if verbose:
        capped = "" if limit is None or limit >= n_available else \
            f"  (CAPPED from {n_available} available -> {len(paths)})"
        print(f"[cohort] {len(paths)} subjects{capped}", flush=True)

    os.makedirs(cache_dir, exist_ok=True)
    cache_path = os.path.join(cache_dir, f"cohort_records_{_config_tag()}.json")
    cached: Dict[str, dict] = {}
    if os.path.isfile(cache_path):
        with open(cache_path) as fh:
            cached = {r["subject"]: r for r in json.load(fh)}

    todo = [p for p in paths if os.path.basename(p) not in cached]
    if verbose:
        print(f"[cohort] {len(cached)} cached, {len(todo)} to label "
              f"(cache={os.path.basename(cache_path)})", flush=True)

    if todo:
        results: List[dict]
        if n_jobs > 1:
            from multiprocessing import Pool
            with Pool(min(n_jobs, len(todo))) as pool:
                results = pool.map(_process_one, todo)
        else:
            results = [_process_one(p) for p in todo]
        for r in results:
            cached[r["subject"]] = r
        # Persist the full cache.
        with open(cache_path, "w") as fh:
            json.dump(list(cached.values()), fh)

    records = []
    for p in paths:
        r = cached[os.path.basename(p)]
        records.append(SubjectRecord(
            subject=r["subject"], shape_family=r["shape_family"],
            n_nodes=r["n_nodes"], features=r["features"],
            inducible=r["inducible"], sustained_ms=r["sustained_ms"],
            reentry_origin=r["reentry_origin"], fib_burden=r["fib_burden"],
        ))
    if verbose:
        n_ind = sum(r.inducible for r in records)
        print(f"[cohort] inducible {n_ind}/{len(records)} "
              f"({n_ind/max(1,len(records)):.2f})", flush=True)
    return records


def records_to_frame(records: List[SubjectRecord]):
    """Stack records into (X DataFrame, y array, groups list) for evaluation."""
    import pandas as pd

    rows = [r.features for r in records]
    X = pd.DataFrame(rows).reindex(
        sorted({k for row in rows for k in row}), axis=1).fillna(0.0)
    y = np.array([int(r.inducible) for r in records], dtype=int)
    groups = [r.shape_family for r in records]
    return X, y, groups
