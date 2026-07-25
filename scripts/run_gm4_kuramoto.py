"""GM4 third medium: Kuramoto oscillator networks — predict + localize battery.

Runs the same GM4 transfer tests (does SFI predict the instability? does |grad phi2|
localize its origin above a spatial null?) on Kuramoto phase-oscillator networks — a
different dynamical class from the two excitable media (cardiac, FHN). Reuses the exact
spectral/SFI feature stack and the gm4_scale aggregators.
"""
import json

import numpy as np

from asb.config import Config
from asb.experiments.gm2 import localizer_fields
from asb.experiments.gm3 import gm3_extra_features
from asb.experiments.gm4_scale import _localize_at, _predict_at
from asb.experiments.realcohort import FROZEN_SFI
from asb.features import subject_features
from asb.transfer.kuramoto import KuramotoConfig, induce_kuramoto
from asb.transfer.network import NetworkConfig, make_excitable_network

CFG = KuramotoConfig()
NET = NetworkConfig(n_nodes=400, radius=0.10, n_region_grid=3)
KEYS = ("grad_phi2", "perron", "combined", "wdegree", "fibrosis")
_BASE = 90000


def _proc(k: int):
    rng = np.random.default_rng(_BASE + k)
    burden = float(rng.uniform(0.05, 0.55))
    G = make_excitable_network(_BASE + k, NET, lesion_burden=burden)
    label = induce_kuramoto(G, CFG)
    feats = subject_features(G, Config(seed=0, sfi=FROZEN_SFI), np.random.default_rng(k))
    feats.update(gm3_extra_features(G, k_dim=2, include_exact=False))
    rec = {"seed": int(_BASE + k), "burden": burden, "inducible": bool(label.inducible),
           "features": {kk: float(v) for kk, v in feats.items()}}
    if label.inducible and label.reentry_origin is not None:
        o = int(label.reentry_origin)
        fields = localizer_fields(G)
        n = G.n_nodes
        r = np.random.default_rng(700000 + k)
        rand = r.integers(0, n, size=min(200, n))
        loc = {}
        for key in KEYS:
            fld = np.asarray(fields[key], float)
            loc[key] = {"rank": float(np.mean(fld >= fld[o])),
                        "perm": [float(np.mean(fld >= fld[int(v)])) for v in rand]}
        rec["loc"] = loc
    return rec


if __name__ == "__main__":
    import time
    from multiprocessing import Pool
    t0 = time.time()
    N = 500
    print(f"[kuramoto] building + labelling {N} oscillator networks", flush=True)
    with Pool(2) as pool:                    # gentle: leave cores for the 100k run
        recs = pool.map(_proc, range(N))
    n_ind = sum(r["inducible"] for r in recs)
    print(f"cohort N={N} unstable={n_ind} ({n_ind / N:.0%})  in {time.time() - t0:.0f}s", flush=True)

    pred = _predict_at(recs, n_boot=5000, n_splits=5, seed=0)
    loc = _localize_at(recs, seed=0)
    out = {"medium": "kuramoto_network", "dynamical_class": "phase_oscillator",
           "N": N, "n_unstable": n_ind, "n_nodes": NET.n_nodes,
           "predict": pred, "localize": loc}
    json.dump(out, open("results/gm4_kuramoto_metrics.json", "w", encoding="utf-8"), indent=2)

    print("[predict — SFI vs competitors]", flush=True)
    for v, r in pred.items():
        for clf in ("lr", "gbt"):
            d = r[clf]
            print(f"  {v} {clf}: dAUC={d['grouped_delta_auc']:+.4f} p={d['delong_p']:.3f} "
                  f"{'MET' if d['endpoint_met'] else 'null'}", flush=True)
    print("[localize — origin vs spatial null]", flush=True)
    for key, r in loc["fields"].items():
        print(f"  {key}: rank={r['mean_origin_rank']:.3f} p={r['perm_p']:.4f} {r['verdict']}",
              flush=True)
    print("KURAMOTO_DONE", flush=True)
