"""Novelty experiment 1 (real cardiac model): source-sink reentry at a healthy isthmus.

Monodomain reaction-diffusion (crossfield S1-S2) on a 2-D sheet with a healthy narrow
isthmus between two scar blocks. Dense scar is treated as non-conducting (edges removed)
in the conduction graph -- physically correct, and it makes the isthmus a sharp Fiedler
bottleneck (|grad phi2| peaks there, ~10x elsewhere, validated). The fibrosis peak sits in
the scar, spatially dissociated from the isthmus. The simulator decides independently
whether reentry nucleates at the isthmus; we then measure whether |grad phi2| localizes the
origin (KEEP) while fibrosis-burden misses it (DELETE) -- the connectivity-beats-substrate
claim, in the one model class (wavefront-curvature source-sink) that can dissociate them.

Idempotent: skips if results/monodomain_isthmus_metrics.json already exists.
"""
import json
import os

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

from asb.experiments.gm2 import localizer_fields
from asb.labels.monodomain import MonodomainConfig, simulate_monodomain
from asb.substrate.mesh import mesh_edges
from asb.substrate.sheet import SheetConfig, make_isthmus_sheet
from asb.types import AtrialGraph

SCAR = 0.5  # fibrosis above this is non-conducting scar (edges removed)
OUT = "results/monodomain_isthmus_metrics.json"


def _scar_removed_graph(mesh):
    """AtrialGraph with dense-scar edges removed, restricted to the largest component."""
    edges = mesh_edges(mesh.faces)
    fib = mesh.fibrosis
    n = mesh.n_points
    keep = ~((fib[edges[:, 0]] > SCAR) | (fib[edges[:, 1]] > SCAR))
    e = edges[keep]
    A = sp.csr_matrix((np.ones(2 * len(e)), (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])),
                      shape=(n, n))
    ncomp, lab = connected_components(A, directed=False)
    big = lab == np.argmax(np.bincount(lab))
    idx = np.flatnonzero(big)
    remap = -np.ones(n, np.int64); remap[idx] = np.arange(idx.size)
    e2 = remap[e[big[e[:, 0]] & big[e[:, 1]]]]
    g = AtrialGraph(coords=mesh.points[idx], edges=e2, weights=np.ones(len(e2)),
                    fibrosis=fib[idx], uac=mesh.uac[idx], region=mesh.region[idx],
                    shape_family=mesh.shape_family, meta=dict(mesh.meta))
    return g, remap, big


def _rank(field, node):
    field = np.asarray(field, float)
    return float(np.mean(field >= field[node]))


def main():
    if os.path.exists(OUT):
        print(f"{OUT} exists -> skip", flush=True); return
    import time
    t0 = time.time()
    widths = [4.0, 6.0, 8.0]
    cis = [160.0, 190.0, 220.0]
    yoff = [-6.0, 0.0, 6.0]     # vary isthmus vertical position -> distinct substrates
    recs = []
    for w in widths:
        for dy in yoff:
            for ci in cis:
                cfg = SheetConfig(lx=60, ly=60, spacing=0.8, isthmus_w=w)
                mesh = make_isthmus_sheet(0, cfg)
                # shift the isthmus channel by dy (rebuild fibrosis around shifted mid-line)
                if dy != 0.0:
                    y = mesh.points[:, 1]
                    dm = np.abs(y - (30.0 + dy))
                    inc = (mesh.points[:, 0] >= 20.0) & (mesh.points[:, 0] <= 40.0)
                    edge = np.clip((dm - w / 2.0) / 2.0, 0.0, 1.0)
                    mesh.fibrosis = np.where(inc, 0.85 * edge, 0.0)
                    mesh.meta["isthmus_xy"] = [30.0, 30.0 + dy]
                    mesh.meta["fibrosis_peak"] = int(np.argmax(mesh.fibrosis))
                g, remap, big = _scar_removed_graph(mesh)
                fields = localizer_fields(g)
                mc = MonodomainConfig(protocol="crossfield", s2_coupling=ci)
                res = simulate_monodomain(mesh, mc, np.array([0]), record_activation=True)
                rec = {"isthmus_w": w, "dy": dy, "s2_ci": ci,
                       "inducible": bool(res["inducible"]),
                       "sustained_ms": float(res.get("sustained_ms", 0.0))}
                if res["inducible"] and res["reentry_origin"] is not None:
                    o_full = int(res["reentry_origin"])
                    oc = mesh.points[o_full, :2]
                    isx = np.array(mesh.meta["isthmus_xy"])
                    fp = mesh.points[mesh.meta["fibrosis_peak"], :2]
                    rec["dist_origin_isthmus"] = float(np.linalg.norm(oc - isx))
                    rec["dist_origin_fibrosis"] = float(np.linalg.norm(oc - fp))
                    o_g = int(remap[o_full]) if big[o_full] else -1
                    if o_g >= 0:
                        rec["grad_phi2_rank"] = _rank(fields["grad_phi2"], o_g)
                        rec["fibrosis_rank"] = _rank(fields["fibrosis"], o_g)
                recs.append(rec)
                print(f"[iw={w} dy={dy} ci={ci}] inducible={rec['inducible']} "
                      f"sustained={rec['sustained_ms']:.0f} "
                      f"gp_rank={rec.get('grad_phi2_rank','-')} fib_rank={rec.get('fibrosis_rank','-')}",
                      flush=True)

    ind = [r for r in recs if r["inducible"] and "grad_phi2_rank" in r]
    summary = {"n_configs": len(recs), "n_inducible_localized": len(ind)}
    if ind:
        summary["grad_phi2_mean_rank"] = float(np.mean([r["grad_phi2_rank"] for r in ind]))
        summary["fibrosis_mean_rank"] = float(np.mean([r["fibrosis_rank"] for r in ind]))
        summary["mean_dist_origin_isthmus"] = float(np.mean([r["dist_origin_isthmus"] for r in ind]))
        summary["mean_dist_origin_fibrosis"] = float(np.mean([r["dist_origin_fibrosis"] for r in ind]))
    out = {"records": recs, "summary": summary, "elapsed_s": time.time() - t0,
           "note": "crossfield monodomain on a healthy-isthmus sheet; scar non-conducting"}
    os.makedirs("results", exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=2)
    print(f"MONODOMAIN_ISTHMUS_DONE  {summary}", flush=True)


if __name__ == "__main__":
    main()
