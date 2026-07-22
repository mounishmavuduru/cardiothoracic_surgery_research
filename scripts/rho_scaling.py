"""Novelty experiment 2: the rho scaling law.

'rho ~ 2000' is a number; this makes it a *law*. The validity radius is
rho = ||dL|| / (lambda3 - lambda2). The gap (lambda3 - lambda2) of a diffusively-
coupled medium collapses with system size by Weyl's law -- for a d-dimensional
domain the low Laplacian eigenvalues scale as lambda_k ~ (k/N)^(2/d), so the gap
scales as N^(-2/d) and rho GROWS as N^(2/d). Consequence: the linear spectral
biomarker fails *more*, and predictably, as the medium is resolved finer -- the
opposite of the usual 'more data helps'. This script measures the gap-scaling
exponent across (i) 2-D and 3-D random-geometric media at fixed mean degree and
(ii) real atrial-surface meshes coarsened to varying N, and tests each against the
Weyl prediction 2/d.

Honest scope: Weyl's law is classical; the contribution is applying it to
characterize biomarker-validity scaling and confirming it empirically, with the
same exponent, across independent media.
"""
import json

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree


def _laplacian_gap(edges, n):
    """(lambda2, lambda3, gap) of the unweighted graph Laplacian, largest component."""
    if edges.shape[0] == 0:
        return np.nan, np.nan, np.nan
    i, j = edges[:, 0], edges[:, 1]
    A = sp.csr_matrix((np.ones(2 * len(i)), (np.concatenate([i, j]),
                       np.concatenate([j, i]))), shape=(n, n))
    ncomp, lab = connected_components(A, directed=False)
    if ncomp > 1:
        keep = lab == int(np.argmax(np.bincount(lab, minlength=ncomp)))
        idx = np.flatnonzero(keep)
        remap = -np.ones(n, np.int64); remap[idx] = np.arange(idx.size)
        m = keep[edges[:, 0]] & keep[edges[:, 1]]
        edges = remap[edges[m]]; n = idx.size
        i, j = edges[:, 0], edges[:, 1]
        A = sp.csr_matrix((np.ones(2 * len(i)), (np.concatenate([i, j]),
                           np.concatenate([j, i]))), shape=(n, n))
    deg = np.asarray(A.sum(1)).ravel()
    L = sp.diags(deg) - A
    k = 4
    vals = eigsh(L.astype(float), k=k, which="SM",
                 v0=np.ones(n) / np.sqrt(n), return_eigenvectors=False, maxiter=5000)
    vals = np.sort(vals)
    return float(vals[1]), float(vals[2]), float(vals[2] - vals[1])


def _rand_geom(n, d, mean_deg, seed):
    """Random-geometric graph in [0,1]^d with radius set for ~constant mean degree."""
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0, 1, size=(n, d))
    # expected neighbours ~ n * V_d * r^d ; solve r for target mean_deg
    from math import gamma, pi
    Vd = pi ** (d / 2) / gamma(d / 2 + 1)
    r = (mean_deg / (n * Vd)) ** (1.0 / d)
    tree = cKDTree(pts)
    pairs = tree.query_pairs(r=r, output_type="ndarray")
    if pairs.shape[0] == 0:
        _, nn = tree.query(pts, k=min(7, n))
        pairs = np.array([[a, b] for a, row in enumerate(nn[:, 1:]) for b in row], np.int64)
    return np.unique(np.sort(pairs, axis=1), axis=0).astype(np.int64)


def _fit_exponent(Ns, gaps):
    """Slope of log(gap) vs log(N): gap ~ N^slope (expect slope = -2/d)."""
    Ns = np.asarray(Ns, float); gaps = np.asarray(gaps, float)
    ok = np.isfinite(gaps) & (gaps > 0)
    if ok.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(Ns[ok]), np.log(gaps[ok]), 1)[0])


def main():
    out = {"random_geometric": {}, "atrial_mesh": {}}
    NS = [300, 600, 1200, 2400, 4800]

    for d in (2, 3):
        rows = []
        for N in NS:
            gs = []
            for s in range(4):
                e = _rand_geom(N, d, mean_deg=12, seed=1000 * d + 10 * s + N)
                gs.append(_laplacian_gap(e, N))
            gap = float(np.nanmedian([g[2] for g in gs]))
            rows.append({"N": N, "gap": gap})
            print(f"[rand-geom d={d}] N={N} gap={gap:.5f}", flush=True)
        slope = _fit_exponent([r["N"] for r in rows], [r["gap"] for r in rows])
        out["random_geometric"][f"d{d}"] = {"rows": rows, "fitted_exponent": slope,
                                             "weyl_prediction": -2.0 / d}
        print(f"[rand-geom d={d}] fitted gap~N^{slope:.3f}  (Weyl 2/d = {-2.0/d:.3f})", flush=True)

    # Real atrial surface (a 2-manifold => expect exponent ~ -1) via coarsening.
    try:
        from asb.substrate.roney import coarsen_mesh
        from asb.substrate.uw_boyle import load_uw_mesh, uw_preablation_paths
        from asb.substrate.mesh import mesh_to_graph
        paths = uw_preablation_paths()
        if paths:
            mesh = load_uw_mesh(paths[0])
            rows = []
            for target in (400, 800, 1600, 3200, 6400):
                cm = coarsen_mesh(mesh, target)
                g = mesh_to_graph(cm)
                e = np.asarray(g.edges, np.int64)
                l2, l3, gap = _laplacian_gap(e, g.n_nodes)
                rows.append({"target": target, "n": int(g.n_nodes), "gap": gap})
                print(f"[atrial mesh] n={g.n_nodes} gap={gap:.5f}", flush=True)
            slope = _fit_exponent([r["n"] for r in rows], [r["gap"] for r in rows])
            out["atrial_mesh"] = {"rows": rows, "fitted_exponent": slope,
                                  "weyl_prediction_2manifold": -1.0}
            print(f"[atrial mesh] fitted gap~N^{slope:.3f}  (2-manifold Weyl = -1.0)", flush=True)
    except Exception as e:  # noqa: BLE001
        out["atrial_mesh"] = {"error": str(e)}
        print(f"[atrial mesh] skipped: {e}", flush=True)

    json.dump(out, open("results/rho_scaling.json", "w"), indent=2)
    print("RHO_SCALING_DONE", flush=True)


if __name__ == "__main__":
    main()
