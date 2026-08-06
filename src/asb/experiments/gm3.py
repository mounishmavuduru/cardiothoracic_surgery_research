r"""GM3 (E4) — BOUND: the validity radius, and does a *correctly computed* SFI recover the signal?

GM1 returned a null: the first-order **single-vector** SFI does not beat the
competitors on real Roney anatomy. The diagnosed cause is that the perioperative
stress is thousands of times larger than the tiny :math:`\lambda_2-\lambda_3` gap,
so the single-Fiedler-vector derivative is applied far outside its perturbative
validity radius. GM3 tests the pre-registered fallback two ways.

**Part A — predictive recovery (three fragility computations, same 7 features each):**

- ``sfi``   : first-order single-vector  ``Σ E[Δw] (φ2_i − φ2_j)²``  (GM1's; invalid when λ2≈λ3);
- ``sfisub``: first-order **subspace/projector** (Davis–Kahan), ``Σ E[Δw] Σ_{c∈cluster}(φ_c,i − φ_c,j)²`` — well-defined under near-degeneracy;
- ``sfiex`` : **exact** per-region Δλ2 under the frozen mean field (recompute λ2 exactly) — valid at any ‖ΔL‖.

Each is a strict add-on to the competitor set, same labels, same nested
GroupKFold + DeLong + 10⁴ bootstrap machinery as GM1. If ``sfiex`` (or ``sfisub``)
recovers ΔAUC≥0.05, p<0.05 where ``sfi`` did not, the null was a *computation*
failure, not a *concept* failure — the honest-core headline.

**Part B — the validity radius (endpoint):** on synthetic graphs with a tunable
gap, sweep ρ = ‖ΔL‖/(λ3−λ2) and measure the relative error of the first-order and
subspace predictions of Δλ2 against the exact recompute; locate the ρ* where
first-order error crosses 10 %, show it matches the Weyl (eigenvalue) /
Davis–Kahan (subspace) scaling, and place the real cohort's ρ on that axis.
"""
from __future__ import annotations

import hashlib
import json
import os
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp

from asb.experiments.gm1 import (
    COMPETITOR_COLS, FIBHET_COLS, _evaluate_pair,
)
from asb.experiments.realcohort import (
    COARSEN_NODES, FROZEN_SFI, _subject_seed, cohort_records, records_to_frame,
)
from asb.sfi import edge_fragility, perturbation_field, subspace_sfi
from asb.spectral import fiedler, smallest_eigpairs
from asb.substrate.mesh import mesh_to_graph
from asb.substrate.roney import coarsen_mesh, load_roney_mesh

__all__ = ["run_gm3", "validity_radius_sweep", "gm3_extra_features"]

#: Below this magnitude an exact eigenvalue drop is at solver noise, so a *relative*
#: error against it carries no information. Rows at or under it report ``None``.
_EXACT_NOISE_FLOOR = 1e-12


def _rel_err(pred: float, exact: float) -> float | None:
    """Relative error of ``pred`` against ``exact``, or ``None`` at solver noise.

    Corrected 2026-08-06. The previous form divided by ``max(abs(exact), 1e-12)``.
    That floor did not guard a division by zero so much as manufacture a small
    error: in the most degenerate row of :func:`degeneracy_sweep` the exact λ2 drop
    is ~7e-15 against a single-vector prediction of ~4e-31 — a true relative error
    of essentially 100 % — which the floor reported as 0.7 %, making the
    single-vector estimator look most accurate exactly where it had failed
    completely. Returning ``None`` says "not measurable here", which is the truth.
    """
    if abs(exact) < _EXACT_NOISE_FLOOR:
        return None
    return abs(pred - exact) / abs(exact)


_SFI_TOP_K = 3


# --------------------------------------------------------------------------- #
# Region assignment (mirrors asb.sfi._edge_regions) + region summary features.
# --------------------------------------------------------------------------- #
def _edge_regions(edges: np.ndarray, region_ids: np.ndarray) -> np.ndarray:
    ri = region_ids[edges[:, 0]]
    rj = region_ids[edges[:, 1]]
    return np.where(ri == rj, ri.astype(np.int64), np.int64(-1))


def _region_sums(per_edge: np.ndarray, er: np.ndarray) -> np.ndarray:
    out = []
    for r in np.unique(er):
        out.append(float(per_edge[er == r].sum()))
    return np.asarray(out, dtype=float)


def _region_summary(values: np.ndarray, total: float, prefix: str) -> Dict[str, float]:
    """The 7 comparable summary features shared by every SFI variant."""
    v = np.asarray(values, dtype=float)
    out = {
        f"{prefix}total": float(total),
        f"{prefix}region_max": float(v.max()) if v.size else 0.0,
        f"{prefix}region_mean": float(v.mean()) if v.size else 0.0,
        f"{prefix}region_std": float(v.std()) if v.size else 0.0,
    }
    ranked = np.sort(v)[::-1] if v.size else np.array([])
    for r in range(_SFI_TOP_K):
        out[f"{prefix}top{r + 1}"] = float(ranked[r]) if r < ranked.size else 0.0
    return out


def _laplacian(G) -> sp.csr_matrix:
    W = G.adjacency()
    d = np.asarray(W.sum(axis=1)).ravel()
    return sp.csr_matrix(sp.diags(d) - W)


def _exact_region_delta_lambda2(
    L: sp.csr_matrix, edges: np.ndarray, expected_dw: np.ndarray,
    region_ids: np.ndarray, lam2_base: float,
) -> Tuple[np.ndarray, float]:
    """Exact per-region drop in λ2 when the region's mean uncoupling is applied.

    For each region, form L' = L - ΔL_region (ΔL_region built from the region's
    edges scaled by their mean expected_dw) and recompute λ2 exactly; the drop is
    ``lam2_base - λ2(L')``. Returns the per-region drops (region order = np.unique)
    and the whole-atrium drop.
    """
    n = L.shape[0]
    er = _edge_regions(edges, region_ids)
    regions = np.unique(er)

    def _dl(mask: np.ndarray) -> float:
        e = edges[mask]
        dw = expected_dw[mask]
        if e.shape[0] == 0:
            return 0.0
        i, j = e[:, 0], e[:, 1]
        rows = np.concatenate([i, j]); cols = np.concatenate([j, i])
        data = np.concatenate([dw, dw])
        Wd = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
        dL = sp.diags(np.asarray(Wd.sum(1)).ravel()) - Wd
        vals, _ = smallest_eigpairs((L - dL).tocsr(), k=2)
        return float(lam2_base - vals[1])

    drops = np.array([_dl(er == r) for r in regions], dtype=float)
    total = _dl(np.ones(edges.shape[0], dtype=bool))
    return drops, total


def gm3_extra_features(graph, cfg=FROZEN_SFI, *, k_dim: int = 2,
                       include_exact: bool = True) -> Dict[str, float]:
    """Subspace-SFI (``sfisub_``) and exact-Δλ2 SFI (``sfiex_``) summary features.

    Both are computed on the *same* frozen expected uncoupling field as the
    single-vector ``sfi_`` block, so the three are directly comparable add-ons.
    ``include_exact=False`` skips the (eigsh-heavy) exact per-region Δλ2 block so
    the cheap subspace features can be extracted for very large cohorts.
    """
    G = graph
    L = _laplacian(G)
    edges = np.asarray(G.edges, np.int64)
    region = np.asarray(G.region, np.int64)
    field = perturbation_field(G, cfg, np.random.default_rng(0))
    expected_dw = field["expected_dw"]

    er = _edge_regions(edges, region)

    # --- subspace SFI (well-defined under near-degeneracy) ---
    sub_edge = subspace_sfi(L, edges, expected_dw, k_dim=k_dim)
    sub_regions = _region_sums(sub_edge, er)
    feats = _region_summary(sub_regions, float(sub_edge.sum()), "sfisub_")

    # --- exact per-region Δλ2 under the mean field (valid at any ‖ΔL‖) ---
    if include_exact:
        lam2_base, _ = fiedler(L)
        ex_regions, ex_total = _exact_region_delta_lambda2(
            L, edges, expected_dw, region, lam2_base)
        feats.update(_region_summary(ex_regions, ex_total, "sfiex_"))
    return feats


# --------------------------------------------------------------------------- #
# Cache of GM3 extra features (recomputing them needs the graph, not just labels).
# --------------------------------------------------------------------------- #
def _gm3_tag(k_dim: int) -> str:
    payload = json.dumps({"coarsen": COARSEN_NODES, "kdim": k_dim,
                          "dw": FROZEN_SFI.delta_w_mean_frac, "v": 1}, sort_keys=True)
    return hashlib.sha1(payload.encode()).hexdigest()[:10]


def _extra_for_path(path: str, k_dim: int) -> Dict[str, object]:
    mesh = coarsen_mesh(load_roney_mesh(path), COARSEN_NODES)
    G = mesh_to_graph(mesh)
    feats = gm3_extra_features(G, k_dim=k_dim)
    return {"subject": os.path.basename(path), "features": feats}


def _cohort_extra_features(
    roney_dir: str, subjects: Sequence[str], *, k_dim: int, n_jobs: int,
    cache_dir: str,
) -> Dict[str, Dict[str, float]]:
    cache_path = os.path.join(cache_dir, f"gm3_features_{_gm3_tag(k_dim)}.json")
    cached: Dict[str, dict] = {}
    if os.path.isfile(cache_path):
        with open(cache_path, encoding="utf-8") as fh:
            cached = {r["subject"]: r["features"] for r in json.load(fh)}
    todo = [os.path.join(roney_dir, s) for s in subjects if s not in cached]
    if todo:
        print(f"[gm3] computing subspace+exact SFI for {len(todo)} subjects", flush=True)
        if n_jobs > 1:
            from multiprocessing import Pool
            with Pool(min(n_jobs, len(todo))) as pool:
                results = pool.map(_extra_for_path_kdim(k_dim), todo)
        else:
            results = [_extra_for_path(p, k_dim) for p in todo]
        for r in results:
            cached[r["subject"]] = r["features"]
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump([{"subject": s, "features": f} for s, f in cached.items()], fh)
    return cached


class _extra_for_path_kdim:
    """Picklable partial of _extra_for_path with a bound k_dim (for Pool.map)."""
    def __init__(self, k_dim: int):
        self.k_dim = k_dim

    def __call__(self, path: str) -> Dict[str, object]:
        return _extra_for_path(path, self.k_dim)


# --------------------------------------------------------------------------- #
# Part B — validity radius sweep on synthetic graphs.
# --------------------------------------------------------------------------- #
def _path_graph(n: int):
    """Path graph P_n (unit weights): simple λ2, analytic gap λ_k=2-2cos(kπ/n)."""
    edges = np.array([[k, k + 1] for k in range(n - 1)], np.int64)
    w = np.ones(edges.shape[0], dtype=float)
    return edges, w


def _lap(n, edges, w):
    i, j = edges[:, 0], edges[:, 1]
    rows = np.concatenate([i, j]); cols = np.concatenate([j, i])
    data = np.concatenate([w, w])
    W = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    return (sp.diags(np.asarray(W.sum(1)).ravel()) - W).tocsr()


def _single_edge_dL(n, edge, g):
    i, j = int(edge[0]), int(edge[1])
    rows = [i, j, i, j]; cols = [j, i, i, j]; data = [g, g, -g, -g]
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n))


def validity_radius_sweep(
    *, n: int = 16, block: int = 3,
    t_grid: Sequence[float] = (0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.55, 0.75, 0.95),
) -> Dict[str, object]:
    """Trace the first-order SFI validity boundary on a simple-λ2 path graph.

    On the path P_n (simple λ2, moderate gap) reduce a *localized block* of the
    first ``block`` edges by fraction ``t`` (a localized perturbation couples to
    the higher modes, so second-order error is visible). For each ``t`` we report
    the achieved ρ = ‖ΔL‖/(λ3−λ2), the exact Δλ2, the first-order prediction
    ``Σ t·w·(φ2_i−φ2_j)²`` and its relative error, the 2-D subspace prediction
    error, and whether the Weyl bound ``|Δλ2| ≤ ‖ΔL‖`` holds. ρ* is where the
    first-order relative error first exceeds 10 %.
    """
    edges, w = _path_graph(n)
    L = _lap(n, edges, w)
    vals, vecs = smallest_eigpairs(L, k=4)
    lam2, lam3 = float(vals[1]), float(vals[2])
    gap = lam3 - lam2
    phi2 = vecs[:, 1]
    frag = edge_fragility(phi2, edges)
    block_mask = np.zeros(edges.shape[0], dtype=bool)
    block_mask[:block] = True

    rows = []
    for t in t_grid:
        dw = np.where(block_mask, t * w, 0.0)
        i, j = edges[:, 0], edges[:, 1]
        rr = np.concatenate([i, j]); cc = np.concatenate([j, i]); dd = np.concatenate([dw, dw])
        Wd = sp.csr_matrix((dd, (rr, cc)), shape=(n, n))
        dL = (sp.diags(np.asarray(Wd.sum(1)).ravel()) - Wd).tocsr()
        dLnorm = float(np.abs(np.linalg.eigvalsh(dL.toarray())).max())
        rho = dLnorm / gap
        pred_first = float(np.sum(dw * frag))
        pred_sub = float(np.sum(subspace_sfi(L, edges, dw, k_dim=2)))
        vex, _ = smallest_eigpairs((L - dL).tocsr(), k=3)
        exact = float(lam2 - vex[1])
        # ``subspace_sfi(k_dim=2)`` predicts the drop in the SUM λ2+λ3, not in λ2
        # alone, so it must be scored against that sum. Corrected 2026-08-06: this
        # row previously compared it to the exact λ2 drop, an apples-to-oranges
        # target that reported a relative error near 12 and read as though the
        # subspace estimator were a thousand times worse than the first-order one.
        # ``degeneracy_sweep`` below always used the correct target.
        exact_sum = float((lam2 + lam3) - (vex[1] + vex[2]))
        rel_first = _rel_err(pred_first, exact)
        rel_sub = _rel_err(pred_sub, exact_sum)
        rows.append({
            "t": float(t), "rho": float(rho), "dL_norm": float(dLnorm),
            "exact_dlam2": exact, "exact_sum_dlam23": exact_sum,
            "pred_first": pred_first, "pred_subspace": pred_sub,
            "rel_err_first": rel_first, "rel_err_subspace": rel_sub,
            "weyl_ok": bool(abs(exact) <= dLnorm + 1e-9),
        })

    rho_star = next((r["rho"] for r in rows
                     if r["rel_err_first"] is not None and r["rel_err_first"] > 0.10), None)
    return {"lam2": lam2, "lam3": lam3, "gap": gap, "rho_star_10pct": rho_star,
            "sweep": rows}


def degeneracy_sweep(
    *, n_clust: int = 8, bridges: Sequence[float] = (1.0, 0.3, 0.1, 0.03, 0.01, 0.003),
    perturb_frac: float = 0.2, seed: int = 0,
) -> Dict[str, object]:
    """Show the subspace SFI beats the single vector as λ2 → degeneracy.

    Three equal cliques on a path of two bridges; shrinking ``bridges`` drives the
    two lowest non-trivial modes together (λ2 → λ3, near-degenerate). At a fixed
    diffuse perturbation we compare the single-vector prediction of the λ2 drop
    against the 2-D subspace prediction of the drop in the *sum* λ2+λ3 (the
    degeneracy-robust invariant), each scored against its own exact target.

    What this sweep shows, measured 2026-08-06 after fixing the two defects below
    (the per-bridge mask redraw, and the relative-error floor in :func:`_rel_err`):
    the subspace error is lower than the single-vector error at every gap tested —
    0.0681 vs 0.0846 at gap 0.167, falling through to 0.0760 vs 0.0906 at gap
    0.00075 — and both series are smooth and monotone. The ``bridge = 1.0`` control
    row reproduces the previously stored value exactly, which validates the setup.

    What it does **not** show, and what was withdrawn from ``SFI_THEORY.md`` §5:
    that the single-vector error blows up as the gap closes. It converges to about
    0.091. The defensible statement is that the subspace form is better conditioned
    and consistently more accurate, not that the single-vector form diverges here.

    ``results/gm3_metrics.json`` still holds the pre-fix values for the two affected
    columns; it is not regenerated because the Roney arm has grown from 62 to 100
    subjects, so a full re-run is not a like-for-like replacement. No manuscript
    number depends on either column.
    """
    rng = np.random.default_rng(seed)
    n = 3 * n_clust
    base_edges = []
    for c in range(3):
        b = c * n_clust
        for a in range(n_clust):
            for d in range(a + 1, n_clust):
                base_edges.append((b + a, b + d))
    bridge_pairs = [(0, n_clust), (n_clust, 2 * n_clust)]
    edges = np.array(base_edges + bridge_pairs, np.int64)
    # Drawn ONCE, outside the sweep. Until 2026-08-06 this was redrawn inside the
    # loop, so the perturbation moved together with the gap and the two error
    # series could not be attributed to the gap at all -- which is the only thing
    # the sweep varies on purpose. The docstring already said "at a fixed diffuse
    # perturbation"; now the code does that.
    mask = rng.random(edges.shape[0]) < perturb_frac
    rows = []
    for br in bridges:
        w = np.ones(edges.shape[0]); w[-2:] = br
        L = _lap(n, edges, w)
        vals, vecs = smallest_eigpairs(L, k=5)
        lam2, lam3, lam4 = float(vals[1]), float(vals[2]), float(vals[3])
        gap23 = lam3 - lam2
        dw = np.where(mask, perturb_frac * w, 0.0)
        i, j = edges[:, 0], edges[:, 1]
        rr = np.concatenate([i, j]); cc = np.concatenate([j, i]); dd = np.concatenate([dw, dw])
        Wd = sp.csr_matrix((dd, (rr, cc)), shape=(n, n))
        dL = (sp.diags(np.asarray(Wd.sum(1)).ravel()) - Wd).tocsr()
        vex, _ = smallest_eigpairs((L - dL).tocsr(), k=3)
        # (i) single-vector predicts the λ2 drop; compare to the exact λ2 drop —
        # the ill-posed target under degeneracy. Whether its error grows with the
        # gap is what the sweep measures, not something to assert here.
        exact_lam2_drop = float(lam2 - vex[1])
        frag2 = edge_fragility(vecs[:, 1], edges)
        pred_single = float(np.sum(dw * frag2))
        err_single = _rel_err(pred_single, exact_lam2_drop)
        # (ii) subspace predicts the degeneracy-invariant λ2+λ3 drop; compare to exact.
        exact_sum_drop = float((lam2 + lam3) - (vex[1] + vex[2]))
        pred_sub = float(np.sum(subspace_sfi(L, edges, dw, k_dim=2)))
        err_sub = _rel_err(pred_sub, exact_sum_drop)
        rows.append({
            "bridge": float(br), "lam2": lam2, "lam3": lam3, "gap23": float(gap23),
            "exact_lam2_drop": exact_lam2_drop, "exact_sum_drop": exact_sum_drop,
            "pred_single": pred_single, "pred_subspace": pred_sub,
            "err_single_vec": err_single, "err_subspace": err_sub,
        })
    return {"sweep": rows}


# --------------------------------------------------------------------------- #
# Orchestrator.
# --------------------------------------------------------------------------- #
def run_gm3(
    roney_dir: str = "data/roney",
    *, limit: int | None = None, n_jobs: int = 4, n_boot: int = 10000,
    n_splits: int = 5, seed: int = 0, k_dim: int = 2, outputs_dir: str = "outputs",
) -> dict:
    """Run GM3 Part A (subspace/exact-SFI head-to-head) + Part B (validity sweep)."""
    records = cohort_records(roney_dir, limit=limit, n_jobs=n_jobs,
                             cache_dir=outputs_dir)
    subjects = [r.subject for r in records]
    extra = _cohort_extra_features(roney_dir, subjects, k_dim=k_dim,
                                   n_jobs=n_jobs, cache_dir=outputs_dir)

    # Merge extra features into each record's feature dict.
    for r in records:
        r.features.update(extra.get(r.subject, {}))

    X, y, groups = records_to_frame(records)
    y = np.asarray(y, int)
    n_ind = int(y.sum())

    sfi_single = sorted(c for c in X.columns if c.startswith("sfi_")
                        and not c.startswith(("sfisub_", "sfiex_")))
    # Restrict single-vector to the same 7 region-summary features for a fair match.
    seven = ("total", "region_max", "region_mean", "region_std", "top1", "top2", "top3")
    sfi_single7 = [f"sfi_{s}" for s in seven if f"sfi_{s}" in X.columns]
    sfi_sub = sorted(c for c in X.columns if c.startswith("sfisub_"))
    sfi_ex = sorted(c for c in X.columns if c.startswith("sfiex_"))

    cfg = {"n_splits": n_splits, "seed": seed, "n_boot": n_boot}
    variants = {"single_vector": sfi_single7, "subspace": sfi_sub, "exact": sfi_ex}

    comparisons: Dict[str, object] = {}
    if min(n_ind, len(y) - n_ind) >= 2:
        for vname, cols in variants.items():
            comparisons[f"competitors_vs_+{vname}"] = _evaluate_pair(
                X, y, groups, COMPETITOR_COLS, cols, cfg)
            comparisons[f"fibrosis_vs_+{vname}"] = _evaluate_pair(
                X, y, groups, FIBHET_COLS, cols, cfg)

    sweep = validity_radius_sweep()
    degen = degeneracy_sweep()

    metrics = {
        "cohort": {"n_subjects": int(len(y)), "n_inducible": n_ind,
                   "inducible_fraction": float(n_ind / len(y)) if len(y) else 0.0},
        "variant_columns": {k: list(v) for k, v in variants.items()},
        "comparisons": comparisons,
        "validity_radius_sweep": sweep,
        "degeneracy_sweep": degen,
        "config": {"n_boot": n_boot, "k_dim": k_dim, "label_source": "monodomain_ms"},
        "label_disclaimer": (
            "Labels are monodomain Mitchell-Schaeffer simulator verdicts on real "
            "Roney LA anatomy. NOT clinical POAF; NOT openCARP (deferred)."),
    }
    os.makedirs(outputs_dir, exist_ok=True)
    with open(os.path.join(outputs_dir, "gm3_metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2)
    _write_report(os.path.join(outputs_dir, "gm3_report.md"), metrics)
    return metrics


def _fmt(x) -> str:
    # ``None`` now reaches here from :func:`_rel_err`, meaning "not measurable at
    # solver noise" rather than "missing"; both render as n/a in the report.
    if x is None:
        return "n/a"
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "n/a" if not np.isfinite(xf) else f"{xf:.4f}"


def _write_report(path: str, m: dict) -> None:
    c = m["cohort"]
    lines = [
        "# GM3 — validity radius + subspace/exact SFI recovery (real Roney, monodomain-MS labels)",
        "",
        "> Labels are a monodomain Mitchell-Schaeffer **simulator verdict**, not clinical POAF, not openCARP.",
        "",
        f"- Subjects **{c['n_subjects']}**, inducible **{c['n_inducible']}** "
        f"({_fmt(c['inducible_fraction'])}), patient-held-out GroupKFold.",
        "",
        "## Part A — does a correctly-computed SFI recover the signal?",
        "",
        "| add-on to competitors | clf | grouped base | grouped +SFI | grouped ΔAUC | DeLong p | boot ΔAUC [95% CI] | endpoint |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for name, comp in m["comparisons"].items():
        if not name.startswith("competitors_vs_"):
            continue
        variant = name.replace("competitors_vs_+", "")
        for clf, r in comp.items():
            b = r["bootstrap"]
            lines.append(
                f"| {variant} | {clf} | {_fmt(r['grouped_auc_base'])} | "
                f"{_fmt(r['grouped_auc_sfi'])} | {_fmt(r['grouped_delta_auc'])} | "
                f"{_fmt(r['delong_p'])} | {_fmt(b['delta_mean'])} "
                f"[{_fmt(b['ci_low'])}, {_fmt(b['ci_high'])}] | "
                f"{'MET' if r['endpoint_met'] else 'not met'} |")
    lines += ["", "## Part B — validity radius (synthetic, exact vs first-order vs subspace)", ""]
    sw = m["validity_radius_sweep"]
    lines += [f"- λ2={_fmt(sw['lam2'])}, λ3={_fmt(sw['lam3'])}, gap={_fmt(sw['gap'])}; "
              f"**ρ\\* (first-order 10% error) = {sw['rho_star_10pct']}**", "",
              "| ρ=‖ΔL‖/gap | exact Δλ2 | first-order err | subspace err | Weyl ok |",
              "| --- | --- | --- | --- | --- |"]
    for r in sw["sweep"]:
        lines.append(f"| {_fmt(r['rho'])} | {_fmt(r['exact_dlam2'])} | "
                     f"{_fmt(r['rel_err_first'])} | {_fmt(r['rel_err_subspace'])} | "
                     f"{r['weyl_ok']} |")
    lines += ["", "> Real-cohort ρ (from GM1 validity scan) has median ≈ 2422 — far right of "
              "this table, i.e. deep in the regime where only the exact recompute is trustworthy.", ""]
    dg = m.get("degeneracy_sweep", {})
    if dg:
        lines += ["## Part B(ii) — near-degeneracy: subspace SFI vs single Fiedler vector", "",
                  "As the bridge shrinks the λ2–λ3 gap closes; error is vs the exact drop in "
                  "the degeneracy-invariant sum λ2+λ3.", "",
                  "| bridge | gap(λ3−λ2) | single-vector err | subspace err |",
                  "| --- | --- | --- | --- |"]
        for r in dg["sweep"]:
            lines.append(f"| {_fmt(r['bridge'])} | {_fmt(r['gap23'])} | "
                         f"{_fmt(r['err_single_vec'])} | {_fmt(r['err_subspace'])} |")
        lines.append("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
