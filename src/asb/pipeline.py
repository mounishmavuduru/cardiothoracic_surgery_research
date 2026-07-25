"""End-to-end AtrialSpectralBench pipeline (config-driven, synthetic data).

Wires the whole engine together on the in-repo synthetic cohort:

``make_cohort`` -> ``mesh_to_graph`` -> spectral + SFI features + ``mock_ep`` labels
-> ``build_design_matrix`` -> grouped/naive cross-validated AUC (fibrosis vs
fibrosis+SFI) -> paired DeLong test -> hotspot/reentry colocalization-vs-null ->
figures.

Outputs (written under ``cfg.outputs_dir``)
-------------------------------------------
- ``metrics.json``        : all scalar results (AUCs, DeLong, colocalization).
- ``results_report.md``   : a human-readable summary.
- ``fig_fiedler_atrium.png``, ``fig_sfi_vs_origin.png``, ``fig_eigen_spectrum.png``,
  ``fig_roc_panel.png``, ``fig_sensitivity_bars.png``.

The ``mock_ep`` labels are a development stand-in for openCARP — **never** clinical
post-operative AF. Everything is deterministic in ``cfg.seed``.

Runnable as ``python -m asb.pipeline --config configs/default.yaml`` (see also the
``asb`` console script in :mod:`asb.cli`).
"""
from __future__ import annotations

import argparse
import json
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp

from asb.config import Config
from asb.evaluation import colocalization_vs_null, delong_test, nested_group_kfold_auc
from asb.features import build_design_matrix, subject_features  # noqa: F401
from asb.figures import (
    fig_eigen_spectrum,
    fig_fiedler_atrium,
    fig_roc_panel,
    fig_sensitivity_bars,
    fig_sfi_vs_origin,
)
from asb.labels.mock_ep import induce
from asb.sfi import edge_fragility, hotspot_map, perturbation_field, sfi_region
from asb.spectral import fiedler, perron, smallest_eigpairs
from asb.substrate.mesh import mesh_to_graph
from asb.substrate.synthetic import make_cohort
from asb.types import AtrialGraph, AtrialMesh, InducibilityLabel

__all__ = ["run", "main"]

_FEATURE_SET_FIBROSIS = "fibrosis"
_FEATURE_SET_FIBROSIS_SFI = "fibrosis+SFI"


def _combinatorial_laplacian(G: AtrialGraph) -> sp.csr_matrix:
    """Combinatorial Laplacian ``L = D - W`` from a graph's adjacency."""
    W = G.adjacency()
    d = np.asarray(W.sum(axis=1)).ravel()
    return sp.csr_matrix(sp.diags(d) - W)


def _grouped_oof_probs(X, y, groups, cols, cfg: Config) -> np.ndarray:
    """Group-out-of-fold positive-class probabilities for a column subset.

    Mirrors the leak-controlled protocol of
    :func:`asb.evaluation.nested_group_kfold_auc` (GroupKFold + standardized logistic
    regression) but returns the per-subject OOF probabilities themselves so the two
    feature sets can be compared with a paired DeLong test on identical labels.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    y = np.asarray(y).astype(int).ravel()
    groups = np.asarray(groups).ravel()
    Xf = np.asarray(X.loc[:, list(cols)].to_numpy(), dtype=float)

    n_groups = len(np.unique(groups))
    n_splits = int(min(cfg.eval.n_splits, n_groups))
    n_splits = max(n_splits, 2)

    oof = np.full(len(y), np.nan, dtype=float)
    gkf = GroupKFold(n_splits=min(n_splits, n_groups))
    for train_idx, test_idx in gkf.split(Xf, y, groups):
        if len(np.unique(y[train_idx])) < 2:
            oof[test_idx] = float(y[train_idx].mean()) if len(train_idx) else 0.5
            continue
        clf = Pipeline(
            [
                ("scale", StandardScaler()),
                ("clf", LogisticRegression(
                    max_iter=1000, solver="liblinear", random_state=cfg.seed)),
            ]
        )
        clf.fit(Xf[train_idx], y[train_idx])
        oof[test_idx] = clf.predict_proba(Xf[test_idx])[:, 1]

    if np.isnan(oof).any():
        oof = np.where(np.isnan(oof), float(y.mean()), oof)
    return oof


def _select_feature_sets(columns) -> Dict[str, List[str]]:
    """Build the {"fibrosis", "fibrosis+SFI"} column selections from ``X.columns``."""
    fibrosis_cols = sorted(c for c in columns if c.startswith("fibhet_"))
    sfi_cols = sorted(c for c in columns if c.startswith("sfi_"))
    return {
        _FEATURE_SET_FIBROSIS: fibrosis_cols,
        _FEATURE_SET_FIBROSIS_SFI: fibrosis_cols + sfi_cols,
    }


def _pick_colocalization_subject(
    graphs, labels
) -> Optional[int]:
    """Index of the first inducible subject with a valid reentry origin, else None."""
    for k, lab in enumerate(labels):
        if lab.inducible and lab.reentry_origin is not None:
            return k
    return None


def run(cfg: Config) -> dict:
    """Run the full synthetic-data benchmark and write outputs.

    Parameters
    ----------
    cfg : Config
        Full configuration. ``cfg.cohort`` sizes the synthetic cohort, ``cfg.sfi``
        parameterizes the Spectral Fragility Index field, ``cfg.labels`` drives the
        ``mock_ep`` inducibility surrogate, and ``cfg.eval`` sets the cross-validation
        and null-resampling counts. Everything is deterministic in ``cfg.seed``.

    Returns
    -------
    dict
        The metrics dictionary that is also written to
        ``<outputs_dir>/metrics.json``.
    """
    outputs_dir = cfg.outputs_dir
    os.makedirs(outputs_dir, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)

    # 1. Synthetic cohort -> conduction graphs.
    meshes: List[AtrialMesh] = make_cohort(cfg.cohort)
    graphs: List[AtrialGraph] = [mesh_to_graph(m) for m in meshes]

    # 2. mock_ep inducibility labels (development stand-in; NOT clinical POAF).
    labels: List[InducibilityLabel] = []
    for k, G in enumerate(graphs):
        labels.append(induce(G, cfg.labels, np.random.default_rng(cfg.seed + 1 + k)))

    # 3. Design matrix (baseline + spectral + SFI), labels, groups.
    X, y, groups = build_design_matrix(graphs, labels, cfg, rng)
    feature_sets = _select_feature_sets(X.columns)

    n_subjects = int(len(y))
    n_inducible = int(y.sum())
    both_classes = 0 < n_inducible < n_subjects

    # 4. Grouped vs naive cross-validated AUC per feature set.
    auc_results = nested_group_kfold_auc(
        X, y, groups, feature_sets, cfg.eval, cfg.seed
    )

    # 5. Paired DeLong: does fibrosis+SFI beat fibrosis alone?
    delong: Dict[str, float] = {}
    if both_classes:
        prob_fib = _grouped_oof_probs(
            X, y, groups, feature_sets[_FEATURE_SET_FIBROSIS], cfg)
        prob_sfi = _grouped_oof_probs(
            X, y, groups, feature_sets[_FEATURE_SET_FIBROSIS_SFI], cfg)
        try:
            delong = delong_test(y, prob_sfi, prob_fib)
        except ValueError:
            delong = {}

    # 6. Colocalization of predicted hotspot vs true reentry origin.
    coloc: Dict[str, float] = {}
    coloc_subject = _pick_colocalization_subject(graphs, labels)
    hotspot_for_fig: Optional[np.ndarray] = None
    origin_for_fig: Optional[int] = None
    mesh_for_fig: AtrialMesh = meshes[0]
    phi2_for_fig: Optional[np.ndarray] = None
    if coloc_subject is not None:
        G = graphs[coloc_subject]
        L = _combinatorial_laplacian(G)
        _, phi2 = fiedler(L)
        _, perron_v = perron(G.adjacency())
        hotspot = hotspot_map(G, phi2, perron_v)
        origin = int(labels[coloc_subject].reentry_origin)
        coloc = colocalization_vs_null(
            hotspot, origin, G.uac, cfg.eval.n_null, cfg.seed
        )
        hotspot_for_fig = hotspot
        origin_for_fig = origin
        mesh_for_fig = meshes[coloc_subject]
        phi2_for_fig = np.asarray(phi2, dtype=float).ravel()

    # Fiedler vector / spectrum for the figure subject (fall back to subject 0).
    fig_subject = coloc_subject if coloc_subject is not None else 0
    G_fig = graphs[fig_subject]
    L_fig = _combinatorial_laplacian(G_fig)
    if phi2_for_fig is None:
        _, phi2_fig = fiedler(L_fig)
        phi2_for_fig = np.asarray(phi2_fig, dtype=float).ravel()
    eig_k = int(min(8, max(3, G_fig.n_nodes)))
    eig_vals, _ = smallest_eigpairs(L_fig, k=eig_k)

    # Per-region SFI of the figure subject, for the sensitivity bars.
    field = perturbation_field(G_fig, cfg.sfi, rng)
    region_sfi = sfi_region(G_fig, phi2_for_fig, G_fig.edges, G_fig.region,
                            field["expected_dw"])
    sfi_bars = {f"region {int(r)}": float(v) for r, v in region_sfi.items()}
    if hotspot_for_fig is None:
        # No inducible subject -> show the figure subject's hotspot map with no origin.
        _, perron_v = perron(G_fig.adjacency())
        hotspot_for_fig = hotspot_map(G_fig, phi2_for_fig, perron_v)
        mesh_for_fig = meshes[fig_subject]

    # 7. Figures.
    figures = {
        "fiedler_atrium": fig_fiedler_atrium(
            meshes[fig_subject], phi2_for_fig,
            os.path.join(outputs_dir, "fig_fiedler_atrium.png")),
        "sfi_vs_origin": fig_sfi_vs_origin(
            mesh_for_fig, hotspot_for_fig, origin_for_fig,
            os.path.join(outputs_dir, "fig_sfi_vs_origin.png")),
        "eigen_spectrum": fig_eigen_spectrum(
            eig_vals, os.path.join(outputs_dir, "fig_eigen_spectrum.png")),
        "roc_panel": fig_roc_panel(
            auc_results, os.path.join(outputs_dir, "fig_roc_panel.png")),
        "sensitivity_bars": fig_sensitivity_bars(
            sfi_bars, os.path.join(outputs_dir, "fig_sensitivity_bars.png")),
    }

    # 8. Assemble + persist metrics.
    metrics = {
        "config": {
            "seed": cfg.seed,
            "n_base": cfg.cohort.n_base,
            "n_variants": cfg.cohort.n_variants,
            "n_nodes": cfg.cohort.n_nodes,
            "protocol": cfg.labels.protocol,
            "label_source": cfg.labels.source,
        },
        "cohort": {
            "n_subjects": n_subjects,
            "n_inducible": n_inducible,
            "inducible_fraction": float(n_inducible / n_subjects) if n_subjects else 0.0,
            "n_shape_families": int(len(set(groups))),
            "n_features": int(X.shape[1]),
            "both_classes_present": bool(both_classes),
        },
        "auc": auc_results,
        "delong_sfi_vs_fibrosis": delong,
        "colocalization": coloc,
        "figures": figures,
        "feature_sets": feature_sets,
        "label_disclaimer": (
            "Labels are mock_ep, an in-repo development stand-in for openCARP. "
            "They are NOT clinical post-operative AF and must not be reported as such."
        ),
    }

    metrics_path = os.path.join(outputs_dir, "metrics.json")
    with open(metrics_path, "w", encoding="utf-8") as fh:
        json.dump(_sanitize(metrics), fh, indent=2, allow_nan=False,
                  default=_json_default)

    report_path = os.path.join(outputs_dir, "results_report.md")
    _write_report(report_path, metrics)

    metrics["metrics_path"] = metrics_path
    metrics["report_path"] = report_path
    return metrics


def _sanitize(obj):
    """Recursively convert NaN/inf floats to ``None`` so metrics.json is valid JSON.

    Non-finite AUCs arise legitimately when the mock_ep labels are a single class
    (e.g. all atria inducible at high mesh resolution); emitting ``null`` keeps the
    file strictly parseable by any JSON consumer rather than the ``NaN`` token.
    """
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def _json_default(obj):
    """JSON fallback: numpy scalars/arrays -> native Python; NaN/inf -> None."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return str(obj)


def _fmt(x) -> str:
    """Format a possibly-NaN/None scalar for the markdown report."""
    if x is None:
        return "n/a"
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(xf):
        return "n/a"
    return f"{xf:.4f}"


def _write_report(path: str, metrics: dict) -> None:
    """Write the human-readable markdown results report."""
    cohort = metrics["cohort"]
    auc = metrics["auc"]
    delong = metrics.get("delong_sfi_vs_fibrosis") or {}
    coloc = metrics.get("colocalization") or {}

    lines: List[str] = []
    lines.append("# AtrialSpectralBench — synthetic-data results report")
    lines.append("")
    lines.append(
        "> **Disclaimer.** Inducibility labels are `mock_ep`, an in-repo development "
        "stand-in for openCARP. They are **not** clinical post-operative AF (POAF) "
        "and must not be interpreted or reported as such."
    )
    lines.append("")
    lines.append("## Cohort")
    lines.append("")
    lines.append(f"- Subjects: **{cohort['n_subjects']}** "
                 f"across **{cohort['n_shape_families']}** shape families")
    lines.append(f"- Inducible (mock_ep): **{cohort['n_inducible']}** "
                 f"({_fmt(cohort['inducible_fraction'])} fraction)")
    lines.append(f"- Features per subject: **{cohort['n_features']}**")
    lines.append(f"- Both label classes present: **{cohort['both_classes_present']}**")
    lines.append("")

    lines.append("## Cross-validated AUC (grouped = leak-controlled)")
    lines.append("")
    lines.append("| feature set | grouped AUC | naive AUC |")
    lines.append("| --- | --- | --- |")
    for name, res in auc.items():
        lines.append(f"| {name} | {_fmt(res.get('grouped_auc'))} "
                     f"| {_fmt(res.get('naive_auc'))} |")
    lines.append("")

    lines.append("## DeLong test — fibrosis+SFI vs fibrosis (grouped OOF)")
    lines.append("")
    if delong:
        lines.append(f"- AUC(fibrosis+SFI) = **{_fmt(delong.get('auc_a'))}**")
        lines.append(f"- AUC(fibrosis)     = **{_fmt(delong.get('auc_b'))}**")
        lines.append(f"- Delta (SFI gain)  = **{_fmt(delong.get('delta'))}**")
        lines.append(f"- z = {_fmt(delong.get('z'))}, p = {_fmt(delong.get('p'))}")
    else:
        lines.append("- Not computed (only one label class present in this cohort).")
    lines.append("")

    lines.append("## Hotspot–reentry colocalization vs spatial null")
    lines.append("")
    if coloc:
        lines.append(f"- Observed hotspot→origin UAC distance: "
                     f"**{_fmt(coloc.get('observed'))}**")
        lines.append(f"- Null 5th percentile: {_fmt(coloc.get('null_p5'))}")
        lines.append(f"- Percentile / p: {_fmt(coloc.get('p'))}")
        lines.append(f"- Passed (obs ≤ null 5th pct): **{coloc.get('passed')}**")
    else:
        lines.append("- Not computed (no inducible subject with a reentry origin).")
    lines.append("")

    lines.append("## Figures")
    lines.append("")
    for name, p in metrics["figures"].items():
        lines.append(f"- `{name}`: `{os.path.basename(p)}`")
    lines.append("")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point: ``asb run --config <yaml>`` (also ``python -m asb.pipeline``).

    Parameters
    ----------
    argv : list of str, optional
        Argument vector (defaults to ``sys.argv[1:]``).

    Returns
    -------
    int
        Process exit code (0 on success).
    """
    parser = argparse.ArgumentParser(
        prog="asb.pipeline",
        description="Run the AtrialSpectralBench synthetic-data benchmark.",
    )
    parser.add_argument(
        "--config", type=str, default=None,
        help="Path to a YAML config (defaults to built-in Config()).")
    args = parser.parse_args(argv)

    cfg = Config.load(args.config) if args.config else Config()
    metrics = run(cfg)

    print(f"Wrote {metrics['metrics_path']}")
    print(f"Wrote {metrics['report_path']}")
    for name, p in metrics["figures"].items():
        print(f"Wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
