"""Assemble the per-subject design matrix for AtrialSpectralBench.

This is the top-of-stack feature layer: it turns each :class:`asb.types.AtrialGraph`
(one atrium) into a flat, *namespaced* feature vector and stacks the cohort into a
:class:`pandas.DataFrame` design matrix ``X`` with a matching binary label vector
``y`` and a ``groups`` list (the ``shape_family`` of each subject, used for grouped
cross-validation so statistically-resampled near-duplicate atria never straddle a
train/test split).

Feature namespaces (column prefixes)
------------------------------------
``fibhet_``
    Fibrosis-*heterogeneity* baseline descriptors (mean burden, spatial entropy,
    patch size). These are the "fibrosis" competitor feature set.
``conn_``
    Classical graph-connectivity baselines (global min-cut, percolation threshold,
    algebraic connectivity :math:`\\lambda_2` alone).
``spec_``
    Global spectral descriptors of the Laplacian (gap, spectral entropy, Cheeger
    estimate, Perron localization, heat-kernel trace, Fiedler-vector statistics).
``sfi_``
    Summary statistics of the per-region **Spectral Fragility Index** — the protected
    novel seed. These are the columns the "fibrosis+SFI" set adds on top of ``fibhet_``.

The ``y`` labels come from :mod:`asb.labels.mock_ep`, a development stand-in for
openCARP. They are **never** clinical post-operative AF and must not be reported as
such; ``build_design_matrix`` merely consumes whatever :class:`InducibilityLabel`
objects it is given.

All compute is pure (no I/O). The only stochastic dependency is the SFI perturbation
field, which is driven entirely by the caller-supplied ``rng``.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from asb.baselines import baseline_feature_vector
from asb.config import Config, SFIConfig
from asb.sfi import edge_fragility, perturbation_field, sfi_region
from asb.spectral import fiedler, inverse_participation_ratio, spectral_features
from asb.types import AtrialGraph, InducibilityLabel

__all__ = [
    "spectral_feature_vector",
    "sfi_feature_vector",
    "subject_features",
    "build_design_matrix",
    "FIBROSIS_HETEROGENEITY_KEYS",
    "SFI_KEY_PREFIX",
]

#: Baseline keys treated as "fibrosis-heterogeneity" competitors (namespaced with the
#: ``fibhet_`` prefix in the design matrix).
FIBROSIS_HETEROGENEITY_KEYS = (
    "fibrosis_burden",
    "fibrosis_spatial_entropy",
    "fibrosis_patch_size",
)
#: Column prefix marking the Spectral Fragility Index feature block.
SFI_KEY_PREFIX = "sfi_"

# Number of top-ranked per-region SFI values exposed as ordered columns.
_SFI_TOP_K = 3


def _combinatorial_laplacian(G: AtrialGraph) -> sp.csr_matrix:
    """Combinatorial Laplacian ``L = D - W`` built from a graph's adjacency."""
    W = G.adjacency()
    d = np.asarray(W.sum(axis=1)).ravel()
    return sp.csr_matrix(sp.diags(d) - W)


def spectral_feature_vector(G: AtrialGraph) -> Dict[str, float]:
    """Global spectral descriptors of one atrial conduction graph.

    Combines :func:`asb.spectral.spectral_features` (spectral gap, spectral entropy,
    Cheeger estimate, near-zero-mode count, Perron localization, heat-kernel trace,
    Fiedler value) with a few statistics of the Fiedler vector itself (its inverse
    participation ratio and absolute-value spread), which summarize how localized the
    slowest global conduction mode is.

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph.

    Returns
    -------
    dict of str -> float
        Named spectral features (un-prefixed; :func:`subject_features` adds the
        ``spec_`` namespace).
    """
    feats = dict(spectral_features(G))

    L = _combinatorial_laplacian(G)
    if G.n_nodes >= 2:
        _, phi2 = fiedler(L)
        phi2 = np.asarray(phi2, dtype=float).ravel()
        feats["fiedler_ipr"] = inverse_participation_ratio(phi2)
        feats["fiedler_abs_mean"] = float(np.mean(np.abs(phi2)))
        feats["fiedler_abs_max"] = float(np.max(np.abs(phi2)))
    else:
        feats["fiedler_ipr"] = 0.0
        feats["fiedler_abs_mean"] = 0.0
        feats["fiedler_abs_max"] = 0.0
    return feats


def sfi_feature_vector(
    G: AtrialGraph, cfg: SFIConfig, rng: np.random.Generator
) -> Dict[str, float]:
    """Per-region Spectral Fragility Index summary statistics for one atrium.

    Computes the Fiedler pair, the fibrosis-weighted diffuse expected uncoupling
    field (:func:`asb.sfi.perturbation_field`), and the first-order analytic
    per-region SFI (:func:`asb.sfi.sfi_region`). The per-region SFI values (the
    expected drop in algebraic connectivity contributed by each region's share of the
    diffuse stress) are then reduced to interpretable scalars: the whole-atrium total,
    the max / mean / std across regions, and the top-``k`` region values in
    descending order.

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph.
    cfg : SFIConfig
        SFI field parameters (mean fractional drop, coefficient of variation,
        fibrosis coupling). Only the analytic expectation is used here — no
        Monte-Carlo — so this stays fast for whole-cohort feature extraction.
    rng : numpy.random.Generator
        Passed to :func:`asb.sfi.perturbation_field` for interface uniformity; the
        analytic expectation it returns is deterministic in ``G``/``cfg``.

    Returns
    -------
    dict of str -> float
        Named SFI summary features (un-prefixed; :func:`subject_features` adds the
        ``sfi_`` namespace):

        - ``total``           : sum of the per-edge expected :math:`\\lambda_2` drop.
        - ``region_max``      : largest per-region SFI.
        - ``region_mean``     : mean per-region SFI.
        - ``region_std``      : standard deviation across regions.
        - ``top1`` .. ``topK``: the ``k`` largest per-region SFI values, descending.
        - ``edge_frag_max``   : largest single-edge Fiedler fragility.
        - ``edge_frag_mean``  : mean single-edge Fiedler fragility.
    """
    n = G.n_nodes
    if n < 2 or G.n_edges == 0:
        out = {
            "total": 0.0,
            "region_max": 0.0,
            "region_mean": 0.0,
            "region_std": 0.0,
            "edge_frag_max": 0.0,
            "edge_frag_mean": 0.0,
        }
        for r in range(_SFI_TOP_K):
            out[f"top{r + 1}"] = 0.0
        return out

    L = _combinatorial_laplacian(G)
    _, phi2 = fiedler(L)
    phi2 = np.asarray(phi2, dtype=float).ravel()

    field = perturbation_field(G, cfg, rng)
    expected_dw = field["expected_dw"]

    region_sfi = sfi_region(G, phi2, G.edges, G.region, expected_dw)
    values = np.asarray(list(region_sfi.values()), dtype=float)

    frag = edge_fragility(phi2, G.edges)
    total = float(np.sum(expected_dw * frag))

    out: Dict[str, float] = {
        "total": total,
        "region_max": float(values.max()) if values.size else 0.0,
        "region_mean": float(values.mean()) if values.size else 0.0,
        "region_std": float(values.std()) if values.size else 0.0,
        "edge_frag_max": float(frag.max()) if frag.size else 0.0,
        "edge_frag_mean": float(frag.mean()) if frag.size else 0.0,
    }

    ranked = np.sort(values)[::-1] if values.size else np.array([], dtype=float)
    for r in range(_SFI_TOP_K):
        out[f"top{r + 1}"] = float(ranked[r]) if r < ranked.size else 0.0
    return out


def subject_features(
    G: AtrialGraph, cfg: Config, rng: np.random.Generator
) -> Dict[str, float]:
    """Merge baseline, spectral and SFI features into one namespaced vector.

    The three sub-vectors are merged with disjoint column prefixes so downstream
    evaluation can select any feature block by name:

    - baseline fibrosis-heterogeneity descriptors -> ``fibhet_``
    - baseline connectivity descriptors            -> ``conn_``
    - spectral descriptors                         -> ``spec_``
    - Spectral Fragility Index summaries           -> ``sfi_``

    Parameters
    ----------
    G : AtrialGraph
        Weighted conduction graph for one subject.
    cfg : Config
        Full configuration; only ``cfg.sfi`` is used here (for the SFI field).
    rng : numpy.random.Generator
        RNG threaded into the SFI perturbation field.

    Returns
    -------
    dict of str -> float
        Flat, namespaced feature dictionary for the subject.
    """
    baseline = baseline_feature_vector(G)
    feats: Dict[str, float] = {}
    for key, val in baseline.items():
        prefix = "fibhet_" if key in FIBROSIS_HETEROGENEITY_KEYS else "conn_"
        feats[f"{prefix}{key}"] = float(val)

    for key, val in spectral_feature_vector(G).items():
        feats[f"spec_{key}"] = float(val)

    for key, val in sfi_feature_vector(G, cfg.sfi, rng).items():
        feats[f"{SFI_KEY_PREFIX}{key}"] = float(val)

    return feats


def build_design_matrix(
    graphs: Sequence[AtrialGraph],
    labels: Sequence[InducibilityLabel],
    cfg: Config,
    rng: np.random.Generator,
) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """Stack a cohort into a design matrix, label vector and group list.

    Parameters
    ----------
    graphs : sequence of AtrialGraph
        One conduction graph per subject.
    labels : sequence of InducibilityLabel
        Matching inducibility labels (``mock_ep`` in development — never clinical
        POAF). ``y[i] = int(labels[i].inducible)``.
    cfg : Config
        Full configuration (``cfg.sfi`` drives the SFI feature block).
    rng : numpy.random.Generator
        RNG threaded per subject into the SFI perturbation field (a distinct child
        stream per subject keeps extraction order-independent and reproducible).

    Returns
    -------
    X : pandas.DataFrame
        Design matrix, one row per subject, namespaced feature columns (sorted for a
        stable column order).
    y : numpy.ndarray
        Binary inducibility flags (int), shape ``(n_subjects,)``.
    groups : list of str
        ``shape_family`` of each subject, for :class:`sklearn.model_selection.GroupKFold`.
    """
    if len(graphs) != len(labels):
        raise ValueError(
            f"graphs and labels length mismatch: {len(graphs)} vs {len(labels)}"
        )

    # Deterministic per-subject child RNGs so extraction is independent of order.
    seeds = rng.integers(0, 2**32 - 1, size=len(graphs))

    rows: List[Dict[str, float]] = []
    y = np.empty(len(graphs), dtype=int)
    groups: List[str] = []
    for k, (G, lab) in enumerate(zip(graphs, labels)):
        child = np.random.default_rng(int(seeds[k]))
        rows.append(subject_features(G, cfg, child))
        y[k] = int(bool(lab.inducible))
        groups.append(str(G.shape_family))

    X = pd.DataFrame(rows)
    # Stable, sorted column order; fill any (defensive) missing keys with 0.
    X = X.reindex(sorted(X.columns), axis=1).fillna(0.0)
    return X, y, groups
