"""Matplotlib figure writers for AtrialSpectralBench.

Every function is a **pure renderer**: it takes already-computed data plus an output
path, writes a PNG, and returns the path. There is no compute here (no eigensolves,
no SFI) — all numbers are produced upstream in the compute modules and merely drawn.

The module uses matplotlib's non-interactive ``Agg`` backend so figures render
headlessly (no display needed) inside the pipeline and CI.
"""
from __future__ import annotations

import os
from typing import Dict, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")  # headless backend; must precede pyplot import.

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from asb.types import AtrialMesh  # noqa: E402

__all__ = [
    "fig_fiedler_atrium",
    "fig_sfi_vs_origin",
    "fig_eigen_spectrum",
    "fig_roc_panel",
    "fig_sensitivity_bars",
]


def _ensure_parent(out: str) -> None:
    """Create the parent directory of ``out`` if needed."""
    parent = os.path.dirname(os.path.abspath(out))
    os.makedirs(parent, exist_ok=True)


def _2d_projection(points: np.ndarray) -> np.ndarray:
    """Project (n, 3) points to a 2D plane via the two dominant PCA axes.

    Gives a stable, orientation-independent flat view of a closed surface for
    scatter rendering (no external 3D dependencies).
    """
    pts = np.asarray(points, dtype=float)
    if pts.shape[1] == 2:
        return pts
    centered = pts - pts.mean(axis=0, keepdims=True)
    # Right singular vectors are the principal axes; take the top two.
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return centered @ vt[:2].T


def fig_fiedler_atrium(mesh: AtrialMesh, phi2: Sequence[float], out: str) -> str:
    """Scatter the atrial surface colored by the Fiedler vector ``phi2``.

    Parameters
    ----------
    mesh : AtrialMesh
        Surface providing vertex ``points``.
    phi2 : (n,) array_like
        Fiedler vector value per vertex.
    out : str
        Output PNG path.

    Returns
    -------
    str
        ``out``.
    """
    _ensure_parent(out)
    xy = _2d_projection(mesh.points)
    phi2 = np.asarray(phi2, dtype=float).ravel()

    fig, ax = plt.subplots(figsize=(5, 4.2))
    sc = ax.scatter(xy[:, 0], xy[:, 1], c=phi2, cmap="RdBu_r", s=12, edgecolors="none")
    fig.colorbar(sc, ax=ax, label=r"$\varphi_2$ (Fiedler vector)")
    ax.set_title("Fiedler mode over atrial surface")
    ax.set_xlabel("PCA axis 1")
    ax.set_ylabel("PCA axis 2")
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def fig_sfi_vs_origin(
    mesh: AtrialMesh, sfi: Sequence[float], origin, out: str
) -> str:
    """Scatter per-node SFI / hotspot scores with the true reentry origin marked.

    Parameters
    ----------
    mesh : AtrialMesh
        Surface providing vertex ``points``.
    sfi : (n,) array_like
        Per-node SFI hotspot score (larger = more fragile / more likely initiation).
    origin : int or None
        Node index of the true reentry origin (marked with a star), or ``None``.
    out : str
        Output PNG path.

    Returns
    -------
    str
        ``out``.
    """
    _ensure_parent(out)
    xy = _2d_projection(mesh.points)
    sfi = np.asarray(sfi, dtype=float).ravel()

    fig, ax = plt.subplots(figsize=(5, 4.2))
    sc = ax.scatter(xy[:, 0], xy[:, 1], c=sfi, cmap="magma", s=12, edgecolors="none")
    fig.colorbar(sc, ax=ax, label="SFI hotspot score")
    if origin is not None and 0 <= int(origin) < xy.shape[0]:
        o = int(origin)
        ax.scatter(
            xy[o, 0], xy[o, 1], marker="*", s=320, c="cyan",
            edgecolors="black", linewidths=1.0, label="reentry origin", zorder=5,
        )
        ax.legend(loc="best")
    ax.set_title("SFI hotspot map vs reentry origin")
    ax.set_xlabel("PCA axis 1")
    ax.set_ylabel("PCA axis 2")
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def fig_eigen_spectrum(vals: Sequence[float], out: str) -> str:
    """Stem plot of the low Laplacian eigenvalues.

    Parameters
    ----------
    vals : (k,) array_like
        Eigenvalues in ascending order.
    out : str
        Output PNG path.

    Returns
    -------
    str
        ``out``.
    """
    _ensure_parent(out)
    vals = np.asarray(vals, dtype=float).ravel()
    idx = np.arange(vals.size)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.stem(idx, vals, basefmt=" ")
    ax.set_title("Low Laplacian spectrum")
    ax.set_xlabel("mode index $k$")
    ax.set_ylabel(r"eigenvalue $\lambda_k$")
    if vals.size >= 3:
        # Annotate the algebraic connectivity and the spectral gap.
        ax.axhline(vals[1], color="tab:green", ls="--", lw=1,
                   label=rf"$\lambda_2={vals[1]:.3g}$")
        ax.axhline(vals[2], color="tab:orange", ls=":", lw=1,
                   label=rf"$\lambda_3={vals[2]:.3g}$")
        ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def fig_roc_panel(results: Mapping[str, Mapping[str, float]], out: str) -> str:
    """Grouped bar chart of grouped vs naive AUC per feature set.

    Parameters
    ----------
    results : dict
        ``{feature_set_name: {"grouped_auc": float, "naive_auc": float}}`` as
        returned by :func:`asb.evaluation.nested_group_kfold_auc`.
    out : str
        Output PNG path.

    Returns
    -------
    str
        ``out``.
    """
    _ensure_parent(out)
    names = list(results.keys())
    grouped = [float(results[n].get("grouped_auc", np.nan)) for n in names]
    naive = [float(results[n].get("naive_auc", np.nan)) for n in names]

    x = np.arange(len(names))
    width = 0.38

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(x - width / 2, grouped, width, label="grouped (leak-controlled)",
           color="tab:blue")
    ax.bar(x + width / 2, naive, width, label="naive (leaky)", color="tab:red",
           alpha=0.7)
    ax.axhline(0.5, color="gray", ls="--", lw=1, label="chance")
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=10, ha="right")
    ax.set_ylabel("ROC-AUC")
    ax.set_title("Cross-validated AUC by feature set")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def fig_sensitivity_bars(sobol: Mapping[str, float], out: str) -> str:
    """Horizontal bar chart of a named sensitivity / importance mapping.

    Parameters
    ----------
    sobol : dict of str -> float
        Named sensitivity scores (e.g. per-region SFI, or feature importances).
        Sorted by magnitude for display.
    out : str
        Output PNG path.

    Returns
    -------
    str
        ``out``.
    """
    _ensure_parent(out)
    items = sorted(sobol.items(), key=lambda kv: float(kv[1]))
    labels = [str(k) for k, _ in items]
    vals = [float(v) for _, v in items]

    fig, ax = plt.subplots(figsize=(6, max(3, 0.4 * len(labels) + 1)))
    ax.barh(np.arange(len(labels)), vals, color="tab:purple")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("sensitivity / importance")
    ax.set_title("Feature sensitivity")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out
