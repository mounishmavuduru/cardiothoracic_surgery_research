"""Streamlit dashboard for AtrialSpectralBench (import-guarded).

Interactive explorer for the synthetic cohort: pick an atrium, view its Fiedler
mode (:math:`\\varphi_2`) coloring, its SFI hotspot map, and its per-subject
``mock_ep`` inducibility + spectral-fragility risk readout.

Import safety
-------------
Streamlit is an **optional** dependency. This module must import cleanly even when
Streamlit is absent, so the ``streamlit`` import is wrapped in a guard: if it is
missing, :data:`STREAMLIT_AVAILABLE` is ``False`` and :func:`main` raises a helpful
error instead of crashing at import time. The core ``asb`` package therefore never
depends on Streamlit.

Run with::

    streamlit run src/asb/dashboard/app.py

.. warning::

    The inducibility labels shown here are ``mock_ep``, a development stand-in for
    openCARP — **not** clinical post-operative AF.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import scipy.sparse as sp

try:  # Import guard: the core package must work without streamlit installed.
    import streamlit as st

    STREAMLIT_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only in streamlit-free envs.
    st = None  # type: ignore[assignment]
    STREAMLIT_AVAILABLE = False

from asb.config import Config
from asb.labels.mock_ep import induce
from asb.sfi import hotspot_map, perturbation_field, sfi_region
from asb.spectral import fiedler, perron
from asb.substrate.mesh import mesh_to_graph
from asb.substrate.synthetic import make_cohort
from asb.types import AtrialGraph

__all__ = ["STREAMLIT_AVAILABLE", "main"]


def _combinatorial_laplacian(G: AtrialGraph) -> sp.csr_matrix:
    """Combinatorial Laplacian ``L = D - W`` from a graph's adjacency."""
    W = G.adjacency()
    d = np.asarray(W.sum(axis=1)).ravel()
    return sp.csr_matrix(sp.diags(d) - W)


def _project_2d(points: np.ndarray) -> np.ndarray:
    """PCA-project (n, 3) points to 2D for scatter display."""
    pts = np.asarray(points, dtype=float)
    centered = pts - pts.mean(axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return centered @ vt[:2].T


def _build_cohort(cfg: Config):
    """Generate the synthetic cohort meshes and graphs (cached in the app)."""
    meshes = make_cohort(cfg.cohort)
    graphs = [mesh_to_graph(m) for m in meshes]
    return meshes, graphs


def main(argv: Optional[list] = None) -> int:
    """Render the Streamlit dashboard.

    Parameters
    ----------
    argv : list, optional
        Unused; present for a uniform entry-point signature.

    Returns
    -------
    int
        ``0`` on success.

    Raises
    ------
    RuntimeError
        If Streamlit is not installed (install the ``dashboard`` extra).
    """
    if not STREAMLIT_AVAILABLE:
        raise RuntimeError(
            "Streamlit is not installed. Install the dashboard extra with "
            "`pip install -e \".[dashboard]\"` and run "
            "`streamlit run src/asb/dashboard/app.py`."
        )

    import matplotlib.pyplot as plt

    st.set_page_config(page_title="AtrialSpectralBench", layout="wide")
    st.title("AtrialSpectralBench — spectral fragility explorer")
    st.warning(
        "Inducibility labels are `mock_ep`, a development stand-in for openCARP. "
        "They are NOT clinical post-operative AF."
    )

    # --- Sidebar controls -------------------------------------------------- #
    st.sidebar.header("Cohort")
    n_base = st.sidebar.slider("Base anatomies", 1, 6, 2)
    n_variants = st.sidebar.slider("Variants per anatomy", 1, 6, 2)
    n_nodes = st.sidebar.select_slider(
        "Approx. nodes", options=[42, 162, 642, 1200], value=162)
    seed = st.sidebar.number_input("Seed", value=0, step=1)

    cfg = Config()
    cfg.cohort.n_base = int(n_base)
    cfg.cohort.n_variants = int(n_variants)
    cfg.cohort.n_nodes = int(n_nodes)
    cfg.cohort.seed = int(seed)
    cfg.seed = int(seed)

    cache = st.cache_data(_build_cohort) if hasattr(st, "cache_data") else _build_cohort
    meshes, graphs = cache(cfg)

    subject = st.sidebar.slider("Subject", 0, len(graphs) - 1, 0)
    mesh = meshes[subject]
    G = graphs[subject]

    # --- Compute spectral + SFI + label for the chosen subject ------------- #
    L = _combinatorial_laplacian(G)
    lam2, phi2 = fiedler(L)
    _, perron_v = perron(G.adjacency())
    hotspot = hotspot_map(G, phi2, perron_v)

    field = perturbation_field(G, cfg.sfi, np.random.default_rng(cfg.seed))
    region_sfi = sfi_region(G, phi2, G.edges, G.region, field["expected_dw"])
    total_sfi = float(sum(region_sfi.values()))

    label = induce(G, cfg.labels, np.random.default_rng(cfg.seed + 1 + subject))

    xy = _project_2d(mesh.points)

    # --- Readout ----------------------------------------------------------- #
    c1, c2, c3 = st.columns(3)
    c1.metric("Shape family", str(G.shape_family))
    c1.metric("Nodes / edges", f"{G.n_nodes} / {G.n_edges}")
    c2.metric("Algebraic connectivity λ₂", f"{lam2:.4g}")
    c2.metric("Total SFI (Σ expected Δλ₂)", f"{total_sfi:.4g}")
    c3.metric("mock_ep inducible", str(label.inducible))
    c3.metric("Reentry origin node", str(label.reentry_origin))

    # --- Figures ----------------------------------------------------------- #
    fcol1, fcol2 = st.columns(2)

    with fcol1:
        st.subheader("Fiedler mode φ₂")
        fig1, ax1 = plt.subplots(figsize=(4.5, 4))
        sc1 = ax1.scatter(xy[:, 0], xy[:, 1], c=phi2, cmap="RdBu_r", s=10)
        fig1.colorbar(sc1, ax=ax1, label="φ₂")
        ax1.set_aspect("equal", adjustable="datalim")
        st.pyplot(fig1)
        plt.close(fig1)

    with fcol2:
        st.subheader("SFI hotspot map")
        fig2, ax2 = plt.subplots(figsize=(4.5, 4))
        sc2 = ax2.scatter(xy[:, 0], xy[:, 1], c=hotspot, cmap="magma", s=10)
        fig2.colorbar(sc2, ax=ax2, label="hotspot score")
        if label.reentry_origin is not None:
            o = int(label.reentry_origin)
            ax2.scatter(xy[o, 0], xy[o, 1], marker="*", s=260, c="cyan",
                        edgecolors="black", label="reentry origin")
            ax2.legend(loc="best")
        ax2.set_aspect("equal", adjustable="datalim")
        st.pyplot(fig2)
        plt.close(fig2)

    st.subheader("Per-region Spectral Fragility Index")
    st.bar_chart({str(int(r)): v for r, v in region_sfi.items()})

    return 0


if __name__ == "__main__":  # pragma: no cover - only under `streamlit run`.
    main()
