r"""Kuramoto oscillator-network instability labeller (GM4 **third** medium).

The first two media in GM4 are *excitable* (cardiac monodomain, FitzHugh--Nagumo):
instability = a self-sustaining propagating pulse (reentry / seizure). This module
adds a genuinely *different dynamical class* — **phase oscillators** (Kuramoto 1975) —
where instability = **loss of synchrony** (a persistent phase-slipping cluster). Testing
the same spectral-fragility calculus across excitable *and* oscillatory dynamics is a
strictly stronger transfer claim than two excitable media.

Why the Fiedler pair is the right object here (mechanistic, not incidental)
--------------------------------------------------------------------------
Linear stability of the fully phase-locked Kuramoto state is governed by the graph
**Laplacian spectrum** (the master-stability-function / Barahona--Pecora synchronizability
picture): the algebraic connectivity :math:`\lambda_2` sets the slowest-decaying transverse
mode, and when coupling weakens the sync manifold destabilizes **along the Fiedler mode**
:math:`\varphi_2`. So a network that is marginally coupled splits into weakly-connected
groups across the :math:`\varphi_2` sign structure, and the desynchronization nucleates
where :math:`|\nabla\varphi_2|` is large (the connectivity bottleneck). The localizer test
— does :math:`|\nabla\varphi_2|` find the desync origin above a spatial null — is therefore
motivated by the *same* spectral mechanism as cardiac reentry, on completely different
dynamics.

Honest positioning (identical to GM4's other media)
---------------------------------------------------
That :math:`\lambda_2` governs synchronizability is established prior art
(Pecora--Carroll PRL 1998; Barahona--Pecora PRL 2002), and the derivative identity
:math:`\partial\lambda_2/\partial w_{ij}=(\varphi_i-\varphi_j)^2` is classical
(Ghosh--Boyd CDC 2006). The contribution is **not** those facts; it is testing whether the
per-edge fragility field spatially localizes the instability origin and whether that
behaviour transfers to a third, non-excitable medium, against a spatial null.

Model (Kuramoto on the weighted conduction graph), explicit forward Euler:

.. math::
    \dot\theta_i = \omega_i + K \sum_{j} W_{ij}\,\sin(\theta_j-\theta_i),

with :math:`W_{ij}` the (lesion-reduced) edge conductances and :math:`\omega_i` heterogeneous
natural frequencies. The lesion field lowers coupling in patches, so those patches lose lock
first — the structural analogue of the low-coupling fibrosis substrate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np

from asb.types import AtrialGraph, InducibilityLabel

__all__ = ["KuramotoConfig", "simulate_kuramoto", "induce_kuramoto", "origin_variants"]


@dataclass
class KuramotoConfig:
    """Kuramoto parameters (calibrated so the lesion-burden sweep yields both classes)."""

    K: float = 3.0                 # global coupling gain (calibrated: ~50% desync)
    sigma_omega: float = 0.8       # natural-frequency spread (std of omega_i)
    dt: float = 0.02
    t_total: float = 120.0         # total integration time
    settle_frac: float = 0.5       # fraction of the run discarded as transient
    r_sync_thresh: float = 0.90    # order parameter above which the net is "synchronized"
    slip_thresh: float = 0.25      # per-node frequency-detuning that counts as "slipping"
    freq_seed_offset: int = 10000  # deterministic omega / initial-phase seed offset


def _node_mean_frequency(theta_hist: np.ndarray, dt: float) -> np.ndarray:
    """Time-averaged phase velocity per node over a phase history (unwrapped)."""
    unwrapped = np.unwrap(theta_hist, axis=0)
    return (unwrapped[-1] - unwrapped[0]) / (dt * (theta_hist.shape[0] - 1))


def simulate_kuramoto(G: AtrialGraph, cfg: KuramotoConfig) -> Dict[str, object]:
    """Run one Kuramoto simulation; detect a persistent desynchronization cluster.

    Deterministic in the graph + config (frequencies/initial phases seeded from the
    graph's ``meta['seed']``). Returns ``inducible`` (True = fails to synchronize / has a
    persistent slipping cluster), ``reentry_origin`` (the most-detuned node = desync core),
    and diagnostics for the origin-sensitivity sweep.
    """
    n = G.n_nodes
    edges = np.asarray(G.edges, np.int64)
    w = np.asarray(G.weights, float)
    ii, jj = edges[:, 0], edges[:, 1]

    seed = int(G.meta.get("seed", 0)) + cfg.freq_seed_offset
    rng = np.random.default_rng(seed)
    omega = cfg.sigma_omega * rng.standard_normal(n)
    theta = rng.uniform(-np.pi, np.pi, n)

    n_steps = int(cfg.t_total / cfg.dt)
    settle = int(cfg.settle_frac * n_steps)
    # store phase history only for the post-settle window (for mean-frequency estimate)
    hist_len = n_steps - settle + 1
    theta_hist = np.empty((hist_len, n), dtype=float)
    r_series = np.empty(n_steps, dtype=float)
    K, dt = cfg.K, cfg.dt

    h = 0
    for step in range(n_steps):
        dtheta = theta[jj] - theta[ii]
        s = w * np.sin(dtheta)
        coup = np.zeros(n)
        np.add.at(coup, ii, s)
        np.add.at(coup, jj, -s)
        theta = theta + dt * (omega + K * coup)
        r_series[step] = np.abs(np.exp(1j * theta).mean())
        if step >= settle:
            theta_hist[h] = theta
            h += 1

    r_final = float(r_series[settle:].mean())
    node_freq = _node_mean_frequency(theta_hist[:h], dt)
    detune = np.abs(node_freq - np.median(node_freq))
    max_slip = float(detune.max())
    slip_frac = float(np.mean(detune > cfg.slip_thresh))

    # Unstable iff the network does not lock into global synchrony AND a genuine
    # slipping cluster exists (rules out trivial near-sync noise).
    inducible = bool(r_final < cfg.r_sync_thresh and max_slip > cfg.slip_thresh)

    # Primary origin = the desync core: the node whose mean frequency is most detuned
    # from the bulk (the anchor of the phase-slipping cluster), analogous to the
    # sustained-activity core in the excitable media.
    origin: Optional[int] = int(np.argmax(detune)) if inducible else None
    return {
        "inducible": inducible,
        "reentry_origin": origin,
        "r_final": r_final,
        "max_slip": max_slip,
        "slip_frac": slip_frac,
        "detune": detune,
        "node_freq": node_freq,
    }


def origin_variants(res: Dict[str, object]) -> Dict[str, Optional[int]]:
    """Alternate desync-origin definitions for the origin-sensitivity transparency sweep."""
    detune = np.asarray(res["detune"], float)
    if detune.max() <= 0:
        return {"detune_core": None, "second_detune": None}
    order = np.argsort(detune)[::-1]
    return {
        "detune_core": int(order[0]),
        "second_detune": int(order[1]) if detune.size > 1 else int(order[0]),
    }


def induce_kuramoto(
    G: AtrialGraph, cfg: KuramotoConfig, rng: Optional[np.random.Generator] = None
) -> InducibilityLabel:
    """Kuramoto desynchronization verdict for one network (``source='kuramoto_network'``).

    Deterministic (frequencies/phases seeded from the graph); ``rng`` accepted for
    interface uniformity but unused. Simulator verdict — a dynamical desync instability,
    never a clinical label.
    """
    del rng
    res = simulate_kuramoto(G, cfg)
    return InducibilityLabel(
        inducible=bool(res["inducible"]),
        reentry_origin=res["reentry_origin"],
        protocol="kuramoto_random_phase",
        source="kuramoto_network",
        meta={"note": "Kuramoto oscillator-network desynchronization; simulator verdict",
              "r_final": float(res["r_final"]), "max_slip": float(res["max_slip"]),
              "slip_frac": float(res["slip_frac"])},
    )
