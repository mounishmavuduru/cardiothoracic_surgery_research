r"""FitzHugh--Nagumo excitable-network instability labeller (GM4 second medium).

A genuine nonlinear excitable-media simulator on an abstract 2-D neural network:
diffusively-coupled FitzHugh--Nagumo (FHN) nodes in the **excitable** regime
(a=0.7, b=0.8, ε=0.08 — the canonical FitzHugh 1961 set). A heterogeneous **lesion**
field does two physiologically-motivated things, both literature-grounded:

1. it **reduces inter-node coupling** (already baked into the graph edge weights) →
   source--sink mismatch / propagation block (the cardiac-reentry mechanism, which
   transfers to disinhibited cortex; Huang et al. *J Neurosci* 2004);
2. it **raises local excitability** via a depolarizing bias ``I_bias·lesion`` that
   pushes lesion tissue toward the Hopf/saddle-node boundary → self-sustained
   "seizure-like" runaway activity — the Epileptor / Wilson--Cowan picture of an
   epileptic focus (Jirsa et al. *Brain* 2014; Meijer et al. *J Math Neurosci* 2015).

A stimulus kick is applied and the network is labelled **unstable** iff self-sustaining
activity persists far beyond a single transit. The instability is driven by
excitability + coupling heterogeneity, **never** by the Fiedler value :math:`\lambda_2`
— so a spectral-fragility link is a genuine finding, exactly as in the cardiac case.

Honest positioning
------------------
That :math:`\lambda_2` governs excitable-network synchronization/stability is
established prior art (Pecora--Carroll master-stability-function, PRL 1998; the Fiedler
value rises at seizure onset, Bomela et al. *Sci Rep* 2020), and the identity
:math:`\partial\lambda_2/\partial w_{ij} = (\varphi_i-\varphi_j)^2` is classical
(Ghosh--Boyd, IEEE CDC 2006). GM4's contribution is **not** those facts; it is testing
whether the per-edge fragility field **spatially localizes** a nonlinear instability
origin and whether that behaviour **transfers** from cardiac tissue to a neural network,
against a spatial null.

Model (excitable FHN, dimensionless), integrated by explicit forward Euler:

.. math::
    \dot v_i &= v_i - v_i^3/3 - w_i + I_i(t) + I_{bias}\,\ell_i + g\,(P v - v)_i, \\
    \dot w_i &= \varepsilon\,(v_i + a - b\,w_i),

with :math:`\ell` the lesion field and :math:`P=D^{-1}W` the row-normalized weighted
adjacency (so ``(Pv-v)=-L_{rw}v`` is a bounded graph diffusion whose conductances are
reduced inside the lesion).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
import scipy.sparse as sp

from asb.types import AtrialGraph, InducibilityLabel

__all__ = ["FHNConfig", "simulate_fhn", "induce_fhn", "origin_variants"]


@dataclass
class FHNConfig:
    """FitzHugh--Nagumo parameters + stimulation protocol (calibrated, excitable)."""

    a: float = 0.7
    b: float = 0.8
    eps: float = 0.08          # slow-variable rate; ε ≪ 1 => excitable, not oscillatory
    g: float = 4.0             # diffusive coupling strength
    # Lesion -> local depolarizing bias (epileptic-focus hyperexcitability). Calibrated
    # so low-lesion nets stay stable and high-lesion nets go self-sustaining.
    i_bias: float = 0.7
    dt: float = 0.02
    # Stimulus: a deterministic planar kick from the left edge (x < stim_x).
    stim_x: float = 0.12
    stim_amp: float = 1.0
    stim_onset: float = 5.0
    stim_dur: float = 3.0
    # Observation window after the stimulus (FHN time units).
    observe_after: float = 300.0
    v_thresh: float = 0.0          # "active" when v > v_thresh
    quiescent_frac: float = 0.02   # active fraction below this => quiescent
    quiescent_win: float = 20.0    # sustained quiescence => terminated (stable)
    # Unstable iff self-sustained activity persists at least this long post-stimulus.
    instability_min_time: float = 200.0


def _row_normalized_adjacency(G: AtrialGraph) -> sp.csr_matrix:
    """P = D^{-1} W (row-stochastic); coupling term is g*(P v - v)."""
    W = G.adjacency().tocsr()
    deg = np.asarray(W.sum(axis=1)).ravel()
    inv = np.zeros_like(deg)
    nz = deg > 0
    inv[nz] = 1.0 / deg[nz]
    return sp.diags(inv) @ W


def simulate_fhn(G: AtrialGraph, cfg: FHNConfig) -> Dict[str, object]:
    """Run one FHN excitable-network simulation; detect self-sustaining instability.

    Deterministic in the graph + config (a fixed planar stimulus). Returns
    ``inducible`` (True = unstable/seizure-like), ``reentry_origin`` (earliest node of
    the sustained post-stimulus activity), ``sustained_time`` and diagnostics.
    """
    n = G.n_nodes
    P = _row_normalized_adjacency(G)
    xy = np.asarray(G.coords, float)[:, :2]
    lesion = np.clip(np.asarray(G.fibrosis, float), 0.0, 1.0)
    i_bias = cfg.i_bias * lesion

    v = np.full(n, -1.2)
    w = np.full(n, -0.6)

    stim_mask = xy[:, 0] < cfg.stim_x
    if not stim_mask.any():
        stim_mask[int(np.argmin(xy[:, 0]))] = True

    total = cfg.stim_onset + cfg.stim_dur + cfg.observe_after
    dt = cfg.dt
    n_steps = int(total / dt)
    post_stim_t = cfg.stim_onset + cfg.stim_dur

    activation_last = np.full(n, -1.0)   # last up-crossing time per node
    activation_first = np.full(n, -1.0)  # first up-crossing time per node
    active_duration = np.zeros(n)        # total post-stimulus supra-threshold time
    prev_active = v > cfg.v_thresh
    quiescent_run = 0.0
    last_active_t = 0.0
    terminated: Optional[float] = None

    for step in range(n_steps):
        t = step * dt
        coup = cfg.g * (P @ v - v)
        stim = np.zeros(n)
        if cfg.stim_onset <= t < cfg.stim_onset + cfg.stim_dur:
            stim[stim_mask] = cfg.stim_amp
        v = v + dt * (v - v ** 3 / 3.0 - w + stim + coup + i_bias)
        w = w + dt * (cfg.eps * (v + cfg.a - cfg.b * w))

        active = v > cfg.v_thresh
        up = active & (~prev_active)
        activation_last[up] = t
        first = up & (activation_first < 0)
        activation_first[first] = t
        prev_active = active

        if t > post_stim_t:
            active_duration[active] += dt
            frac = float(active.mean())
            if frac >= cfg.quiescent_frac:
                last_active_t = t
                quiescent_run = 0.0
            else:
                quiescent_run += dt
                if quiescent_run >= cfg.quiescent_win and terminated is None:
                    terminated = t
                    break

    sustained_time = max(0.0, last_active_t - post_stim_t)
    inducible = bool(sustained_time >= cfg.instability_min_time)

    # Primary origin = the SUSTAINED-ACTIVITY CORE (node active longest in the
    # post-stimulus window) — the anchor of the self-sustaining instability, matching
    # the physical "where the seizure lives" intent. Alternate definitions are exposed
    # for the origin-sensitivity transparency sweep (GM4).
    origin: Optional[int] = None
    if inducible:
        origin = int(np.argmax(active_duration)) if active_duration.max() > 0 else int(np.argmax(v))
    return {
        "inducible": inducible,
        "reentry_origin": origin,
        "sustained_time": float(sustained_time),
        "terminated_at": terminated,
        "activation_last": activation_last,
        "activation_first": activation_first,
        "active_duration": active_duration,
        "post_stim_t": float(post_stim_t),
    }


def origin_variants(res: Dict[str, object]) -> Dict[str, Optional[int]]:
    """Alternate instability-origin definitions from a :func:`simulate_fhn` result.

    Used by GM4 to report how the localizer keep/delete verdicts depend on the
    (inherently ambiguous) definition of an excitable-instability "origin".
    """
    ad = np.asarray(res["active_duration"], float)
    af = np.asarray(res["activation_first"], float)
    al = np.asarray(res["activation_last"], float)
    post = float(res["post_stim_t"])
    if ad.max() <= 0:
        return {"sustained_core": None, "first_activation": None, "earliest_last": None}
    post_mask = al > post
    return {
        "sustained_core": int(np.argmax(ad)),
        "first_activation": (int(np.argmin(np.where(post_mask, af, np.inf)))
                             if post_mask.any() else None),
        "earliest_last": (int(np.argmin(np.where(post_mask, al, np.inf)))
                          if post_mask.any() else None),
    }


def induce_fhn(
    G: AtrialGraph, cfg: FHNConfig, rng: Optional[np.random.Generator] = None
) -> InducibilityLabel:
    """FHN-network instability verdict for one network (``source='fhn_network'``).

    Deterministic (a fixed planar stimulus); ``rng`` is accepted for interface
    uniformity but unused. Simulator verdict — a nonlinear dynamical instability,
    never a clinical label.
    """
    del rng
    res = simulate_fhn(G, cfg)
    return InducibilityLabel(
        inducible=bool(res["inducible"]),
        reentry_origin=res["reentry_origin"],
        protocol="fhn_planar_S1",
        source="fhn_network",
        meta={"note": "FitzHugh-Nagumo excitable-network instability; simulator verdict",
              "sustained_time": float(res["sustained_time"]),
              "terminated_at": res["terminated_at"]},
    )
