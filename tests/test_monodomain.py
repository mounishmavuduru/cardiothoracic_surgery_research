"""Tests for the monodomain Mitchell--Schaeffer inducibility labeller.

Kept fast: small strips / short observation windows. Verifies the physics gates
(planar CV in range, APD shortens with fibrosis, determinism, operator shapes)
that the E0 calculation-verification relies on.
"""
from __future__ import annotations

import numpy as np

from asb.labels.monodomain import (
    MonodomainConfig,
    cotangent_operator,
    induce_monodomain,
    measure_planar_cv,
    simulate_monodomain,
    single_cell_apd,
)
from asb.substrate.synthetic import make_base_atrium, paint_fibrosis


def _small_cfg(**kw):
    base = dict(d0=0.2, tau_close=110.0, stim_amp=0.15, dt=0.05,
                observe_after=200.0, reentry_min_ms=600.0)
    base.update(kw)
    return MonodomainConfig(**base)


def test_planar_cv_physiological():
    cv = measure_planar_cv(_small_cfg())["cv_m_per_s"]
    # Human LA CV band 0.3-1.2 m/s.
    assert 0.2 < cv < 1.4


def test_apd_shortens_with_fibrosis():
    cfg = _small_cfg()
    apd_healthy = single_cell_apd(cfg)["apd90"]
    apd_fibrotic = single_cell_apd(cfg, tau_close=cfg.tau_close * 0.5)["apd90"]
    assert apd_healthy > apd_fibrotic > 0
    # Atrial APD is a few hundred ms.
    assert 100 < apd_healthy < 350


def test_cotangent_operator_shapes():
    mesh = paint_fibrosis(make_base_atrium(seed=0, n_nodes=500), seed=1, burden=0.2)
    from asb.substrate.roney import coarsen_mesh
    coarse = coarsen_mesh(mesh, 300)
    cfg = _small_cfg()
    edges, T, mass, length = cotangent_operator(coarse, cfg)
    assert edges.shape[1] == 2
    assert np.all(T >= 0.0)          # clamped non-negative (M-matrix)
    assert np.all(mass > 0.0)        # positive lumped mass
    assert np.all(length > 0.0)


def test_planar_wave_propagates():
    # A single beat's activation times increase away from the paced end.
    cfg = _small_cfg(protocol="S1S2", n_s1=1, observe_after=120.0,
                     s1_cycle_length=400.0, stim_radius_mm=1.0)
    res = measure_planar_cv(cfg)
    assert res["t2"] > res["t1"] > 0   # farther probe activates later


def test_induce_deterministic():
    mesh = paint_fibrosis(make_base_atrium(seed=3, n_nodes=800), seed=4, burden=0.4)
    from asb.substrate.roney import coarsen_mesh
    coarse = coarsen_mesh(mesh, 500)
    cfg = _small_cfg(observe_after=700.0, n_pacing_sites=1)
    a = induce_monodomain(coarse, cfg, np.random.default_rng(0), burst_cls=(150.0,))
    b = induce_monodomain(coarse, cfg, np.random.default_rng(0), burst_cls=(150.0,))
    assert a.inducible == b.inducible
    assert a.source == "monodomain_ms"
    assert abs(a.meta["best_sustained_ms"] - b.meta["best_sustained_ms"]) < 1e-6


def test_simulate_stable_bounded():
    # Implicit diffusion keeps V bounded (no CFL blow-up) on a coarse mesh.
    mesh = paint_fibrosis(make_base_atrium(seed=5, n_nodes=600), seed=6, burden=0.3)
    from asb.substrate.roney import coarsen_mesh
    coarse = coarsen_mesh(mesh, 400)
    cfg = _small_cfg(protocol="burst", burst_cycle_length=150.0, n_burst=3,
                     observe_after=200.0)
    res = simulate_monodomain(coarse, cfg, np.array([0]), record_activation=True)
    assert float(np.max(res["maxV_trace"])) < 1.2  # no numerical blow-up
