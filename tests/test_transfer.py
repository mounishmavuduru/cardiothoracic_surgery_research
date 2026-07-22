r"""Tests for the GM4 second excitable medium (network generator + FHN labeller).

Verifies (a) the network is a valid weighted graph the spectral/SFI stack accepts,
(b) the FHN labeller is deterministic and its instability grows with lesion burden
(the calibrated substrate → instability link), and (c) — the load-bearing GM4 claim —
that the cardiac SFI/spectral code runs on the neural network **verbatim**.
"""
from __future__ import annotations

import numpy as np

from asb.transfer.network import NetworkConfig, make_excitable_network, make_network_cohort
from asb.transfer.fhn import FHNConfig, simulate_fhn, induce_fhn


def _small_cfg():
    return NetworkConfig(n_nodes=400, radius=0.10, n_region_grid=3)


def test_network_is_valid_weighted_graph():
    G = make_excitable_network(1, _small_cfg())
    assert G.n_nodes > 0 and G.n_edges > 0
    # Connected (single component): every node reachable via the edge set.
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    ncomp, _ = connected_components(G.adjacency(), directed=False)
    assert ncomp == 1
    assert np.all(G.weights > 0)
    assert np.all((G.fibrosis >= 0) & (G.fibrosis <= 1))
    assert G.uac.min() >= 0.0 and G.uac.max() <= 1.0
    assert G.coords.shape[1] == 3


def test_network_deterministic_in_seed():
    a = make_excitable_network(7, _small_cfg())
    b = make_excitable_network(7, _small_cfg())
    assert a.n_nodes == b.n_nodes
    assert np.allclose(a.coords, b.coords)
    assert np.allclose(a.weights, b.weights)


def test_cardiac_spectral_sfi_stack_runs_on_network_verbatim():
    """The GM4 thesis: the fragility calculus is graph-universal (literal reuse)."""
    from asb.spectral import fiedler, perron
    from asb.sfi import edge_fragility, perturbation_field, sfi_region, subspace_sfi
    from asb.config import SFIConfig
    import scipy.sparse as sp

    G = make_excitable_network(3, _small_cfg())
    W = G.adjacency()
    L = sp.csr_matrix(sp.diags(np.asarray(W.sum(1)).ravel()) - W)
    lam2, phi2 = fiedler(L)
    assert np.isfinite(lam2)
    frag = edge_fragility(phi2, G.edges)
    assert np.all(frag >= 0)
    field = perturbation_field(G, SFIConfig(), np.random.default_rng(0))
    reg = sfi_region(G, phi2, G.edges, G.region, field["expected_dw"])
    assert all(v >= -1e-12 for v in reg.values())
    sub = subspace_sfi(L, G.edges, field["expected_dw"], k_dim=2)
    assert np.all(np.isfinite(sub)) and np.all(sub >= -1e-9)
    rho, v = perron(W)
    assert rho > 0 and np.all(v >= -1e-9)


def test_fhn_label_deterministic():
    G = make_excitable_network(5, _small_cfg(), lesion_burden=0.3)
    cfg = FHNConfig(observe_after=120.0, instability_min_time=80.0)
    a = induce_fhn(G, cfg)
    b = induce_fhn(G, cfg)
    assert a.inducible == b.inducible
    assert a.reentry_origin == b.reentry_origin
    assert a.source == "fhn_network"


def test_fhn_instability_grows_with_lesion_burden():
    """Population-level: mean sustained (seizure-like) time rises with lesion burden.

    The substrate→instability link is a cohort property (single networks vary with
    topology), so we average over several seeds rather than assert per-network.
    """
    cfg = FHNConfig(observe_after=200.0)
    lo = [simulate_fhn(make_excitable_network(s, _small_cfg(), lesion_burden=0.05), cfg)["sustained_time"]
          for s in range(6)]
    hi = [simulate_fhn(make_excitable_network(s, _small_cfg(), lesion_burden=0.45), cfg)["sustained_time"]
          for s in range(6)]
    assert np.mean(hi) > np.mean(lo)
    # A well-formed inducible label carries an in-range origin.
    for s in range(6):
        G = make_excitable_network(s, _small_cfg(), lesion_burden=0.45)
        lab = induce_fhn(G, cfg)
        if lab.inducible:
            assert lab.reentry_origin is not None and 0 <= lab.reentry_origin < G.n_nodes


def test_network_cohort_has_variety():
    cohort = make_network_cohort(6, _small_cfg(), seed=0)
    assert len(cohort) == 6
    burdens = [G.meta["lesion_burden"] for G in cohort]
    assert max(burdens) > min(burdens)  # spans a gradient
