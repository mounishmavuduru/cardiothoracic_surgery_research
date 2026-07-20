"""Tests for the labels layer: mock_ep surrogate and openCARP interface.

These tests depend only on ``asb.types``, ``asb.config``, numpy, and the two
label modules under test. Inputs (a small synthetic ``AtrialGraph`` / a small
``AtrialMesh``) are constructed locally so the tests are hermetic and do not
import sibling ``asb.*`` modules being written concurrently.
"""
from __future__ import annotations

import numpy as np
import pytest

from asb.config import LabelConfig
from asb.labels.mock_ep import induce, paced_activation
from asb.labels.opencarp import (
    OpenCARPUnavailableError,
    make_pacing_protocol,
    reproduce_niederer_benchmark,
    run_opencarp,
    write_carp_mesh,
)
from asb.types import AtrialGraph, AtrialMesh, InducibilityLabel


# --------------------------------------------------------------------------- #
# Local synthetic fixtures
# --------------------------------------------------------------------------- #
def _grid_graph(nx: int = 6, ny: int = 5, *, fibrotic: bool = True) -> AtrialGraph:
    """A small 2D grid conduction graph with an optional fibrotic slow region.

    A grid gives clear source-sink expansions and a heterogeneous wavelength
    field when part of it is made low-conductance / fibrotic.
    """
    n = nx * ny
    xs, ys = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    coords = np.column_stack([xs.ravel(), ys.ravel(), np.zeros(n)]).astype(float)

    def idx(ix, iy):
        return ix * ny + iy

    edges = []
    for ix in range(nx):
        for iy in range(ny):
            if ix + 1 < nx:
                edges.append((idx(ix, iy), idx(ix + 1, iy)))
            if iy + 1 < ny:
                edges.append((idx(ix, iy), idx(ix, iy + 1)))
    edges = np.array([(min(a, b), max(a, b)) for a, b in edges], dtype=np.int64)

    fibrosis = np.zeros(n)
    weights = np.ones(edges.shape[0])
    if fibrotic:
        # A patch of fibrosis in the middle columns lowers local conductance.
        patch = (xs.ravel() >= nx // 3) & (xs.ravel() <= nx // 2)
        fibrosis[patch] = 0.8
        fib_edge = 0.5 * (fibrosis[edges[:, 0]] + fibrosis[edges[:, 1]])
        weights = np.clip(1.0 - fib_edge, 0.05, 1.0)

    uac = np.column_stack(
        [xs.ravel() / max(nx - 1, 1), ys.ravel() / max(ny - 1, 1)]
    ).astype(float)
    region = np.zeros(n, dtype=np.int64)
    return AtrialGraph(
        coords=coords,
        edges=edges,
        weights=weights,
        fibrosis=fibrosis,
        uac=uac,
        region=region,
        shape_family="grid_test",
    )


def _small_mesh() -> AtrialMesh:
    """A tiny two-triangle quad mesh (all vertices used by >= 1 face)."""
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float
    )
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    fibres = np.tile([1.0, 0.0, 0.0], (4, 1))
    uac = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    fibrosis = np.zeros(4)
    region = np.array([0, 0, 1, 1], dtype=np.int64)
    return AtrialMesh(
        points=points,
        faces=faces,
        fibres=fibres,
        uac=uac,
        fibrosis=fibrosis,
        region=region,
        shape_family="mesh_test",
    )


# --------------------------------------------------------------------------- #
# mock_ep.induce / paced_activation
# --------------------------------------------------------------------------- #
def test_induce_returns_valid_label():
    G = _grid_graph()
    cfg = LabelConfig(protocol="S1S2", n_pacing_sites=3)
    label = induce(G, cfg, np.random.default_rng(0))

    assert isinstance(label, InducibilityLabel)
    assert isinstance(label.inducible, bool)
    assert label.source == "mock_ep"
    assert label.protocol == "S1S2"
    # reentry_origin is an int-or-None, consistent with inducible.
    if label.inducible:
        assert isinstance(label.reentry_origin, int)
        assert 0 <= label.reentry_origin < G.n_nodes
    else:
        assert label.reentry_origin is None


def test_induce_deterministic_for_fixed_rng():
    G = _grid_graph()
    cfg = LabelConfig(protocol="S1S2", n_pacing_sites=3)
    a = induce(G, cfg, np.random.default_rng(42))
    b = induce(G, cfg, np.random.default_rng(42))
    assert a.inducible == b.inducible
    assert a.reentry_origin == b.reentry_origin
    assert a.protocol == b.protocol
    assert a.source == b.source


def test_induce_source_is_mock_ep_for_any_protocol():
    G = _grid_graph()
    for protocol in ("S1S2", "burst"):
        label = induce(
            G, LabelConfig(protocol=protocol), np.random.default_rng(1)
        )
        assert label.source == "mock_ep"
        assert label.protocol == protocol


def test_fibrosis_tends_to_increase_inducibility_signal():
    """High-fibrosis / low-conductance atria score a higher reentry index."""
    cfg = LabelConfig(protocol="burst", n_pacing_sites=3)
    healthy = induce(_grid_graph(fibrotic=False), cfg, np.random.default_rng(7))
    fibrotic = induce(_grid_graph(fibrotic=True), cfg, np.random.default_rng(7))
    assert (
        fibrotic.meta["reentry_index"] >= healthy.meta["reentry_index"]
    )


def test_paced_activation_shape_and_values():
    G = _grid_graph()
    cfg = LabelConfig()
    times = paced_activation(G, site=0, cfg=cfg)
    assert times.shape == (G.n_nodes,)
    assert times[0] == 0.0
    # Grid is connected -> all finite; activation times are nonnegative.
    assert np.all(np.isfinite(times))
    assert np.all(times >= 0.0)


# --------------------------------------------------------------------------- #
# opencarp interface
# --------------------------------------------------------------------------- #
def test_write_carp_mesh_row_counts(tmp_path):
    mesh = _small_mesh()
    out = write_carp_mesh(mesh, str(tmp_path))

    for key in ("pts", "elem", "lon"):
        assert key in out
        assert out[key].endswith(f"mesh.{key}") or key == "lon"
    assert out["n_points"] == mesh.n_points
    assert out["n_elems"] == mesh.faces.shape[0]

    pts_lines = open(out["pts"]).read().splitlines()
    elem_lines = open(out["elem"]).read().splitlines()
    lon_lines = open(out["lon"]).read().splitlines()

    # Header line + one row per vertex / element.
    assert int(pts_lines[0]) == mesh.n_points
    assert len(pts_lines) == mesh.n_points + 1
    assert int(elem_lines[0]) == mesh.faces.shape[0]
    assert len(elem_lines) == mesh.faces.shape[0] + 1
    # .lon: header "1" (one fibre per element) + one row per element.
    assert lon_lines[0].strip() == "1"
    assert len(lon_lines) == mesh.faces.shape[0] + 1
    # Element rows carry the triangle tag.
    assert elem_lines[1].startswith("Tr ")


def test_make_pacing_protocol_fields():
    s1s2 = make_pacing_protocol(LabelConfig(protocol="S1S2", n_pacing_sites=4))
    assert s1s2["name"] == "S1S2"
    assert s1s2["n_sites"] == 4
    assert "s2_coupling_ms" in s1s2

    burst = make_pacing_protocol(LabelConfig(protocol="burst"))
    assert burst["name"] == "burst"
    assert "burst_cycle_length_ms" in burst


def test_run_opencarp_raises_when_binary_absent(tmp_path):
    mesh_dir = tmp_path / "carp"
    mesh_dir.mkdir()
    protocol = make_pacing_protocol(LabelConfig())
    with pytest.raises(OpenCARPUnavailableError) as exc:
        run_opencarp(str(mesh_dir), protocol)
    # Message carries an install pointer.
    assert "opencarp.org" in str(exc.value)


def test_run_opencarp_explicit_missing_binary_raises(tmp_path):
    mesh_dir = tmp_path / "carp"
    mesh_dir.mkdir()
    protocol = make_pacing_protocol(LabelConfig())
    with pytest.raises(OpenCARPUnavailableError):
        run_opencarp(
            str(mesh_dir), protocol, opencarp_bin="/nonexistent/openCARP"
        )


def test_reproduce_niederer_benchmark_deferred():
    with pytest.raises(OpenCARPUnavailableError):
        reproduce_niederer_benchmark()
