"""Tests for the substrate layer: mesh<->graph, synthetic anatomy, loaders.

These tests depend only on ``asb.types``, ``asb.config``, numpy, and the
substrate modules under test (plus ``asb.graph`` indirectly via
``mesh_to_graph``). Inputs are constructed locally.
"""
from __future__ import annotations

import numpy as np
import pytest

from asb.config import CohortConfig
from asb.substrate.loaders import (
    DataUnavailableError,
    load_fibre_atlas,
    load_rodero_cohort,
    load_roney_meshes,
)
from asb.substrate.mesh import geodesic_edge_length, mesh_edges, mesh_to_graph
from asb.substrate.synthetic import (
    make_base_atrium,
    make_cohort,
    paint_fibrosis,
    shape_variant,
)


# --------------------------------------------------------------------------- #
# mesh.py
# --------------------------------------------------------------------------- #
def _tetra_mesh():
    """A tiny closed tetrahedron surface for local mesh tests."""
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float
    )
    faces = np.array(
        [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=np.int64
    )
    return points, faces


def test_mesh_edges_unique_and_sorted():
    _points, faces = _tetra_mesh()
    edges = mesh_edges(faces)
    # Tetrahedron has 6 undirected edges.
    assert edges.shape == (6, 2)
    assert np.all(edges[:, 0] < edges[:, 1])
    # Rows are unique.
    assert len(np.unique(edges, axis=0)) == edges.shape[0]


def test_geodesic_edge_length():
    points, faces = _tetra_mesh()
    edges = mesh_edges(faces)
    lengths = geodesic_edge_length(points, edges)
    assert lengths.shape == (edges.shape[0],)
    assert np.all(lengths > 0)


def test_mesh_to_graph_returns_valid_graph():
    mesh = make_base_atrium(seed=1, n_nodes=200)
    G = mesh_to_graph(mesh, along=1.0, cross=0.3)
    # i < j on every edge.
    assert np.all(G.edges[:, 0] < G.edges[:, 1])
    # Positive weights.
    assert np.all(G.weights > 0)
    # Edge indices valid.
    assert G.edges.min() >= 0
    assert G.edges.max() < G.n_nodes
    # Carried-through fields.
    assert G.shape_family == mesh.shape_family
    assert G.coords.shape == mesh.points.shape


# --------------------------------------------------------------------------- #
# synthetic.py
# --------------------------------------------------------------------------- #
def _assert_valid_surface(mesh):
    n = mesh.points.shape[0]
    faces = mesh.faces
    # Faces index in range.
    assert faces.min() >= 0
    assert faces.max() < n
    # Every vertex used by >= 1 face.
    used = np.unique(faces)
    assert used.size == n
    assert set(used.tolist()) == set(range(n))


def test_base_atrium_is_valid_surface():
    mesh = make_base_atrium(seed=3, n_nodes=200)
    _assert_valid_surface(mesh)
    # Unit fibre vectors.
    norms = np.linalg.norm(mesh.fibres, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-6)
    # UAC in [0, 1].
    assert mesh.uac.min() >= 0.0
    assert mesh.uac.max() <= 1.0


def test_make_base_atrium_deterministic():
    a = make_base_atrium(seed=5, n_nodes=200)
    b = make_base_atrium(seed=5, n_nodes=200)
    assert np.array_equal(a.points, b.points)
    assert np.array_equal(a.faces, b.faces)


def test_paint_fibrosis_in_range():
    mesh = make_base_atrium(seed=2, n_nodes=200)
    painted = paint_fibrosis(mesh, seed=9, burden=0.2, n_patches=5)
    assert painted.fibrosis.min() >= 0.0
    assert painted.fibrosis.max() <= 1.0
    # Original untouched (returns a copy).
    assert np.all(mesh.fibrosis == 0.0)
    # Deterministic.
    again = paint_fibrosis(mesh, seed=9, burden=0.2, n_patches=5)
    assert np.array_equal(painted.fibrosis, again.fibrosis)


def test_shape_variant_preserves_family_and_validity():
    base = make_base_atrium(seed=4, n_nodes=200, family="fam_X")
    var = shape_variant(base, seed=11)
    assert var.shape_family == base.shape_family == "fam_X"
    _assert_valid_surface(var)
    # Same topology, perturbed geometry.
    assert np.array_equal(var.faces, base.faces)
    assert not np.array_equal(var.points, base.points)


def test_make_cohort_shapes_and_families():
    cfg = CohortConfig(n_base=2, n_variants=3, n_nodes=200, seed=0)
    cohort = make_cohort(cfg)
    assert len(cohort) == cfg.n_base * cfg.n_variants
    families = {m.shape_family for m in cohort}
    assert families == {"family_0", "family_1"}
    for m in cohort:
        _assert_valid_surface(m)
        assert m.fibrosis.min() >= 0.0
        assert m.fibrosis.max() <= 1.0


# --------------------------------------------------------------------------- #
# loaders.py
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "loader,record",
    [
        (load_rodero_cohort, "4506930"),
        (load_roney_meshes, "5801337"),
        (load_fibre_atlas, "3764917"),
    ],
)
def test_loaders_raise_on_missing_path(loader, record, tmp_path):
    missing = str(tmp_path / "does_not_exist")
    with pytest.raises(DataUnavailableError) as exc:
        loader(missing)
    assert record in str(exc.value)


def test_loaders_raise_on_empty_dir(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(DataUnavailableError):
        load_rodero_cohort(str(empty))
