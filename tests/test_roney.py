"""Tests for the real Roney loader + field-preserving coarsening."""
from __future__ import annotations

import numpy as np

from asb.substrate.roney import (
    IIR_DENSE,
    IIR_HEALTHY,
    coarsen_mesh,
    iir_to_fibrosis,
    load_roney_mesh,
    read_vtk_polydata_full,
    uac_regions,
)


def _write_tiny_vtk(path):
    """A tiny closed-ish VTK PolyData with UAC1/UAC2/IIR + fiber_endo arrays."""
    # 4 points forming 2 triangles (a small quad).
    pts = np.array([[0, 0, 0], [1000, 0, 0], [1000, 1000, 0], [0, 1000, 0]], float)
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    uac1 = np.array([0.0, 1.0, 1.0, 0.0])
    uac2 = np.array([0.0, 0.0, 1.0, 1.0])
    iir = np.array([1.0, 1.16, 1.32, 1.5])  # -> 0, 0.5, 1, 1
    lines = ["# vtk DataFile Version 3.0", "tiny", "ASCII", "DATASET POLYDATA",
             f"POINTS {len(pts)} float"]
    lines += [" ".join(f"{v:.1f}" for v in p) for p in pts]
    lines.append(f"POLYGONS {len(faces)} {len(faces)*4}")
    lines += ["3 " + " ".join(str(i) for i in f) for f in faces]
    lines.append(f"POINT_DATA {len(pts)}")
    lines += ["SCALARS UAC1 float", "LOOKUP_TABLE default"] + [f"{v}" for v in uac1]
    lines += ["SCALARS UAC2 float", "LOOKUP_TABLE default"] + [f"{v}" for v in uac2]
    lines += ["SCALARS IIR float", "LOOKUP_TABLE default"] + [f"{v}" for v in iir]
    lines.append(f"CELL_DATA {len(faces)}")
    lines += ["VECTORS fiber_endo float", "1 0 0", "0 1 0"]
    path.write_text("\n".join(lines) + "\n")


def test_iir_to_fibrosis_ramp():
    f = iir_to_fibrosis(np.array([0.5, IIR_HEALTHY, 1.16, IIR_DENSE, 1.6]))
    assert f[0] == 0.0  # below healthy
    assert f[1] == 0.0  # at healthy
    assert abs(f[2] - 0.5) < 0.02  # midpoint
    assert f[3] == 1.0  # at dense
    assert f[4] == 1.0  # above dense clipped


def test_read_and_load(tmp_path):
    p = tmp_path / "Mesh_test.vtk"
    _write_tiny_vtk(p)
    parsed = read_vtk_polydata_full(str(p))
    assert parsed["points"].shape == (4, 3)
    assert parsed["faces"].shape == (2, 3)
    assert set(parsed["point_scalars"]) == {"UAC1", "UAC2", "IIR"}
    assert "fiber_endo" in parsed["cell_vectors"]

    mesh = load_roney_mesh(str(p))
    assert mesh.n_points == 4
    # microns -> mm.
    assert np.allclose(mesh.points.max(), 1.0)
    # UAC carried through, in [0, 1].
    assert mesh.uac.min() >= 0.0 and mesh.uac.max() <= 1.0
    # Fibrosis from IIR ramp.
    assert np.allclose(mesh.fibrosis, [0.0, 0.5, 1.0, 1.0], atol=0.02)
    # Fibres unit-normalized.
    norms = np.linalg.norm(mesh.fibres, axis=1)
    assert np.allclose(norms, 1.0)


def test_uac_regions_range():
    uac = np.array([[0.0, 0.0], [0.99, 0.99], [0.5, 0.5]])
    r = uac_regions(uac, n_alpha=3, n_beta=3)
    assert r.min() >= 0 and r.max() <= 8


def test_coarsen_preserves_fields():
    # Build a moderately fine sphere-ish mesh via the synthetic generator.
    from asb.substrate.synthetic import make_base_atrium, paint_fibrosis

    mesh = paint_fibrosis(make_base_atrium(seed=1, n_nodes=2500), seed=2, burden=0.3)
    coarse = coarsen_mesh(mesh, 600)
    assert coarse.n_points < mesh.n_points
    assert coarse.n_points > 100
    # Fields stay valid.
    assert coarse.fibrosis.min() >= 0.0 and coarse.fibrosis.max() <= 1.0
    assert coarse.uac.min() >= 0.0 and coarse.uac.max() <= 1.0
    assert np.allclose(np.linalg.norm(coarse.fibres, axis=1), 1.0, atol=1e-6)
    # Coarse mesh is a single connected component (mesh_to_graph works).
    from asb.substrate.mesh import mesh_to_graph
    G = mesh_to_graph(coarse)
    assert G.n_edges > 0
    # Mean fibrosis roughly preserved (coarsening averages).
    assert abs(coarse.fibrosis.mean() - mesh.fibrosis.mean()) < 0.15
