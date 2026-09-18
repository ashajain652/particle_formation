"""mesh.py: gmsh sphere generation and cache, boundary faces, outward normals, patch angles."""
import math
import os

import numpy as np
import pytest

from reentry_model import mesh

R = 0.05


def test_coarse_sphere_geometry(coarse_sphere_mesh):
    m = coarse_sphere_mesh
    assert 2000 < m.n_nodes < 5000 and m.n_elements > 4 * m.n_nodes // 2
    assert m.volume() == pytest.approx(4.0 / 3.0 * math.pi * R ** 3, rel=3e-3)
    s = m.surface()
    assert s.area == pytest.approx(4.0 * math.pi * R ** 2, rel=2e-3)
    assert np.allclose(np.linalg.norm(s.centroids, axis=1), R, atol=3e-4)      # centroids sag by ~h^2/(6R) below the sphere
    assert np.all(np.einsum("ij,ij->i", s.normals, s.centroids) > 0.95 * np.linalg.norm(s.centroids, axis=1))   # outward, within 18 deg of radial
    assert np.allclose(np.linalg.norm(s.normals, axis=1), 1.0)
    assert np.linalg.norm(m.points[m.centre_node()]) < 1e-9
    assert m.boundary_nodes().size == np.unique(s.faces).size and m.params["radius_m"] == R


def test_angles_and_patch_lookup(coarse_sphere_mesh):
    s = coarse_sphere_mesh.surface()
    theta = s.angles_to([1.0, 0.0, 0.0])
    assert theta.min() < math.radians(5.0) and theta.max() > math.radians(175.0)
    windward = theta < math.pi / 2
    assert s.areas[windward].sum() == pytest.approx(0.5 * s.area, rel=2e-2)
    assert s.normals[s.patch_toward([1.0, 0.0, 0.0])][0] > 0.99
    assert s.normals[s.patch_toward([-1.0, 0.0, 0.0])][0] < -0.99
    x = coarse_sphere_mesh.points[:, 0]
    assert np.allclose(s.facet_mean(x), s.centroids[:, 0])


def test_boundary_faces_are_single_use_faces():
    tets = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
    faces, opposite = mesh.boundary_faces(tets)
    assert faces.shape == (6, 3) and not any(set(f) == {1, 2, 3} for f in faces)
    assert opposite.shape == (6,)


def test_generation_is_cached(tmp_path):
    p1 = mesh.generate_sphere_mesh(0.01, 2e-3, 4e-3, str(tmp_path))
    mtime = os.path.getmtime(p1)
    p2 = mesh.generate_sphere_mesh(0.01, 2e-3, 4e-3, str(tmp_path))
    assert p1 == p2 and os.path.getmtime(p2) == mtime
    assert os.path.basename(p1) == "sphere_R10.000mm_hs2.000mm_hc4.000mm.msh"
    m = mesh.load_mesh(p1)
    assert m.volume() == pytest.approx(4.0 / 3.0 * math.pi * 0.01 ** 3, rel=3e-2)        # crude 2 mm elements on a 10 mm sphere


def test_points_are_mutable_and_surface_follows(coarse_sphere_mesh):
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets)
    area0 = m.surface().area
    m.points *= 0.5
    assert m.surface().area == pytest.approx(0.25 * area0) and m.volume() == pytest.approx(0.125 * coarse_sphere_mesh.volume())
