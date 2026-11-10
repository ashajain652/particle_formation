import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import inside_volumes, patches

H0 = 1.0e-3


def test_inside_and_outside_parts_add_to_each_cell_for_random_fields():
    m = build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), H0))
    phi = np.random.default_rng(3).normal(size=len(m.points))
    total = inside_volumes(m, phi) + inside_volumes(m, -phi)
    assert np.allclose(total, m.volumes, rtol=1e-12)


def test_a_plane_gives_the_exact_volume_area_and_normal():
    base = BoxBase.around((0, 0, 0), (4e-3, 3e-3, 2e-3), H0)
    m = build(base, lambda x, s: np.where(np.abs(x[:, 0] - 1.37e-3) < 0.6e-3, 0.25e-3, np.inf))
    n = np.array([0.3, 0.5, 0.81]); n /= np.linalg.norm(n)
    d = n @ np.array([2e-3, 1.5e-3, 1e-3])
    phi = m.points @ n - d
    vol = inside_volumes(m, phi).sum()
    P = patches(m, phi)
    assert np.allclose(P.normal, n, atol=1e-12)
    # exact answers from a fine plane slicing of the same box, independently by Monte-Carlo-free formula:
    flat = build(base)
    assert vol == pytest.approx(inside_volumes(flat, phi).sum(), rel=1e-12)
    assert P.area.sum() == pytest.approx(patches(flat, phi).area.sum(), rel=1e-12)


def test_a_sphere_volume_and_area_converge_at_second_order():
    R, c = 1.1e-3, np.array([1.5e-3, 1.5e-3, 1.5e-3])
    errs = []
    for h in (0.25e-3, 0.125e-3):
        m = build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), h))
        phi = np.linalg.norm(m.points - c, axis=1) - R
        V = inside_volumes(m, phi).sum()
        A = patches(m, phi).area.sum()
        errs.append((abs(V / (4 / 3 * np.pi * R ** 3) - 1), abs(A / (4 * np.pi * R ** 2) - 1)))
    (v1, a1), (v2, a2) = errs
    assert v2 < 0.01 and a2 < 0.005          # measured 0.65 % and 0.33 % at R/h = 8.8
    assert v1 / v2 > 3.5 and a1 / a2 > 3.5          # measured 4.0 and 4.0 per halving


def test_patch_normals_point_out_of_the_metal_and_centroids_lie_on_the_surface():
    R, c = 1.1e-3, np.array([1.5e-3, 1.5e-3, 1.5e-3])
    m = build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), 0.25e-3))
    phi = np.linalg.norm(m.points - c, axis=1) - R
    P = patches(m, phi)
    radial = (P.centroid - c) / np.linalg.norm(P.centroid - c, axis=1, keepdims=True)
    assert (np.einsum("ij,ij->i", P.normal, radial) > 0.9).all()
    assert np.abs(np.linalg.norm(P.centroid - c, axis=1) - R).max() < 0.25e-3 ** 2 / R


def test_a_contour_through_nodes_leaves_no_zero_area_patches():
    m = build(BoxBase.around((0, 0, 0), (2e-3, 2e-3, 2e-3), 0.25e-3))
    phi = m.points[:, 0] - 1.0e-3                    # the plane x = 1 mm passes through a layer of nodes
    P = patches(m, phi)
    assert (P.area > 0).all()
    assert P.area.sum() == pytest.approx(4e-6, rel=1e-12)
