import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.criteria import kappa_max_nodal, make_target, principal_curvatures
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import redistance


C = np.array([1.5e-3, 1.5e-3, 1.5e-3])


def _mesh(h=0.0625e-3):
    return build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), h))


def test_sphere_curvature_is_recovered():
    m = _mesh()
    r = np.linalg.norm(m.points - C, axis=1)
    phi = redistance(m, r - 1.1e-3)
    P = patches(m, phi)
    k1, k2 = principal_curvatures(m, phi, P.centroid)
    assert np.median(k1) == pytest.approx(1 / 1.1e-3, rel=0.05)
    assert np.median(k2) == pytest.approx(1 / 1.1e-3, rel=0.05)


def test_a_cylinder_has_one_curved_and_one_flat_direction():
    m = _mesh()
    phi = redistance(m, np.linalg.norm(m.points[:, :2] - C[:2], axis=1) - 0.8e-3)
    P = patches(m, phi)
    mid = np.abs(P.centroid[:, 2] - C[2]) < 0.8e-3
    k1, k2 = principal_curvatures(m, phi, P.centroid[mid])
    assert np.median(k1) == pytest.approx(1 / 0.8e-3, rel=0.05)
    assert np.abs(np.median(k2)) < 0.05 / 0.8e-3


def test_the_fine_region_reaches_the_margin_below_the_floor():
    m = _mesh(0.25e-3)
    r = np.linalg.norm(m.points - C, axis=1)
    phi = redistance(m, r - 1.3e-3)
    psi = redistance(m, r - 1.0e-3)
    k = np.zeros(len(r))
    tgt = make_target(m, phi, psi, k, k, fine=0.125e-3, margin=0.3e-3, lead=0.1e-3, floor=0.03e-3)
    new = build(m.base, tgt)
    rc = np.linalg.norm(new.centroids - C, axis=1)
    band = (rc > 0.65e-3) & (rc < 1.25e-3)
    assert (new.size[band] <= 0.125e-3 + 1e-15).all()
    assert (new.size[rc < 0.1e-3] == 0.25e-3).all()   # psi < -0.9 mm: beyond margin + lead + reach


def test_curvature_refines_a_tight_rim_down_to_the_floor():
    m = _mesh(0.125e-3)
    phi = redistance(m, np.linalg.norm(m.points[:, :2] - C[:2], axis=1) - 0.2e-3)   # 0.2 mm wire
    km = kappa_max_nodal(m, phi)
    tgt = make_target(m, phi, phi, km, km, fine=np.inf, margin=0.0, lead=0.0, floor=0.0625e-3)
    new = build(m.base, tgt)
    near = np.abs(np.linalg.norm(new.centroids[:, :2] - C[:2], axis=1) - 0.2e-3) < 0.05e-3
    assert (new.size[near] <= 0.0625e-3 + 1e-15).all()


def test_a_fully_molten_sphere_is_fine_everywhere_inside_the_metal_even_on_faces_of_the_old_mesh():
    base = BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), 0.5e-3)
    r0 = lambda x: np.linalg.norm(x - C, axis=1) - 1.2e-3
    old = build(base, lambda x, s: np.where(r0(x) < 0.94 * s, 0.125e-3, np.inf))     # refined inside: faces at 0.5 mm
    phi = redistance(old, r0(old.points))
    psi = np.full(len(old.points), np.inf)                                             # nothing rigid
    k = kappa_max_nodal(old, phi, band=1e-3)
    tgt = make_target(old, phi, psi, k, k, fine=0.125e-3, margin=0.3e-3, lead=0.1e-3, floor=0.03e-3)
    base_centroids = build(base).centroids
    inside = r0(base_centroids) < -0.1e-3
    got = tgt(base_centroids, base.h0 * np.ones(len(base_centroids)))
    assert inside.sum() > 20 and (got[inside] <= 0.125e-3).all() and not np.isnan(got).any()
    new = build(base, tgt)
    rc = np.linalg.norm(new.centroids - C, axis=1)
    assert (new.size[rc < 1.0e-3] <= 0.125e-3 + 1e-15).all()


def test_the_curvature_of_a_field_with_no_finite_value_is_zero():
    m = _mesh(0.25e-3)
    out = kappa_max_nodal(m, np.full(len(m.points), np.inf))
    assert (out == 0).all() and out.shape == (len(m.points),)
