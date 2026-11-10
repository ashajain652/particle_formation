import numpy as np
import pytest

from lsproto.boxmesh import BoxBase, build
from lsproto.criteria import kappa_max_nodal, make_target, principal_curvatures
from lsproto.cutcell import patches
from lsproto.levelset import redistance
from lsproto.surfaces import T_RIGID, molten_rigid_surface, outer_surface

C = np.array([1.5e-3, 1.5e-3, 1.5e-3])


def _mesh(h=0.0625e-3):
    return build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), h))


def test_the_outer_surface_stands_the_skin_thickness_outside_the_base():
    m = _mesh(0.125e-3)
    r = np.linalg.norm(m.points - C, axis=1)
    phi = outer_surface(m, redistance(m, r - 0.9e-3), np.full(len(r), 0.2e-3))
    radius = np.linalg.norm(patches(m, phi).centroid - C, axis=1)
    assert np.abs(radius - 1.1e-3).max() < 0.02e-3


def test_a_radial_temperature_gives_the_rigid_sphere_of_the_analytic_radius():
    m = _mesh(0.125e-3)
    r = np.linalg.norm(m.points - C, axis=1)
    phi = redistance(m, r - 1.2e-3)
    T = 850.0 + (r / 1e-3) * 60.0                     # T_rigid at r = 0.7517 mm
    psi, crust = molten_rigid_surface(m, phi, T)
    radius = np.linalg.norm(patches(m, psi).centroid - C, axis=1)
    assert np.abs(radius - (T_RIGID - 850.0) / 60.0 * 1e-3).max() < 0.02e-3
    assert not crust.any()
    assert (psi >= phi - 1e-15).all()


def test_a_cold_crust_floating_on_melt_is_not_part_of_the_rigid_region():
    m = _mesh(0.125e-3)
    r = np.linalg.norm(m.points - C, axis=1)
    phi = redistance(m, r - 1.2e-3)
    T = np.where(r < 0.6e-3, 850.0, 920.0)            # rigid core, molten shell
    T = np.where(r > 1.0e-3, 850.0, T)                # a frozen skin on top of the melt
    psi, crust = molten_rigid_surface(m, phi, T)
    radius = np.linalg.norm(patches(m, psi).centroid - C, axis=1)
    assert crust.any()
    assert np.abs(radius - 0.6e-3).max() < 0.1e-3     # psi follows the core only
    assert (psi >= phi - 1e-15).all()


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


def test_fully_molten_metal_has_no_rigid_region():
    m = _mesh(0.25e-3)
    r = np.linalg.norm(m.points - C, axis=1)
    phi = redistance(m, r - 1.2e-3)
    psi, crust = molten_rigid_surface(m, phi, np.full(len(r), 950.0))
    assert np.isinf(psi).all() and (psi > 0).all() and not crust.any()


def test_local_curvature_equals_the_whole_mesh_recovery_at_the_surface():
    from lsproto.criteria import hessian_field
    m = _mesh(0.125e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.1e-3)
    P = patches(m, phi)
    k_local = principal_curvatures(m, phi, P.centroid)
    k_global = principal_curvatures(m, phi, P.centroid, hessian_field(m, phi))
    assert np.allclose(k_local, k_global, rtol=1e-12)
    km = kappa_max_nodal(m, phi, band=0.3e-3)
    assert (km[np.abs(phi) >= 0.3e-3] == 0).all() and (km[np.abs(phi) < 0.2e-3] > 0).all()
