import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import redistance
from reentry_model.levelset.surfaces import T_RIGID, molten_rigid_surface, outer_surface


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


def test_fully_molten_metal_has_no_rigid_region():
    m = _mesh(0.25e-3)
    r = np.linalg.norm(m.points - C, axis=1)
    phi = redistance(m, r - 1.2e-3)
    psi, crust = molten_rigid_surface(m, phi, np.full(len(r), 950.0))
    assert np.isinf(psi).all() and (psi > 0).all() and not crust.any()
