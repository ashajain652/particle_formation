import numpy as np

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.criteria import kappa_max_nodal, principal_curvatures
from reentry_model.levelset.cutcell import patches
from reentry_model.levelset.levelset import redistance

C = np.array([1.5e-3, 1.5e-3, 1.5e-3])


def _sphere_mesh(h=0.125e-3):
    return build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), h))


def test_a_banded_reset_matches_the_full_one_inside_the_band_and_never_exceeds_it_outside():
    m = _sphere_mesh(0.125e-3)
    phi0 = 3.0 * (np.linalg.norm(m.points - C, axis=1) - 0.9e-3)
    full = redistance(m, phi0)
    banded = redistance(m, phi0, band=0.3e-3)
    inside = np.abs(full) < 0.25e-3
    assert np.array_equal(banded[inside], full[inside])
    assert (np.abs(banded) <= np.abs(full) + 1e-15).all() and np.abs(banded).max() == 0.3e-3


def test_local_curvature_equals_the_whole_mesh_recovery_at_the_surface():
    from reentry_model.levelset.criteria import hessian_field
    m = _sphere_mesh(0.125e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.1e-3)
    P = patches(m, phi)
    k_local = principal_curvatures(m, phi, P.centroid)
    k_global = principal_curvatures(m, phi, P.centroid, hessian_field(m, phi))
    assert np.allclose(k_local, k_global, rtol=1e-12)
    km = kappa_max_nodal(m, phi, band=0.3e-3)
    assert (km[np.abs(phi) >= 0.3e-3] == 0).all() and (km[np.abs(phi) < 0.2e-3] > 0).all()
