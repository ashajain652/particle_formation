import numpy as np
import pytest

from reentry_model.levelset.boxmesh import BoxBase, build
from reentry_model.levelset.cutcell import inside_volumes, patches
from reentry_model.levelset.levelset import redistance
from reentry_model.levelset.surfaces import molten_rigid_surface, outer_surface, skin_base
from reentry_model.shared.material import Material

from cases import box_corners, box_sdf, edge_length, rounding_volume, stage3_box, volume_between


C = np.array([1.5e-3, 1.5e-3, 1.5e-3])
T_RIGID = Material.from_drama_json("AA7075_scheil").T_rigid     # K, 895.10 (the material these fields mean)
H = 1.2e-3                 # m: the skin depth (design section 21.3, SkinSettings.depth)
SHARP = 0.1                # of the rounding volume: an edge is sharp within this (plan 03b, planning ruling 11)
BASE_BAND = 3.0e-3         # m: the distance reset's band in the skin base's test, beyond H by more than a fine cell


def _mesh(h=0.0625e-3):
    return build(BoxBase.around((0, 0, 0), (3e-3, 3e-3, 3e-3), h))


def test_the_outer_surface_stands_the_skin_thickness_outside_the_base():
    """`outer_surface`, kept for stage 1's timing script only (plan 03b Task 5 retired its use by the body)."""
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
    psi, crust = molten_rigid_surface(m, phi, T, T_RIGID)
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
    psi, crust = molten_rigid_surface(m, phi, T, T_RIGID)
    radius = np.linalg.norm(patches(m, psi).centroid - C, axis=1)
    assert crust.any()
    assert np.abs(radius - 0.6e-3).max() < 0.1e-3     # psi follows the core only
    assert (psi >= phi - 1e-15).all()


def test_fully_molten_metal_has_no_rigid_region():
    m = _mesh(0.25e-3)
    r = np.linalg.norm(m.points - C, axis=1)
    phi = redistance(m, r - 1.2e-3)
    psi, crust = molten_rigid_surface(m, phi, np.full(len(r), 950.0), T_RIGID)
    assert np.isinf(psi).all() and (psi > 0).all() and not crust.any()


def test_the_skin_base_of_a_box_is_the_box_shrunk_by_the_depth_with_sharp_edges():
    """Design section 21.4: beta, the zero contour of phi + H with phi a signed distance, is the box shrunk by H with its
    convex edges sharp, its contour exactly phi + H's (the reset keeps it; plan 03b Task 5). The stage 3 box (8 x 10 x 10 mm, 0.13 mm off the grid), phi its exact signed distance, on the
    mesh the body's refinement rule gives it (`cases.stage3_box`). The reference is the sharp shrunk box's exact signed distance on the same nodes (the controller's ruling
    P1), so that the P1 grid's own cut at every edge is on both sides; within a tenth of what rounding the edges to H
    would change (ruling 11), both in volume and in the volume between the two contours."""
    box = stage3_box()
    lo, hi = box_corners()
    beta = skin_base(box.mesh, box.phi, H, band=BASE_BAND)
    ideal = box_sdf(lo + H, hi - H)(box.mesh.points)
    V = inside_volumes(box.mesh, beta).sum()
    ref = inside_volumes(box.mesh, ideal).sum()
    analytic = np.prod(hi - lo - 2.0 * H)
    tol = SHARP * rounding_volume(H, edge_length(hi - lo - 2.0 * H))
    between = volume_between(box.mesh, beta, ideal)
    print("skin base: V {:.4f} mm3, grid reference {:.4f}, analytic {:.4f}, tolerance {:.4f}; between {:.4f} mm3".format(
        V * 1e9, ref * 1e9, analytic * 1e9, tol * 1e9, between * 1e9))
    assert abs(V - ref) <= tol
    assert between <= tol
    # the contour is phi + H's own (plan 03b Task 5): the volume's metal and the division's layer agree cell by cell
    assert np.array_equal(inside_volumes(box.mesh, beta), inside_volumes(box.mesh, box.phi + H))


def test_the_skin_base_refuses_a_band_that_does_not_reach_the_depth():
    """The skin base needs phi a distance at least H into the metal: a reset band no deeper than H is refused."""
    m = _mesh(0.25e-3)
    phi = redistance(m, np.linalg.norm(m.points - C, axis=1) - 1.4e-3)
    for band in (0.5 * H, H):
        with pytest.raises(ValueError, match="band"):
            skin_base(m, phi, H, band=band)
    beta = skin_base(m, phi, 0.4e-3, band=0.8e-3)
    assert np.abs(np.linalg.norm(patches(m, beta).centroid - C, axis=1) - 1.0e-3).max() < 0.05e-3
