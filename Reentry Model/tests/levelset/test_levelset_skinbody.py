"""Tests of the level set's skinned body (design sections 10.1, 12.2 and 20; level-set stage 3 Task 9): one macro step
coupling the skins on the gas-surface markers to the cut-cell 3D solve -- the books of mass and enthalpy, a cold body left
alone, the deferred load reaching the volume, the markers on phi after the surface moves, and a skin deeper than its
radius of curvature held back on its own.

The cube is the AA7075 Scheil alloy's (the shared Material). The test meshes are coarser than the design's defaults
(fine 1.4 mm instead of 0.7 mm, a 0.7 mm curvature floor instead of 0.25 mm, a 2 mm fine margin instead of 4 mm) so that a
macro step takes tens of seconds; the skins are the development preset's (20 um cells)."""
import time

import numpy as np
import pytest

from reentry_model.levelset import skinbody
from reentry_model.levelset.boxmesh import BoxBase
from reentry_model.levelset.skin import SKIN_PRESETS, Skins, SkinSettings
from reentry_model.levelset.transfer import interpolate
from reentry_model.shared.material import Material

HALF = 5e-3            # m: the cube's half side (a 10 mm cube)
OFF = 0.13e-3          # m: the cube's centre off the grid, so that its faces cut the cells
SETTINGS = SkinSettings(cell=SKIN_PRESETS["dev"])
MESH = dict(fine_margin=2e-3, fine_size=1.4e-3, cell_floor=0.7e-3)


def aa7075():
    return Material.from_drama_json("AA7075_scheil")


def cube(x):
    return np.max(np.abs(x - OFF) - HALF, axis=1)


def cube_body(**kw):
    base = BoxBase.around((-9e-3,) * 3, (9e-3,) * 3, 2.8e-3)
    args = dict(MESH, emissivity=0.0)
    args.update(kw)
    return skinbody.LevelSetSkinBody(aa7075(), base, cube, 850.0, SETTINGS, **args)


def one_face(q):
    """q [W/m2] on the markers facing +x (the cube's +x face, its rounded edges by their normals' share), none elsewhere."""
    return lambda markers: q * np.clip((markers.normal[:, 0] - 0.5) / 0.5, 0.0, 1.0)


@pytest.fixture(scope="module")
def heated():
    body = cube_body()
    before = body.books()
    started = time.perf_counter()
    book = body.advance(0.5, q_prescribed=one_face(2e6))
    return body, before, book, time.perf_counter() - started


def test_a_macro_step_closes_mass_and_enthalpy(heated):
    """The 10 mm cube at 850 K heated 2 MW/m2 on its +x face for one 0.5 s macro step: the mass closes to 1e-12 and the
    enthalpy to 1e-8 of the stored enthalpy plus the heat in, the volume and the skins both gained heat, the deferred
    interface load is a small part of the heat in, and the heat in is what the face took."""
    body, before, book, seconds = heated
    print("heated cube: {:.1f} s for the step, {} parts, at most {} passes; heat in {:.4f} J, radiated {:.3e} J; "
          "stored volume {:.4f} J + skins {:.4f} J; pending {:.3e} J; residuals mass {:.2e}, enthalpy {:.2e}; "
          "markers {}; top cell at most {:.2f} K".format(
              seconds, book.parts, book.passes, book.enthalpy_in, book.enthalpy_radiated, book.enthalpy_volume,
              book.enthalpy_skins, book.pending_J, book.mass_residual, book.enthalpy_residual, book.n_markers,
              np.nanmax(body.skins.T()[:, 0])))
    assert abs(book.mass_residual) <= 1e-12
    assert abs(book.enthalpy_residual) <= 1e-8
    area = (body.markers.area * np.clip((body.markers.normal[:, 0] - 0.5) / 0.5, 0.0, 1.0)).sum()
    assert book.enthalpy_in == pytest.approx(2e6 * area * 0.5, rel=0.05)
    assert book.enthalpy_volume > before["enthalpy_volume"] and book.enthalpy_skins > before["enthalpy_skins"]
    assert abs(book.pending_J) < 0.01 * book.enthalpy_in
    for key in ("n_markers", "marker_splits", "marker_merges", "markers_trapped", "interface_mismatch_J",
                "volume_correction_max", "n_pieces", "frontal_area_m2"):
        assert key in book.row()
    assert book.n_pieces == 1 and book.frontal_area_m2 == pytest.approx((2 * HALF + 1.6e-3) ** 2, rel=0.1)
    for _, key in skinbody.SKIN_SUMMARY_COLUMNS:
        assert key in book.row()


def test_markers_follow_the_recession():
    """The cube's +x face at 6 MW/m2 for 0.5 s with the liquid removed as it forms (remove="feed"): the face recedes
    (metal leaves, the skins draw from the volume and the skin base beta moves inward), the books close, and after the
    step every marker lies within 1e-9 m of phi = 0, the outer surface rebuilt from beta and the skins."""
    body = cube_body(remove="feed")
    x0 = body.markers.x[body.markers.normal[:, 0] > 0.9, 0].mean()
    book = body.advance(0.5, q_prescribed=one_face(6e6))
    mk = body.markers
    face = mk.normal[:, 0] > 0.9
    off = np.abs(interpolate(body.mesh, body.phi, mk.x))
    print("receding face: removed {:.3e} kg, drawn {:.3e} kg; face moved {:.3f} mm; markers off phi at most {:.2e} m; "
          "residuals mass {:.2e}, enthalpy {:.2e}".format(book.mass_removed, book.mass_drawn,
                                                          (x0 - mk.x[face, 0].mean()) * 1e3, off.max(),
                                                          book.mass_residual, book.enthalpy_residual))
    assert book.mass_removed > 0.0 and book.mass_drawn > 0.0
    assert mk.x[face, 0].mean() < x0
    assert off.max() <= 1e-9
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8


def test_a_cold_body_is_unchanged_by_the_skins():
    """With no heat, no flow and no radiation, a macro step leaves the volume's temperatures, the skin base beta and the
    skins as they were."""
    body = cube_body()
    T0, beta0 = body.solver.T.copy(), body.beta.copy()
    sk = body.skins
    m0, H0, c0 = sk.m.copy(), sk.H.copy(), sk.count.copy()
    book = body.advance(0.5)
    active = body.solver.active
    assert np.abs(body.solver.T[active] - T0[active]).max() <= 1e-9
    assert np.array_equal(body.beta, beta0)
    assert np.array_equal(sk.count, c0)
    np.testing.assert_allclose(sk.m, m0, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(sk.H, H0, rtol=1e-12, atol=0.0)
    assert book.enthalpy_in == 0.0 and abs(book.pending_J) <= 1e-9 * abs(book.enthalpy_volume)


def test_pending_load_reaches_the_volume():
    """A deferred load of 5 J at a point under the +x face is delivered to the volume, at the patch of beta nearest the
    point, in the next part (a 0.05 s step): the body's stored enthalpy rises by the 5 J less the load still deferred
    (to 1e-9), the volume keeps most of it and the skins above the warmed base take the rest through their base (they
    have no other source here), the load still deferred is small, and the volume's hottest node lies near the point."""
    body = cube_body()
    point = np.array([HALF + OFF - 1.0e-3, OFF, OFF])
    body.add_pending(point, 5.0)
    before = body.books()
    book = body.advance(0.05)
    gained = book.enthalpy_volume + book.enthalpy_skins - before["enthalpy_volume"] - before["enthalpy_skins"]
    print("pending load: volume gained {:.6f} J, skins {:.6f} J, still deferred {:.3e} J".format(
        book.enthalpy_volume - before["enthalpy_volume"], book.enthalpy_skins - before["enthalpy_skins"], book.pending_J))
    assert gained + book.pending_J == pytest.approx(5.0, rel=1e-9)
    assert book.enthalpy_volume - before["enthalpy_volume"] > 2.5
    assert book.enthalpy_skins - before["enthalpy_skins"] >= 0.0
    assert abs(book.pending_J) < 0.05
    T = body.solver.T
    hot = np.argmax(np.where(body.solver.active, T, -np.inf))
    assert np.linalg.norm(body.mesh.points[hot] - point) < 3e-3


def test_a_skin_deeper_than_its_radius_of_curvature_is_held_back_alone():
    """Skins.dz raises for every column if one column is deeper than its radius of curvature (a skin of 0.8 mm on the
    cube's rounded edge, about 0.8 mm in radius). The body caps each column's principal curvatures at
    KAPPA_DEPTH_SHARE over its flat depth (a sub-grid edge, design 20.6's known limitation 1): only the too-deep column is
    changed, and its thickness is then defined."""
    mat = aa7075()
    kappa = np.array([[40.0, 40.0], [1250.0, 0.0], [-300.0, 20.0]])
    sk = Skins(mat, np.full(3, 1e-6), kappa, SETTINGS)
    sk.fill(mat.rho * 1e-6 * 0.8e-3, 850.0)
    with pytest.raises(ValueError, match="deeper than its radius of curvature"):
        sk.dz()
    capped = skinbody.cap_curvature(sk, kappa)
    assert np.array_equal(capped[[0, 2]], kappa[[0, 2]])
    assert capped[1, 0] == pytest.approx(skinbody.KAPPA_DEPTH_SHARE / 0.8e-3, rel=0.05) and capped[1, 1] == 0.0
    sk.kappa = capped
    assert np.all(np.isfinite(sk.thickness()))


def test_a_step_at_a_gas_state_with_the_flow_and_the_spray_on_closes_its_books():
    """The cube at a gas state (77.5 km, 7.5 km/s, the freestream density raised to 1e-4 kg/m3 so that it melts within
    the step) with the shared heating and surface flow (`loads.GasLoads`), the runoff driven by the surface flow and the
    deceleration, radiation, and Girin's plan (the shared SprayModel): one 0.1 s macro step closes mass and enthalpy,
    the plan's books hold exactly what the step booked as sprayed, the history row carries the droplet columns, and
    what a trial pass is undone to restores the cells, the flow memory and the plan's books. (At this state the top
    does not reach the liquidus within 0.1 s, so nothing sprays: the release itself is tested in Task 8's tests.)"""
    from test_levelset_loads import _state
    from reentry_model.levelset import loads
    from reentry_model.shared import heating, spray, surface_flow
    from reentry_model.shared.material import LiquidProperties
    mat = aa7075()
    gl = loads.GasLoads(heating.PhysicsHeating(), surface_flow.SurfaceFlow(), mat.liquid.rho)
    model = spray.SprayModel(LiquidProperties(mat.liquid.rho, mat.liquid.mu, mat.liquid.sigma), spray_gates="2026-11-10")
    body = cube_body(loads=gl, spray_model=model, T_ambient=200.0, emissivity=0.2)
    started = time.perf_counter()
    book = body.advance(0.1, state=_state(rho=1e-4))
    print("gas state: {:.1f} s, {} parts, {} passes; heat in {:.3f} J, radiated {:.3e} J, sprayed {:.3e} kg; residuals "
          "mass {:.2e}, enthalpy {:.2e}; momentum to the wall {}".format(
              time.perf_counter() - started, book.parts, book.passes, book.enthalpy_in, book.enthalpy_radiated,
              book.mass_sprayed, book.mass_residual, book.enthalpy_residual, book.momentum_to_wall))
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8
    assert book.enthalpy_in > 0.0 and book.enthalpy_radiated > 0.0
    plan = body.plan
    assert plan is not None
    assert plan.mass.sum() + plan.mass_slurry.sum() == pytest.approx(book.mass_sprayed, rel=1e-12, abs=1e-30)
    row = book.row()
    assert row["altitude_km"] == pytest.approx(77.5) and np.isfinite(row["flowing_depth_mean_mm"])
    # the state a trial pass is undone to: the cells, the flow memory and the plan's books (no release spray here, so
    # the restore is checked directly)
    ctx = {"plan": plan}
    saved = body._save(ctx)
    m, memory, mass = body.skins.m.copy(), body.skins.flow_memory.copy(), plan.mass.copy()
    plan.mass += 1.0
    body.skins.flow_memory = memory + 1.0
    body.skins.m[0, 0] *= 2.0
    body._restore(saved, ctx)
    assert np.array_equal(plan.mass, mass) and np.array_equal(body.skins.flow_memory, memory)
    assert np.array_equal(body.skins.m, m)
