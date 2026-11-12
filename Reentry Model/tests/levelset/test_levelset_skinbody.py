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

from cases import CAPPED_SHARE, DESIGN_FINE, DESIGN_SPACING, OWNERSHIP_FACTOR, gas_state, nearest_ratio, seeded

HALF = 5e-3            # m: the cube's half side (a 10 mm cube)
OFF = 0.13e-3          # m: the cube's centre off the grid, so that its faces cut the cells
SETTINGS = SkinSettings(cell=SKIN_PRESETS["dev"])
MESH = dict(fine_margin=2e-3, fine_size=1.4e-3, cell_floor=0.7e-3)
OWED_SHARE = 0.05      # the heat still owed to the volume at a step's end, summed in absolute value per point, and the last
                       # part's interface mismatch, each at most this share of the step's heat in (fix round 1, Important 2).
                       # Measured once each: 4e-4 on the heated face; 1.2e-2 on the receding face; 2.9e-2 on the spraying
                       # face's 0.05 s step. What is owed is almost all the base moves' deferral -- the draw values the
                       # metal it takes on the profile capped at T_rigid, the volume removes it at its own, hotter,
                       # temperatures -- about 4 % of the heat in over the receding step, most of it delivered within it
                       # (the mismatch itself is 1e-4 of the heat in or less)
OWED_SHARE_HEATED = 0.01    # the same on the heated face, whose step owes little (measured 5.3e-4; fix round 5)
PENDING_FLUX_SHARE = 10.0   # the pending heat's flux on a patch at most this times the largest surface flux (fix round 3)


def owed(body, book):
    """(the heat still owed to the volume, summed in absolute value [J], the last part's sum |mismatch| [J]) as shares of
    the step's heat in, printed with the step's deferred residuals by source."""
    pend, mism = np.abs(body.pending_J).sum() / book.enthalpy_in, book.interface_mismatch_J / book.enthalpy_in
    lo, hi = skinbody.GROUP_AREA_RANGE
    # a group formed again stops at the first patch that takes it past the base area, so it may overshoot the range's
    # top by that one patch (accepted, fix round 5: 2.17 on the full-size cube's second step): less its largest patch,
    # every group is within the range
    assert lo <= book.group_area_ratio_min and book.group_area_ratio_trimmed_max <= hi
    print("  pending flux at most {:.3e} W/m2 against the surface's {:.3e} W/m2".format(book.pending_flux_max,
                                                                                    book.surface_flux_max))
    assert book.pending_flux_max <= PENDING_FLUX_SHARE * book.surface_flux_max
    print("  owed |pending| {:.2e} and mismatch {:.2e} of the heat in; deferred Newton {:.2e} J, moves {:.2e} J, "
          "remesh {:.2e} J; groups {:.2f} to {:.2f} of the base areas".format(
              pend, mism, book.residual_newton_J, book.residual_move_J, book.residual_remesh_J,
              book.group_area_ratio_min, book.group_area_ratio_max))
    return pend, mism


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
    assert max(owed(body, book)) <= OWED_SHARE_HEATED
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
    assert max(owed(body, book)) <= OWED_SHARE
    # the moves' deferral (booked as pending at the bases) is mostly delivered to the volume within the step: what is
    # still owed at its end is under half of what the moves deferred over it, and the books close with it
    assert np.abs(body.pending_J).sum() < 0.5 * book.residual_move_J


COLD_FACTOR = 2.0         # a skin's mass and enthalpy over a cold step change by at most this times its own area change
                          # plus its stretch's departure from one (fix round 1 of Task 9b; measured: 1.08 times at most)
COLD_STRETCH = 1e-5       # the stretch of a cold step within this of one (measured 2.4e-6: phi itself moves by 1.9e-8 m
                          # on a cold step, the step end not settling, Task 11)


def test_a_cold_body_is_unchanged_by_the_skins():
    """With no heat, no flow and no radiation, a macro step leaves the volume's temperatures and the skin base beta as
    they were, and the skins at their temperatures with their books. Since Task 9b each column's metal follows the area
    the step end's tessellation re-partitions, and on a cold step that is round-off of the snap: phi moves by up to
    1.9e-8 m, the markers by 13 nm, their areas by up to 3.0e-6 (measured). Each column's mass and enthalpy change by at
    most COLD_FACTOR times its own area change plus its stretch's departure from one; the stretch, as the step's book
    records it, is one within COLD_STRETCH; the book's share of area passed between markers is the transfer's; the skins'
    totals hold to 1e-14. (Fix round 1: a column's change bounded by its own area change alone failed where its
    neighbours' areas moved and its own hardly did, by up to 40 times.)"""
    body = cube_body()
    T0, beta0 = body.solver.T.copy(), body.beta.copy()
    sk, mk = body.skins, body.markers
    M0, E0 = (t[mk.skin_index] for t in sk.totals()[:2])
    A0, n0, T_skin0 = sk.area[mk.skin_index].copy(), mk.n, np.nanmax(np.abs(sk.T() - 850.0))
    book = body.advance(0.5)
    active = body.solver.active
    assert np.abs(body.solver.T[active] - T0[active]).max() <= 1e-9
    assert np.array_equal(body.beta, beta0)
    assert mk.n == n0 and len(body.last_settle) == 1
    rep = body.last_settle[0]
    M1, E1 = (t[mk.skin_index] for t in sk.totals()[:2])
    bound = COLD_FACTOR * (np.abs(mk.area / A0 - 1.0) + np.abs(rep.stretch - 1.0))
    T = rep.transfer
    print("cold step: columns' mass changed by at most {:.1e} ({:.2f} of the bound), areas by {:.1e}; stretch {:.7f} to "
          "{:.7f}; {:.1e} of the area passed between markers".format(
              np.abs(M1 / M0 - 1).max(), (np.abs(M1 / M0 - 1) / bound).max(), np.abs(mk.area / A0 - 1).max(),
              book.stretch_min, book.stretch_max, book.area_repartitioned_share))
    assert np.all(np.abs(M1 / M0 - 1) <= bound) and np.all(np.abs(E1 / E0 - 1) <= bound)
    assert abs(M1.sum() / M0.sum() - 1) <= BOOKS_RTOL and abs(E1.sum() / E0.sum() - 1) <= BOOKS_RTOL
    assert book.stretch_min == rep.stretch.min() and book.stretch_max == rep.stretch.max()
    assert 1.0 - COLD_STRETCH <= book.stretch_min <= book.stretch_max <= 1.0 + COLD_STRETCH
    assert book.area_repartitioned_share == pytest.approx((T.sum() - T.diagonal().sum()) / mk.area.sum(), rel=1e-12)
    assert np.nanmax(np.abs(sk.T() - 850.0)) <= max(T_skin0, 1e-9)
    assert book.enthalpy_in == 0.0 and abs(book.pending_J) <= 1e-9 * abs(book.enthalpy_volume)
    for key in ("area_repartitioned_share", "stretch_min", "stretch_max"):
        assert key in book.row()


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

    # the load is spread over the group of the marker whose base is nearest the point, never on one patch alone
    body.add_pending(point, 5.0)
    S, _, base, _ = body._interface()
    spread = body._pending_per_patch(S, base)
    assert spread.sum() == pytest.approx(body.pending_J.sum(), rel=1e-12)
    assert (spread != 0.0).sum() > 1


def test_a_skin_deeper_than_its_radius_of_curvature_is_held_back_alone():
    """Skins.dz raises for every column if one column is deeper than its radius of curvature (a skin of 0.8 mm on the
    cube's rounded edge, about 0.8 mm in radius). The body caps each column's principal curvatures at
    KAPPA_DEPTH_SHARE over its flat depth either way (a sub-grid edge, design 20.6's known limitation 1): only the
    columns over it change (the too-deep one, and a sharply concave one), and the thickness is then defined."""
    mat = aa7075()
    kappa = np.array([[40.0, 40.0], [1250.0, 0.0], [-100.0, 20.0], [-1500.0, 0.0]])
    sk = Skins(mat, np.full(4, 1e-6), kappa, SETTINGS)
    sk.fill(mat.rho * 1e-6 * 0.8e-3, 850.0)
    with pytest.raises(ValueError, match="deeper than its radius of curvature"):
        sk.dz()
    capped = skinbody.cap_curvature(sk, kappa)
    assert np.array_equal(capped[[0, 2]], kappa[[0, 2]])
    assert capped[1, 0] == pytest.approx(skinbody.KAPPA_DEPTH_SHARE / 0.8e-3, rel=0.05) and capped[1, 1] == 0.0
    # a concave column is held too (its widening base would raise its minimum mass faster than its mass)
    assert capped[3, 0] == pytest.approx(-skinbody.KAPPA_DEPTH_SHARE / 0.8e-3, rel=0.05) and capped[3, 1] == 0.0
    sk.kappa = capped
    assert np.all(np.isfinite(sk.thickness()))


def test_a_step_at_a_gas_state_with_the_flow_and_the_spray_on_closes_its_books():
    """The cube at a gas state (77.5 km, 7.5 km/s, the freestream density raised to 1e-4 kg/m3 so that it melts within
    the step) with the shared heating and surface flow (`loads.GasLoads`), the runoff driven by the surface flow and the
    deceleration, radiation, and Girin's plan (the shared SprayModel): one 0.1 s macro step closes mass and enthalpy,
    the plan's books hold exactly what the step booked as sprayed, the history row carries the droplet columns, and
    what a trial pass is undone to restores the cells, the flow memory and the plan's books. (At this state the top
    does not reach the liquidus within 0.1 s, so nothing sprays: the release itself is tested in Task 8's tests.)"""
    from reentry_model.levelset import loads
    from reentry_model.shared import heating, spray, surface_flow
    from reentry_model.shared.material import LiquidProperties
    mat = aa7075()
    gl = loads.GasLoads(heating.PhysicsHeating(), surface_flow.SurfaceFlow(), mat.liquid.rho)
    model = spray.SprayModel(LiquidProperties(mat.liquid.rho, mat.liquid.mu, mat.liquid.sigma), spray_gates="2026-11-10")
    body = cube_body(loads=gl, spray_model=model, T_ambient=200.0, emissivity=0.2)
    started = time.perf_counter()
    book = body.advance(0.1, state=gas_state(rho=1e-4))
    print("gas state: {:.1f} s, {} parts, {} passes; heat in {:.3f} J, radiated {:.3e} J, sprayed {:.3e} kg; residuals "
          "mass {:.2e}, enthalpy {:.2e}; momentum to the wall {}".format(
              time.perf_counter() - started, book.parts, book.passes, book.enthalpy_in, book.enthalpy_radiated,
              book.mass_sprayed, book.mass_residual, book.enthalpy_residual, book.momentum_to_wall))
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8
    assert book.enthalpy_in > 0.0 and book.enthalpy_radiated > 0.0
    assert max(owed(body, book)) <= OWED_SHARE
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


def test_a_trial_pass_leaves_the_curvatures_as_they_were():
    """Fix round 1 of Task 9 (review, Important 1): the curvature cap runs in every heat and flow sub-step, so a trial
    pass may lower a column's curvature as its skin swells; the trial is undone with the curvatures too. The heated face's
    columns start exactly at their cap; a trial part of 0.0625 s at 6 MW/m2 melts and swells them, lowering their cap;
    after the restore every curvature is bit for bit the one before."""
    body = cube_body()
    sk, mk = body.skins, body.markers
    face = np.zeros(sk.n, bool)
    face[mk.skin_index[mk.normal[:, 0] > 0.9]] = True
    T, valid = sk.T(), sk.valid
    vol = np.zeros(sk.m.shape)
    vol[valid] = sk.m[valid] / sk.density(T[valid])
    depth = vol.sum(axis=1) / sk.area
    sk.kappa[face, 0] = skinbody.KAPPA_DEPTH_SHARE / depth[face]          # at the cap
    ctx = {"q": np.where(face, 6e6, 0.0), "forcing": None, "plan": None}
    S, _, _, g = body._interface()
    T_b = S @ body.solver.patch_T()
    kappa0 = sk.kappa.copy()
    saved = body._save(ctx)
    body._pass(0.0625, T_b, g, sk.minimum_mass(), ctx, trial=True)
    lowered = sk.kappa[face, 0] < kappa0[face, 0]
    body._restore(saved, ctx, solver=False)
    print("trial pass: {} of {} face columns had their curvature lowered".format(int(lowered.sum()), int(face.sum())))
    assert lowered.any()
    assert np.array_equal(sk.kappa, kappa0)


def gas_body():
    """The cube with the shared heating and surface flow (`loads.GasLoads`) and Girin's spray model, radiating."""
    from reentry_model.levelset import loads
    from reentry_model.shared import heating, spray, surface_flow
    from reentry_model.shared.material import LiquidProperties
    mat = aa7075()
    gl = loads.GasLoads(heating.PhysicsHeating(), surface_flow.SurfaceFlow(), mat.liquid.rho)
    model = spray.SprayModel(LiquidProperties(mat.liquid.rho, mat.liquid.mu, mat.liquid.sigma), spray_gates="2026-11-10")
    return cube_body(loads=gl, spray_model=model, T_ambient=200.0, emissivity=0.2)


PART_BOUND = 100      # parts a flow-on step of the test cube may take at most: a tenth of PARTS_MAX (measured 3 and 18
                      # parts; the skyscrapers' rule asked 5 and then 6 448)
CONTINUUM = dict(rho=1e-4, kn=1e-3)   # 7.5 km/s at a density of 1e-4 kg/m3 with a body Knudsen number of 1e-3: the shock
                                      # layer's branch, where Girin's closure gives a conjugate depth (at the merged
                                      # branch's Knudsen number of 0.03 the gates hold every marker: measured)


@pytest.fixture(scope="module")
def sprayed():
    """The cube with liquid formed on its +x face by 0.25 s of 6 MW/m2 (no flow, nothing removed), then one 0.05 s macro
    step at the continuum gas state with the runoff and Girin's plan on."""
    body = gas_body()
    melt = body.advance(0.25, q_prescribed=one_face(6e6))
    started = time.perf_counter()
    book = body.advance(0.05, state=gas_state(**CONTINUUM))
    return body, melt, book, time.perf_counter() - started


def test_liquid_formed_on_the_face_sprays_and_the_books_close(sprayed):
    """Fix round 1 of Task 9 (ruling b): liquid formed on the +x face sprays in the next 0.05 s macro step at the
    continuum gas state: the released mass is positive and equals the plan's books and the history's, the mass and
    enthalpy close, and the heat still owed and the interface mismatch stay within OWED_SHARE of the heat in."""
    body, melt, book, seconds = sprayed
    plan = body.plan
    liquid = melt.summary["liquid_depth_mean_mm"]
    print("spraying face: liquid {:.3f} mm mean before; {:.1f} s, {} parts; sprayed {:.3e} kg from {} markers; heat in "
          "{:.3f} J; residuals mass {:.2e}, enthalpy {:.2e}; droplets {:.0f}, median radius {:.1f} um".format(
              liquid, seconds, book.parts, book.mass_sprayed, int((plan.rate > 0).sum()), book.enthalpy_in,
              book.mass_residual, book.enthalpy_residual, book.summary["n_released"], book.summary["r_median_um"]))
    assert liquid > 0.0 and book.mass_sprayed > 0.0
    assert plan.mass.sum() + plan.mass_slurry.sum() == pytest.approx(book.mass_sprayed, rel=1e-12)
    assert book.summary["released_step_mass_kg"] == pytest.approx(book.mass_sprayed, rel=1e-12)
    assert abs(book.mass_residual) <= 1e-12 and abs(book.enthalpy_residual) <= 1e-8
    assert max(owed(body, book)) <= OWED_SHARE


def test_two_flow_on_steps_keep_their_part_count_bounded(sprayed):
    """Fix round 1 of Task 9 (ruling a), the smallest case found that reproduces the skyscrapers' runaway: the spraying
    step above and a second 0.05 s step at the same state. Measured once each: with the nearest-base patch groups 6 and
    24 parts, where the skyscrapers' rule (every accepted part's move ratio, scaled to the macro step) asked 31 and then
    31 521, far past PARTS_MAX; with the groups of fix round 2, 3 and 18 parts, their rule asking 5 and then 6 448. A
    part that meets its allowance only after splits, or never, carries the splits' shortness into their estimate. Here
    a step's hint is at most 1.25 times the parts it took: both steps stay within PART_BOUND parts, the second step's
    hint is below the skyscrapers' for it, both hints are 1.25 times their step's parts at most, and the books
    close."""
    body, _, first, _ = sprayed
    assert first.parts_next <= max(1, int(np.ceil(skinbody.MOVE_HINT_MARGIN * first.parts)))
    second = body.advance(0.05, state=gas_state(**CONTINUUM))
    print("two flow-on steps: parts {} then {}; hints {} then {}; the skyscrapers' hints {} then {}".format(
        first.parts, second.parts, first.parts_next, second.parts_next, first.parts_skyscrapers,
        second.parts_skyscrapers))
    assert first.parts <= PART_BOUND and second.parts <= PART_BOUND and second.parts_next < second.parts_skyscrapers
    assert second.parts_next <= max(1, int(np.ceil(skinbody.MOVE_HINT_MARGIN * second.parts)))
    assert abs(second.mass_residual) <= 1e-12 and abs(second.enthalpy_residual) <= 1e-8
    assert max(owed(body, second)) <= OWED_SHARE



def test_a_marker_whose_nearest_patches_are_few_gets_a_grown_group_and_the_law_stays_exact():
    """Fix round 2 of Task 9 (the ruling on concern 1): a marker's base law is spread over a group of patches of beta
    whose area is within GROUP_AREA_RANGE of its skin's base area. Flat patches of 0.1 x 0.1 mm on a plane carry five
    base points of base area 1 mm2 each, one at the origin ringed by four 0.02 mm away, so that by nearest base it owns
    almost nothing (the full-size failure of fix round 1: a whole law on a 0.1 mm2 group, a Biot number of about 800). Every group ends within range, the crowded marker's grown from its nearest patches, and the law each patch
    receives sums back, over the patches, to every marker's a and b exactly (to round-off)."""
    h = 0.1e-3
    ix, iy = np.meshgrid(np.arange(-15, 26), np.arange(-15, 16), indexing="ij")
    centroid = np.column_stack([(ix.ravel() + 0.5) * h, (iy.ravel() + 0.5) * h, np.zeros(ix.size)])
    area = np.full(len(centroid), h * h)
    d = 2e-5
    base = np.array([[0.0, 0.0, 0.0], [d, 0.0, 0.0], [-d, 0.0, 0.0], [0.0, d, 0.0], [0.0, -d, 0.0]])
    A_base = np.full(len(base), 1e-6)
    near = skinbody.patch_groups(centroid, area, base, np.zeros(len(base)))[1]      # nearest base alone, no target
    S, A_group = skinbody.patch_groups(centroid, area, base, A_base)
    ratio = A_group / A_base
    print("patch groups: nearest-base areas {} mm2, grown {} mm2".format(np.round(near * 1e6, 3),
                                                                          np.round(A_group * 1e6, 3)))
    assert near[0] < 0.1 * A_base[0]
    lo, hi = skinbody.GROUP_AREA_RANGE
    assert np.all(ratio >= lo) and np.all(ratio <= hi)
    np.testing.assert_allclose(np.asarray(S.sum(axis=1)).ravel(), 1.0, rtol=1e-14)
    a, b = np.array([3.0, -2.0, 5.0, 1.0, 0.5]), np.array([-16.0, -0.4, -0.1, -2.0, -1.0])
    q0, q1 = (S.T @ a) / area, (S.T @ b) / area
    assert (q0 * area).sum() == pytest.approx(a.sum(), rel=1e-14)
    assert (q1 * area).sum() == pytest.approx(b.sum(), rel=1e-14)


def test_the_concave_cap_keeps_a_concave_column_from_drawing_without_end():
    """Fix round 2 of Task 9 (the ruling on concern 3): a column concave enough that |kappa| stub > 1 (kappa = -2000 1/m
    under a 0.8 mm stub, |kappa| stub = 1.6) widens with depth so fast that its minimum mass, rho_s A(z_bottom) stub, grows
    faster than its mass: each draw deepens it and raises its minimum again. Without the cap, twenty draws more than
    double its mass; with its curvatures held by `cap_curvature` (re-capped before every draw, as the body does), its
    minimum mass stays below its mass after the first fill and the twenty draws add under a percent."""
    mat = aa7075()
    stub = SETTINGS.stub

    def column():
        sk = Skins(mat, np.full(1, 1e-6), np.array([[-2000.0, 0.0]]), SETTINGS)
        sk.fill(mat.rho * 1e-6 * stub, 850.0)
        return sk

    free = column()
    m0 = free.m.sum()
    for _ in range(20):
        free.draw(850.0, 0.0, free.minimum_mass())
        free.regrid()
    held = column()
    held.kappa = skinbody.cap_curvature(held, held.kappa)
    held.draw(850.0, 0.0, held.minimum_mass())
    held.regrid()
    m1, added = held.m.sum(), 0.0
    for _ in range(20):
        held.kappa = skinbody.cap_curvature(held, held.kappa)
        dm, _ = held.draw(850.0, 0.0, held.minimum_mass())
        added += float(dm.sum())
        held.regrid()
    print("concave column: without the cap {:.2f} times its mass after twenty draws; with it {:.2e} of it added, "
          "minimum mass {:.4f} of the mass".format(free.m.sum() / m0, added / m1, held.minimum_mass()[0] / held.m.sum()))
    assert abs(-2000.0) * stub > 1.0
    assert free.m.sum() > 2.0 * m0
    assert added < 0.01 * m1 and held.minimum_mass()[0] <= held.m.sum() * (1 + 1e-12)


def _hints(move, hint, steps=5, dt=0.5):
    """The part-count hints over consecutive macro steps of length dt [s] driven by the body's own splitting and hint
    rules (`skinbody.split_part`, `skinbody.next_hint`) with a synthetic base move: `move(t0, t1)` [-] the largest move
    ratio of a part from t0 to t1 [s] into the step. Returns the hints, the first being the one the steps start with."""
    out = [hint]
    for _ in range(steps):
        clock = {"t": 0.0}

        def run(h):
            t0 = clock["t"]

            def accept():
                clock["t"] = t0 + h
            return move(t0, t0 + h), accept, lambda: None
        counts = [skinbody.split_part(run, dt / hint, dt, 0 if hint == 1 else 1) for _ in range(hint)]
        assert clock["t"] == pytest.approx(dt, rel=1e-12)
        hint = skinbody.next_hint(counts)
        out.append(hint)
    return out


def test_a_move_that_does_not_shrink_with_the_part_does_not_raise_the_part_count():
    """A jump of three allowances at each step's start and a refill of a twentieth of the allowance in every part, neither
    shrinking with the part's length (the full-size cube's hand-down of a backlog, fix round 5): splitting cannot help,
    the parts it produced do not count, and the hint falls to one and stays there. (Counting every accepted part, as
    before, raised it by about 1.25 each step: the skyscrapers' stall.)"""
    hints = _hints(lambda t0, t1: (3.0 if t0 == 0.0 else 0.0) + 0.05, 4)
    print("fixed jump: hints", hints)
    assert all(b <= a for a, b in zip(hints, hints[1:]))
    assert hints[-1] == 1


def test_a_move_that_scales_with_the_part_sets_the_part_count_and_holds_it():
    """5.5 allowances of base move per macro step, spread evenly: from one part the hint grows to ceil(1.25 x 5.5) = 7,
    what the move asks for, and holds there."""
    hints = _hints(lambda t0, t1: 5.5 * (t1 - t0) / 0.5, 1)
    print("scaling move: hints", hints)
    assert hints[1:] == [7] * 5


def test_small_moves_let_the_part_count_fall():
    """A fifth of an allowance per macro step: from 20 parts the hint falls to one and stays."""
    hints = _hints(lambda t0, t1: 0.2 * (t1 - t0) / 0.5, 20)
    print("small moves: hints", hints)
    assert hints[1:] == [1] * 5


def test_a_piece_the_draws_use_up_leaves_and_the_largest_piece_gives_the_rest():
    """Two pieces of beta, a 1.2 mm sphere beside a 0.3 mm one, and draws at a point in the small one asking twice its
    volume (as on the full-size cube's third step, where three markers over a 1.9e-10 m3 sliver drew 2.0e-9 m3): the
    small piece leaves, the largest gives the rest, the total volume changes by exactly the draws and the stored
    enthalpy by what the move returns (fix round 5)."""
    from materials import ConstMaterial
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.cutcell import inside_volumes
    from reentry_model.levelset.cutthermal import CutThermal
    from reentry_model.levelset.levelset import pieces, redistance
    from reentry_model.levelset.stepend import StepEnd

    class Host(StepEnd):
        pass
    mesh = build(BoxBase.around((0.0, 0.0, 0.0), (4e-3, 3e-3, 3e-3), 0.125e-3))
    x = mesh.points
    big, small = np.array([1.4e-3, 1.5e-3, 1.5e-3]), np.array([3.1e-3, 1.5e-3, 1.5e-3])
    beta = redistance(mesh, np.minimum(np.linalg.norm(x - big, axis=1) - 1.2e-3, np.linalg.norm(x - small, axis=1) - 0.3e-3))
    label, n = pieces(mesh, beta)
    v = inside_volumes(mesh, beta)
    assert n == 2
    host = Host()
    host.mesh, host.beta, host.band, host.fine_margin = mesh, beta, 1e-3, 1e-3
    host.solver = CutThermal(ConstMaterial())
    host.solver.setup(mesh, beta)
    host.solver.set_temperature(500.0)
    tally = skinbody._Tally(0)
    E0, V0 = host.solver.energy(), v.sum()
    dV = np.array([-2.0 * v[label == 1].sum()])
    removed = host._move_beta(small[None, :], dV, np.zeros(1), tally)
    after = inside_volumes(mesh, host.beta)
    assert pieces(mesh, host.beta)[1] == 1 and tally.pieces_used == 1
    assert after.sum() == pytest.approx(V0 + dV.sum(), rel=1e-12, abs=0.0)
    assert removed == pytest.approx(E0 - host.solver.energy(), rel=1e-12)


def test_draws_asking_more_than_the_largest_piece_holds_raise():
    """Task 9b (minor, left from Task 9): draws at a point inside the only piece of beta, a 1.2 mm sphere, asking twice
    its volume: no piece is left to give the rest, so the move raises, naming the largest piece, and leaves beta as it
    was."""
    from materials import ConstMaterial
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.cutcell import inside_volumes
    from reentry_model.levelset.cutthermal import CutThermal
    from reentry_model.levelset.levelset import redistance
    from reentry_model.levelset.stepend import StepEnd

    class Host(StepEnd):
        pass
    mesh = build(BoxBase.around((0.0, 0.0, 0.0), (3e-3, 3e-3, 3e-3), 0.125e-3))
    centre = np.array([1.5e-3, 1.5e-3, 1.5e-3])
    beta = redistance(mesh, np.linalg.norm(mesh.points - centre, axis=1) - 1.2e-3)
    host = Host()
    host.mesh, host.beta, host.band, host.fine_margin = mesh, beta, 1e-3, 1e-3
    host.solver = CutThermal(ConstMaterial())
    host.solver.setup(mesh, beta)
    host.solver.set_temperature(500.0)
    dV = np.array([-2.0 * inside_volumes(mesh, beta).sum()])
    with pytest.raises(RuntimeError, match="largest piece"):
        host._move_beta(centre[None, :], dV, np.zeros(1), skinbody._Tally(0))
    assert np.array_equal(host.beta, beta)


def test_the_snap_books_the_markers_its_last_round_brings_from_afar_and_their_slide():
    """Task 9b (minor, left from Task 9): phi = 0.05 (|x| - 2 mm), rising only 0.05 per metre along its normal, and
    three markers 1 mm outside its zero contour: each recovered-normal step removes a twentieth of the distance, so 19
    rounds leave them 1 mm x 0.95^38 = 0.14 mm off, where |phi| = 7e-6 is beyond the Newton steps' reach (1e-6); the last
    round brings all three on by Newton steps from afar. The step's tally counts the three (`snap_fallback`) and holds
    the largest sideways part of those Newton steps (`snap_slide`), as a replay of the rounds measures them, and every
    marker ends on phi."""
    from reentry_model.levelset import markers as M
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.markers import Markers
    from reentry_model.levelset.stepend import SNAP_ROUNDS, SNAP_TOL, StepEnd

    class Host(StepEnd):
        pass
    mesh = build(BoxBase.around((-4e-3,) * 3, (4e-3,) * 3, 1e-3))
    phi = 0.05 * (np.linalg.norm(mesh.points, axis=1) - 2e-3)
    u = np.array([[1.0, 0.2, 0.1], [-0.3, 1.0, 0.4], [0.2, -0.5, -1.0]])
    u /= np.linalg.norm(u, axis=1)[:, None]
    x0 = 3e-3 * u

    def markers():
        return Markers(x0, u, np.ones(3), np.zeros((3, 2)), np.arange(3), np.full(3, -1), 1e-3)
    host = Host()
    host.mesh, host.phi, host.markers = mesh, phi, markers()
    tally = skinbody._Tally(0)
    host._snap(tally)
    # the replay: the first rounds within the reach, then the last round's recovered passes and Newton steps by hand
    mk = markers()
    for _ in range(SNAP_ROUNDS - 1):
        mk.snap(mesh, phi)
    x = mk.x
    for _ in range(M.SNAP_PASSES):
        val, n = M._field(mesh, phi, x)
        x = x - val[:, None] * n
    val, _ = M._field(mesh, phi, x)
    far = int((np.abs(val) > M.NEWTON_REACH).sum())
    after = M._cell_newton(mesh, phi, x)
    _, n = M._field(mesh, phi, after)
    d = after - x
    slide = np.linalg.norm(d - np.einsum("ij,ij->i", d, n)[:, None] * n, axis=1).max()
    off = np.abs(M._field(mesh, phi, host.markers.x)[0]).max()
    print("snap from afar: {} markers by the fallback, largest slide {:.3e} m (replay {:.3e} m), off phi {:.1e}".format(
        tally.snap_fallback, tally.snap_slide, slide, off))
    assert far == 3 and tally.snap_fallback == 3
    assert tally.snap_slide == slide and slide > 0.0
    assert off <= SNAP_TOL


def test_the_test_cube_s_markers_keep_their_neighbours():
    """Task 9b (the ruling on the neighbour test): on the test cube as the body builds it (1.4 mm fine cells, 0.7 mm at
    the rounded edges, markers one per 2 mm voxel) under 1 % of the cells are bounded by the cap and every area is within
    a factor of 1.5 of the area nearest the marker. Measured before the ruling: 123 of 219 cells capped, areas 0.38 to
    2.67 times the nearest-marker areas."""
    body = cube_body()
    mk = body.markers
    mk.tessellate(body.mesh, body.phi)
    r = nearest_ratio(mk, body.mesh, body.phi)
    print("test cube: {} markers, {} capped; areas {:.3f} to {:.3f} times the nearest-marker areas".format(
        mk.n, int(mk.capped.sum()), r.min(), r.max()))
    assert mk.capped.mean() < CAPPED_SHARE
    assert r.min() >= 1.0 / OWNERSHIP_FACTOR and r.max() <= OWNERSHIP_FACTOR



# -- the step end's re-partition (Task 9b) -----------------------------------------------------------------------------
BOOKS_RTOL = 1e-14    # the skins' mass and enthalpy over a step end, relative (Task 9b)
SLAB_HALF = (6e-3, 6e-3, 2e-3)    # m: the flat box's half sides; its top face is the flat surface the markers move on


def _slab_host(mat, T=850.0):
    """A host of the step end (`stepend.StepEnd`) on a closed flat box (a slab 12 x 12 x 4 mm, its centre OFF off the
    grid; 1 mm cells): on its top face 36 markers on a regular grid 2 mm apart (x, y = OFF -5, -3, ..., 5 mm), on its
    other faces markers seeded one per 2 mm voxel as the body seeds them. Each skin starts a stub of solid metal at T
    [K]; the host is settled until a step end splits and merges none (at most five times). beta lies the stub below
    phi. Returns (host, the top grid's markers (6, 6) by x then y)."""
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.cutcell import patches
    from reentry_model.levelset.markers import Markers
    from reentry_model.levelset.stepend import StepEnd, cap_curvature

    class Host(StepEnd):
        pass
    h = 2e-3
    mesh = build(BoxBase.around((-8e-3, -8e-3, -4e-3), (8e-3, 8e-3, 4e-3), 1e-3))
    phi = np.max(np.abs(mesh.points - OFF) - np.array(SLAB_HALF), axis=1)
    P = patches(mesh, phi)
    side = np.flatnonzero(P.normal[:, 2] < 0.5)
    key = np.floor((P.centroid[side] + 8e-3) / h).astype(np.int64)
    order = np.lexsort((-P.area[side], key[:, 2], key[:, 1], key[:, 0]))
    k = key[order]
    seeds = side[np.sort(order[np.r_[True, np.any(k[1:] != k[:-1], axis=1)]])]
    g = OFF + (np.arange(6) - 2.5) * h
    X, Y = np.meshgrid(g, g, indexing="ij")
    grid = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, OFF + SLAB_HALF[2])])
    x = np.concatenate([grid, P.centroid[seeds]])
    n = len(x)
    normal = np.concatenate([np.tile([0.0, 0.0, 1.0], (len(grid), 1)), P.normal[seeds]])
    sk = Skins(mat, np.ones(n), np.zeros((n, 2)), SETTINGS)
    mk = Markers(x, normal, np.ones(n), np.zeros((n, 2)), np.arange(n), np.full(n, -1), h, skins=sk)
    mk.tessellate(mesh, phi)
    sk.area[mk.skin_index] = mk.area
    sk.kappa[mk.skin_index] = cap_curvature(sk, mk.kappa, np.full(n, SETTINGS.stub))
    sk.fill(sk.minimum_mass(), T)
    host = Host()
    host.mesh, host.phi, host.beta, host.markers, host.skins, host.spacing = mesh, phi, phi + SETTINGS.stub, mk, sk, h
    host.band, host.fine_margin, host.material, host.rho = 4e-3, 2e-3, mat, mat.rho
    for _ in range(5):
        changes = mk.splits + mk.merges
        host._settle(skinbody._Tally(0))
        if mk.splits + mk.merges == changes:
            break
    else:
        raise RuntimeError("the slab's markers did not settle")
    top = np.array([[np.flatnonzero(np.all(np.isclose(mk.x, p, rtol=0.0, atol=1e-12), axis=1))[0] for p in row]
                    for row in grid.reshape(6, 6, 3)])
    return host, top


def _flat_depth(sk):
    """Each skin row's flat depth [m]: its volume, the cells' m / rho(T), over its top area (what mixing columns at one
    density averages by area; the thickness under a curved cross-section is not)."""
    T, valid = sk.T(), sk.valid
    vol = np.zeros(sk.m.shape)
    vol[valid] = sk.m[valid] / sk.density(T[valid])
    return vol.sum(axis=1) / sk.area


def _through_base(host):
    """The volume [m3] the next part's draws and hand-downs would move through the skins' base (`skinbody`'s
    `_part_once`: draws to the minimum mass less DRAW_SLACK, hand-downs keeping the minimum mass), on a copy."""
    import copy
    sk = copy.deepcopy(host.skins)
    drawn, _ = sk.draw(850.0, 0.0, sk.minimum_mass() * (1.0 - skinbody.DRAW_SLACK))
    handed, _, _ = sk.hand_down(SETTINGS.stub, m_keep=sk.minimum_mass())
    return float(np.abs(handed - drawn).sum()) / host.rho, handed, drawn


def test_a_repartition_moves_no_metal_through_the_base():
    """Task 9b: liquid skins of one to three stubs' mass on a closed flat box, the 16 markers in the middle of its top
    face shifted along it by up to 0.4 mm each, the surface unchanged: the step end's tessellation re-partitions the area
    (the shifted cells' areas change by 0.6 to 1.13 times) and each column's metal follows its area by the hand-over
    rule. A third of the shifted markers' skins (flat columns) start exactly at their minimum mass, the others above it. The skins' mass
    and enthalpy close to 1e-14 and their momentum to round-off; no column ends below its minimum mass, so the next
    part's draws and hand-downs move no metal through the base (at most 1e-14 of beta's volume); every new column's
    depth (its flat depth: its volume over its top area) lies within its donors'; the stretch is one, since the surface
    did not change. (Before Task 9b each skin
    kept its metal on its new area: the columns that grew fell below their minimum and drew.)"""
    from materials import LinearMeltMaterial
    from reentry_model.levelset.cutcell import inside_volumes
    mat = LinearMeltMaterial()
    host, top = _slab_host(mat)
    mk, sk = host.markers, host.skins
    middle = top[1:5, 1:5].ravel()
    rng = np.random.default_rng(7)
    factor = rng.uniform(1.5, 3.0, mk.n)
    factor[middle] = np.where(rng.random(len(middle)) < 1.0 / 3.0, 1.0, rng.uniform(1.0, 3.0, len(middle)))
    row = np.argsort(mk.skin_index)                                    # each skin row's marker
    sk.fill(factor[row] * sk.minimum_mass(), 950.0)                    # liquid: above the liquidus, 908 K
    sk.U = np.where(sk.valid[..., None], np.concatenate([rng.normal(0.0, 1.0, (sk.n, 1, 2)),
                                                         np.zeros((sk.n, 1, 1))], axis=2), 0.0)
    M0, E0, P0 = (np.sum(t, axis=0) for t in sk.totals())
    P_size = np.abs(sk.m[..., None] * sk.U).sum()
    thick0, A0, n0 = _flat_depth(sk)[mk.skin_index], mk.area.copy(), mk.n
    mk.x[middle, :2] += rng.uniform(-0.4e-3, 0.4e-3, (len(middle), 2))
    host._settle(skinbody._Tally(0))
    assert mk.n == n0
    M1, E1, P1 = (np.sum(t, axis=0) for t in sk.totals())
    moved, handed, drawn = _through_base(host)
    V = inside_volumes(host.mesh, host.beta).sum()
    ratio = mk.area[middle] / A0[middle]
    print("re-partition: shifted areas {:.3f} to {:.3f} of their old; books mass {:.1e} enthalpy {:.1e} momentum {:.1e}; "
          "through the base {:.1e} of beta's volume; least mass over minimum {:.15f}".format(
              ratio.min(), ratio.max(), M1 / M0 - 1, E1 / E0 - 1, np.abs(P1 - P0).max() / P_size, moved / V,
              (sk.m.sum(axis=1) / sk.minimum_mass()).min()))
    assert ratio.min() < 0.8 and ratio.max() > 1.1
    assert np.all(sk.m.sum(axis=1) >= sk.minimum_mass() * (1.0 - 1e-13))
    assert moved <= 1e-14 * V
    assert abs(M1 / M0 - 1) <= BOOKS_RTOL and abs(E1 / E0 - 1) <= BOOKS_RTOL
    assert np.abs(P1 - P0).max() <= BOOKS_RTOL * P_size
    T, stretch = host.last_settle[0].transfer, host.last_settle[0].stretch
    assert np.abs(stretch - 1.0).max() <= 1e-14
    thick1 = _flat_depth(sk)[mk.skin_index]
    for i in range(mk.n):
        donors = T[:, i].nonzero()[0]
        assert thick0[donors].min() * (1 - 1e-12) <= thick1[i] <= thick0[donors].max() * (1 + 1e-12)


def test_a_uniform_stretch_still_draws(monkeypatch):
    """Task 9b: the markers of a 10 mm sphere on the design's 0.7 mm cells, one per 2 mm voxel, their skins solid columns
    exactly at their minimum mass, carried radially onto a 10.4 mm sphere: the step end finds the surface stretched by
    (10.4/10)^2 (the stretch, sampled from the old regions and smoothed over each marker's neighbours, within 1 % of it
    for every marker) and re-partitions under 2 % of the area; every column's flat depth falls by its stretch, every
    column falls short of its new minimum, and the next part draws, in all, the new minimum masses less the metal there
    is (section 9.2: a stretched skin thins and draws). The skins' books close to 1e-14. (Split and merge are switched off:
    the stretch takes some cells past twice the target area, and the markers are compared one by one.)"""
    from materials import ConstMaterial
    from reentry_model.levelset import markers
    monkeypatch.setattr(markers, "SPLIT_FACTOR", np.inf)
    monkeypatch.setattr(markers, "MERGE_FRACTION", 0.0)
    from reentry_model.levelset.boxmesh import build
    from reentry_model.levelset.repartition import surface_samples
    from reentry_model.levelset.stepend import StepEnd, cap_curvature

    class Host(StepEnd):
        pass
    R0, R1 = 10e-3, 10.4e-3
    mat = ConstMaterial()
    base = BoxBase.around((-R1 - 5.6e-3,) * 3, (R1 + 5.6e-3,) * 3, 2.8e-3)
    r = lambda x: np.linalg.norm(x, axis=1)
    mesh = build(base, lambda x, s: np.where((r(x) > R0 - 1.4e-3 - 0.87 * s) & (r(x) < R1 + 1.4e-3 + 0.87 * s),
                                             DESIGN_FINE, np.inf))
    phi0, phi1 = r(mesh.points) - R0, r(mesh.points) - R1
    mk = seeded(mesh, phi0, DESIGN_SPACING, base.lo, material=mat)
    mk.tessellate(mesh, phi0)
    sk = mk.skins
    sk.area[mk.skin_index] = mk.area
    sk.kappa[mk.skin_index] = cap_curvature(sk, mk.kappa, np.full(mk.n, SETTINGS.stub))
    sk.fill(mat.rho * sk.area * SETTINGS.stub, 600.0)
    for _ in range(60):                     # each column at its minimum mass: the fixed point of M = rho A(z_bottom(M)) stub
        sk.fill(sk.minimum_mass(), 600.0)
    assert np.allclose(sk.m.sum(axis=1), sk.minimum_mass(), rtol=1e-14, atol=0.0)
    old = surface_samples(mesh, phi0, DESIGN_SPACING / 16.0)
    M0, E0 = sk.m.sum(), sk.H.sum()
    A0, min0, thick0 = mk.area.copy(), sk.minimum_mass().sum(), _flat_depth(sk)[mk.skin_index]
    mk.x = mk.x * (R1 / R0)
    host = Host()
    host.mesh, host.phi, host.beta, host.markers, host.skins, host.spacing = mesh, phi1, phi1 + SETTINGS.stub, mk, sk, \
        DESIGN_SPACING
    host.band, host.fine_margin, host.material, host.rho = 4e-3, 2e-3, mat, mat.rho
    host._settle(skinbody._Tally(0), old)
    s2 = (R1 / R0) ** 2
    log = host.last_settle[0]
    T, stretch = log.transfer, log.stretch
    moved = (T.sum() - T.diagonal().sum()) / T.sum()
    m, m_min, E1, thick1 = sk.m.sum(axis=1), sk.minimum_mass(), sk.H.sum(), _flat_depth(sk)[mk.skin_index]
    drawn, _ = sk.draw(600.0, 0.0, m_min)
    print("uniform stretch: {} markers, stretch {:.4f} to {:.4f} against {:.4f}; re-partitioned {:.2%} of the area; "
          "drew {:.4e} kg against the minimum masses' rise {:.4e} kg; books {:.1e}".format(
              mk.n, stretch.min(), stretch.max(), s2, moved, drawn.sum(), m_min.sum() - min0, M0 / m.sum() - 1))
    assert np.abs(stretch / s2 - 1.0).max() < 0.01 and moved < 0.02
    assert abs(m.sum() / M0 - 1) <= BOOKS_RTOL and abs(E1 / E0 - 1) <= BOOKS_RTOL
    assert np.all(m < m_min)
    assert drawn.sum() == pytest.approx(m_min.sum() - M0, rel=1e-12)
    assert np.median(thick1 / thick0) == pytest.approx(1.0 / s2, rel=0.01)


def test_a_thick_skin_losing_area_gives_its_metal_to_its_neighbour():
    """Task 9b: on the closed flat box every skin is a solid stub at 850 K but one in the middle of the top face, 2 mm
    deep, its stub under 1.2 mm of liquid at 950 K. Its two neighbours along x close in on it, each by half a spacing,
    so its area halves and a quarter of its old area passes to each. Its metal follows: its own column keeps its depth
    (it does not thicken, so it hands nothing down), and the heat its liquid holds above the stubs' (its enthalpy less
    that of the same mass at 850 K, which adds when columns mix) goes, column by column, exactly where the transfer
    sends its area: most of it to the two neighbours, none to beta. (Before Task 9b the column kept its metal on half
    the area and stood twice as deep, its rigid metal below the stub handed down to beta at the next part.)"""
    from materials import LinearMeltMaterial
    mat = LinearMeltMaterial()
    host, top = _slab_host(mat)
    mk, sk = host.markers, host.skins
    h, stub = host.spacing, SETTINGS.stub
    thick, left, right = top[2, 3], top[1, 3], top[3, 3]
    row = mk.skin_index[thick]
    n_l = 60
    K = int(sk.count.max())
    m, H = np.zeros((sk.n, K + n_l)), np.zeros((sk.n, K + n_l))
    m[:, :K], H[:, :K] = sk.m, sk.H
    count = sk.count.copy()
    c = count[row]
    m[row, n_l:n_l + c], H[row, n_l:n_l + c] = sk.m[row, :c], sk.H[row, :c]
    m[row, :n_l] = mat.liquid.rho * sk.area[row] * 1.2e-3 / n_l
    H[row, :n_l] = m[row, :n_l] * mat.enthalpy(950.0)
    count[row] = c + n_l
    sk.install(m, H, count)

    def excess():
        """Each marker's enthalpy above the same mass at 850 K [J] (additive when columns mix)."""
        return (sk.H.sum(axis=1) - sk.m.sum(axis=1) * mat.enthalpy(850.0))[mk.skin_index]
    E_ex, M0, depth0, A0, n0 = excess(), sk.m.sum(), sk.thickness()[row], mk.area[thick], mk.n
    mk.x[left, 0] += 0.5 * h
    mk.x[right, 0] -= 0.5 * h
    host._settle(skinbody._Tally(0))
    assert mk.n == n0
    depth1, after = sk.thickness()[row], excess()
    _, handed, _ = _through_base(host)
    print("thick skin: area {:.3f} of its old, depth {:.4f} of its old, handed down {:.2e} kg; its liquid's heat now "
          "{:.3f} in itself, {:.3f} and {:.3f} in its neighbours".format(
              mk.area[thick] / A0, depth1 / depth0, handed[row], after[thick] / E_ex[thick],
              after[left] / E_ex[thick], after[right] / E_ex[thick]))
    assert mk.area[thick] / A0 == pytest.approx(0.5, abs=0.02)
    assert depth1 <= depth0 * (1.0 + 1e-12) and depth1 >= 0.95 * depth0
    assert handed[row] == 0.0
    assert after[left] > 0.2 * E_ex[thick] and after[right] > 0.2 * E_ex[thick]
    assert abs(sk.m.sum() / M0 - 1.0) <= BOOKS_RTOL
    tr = host.last_settle[0].transfer
    share = np.asarray(tr[thick].todense()).ravel() / tr[thick].sum()
    np.testing.assert_allclose(after, E_ex[thick] * share, rtol=0.0, atol=1e-12 * E_ex[thick])
    assert share[left] + share[right] >= 0.9 * (1.0 - share[thick])


def test_a_marker_the_second_tessellation_leaves_without_area_is_released(monkeypatch):
    """Fix round 1 (minor 1): the step end's tessellation after split and merge may leave a marker no area, as the first
    may; it is released (its skin's metal to the volume) before its skins are re-partitioned, as after the first. Here
    the test cube's step end reports a merge, so that it tessellates a second time, and that tessellation leaves one
    marker no area: the marker leaves, counted as released, and the body's mass closes (before, the re-partition ran
    with the zero-area column and the marker stayed, with no area)."""
    from reentry_model.levelset.markers import Markers
    body = cube_body()
    mk = body.markers
    calls = {"n": 0}
    tessellate = Markers.tessellate

    def tess(self, mesh, phi):
        faces = tessellate(self, mesh, phi)
        calls["n"] += 1
        if calls["n"] == 2:                               # the tessellation after split and merge
            self.area[0] = 0.0
        return faces
    monkeypatch.setattr(Markers, "tessellate", tess)
    monkeypatch.setattr(Markers, "merge", lambda self, spacing: 1)
    before = body.books()
    n0, released0 = mk.n, mk.trapped
    body._settle(skinbody._Tally(0))
    after = body.books()
    M0, M1 = before["mass_volume"] + before["mass_skins"], after["mass_volume"] + after["mass_skins"]
    print("second tessellation: {} markers then {}, released {}; mass {:.1e}".format(n0, mk.n, mk.trapped - released0,
                                                                                    M1 / M0 - 1))
    assert mk.n == n0 - 1 and mk.trapped == released0 + 1
    assert np.all(mk.area > 0.0)
    assert abs(M1 / M0 - 1) <= 1e-12
