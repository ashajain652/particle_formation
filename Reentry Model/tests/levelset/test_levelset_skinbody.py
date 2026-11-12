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
OWED_SHARE = 0.05      # the heat still owed to the volume at a step's end, summed in absolute value per point, and the last
                       # part's interface mismatch, each at most this share of the step's heat in (fix round 1, Important 2).
                       # Measured once each: 4e-4 on the heated face; 1.2e-2 on the receding face; 2.9e-2 on the spraying
                       # face's 0.05 s step. What is owed is almost all the base moves' deferral -- the draw values the
                       # metal it takes on the profile capped at T_rigid, the volume removes it at its own, hotter,
                       # temperatures -- about 4 % of the heat in over the receding step, most of it delivered within it
                       # (the mismatch itself is 1e-4 of the heat in or less)


def owed(body, book):
    """(the heat still owed to the volume, summed in absolute value [J], the last part's sum |mismatch| [J]) as shares of
    the step's heat in, printed with the step's deferred residuals by source."""
    pend, mism = np.abs(body.pending_J).sum() / book.enthalpy_in, book.interface_mismatch_J / book.enthalpy_in
    lo, hi = skinbody.GROUP_AREA_RANGE
    assert lo <= book.group_area_ratio_min and book.group_area_ratio_max <= hi
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
    assert max(owed(body, book)) <= OWED_SHARE
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
    from test_levelset_loads import _state
    body = gas_body()
    melt = body.advance(0.25, q_prescribed=one_face(6e6))
    started = time.perf_counter()
    book = body.advance(0.05, state=_state(**CONTINUUM))
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
    from test_levelset_loads import _state
    body, _, first, _ = sprayed
    assert first.parts_next <= max(1, int(np.ceil(skinbody.MOVE_HINT_MARGIN * first.parts)))
    second = body.advance(0.05, state=_state(**CONTINUUM))
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
