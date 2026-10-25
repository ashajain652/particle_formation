"""The skin body (layered-skin design sections 6.6 and 7, sub-plan 20 Task 7): skins on every corner from time zero, the
trial and real passes around the 3D step, the base moving the 3D mesh, and the books."""
import numpy as np
import pytest

from reentry_model import heating, material, mesh, movemesh, skin, skinbody
from reentry_model.thermal import thermal_solver

DEV = skin.SkinSettings(cell=20e-6)


def make(m, T0=880.0, name="skfem", settings=DEV, **kw):
    copy = mesh.VolumeMesh(m.points.copy(), m.tets.copy(), dict(m.params))
    mat = material.Material.from_drama_json("AA7075_scheil")
    return skinbody.SkinBody(copy, mat, thermal_solver(name), 1.47, settings, T0=T0, v_hat=(-1.0, 0.0, 0.0), **kw)


def front_loads(b, q0):
    return heating.HeatingResult(q0 * np.clip(np.cos(b.theta), 0.0, None), q0, 0.0, 0.0, 0.0)


def test_every_corner_starts_with_its_minimum_taken_from_the_3d_model(coarse_sphere_mesh):
    b = make(coarse_sphere_mesh)
    m_min = skin.minimum_mass(b.material, b.columns.area, DEV)
    V0 = coarse_sphere_mesh.volume()
    assert b.skins.n == b.columns.n == len(np.unique(coarse_sphere_mesh.surface().faces))
    assert b.mass() == pytest.approx(b.material.rho * V0, rel=1e-14)                  # nothing made or lost
    assert b.energy() == pytest.approx(b.material.rho * V0 * float(b.material.enthalpy(880.0)), rel=1e-12)
    # each column holds what its corner swept moving in by its minimum: within 3 % of it (prisms against dual areas
    # on 4 mm triangles: 2.4 % for the 0.8 mm minimum of option B', under 2 % for 0.4 mm), so the first sub-step draws
    # any shortfall
    assert np.abs(b.skins.column_mass() / m_min - 1.0).max() < 0.03
    assert np.allclose(b.skins.T[b.skins.valid], 880.0)
    assert b.surface.n_patches == b.columns.n and b.theta[b.i_stag] < np.radians(5.0)


def test_a_body_with_no_heat_stays_uniform_and_its_mesh_and_skins_keep_their_books(coarse_sphere_mesh):
    b = make(coarse_sphere_mesh, T0=600.0)
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.emissivity = 0.0
    M0, E0 = b.mass(), b.energy()
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, zero)
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert np.abs(b.solver.temperature() - 600.0).max() < 1e-6
    assert np.abs(b.skins.T[b.skins.valid] - 600.0).max() < 1e-6
    assert b.removed_mass == 0.0 and b.jump_faults == 0


def test_heating_with_removal_keeps_mass_and_energy_exact(coarse_sphere_mesh):
    """3 MW/m2 on the front for 1.5 s from 880 K: liquid forms and leaves, the skins absorb the non-rigid metal below
    them, the base recedes; the mass removed plus the mass left is the mass at start, and the energy books close to 1e-8
    of the heat absorbed. The three 0.5 s steps run as 1, 8 and 8 parts and 53 g is removed (measured 2026-10-14).
    3 MW/m2, not 6: on these 4 mm elements the rim of a strongly melting region is one edge wide, a corner there absorbs
    a 0.7 mm lens its neighbours do not have, and the mesh refuses the move; the 2 mm production mesh follows it
    (`test_strong_melting_on_the_production_mesh`, reference tier)."""
    b = make(coarse_sphere_mesh)
    M0 = b.mass()
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, front_loads(b, 3e6))
    assert b.removed_mass > 0.0 and b.drawn_mass > 0.0
    assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8                  # 2e-11 measured (2026-10-14)
    stats = b.skin_stats()
    assert stats["base_recession_max_mm"] > 0.5                       # 4.5 mm at the nose (measured, 2026-10-14)
    assert b.surface_temperature().max() <= b.material.T_feed + 1.0   # liquid left as it formed
    # absorption keeps up: no base is left more than 2 K above T_rigid, no non-rigid metal stranded in the 3D model
    # (0.3 K measured, 2026-10-14; 7 to 8 K under the old quarter-edge cap per macro step)
    assert b.solver.temperature()[b.moving.corners].max() - b.material.T_rigid < 2.0
    assert stats["skin_cells"] == b.skins.count.sum() and b.last.Q_conv > 0.0


def test_hand_down_passes_rigid_bottom_metal_beyond_the_minimum():
    mat = material.Material.from_drama_json("AA7075_scheil")
    sk = skin.Skins(mat, np.full(3, 4e-6), DEV)
    m_min = skin.minimum_mass(mat, sk.area, DEV)
    sk.fill(1.5 * m_min, np.array([700.0, 700.0, 905.0]))            # two cold, one slurry (above T_rigid)
    sk.take_bottom(np.array([0.0, 0.6, 0.0]) * m_min)                 # the second below its minimum
    M0, E0 = sk.column_mass(), sk.H.sum(axis=1)
    handed, H = skin.hand_down(sk, m_min, mat.T_rigid)
    assert handed == pytest.approx(np.array([0.5 * m_min[0], 0.0, 0.0]), rel=1e-12)
    assert sk.column_mass() == pytest.approx(M0 - handed, rel=1e-13)
    assert sk.H.sum(axis=1) == pytest.approx(E0 - H, rel=1e-13)
    assert H[0] == pytest.approx(handed[0] * float(mat.enthalpy(700.0)), rel=1e-9)


def test_handing_down_moves_the_base_out_and_the_3d_model_takes_the_metal(coarse_sphere_mesh):
    """A cold skin above its minimum hands its bottom metal to the 3D model: the corners move out by its volume, the 3D
    model's volume grows by it, and mass and energy are exact."""
    b = make(coarse_sphere_mesh, T0=600.0)
    b.emissivity = 0.0
    extra = 0.5 * skin.minimum_mass(b.material, b.columns.area, DEV)
    b.skins.add_bottom(extra, float(b.material.enthalpy(600.0)))
    b.skins.regrid()
    M0, E0, V0 = b.mass(), b.energy(), b.mesh.volume()
    # every column's excess over its minimum, the added metal and what it held above its minimum from the start
    excess = np.maximum(b.skins.column_mass() - skin.minimum_mass(b.material, b.columns.area, DEV), 0.0).sum()
    b.advance(0.5, 0.5, heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0))
    assert b.handed_mass == pytest.approx(excess, rel=1e-9) and excess > extra.sum()
    assert b.mesh.volume() - V0 == pytest.approx(excess / b.material.rho, rel=0.02)
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert np.abs(b.solver.temperature() - 600.0).max() < 1e-6


# -- how far a hot base absorbs (spec 6.6) ------------------------------------------------------------------------------
def test_the_absorption_depth_is_the_excess_over_the_gradient_along_the_normal_and_is_capped(coarse_sphere_mesh, monkeypatch):
    """The metal a base absorbs is the depth, along the column, from the base node to where the 3D temperature falls to
    T_rigid. The field imposed is exactly linear, T = T_rigid + G (x_c - x) with G = 1000 K/m, hotter towards -x (the front,
    where `make`'s v_hat points) with its T_rigid contour at x_c, 0.3 of the radius behind the centre: the corners beyond
    the contour are below T_rigid, the front hot corners are hotter than the interior below them (the gradient along the
    outward normal is positive, G times the normal's -x component), and the rear hot corners (the field falls outward
    there) have a gradient below GRADIENT_FLOOR. The 3D gradient is exact for a linear field, so each depth is known
    in closed form."""
    G = 1000.0
    b = make(coarse_sphere_mesh)
    mat, corners, n = b.material, b.moving.corners, b.columns.normal
    P = b.mesh.points
    radial = P[corners] - P.mean(axis=0)
    assert np.all(np.einsum("ij,ij->i", n, radial) > 0.99 * np.linalg.norm(radial, axis=1))      # the normals point outward
    x_c = 0.3 * np.linalg.norm(radial, axis=1).mean()
    T = mat.T_rigid + G * (x_c - P[:, 0])
    b.solver.set_temperature(T)
    excess = T[corners] - mat.T_rigid
    g_n = n @ np.array([-G, 0.0, 0.0])                       # the gradient of T along the outward normal, in closed form
    hot = excess > 0.0
    usable = hot & (g_n >= skinbody.GRADIENT_FLOOR)
    raw = np.where(hot, excess / np.maximum(g_n, skinbody.GRADIENT_FLOOR), 0.0)       # m, before the cap
    assert hot.sum() > 100 and (~hot).sum() > 100 and (hot & ~usable).sum() > 100
    # the sign: at a front hot corner the surface is hotter than the interior one millimetre below it, and that is a
    # positive gradient along the outward normal -- with the normal or the gradient reversed the depths below fail
    front = usable & (n[:, 0] < -0.5)
    inner = mat.T_rigid + G * (x_c - (P[corners] - 1e-3 * n)[:, 0])
    assert front.sum() > 100 and np.all(T[corners][front] > inner[front]) and np.all(g_n[front] > 0.0)

    # a whole part's allowance, a quarter of the mean surface edge (about 1.0 mm of 4.0 mm): every hot corner needs more
    cap = skinbody.MOVE_MAX * b.columns.edge_length.mean()
    assert np.all(raw[hot] > cap)
    depth = b._absorb_depth(cap)
    assert np.all(depth[~hot] == 0.0)                                      # below T_rigid: nothing to absorb
    assert depth[hot] == pytest.approx(cap, rel=1e-12) and b.absorb_capped == hot.sum()
    # per corner: what a corner's draws already used of the allowance is not absorbed again
    used = np.where(np.arange(len(raw)) % 2 == 0, 0.6 * cap, 0.0)
    b.absorb_capped = 0
    assert b._absorb_depth(cap - used, count=False)[hot] == pytest.approx((cap - used)[hot], rel=1e-12)
    assert b.absorb_capped == 0                                            # not counted when asked not to

    # an allowance too large to bind for the corners with a real gradient: their depth is the excess over the gradient
    cap = 0.4                                                              # m
    b.absorb_capped = 0
    depth = b._absorb_depth(cap)
    below = hot & (raw < cap)
    capped = hot & ~below
    assert below.sum() > 100 and capped.sum() > 100 and (capped & ~usable).sum() > 100     # the last: gradient under the floor
    assert np.all(depth[~hot] == 0.0)
    assert depth[below] == pytest.approx(raw[below], rel=1e-9) and np.all(depth[below] > 0.0)
    assert depth[capped] == pytest.approx(cap, rel=1e-12)                  # a gradient under the floor takes its allowance
    assert b.absorb_capped == capped.sum()


# -- what triggers the interface repeat (Task 7, measured in Task 10) --------------------------------------------------
def test_every_repeat_rule_keeps_the_books_and_records_what_triggered_it(coarse_sphere_mesh):
    """Each undamped repeat rule (measured in Task 10; the damped rule is the default) leaves the books exact, and the
    first pass's mismatch is recorded whatever the rule. 2 MW/m2 from 880 K on 4 mm elements: without a repeat, 3 MW/m2
    already runs away -- the 3D base overshoots T_feed and every sub-step removes and redraws a whole skin (measured
    while planning)."""
    results = {}
    for rule in ("column", "total", "never"):
        b = make(coarse_sphere_mesh, repeat_rule=rule, part_max=None)      # Task 10's conditions: whole 0.5 s parts
        M0 = b.mass()
        for k in range(3):
            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 2e6))
        assert b.removed_mass > 0.0
        assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
        assert abs(b.energy_balance_residual()) < 1e-8
        stats = b.skin_stats()
        assert stats["mismatch_total_ratio"] > 0.0 and stats["mismatch_worst_ratio"] >= stats["mismatch_total_ratio"]
        results[rule] = b.interface_repeats
    assert results["never"] == 0 and results["column"] >= results["total"] > 0
    with pytest.raises(ValueError):
        make(coarse_sphere_mesh, repeat_rule="sometimes")
    assert make(coarse_sphere_mesh).repeat_rule == skinbody.DEFAULT_REPEAT_RULE == "damped"     # Asha, 2026-10-12


def test_damped_repetition_converges_where_the_undamped_repeat_parks_the_heat(coarse_sphere_mesh):
    """3 MW/m2 on 4 mm elements from 880 K, one 0.5 s step: the skins' base heat falls off a cliff a little above the
    start-of-step base temperature, so the per-column rule's one repeat leaves heat in the deferred mismatch (7.7 % of
    the step's absorbed heat, measured) while the damped rule repeats with Aitken's relaxation until the mismatch is
    within INTERFACE_TOLERANCE (9 passes, 0.14 % left, measured). Mass is exact under both; the energy residual is below
    1e-8 under the damped rule (the column rule's, 1.2e-8 after this one step, is the parked heat's book-keeping, not
    asserted here)."""
    parked = {}
    for rule in ("column", "damped"):
        b = make(coarse_sphere_mesh, repeat_rule=rule, part_max=None)      # a whole 0.5 s part, where the cliff bites
        M0 = b.mass()
        b.advance(0.5, 0.5, front_loads(b, 3e6))
        assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
        parked[rule] = np.abs(b.pending).sum() / b.absorbed_heat
    assert abs(b.energy_balance_residual()) < 1e-8
    assert parked["column"] > 0.05
    assert parked["damped"] < 0.005
    assert 2 < b.last_passes < skinbody.DAMPED_MAX_PASSES and b.last_total <= skin.INTERFACE_TOLERANCE


# -- the move limit (2026-10-12) ---------------------------------------------------------------------------------------
def test_a_step_that_would_move_a_base_too_far_runs_as_parts_and_matches_the_shorter_steps(coarse_sphere_mesh, monkeypatch):
    """6 MW/m2 from 880 K on 4 mm elements: a 0.5 s step would draw and absorb more than MOVE_MAX of the mean edge
    (1 mm), so it is run again as shorter parts, every move's inward part within the limit (a mesh refused moves of
    1.65 mm here without it, measured; the tangential smoothing slides corners up to 1.4 mm along the surface besides);
    the parts give what the same shorter steps give when taken one by one with the same heating."""
    b = make(coarse_sphere_mesh, part_max=None)                          # the move limit alone, not the melting cap
    loads = front_loads(b, 6e6)
    moves = []
    policy = b.moving.move_with_policy

    def recorded(u, h_inward=None):
        moves.append(float(np.abs(np.einsum("ij,ij->i", u, b.columns.normal)).max()))
        return policy(u, h_inward)
    b.moving.move_with_policy = recorded
    M0 = b.mass()
    b.advance(0.5, 0.5, loads)
    n = b.step_parts
    limit = skinbody.MOVE_MAX * b.columns.edge_length.mean()
    assert n > 1 and len(moves) == n and max(moves) < 1.01 * limit       # the smoothing's own normal correction on top
    assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8
    monkeypatch.setattr(skinbody, "MOVE_HINT_MARGIN", 1.0)              # one step at a time, as the parts were taken
    c = make(coarse_sphere_mesh, part_max=None)
    for k in range(n):
        c.advance(0.5 * (k + 1) / n, 0.5 / n, loads)
    assert c.removed_mass == pytest.approx(b.removed_mass, rel=1e-6)
    assert c.energy() == pytest.approx(b.energy(), rel=1e-9)
    assert np.abs(c.moving.mesh.points[c.moving.corners] - b.moving.mesh.points[b.moving.corners]).max() < 1e-9


def test_a_base_fault_is_a_bottom_cell_outside_the_temperatures_of_the_metal_that_crossed():
    """The fault counter (spec 6.6, redefined 2026-10-12): a bottom cell within the tolerance of the range between the
    node's temperature and the 3D profile's at the depth swept is no fault, however far from the node; outside it, or
    above the node where nothing was drawn, it is; a column nothing crossed never counts."""
    tol = skin.JUMP_TOLERANCE
    T_node = np.full(6, 895.0)
    T_deep = np.array([884.0, 884.0, 884.0, 895.0, 895.0, 884.0])        # 0.5 mm drawn at 22 K/mm; nothing drawn
    T_bottom = np.array([886.0, 884.0 - tol + 0.1, 884.0 - tol - 0.1, 895.0 + tol + 0.1, 895.0 - tol + 0.1, 870.0])
    crossed = np.array([True, True, True, True, True, False])
    assert skinbody.jump_faults(T_bottom, T_node, T_deep, crossed).tolist() == [False, False, True, True, False, False]


# -- the part length while melting and the heating in time (2026-10-12) ------------------------------------------------
def test_while_melting_a_step_runs_as_short_parts_and_each_part_takes_the_heating_at_its_middle(coarse_sphere_mesh):
    """No part is longer than PART_MAX while some skin's top is non-rigid, none is capped before (from 880 K the first
    0.5 s step is one part; with the tops melting the second is 8, measured), and every pass of a part takes the
    heating interpolated linearly between the last macro step's and this one's at the part's middle, so a macro step
    receives the trapezoid of the two."""
    b = make(coarse_sphere_mesh)
    seen = []
    advance = skin.advance

    def recorded(sk, dt, q_top, *a, **k):
        seen.append((dt, np.array(q_top, dtype=float)))
        return advance(sk, dt, q_top, *a, **k)
    parts = []
    part = b._advance_part

    def recorded_part(t, dt, dt_macro, loads, depth):
        parts.append((t, dt))
        return part(t, dt, dt_macro, loads, depth)
    b._advance_part = recorded_part
    skinbody.skin_mod.advance = recorded
    try:
        q1 = front_loads(b, 2e6)
        b.advance(0.5, 0.5, q1)
        assert b.step_parts == 1 and all(np.array_equal(q, q1.q_conv) for _, q in seen)   # no earlier step: held
        assert np.any(b.skins.T[:, 0] >= b.material.T_rigid)
        seen.clear(), parts.clear()
        q2 = front_loads(b, 3e6)
        b.advance(1.0, 0.5, q2)
    finally:
        skinbody.skin_mod.advance = advance
    assert b.step_parts >= int(np.ceil(0.5 / skinbody.PART_MAX)) and max(dt for _, dt in parts) <= skinbody.PART_MAX + 1e-12
    expected = {round(dt, 12): [] for _, dt in parts}
    for t, dt in parts:
        frac = ((t - 0.5 * dt) - 0.5) / 0.5
        expected[round(dt, 12)].append(q1.q_conv + (q2.q_conv - q1.q_conv) * frac)
    for dt, q in seen:                                                     # every pass of every part
        assert any(np.allclose(q, e, rtol=1e-12, atol=0.0) for e in expected[round(dt, 12)])


def test_the_step_in_which_a_top_becomes_non_rigid_is_capped_too(coarse_sphere_mesh):
    """The part cap's look-ahead (Asha, 2026-10-16): a macro step whose tops are all below T_rigid at its start still runs
    in parts of at most PART_MAX when a top's rise over the last step, carried on over this one, would reach T_rigid;
    a smaller rise leaves the step whole."""
    mat = material.Material.from_drama_json("AA7075_scheil")
    for rise, capped in ((10.0, True), (4.0, False)):
        b = make(coarse_sphere_mesh)
        no_heat = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
        b.advance(0.5, 0.5, no_heat)                                     # the tops' start temperature recorded: 880 K
        sk = b.skins
        sk.T[:, 0] = 880.0 + rise                                        # as if the last step had raised them
        sk.H[:, 0] = sk.m[:, 0] * mat.enthalpy(sk.T[:, 0])
        b._top_start = np.full(sk.n, 880.0)
        b.advance(1.0, 0.5, no_heat)
        assert (b.step_parts == 8) is capped and b.capped_ahead == int(capped)   # 0.5 s / 0.0625 s
        assert np.all(sk.T[sk.valid] < mat.T_rigid)


# -- the reference tier: strong melting on the production mesh, and level 2b (2026-10-14) -----------------------------
@pytest.mark.reference
def test_strong_melting_on_the_production_mesh():
    """6 MW/m2 on the front of the 100 mm production sphere (2 mm surface, dense band, the skin runs' fine band of
    0.7 mm to 4 mm) for 1.5 s from 880 K: absorption takes the non-rigid metal as it forms, the base stays within 2 K of
    T_rigid, no skin cell is left buried above T_feed, and the books are exact. The 4 mm test sphere is refused at the
    rim of the melting region here; the 2 mm mesh follows it, moving no corner more than its quarter-edge allowance
    (0.49 mm). Measured 2026-10-14 without the fine band: 172.9 g removed, the base 1.2 K above T_rigid at most, one
    rebuild."""
    m = movemesh.refined_mesh(mesh.sphere_mesh(0.05, band=mesh.DEFAULT_BAND), 4e-3, 0.7e-3)
    b = make(m)
    M0 = b.mass()
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, front_loads(b, 6e6))
    assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8
    assert b.solver.temperature()[b.moving.corners].max() - b.material.T_rigid < 2.0
    assert not np.any(b.skins.valid & (b.skins.T > b.material.T_feed))
    assert 0.16 < b.removed_mass < 0.19                                    # kg; the heat could melt at most about 0.19


@pytest.mark.reference
def test_the_hot_plate_passes_level_2b():
    """The level-2 plate started at 880 K under 2 MW/m2 (`analysis/skin_plate.py --case hot`) at the flight's 0.5 s macro
    step, the skins at the solid's density as the reference holds its metal: the centre's recession within 5 % of the
    one-dimensional reference at 1.5 and 2.0 s, the skins' flowing depth within 10 % of its exact non-rigid depth
    (1.49-1.67 mm) and the centre's base within 1 K of T_rigid: the skins hold the non-rigid layer, none is left in the
    3D model."""
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
    import skin_plate
    T_start, q0 = skin_plate.CASES["hot"]
    body, rows, n_centre, failure = skin_plate.run(dt=0.5, t_end=2.0, verbose=False, T_start=T_start, q0=q0, swell=False)
    assert failure is None
    ref = skin_plate.hot_reference()
    at = lambda t_, key: float(np.interp(t_, [x["t"] for x in rows], [x[key] for x in rows]))
    for t_ in (1.5, 2.0):
        assert at(t_, "recession_mm") == pytest.approx(ref[t_][0], rel=0.05)
    for t_ in (1.0, 1.5, 2.0):
        assert at(t_, "flowing_mm") == pytest.approx(ref[t_][1], rel=skin_plate.HOT_FLOWING_TOL)
        assert at(t_, "base_excess_K") <= skin_plate.HOT_BASE_TOL_K
    assert max(abs(x["mass_error"]) for x in rows) < 1e-12 and max(abs(x["energy_residual"]) for x in rows) < 1e-8



def test_the_mass_is_measured_again_after_every_move(coarse_sphere_mesh):
    """Asha, 2026-10-17: the body measures its 3D volume once per move, not at every call -- the trajectory asks at
    every evaluation, and on the 100 mm sphere's fine band the measurements took 12 % of a flight's first steps. After
    steps that move the corners the mass is still exactly the 3D model's volume and the skins' metal, as measured
    afresh."""
    b = make(coarse_sphere_mesh)
    fresh = lambda: b.material.rho * b.mesh.volume() + b.skins.mass()
    assert b.mass() == fresh()
    V0 = b.mesh.volume()
    for k in range(2):
        b.advance(0.5 * (k + 1), 0.5, front_loads(b, 3e6))
        assert b.mass() == fresh()
    assert b.mesh.volume() != V0                                     # the corners moved


# -- the flow in the skin body (sub-plan 21 Task 6c) -------------------------------------------------------------------
def newtonian(b, p_stag=11e3, tau=50.0, decel=30.0):
    """A prescribed gas-side forcing on the body's corners: a Newtonian pressure falling away from the stagnation point,
    `tau` Pa of gas shear on the windward face and the deceleration along the flight direction."""
    cos = np.clip(b.columns.normal @ b.v_hat, 0.0, 1.0)
    return skinbody.FlowForcing(-2.0 * p_stag * np.sqrt(1.0 - cos ** 2) * cos / 0.05, np.where(cos > 0.0, tau, 0.0),
                                decel * b.v_hat)


def same_body(a, b):
    """Whether two skin bodies hold the same state, to the last bit."""
    return (np.array_equal(a.solver.temperature(), b.solver.temperature()) and np.array_equal(a.mesh.points, b.mesh.points)
            and all(np.array_equal(x, y) for x, y in zip(a.skins.state(), b.skins.state()))
            and np.array_equal(a.pending, b.pending))


def test_the_flow_is_off_by_default_and_turned_on_where_nothing_flows_it_changes_nothing(coarse_sphere_mesh):
    """The flow is off by default (the ruling of 2026-10-22: until stage 4's spraying gives the liquid a sink), and the
    heat-only skins are then unchanged: their tests above pass as they did, and the 3 MW/m2 removal run of
    `test_heating_with_removal_keeps_mass_and_energy_exact` reproduces Task 6b's code to the last bit, every array of
    the 3D model and the skins (measured 2026-10-22, numpy's generator seeded as a run seeds it: pyamg draws its
    starting vectors from it, so two unseeded runs differ in round-off). Turned on in a body where nothing flows -- 1
    MW/m2 on a body at 600 K for two 0.5 s steps -- the flow is driven and its books kept but no cell moves, and the
    body is the flow-off body to the last bit."""
    bodies = []
    for flow in (False, True):
        np.random.seed(12345)
        b = make(coarse_sphere_mesh, T0=600.0, flow=flow, remove=not flow)
        for k in range(2):
            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 1e6), forcing=newtonian(b))
        bodies.append(b)
    off, on = bodies
    assert off.flow_on is False and skinbody.FLOW_DEFAULT is False and on.flow_on is True
    assert make(coarse_sphere_mesh, T0=600.0).flow_on is False
    assert same_body(off, on) and np.all(on.skins.T[on.skins.valid] < on.material.T_rigid)
    a, c = off.skin_stats(), on.skin_stats()
    assert all(a[k] == c[k] for k in a if not k.startswith("wall_"))
    assert c["speed_max_ms"] == 0.0 and c["momentum_residual"] == 0.0 and c["flow_courant_max"] == 0.0


def test_a_skin_body_with_flow_keeps_its_mass_and_energy_closes_its_momentum_and_reports_the_flow(coarse_sphere_mesh):
    """3 MW/m2 on the front of the 4 mm test sphere from 880 K for two 0.25 s steps, flow on with a Newtonian pressure
    gradient, 50 Pa of gas shear and 30 m/s2 of deceleration, the liquid kept, one flow sub-step a heat sub-step (the
    books, not the flow's accuracy, are tested; the default's ten cost ten times the flow): the mass is the starting
    mass to 1e-13, the energy books close to 1e-8 of the heat absorbed, and the momentum balance -- the skins' momentum
    against the impulses of the push, the shear and the walls, and what left with metal handed down -- closes to
    round-off. The history's flow columns are the state's: the largest speed in any cell, the corner-to-corner roughness
    of the top, the flow sub-steps' largest Courant number, Rayleigh-Taylor growth and linearisations. Measured
    2026-10-22: mass to 1.1e-16, energy 9.8e-13, momentum 2.8e-17; the fastest cell 1.6 m/s, the flowing layer up to
    1.38 mm deep, Courant numbers to 6.9, growth 2.1e-3, roughness 0.14."""
    b = make(coarse_sphere_mesh, flow=True, remove=False, flow_cycles=1)
    M0 = b.mass()
    for k in range(2):
        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), forcing=newtonian(b))
    stats = b.skin_stats()
    assert b.removed_mass == 0.0 and b.mass() == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8
    assert stats["momentum_residual"] < 1e-12 and b.impulse_size > 0.0
    sk = b.skins
    speed = np.linalg.norm(sk.U, axis=2)
    assert stats["speed_max_ms"] == float(speed.max()) > 0.0
    i, j = b.columns.edges[:, 0], b.columns.edges[:, 1]
    H = sk.thickness()
    assert stats["roughness_max"] == float((np.abs(H[i] - H[j]) / b.columns.edge_length).max()) > 0.0
    assert stats["flow_courant_max"] > 0.0 and 0.0 <= stats["flow_growth_max"] < 1.0 and stats["flow_passes_max"] >= 1
    assert stats["flow_cycles_max"] == 1
    assert stats["wall_flow_s"] > 0.0 and stats["flowing_depth_max_mm"] > 0.0


def test_the_flow_is_driven_by_the_surface_flow_on_the_corners_taken_at_each_part_s_middle(coarse_sphere_mesh, monkeypatch):
    """Every macro step the skin body evaluates the surface flow on its corners, as the heating is, from the step's aero
    state and the tops' temperatures: the pressure gradient along the meridian, the gas shear along it on the top cells,
    and the deceleration, which the metal feels along the flight direction. Each part takes them interpolated linearly
    between the last macro step's and this one's at its middle, as it takes the heating."""
    pytest.importorskip("cantera")
    from test_reentry_model_surface_flow import aero_state_at
    from reentry_model import surface_flow
    b = make(coarse_sphere_mesh, T0=600.0, flow=True, remove=False)
    seen = []
    advance = skin.advance

    def recorded(sk, dt, *a, **k):
        seen.append((b.columns.normal.copy(), k["flow"]))              # the part's normals: the corners move after it
        return advance(sk, dt, *a, **k)
    monkeypatch.setattr(skinbody.skin_mod, "advance", recorded)
    evaluate = lambda a: surface_flow.SurfaceFlow().evaluate(a, b.theta, b.nose_radius(), b.material.liquid.rho,
                                                             b.surface_temperature())
    states = (aero_state_at(90.0), aero_state_at(92.0))
    first = evaluate(states[0])
    b.advance(0.5, 0.5, front_loads(b, 1e6), state=states[0])
    assert all(np.array_equal(f.dp_ds, first.dp_ds) for _, f in seen) and len(seen) >= 2
    normal, f = seen[0]
    np.testing.assert_allclose(f.shear, skin.along_meridian(normal, b.v_hat, first.tau), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(f.g, first.deceleration * b.v_hat, rtol=1e-12)
    assert f.rho_l == b.skins.rho_l and f.sigma == b.material.liquid.sigma and f.cycles is None    # chosen (Task 6d)
    second = evaluate(states[1])
    seen.clear()
    b._hint = 2                                                          # two parts, to see the interpolation
    b.advance(1.0, 0.5, front_loads(b, 1e6), state=states[1])
    assert b.step_parts == 2
    for frac in (0.25, 0.75):                                           # the two parts' middles
        dp = first.dp_ds + (second.dp_ds - first.dp_ds) * frac
        assert any(np.allclose(f.dp_ds, dp, rtol=1e-12, atol=1e-9) for _, f in seen)


def test_a_skin_body_s_flow_takes_the_counts_each_trial_pass_chose_in_its_real_pass(coarse_sphere_mesh, monkeypatch):
    """Sub-plan 21 Task 6d: with the flow on, the flow sub-steps in each heat sub-step are chosen from the last one's
    Courant number and growth (the default), and the real pass of every part takes the counts its trial pass chose, so
    that the base law is made by the pass it predicts; the history reports the most a heat sub-step took. 3 MW/m2 on the
    front of the 4 mm test sphere from 880 K for a 0.25 s step, as above: the books close as with a fixed count."""
    calls = []
    advance, substeps = skin.advance, skin.flow_substeps

    def recorded_advance(sk, dt, *a, **k):
        calls.append(("trial" if k.get("T_b1") is None else "real", []))
        return advance(sk, dt, *a, **k)

    def recorded_substeps(sk, flow, dt, cycles=None):
        book = substeps(sk, flow, dt, cycles)
        calls[-1][1].append((cycles, book.cycles))                    # the count handed down, and the one used
        return book
    monkeypatch.setattr(skin, "advance", recorded_advance)
    monkeypatch.setattr(skin, "flow_substeps", recorded_substeps)
    b = make(coarse_sphere_mesh, flow=True, remove=False)
    assert b.flow_cycles is None
    M0 = b.mass()
    b.advance(0.25, 0.25, front_loads(b, 3e6), forcing=newtonian(b))
    kinds = [k for k, _ in calls]
    assert kinds == ["trial", "real"] * (len(calls) // 2) and len(calls) >= 2
    chosen = [[c for c, _ in cs] for _, cs in calls]
    assert all(chosen[i] == chosen[i + 1] for i in range(0, len(calls), 2))
    used = [u for _, cs in calls for _, u in cs if u]                 # where something flowed
    assert used and max(used) <= skin.FLOW_CYCLES_MAX and min(used) >= skin.FLOW_CYCLES_MIN
    stats = b.skin_stats()
    assert max(u for _, u in calls[-1][1]) <= stats["flow_cycles_max"] <= max(used) and stats["flow_courant_max"] > 0.0
    assert b.mass() == pytest.approx(M0, rel=1e-13) and abs(b.energy_balance_residual()) < 1e-8
    assert stats["momentum_residual"] < 1e-12


# -- the droplet momentum sink (sub-plan 21 Task 7) ---------------------------------------------------------------------
def test_a_skin_body_with_a_prescribed_release_keeps_its_books_and_reports_the_droplet_momentum_ratio(coarse_sphere_mesh):
    """Sub-plan 21 Task 7, the placeholder for stage 4's spraying: a release of 2 kg/(m2 s) from every column's top, liquid
    first, as sub-plan 22's Girin rate will be given, on the 4 mm test sphere heated at 3 MW/m2 from 880 K with the flow on
    (as above: two 0.25 s steps, one flow sub-step a heat sub-step). The tops melt in the second step and the release
    takes the liquid as it forms, 0.53 g. It leaves the body: the mass left and the mass released are the starting mass
    to the last bit, the energy books close with the released enthalpy to 9.8e-13 of the heat absorbed, and the
    momentum balance, the droplets' momentum lost with the metal, to 2.3e-17. The history reports the released mass,
    fact 29's droplet-momentum ratio Pi over the macro step, which is its definition on the step's books (0.0048 here,
    with the tops barely liquid and slow), and the columns flagged above 0.2, none (measured 2026-10-24). A release
    needs the flow: with it off the heat-only skins remove the liquid as it forms."""
    with pytest.raises(ValueError):
        make(coarse_sphere_mesh, release=skin.prescribed_release(2.0))
    b = make(coarse_sphere_mesh, flow=True, remove=False, flow_cycles=1, release=skin.prescribed_release(2.0))
    M0 = b.mass()
    for k in range(2):
        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), forcing=newtonian(b))
    stats = b.skin_stats()
    assert b.removed_mass == 0.0 and b.released_mass > 0.0 and stats["released_mass_kg"] == b.released_mass
    assert b.mass() + b.released_mass == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8 and b.released_enthalpy > 0.0
    assert stats["momentum_residual"] < 1e-12 and np.abs(b.droplet_momentum).max() > 0.0
    Pi, _, flagged = skin.droplet_momentum_ratio(b.droplet_momentum, b.droplet_shear)
    assert stats["droplet_momentum_ratio"] == Pi > 0.0 and stats["droplet_momentum_flagged"] == flagged
    off = make(coarse_sphere_mesh, flow=True, remove=False, flow_cycles=1)
    assert off.skin_stats()["droplet_momentum_ratio"] == 0.0 and off.skin_stats()["released_mass_kg"] == 0.0


# -- level 2 with the flow on: the plate (sub-plan 21 Task 8) -----------------------------------------------------------
def skin_plate_module():
    import os
    import sys
    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis")
    if here not in sys.path:
        sys.path.insert(0, here)
    import skin_plate
    return skin_plate


# a reduced plate for the unit tier: 24 mm across, heated within 7 mm of its centre, the flux gone at 10 mm; 2 mm
# elements and no fine band (box_mesh is structured: no meshing), the skins' development cells
SMALL_PLATE = dict(verbose=False, width=24e-3, r_flat=7e-3, r_zero=10e-3, fine=None)


def test_the_plate_with_the_flow_off_is_sub_plan_20_s_plate():
    """`analysis/skin_plate.py` with the flow off builds and steps the plate as sub-plan 20's script did, so its level-2
    results stand: a reduced plate from 880 K under 6 MW/m2 for two 0.1 s steps, which melts, removes liquid and hands
    metal down, is the body built and advanced by hand to the last bit (numpy's generator seeded: pyamg draws from it).
    Sub-plan 20's script itself, run on the planning copy's package, gives Task 8's every array and history column to the
    last bit on such a plate over three steps (measured 2026-10-24, measure/task8/bitid_flowoff.py). A push, a release or
    a flow sub-step count with the flow off is refused."""
    sp_ = skin_plate_module()
    np.random.seed(12345)
    body, rows, _, failure = sp_.run(dt=0.1, t_end=0.2, T_start=880.0, q0=6e6, **SMALL_PLATE)
    np.random.seed(12345)
    mat = material.Material.from_drama_json(material.MATERIAL_NAMES["AA7075_scheil"])
    m = mesh.box_mesh(sp_.LENGTH, 24e-3, 24e-3, 2e-3)
    m.points[:, 1:] -= 12e-3
    hand = skinbody.SkinBody(m, mat, thermal_solver("skfem"), mat.rho * m.volume(),
                             skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"], swell=False), T0=880.0, emissivity=0.0,
                             v_hat=(-1.0, 0.0, 0.0), deform_depth=4e-3, repeat_rule=skinbody.DEFAULT_REPEAT_RULE)
    p = hand._top0
    radius = np.where(np.abs(p[:, 0]) < 0.05e-3, np.hypot(p[:, 1], p[:, 2]), np.inf)
    for k in range(2):
        q = sp_.projected_flux(hand, radius, 7e-3, 10e-3, 6e6)
        hand.advance(0.1 * (k + 1), 0.1, heating.HeatingResult(q, 6e6, 0.0, 0.0, 0.0))
    assert failure is None and body.removed_mass > 0.0 and same_body(body, hand)
    for kw in ({"push": sp_.PUSHES["flight"]}, {"release": 1.0}, {"flow_cycles": 2}):
        with pytest.raises(ValueError):
            sp_.run(t_end=0.0, **SMALL_PLATE, **kw)


def test_a_plate_pushed_along_its_face_closes_its_books_and_its_metal_runs_the_way_the_push_points():
    """Spec 11.1's level 2 with the flow on: the reduced plate from 880 K under 6 MW/m2 for two 0.1 s steps, its face a
    patch of the 100 mm sphere's shoulder 45 degrees from the stagnation point at 93 s (`skin_plate.PUSHES['flight']`:
    dp/ds -3.42e5 Pa/m and 80 Pa of shear down the face's meridian, -y, and 29.8 m/s2 of deceleration along the tilted
    flight direction), the flow sub-steps chosen, the liquid kept. The books close to round-off -- the mass to 1e-13, the
    energy to 1e-8 of the heat absorbed, the momentum to 1e-12 of the impulses -- and the metal runs down the face: on the
    face the columns downstream of the heat (y < 0) gain more than a hundred times what those upstream do, and the face's
    volume flux points within 10 degrees of -y. Measured 2026-10-24 (measure/task8/measure_tests.py): mass 2.2e-16,
    energy 1.3e-12, momentum 3.1e-17 (4.7e-17 in a second run: pyamg's starting vectors are not seeded here); gained
    downstream 3.3e-5 kg, upstream 7.2e-8 kg (464 times); the face's flux 5.0 degrees from -y (on the plate's crater
    wall the face's columns tilt), the fastest cell 2.4 m/s, 10 and 4 flow sub-steps a heat sub-step at most. On this
    plate the runoff reaches the edge 2 mm beyond the heat within the first step and runs over it: with the side's
    columns the whole flux is 19 degrees from -y."""
    sp_ = skin_plate_module()
    body, rows, _, failure, records = sp_.run(dt=0.1, t_end=0.2, T_start=880.0, q0=6e6, flow=True,
                                              push=sp_.PUSHES["flight"], record=True, **SMALL_PLATE)
    assert failure is None and body.flow_on and body.removed_mass == 0.0 and body.released_mass == 0.0
    assert max(abs(x["mass_error"]) for x in rows) < 1e-13
    assert max(abs(x["energy_residual"]) for x in rows) < 1e-8
    assert max(x["momentum_residual"] for x in rows) < 1e-12 and body.impulse_size > 0.0
    face = sp_.face_columns(body)
    y = body._top0[:, 1]
    gained = np.maximum(records[-1]["mass"] - records[0]["mass"], 0.0)
    assert gained[face & (y < 0.0)].sum() > 100.0 * gained[face & (y > 0.0)].sum() > 0.0
    flux = records[-1]["flux"][face].sum(axis=0)
    assert -flux[1] / np.linalg.norm(flux) > np.cos(np.radians(10.0))
    assert rows[-1]["speed_max_ms"] > 0.0 and rows[-1]["flow_cycles_max"] >= 1


def test_a_plate_with_the_flow_on_and_nothing_flowing_is_the_flow_off_plate():
    """With no push and no metal above T_rigid -- the reduced plate from 700 K under 1 MW/m2 for two 0.1 s steps -- the
    flow-on plate is the flow-off plate to the last bit (numpy's generator seeded): no cell flows, so no flow sub-step
    is taken, nothing moves and nothing is released. Once the top turns to slurry the two part, even with no push and
    nothing liquid: surface tension acts where a top flows (6b's ruling) and moves the slurry wherever the tops' heights
    differ -- on the plate where the draws and hand-downs have left neighbouring columns a little apart -- so that
    starting at 880 K under 2 MW/m2, 45 tops of slurry and none liquid, the flow-on plate's skins differ from the
    flow-off plate's by 6.8e-7 of their largest thickness after 0.2 s and its 3D field by 6.3e-4 K, at speeds of at most
    1.6e-5 m/s (measured 2026-10-24, measure/task8/measure_tests.py)."""
    sp_ = skin_plate_module()
    bodies = []
    for flow in (False, True):
        np.random.seed(12345)
        bodies.append(sp_.run(dt=0.1, t_end=0.2, T_start=700.0, q0=1e6, flow=flow, **SMALL_PLATE)[0])
    off, on = bodies
    assert on.flow_on and same_body(off, on) and np.all(on.skins.T[on.skins.valid] < on.material.T_rigid)
    assert on.skin_stats()["speed_max_ms"] == 0.0 and on.momentum_residual() == 0.0
