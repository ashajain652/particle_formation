"""The skin body (layered-skin design sections 6.6 and 7, sub-plan 20 Task 7): skins on every corner from time zero, the
trial and real passes around the 3D step, the base moving the 3D mesh, and the books."""
import numpy as np
import pytest

from reentry_model import heating, material, mesh, skin, skinbody
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
    # each column holds what its corner swept moving in by its minimum: within 2 % of it (prisms against dual areas
    # on 4 mm triangles), so the first sub-step draws any shortfall
    assert np.abs(b.skins.column_mass() / m_min - 1.0).max() < 0.02
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


# -- the reference tier: strong melting on the production mesh, and level 2b (2026-10-14) -----------------------------
@pytest.mark.reference
def test_strong_melting_on_the_production_mesh():
    """6 MW/m2 on the front of the 100 mm production sphere (2 mm surface, dense band) for 1.5 s from 880 K: absorption
    takes the non-rigid metal as it forms, the base stays within 2 K of T_rigid, no skin cell is left buried above
    T_feed, and the books are exact. The 4 mm test sphere is refused at the rim of the melting region here; the 2 mm
    mesh follows it, moving no corner more than its quarter-edge allowance (0.49 mm). Measured 2026-10-14: 172.9 g
    removed, the base 1.2 K above T_rigid at most, one rebuild."""
    m = mesh.sphere_mesh(0.05, band=mesh.DEFAULT_BAND)
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
    step: the centre's recession within 5 % of the one-dimensional reference at 1.5 and 2.0 s, the skins' flowing depth
    within 10 % of its exact non-rigid depth (1.49-1.67 mm) and the centre's base within 1 K of T_rigid: the skins hold
    the non-rigid layer, none is left in the 3D model."""
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
    import skin_plate
    T_start, q0 = skin_plate.CASES["hot"]
    body, rows, n_centre, failure = skin_plate.run(dt=0.5, t_end=2.0, verbose=False, T_start=T_start, q0=q0)
    assert failure is None
    ref = skin_plate.hot_reference()
    at = lambda t_, key: float(np.interp(t_, [x["t"] for x in rows], [x[key] for x in rows]))
    for t_ in (1.5, 2.0):
        assert at(t_, "recession_mm") == pytest.approx(ref[t_][0], rel=0.05)
    for t_ in (1.0, 1.5, 2.0):
        assert at(t_, "flowing_mm") == pytest.approx(ref[t_][1], rel=skin_plate.HOT_FLOWING_TOL)
        assert at(t_, "base_excess_K") <= skin_plate.HOT_BASE_TOL_K
    assert max(abs(x["mass_error"]) for x in rows) < 1e-12 and max(abs(x["energy_residual"]) for x in rows) < 1e-8

