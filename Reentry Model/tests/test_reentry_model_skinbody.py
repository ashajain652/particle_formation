"""The skin body (layered-skin design sections 6.6 and 7, sub-plan 20 Task 7): skins on every corner from time zero, the
trial and real passes around the 3D step, the base moving the 3D mesh, and the books."""
import numpy as np
import pytest

from reentry_model import coupled, heating, material, mesh, movemesh, skin, skinbody
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


def test_hand_down_passes_only_bottom_metal_within_the_fault_tolerance_of_the_base_node():
    """Spec 6.6, amended 2026-10-26 (sub-plan 21 Task 9): rigid metal at a skin's bottom passes to the 3D model only while
    no cell handed is warmer than the base node by more than the fault tolerance (5 K, `T_limit` the node's temperature
    plus it); warmer metal stays in the skin until it has cooled. Three skins half their minimum above it on nodes at
    700 K: at the node's temperature throughout (the whole excess goes, as without a limit), rising 50 K from bottom to
    top (only the bottom cells within 705 K go, whole), and with a bottom cell 10 K warm (nothing goes). With `span` the
    coldest and warmest of the cells it took from come back, for the fault counter (nan where none went). `handable`
    gives what a hand-down would take, taking nothing: with no limit, less what went, it is the metal the tolerance holds
    in the skin (the skin body's held mass)."""
    mat = material.Material.from_drama_json("AA7075_scheil")
    sk = skin.Skins(mat, np.full(3, 4e-6), DEV)
    m_min = skin.minimum_mass(mat, sk.area, DEV)
    sk.fill(1.5 * m_min, 700.0)
    k = sk.count
    T = np.array([np.full(sk.m.shape[1], 700.0), np.linspace(750.0, 700.0, k[1]), np.linspace(760.0, 710.0, k[2])])
    sk._install(sk.m, sk.m * mat.enthalpy(T), k)
    m0, H0 = sk.m.copy(), sk.H.copy()
    able, limited = skin.handable(sk, m_min, mat.T_rigid), skin.handable(sk, m_min, mat.T_rigid, np.full(3, 705.0))
    assert np.array_equal(sk.m, m0) and np.array_equal(sk.H, H0)            # it takes nothing
    free = skin.hand_down(sk.rows(np.arange(3)), m_min, mat.T_rigid)[0]
    handed, H, lo, hi = skin.hand_down(sk, m_min, mat.T_rigid, T_limit=np.full(3, 705.0), span=True)
    np.testing.assert_array_equal(able, free)
    np.testing.assert_array_equal(limited, handed)
    assert free == pytest.approx(0.5 * m_min, rel=1e-12) and handed[0] == pytest.approx(free[0], rel=1e-12)
    assert lo[0] == pytest.approx(700.0, abs=1e-9) and hi[0] == pytest.approx(700.0, abs=1e-9)
    went = int(np.sum(T[1] <= 705.0))                                  # the bottom cells within 5 K of the node
    assert 0 < went < k[1] and handed[1] == pytest.approx(sk.m[1, 0] * went, rel=1e-12) and handed[1] < free[1]
    assert lo[1] == pytest.approx(700.0, abs=1e-9) and 700.0 < hi[1] <= 705.0 and sk.T[1, sk.count[1] - 1] > 705.0
    assert handed[2] == 0.0 and H[2] == 0.0 and np.isnan(lo[2]) and np.isnan(hi[2])


def test_hand_down_takes_at_most_its_cap_a_cell_going_in_part():
    """Spec 6.6 and 7.1, amended 2026-10-27 (Asha, option 1 of Task 9's third hand-back): a part's hand-down at a column
    takes at most `cap` [kg], the metal that moves its base by the move allowance, a cell going in part where the cap
    ends inside it; the rest waits for later parts. On the full pushed plate at 1.5 s the lowest held cells were 0.18 mm
    of solid metal (median) and up to 0.79 mm against the 0.556 mm allowance, each released whole the moment it came
    within the fault tolerance of its node, so a release did not shrink with the part: the part was split six times in
    vain and kept over the allowance, and the 12.4 g held came down in 0.039 s (measure/task9, item 4). Three skins at
    their minimum plus 0.8 mm, graded from 704 K at the top to 700 K at the bottom, all within a 705 K limit: capped
    inside a cell (0.31 of the minimum), inside another (0.555), and not at all (twice the minimum). Mass and enthalpy
    are conserved, the cell taken in part keeps its temperature and stays the bottom cell, the coldest and warmest cells
    taken include it, and a cap that does not bind takes to the last bit what no cap takes."""
    mat = material.Material.from_drama_json("AA7075_scheil")
    sk = skin.Skins(mat, np.full(3, 4e-6), DEV)
    m_min = skin.minimum_mass(mat, sk.area, DEV)
    sk.fill(2.0 * m_min, 700.0)
    k = sk.count
    T = np.array([np.linspace(704.0, 700.0, k[i]) for i in range(3)])
    sk._install(sk.m, sk.m * mat.enthalpy(T), k)
    limit, cap = np.full(3, 705.0), np.array([0.31, 0.555, 2.0]) * m_min
    free, free_H = skin.hand_down(sk.rows(np.arange(3)), m_min, mat.T_rigid, limit)
    unlimited = sk.rows(np.arange(3))
    cell = sk.m[0, 0]                                                    # fill makes equal cells
    M0, E0 = sk.column_mass(), sk.H.sum(axis=1)
    handed, H, lo, hi = skin.hand_down(sk, m_min, mat.T_rigid, limit, span=True, cap=cap)
    assert free == pytest.approx(m_min, rel=1e-12)                      # every cell within the limit: the whole excess
    np.testing.assert_allclose(handed[:2], cap[:2], rtol=1e-12)
    assert handed[2] == free[2] and H[2] == free_H[2]                    # not binding: the same to the last bit
    assert sk.column_mass() == pytest.approx(M0 - handed, rel=1e-13) and sk.H.sum(axis=1) == pytest.approx(E0 - H, rel=1e-13)
    for c in (0, 1):                                                     # the cell the cap ended in, left in part
        whole = int(handed[c] // cell)
        rest = sk.m[c, sk.count[c] - 1]
        assert 0.0 < rest < cell and rest == pytest.approx((whole + 1) * cell - handed[c], rel=1e-9)
        assert sk.T[c, sk.count[c] - 1] == pytest.approx(T[c, k[c] - 1 - whole], abs=1e-9)
        assert lo[c] == pytest.approx(700.0, abs=1e-9) and hi[c] == pytest.approx(T[c, k[c] - 1 - whole], abs=1e-9)
    big = skin.hand_down(unlimited, m_min, mat.T_rigid, limit, cap=np.full(3, 1.0))
    np.testing.assert_array_equal(big[0], free)
    np.testing.assert_array_equal(big[1], free_H)


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


def deposit_on(body, cols, T_top_old, depth=0.5e-3, T_dep=890.0):
    """The columns `cols` given a rigid deposit `depth` [m of solid] at T_dep on top of their skins, whose own metal is
    warmed linearly from the base node's temperature at its bottom to T_top_old [K, per column] at its top -- the
    deposit's heat conducting down, as where runoff froze on Task 8's pushed plate -- the cells then split and merged."""
    sk, mat = body.skins, body.material
    n, K = sk.m.shape
    nd = int(round(depth / DEV.cell))
    T_node = body.solver.temperature()[body.moving.corners]
    m, H, count = np.zeros((n, K + nd)), np.zeros((n, K + nd)), sk.count.copy()
    m[:, :K], H[:, :K] = sk.m, sk.H
    for c, T_top in zip(cols, T_top_old):
        k = sk.count[c]
        m[c, :nd + k] = np.r_[np.full(nd, mat.rho * sk.area[c] * depth / nd), sk.m[c, :k]]
        H[c, :nd + k] = m[c, :nd + k] * mat.enthalpy(np.r_[np.full(nd, T_dep), np.linspace(T_top, T_node[c], k)])
        count[c] = nd + k
    sk._install(m, H, count)
    sk.regrid()


def test_a_warm_deposit_goes_down_only_near_the_node_s_temperature_at_its_own_enthalpy_and_heats_no_cell(
        coarse_sphere_mesh):
    """Sub-plan 21 Task 9's cause and its fix (spec 6.6, amended 2026-10-26). Thirty front columns of the 4 mm sphere at
    700 K hold a 0.5 mm deposit at 890 K over their 0.8 mm skin, which the deposit has warmed linearly to 880 K at its
    top (fifteen columns; 225 K/mm, as steep as the deposit's heat front gets) or to 706 K (fifteen; 7.5 K/mm), and the
    base's operations run (`SkinBody._move_base`, no pass before them). Before Task 9 the whole 0.5 mm went down, the
    moving mesh booked it at the base node's enthalpy -- 420 429 J/kg against the steep slab's own 537 279 J/kg -- and
    the reconciliation put the difference into the bottom cell the hand-down left: 1 614 K in a column whose hottest
    metal was 890 K (measured 2026-10-25, measure/task9/test_task9_cause.py). Now only cells within the fault tolerance
    of the node go down: none of the steep skins' -- the regrid has merged their small bottom cells into graded ones
    85.6 K above the node -- and the gentle skins' whole excess, 2.8 K above it, which the 3D model takes at its own
    423 923 J/kg rather than the node's 420 429 J/kg; the steep skins' deposit waits above their minimum, no cell ends
    hotter than the hottest metal of its column before, and mass and energy are exact (measured 2026-10-26). The body
    records what waits (Asha, 2026-10-26: the hold is unbounded, and sub-plan 22's flow-on flight reports it): the held
    mass, the steep skins' whole excess, and its deepest column's depth."""
    b = make(coarse_sphere_mesh, T0=700.0)
    mat, sk = b.material, b.skins
    cols = np.argsort(b.theta)[:30]
    steep, gentle = cols[:15], cols[15:]
    deposit_on(b, cols, np.r_[np.full(15, 880.0), np.full(15, 706.0)])
    m_min = skin.minimum_mass(mat, sk.area, DEV)
    before = sk.rows(np.arange(sk.n))                                  # a copy, for what the hand-down should take
    T_hot = np.where(sk.valid, sk.T, 0.0).max(axis=1)
    T_node = b.solver.temperature()[b.moving.corners]
    M0, E0, area = b.mass(), b.energy(), b.columns.area.copy()          # the areas the base move measures depths on
    b._move_base(skin.SkinBook(*(np.zeros(sk.n) for _ in range(8)), momentum_out=np.zeros((sk.n, 3))), m_min,
                 b._normal_gradient())
    T_after = np.where(sk.valid, sk.T, 0.0).max(axis=1)
    worst = int(np.argmax(T_after - T_hot))
    assert np.all(T_after <= T_hot + 0.5), "{:.1f} K in a column whose hottest metal was {:.1f} K".format(
        T_after[worst], T_hot[worst])
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    excess = before.column_mass() - m_min
    handed, handed_H, lo, hi = skin.hand_down(before, m_min, mat.T_rigid, T_node + skin.JUMP_TOLERANCE, span=True)
    assert np.all(handed[steep] == 0.0) and np.all(handed[gentle] > 0.0)
    assert np.all(hi[gentle] <= T_node[gentle] + skin.JUMP_TOLERANCE) and np.all(sk.column_mass()[steep] > 1.5 * m_min[steep])
    np.testing.assert_allclose(handed[gentle], excess[gentle], rtol=1e-12)
    stats = b.skin_stats()
    assert b.held_mass == pytest.approx(excess[steep].sum(), rel=1e-12) and stats["held_mass_kg"] == b.held_mass
    assert stats["held_depth_max_mm"] == pytest.approx((excess[steep] / (mat.rho * area[steep])).max() * 1e3, rel=1e-12)
    r = b.last_move
    np.testing.assert_allclose(r.swept_enthalpy[gentle] / (mat.rho * r.swept_volume[gentle]),
                               handed_H[gentle] / handed[gentle], rtol=1e-12)


def test_a_deposit_thicker_than_the_move_allowance_comes_down_over_several_parts_and_splits_none(coarse_sphere_mesh):
    """Spec 6.6 and 7.1, amended 2026-10-27 (Asha, option 1 of Task 9's third hand-back): a part's hand-down moves no
    base further than the move allowance, a cell going in part, so a hand-down never splits a part; draws still can.
    Thirty front columns of the 4 mm sphere at 700 K are given 2.5 move allowances of metal at the node's temperature
    at their bottoms (the allowance, a quarter of the mean edge, is about 1 mm here) and three 0.5 s macro steps run
    with no heat. Under the rule before, the first part asked for 2.51 allowances at every length from 0.5 s to 0.69
    ms, was split six times in vain and kept the last with a move of 2.51 allowances (measured 2026-10-27; on the full
    pushed plate 8 of the last 15 parts, 1.015-1.19 allowances, measure/task9 item 4). Now no part asks for the
    allowance and none is split, every base moves out by
    at most the allowance (the smoothing's own normal correction on top), the metal comes down an allowance a part
    while the body's held mass records what waits -- 1.5 allowances after the first part, nothing once it is down --
    no fault is counted where a cell went down in part, and mass and energy are exact."""
    b = make(coarse_sphere_mesh, T0=700.0)
    b.emissivity = 0.0
    mat, sk, col = b.material, b.skins, b.columns
    cols = np.argsort(b.theta)[:30]
    allow, area0 = skinbody.MOVE_MAX * col.edge_length.mean(), col.area.copy()   # the areas grow as the bases move out
    extra = np.zeros(sk.n)
    extra[cols] = 2.5 * allow * mat.rho * area0[cols]
    sk.add_bottom(extra, float(mat.enthalpy(700.0)))
    sk.regrid()
    ratios, outward, held = [], [], []
    ratio, policy = b._move_ratio, b.moving.move_with_policy

    def recorded_ratio(book, m_min):
        ratios.append(ratio(book, m_min))
        return ratios[-1]

    def recorded_move(u, h_inward=None, h_outward=None):
        outward.append(float(np.einsum("ij,ij->i", u, b.columns.normal).max()))
        return policy(u, h_inward, h_outward)
    b._move_ratio, b.moving.move_with_policy = recorded_ratio, recorded_move
    M0, E0 = b.mass(), b.energy()
    zero = heating.HeatingResult(np.zeros(col.n), 0.0, 0.0, 0.0, 0.0)
    parts = 0
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, zero)
        parts += b.step_parts
        held.append(b.held_mass)
    assert max(ratios) < 1.0 and len(ratios) == parts                    # none asked for more, none was split
    assert max(outward) < 1.01 * allow
    assert b.handed_mass > extra.sum()                                    # all of it down, and the skins' own excess
    assert held[0] == pytest.approx(1.5 * allow * mat.rho * area0[cols].sum(), rel=0.02)   # 1.51 measured (2026-10-27)
    assert held[1] < 1e-3 * held[0] and held[2] < 1e-3 * held[0]
    assert b.jump_faults == 0
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert np.abs(b.solver.temperature() - 700.0).max() < 1e-6


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

    def recorded(u, h_inward=None, h_outward=None):
        moves.append(float(np.abs(np.einsum("ij,ij->i", u, b.columns.normal)).max()))
        return policy(u, h_inward, h_outward)
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


def test_a_hand_down_s_bottom_cell_is_judged_against_the_handed_cells_and_the_cell_it_uncovered():
    """The fault counter as amended 2026-10-26, and after the build (spec 6.6; sub-plan 21 Task 9): the bottom cell, read
    after the reconciliation, is judged against the metal that crossed -- where metal was drawn or absorbed, the range
    between the node's temperature and the 3D profile's at the depth swept, as before; where metal was handed down, the
    range from the handed cells' coldest and warmest temperatures to that of the cell the hand-down uncovered, as the
    hand-down left it; a column that did both may lie within either; a column nothing crossed never counts. So the skin's
    own warmer interior that a hand-down uncovers is no fault (column 2: 28 K above the warmest cell handed), while a
    bottom cell the reconciliation then heats beyond that range still is. Judged by the handed cells alone, the counter
    counted the uncovered cell, 1 900 of the 1 942 faults on Task 8's pushed plate to 1.5 s; by the node alone, before
    2026-10-26, it counted it too and never saw the superheating."""
    tol = skin.JUMP_TOLERANCE
    T_node = np.full(6, 800.0)
    T_deep = np.array([790.0, 800.0, 800.0, 800.0, 790.0, 800.0])
    inward = np.array([True, False, False, False, True, False])        # drawn or absorbed
    handed = np.array([False, True, True, True, True, False])
    lo = np.array([np.nan, 800.0, 800.0, 800.0, 820.0, np.nan])        # the handed cells' coldest and warmest
    hi = np.array([np.nan, 812.0, 812.0, 812.0, 825.0, np.nan])
    uncovered = np.array([np.nan, 812.0, 840.0, 805.0, 826.0, np.nan])  # the bottom cell the hand-down left
    T_bottom = np.array([791.0, 812.0 + tol - 0.1, 840.0, 800.0 - tol - 0.1, 822.0, 2000.0])
    faults = skinbody.jump_faults(T_bottom, T_node, T_deep, inward, handed=(handed, lo, hi, uncovered))
    assert faults.tolist() == [False, False, False, True, False, False]
    heated = T_bottom.copy()
    heated[2] = 840.0 + tol + 0.1                                      # the reconciliation heated it beyond the range
    assert skinbody.jump_faults(heated, T_node, T_deep, inward, handed=(handed, lo, hi, uncovered))[2]
    assert skinbody.jump_faults(T_bottom, T_node, T_deep, inward, handed=(handed, lo, hi)).tolist() == [
        False, False, True, True, False, False]                       # the handed cells alone
    assert skinbody.jump_faults(T_bottom, T_node, T_deep, inward).tolist() == [False, False, False, False, True, False]


def graded_column(bottom_excess, dz=(0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3), T_base=895.0, slope=8000.0):
    """One column of graded cells (top first): a linear profile rising from `T_base` at the base by `slope` K/m, the bottom
    cell `bottom_excess` K above the line at its centre (Task 5k's measured case: a bottom cell about 2 mm thick and 8 K
    warmer at its centre than the base, on a profile rising 8 K/mm)."""
    dz = np.array(dz)[None, :]
    height = np.cumsum(dz[:, ::-1], axis=1)[:, ::-1] - 0.5 * dz            # each centre's height above the base [m]
    T = T_base + slope * height
    T[0, -1] += bottom_excess
    return T, dz, np.array([dz.shape[1]])


def test_a_thick_bottom_cell_on_a_smooth_profile_is_not_a_jump_fault():
    """Task 5k, Asha's option (a), 2026-11-10: a bottom cell 2 mm thick, 8 K warmer at its centre than the base on a linear
    profile (so the base itself is at the node's temperature) is no fault under the new rule, as the 50 mm flight's 336
    faults a step were not; the bottom cell's rule, kept as `jump_faults_bottom_cell`, counts it (8 K against 5 K). Fails on the
    old rule (its `[False]` assertion): the one of the four that does; the other three guard the "only removes faults" property."""
    T, dz, count = graded_column(0.0)
    T_base = skinbody.base_temperature(T, dz, count)
    assert T_base[0] == pytest.approx(895.0, abs=1e-9) and T[0, -1] - 895.0 == pytest.approx(8.0)
    node, crossed = np.array([895.0]), np.array([True])
    args = (T[:, -1], node, node, crossed)
    assert skinbody.jump_faults(*args).tolist() == [True]                          # today's rule: the bottom cell alone
    assert skinbody.jump_faults(*args, T_base=T_base).tolist() == [False]


def test_a_genuine_jump_in_a_thick_bottom_cell_is_still_a_fault():
    """The same column with the bottom cell 8 K above the line through the cells above it (16 K above the base node): the
    extrapolation to the base does not hide it (it lands 12.6 K above the node, above the 5 K tolerance). A guard: passes on the old
    rule too, and fails on today's code only because `base_temperature` is missing."""
    T, dz, count = graded_column(8.0)
    T_base = skinbody.base_temperature(T, dz, count)
    assert T_base[0] > 895.0 + skin.JUMP_TOLERANCE
    node, crossed = np.array([895.0]), np.array([True])
    assert skinbody.jump_faults(T[:, -1], node, node, crossed, T_base=T_base).tolist() == [True]


def test_fresh_metal_at_the_node_under_a_steep_step_is_still_not_a_fault():
    """A thin bottom cell (0.1 mm) at the node's temperature below a cell 205 K warmer: the extrapolation would put the base
    18.6 K below the node, outside the range, but the bottom cell is in it, so the column passes on its own temperature: the
    change can only remove faults the profile explains, never add one. A guard: passes on the old rule too."""
    T, dz = np.array([[1100.0, 895.0]]), np.array([[1.0e-3, 0.1e-3]])
    T_base = skinbody.base_temperature(T, dz, np.array([2]))
    node, crossed = np.array([895.0]), np.array([True])
    assert T_base[0] < 895.0 - skin.JUMP_TOLERANCE                                 # the extrapolation alone would fault it
    assert skinbody.jump_faults(T[:, -1], node, node, crossed, T_base=T_base).tolist() == [False]
    assert skinbody.jump_faults(T[:, -1], node, node, crossed).tolist() == [False]


def test_a_one_cell_column_is_judged_by_its_own_temperature():
    """With one cell there is no profile to extrapolate: the base temperature is the cell's, so the two rules agree, in both
    directions (897 K passes, 903 K faults, against a node at 895 K and a tolerance of 5 K). A guard: passes on the old rule too."""
    T, dz = np.array([[897.0, 0.0], [903.0, 0.0]]), np.array([[2.0e-3, 0.0], [2.0e-3, 0.0]])
    count = np.array([1, 1])
    T_base = skinbody.base_temperature(T, dz, count)
    assert T_base.tolist() == [897.0, 903.0]
    node, crossed = np.full(2, 895.0), np.full(2, True)
    assert skinbody.jump_faults(T[:, 0], node, node, crossed, T_base=T_base).tolist() == [False, True]
    assert skinbody.jump_faults(T[:, 0], node, node, crossed).tolist() == [False, True]


def test_the_history_and_the_results_carry_both_fault_counts(coarse_sphere_mesh):
    """`jump_faults` (the new rule) and `jump_faults_bottom_cell` (today's, kept so that earlier flights stay comparable) are
    both skin columns, next to each other, and in the skins' stats; the new rule never counts a column the old does not."""
    assert "jump_faults_bottom_cell" in coupled.SKINBODY_COLUMNS
    i = coupled.SKINBODY_COLUMNS.index("jump_faults")
    assert coupled.SKINBODY_COLUMNS[i + 1] == "jump_faults_bottom_cell"
    b = make(coarse_sphere_mesh, T0=890.0)
    s = b.skin_stats()
    assert s["jump_faults"] == 0 and s["jump_faults_bottom_cell"] == 0
    b.jump_faults_bottom_cell = 3
    assert b.skin_stats()["jump_faults_bottom_cell"] == 3 and b.skin_stats()["jump_faults"] <= 3


def stepped_deposit(body, cols, cold=4, T_warm=760.0, depth=0.5e-3, T_dep=890.0):
    """The columns `cols` given a rigid deposit `depth` [m of solid] at T_dep over their skins, whose bottom `cold` cells
    sit at the base node's temperature and the rest at T_warm -- a skin steeply warmer above its base -- the cells left
    as they are (no regrid), so that a hand-down stops at the first warm cell."""
    sk, mat = body.skins, body.material
    n, K = sk.m.shape
    nd = int(round(depth / DEV.cell))
    T_node = body.solver.temperature()[body.moving.corners]
    m, H, count = np.zeros((n, K + nd)), np.zeros((n, K + nd)), sk.count.copy()
    m[:, :K], H[:, :K] = sk.m, sk.H
    for c in cols:
        k = sk.count[c]
        T_old = np.r_[np.full(k - cold, T_warm), np.full(cold, T_node[c])]
        m[c, :nd + k] = np.r_[np.full(nd, mat.rho * sk.area[c] * depth / nd), sk.m[c, :k]]
        H[c, :nd + k] = m[c, :nd + k] * mat.enthalpy(np.r_[np.full(nd, T_dep), T_old])
        count[c] = nd + k
    sk._install(m, H, count)


def test_the_cell_a_hand_down_uncovers_is_no_fault_unless_the_reconciliation_heats_it(coarse_sphere_mesh, monkeypatch):
    """Spec 6.6's fault as amended after the build (2026-10-26): thirty front columns of the 4 mm sphere at 700 K hold a
    0.5 mm deposit over a skin whose bottom four cells are at the node's 700 K and the rest at 760 K; the base's
    operations (`SkinBody._move_base`) hand the four down and uncover a cell 60 K warmer than any handed, which the
    reconciliation leaves as it is: no fault, where judged by the handed cells alone all thirty were. With the move's
    sweep made to return 0.05 J more to each of those columns' bottoms -- a reconciliation heating the cell it uncovered --
    all thirty are faults."""
    counts = {}
    for heat in (False, True):
        b = make(coarse_sphere_mesh, T0=700.0)
        cols = np.argsort(b.theta)[:30]
        stepped_deposit(b, cols)
        sk = b.skins
        if heat:
            policy = b.moving.move_with_policy

            def heated(u, h_inward=None, h_outward=None):
                r = policy(u, h_inward, h_outward)
                r.swept_enthalpy[cols] += 0.05
                return r
            monkeypatch.setattr(b.moving, "move_with_policy", heated)
        m_min = skin.minimum_mass(b.material, sk.area, DEV)
        b._move_base(skin.SkinBook(*(np.zeros(sk.n) for _ in range(8)), momentum_out=np.zeros((sk.n, 3))), m_min,
                     b._normal_gradient())
        counts[heat] = b.jump_faults
        if not heat:
            bottom = sk.T[cols, sk.count[cols] - 1]
            assert np.abs(bottom - 760.0).max() < 0.5 and b.handed_mass > 0.0
    assert counts == {False: 0, True: 30}


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
    same = lambda x, y: x == y or (np.isnan(x) and np.isnan(y))     # a column empty in both (sub-plan 22 Task 4's)
    assert all(same(a[k], c[k]) for k in a if not k.startswith("wall_"))
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


def test_a_skin_body_s_real_pass_takes_at_least_its_trial_pass_s_counts_and_its_own_where_the_trial_flowed_nothing(
        coarse_sphere_mesh, monkeypatch):
    """Sub-plan 21 Task 6d, amended by Task 9 (spec 16, 2026-10-26): with the flow on, the flow sub-steps in each heat
    sub-step are chosen from the last one's Courant number and growth (the default), and the real pass of every part
    takes, in each heat sub-step, the larger of the count its trial pass chose and the count its own memory chooses,
    and its own where the trial pass took no flow sub-step -- so that a real pass whose base melts more than its trial
    pass's runs at the Courant target rather than above it (18.8 against 6 on Task 8's pushed plate), and does not replay
    counts chosen where nothing flowed (380 empty flow sub-steps there) -- while wherever the two passes flow alike the
    base law is made by the pass it predicts (6d). The history reports the most a heat sub-step took. 3 MW/m2 on the
    front of the 4 mm test sphere from 880 K for a 0.25 s step, as above: the books close as with a fixed count."""
    calls = []
    advance, substeps, choose = skin.advance, skin.flow_substeps, skin.flow_cycles

    def recorded_advance(sk, dt, *a, **k):
        calls.append(("trial" if k.get("T_b1") is None else "real", [], []))
        return advance(sk, dt, *a, **k)

    def recorded_substeps(sk, flow, dt, cycles=None):
        book = substeps(sk, flow, dt, cycles)
        calls[-1][1].append((cycles, book.cycles))                    # the count handed down, and the one used
        return book

    def recorded_choice(memory, dt):
        count = choose(memory, dt)
        calls[-1][2].append(count)                                    # what each pass's own memory chose
        return count
    monkeypatch.setattr(skin, "advance", recorded_advance)
    monkeypatch.setattr(skin, "flow_substeps", recorded_substeps)
    monkeypatch.setattr(skin, "flow_cycles", recorded_choice)
    b = make(coarse_sphere_mesh, flow=True, remove=False)
    assert b.flow_cycles is None
    M0 = b.mass()
    b.advance(0.25, 0.25, front_loads(b, 3e6), forcing=newtonian(b))
    kinds = [k for k, _, _ in calls]
    assert kinds == ["trial", "real"] * (len(calls) // 2) and len(calls) >= 2
    for (_, trial, mine_t), (_, real, mine_r) in zip(calls[0::2], calls[1::2]):
        assert [c for c, _ in trial] == mine_t                        # the trial pass takes its own choice
        took = [u for _, u in trial]
        assert [c for c, _ in real] == [r if t == 0 else max(c, r) for t, c, r in zip(took, mine_t, mine_r)]
    used = [u for _, cs, _ in calls for _, u in cs if u]              # where something flowed
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
    face the columns downstream of the heat (y < 0) gain more than ten times what those upstream do, and the face's
    volume flux points within 10 degrees of -y. Measured 2026-10-24 (measure/task8/measure_tests.py): mass 2.2e-16,
    energy 1.3e-12, momentum 3.1e-17 (4.7e-17 in a second run: pyamg's starting vectors are not seeded here); gained
    downstream 3.3e-5 kg, upstream 7.2e-8 kg (464 times); the face's flux 5.0 degrees from -y (on the plate's crater
    wall the face's columns tilt), the fastest cell 2.4 m/s, 10 and 4 flow sub-steps a heat sub-step at most. On this
    plate the runoff reaches the edge 2 mm beyond the heat within the first step and runs over it: with the side's
    columns the whole flux is 19 degrees from -y. Since Task 9 (2026-10-26) no part is run again for its moves -- three
    parts, where 14 were tried and 10 kept -- and the skins' mass, which this compares, ends otherwise shared between the
    skins and the 3D model at the heat's upstream rim: 19 upstream columns, flowing 0.39-0.51 mm deep, end at 0.97-1.03
    of their minimum and 1.86e-6 kg above where they started, against 3.42e-5 kg gained downstream (18 times, against 464;
    the threshold is Task 8's ruling 6's ten), while the upstream skins lost 7.0e-6 kg in all; the flux 5.0 degrees from
    -y as before; books 2.2e-16, 3.1e-12, 1.1e-17."""
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
    assert gained[face & (y < 0.0)].sum() > 10.0 * gained[face & (y > 0.0)].sum()
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


def test_on_the_flow_off_plate_the_fault_tolerance_holds_no_hand_down_back(monkeypatch):
    """Of Task 9's changes to the base (spec 6.6, amended 2026-10-26) only the booking of handed-down metal at its own
    enthalpy moves a flow-off result; the limit that keeps metal warmer than the node by more than the fault tolerance in
    the skin holds nothing back where the skins' bottoms are at their nodes' temperature. On the reduced flow-off plate
    from 880 K under 6 MW/m2 for two 0.1 s steps, which melts, removes liquid and hands metal down, every hand-down takes
    to the last bit what it would with no limit. Measured 2026-10-26 (measure/task9/flowoff_state.py, TASK9_NO_A=1): with
    the limit alone the level-2 plates, cold and hot, are the same to 1.0 s, every array and history column; with the
    booking too this plate's cells move by at most 0.048 K and its 3D field by 0.0042 K, its masses and recession not at
    all, and the level-2 plates' centres by about 1e-8 mm at 1.0 s. Nor does the cap on a part's hand-down at the move
    allowance (amended 2026-10-27) bind there: every hand-down takes what it would with neither limit nor cap."""
    sp_ = skin_plate_module()
    held = []
    hand_down = skin.hand_down

    def compared(sk, m_min, T_rigid, T_limit=None, span=False, cap=None):
        free = hand_down(sk.rows(np.arange(sk.n)), m_min, T_rigid)[0]
        out = hand_down(sk, m_min, T_rigid, T_limit, span, cap=cap)
        held.append((T_limit is not None or cap is not None) and not np.array_equal(out[0], free))
        return out
    monkeypatch.setattr(skin, "hand_down", compared)
    body, rows, _, failure = sp_.run(dt=0.1, t_end=0.2, T_start=880.0, q0=6e6, **SMALL_PLATE)
    assert failure is None and body.handed_mass > 0.0 and held and not any(held)
    assert body.skin_stats()["held_mass_kg"] == 0.0 and body.skin_stats()["held_depth_max_mm"] == 0.0


# -- Girin's spraying from the top cells (sub-plan 22 Task 1c) ----------------------------------------------------------
def girin_flow(n, closure, u, tau, decel=106.0, branch=None):
    """The fields Girin's evaluation reads from the surface flow, per column: the spray tests' FakeFlow (edge state,
    closure, deceleration) with the gas shear and push the lubrication closure reads."""
    from reentry_model import surface_flow
    from test_reentry_model_spray import FakeFlow
    f = FakeFlow(n, branch=surface_flow.BRANCH_SHOCK_LAYER if branch is None else branch, decel=decel)
    f.closure, f.u_eff = np.asarray(closure), np.asarray(u, dtype=float)
    f.u_e = f.u_eff
    f.tau, f.G = np.asarray(tau, dtype=float), np.zeros(n)
    return f


def plan_on(sk, col, flow, slurry=False, theta=None):
    from reentry_model import spray
    from test_reentry_model_spray import FakeState
    return skinbody.girin_plan(sk, flow, FakeState(), spray.SprayModel(sk.material.liquid), np.zeros(sk.n) if theta is None
                               else theta, col.area, col.edges, 0.05, 0.0125, slurry=slurry)


def crest(lam, h):
    return (3.0 * lam ** 2 * h / (4.0 * np.pi)) ** (1.0 / 3.0)


def test_girin_gives_each_prescribed_column_its_mode_from_the_skin_s_depths_once_for_the_macro_step():
    """Spec 8 and 6.7 with Asha's decision of 2026-10-29 (option B): b and the regime test read the flowing depth, the
    viscous film's floor too; the Rayleigh-Taylor mode reads the liquid pool and the molten region the wave must fit in
    is the liquid's; off Girin's closure a liquid on slurry is held from the shear modes; the cap reads the liquid with
    slurry spraying off. Six columns 20 mm apart at the nose (one molten region about 30 mm across): 1 mm of liquid over
    rigid metal under Girin's closure (thick); 50 um of liquid on the rigid floor (thin); 5 mm of liquid with no shear
    mode at 5 m/s (Rayleigh-Taylor, lambda* 2 cm at 106 m/s2); 0.3 mm of liquid over 4 mm of slurry the same (held: read
    as one flowing film the Rayleigh-Taylor mode would take it); 50 um of liquid over slurry under the Couette closure
    (held); exposed slurry (nothing: no liquid to take). The plan's rates are the shared spray model's on those depths,
    and its release gives them whatever the skins hold until the next macro step."""
    from reentry_model import spray, surface_flow
    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
    from test_reentry_model_spray import FakeState
    col, _ = strip_grid(3, 2, h=20e-3)
    temps = [[LIQUID] * 86 + [RIGID] * 10, [LIQUID] * 4 + [RIGID] * 10, [LIQUID] * 428 + [RIGID] * 10,
             [LIQUID] * 26 + [SLURRY] * 340 + [RIGID] * 10, [LIQUID] * 4 + [SLURRY] * 50 + [RIGID] * 10,
             [SLURRY] * 50 + [RIGID] * 10]
    sk = stacked(col, temps)
    G, C = surface_flow.CLOSURE_GIRIN, surface_flow.CLOSURE_COUETTE
    flow = girin_flow(6, [G, G, C, C, C, G], [2000.0, 2000.0, 5.0, 5.0, 2000.0, 2000.0], [50.0, 50.0, 0.01, 0.01, 50.0, 50.0])
    plan = plan_on(sk, col, flow)
    mat, liq = sk.material, sk.material.liquid
    F, L = sk.depth_above(mat.T_rigid), sk.depth_above(mat.T_liquidus)
    assert plan.result.branch.tolist() == [spray.BRANCH_THICK, spray.BRANCH_THIN, spray.BRANCH_RT, -1, -1,
                                           spray.BRANCH_THICK]
    assert (plan.rate > 0.0).tolist() == [True, True, True, False, False, False] and not plan.rate_slurry.any()
    np.testing.assert_array_equal(plan.flowing, F)
    np.testing.assert_array_equal(plan.liquid, L)
    np.testing.assert_array_equal(plan.cap_depth, L)
    assert plan.on_slurry.tolist() == [False, False, False, True, True, True]
    # the shared model on the skin's depths, by hand
    delta_m, _ = spray.melt_layer(flow, liq)
    from reentry_model import film
    v_s = film.lubrication(flow.tau, flow.G, F, delta_m, liq.mu, b_layer=F)[0]
    extent = np.where(L > 0.0, 2.0 * np.sqrt(col.area[L > 0.0].sum() / np.pi), 0.0)     # one liquid region
    run = np.cumprod(sk.valid & (sk.T >= mat.T_liquidus), axis=1).astype(bool)
    res = spray.SprayModel(liq).evaluate(flow, FakeState(), F, delta_m, v_s, np.ones(6, bool), 0.0125, col.area,
                                         (sk.m * run).sum(axis=1), 0.05, b_layer=F, extent=extent,
                                         deceleration_n=np.full(6, 106.0), regime_layer=F, on_slurry=plan.on_slurry,
                                         rt_layer=L, cap_depth=L)
    np.testing.assert_allclose(plan.rate, res.mdot, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(plan.result.r, res.r, rtol=1e-12, atol=0.0)
    assert plan.rate[2] == pytest.approx(liq.rho * L[2] / res.growth[2], rel=1e-12)     # the pool shed in a growth time
    # read as one flowing film the fourth column would be the Rayleigh-Taylor mode's (option A, not taken)
    one = spray.SprayModel(liq).evaluate(flow, FakeState(), F, delta_m, v_s, np.ones(6, bool), 0.0125, col.area,
                                         (sk.m * run).sum(axis=1), 0.05, b_layer=F, extent=np.full(6, extent.max()),
                                         deceleration_n=np.full(6, 106.0), regime_layer=F, cap_depth=L)
    assert one.branch[3] == spray.BRANCH_RT and one.mdot[3] > 0.0
    # held for the macro step: the release gives the plan's rates whatever the skins then hold
    before = plan.release(sk)
    skin.release_top(sk, sk.column_mass())
    assert plan.release(sk) is before and np.array_equal(before.rate, plan.rate)
    assert np.array_equal(before.rate_slurry, plan.rate_slurry) and np.isinf(before.T_slurry).all()


def test_the_cap_reads_the_liquid_with_slurry_spraying_off_the_flowing_depth_with_it_on_and_never_the_area():
    """The decision of 2026-10-24: a droplet holds at most one crest of its own wave on the depth the release may take --
    the liquid with slurry spraying off, the whole flowing depth with it on -- and the column's area does not enter. Six
    columns of one liquid cell over 0.5 mm of slurry under Girin's closure, sheared near the thick branch's threshold so
    that its wavelength is millimetres: the liquid's 15 um binds, the flowing 0.6 mm does not; the corner columns of the
    strip have smaller areas than the others and the same droplet."""
    from reentry_model import surface_flow
    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
    col, _ = strip_grid(3, 2, h=20e-3)
    sk = stacked(col, [[LIQUID] + [SLURRY] * 50 + [RIGID] * 10] * 6)
    flow = girin_flow(6, [surface_flow.CLOSURE_GIRIN] * 6, [2000.0] * 6, [20.0] * 6)
    off, on = plan_on(sk, col, flow), plan_on(sk, col, flow, slurry=True)
    mat = sk.material
    L, F = sk.depth_above(mat.T_liquidus), sk.depth_above(mat.T_rigid)
    lam = off.result.wavelength
    assert np.isfinite(lam).all() and (off.rate > 0.0).all() and np.ptp(col.area) > 0.1 * col.area.max()
    np.testing.assert_allclose(off.result.r, np.minimum(0.17 * lam, crest(lam, L)), rtol=1e-12)
    np.testing.assert_allclose(on.result.r, np.minimum(0.17 * lam, crest(lam, F)), rtol=1e-12)
    assert (off.result.r < 0.17 * lam).all() and (on.result.r == 0.17 * lam).all()
    np.testing.assert_array_equal(on.cap_depth, F)
    assert np.ptp(off.result.r) <= 1e-12 * off.result.r.max() and np.ptp(on.result.r) <= 1e-12 * on.result.r.max()


def test_the_slurry_cascade_sprays_exposed_slurry_by_its_own_properties_only_with_the_flag():
    """Spec 8's cascade (the 2026-10-07 design's 8.1-8.2, Asha 2026-10-29, option i): with `--slurry-spray on` the slurry
    the release exposes sprays by Girin's model with its own properties -- Li et al.'s viscosity at its solid fraction and
    the shear rate of its own sheared layer, found by fixed point (five passes at most, the rate held at 367 1/s), the
    cell's mixture density and 0.80 N/m -- its layer the slurry from its top down to T_rigid; the Rayleigh-Taylor mode is
    the liquid's. Five columns of slurry under 15 um of liquid, the slurry's top cell at 902 K (30 % solid), 906.0 K, 907.0
    K, 907.5 K and 907.9 K (the last three in the bridge to the liquid). Off, nothing; on, the slurry's viscous film blocks
    its thin mode wherever the slurry is more than a few times the liquid's viscosity, and it sprays only near the
    liquidus."""
    from reentry_model import rheology, surface_flow
    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
    col, _ = strip_grid(5, 2, h=20e-3)
    tops = [902.0, 906.0, 907.0, 907.5, 907.9]
    temps = [[LIQUID] + [T] + [SLURRY] * 48 + [RIGID] * 10 for T in tops] * 2
    sk = stacked(col, temps)
    flow = girin_flow(col.n, [surface_flow.CLOSURE_GIRIN] * col.n, [2000.0] * col.n, [50.0] * col.n)
    off, on = plan_on(sk, col, flow), plan_on(sk, col, flow, slurry=True)
    assert not off.rate_slurry.any() and np.isnan(off.slurry["mu"]).all()
    mat = sk.material
    s = on.slurry
    f_s = 1.0 - sk.props.liquid_fraction(sk.T[:, 1])                    # the cells' own, as their viscosity's
    f_b = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])
    np.testing.assert_allclose(s["mu"], rheology.viscosity(sk.T[:, 1], f_s, s["gamma"], mat.liquid.mu, f_b), rtol=1e-12)
    np.testing.assert_allclose(s["rho"], sk.density(sk.T[:, 1]), rtol=1e-12)
    assert (s["gamma"] <= rheology.GAMMA_MAX).all() and (s["mu"] > mat.liquid.mu).all()
    dz = sk.dz()
    np.testing.assert_allclose(s["depth"], on.flowing - dz[:, 0], rtol=1e-12)          # the slurry, top to T_rigid
    blocked = s["gamma_tau"] >= 1.0
    assert (on.rate_slurry[blocked] == 0.0).all() and blocked[[0, 5]].all()             # 30 % solid: blocked
    assert (on.rate_slurry[~blocked] > 0.0).all() and (~blocked[[4, 9]]).all()          # 907.9 K: sprays
    assert (np.diff(s["mu"][:5]) < 0.0).all()                                           # runnier toward the liquidus
    np.testing.assert_array_equal(s["T"], sk.T[:, 1])
    np.testing.assert_array_equal(on.release(sk).T_slurry, sk.T[:, 1])     # its rate for slurry at least as warm
    np.testing.assert_array_equal(on.rate, off.rate)            # the liquid's own mode: the flag moves only its cap


def test_slurry_modes_composite():
    """The slurry's mode under the composite law: its viscosity is `composite_viscosity` of the exposed cell's solid
    fraction to 1e-12, whatever the shear rate, and the recorded rate is not clipped at 367 1/s (the bridge law's
    rate cap). The cascade test's five columns at a high gas shear."""
    import dataclasses
    from reentry_model import rheology, surface_flow
    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
    col, _ = strip_grid(5, 2, h=20e-3)
    tops = [902.0, 906.0, 907.0, 907.5, 907.9]
    sk = stacked(col, [[LIQUID] + [T] + [SLURRY] * 48 + [RIGID] * 10 for T in tops] * 2)
    sk.settings = dataclasses.replace(sk.settings, rheology="composite")
    flow = girin_flow(col.n, [surface_flow.CLOSURE_GIRIN] * col.n, [2000.0] * col.n, [500.0] * col.n)
    s = plan_on(sk, col, flow, slurry=True).slurry
    mat = sk.material
    f_s = 1.0 - sk.props.liquid_fraction(sk.T[:, 1])
    np.testing.assert_allclose(s["mu"], rheology.composite_viscosity(f_s, mat.liquid.mu), rtol=1e-12)
    assert (s["gamma"] > rheology.GAMMA_MAX).any()


def test_a_skin_body_sprays_by_girin_from_its_top_cells_once_a_macro_step_and_keeps_its_books(coarse_sphere_mesh,
                                                                                              monkeypatch):
    """The skin body with the flow on and a spray model (sub-plan 22 Task 1c): at the start of every macro step Girin is
    evaluated once on the corners, from the surface flow the flow's forcing evaluates there and the skins' depths then;
    every part's Flow releases the plan's two rates, the same in every call, through the cascade in every flow sub-step.
    The 4 mm test sphere heated at 3 MW/m2 from 890 K under the 100 mm flight's gas state at 93 s, two 0.25 s steps, one
    flow sub-step a heat sub-step: its tops melt in the first, and in the second the liquid is sprayed; the released
    metal leaves the body, the mass left and released is the starting mass to round-off, the energy books close with the
    released enthalpy and the momentum balance with the droplets' momentum. The plan's record counts what each column's
    mode released in the step's accepted passes, and the droplets: their mass over one droplet's at the plan's radius.
    A spray needs the flow and is not a prescribed release."""
    pytest.importorskip("cantera")
    from test_reentry_model_surface_flow import aero_state_at
    from reentry_model import spray
    mat = material.Material.from_drama_json("AA7075_scheil")
    model = spray.SprayModel(mat.liquid)
    with pytest.raises(ValueError):
        make(coarse_sphere_mesh, spray=model)
    with pytest.raises(ValueError):
        make(coarse_sphere_mesh, flow=True, remove=False, spray=model, release=skin.prescribed_release(1.0))
    b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=model)
    plans, flows = [], []
    girin, advance = skinbody.girin_plan, skin.advance

    def planned(*a, **k):
        plans.append(girin(*a, **k))
        return plans[-1]

    def recorded(sk, dt, *a, **k):
        flows.append(k["flow"])
        return advance(sk, dt, *a, **k)
    monkeypatch.setattr(skinbody, "girin_plan", planned)
    monkeypatch.setattr(skinbody.skin_mod, "advance", recorded)
    M0 = b.mass()
    state = aero_state_at(93.0)
    for k in range(2):
        flows.clear()
        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), state=state)
        assert len(plans) == k + 1 and b.spray_plan is plans[-1]
        assert all(f.release == plans[-1].release for f in flows) and len(flows) >= 2
    p = b.spray_plan
    assert b.released_mass > 0.0 and p.mass.sum() > 0.0 and not p.mass_slurry.any()
    assert b.mass() + b.released_mass == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8 and b.skin_stats()["momentum_residual"] < 1e-12
    ok = p.mass > 0.0
    assert (p.rate[ok] > 0.0).all() and np.isfinite(p.result.r[ok]).all()
    np.testing.assert_allclose(p.number[ok], p.mass[ok] / (4.0 / 3.0 * np.pi * mat.liquid.rho * p.result.r[ok] ** 3),
                               rtol=1e-12)
    assert b.droplet_number >= p.number.sum() > 0.0


# -- the droplet records (sub-plan 22 Task 2) ---------------------------------------------------------------------------
def test_a_plan_s_records_are_one_per_column_and_mode_with_the_released_cells_solid_fraction_enthalpy_and_velocity():
    """Spec 8 and 10.2: the droplet records are combined per column, per mode, per macro step, summing mass and number
    over the sub-steps; a column's radius is its plan's, so combining changes none; the liquid's mode and the slurry's
    are separate records; `solid_fraction` and `h_kg` are mass-weighted over the cells each record released, and the
    droplets leave at the mean velocity of the metal they are made of. Three columns of 10 um cells under Girin's
    closure, slurry spraying on: liquid alone (the thin branch), exposed slurry three cells of it in the bridge at 907.9 K
    (the slurry's thick mode), and one liquid cell over that slurry -- the liquid moving at 2 m/s along x, the slurry at
    0.5 m/s along y. Three flow sub-steps of 12 us: the third column gives its liquid cell in the first two and slurry
    after it. Four rows, the liquid's mode's first; the size histograms accumulated each sub-step are the records'."""
    from reentry_model import spray, surface_flow
    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
    from test_reentry_model_spray import FakeState
    col, pts = strip_grid(3, 2, h=20e-3)
    WARM = 907.9
    sk = stacked(col, [[LIQUID] * 5 + [RIGID] * 10, [WARM] * 3 + [SLURRY] * 47 + [RIGID] * 10,
                       [LIQUID] + [WARM] * 3 + [SLURRY] * 46 + [RIGID] * 10, [RIGID] * 10, [RIGID] * 10, [RIGID] * 10])
    mat = sk.material
    liquid = sk.valid & (sk.T >= mat.T_liquidus)
    sk.U = np.where(liquid[..., None], [2.0, 0.0, 0.0], [0.0, 0.5, 0.0]) * sk.valid[..., None]
    flow = girin_flow(6, [surface_flow.CLOSURE_GIRIN] * 6, [2000.0] * 6, [50.0] * 6)

    class State(FakeState):
        h = 71e3
    theta = np.radians([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    plan = skinbody.girin_plan(sk, flow, State(), spray.SprayModel(mat.liquid), theta, col.area, col.edges, 0.05, 0.0125,
                               slurry=True, points=pts)
    assert plan.result.branch[:3].tolist() == [spray.BRANCH_THIN, spray.BRANCH_THICK, spray.BRANCH_THICK]
    assert (plan.rate > 0.0).tolist() == [True, False, True] + [False] * 3
    assert (plan.rate_slurry > 0.0).tolist() == [False, True, True] + [False] * 3
    hist = np.zeros((2, spray.N_BINS))
    taken, m0 = np.zeros((2, 6)), sk.m[:, 0].copy()
    for _ in range(3):
        cascade = plan.release(sk)
        dm, _ = skin.cascade_mass(sk, cascade.rate, cascade.rate_slurry, col.area, 12e-6, True, cascade.T_slurry)
        ph = skin.release_top(sk, dm, slurry=True).by_phase
        fb = skin.FlowBook.zero(6)
        fb.phase_mass, fb.phase_H, fb.phase_p, fb.phase_solid = ph.mass, ph.enthalpy, ph.momentum, ph.solid
        plan.book(fb)
        hist += cascade.histogram(ph.mass)
        taken += ph.mass
    assert taken[0, 2] == pytest.approx(m0[2], rel=1e-12) and 0.0 < taken[1, 2] < m0[2]        # all its liquid
    plan.count(mat.liquid.rho)
    rows = plan.source_block(0.0125)
    C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
    i, mode = np.array([0, 2, 1, 2]), np.array([0, 0, 1, 1])
    assert rows.shape == (4, 27) and np.isnan(rows[:, C["release_angle_register_deg"]]).all()    # no centroid given
    np.testing.assert_array_equal(rows[:, C["x_m"]:C["z_m"] + 1], pts[i])
    np.testing.assert_array_equal(rows[:, C["theta_deg"]], np.degrees(theta[i]))
    assert (rows[:, C["time_s"]] == 0.0125).all() and (rows[:, C["altitude_km"]] == 71.0).all()
    s = plan.slurry
    r = np.where(mode == 0, plan.result.r[i], s["r"][i])
    np.testing.assert_array_equal(rows[:, C["r_m"]], r)                                  # combining changes no radius
    np.testing.assert_array_equal(rows[:, C["dm_kg"]], taken[mode, i])
    rho = np.where(mode == 0, mat.liquid.rho, s["rho"][i])
    np.testing.assert_allclose(rows[:, C["dn"]], taken[mode, i] / (4.0 / 3.0 * np.pi * rho * r ** 3), rtol=1e-14)
    for key, liq, sl in (("branch", plan.result.branch, s["branch"]), ("delta_m_m", plan.result.delta_m, s["delta_m"]),
                         ("we_s", plan.result.we_s, s["we_s"]), ("b_m", plan.flowing, s["depth"]),
                         ("h_floor_m", plan.result.h_floor, s["h_floor"]), ("gamma_tau", plan.result.gamma_tau,
                                                                             s["gamma_tau"])):
        np.testing.assert_array_equal(rows[:, C[key]], np.where(mode == 0, liq[i], sl[i]))
    # the two new columns: liquid at 950 K carries no solid and its enthalpy; the slurry its cells' at 907.9 K
    f_s = 1.0 - float(sk.props.liquid_fraction(np.array([WARM]))[0])
    assert 0.0 < f_s < 0.1
    np.testing.assert_allclose(rows[:, C["solid_fraction"]], np.where(mode == 0, 0.0, f_s), rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(rows[:, C["h_kg"]], mat.enthalpy(np.where(mode == 0, LIQUID, WARM)), rtol=1e-12)
    # leaving at the metal's own velocity, and each phase's Weber and Ohnesorge numbers at its radius
    np.testing.assert_allclose(rows[:, C["v_s_ms"]], np.where(mode == 0, 2.0, 0.5), rtol=1e-14)
    np.testing.assert_allclose(rows[:, C["tx"]:C["tz"] + 1], np.where(mode[:, None] == 0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                               atol=1e-15)
    sigma = np.where(mode == 0, mat.liquid.sigma, skinbody.SLURRY_SIGMA)
    mu = np.where(mode == 0, mat.liquid.mu, s["mu"][i])
    np.testing.assert_allclose(rows[:, C["we_d"]], State.freestream.rho * State.V ** 2 * 2.0 * r / sigma, rtol=1e-14)
    np.testing.assert_allclose(rows[:, C["oh"]], mu / np.sqrt(rho * sigma * 2.0 * r), rtol=1e-14)
    assert s["mu"][2] > 1.2 * mat.liquid.mu and s["rho"][2] > mat.liquid.rho
    # the histograms accumulated each sub-step are the records' (in the same bins, to round-off)
    n_rec, m_rec = spray.histogram(rows[:, C["r_m"]], rows[:, C["dn"]], rows[:, C["dm_kg"]])
    np.testing.assert_array_equal(hist[0] > 0.0, n_rec > 0.0)
    np.testing.assert_allclose(hist[0], n_rec, rtol=1e-13, atol=0.0)
    np.testing.assert_allclose(hist[1], m_rec, rtol=1e-13, atol=0.0)


def test_a_spraying_skin_body_s_records_close_on_its_books_and_the_coupled_writer_writes_them(coarse_sphere_mesh,
                                                                                            tmp_path, monkeypatch):
    """The skin body's droplet records (sub-plan 22 Task 2), on Task 1c's spraying test with slurry spraying on: the 4 mm
    test sphere heated at 3 MW/m2 from 890 K under the 100 mm flight's gas state at 93 s, two 0.25 s steps, one flow
    sub-step a heat sub-step. The second step sprays the liquid and the slurry under it, and the body keeps one block of
    rows for it: one row per column and mode, the liquid's mode's first, each at its plan's radius; the rows' mass is
    what the accepted passes released from each column (the slurry's part what the cascade took at the slurry's rate)
    and the books' released mass, their enthalpy the books' released enthalpy, and their momentum, mass times speed
    along their direction, the momentum the droplets took from each column, all to round-off; the size histograms
    accumulated each flow sub-step are the rows'. The coupled writer writes the rows and their size distribution and
    reads them back to the last bit.

    Sub-plan 22 Task 4: the rows go to a spool on disk as they come, which gives back the block it was given to the
    last bit; each carries Register et al.'s release angle of its column's top corner from the centroid of the volume
    the top surface enclosed where the plan was evaluated; the history reports the step's released mass, the droplets
    since time zero, the radius median and largest over the step's rows, the smallest angle and the area of the
    columns whose plan sprays, the release angle's mean by mass -- the mean recomputed from the step's rows -- least and
    largest, and b and the resolved sheared depth over delta_m on the wet windward columns; and from those columns the
    coupled writer writes the summary, whose released mass is the step's rows'."""
    pytest.importorskip("cantera")
    from test_reentry_model_surface_flow import aero_state_at
    from reentry_model import coupled, spray
    from reentry_model import trajectory as tj
    mat = material.Material.from_drama_json("AA7075_scheil")
    np.random.seed(12345)              # pyamg's draws, as a run seeds them: the histograms' round-off is the path's
    b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=spray.SprayModel(mat.liquid),
             slurry_spray=True, spool=str(tmp_path / "spool"))
    plans, books, spooled = [], [], []
    girin, book, append = skinbody.girin_plan, skinbody.SprayPlan.book, spray.SourceSpool.append

    def planned(*a, **k):
        plans.append(girin(*a, **k))
        return plans[-1]

    def booked(plan, fb):
        books.append((plan, fb.released.copy(), fb.released_slurry.copy()))
        return book(plan, fb)
    def kept(spool, block):
        spooled.append(block.copy())
        return append(spool, block)
    monkeypatch.setattr(skinbody, "girin_plan", planned)
    monkeypatch.setattr(skinbody.SprayPlan, "book", booked)
    monkeypatch.setattr(spray.SourceSpool, "append", kept)
    state = aero_state_at(93.0)
    stats = [b.skin_stats()]
    for k in range(2):
        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), state=state)
        stats.append(b.skin_stats())
    assert len(b.source_rows) == 1 and b.particles and isinstance(b.source_rows, spray.SourceSpool)   # the first step melts
    rows = spray.source_table(b.source_rows)
    assert len(spooled) == 1 and np.array_equal(rows, spooled[0], equal_nan=True)     # the spool's rows: the block's
    C = {name: i for i, name in enumerate(spray.SOURCE_COLUMNS)}
    p = plans[-1]
    index = {tuple(x): i for i, x in enumerate(p.points)}
    i = np.array([index[tuple(x)] for x in rows[:, C["x_m"]:C["z_m"] + 1]])
    n_liq, n_sl = int((p.mass > 0.0).sum()), int((p.mass_slurry > 0.0).sum())
    mode = (np.arange(len(rows)) >= n_liq).astype(int)
    assert len(rows) == n_liq + n_sl and n_liq > 100 and n_sl > 100 and (rows[:, C["time_s"]] == 0.5).all()
    assert len(set(i[mode == 0])) == n_liq and len(set(i[mode == 1])) == n_sl
    assert ((p.mass > 0.0) & (p.mass_slurry > 0.0)).sum() > 10                    # columns with both records
    np.testing.assert_array_equal(rows[:, C["r_m"]], np.where(mode == 0, p.result.r[i], p.slurry["r"][i]))
    # the books: the mass each column released and what the cascade took at the slurry's rate, the enthalpy, the momentum
    mine = [x for x in books if x[0] is p]
    released, from_slurry = sum(x[1] for x in mine), sum(x[2] for x in mine)
    by_column = np.bincount(i, rows[:, C["dm_kg"]], b.columns.n)
    np.testing.assert_allclose(by_column, released, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(np.bincount(i[mode == 1], rows[mode == 1, C["dm_kg"]], b.columns.n), from_slurry,
                               rtol=1e-12, atol=1e-21)
    assert rows[:, C["dm_kg"]].sum() == pytest.approx(b.released_mass, rel=1e-12)
    assert (rows[:, C["dm_kg"]] * rows[:, C["h_kg"]]).sum() == pytest.approx(b.released_enthalpy, rel=1e-12)
    velocity = rows[:, C["v_s_ms"], None] * rows[:, C["tx"]:C["tz"] + 1]
    momentum = np.stack([np.bincount(i, rows[:, C["dm_kg"]] * velocity[:, k], b.columns.n) for k in range(3)], axis=1)
    np.testing.assert_allclose(momentum, b.droplet_momentum, rtol=0.0, atol=1e-12 * np.abs(b.droplet_momentum).max())
    # the two new columns: the liquid carries no solid (to the property table's interpolation at the liquidus), the
    # slurry less than half its mass, each at its cells' own enthalpy
    s_f, h = rows[:, C["solid_fraction"]], rows[:, C["h_kg"]]
    h_rigid, h_liq = mat.enthalpy(np.array([mat.T_rigid, mat.T_liquidus]))
    assert (s_f[mode == 0] < 1e-9).all() and (s_f[mode == 1] > 0.0).all() and (s_f[mode == 1] < 0.5).all()
    assert (h[mode == 0] > h_liq * (1.0 - 1e-6)).all() and (h_rigid < h[mode == 1]).all() and (h[mode == 1] < h_liq).all()
    # the size histograms, accumulated each flow sub-step, are the rows': in the same bins, to the round-off of numpy's
    # weighted histogram, which takes each bin as a difference of running sums over all the droplets (measured: 1.6e-12
    # droplets in a bin of 0.065, 2e-17 of the 71 000 in all; with the draws seeded 5.3e-16 of the 89 173 in all, and
    # unseeded, after other tests, up to 1.04e-15, sub-plan 22 Task 4)
    n_rec, m_rec = spray.histogram(rows[:, C["r_m"]], rows[:, C["dn"]], rows[:, C["dm_kg"]])
    np.testing.assert_array_equal(b.hist_n > 0.0, n_rec > 0.0)
    np.testing.assert_allclose(b.hist_n, n_rec, rtol=0.0, atol=1e-15 * n_rec.sum())
    np.testing.assert_allclose(b.hist_m, m_rec, rtol=0.0, atol=1e-15 * m_rec.sum())
    # Task 4: the release angle of each row's column top from the plan's centroid, the top surface's where it was
    # evaluated: 111.5-169.9 degrees, the columns within 10.9 degrees of the nose not spraying, none behind the equator,
    # and within 7.3 degrees of 180 less the corner normal's angle theta_deg (measured)
    angle = rows[:, C["release_angle_register_deg"]]
    np.testing.assert_array_equal(p.centre, spray.enclosed_centroid(p.points, b.surface.faces))
    np.testing.assert_array_equal(angle, spray.register_angle(p.points, p.centre, b.v_hat)[i])
    assert 160.0 < angle.max() < 180.0 and angle.min() > 90.0
    assert (np.abs(angle - (180.0 - rows[:, C["theta_deg"]])) < 10.0).all()
    # the history's droplet columns: nothing in the first step, the second's from its rows and its plan
    s0, s1, s2 = stats
    assert s1["released_step_mass_kg"] == 0.0 and np.isnan(s1["r_median_um"]) and np.isnan(s1["release_angle_register_mean_deg"])
    assert np.isnan(s0["b_over_delta_m_median"]) and s0["n_released"] == 0.0 and s0["spraying_area_m2"] == 0.0
    dm = rows[:, C["dm_kg"]]
    assert s2["released_step_mass_kg"] == pytest.approx(dm.sum(), rel=1e-12) and s2["n_released"] == b.droplet_number
    assert s2["r_median_um"] == np.median(rows[:, C["r_m"]]) * 1e6 and s2["r_max_um"] == rows[:, C["r_m"]].max() * 1e6
    on = p.spraying()
    assert s2["theta_cr_deg"] == np.degrees(p.theta[on].min()) and s2["spraying_area_m2"] == p.area[on].sum() > 0.0
    assert s2["release_angle_register_mean_deg"] == pytest.approx((dm * angle).sum() / dm.sum(), rel=1e-14)
    assert (s2["release_angle_register_min_deg"], s2["release_angle_register_max_deg"]) == (angle.min(), angle.max())
    assert b.sprayed_mass == pytest.approx(b.released_mass, rel=1e-12) and b.n_source_rows == len(rows)
    # b and the resolved sheared depth over delta_m on the wet windward columns, the step's end
    F = b.skins.depth_above(mat.T_rigid)
    delta_m, _ = spray.melt_layer(b.last_flow, mat.liquid)
    wet = (b.theta <= 0.5 * np.pi) & (F > 0.0) & np.isfinite(delta_m) & (delta_m > 0.0)
    assert wet.sum() > 100 and s2["b_over_delta_m_median"] == np.median(F[wet] / delta_m[wet])
    sheared = skinbody.sheared_depth(b.skins)
    moving = wet & np.isfinite(sheared)
    assert s2["sheared_depth_over_delta_m_median"] == np.median(sheared[moving] / delta_m[moving]) > 0.0
    # the coupled writer: the same files as the element model's, read back to the last bit, and the summary from the
    # history's droplet columns
    keys = {k for _, k in coupled.SKIN_SUMMARY_COLUMNS} - {"time_s", "altitude_km", "velocity_kms"}
    columns = {k: np.array([x[k] for x in stats], dtype=float) for k in keys}
    columns.update({"time_s": np.array([0.0, 0.25, 0.5]), "altitude_km": np.full(3, state.h / 1e3),
                    "velocity_kms": np.full(3, state.V / 1e3)})
    history = tj.History(columns, np.zeros((3, 6)), None)
    files = coupled.write_particles(str(tmp_path), b, history)
    assert set(files) == {"particles", "particles_summary", "size_distribution"}
    with open(files["particles_summary"]) as fh:
        summary = [line.split(",") for line in fh.read().splitlines()]
    assert summary[0] == [h for h, _ in coupled.SKIN_SUMMARY_COLUMNS] and len(summary) == 4
    assert float(summary[3][3]) == pytest.approx(dm.sum(), rel=1e-8) and float(summary[2][3]) == 0.0
    assert float(summary[3][4]) == pytest.approx(rows[:, C["dn"]].sum(), rel=1e-8)
    npz = np.load(files["particles"])
    assert npz.files == spray.SOURCE_COLUMNS
    assert all(np.array_equal(npz[name], rows[:, k], equal_nan=True) for k, name in enumerate(spray.SOURCE_COLUMNS))
    with open(files["size_distribution"]) as fh:
        total = [line.split(",") for line in fh.read().splitlines()[1:]][-spray.N_BINS:]
    np.testing.assert_allclose([float(x[4]) for x in total], n_rec, rtol=1e-8)
    np.testing.assert_allclose([float(x[5]) for x in total], m_rec, rtol=1e-8)


def test_the_skins_hold_liquid_once_a_top_cell_is_above_the_liquidus(coarse_sphere_mesh):
    """Sub-plan 22 Task 3b: `liquid_present`, the melting onset at which a run's continuum step latches, is true once any
    column's top cell is above the liquidus -- where the liquid depth, the depth option B reads, is not zero -- and not
    for slurry however warm, nor for a material that does not melt."""
    b = make(coarse_sphere_mesh, T0=890.0)
    mat = b.material
    assert not b.liquid_present() and not (b.skins.depth_above(mat.T_liquidus) > 0.0).any()
    b.skins.T[3, 0] = mat.T_liquidus - 1e-6                      # slurry a micro-kelvin under the liquidus
    assert not b.liquid_present()
    b.skins.T[3, 0] = mat.T_liquidus + 1e-6
    assert b.liquid_present() and (b.skins.depth_above(mat.T_liquidus) > 0.0).sum() == 1
    copy = mesh.VolumeMesh(coarse_sphere_mesh.points.copy(), coarse_sphere_mesh.tets.copy(), dict(coarse_sphere_mesh.params))
    cold = skinbody.SkinBody(copy, material.Material.from_drama_json("AA7075_nomelt"), thermal_solver("skfem"), 1.47, DEV,
                             T0=1200.0, v_hat=(-1.0, 0.0, 0.0))
    assert not cold.liquid_present()                             # a material that does not melt holds no liquid


# -- the outputs of a spraying skin flight (sub-plan 22 Task 4) ----------------------------------------------------------
def test_the_resolved_sheared_depth_is_the_straight_line_depth_carrying_the_computed_flow():
    """Spec 6.7's resolved sheared depth (`skinbody.sheared_depth`, Task 4b, Asha 2026-10-31): twice the flow per unit
    width over the top cell's speed -- the depth of the straight line from the top cell's speed to rest that carries the
    computed flow -- with no cut-off. Columns of 11.7 um cells: a profile falling linearly to rest 0.2 mm down, sampled at
    the cells' centres, D / (1 - z0 / D) with z0 the top cell centre's depth, 5.9 um (on Girin's linear profile over
    delta_m, delta_m and about half a top cell), within the sampling of the cell where it comes to rest; a skin moving as
    one piece, twice its thickness; a fast top of five cells over a creep at a hundredth of its speed, the creep counted
    by the flow it carries, 10.7 cells where a tenth of the top's speed read 5.4; a top at rest, nan."""
    from test_reentry_model_skin import LIQUID, RIGID, stacked, strip_grid
    col, _ = strip_grid(3, 2, h=20e-3)
    sk = stacked(col, [[LIQUID] * 30 + [RIGID] * 10, [LIQUID] * 40, [LIQUID] * 40] + [[LIQUID] * 40] * 3)
    dz = sk.dz()
    z = np.cumsum(dz, axis=1) - 0.5 * dz
    D = 2e-4
    U = np.zeros(sk.U.shape)
    U[0, :, 0] = 3.0 * np.clip(1.0 - z[0] / D, 0.0, None)
    U[1, :, 1] = 0.5
    U[2, :, 0] = np.where(np.arange(40) < 5, 3.0, 0.03)
    sk.U = U * sk.valid[..., None]
    assert not hasattr(skinbody, "SHEARED_FRACTION")
    got = skinbody.sheared_depth(sk)
    assert got[0] == pytest.approx(D / (1.0 - z[0, 0] / D), rel=1e-4) and 5e-6 < z[0, 0] < 7e-6
    assert got[0] - D == pytest.approx(z[0, 0], rel=0.05)
    assert got[1] == pytest.approx(2.0 * sk.thickness()[1], rel=1e-14)
    assert np.allclose(dz[2], dz[2, 0]) and got[2] == pytest.approx(10.7 * dz[2, 0], rel=1e-12)
    assert np.isnan(got[3:]).all()


# -- spreading the deposit (sub-plan 22 Task 5c; Asha 2026-11-06) ----------------------------------------------------------
def melting_body(mesh, **kw):
    """The melting test body of the spraying tests above, 890 K under the 100 mm flight's gas state at 93 s, one flow
    sub-step a heat sub-step, flown for two 0.125 s macro steps at 3 MW/m2, numpy's generator seeded as a run seeds it: it
    hands down 3.8 g. Returns the body and its starting mass."""
    pytest.importorskip("cantera")
    from test_reentry_model_surface_flow import aero_state_at
    from reentry_model import spray
    np.random.seed(12345)
    mat = material.Material.from_drama_json("AA7075_scheil")
    b = make(mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=spray.SprayModel(mat.liquid), **kw)
    M0 = b.mass()
    state = aero_state_at(93.0)
    for k in range(2):
        b.advance(0.125 * (k + 1), 0.125, front_loads(b, 3e6), state=state)
    return b, M0


def zero_book(sk):
    """What no draw left: the base's operations alone, as the deposit tests above run them."""
    return skin.SkinBook(*(np.zeros(sk.n) for _ in range(8)), momentum_out=np.zeros((sk.n, 3)))


def bump(b, peak, share=0.9, T=700.0):
    """A deposit at `T` on the bottoms of the columns `peak`, `share` of a part's allowance deep."""
    mat, sk, col = b.material, b.skins, b.columns
    extra = np.zeros(sk.n)
    extra[peak] = share * skinbody.MOVE_MAX * col.edge_length.mean() * mat.rho * col.area[peak]
    sk.add_bottom(extra, float(mat.enthalpy(T)))
    sk.regrid()


def spied(b, monkeypatch):
    """The arguments of the base move's spreading and of its corner move, as they come: ([(hand-down given, Spread)], [the
    enthalpy each outward sweep is to be booked at])."""
    spreads, outward = [], []
    spread, policy = skinbody.spreading.spread_hand_down, b.moving.move_with_policy

    def recorded(*a, **k):
        spreads.append((a[6], spread(*a, **k)))
        return spreads[-1][1]

    def recorded_move(u, h_inward=None, h_outward=None):
        outward.append(h_outward)
        return policy(u, h_inward, h_outward)
    monkeypatch.setattr(skinbody.spreading, "spread_hand_down", recorded)
    b.moving.move_with_policy = recorded_move
    return spreads, outward


def neighbours(b, corners):
    """The corners sharing an edge with any of `corners`, those left out."""
    e = b.columns.edges
    near = np.unique(np.r_[e[np.isin(e[:, 0], corners), 1], e[np.isin(e[:, 1], corners), 0]])
    return near[~np.isin(near, corners)]


def test_a_limit_that_does_not_bind_changes_nothing_to_the_last_bit(coarse_sphere_mesh):
    """Nothing folds near 179.9 degrees, so the melting body flown with that limit is the body flown without one in every
    array and every shared history column, and only the four spreading columns are added, recording nothing moved. Without
    a limit there are no spreading columns: every run without the option keeps its column set."""
    off, _ = melting_body(coarse_sphere_mesh)
    on, _ = melting_body(coarse_sphere_mesh, spread_limit=179.9)
    assert off.handed_mass > 0.0 and on.spread_limit == 179.9 and off.spread_limit is None
    assert same_body(off, on)
    a, c = off.skin_stats(), on.skin_stats()
    same = lambda x, y: x == y or (np.isnan(x) and np.isnan(y))
    assert all(same(a[k], c[k]) for k in a if not k.startswith("wall_"))
    assert set(c) - set(a) == set(skinbody.SPREAD_COLUMNS) and not set(a) & set(skinbody.SPREAD_COLUMNS)
    assert [c[k] for k in skinbody.SPREAD_COLUMNS] == [0.0, 0.0, 0, 0.0]


def test_a_limit_outside_a_fold_s_range_is_refused(coarse_sphere_mesh):
    for bad in (0.0, 180.0, -1.0, float("nan")):
        with pytest.raises(ValueError):
            make(coarse_sphere_mesh, spread_limit=bad)


def test_spreading_keeps_the_books(coarse_sphere_mesh):
    """The melting body with a limit low enough to bind on the 4 mm mesh, whose sphere folds 14 degrees before anything happens:
    0.6 g of the 3.8 g handed down moves between corners, and the mass closes to 1e-12, the energy balance residual is below
    1e-8 and the momentum books close, the hot landing's heat in transit counted (measured: mass exact, energy 1.7e-11,
    momentum 1.4e-17)."""
    b, M0 = melting_body(coarse_sphere_mesh, spread_limit=6.0)
    stats = b.skin_stats()
    assert stats["spread_mass_kg"] == b.spread_mass > 0.0 and stats["spread_edges_over"] > 0
    assert b.mass() + b.released_mass == pytest.approx(M0, rel=1e-12)
    assert abs(b.energy_balance_residual()) < 1e-8
    assert stats["momentum_residual"] < 1e-12


def test_metal_landing_hot_takes_its_heat_to_the_receiving_corner_and_counts_no_fault(coarse_sphere_mesh, monkeypatch):
    """Asha's decision of 2026-11-06: spread metal lands hot, because the lower neighbours of a toe are the cold corners ahead
    of it. A proud corner of the 4 mm sphere at 700 K with its neighbours' nodes 40 K colder, past the 5 K fault tolerance:
    the metal crosses each receiving base at that node's enthalpy, so the moving mesh shares no correction over the zone for
    it, and the rest of its heat, the excess, waits in `pending` at that corner alone, which the books count and the next
    part's 3D solve takes as a load there. No fault is counted, the largest landing is recorded, and the energy is exact
    before and after the 3D model takes the heat."""
    b = make(coarse_sphere_mesh, T0=700.0, spread_limit=10.0)
    b.emissivity = 0.0
    mat, sk, col, nodes = b.material, b.skins, b.columns, b.moving.corners
    peak = np.argsort(b.theta)[:1]
    ring = neighbours(b, peak)
    T = b.solver.temperature().copy()
    T[nodes[ring]] -= 40.0
    b.solver.set_temperature(T)
    bump(b, peak)
    spreads, outward = spied(b, monkeypatch)
    E0, T_node = b.energy(), T[nodes]
    b._move_base(zero_book(sk), skin.minimum_mass(mat, sk.area, DEV), b._normal_gradient())
    handed, s = spreads[0]
    rose = s.landed
    cold = rose & np.isin(np.arange(sk.n), ring)
    assert cold.any() and s.moved > 0.0
    excess = s.enthalpy[rose] - s.handed[rose] * mat.enthalpy(T_node[rose])
    assert excess[np.isin(np.flatnonzero(rose), ring)].min() > 0.0
    np.testing.assert_allclose(b.pending[rose], excess, rtol=1e-12)
    assert (b.pending[~rose] == 0.0).all()
    np.testing.assert_allclose(outward[0][cold], mat.enthalpy(T_node[cold]), rtol=1e-12)    # at the node's own enthalpy
    assert b.jump_faults == 0
    stats = b.skin_stats()
    assert 30.0 < stats["spread_jump_max_K"] <= 40.0 + skin.JUMP_TOLERANCE
    assert b.energy() == pytest.approx(E0, rel=1e-12)
    T_before = b.solver.temperature()[nodes].copy()
    zero = heating.HeatingResult(np.zeros(col.n), 0.0, 0.0, 0.0, 0.0)
    b.advance(0.01, 0.01, zero)
    assert b.energy() == pytest.approx(E0, rel=1e-12)
    assert (b.solver.temperature()[nodes[cold]] > T_before[cold]).all()
    assert b.skin_stats()["spread_jump_max_K"] < stats["spread_jump_max_K"]               # since the last row


def test_a_corner_that_takes_hot_metal_and_gives_more_onward_still_sends_its_landing_heat_to_its_node(coarse_sphere_mesh,
                                                                                                       monkeypatch):
    """A two-step cascade on the 4 mm sphere at 700 K (Finding 1 of Task 5c's review): a peak with its node and metal 60 K
    hotter stands above a ring of corners, themselves above the rest. The ring takes the peak's hot metal and gives more onward than
    it took, so each of its corners that keeps any ends handing down less than it started with. The heat it took is still its
    landing, as Asha decided: its metal crosses at the node's enthalpy, the rest goes to that node through `pending`, not through
    `h_outward` and so over the whole moving zone, and the landing counts in the largest recorded. Energy closes and no fault is
    counted."""
    b = make(coarse_sphere_mesh, T0=700.0, spread_limit=15.0)
    b.emissivity = 0.0
    mat, sk, col, nodes = b.material, b.skins, b.columns, b.moving.corners
    peak = np.argsort(b.theta)[:1]
    ring = neighbours(b, peak)
    u = np.zeros((col.n, 3))
    u[peak], u[ring] = 2.4e-3 * col.normal[peak], 1.2e-3 * col.normal[ring]
    b.moving.move_with_policy(u)
    b._refresh_geometry()
    nodes = b.moving.corners
    T = b.solver.temperature().copy()
    T[nodes[peak]] += 60.0
    b.solver.set_temperature(T)
    bump(b, peak, share=0.9, T=760.0)
    bump(b, ring, share=0.5)
    spreads, outward = spied(b, monkeypatch)
    E0, T_node = b.energy(), b.solver.temperature()[nodes]
    b._move_base(zero_book(sk), skin.minimum_mass(mat, sk.area, DEV), b._normal_gradient())
    handed, s = spreads[0]
    on = ring[(s.handed[ring] < handed[ring]) & (s.handed[ring] > 0.0)]
    assert on.size > 0
    jump = mat.temperature_from_enthalpy(s.enthalpy[on] / s.handed[on]) - T_node[on]
    assert jump.max() > skin.JUMP_TOLERANCE
    excess = s.enthalpy[on] - s.handed[on] * mat.enthalpy(T_node[on])
    np.testing.assert_allclose(b.pending[on], excess, rtol=1e-12)
    np.testing.assert_allclose(outward[0][on], mat.enthalpy(T_node[on]), rtol=1e-12)
    assert b.skin_stats()["spread_jump_max_K"] >= jump.max() and b.jump_faults == 0
    assert b.energy() == pytest.approx(E0, rel=1e-12)


def test_a_donor_that_gives_away_all_its_hand_down_keeps_the_books(coarse_sphere_mesh, monkeypatch):
    """A corner 1.5 mm proud of its neighbours with little to hand down -- 0.2 of a part's allowance, some of which goes to the
    minimum its larger dual cell now needs -- gives all it has and nothing lands on it, being the highest: its hand-down is exactly nothing, so the moving mesh has no enthalpy for it (`h_outward` nan) and books any sweep
    there at the corner's own, and the reconciliation and the fault count hold: mass and energy exact, nothing not finite, and
    no fault."""
    b = make(coarse_sphere_mesh, T0=700.0, spread_limit=15.0)
    b.emissivity = 0.0
    mat, sk, col = b.material, b.skins, b.columns
    peak = np.argsort(b.theta)[:1]
    u = np.zeros((col.n, 3))
    u[peak] = 1.5e-3 * col.normal[peak]
    b.moving.move_with_policy(u)
    b._refresh_geometry()
    bump(b, peak, share=0.2)
    spreads, outward = spied(b, monkeypatch)
    M0, E0 = b.mass(), b.energy()
    b._move_base(zero_book(sk), skin.minimum_mass(mat, sk.area, DEV), b._normal_gradient())
    handed, s = spreads[0]
    assert handed[peak] > 0.0 and s.handed[peak] == 0.0
    gone = (handed > 0.0) & (s.handed == 0.0)
    assert gone.any() and (s.enthalpy[gone] == 0.0).all() and np.isnan(outward[0][gone]).all()
    assert b.mass() == pytest.approx(M0, rel=1e-12) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert np.isfinite(b.solver.temperature()).all() and np.isfinite(sk.T[sk.valid]).all() and np.isfinite(b.pending).all()
    assert b.jump_faults == 0


def test_a_body_saved_without_spreading_attributes_flies_on(coarse_sphere_mesh):
    """The two real states were pickled before spreading existed, so a body from them has none of its attributes on the
    instance: the class's defaults give them, it flies on as a body without a limit, and it takes a limit a resume gives it."""
    b = make(coarse_sphere_mesh, T0=700.0)
    b.emissivity = 0.0
    for key in [k for k in vars(b) if "spread" in k]:
        del vars(b)[key]
    assert not [k for k in vars(b) if "spread" in k]
    assert b.spread_limit is None
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.advance(0.5, 0.5, zero)
    assert not set(b.skin_stats()) & set(skinbody.SPREAD_COLUMNS) and b.spread_stats() == {}
    b.spread_limit = 60.0
    b.advance(1.0, 0.5, zero)
    assert set(b.spread_stats()) == set(skinbody.SPREAD_COLUMNS) and set(skinbody.SPREAD_COLUMNS) <= set(b.skin_stats())


# -- what sets the next macro step's part count (sub-plan 22 Task 5d; Asha 2026-11-08) ---------------------------------------
def held_back_body(mesh, depth_allowances=9.0, **kw):
    """The 4 mm sphere at 700 K with thirty front columns holding `depth_allowances` move allowances of metal at the node's
    temperature at their bottoms and no heat: every part's hand-down is capped at just under one allowance (a deposit that
    deep is not down in the macro steps flown), and nothing is drawn. Returns the body and the columns' deposit [kg]."""
    b = make(mesh, T0=700.0, **kw)
    b.emissivity = 0.0
    mat, sk, col = b.material, b.skins, b.columns
    cols = np.argsort(b.theta)[:30]
    extra = np.zeros(sk.n)
    extra[cols] = depth_allowances * skinbody.MOVE_MAX * col.edge_length.mean() * mat.rho * col.area[cols]
    sk.add_bottom(extra, float(mat.enthalpy(700.0)))
    sk.regrid()
    return b, extra.sum()


def timed_out(stats):
    """A history row without its wall-time columns, which time the run and not the flight."""
    return {k: v for k, v in stats.items() if not k.startswith("wall_")}


def same_rows(a, b):
    """Two lists of (parts, row) the same to the last bit, a column of not-a-number (nothing sprayed yet) the same as another."""
    return len(a) == len(b) and all(
        pa == pb and ra.keys() == rb.keys() and all(ra[k] == rb[k] or (ra[k] != ra[k] and rb[k] != rb[k]) for k in ra)
        for (pa, ra), (pb, rb) in zip(a, b))


def fly_three_macro_steps(b, dt=0.05):
    """Three macro steps of `dt` s (not longer than the melting part length) with no heat; the parts each ran as."""
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    parts = []
    for k in range(3):
        b.advance(dt * (k + 1), dt, zero)
        parts.append(b.step_parts)
    return parts


def test_under_the_moves_rule_a_capped_hand_down_grows_the_next_steps_parts_by_a_quarter(coarse_sphere_mesh):
    """The defect of the spreading flights, pinned (sub-plan 22 Task 5d): a hand-down capped just under the move allowance
    gives a move ratio just under 1 in every part, and the hint is 1.25 times the parts of this step, so the next step
    runs a quarter more parts whose hand-downs are capped again -- 837, 1 047, ..., 9 765 parts on the 50 mm spreading flight,
    each the previous times 1.25 rounded up (measured 2026-11-08). Here a deposit nine allowances deep, three 0.05 s
    macro steps (under the 0.0625 s part length, so that the length does not set the count): 1, 2 and 3 parts, each the
    previous times 1.25 rounded up, the allowance of HAND_DOWN_SHARE's 1e-6 included."""
    assert skinbody.PART_MAX >= 0.05
    b, _ = held_back_body(coarse_sphere_mesh)
    assert b.part_count_from == "moves"
    parts = fly_three_macro_steps(b)
    assert parts[0] == 1
    for before, after in zip(parts, parts[1:]):
        assert after == int(np.ceil(1.25 * before * skinbody.HAND_DOWN_SHARE)) and after > before
    assert b.held_mass > 0.0                                        # the deposit still waits: every hand-down was capped


def test_under_the_melting_rule_a_capped_hand_down_leaves_the_part_count_alone(coarse_sphere_mesh):
    """The same flight under `melting`: the hint reads the draws alone, and nothing is drawn, so every macro step runs as the
    one part the draws and the part length call for, while the same metal comes down (a capped hand-down an allowance a
    part, the rest held) and mass and energy close (sub-plan 22 Task 5d)."""
    b, extra = held_back_body(coarse_sphere_mesh, part_count_from="melting")
    M0, E0 = b.mass(), b.energy()
    assert fly_three_macro_steps(b) == [1, 1, 1] and b.part_count_from == "melting"
    assert b.handed_mass > 0.0 and b.held_mass > 0.0
    assert b.mass() == pytest.approx(M0, rel=1e-13) and b.energy() == pytest.approx(E0, rel=1e-12)
    assert abs(b.energy_balance_residual()) < 1e-8 and b.jump_faults == 0


def test_a_draw_over_the_allowance_splits_a_part_under_either_rule(coarse_sphere_mesh, monkeypatch):
    """The split test is the same under both rules (sub-plan 22 Task 5d): 6 MW/m2 from 880 K for three 0.5 s macro steps on
    4 mm elements, with no melting cap on the parts, draws more than the allowance in a part, which is run again as shorter
    parts. The two rules differ only through the hand-down, so it is taken out of the comparison (the cap given to
    `skin.hand_down` is zero: the skins only grow, the first step of such a flight is the one the unbroken run splits
    into 8 parts), and then the two give the same parts and the same history rows in every step. numpy's generator is
    seeded as a run seeds it, since the heat solver's multigrid hierarchy draws from it."""
    hand_down = skin.hand_down
    monkeypatch.setattr(skin, "hand_down", lambda *a, cap=None, **k: hand_down(*a, cap=np.zeros_like(cap), **k))
    runs = {}
    for rule in skinbody.PART_COUNT_RULES:
        np.random.seed(12345)
        b = make(coarse_sphere_mesh, part_count_from=rule, part_max=None)
        rows = []
        for k in range(3):
            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 6e6))
            rows.append((b.step_parts, timed_out(b.skin_stats())))
        runs[rule] = (b, rows)
    moves, melting = runs["moves"], runs["melting"]
    assert moves[0].handed_mass == 0.0 and moves[0].drawn_mass > 0.0
    assert moves[1][0][0] > 1                                          # the draws split the first step's part
    assert same_rows(moves[1], melting[1])
    assert moves[0].mass() == melting[0].mass() and moves[0].energy() == melting[0].energy()


def test_the_moves_rule_by_default_changes_nothing_to_the_last_bit(coarse_sphere_mesh):
    """`part_count_from="moves"` and the argument left out are one flight, every column of every row and the final mass and
    energy equal to the last bit, the wall-time columns that time the run left out (sub-plan 22 Task 5d), on a flight that melts, hands down and splits parts (numpy's generator
    seeded as a run seeds it before each)."""
    runs = []
    for kw in ({}, {"part_count_from": "moves"}):
        np.random.seed(12345)
        b = make(coarse_sphere_mesh, **kw)
        rows = []
        for k in range(3):
            b.advance(0.5 * (k + 1), 0.5, front_loads(b, 3e6))
            rows.append((b.step_parts, timed_out(b.skin_stats())))
        runs.append((b, rows))
    assert same_rows(runs[0][1], runs[1][1]) and runs[0][1][1][0] > 1
    assert runs[0][0].mass() == runs[1][0].mass() and runs[0][0].energy() == runs[1][0].energy()


def test_a_macro_step_that_would_need_more_than_PARTS_MAX_parts_stops_the_run(coarse_sphere_mesh, monkeypatch):
    """The stall guard (sub-plan 22 Task 5d; the controller, 2026-11-08): a macro step that would start with more than
    PARTS_MAX parts raises RuntimeError naming the step, the count and the limit, before any part runs, under either rule;
    PARTS_MAX parts themselves run. The limit is lowered to 3 so that the test does not run a thousand parts."""
    monkeypatch.setattr(skinbody, "PARTS_MAX", 3)
    for rule in skinbody.PART_COUNT_RULES:
        b = make(coarse_sphere_mesh, T0=600.0, part_count_from=rule)
        b.emissivity = 0.0
        zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
        M0, E0 = b.mass(), b.energy()
        b._hint = 4
        with pytest.raises(RuntimeError, match=r"from 0\.5 s to 1 s would run as 4 parts, more than PARTS_MAX = 3"):
            b.advance(1.0, 0.5, zero)
        assert b.step_parts == 0 and b.mass() == M0 and b.energy() == E0       # nothing ran
        b._hint = 3
        b.advance(1.0, 0.5, zero)
        assert b.step_parts == 3
    assert skinbody.PARTS_MAX == 3 and skinbody.PART_COUNT_RULES == ("moves", "melting")


def test_a_body_saved_without_the_part_count_rule_flies_on(coarse_sphere_mesh):
    """The states saved before the rule existed have no `part_count_from` on the instance: the class's default gives them
    `moves`, they fly on as before, and a resume can give them `melting` (sub-plan 22 Task 5d)."""
    b, _ = held_back_body(coarse_sphere_mesh)
    del vars(b)["part_count_from"]
    assert "part_count_from" not in vars(b) and b.part_count_from == "moves"
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.advance(0.05, 0.05, zero)
    b.advance(0.1, 0.05, zero)
    assert b.step_parts == 2                                         # the moves rule: the hand-down fed the hint
    b.part_count_from = "melting"
    b.advance(0.15, 0.05, zero)
    b.advance(0.2, 0.05, zero)
    assert b.step_parts == 1


def test_an_unknown_part_count_rule_is_refused(coarse_sphere_mesh):
    with pytest.raises(ValueError, match="part_count_from"):
        make(coarse_sphere_mesh, part_count_from="draws")


# -- the hand-down switched off (sub-plan 22 Task 5e; Asha 2026-11-08) ------------------------------------------------------
def fly_macro_steps(b, steps=2, dt=0.05, loads=None):
    """`steps` macro steps of `dt` s (not longer than the melting part length), with no heat unless `loads` is given; the 3D
    model's enclosed volume after each [m3]."""
    loads = loads if loads is not None else heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    volumes = []
    for k in range(steps):
        b.advance(dt * (k + 1), dt, loads)
        volumes.append(b.mesh.volume())
    return volumes


def test_with_the_hand_down_off_no_metal_crosses_the_base_outward(coarse_sphere_mesh):
    """Asha's decision of 2026-11-08 (sub-plan 22 Task 5e): with `hand_down=False` nothing is handed down, so no metal is
    added to the 3D model. On the deposit of `held_back_body` (nine move allowances of metal at the node's temperature,
    no heat, nothing drawn), two 0.05 s macro steps: the default hands down and the enclosed volume of the 3D model grows
    in each; with the hand-down off `handed_mass` stays 0 and the volume does not grow from one step to the next."""
    on, _ = held_back_body(coarse_sphere_mesh)
    off, _ = held_back_body(coarse_sphere_mesh, hand_down=False)
    start = off.mesh.volume()
    assert on.mesh.volume() == start
    v_on, v_off = fly_macro_steps(on), fly_macro_steps(off)
    assert on.handed_mass > 0.0 and v_on[0] > start and v_on[1] > v_on[0]
    assert off.handed_mass == 0.0 and v_off[0] <= start * (1.0 + 1e-12) and v_off[1] <= v_off[0] * (1.0 + 1e-12)


def test_with_the_hand_down_off_the_rigid_metal_stays_in_its_skin_and_is_counted_as_held(coarse_sphere_mesh):
    """The metal the default hands down stays in its skin with the hand-down off, so the skins' mass above their minimum is at
    least the default's after the same steps, and `held_mass` is every kilogram `skin.handable` finds, with no limit or cap
    (it is found on the skins as they stand before the first step: nothing is drawn or heated, so they do not change)."""
    on, extra = held_back_body(coarse_sphere_mesh)
    off, _ = held_back_body(coarse_sphere_mesh, hand_down=False)
    m_min = skin.minimum_mass(off.material, off.columns.area, off.settings)
    able = skin.handable(off.skins, m_min, off.material.T_rigid)
    area = off.columns.area.copy()               # the move changes the areas, and the held depth reads those before it
    assert able.sum() > 0.0
    fly_macro_steps(on, steps=1)
    fly_macro_steps(off, steps=1)
    above = lambda b: b.skins.mass() - m_min.sum()
    assert above(off) >= above(on) and above(off) > above(on) + 0.5 * on.handed_mass
    assert off.held_mass == pytest.approx(able.sum(), rel=1e-9) and off.held_mass > on.held_mass
    assert off.held_depth_max == pytest.approx((able / (off.material.rho * area)).max(), rel=1e-9)


def test_with_the_hand_down_off_the_books_close(coarse_sphere_mesh):
    """Skins that only grow still close the books (mass 1e-12, energy 1e-8, momentum) on the held-back deposit flown for two
    0.05 s macro steps, the second heated at 3 MW/m2 so that draws and the move run too."""
    b, _ = held_back_body(coarse_sphere_mesh, hand_down=False)
    M0 = b.mass()
    b.energy0 = E0 = b.energy()                  # the deposit was added after the body was built: its heat is the start's
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.advance(0.05, 0.05, zero)
    b.advance(0.1, 0.05, front_loads(b, 3e6))
    assert b.handed_mass == 0.0 and b.drawn_mass >= 0.0
    assert b.mass() == pytest.approx(M0, rel=1e-12) and b.jump_faults == 0
    assert abs(b.energy_balance_residual()) < 1e-8 and abs(b.momentum_residual()) < 1e-8
    assert b.energy() == pytest.approx(E0 + b.absorbed_heat - b.removed_enthalpy - b.released_enthalpy, rel=1e-8)


def test_the_hand_down_on_by_default_changes_nothing_to_the_last_bit(coarse_sphere_mesh):
    """`hand_down=True` and the argument left out are one flight: every column of every row and the final mass and energy
    equal to the last bit on the held-back deposit, which hands down (sub-plan 22 Task 5e)."""
    runs = []
    for kw in ({}, {"hand_down": True}):
        b, _ = held_back_body(coarse_sphere_mesh, **kw)
        zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
        rows = []
        for k in range(2):
            b.advance(0.05 * (k + 1), 0.05, zero)
            rows.append((b.step_parts, timed_out(b.skin_stats())))
        runs.append((b, rows))
    assert same_rows(runs[0][1], runs[1][1]) and runs[0][0].handed_mass > 0.0 and runs[0][0].hand_down is True
    assert runs[0][0].mass() == runs[1][0].mass() and runs[0][0].energy() == runs[1][0].energy()


def test_the_hand_down_off_refuses_a_spread_limit(coarse_sphere_mesh):
    """Spreading moves only metal being handed down (Asha, 2026-11-08), so a body with neither is refused."""
    with pytest.raises(ValueError, match="hand_down"):
        make(coarse_sphere_mesh, hand_down=False, spread_limit=60.0)
    assert make(coarse_sphere_mesh, hand_down=False).spread_limit is None
    assert make(coarse_sphere_mesh, spread_limit=60.0).hand_down is True


def test_a_body_saved_without_the_hand_down_switch_flies_on(coarse_sphere_mesh):
    """The states saved before the switch existed have no `hand_down` on the instance: the class's default gives them the
    hand-down on, they fly on as before, and a resume can switch it off (sub-plan 22 Task 5e)."""
    b, _ = held_back_body(coarse_sphere_mesh)
    del vars(b)["hand_down"]
    assert "hand_down" not in vars(b) and b.hand_down is True
    fly_macro_steps(b, steps=1)
    assert b.handed_mass > 0.0
    handed = b.handed_mass
    b.hand_down = False
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.advance(0.1, 0.05, zero)
    assert b.handed_mass == handed


# -- sideways heat between cells that touch (sub-plan 22 Task 5f; Asha 2026-11-08) ------------------------------------------
def maps_made(monkeypatch):
    """The names of the sideways maps the body makes, as it makes them."""
    made = []

    def counted(cls):
        def make_map(*args, **kwargs):
            made.append(cls.__name__)
            return cls(*args, **kwargs)
        return make_map
    for name in ("Sideways", "SidewaysByHeight"):
        monkeypatch.setattr(skin, name, counted(getattr(skin, name)))
    return made


def test_height_pairing_keeps_the_books(coarse_sphere_mesh, monkeypatch):
    """`sideways_pairing="height"` (sub-plan 22 Task 5f) on a melting body that flows: 3 MW/m2 on the front of the 4 mm
    sphere from 893 K, just under T_rigid, for two 0.0625 s macro steps (the melting part length), the liquid kept and
    flowing under a Newtonian push, one flow sub-step a heat sub-step, so the cells' thicknesses change every sub-step and
    metal is drawn and handed down. Every sideways map is the height map; the mass closes to 1e-12, the energy to 1e-8 of
    the heat absorbed and the momentum to round-off."""
    made = maps_made(monkeypatch)
    b = make(coarse_sphere_mesh, T0=893.0, flow=True, remove=False, flow_cycles=1, sideways_pairing="height")
    M0 = b.mass()
    for k in range(2):
        b.advance(0.0625 * (k + 1), 0.0625, front_loads(b, 3e6), forcing=newtonian(b))
    assert b.sideways_pairing == "height" and made and set(made) == {"SidewaysByHeight"}
    assert b.drawn_mass > 0.0 and b.handed_mass > 0.0 and b.skin_stats()["flowing_depth_max_mm"] > 0.0
    assert b.mass() == pytest.approx(M0, rel=1e-12)
    assert abs(b.energy_balance_residual()) < 1e-8
    assert b.momentum_residual() < 1e-12 and b.impulse_size > 0.0


def test_the_depth_pairing_by_default_changes_nothing_to_the_last_bit(coarse_sphere_mesh, monkeypatch):
    """`sideways_pairing="depth"` and the argument left out are one flight, every column of every row and the final mass
    and energy equal to the last bit, every map the depth map (sub-plan 22 Task 5f): 3 MW/m2 on the front from 893 K for two
    0.0625 s macro steps, the liquid removed as it forms, numpy's generator seeded as a run seeds it before each."""
    made = maps_made(monkeypatch)
    runs = []
    for kw in ({}, {"sideways_pairing": "depth"}):
        np.random.seed(12345)
        b = make(coarse_sphere_mesh, T0=893.0, **kw)
        rows = []
        for k in range(2):
            b.advance(0.0625 * (k + 1), 0.0625, front_loads(b, 3e6))
            rows.append((b.step_parts, timed_out(b.skin_stats())))
        runs.append((b, rows))
    assert same_rows(runs[0][1], runs[1][1]) and runs[0][0].removed_mass > 0.0 and runs[0][0].sideways_pairing == "depth"
    assert runs[0][0].mass() == runs[1][0].mass() and runs[0][0].energy() == runs[1][0].energy()
    assert made and set(made) == {"Sideways"} and skinbody.SkinBody.sideways_pairing == "depth"


def test_a_body_saved_without_the_sideways_pairing_flies_on(coarse_sphere_mesh, monkeypatch):
    """The states saved before the pairing existed have no `sideways_pairing` on the instance: the class's default gives
    them the depth pairing, they fly on as before, and a resume can give them the height pairing (sub-plan 22 Task 5f). An
    unknown pairing is refused."""
    made = maps_made(monkeypatch)
    b = make(coarse_sphere_mesh, T0=700.0)
    b.emissivity = 0.0
    del vars(b)["sideways_pairing"]
    assert "sideways_pairing" not in vars(b) and b.sideways_pairing == "depth"
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.advance(0.05, 0.05, zero)
    assert made and set(made) == {"Sideways"}
    del made[:]
    b.sideways_pairing = "height"
    b.advance(0.1, 0.05, zero)
    assert made and set(made) == {"SidewaysByHeight"}
    with pytest.raises(ValueError, match="sideways_pairing"):
        make(coarse_sphere_mesh, sideways_pairing="width")


# -- the windward loads from a sphere fitted to the skins' tops (sub-plan 22 Task 5g, Asha 2026-11-10) --------------------


def with_tops(b, points=None, theta=None):
    """Give the body the top surface of `points` (the corners' tops as they are by default) and `theta` (recomputed from the
    surface by default) and refit the nose as `_refresh_geometry` does: a way to shape the windward cap without flying to it."""
    from reentry_model import columns as columns_mod
    pts = b.surface.points.copy() if points is None else points
    b.surface = columns_mod.CornerSurface.from_points(pts, b.columns.faces)
    b.theta = b.surface.angles_to(b.v_hat) if theta is None else theta
    b._leeward_raw = None                           # the centroid test of the new surface, made at the first read
    b._nose_raw = b._fit_nose()
    return b


def flatten_front(b, factor=2.0, rho_max=0.6):
    """The tops of the front half within `rho_max` of the sphere's radius from the flight axis moved onto a sphere of `factor`
    times its radius touching it at the nose; returns the new points and the sphere's radius [m]."""
    R = float(b.mesh.params["radius_m"])
    pts = b.surface.points.copy()
    c0, v = b.mesh.points.mean(axis=0), b.v_hat
    X = pts - c0
    ax = X @ v
    rho = np.sqrt(np.maximum((X * X).sum(axis=1) - ax ** 2, 0.0))
    front = (ax > 0.0) & (rho < rho_max * R)
    new_ax = -R * (factor - 1.0) + np.sqrt((factor * R) ** 2 - rho[front] ** 2)
    pts[front] += np.outer(new_ax - ax[front], v)
    return pts, factor * R


def test_the_fitted_nose_radius_of_an_untouched_sphere_is_its_own_radius(coarse_sphere_mesh):
    """The 4 mm sphere's tops are the bases moved out by the 0.8 mm skins: the cap within 30 degrees of the flow (about 140
    corners) fits a sphere within 2 % of the mesh's 50 mm radius, which `nose_radius()` returns under the default
    `current`; the bound is 1.67 times the windward transverse radius and the cap is convex (sub-plan 22 Task 5g)."""
    b = make(coarse_sphere_mesh)
    R, f = 0.05, b.nose_fit
    assert f.points >= 50 and not f.concave and abs(f.fitted / R - 1.0) < 0.02
    assert b.size_feedback == "current" and b.nose_radius() == f.used == f.fitted
    assert f.bound == pytest.approx(1.67 * R, rel=0.03)


def test_a_flattened_windward_cap_gives_the_radius_of_the_flattened_cap(coarse_sphere_mesh):
    """Tops on a sphere of twice the radius over the cap (100 mm against 50): the fit is within 5 % of it, and
    `nose_radius()` returns the mesh's 50 mm under `initial`. At twice the radius the fit (100 mm) is flatter than the bound
    (83.5 mm), so under `current` the bound is used (Asha, 2026-11-10, option B); at 1.5 times (75 mm) the fit itself is."""
    b = make(coarse_sphere_mesh)
    pts, R2 = flatten_front(b)
    with_tops(b, pts)
    assert abs(b.nose_fit.fitted / R2 - 1.0) < 0.05 and not b.nose_fit.concave
    assert b.nose_fit.fitted > b.nose_fit.bound and b.nose_radius() == b.nose_fit.bound
    pts15, R15 = flatten_front(b, factor=1.5)
    with_tops(b, pts15)
    assert abs(b.nose_fit.fitted / R15 - 1.0) < 0.05 and b.nose_fit.fitted < b.nose_fit.bound
    assert b.nose_radius() == b.nose_fit.fitted
    b.size_feedback = "initial"
    assert b.nose_radius() == 0.05 and b.nose_fit.fitted > 0.07


def test_a_cap_flatter_than_the_bound_uses_the_bound_and_reports_both(coarse_sphere_mesh):
    """A cap flatter than the flat-face value (the 2R sphere's 100 mm against 1.67 times the 50 mm transverse radius, 83.5 mm)
    is used at the bound, a flat face heating like a sphere of 1.67 times its radius (Asha, 2026-11-10, option B), while the
    fitted radius is still reported above it; a sharper cap has no lower bound."""
    b = make(coarse_sphere_mesh)
    with_tops(b, flatten_front(b)[0])
    f = b.nose_fit
    assert f.fitted > f.bound and f.used == f.bound == b.nose_radius()
    s = b.skin_stats()
    assert s["fitted_nose_radius_mm"] > s["nose_bound_mm"] and s["nose_radius_mm"] == s["nose_bound_mm"] == f.bound * 1e3
    assert s["fitted_nose_radius_mm"] == f.fitted * 1e3 and s["nose_cap_concave"] == 0.0
    f = with_tops(b, flatten_front(b, factor=0.6)[0]).nose_fit
    assert f.fitted < 0.05 and f.used == f.fitted


def test_the_heating_and_the_surface_flow_both_receive_the_fitted_radius(coarse_sphere_mesh):
    """The radius each load is called with, not only what `nose_radius()` returns: the stagnation heating through
    `CoupledRun.loads_at` and the surface flow through `flow_forcing`."""
    from types import SimpleNamespace
    from reentry_model import coupled
    b = make(coarse_sphere_mesh, flow=True)
    with_tops(b, flatten_front(b, factor=1.5)[0])
    seen = {}
    ran = SimpleNamespace(sim=SimpleNamespace(aero_state=lambda t, r, v: "state"), body=b,
                          heating=SimpleNamespace(evaluate=lambda a, theta, T, R, T_mean=None: seen.setdefault("heating", R)))
    coupled.CoupledRun.loads_at(ran, 0.0, np.zeros(6))

    def flow_evaluate(state, theta, R, rho_l, T):
        seen["flow"] = R
        return SimpleNamespace(dp_ds=np.zeros(len(theta)), tau=np.zeros(len(theta)), deceleration=0.0)
    b.surface_flow = SimpleNamespace(evaluate=flow_evaluate)
    b.flow_forcing("state")
    assert seen["heating"] == seen["flow"] == b.nose_fit.fitted > 0.07 and seen["flow"] != 0.05 and b.nose_fit.fitted < b.nose_fit.bound


def test_a_cap_with_fewer_than_four_corners_falls_back_to_the_original_radius(coarse_sphere_mesh):
    """Three corners within 30 degrees of the flow cannot fix a sphere: no fit (not-a-number), the mesh's 50 mm used, the
    cap's three corners counted, the row says so and the bound is still reported."""
    b = make(coarse_sphere_mesh)
    theta = np.full(b.columns.n, 1.0)
    theta[:3] = 0.1
    with_tops(b, theta=theta)
    f = b.nose_fit
    assert f.points == 3 and np.isnan(f.fitted) and f.used == 0.05 == b.nose_radius() and not f.concave
    s = b.skin_stats()
    assert np.isnan(s["fitted_nose_radius_mm"]) and s["nose_radius_mm"] == 50.0 and s["nose_cap_points"] == 3
    assert s["nose_bound_mm"] > 0.0


def test_a_concave_cap_uses_the_flat_face_bound_and_says_so(coarse_sphere_mesh):
    """Tops dished toward the flow (the cap's edge nearer the flow than its middle): the fitted centre is on the flow side, no
    sphere of positive curvature describes it, so the radius used is the bound and the row records a concave cap."""
    b = make(coarse_sphere_mesh)
    R = 0.05
    c0, v = b.mesh.points.mean(axis=0), b.v_hat
    pts = b.surface.points.copy()
    X = pts - c0
    ax = X @ v
    rho = np.sqrt(np.maximum((X * X).sum(axis=1) - ax ** 2, 0.0))
    cap = (ax > 0.0) & (rho < 0.4 * R)
    pts[cap] += np.outer(rho[cap] ** 2 / R, v)                   # dished: twice the sphere's sag the other way, so the edge is nearer the flow than the middle
    theta = np.where(cap, 0.2, 1.0)
    with_tops(b, pts, theta=theta)
    f = b.nose_fit
    assert f.concave and f.used == f.bound == b.nose_radius() and np.isfinite(f.fitted) and f.points == cap.sum()
    assert b.skin_stats()["nose_cap_concave"] == 1.0


def test_a_pool_behind_the_body_does_not_move_the_bound(coarse_sphere_mesh):
    """Thicker skins on the leeward columns, up to the equator, with the angles recomputed by the real refresh: the windward
    columns are those whose base normals face the flow, which a pool lying on top cannot tilt, so the bound is unchanged. (Selecting
    by the tops' own angles let a 1 mm layer raise it from 83.83 to 85.49 mm and a 5 mm one to 91.90 mm: the step at the equator
    tilts leeward corners forward; review of 2026-11-10.)"""
    b = make(coarse_sphere_mesh)
    before = b.nose_fit.bound
    ax = (b.columns.points - b.mesh.points.mean(axis=0)) @ b.v_hat
    lee = b.columns.normal @ b.v_hat < 0.0
    assert lee.sum() > 100 and ax[lee].max() > 0.0                  # the layer reaches the equator and a little beyond
    for depth in (1e-3, 5e-3):
        b.skins.fill(b.skins.column_mass() + np.where(lee, b.material.rho * depth * b.skins.area, 0.0), 880.0)
        b._refresh_geometry()
        assert b.nose_fit.bound == before


def test_a_cap_coplanar_to_rounding_is_a_flat_face_and_reads_the_bound(coarse_sphere_mesh):
    """Cap corners set exactly in a plane have no sphere (the least-squares system is rank-deficient and its minimum-norm answer
    depends on where the plane sits): they are read as concave, no fitted radius, the bound used (review of 2026-11-10)."""
    b = make(coarse_sphere_mesh)
    c0, v = b.mesh.points.mean(axis=0), b.v_hat
    pts = b.surface.points.copy()
    X = pts - c0
    ax = X @ v
    rho = np.sqrt(np.maximum((X * X).sum(axis=1) - ax ** 2, 0.0))
    cap = (ax > 0.0) & (rho < 0.4 * 0.05)
    pts[cap] += np.outer(ax[cap].max() - ax[cap], v)                 # every cap corner on the plane of the nose
    with_tops(b, pts, theta=np.where(cap, 0.2, 1.0))
    f = b.nose_fit
    assert f.concave and np.isnan(f.fitted) and f.used == f.bound == b.nose_radius() and f.points == cap.sum()


def test_size_feedback_initial_changes_nothing_to_the_last_bit(coarse_sphere_mesh):
    """`size_feedback="initial"`, the class default a body saved before Task 5g flies on and the old `nose_radius` (the
    mesh's radius) are one flight: every pre-existing column of every row, and the final mass and energy, equal to the last
    bit. 3 MW/m2 on the front from 893 K for two 0.0625 s macro steps with the flow on, its loads computed through the body's
    own `nose_radius()` (a heating that falls as one over its square root and a pressure gradient as one over it), so that the radius is really read; a fourth arm under `current` flies a different flight, which shows
    that it is. The five nose columns are there, filled from the fit (sub-plan 22 Task 5g)."""
    legacy = skinbody.SkinBody.nose_radius
    runs = []
    for how in ("initial", "default", "legacy", "current"):
        np.random.seed(12345)
        if how == "legacy":
            skinbody.SkinBody.nose_radius = lambda self: float(self.mesh.params.get("radius_m", 0.05))
        try:
            kw = {"size_feedback": "initial"} if how == "initial" else {"size_feedback": "current"} if how == "current" else {}
            b = make(coarse_sphere_mesh, T0=893.0, flow=True, remove=False, flow_cycles=1, **kw)
            if how == "default":
                del vars(b)["size_feedback"]                     # as a state saved before Task 5g
                assert b.size_feedback == "initial"
            rows = []
            for k in range(2):
                R = b.nose_radius()
                loads = front_loads(b, 3e6 * np.sqrt(0.05 / R))
                wind = b.theta < 0.5 * np.pi
                push = skinbody.FlowForcing(np.where(wind, -5e2 / R * np.sin(b.theta), 0.0), np.where(wind, 50.0, 0.0),
                                            30.0 * b.v_hat)
                b.advance(0.0625 * (k + 1), 0.0625, loads, forcing=push)
                rows.append((b.step_parts, timed_out(b.skin_stats())))
        finally:
            skinbody.SkinBody.nose_radius = legacy
        runs.append((b, rows))
    assert all(k in rows[0][1] for _, rows in runs for k in skinbody.NOSE_COLUMNS)
    old = lambda rows: [(p, {k: v for k, v in r.items() if k not in skinbody.NOSE_COLUMNS}) for p, r in rows]
    assert same_rows(old(runs[0][1]), old(runs[1][1])) and same_rows(old(runs[0][1]), old(runs[2][1]))
    assert not same_rows(old(runs[0][1]), old(runs[3][1]))               # the radius is read: the fitted one flies another flight
    assert all(r["nose_radius_mm"] == 50.0 for _, r in runs[0][1]) and runs[0][0].handed_mass > 0.0
    assert runs[0][0].mass() == runs[1][0].mass() == runs[2][0].mass()
    assert runs[0][0].energy() == runs[1][0].energy() == runs[2][0].energy()
    assert skinbody.SkinBody.size_feedback == "initial" and all(np.isfinite(r["fitted_nose_radius_mm"]) for _, r in runs[0][1])
    with pytest.raises(ValueError, match="size_feedback"):
        make(coarse_sphere_mesh, size_feedback="shrinking")


def test_a_body_saved_without_size_feedback_flies_on_as_initial(coarse_sphere_mesh):
    """The states saved before Task 5g have neither `size_feedback` nor a fit on the instance: the class's default gives them
    the mesh's radius and they advance as before, the fit made at the first ask; a resume can give them `current`, which
    reads the fit of the top surface as it stands."""
    b = make(coarse_sphere_mesh, T0=700.0)
    b.emissivity = 0.0
    del vars(b)["size_feedback"], vars(b)["_nose_raw"]
    assert "size_feedback" not in vars(b) and b.size_feedback == "initial" and b._nose_raw is None
    zero = heating.HeatingResult(np.zeros(b.columns.n), 0.0, 0.0, 0.0, 0.0)
    b.advance(0.05, 0.05, zero)
    assert b.nose_radius() == 0.05 and b.skin_stats()["nose_radius_mm"] == 50.0 and np.isfinite(b.nose_fit.fitted)
    b.size_feedback = "current"
    assert b.nose_radius() == b.nose_fit.fitted != 0.05


def test_with_the_hand_down_off_a_release_s_momentum_is_not_counted_twice(coarse_sphere_mesh):
    """Found by sub-plan 22 Task 5h on 2026-11-10 (a defect of Task 5e): `_move_base` added `momentum_out` to the momentum lost
    as "what the handed-down cells held" even where no hand-down ran, when it still held what the last take (the spray
    release's) took, which the release book had already counted. The command line's spraying flight missed its momentum books
    by 4.4e-6 of the impulses with the hand-down off, and the 50 mm comparison flight by up to 1.4e-5 after spraying began.
    The prescribed-release body of sub-plan 21 Task 7, with the hand-down off, closes its momentum books to round-off."""
    b = make(coarse_sphere_mesh, flow=True, remove=False, flow_cycles=1, release=skin.prescribed_release(2.0), hand_down=False)
    for k in range(2):
        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), forcing=newtonian(b))
    assert b.released_mass > 0.0 and np.abs(b.droplet_momentum).max() > 0.0 and b.handed_mass == 0.0
    assert b.skin_stats()["momentum_residual"] < 1e-12


# -- spraying in the slip band (sub-plan 22 Task 5i; the plan of 2026-11-06, Asha 2026-11-10) -------------------------------
def band_plan(sk, col, flow, band, slurry=False, theta=0.6, gates="2026-11-10", nose_radius=None):
    """`girin_plan` on columns at one angle `theta` with a spray model under `gates` and the band `band`."""
    from reentry_model import spray
    from test_reentry_model_spray import FakeState
    model = spray.SprayModel(sk.material.liquid, spray_gates=gates, slip_band=band)
    return skinbody.girin_plan(sk, flow, FakeState(), model, np.full(sk.n, theta), col.area, col.edges, 0.05, 0.0125,
                               slurry=slurry, nose_radius=nose_radius)


TAU = 200.0                    # Pa, the bridged shear on the band cases


def band_case(kn_body=0.04, kn_local=0.05):
    """Eight columns (four kinds, twice) 20 mm apart at one angle of 0.6 rad, on the merged branch (Couette closure, no Girin
    delta_m), the bridged shear TAU: 50 um of liquid over 500 um of slurry (at 902 K) over rigid metal; 50 um of liquid on the
    rigid floor; exposed slurry in the bridge (907.9 K, within a kelvin of the liquidus, the only slurry that sprays); 15 um
    of liquid on the rigid floor, one cell. The edge speed is 1 500 m/s."""
    from test_reentry_model_skin import LIQUID, RIGID, SLURRY, stacked, strip_grid
    from test_reentry_model_spray import slip_flow
    col, _ = strip_grid(4, 2, h=20e-3)
    n = col.n
    temps = [[LIQUID] * 4 + [SLURRY] * 50 + [RIGID] * 10, [LIQUID] * 4 + [RIGID] * 10, [907.9] * 50 + [RIGID] * 10,
             [LIQUID] + [RIGID] * 10] * 2
    sk = stacked(col, temps)
    flow, _ = slip_flow(kn_body, (kn_local,) * n, (0.6,) * n, tau=TAU)
    flow.u_eff = flow.u_e = np.full(n, 1500.0)
    flow.G = np.zeros(n)
    return sk, col, flow


def test_band_off_leaves_a_slipping_film_held():
    """Without the band the merged branch has no delta_m: a liquid on slurry is held from every shear mode (the element
    model's reading, 2026-10-29), and the liquid on the rigid floor, with the gate of 2026-11-10 on, releases nothing."""
    sk, col, flow = band_case()
    plan = band_plan(sk, col, flow, band=False)
    assert plan.rate.sum() == 0.0 and plan.rate_slurry.sum() == 0.0
    assert plan.on_slurry.tolist() == [True, False, True, False] * 2 and (plan.result.branch == -1).all()


def test_band_on_sprays_the_liquid_and_reads_the_flowing_depth():
    """Band on (body Knudsen number 0.04, wall 0.05 everywhere, 200 Pa of shear): the momentum-balance depth is finite on every
    windward patch (54 um), a liquid on slurry is no longer held, its flowing depth (0.6 mm) is deeper than that depth so
    the film is on Girin's thick branch, and 50 um of liquid on the rigid floor is shallower than it: the thin branch.
    Both release; the release is capped by the liquid (slurry spraying off), so the exposed slurry gives nothing."""
    from reentry_model import spray
    sk, col, flow = band_case()
    plan = band_plan(sk, col, flow, band=True)
    delta = plan.result.delta_m
    assert np.isfinite(delta).all() and (40e-6 < delta).all() and (delta < 120e-6).all()
    assert plan.result.branch[:2].tolist() == [spray.BRANCH_THICK, spray.BRANCH_THIN]
    assert plan.rate[0] > 0.0 and plan.rate[1] > 0.0 and plan.rate[2] == 0.0      # no liquid on the exposed slurry
    assert not plan.on_slurry[0] or plan.rate[0] > 0.0
    np.testing.assert_array_equal(plan.cap_depth, sk.depth_above(sk.material.T_liquidus))
    np.testing.assert_array_equal(delta, spray.spray_layer(flow, sk.material.liquid, True, np.full(8, 0.6), 0.05))


def test_band_release_is_capped_by_the_liquid_unless_the_slurry_is_on():
    """Decision 4 of 2026-11-06 and item 3 of 2026-11-10: with --slurry-spray on the exposed slurry in the bridge sprays in the
    band by its own mode (`slurry_modes`, its depth the momentum balance's with its own viscosity) and the cap is the flowing
    depth; with it off the cap stays the liquid's depth and the slurry sprays nothing. The liquid's mode is the same."""
    sk, col, flow = band_case()
    off, on = band_plan(sk, col, flow, band=True), band_plan(sk, col, flow, band=True, slurry=True)
    assert off.rate_slurry.sum() == 0.0
    np.testing.assert_array_equal(on.cap_depth, sk.depth_above(sk.material.T_rigid))
    assert on.rate_slurry[2] > 0.0 and np.isfinite(on.slurry["delta_m"][2]) and on.slurry["delta_m"][2] < on.slurry["depth"][2]
    np.testing.assert_array_equal(on.rate[:2], off.rate[:2])


def test_the_band_s_thin_branch_keeps_the_viscous_floor():
    """Decision 6 of 2026-11-06: the viscous film's floor already applies on the liquid thin branch, and the band changes no
    gate there, at an edge speed of 1 500 m/s: 50 um of liquid on the rigid floor sprays (gamma tau 0.47, below one), 15 um is
    damped by the floor (gamma tau 1.19, above one) and releases nothing; both record the floor."""
    from reentry_model import spray
    sk, col, flow = band_case()
    res = band_plan(sk, col, flow, band=True).result
    assert res.branch[[1, 3]].tolist() == [spray.BRANCH_THIN, spray.BRANCH_THIN]
    assert np.isfinite(res.h_floor[[1, 3]]).all() and res.gamma_tau[1] < 1.0 and res.mdot[1] > 0.0
    assert res.mdot[3] == 0.0 and res.gamma_tau[3] > 1.0


def test_the_band_above_a_body_knudsen_number_of_a_tenth_or_a_wall_knudsen_number_of_a_tenth_does_nothing():
    sk, col, flow = band_case(kn_body=0.2)
    assert band_plan(sk, col, flow, band=True).rate.sum() == 0.0
    sk, col, flow = band_case(kn_body=0.04, kn_local=0.12)             # Asha, 2026-11-10: the wall's own number too
    assert band_plan(sk, col, flow, band=True).rate.sum() == 0.0
    sk, col, flow = band_case()
    assert band_plan(sk, col, flow, band=True).rate.sum() > 0.0         # the same columns inside both limits


def test_the_band_depth_reads_the_nose_radius_it_is_given():
    """The momentum balance's r is the nose radius times sin(theta): half the radius, half the r, and a depth larger by
    2^(1/3), the plan's delta_m^3 = 3 mu^2 I / (rho tau^2 r) with the integral I unchanged."""
    sk, col, flow = band_case()
    a = band_plan(sk, col, flow, band=True, nose_radius=0.05).result.delta_m
    b = band_plan(sk, col, flow, band=True, nose_radius=0.025).result.delta_m
    np.testing.assert_allclose(b / a, 2.0 ** (1.0 / 3.0), rtol=1e-12)
    np.testing.assert_array_equal(band_plan(sk, col, flow, band=True).result.delta_m, a)      # the default: the transverse radius, 0.05


def test_the_thin_gate_leaves_a_slipping_liquid_on_the_floor_unsprayed_with_the_band_off_and_before_sprays_it():
    """Task 5i's gate on the skins: with the band off the 50 um film on the floor has no conjugate depth and sprays nothing under
    2026-11-10; under `before` the thin branch sprays it as it did."""
    sk, col, flow = band_case()
    assert band_plan(sk, col, flow, band=False).rate[1] == 0.0
    from reentry_model import spray
    old = band_plan(sk, col, flow, band=False, gates="before")
    assert old.result.branch[1] == spray.BRANCH_THIN and old.rate[1] > 0.0                  # as before this task


def test_the_skins_hand_girin_plan_their_nose_radius(coarse_sphere_mesh, monkeypatch):
    """`SkinBody.advance` gives `girin_plan` the nose radius the loads read (`nose_radius()`, Task 5g's fit), for the band's
    r = R sin(theta), where the transverse radius is the cap's size limit."""
    pytest.importorskip("cantera")
    from test_reentry_model_surface_flow import aero_state_at
    from reentry_model import spray
    mat = material.Material.from_drama_json("AA7075_scheil")
    b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=spray.SprayModel(mat.liquid))
    seen = []
    girin = skinbody.girin_plan

    def planned(*a, **k):
        seen.append((k.get("nose_radius"), b.nose_radius()))
        return girin(*a, **k)
    monkeypatch.setattr(skinbody, "girin_plan", planned)
    b.advance(0.25, 0.25, front_loads(b, 3e6), state=aero_state_at(93.0))
    assert len(seen) == 1 and seen[0][0] is not None and seen[0][0] == seen[0][1]


# -- the leeward test by the body's centroid (sub-plan 22 Task 5j; Asha, 2026-11-10) ----------------------------------------
def centroid_body(m, **kw):
    return make(m, leeward_test="centroid", **kw)


def put_behind(b, back=0.04):
    """A corner facing the flow (its base normal within 60 degrees of the flight direction, its top angle 0.5 rad) whose top is
    moved `back` m behind the plane through the mass centroid at right angles to the flight direction, as a pool standing
    off the back leaves one: returns its index. The tops, angles and fits are refreshed as `_refresh_geometry` does."""
    c = b.mass_centroid()
    front = np.flatnonzero(b.columns.normal @ b.v_hat > 0.5)
    k = int(front[np.argmin((b.surface.points[front] - c) @ b.v_hat)])      # the one nearest the plane
    pts, theta = b.surface.points.copy(), b.theta.copy()
    pts[k] -= ((pts[k] - c) @ b.v_hat + back) * b.v_hat
    theta[k] = 0.5
    with_tops(b, pts, theta)
    return k


def test_the_mass_centroid_of_an_untouched_sphere_is_its_centre(coarse_sphere_mesh):
    """The mass-weighted centre of the 3D model's elements and the skins' cells (sub-plan 22 Task 5j) of the 4 mm test sphere
    sits within a fifth of a millimetre of its centre, the origin, on a 50 mm radius, and `centroid_offset_mm` (its position
    along the flight direction relative to the 3D mesh's mean point, positive forward) is that small too."""
    b = centroid_body(coarse_sphere_mesh)
    c = b.mass_centroid()
    assert c.shape == (3,) and np.linalg.norm(c) < 2e-4
    assert abs(b.skin_stats()["centroid_offset_mm"]) < 0.2 and b.skin_stats()["leeward_corners"] < 10   # a few equator corners


def test_a_pool_at_the_back_moves_the_centroid_backward_by_its_mass_moment(coarse_sphere_mesh):
    """A 5 mm pool on the leeward columns (base normals facing away from the flow): the centroid moves backward, and to the
    place the mass moment gives, the 3D model's moment (its elements' volumes times the density, at their centroids) plus each
    column's mass at the middle of its skin (equal cells of one density: base plus normal times half the thickness)."""
    b = centroid_body(coarse_sphere_mesh)
    front = b.mass_centroid()
    lee = b.columns.normal @ b.v_hat < 0.0
    b.skins.fill(b.skins.column_mass() + np.where(lee, b.material.rho * 5e-3 * b.skins.area, 0.0), 880.0)
    b._refresh_geometry()
    c = b.mass_centroid()
    assert (c - front) @ b.v_hat < -1e-4                                    # backward: against the flight direction
    mm = b.mesh
    live = mm.active
    vol, at = mm.element_volumes()[live], mm.points[mm.tets[live]].mean(axis=1)
    m3 = b.material.rho * vol
    mid = b.columns.points + b.columns.normal * (0.5 * b.skins.thickness())[:, None]
    m = b.skins.column_mass()
    expected = ((m3[:, None] * at).sum(axis=0) + (m[:, None] * mid).sum(axis=0)) / (m3.sum() + m.sum())
    np.testing.assert_allclose(c, expected, atol=1e-9)
    assert b.skin_stats()["centroid_offset_mm"] == pytest.approx(1e3 * (c - b.mesh.points.mean(axis=0)) @ b.v_hat, abs=1e-9)


def test_the_effective_angle_leaves_a_corner_ahead_of_the_plane_as_it_is_and_turns_one_behind_it_leeward(coarse_sphere_mesh):
    """Under `centroid` a corner is windward only if it is ahead of the plane through the centroid AND its normal faces the flow;
    the effective angle is theta where it is and the larger of theta and 90 degrees plus a microradian where not. `facing`
    returns `theta` itself. The count of corners that face the flow from behind the plane is the history's `leeward_corners`,
    with either setting."""
    b = centroid_body(coarse_sphere_mesh)
    moved = b.load_angles() != b.theta                                       # only an equator corner or two, a hair behind the plane
    assert moved.sum() < 10 and (b.theta[moved] > 1.5).all()
    k = put_behind(b)
    eff = b.load_angles()
    other = np.arange(b.theta.size) != k
    assert b.theta[k] == 0.5 and eff[k] == 0.5 * np.pi + 1e-6
    assert (eff[other & (b.theta < 1.5)] == b.theta[other & (b.theta < 1.5)]).all()      # every other windward corner is as it was
    assert (eff[other] >= b.theta[other]).all()                                         # an equator corner a hair behind the plane is only turned leeward
    away = int([i for i in np.flatnonzero(b.columns.normal @ b.v_hat > 0.5) if i != k][0])
    theta = b.theta.copy()
    theta[away] = 2.0                                                    # ahead of the plane but facing away
    with_tops(b, theta=theta)
    assert b.load_angles()[away] == 2.0 and b.load_angles()[k] == 0.5 * np.pi + 1e-6
    assert b.skin_stats()["leeward_corners"] >= 1
    f = make(coarse_sphere_mesh, leeward_test="facing")
    put_behind(f)
    assert f.load_angles() is f.theta and f.skin_stats()["leeward_corners"] >= 1


def test_a_corner_behind_the_centroid_that_faces_the_flow_gets_no_shear_no_spray_and_no_rayleigh_taylor(coarse_sphere_mesh,
                                                                                                         monkeypatch):
    """The wiring through `SkinBody.advance`: the surface flow and the spray plan both receive the effective angle, so the
    corner put behind the plane (theta 0.5 rad, its normal facing the flow) has no shear and a non-windward patch; the spray
    model is told it is not windward and its Rayleigh-Taylor deceleration is zero; and no droplet leaves it. A corner ahead of
    the plane that faces away gets no shear or spray either."""
    pytest.importorskip("cantera")
    from test_reentry_model_surface_flow import aero_state_at
    from reentry_model import spray, surface_flow
    mat = material.Material.from_drama_json("AA7075_scheil")
    seen = {}

    class Model(spray.SprayModel):
        def evaluate(self, flow, state, b, delta_m, v_s, windward, dt, areas, m_f, *a, **k):
            seen["windward"], seen["w_n"] = windward.copy(), np.array(k["deceleration_n"])
            return super().evaluate(flow, state, b, delta_m, v_s, windward, dt, areas, m_f, *a, **k)

    class Flow(surface_flow.SurfaceFlow):
        def evaluate(self, state, theta, *a, **k):
            seen["flow_theta"] = np.array(theta)
            seen["result"] = super().evaluate(state, theta, *a, **k)
            return seen["result"]
    b = centroid_body(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=Model(mat.liquid),
                      surface_flow=Flow())
    k = put_behind(b)
    away = int([i for i in np.flatnonzero(b.columns.normal @ b.v_hat > 0.5) if i != k][3])
    theta = b.theta.copy()
    theta[away] = 2.0
    with_tops(b, theta=theta)
    eff = b.load_angles().copy()
    plans, girin = [], skinbody.girin_plan
    monkeypatch.setattr(skinbody, "girin_plan", lambda *a, **kw: plans.append((np.array(a[4]), girin(*a, **kw))) or plans[-1][1])
    b.advance(0.25, 0.25, front_loads(b, 3e6), state=aero_state_at(93.0))
    np.testing.assert_array_equal(seen["flow_theta"], eff)
    np.testing.assert_array_equal(plans[0][0], eff)
    assert eff[k] > 0.5 * np.pi and not seen["windward"][k] and seen["w_n"][k] == 0.0 and seen["w_n"][away] == 0.0
    res = seen["result"]
    for j in (k, away):
        assert res.tau[j] == 0.0 and res.dp_ds[j] == 0.0 and plans[0][1].mass[j] == 0.0 and plans[0][1].result.branch[j] == -1


def test_the_heating_reads_the_effective_angle(coarse_sphere_mesh):
    """`CoupledRun.loads_at` hands the heating model the effective angles: under `centroid` the corner behind the plane is past
    90 degrees, under `facing` the heating gets `theta` as it always did (sub-plan 22 Task 5j)."""
    from types import SimpleNamespace
    from reentry_model import coupled
    got = {}
    for how in ("centroid", "facing"):
        b = make(coarse_sphere_mesh, leeward_test=how)
        k = put_behind(b)
        ran = SimpleNamespace(sim=SimpleNamespace(aero_state=lambda t, r, v: "state"), body=b,
                              heating=SimpleNamespace(evaluate=lambda a, theta, T, R, T_mean=None: got.setdefault(how, np.array(theta))))
        coupled.CoupledRun.loads_at(ran, 0.0, np.zeros(6))
        if how == "centroid":
            assert got[how][k] > 0.5 * np.pi and got[how][k] != b.theta[k]
        else:
            np.testing.assert_array_equal(got[how], b.theta)


def test_the_leeward_test_facing_changes_nothing_to_the_last_bit(coarse_sphere_mesh):
    """A body on the class default `facing` (what a body saved before Task 5j flies on) and one whose angles are read as `theta` (the
    method replaced by the old attribute) are one flight: every column of every row but the two new ones, the surface flow's
    shear and pressure gradient, the spray plan's depths, the mass released, and the final mass and energy, equal to the last
    bit. The flight sprays with the surface flow and the spray plan reading the angles (the gas state of the 100 mm flight at
    93 s, 3e6 W/m2 on the front from 890 K, the slip band on, one 0.25 s macro step, as the books test with the band on), so every consumer of the effective
    angle is called (review of Task 5j, Minor 1)."""
    pytest.importorskip("cantera")
    from test_reentry_model_surface_flow import aero_state_at
    from reentry_model import spray
    mat = material.Material.from_drama_json("AA7075_scheil")
    legacy = skinbody.SkinBody.load_angles
    runs = []
    for how in ("default", "legacy"):
        np.random.seed(12345)
        if how == "legacy":
            skinbody.SkinBody.load_angles = lambda self: self.theta
        try:
            kw = {"leeward_test": "facing"} if how == "legacy" else {}
            model = spray.SprayModel(mat.liquid, spray_gates="2026-11-10", slip_band=True)
            b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=model, **kw)
            if how == "default":
                del vars(b)["leeward_test"]                 # as a state saved before Task 5j
                assert b.leeward_test == "facing"
            state = aero_state_at(93.0)
            state.kn = 0.04
            rows = []
            for k in range(1):
                b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), state=state)
                rows.append((b.step_parts, timed_out(b.skin_stats())))
        finally:
            skinbody.SkinBody.load_angles = legacy
        runs.append((b, rows))
    assert all(k in rows[0][1] for _, rows in runs for k in skinbody.LEEWARD_COLUMNS)
    assert runs[0][0].last_flow.in_slip_band and np.isfinite(runs[0][0].spray_plan.result.delta_m).any()   # the consumers really ran
    old = lambda rows: [(p, {k: v for k, v in r.items() if k not in skinbody.LEEWARD_COLUMNS}) for p, r in rows]
    assert same_rows(old(runs[0][1]), old(runs[1][1]))
    for b, _ in runs[1:]:
        for name in ("tau", "dp_ds"):
            np.testing.assert_array_equal(getattr(b.last_flow, name), getattr(runs[0][0].last_flow, name))
        np.testing.assert_array_equal(b.spray_plan.result.delta_m, runs[0][0].spray_plan.result.delta_m)
        assert b.released_mass == runs[0][0].released_mass and b.mass() == runs[0][0].mass() and b.energy() == runs[0][0].energy()
    # a third arm with teeth (final review, Minor 3; Task 5k, made cheap in its fix round 1): the two arms above are one code path, so
    # they cannot fail on the angle plumbing. Under `centroid` a corner facing the flow from behind the plane gets no shear and no
    # pressure gradient from the surface flow; under `facing` the same corner does. The corner is put there synthetically, 4 cm behind
    # the plane (not the ruling's case of a hair behind it: 3 of 2254 on the 4 mm sphere), and the flow is evaluated directly
    # (`flow_forcing`, the routine `advance` calls), with no macro step run and no spray model
    forcing = {}
    for how in ("centroid", "facing"):
        b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, leeward_test=how)
        k = put_behind(b)
        forcing[how] = b.flow_forcing(aero_state_at(93.0))
    assert forcing["centroid"].tau[k] == 0.0 and forcing["centroid"].dp_ds[k] == 0.0
    assert forcing["facing"].tau[k] != 0.0 or forcing["facing"].dp_ds[k] != 0.0
    assert (forcing["centroid"].tau != 0.0).any()                       # the flow is not zero everywhere for another reason
    assert skinbody.SkinBody.leeward_test == "facing" and skinbody.LEEWARD_TESTS == ("facing", "centroid")
    with pytest.raises(ValueError, match="leeward_test"):
        make(coarse_sphere_mesh, leeward_test="behind")


def test_the_fit_takes_only_windward_corners_within_the_cap(coarse_sphere_mesh):
    """Under `centroid` the sphere's cap is the corners windward by the new test within 30 degrees of the flow, and the
    windward transverse radius for the bound takes the windward corners likewise: ten cap corners put behind the plane drop out
    of the fit (the count falls by ten) and a faraway corner behind it, with a base normal facing the flow, leaves the bound as
    it was; under `facing` Task 5g's selection is unchanged (the count stays, the faraway corner raises the bound)."""
    for how in ("centroid", "facing"):
        b = make(coarse_sphere_mesh, leeward_test=how)
        n0, bound0 = b.nose_fit.points, b.nose_fit.bound
        c = b.mass_centroid()
        cap = np.flatnonzero(b.theta <= np.radians(30.0))
        far = int(np.flatnonzero((b.columns.normal @ b.v_hat > 0.5) & (b.theta > np.radians(40.0)))[-1])
        pts = b.surface.points.copy()
        for k in cap[:10]:
            pts[k] -= ((pts[k] - c) @ b.v_hat + 0.04) * b.v_hat
        pts[far] = c - 0.04 * b.v_hat + 0.3 * np.array([0.0, 1.0, 0.0])
        theta = b.theta.copy()
        with_tops(b, pts, theta)
        f = b.nose_fit
        if how == "centroid":
            assert f.points == n0 - 10 and f.bound == pytest.approx(bound0, rel=0.02)
        else:
            assert f.points == n0 and f.bound > 1.5 * bound0


def test_the_history_has_the_centroid_offset_and_the_leeward_corners_with_either_setting(coarse_sphere_mesh):
    """Two new skin columns, `centroid_offset_mm` and `leeward_corners`, with either setting; the corners that face the flow
    from behind the plane through the centroid are counted whatever the setting decides to do with them."""
    for how in ("centroid", "facing"):
        b = make(coarse_sphere_mesh, leeward_test=how)
        s0 = b.skin_stats()
        assert s0["leeward_corners"] < 10 and abs(s0["centroid_offset_mm"]) < 0.2
        put_behind(b)
        s = b.skin_stats()
        assert s["leeward_corners"] == s0["leeward_corners"] + 1 and all(k in s for k in skinbody.LEEWARD_COLUMNS)


# -- the folded-in fixes of the review of Task 5i (sub-plan 22 Task 5j) -----------------------------------------------------
def test_with_the_band_on_the_depth_ratios_read_the_depth_the_plan_sprayed_on(coarse_sphere_mesh):
    """Review of Task 5i, Minor 3: with the band on and a plan for the step, `b_over_delta_m_median` is the flowing depth over
    the plan's own `delta_m` (here 100 um on every column, where a recomputation from the start-of-step flow with the
    end-of-step angles and radius would give something else), not a depth recomputed in `_depth_ratios`."""
    from types import SimpleNamespace
    from reentry_model import spray
    mat = material.Material.from_drama_json("AA7075_scheil")
    b = make(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=spray.SprayModel(mat.liquid, slip_band=True))
    b.last_flow = object()                       # never read: the plan has the depth
    b.spray_plan = SimpleNamespace(result=SimpleNamespace(delta_m=np.full(b.skins.n, 1e-4)))
    flowing = np.full(b.skins.n, 5e-5)           # a flowing depth on every column (880 K metal has none to measure)
    wet = b.windward_corners() & (flowing > 0.0)
    assert wet.sum() > 50
    got = b._depth_ratios(flowing)["b_over_delta_m_median"]
    assert got == pytest.approx(float(np.median(flowing[wet] / 1e-4)), rel=1e-14)


def test_a_spraying_body_in_the_slip_band_keeps_the_three_books(coarse_sphere_mesh):
    """Review of Task 5i, item 4: melting skins advanced with the slip band on (the gas state of the 100 mm flight at 93 s with
    its body Knudsen number set to 0.04, inside the band, so the spraying reads the momentum balance's depth) keep the mass
    (the mass left and released is the starting mass to 1e-12 of it), the energy (the balance's residual below 1e-8 of the heat
    absorbed) and the momentum (the books' residual below 1e-12 of the impulses), with the leeward test by the centroid."""
    pytest.importorskip("cantera")
    from test_reentry_model_surface_flow import aero_state_at
    from reentry_model import spray
    mat = material.Material.from_drama_json("AA7075_scheil")
    model = spray.SprayModel(mat.liquid, spray_gates="2026-11-10", slip_band=True)
    b = centroid_body(coarse_sphere_mesh, T0=890.0, flow=True, remove=False, flow_cycles=1, spray=model)
    M0, state = b.mass(), aero_state_at(93.0)
    state.kn = 0.04
    for k in range(2):
        b.advance(0.25 * (k + 1), 0.25, front_loads(b, 3e6), state=state)
    assert b.last_flow.in_slip_band and np.isfinite(b.spray_plan.result.delta_m).any()      # the band's depth is in use
    s = b.skin_stats()
    assert b.mass() + b.released_mass == pytest.approx(M0, rel=1e-12)
    assert abs(b.energy_balance_residual()) < 1e-8 and s["momentum_residual"] < 1e-12 and b.impulse_size > 0.0
