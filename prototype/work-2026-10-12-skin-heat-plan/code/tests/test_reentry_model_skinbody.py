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
    """6 MW/m2 on the front for 1.5 s from 880 K: liquid forms and leaves, columns draw, the base recedes; the mass
    removed plus the mass left is the mass at start, and the energy books close to 1e-8 of the heat absorbed. Under the
    damped repetition with the move limit (2026-10-12) each 0.5 s step runs as 5 to 7 parts and 169 g is removed, near
    the 190 g the heat could melt from 880 K (measured; 8 g under the one-repeat rule, which left the heat parked)."""
    b = make(coarse_sphere_mesh)
    M0 = b.mass()
    for k in range(3):
        b.advance(0.5 * (k + 1), 0.5, front_loads(b, 6e6))
    assert b.removed_mass > 0.0 and b.drawn_mass > 0.0
    assert b.mass() + b.removed_mass == pytest.approx(M0, rel=1e-13)
    assert abs(b.energy_balance_residual()) < 1e-8                  # 4.7e-9 measured (2026-10-12)
    stats = b.skin_stats()
    assert stats["base_recession_max_mm"] > 0.5                       # 10.2 mm at the nose (measured, 2026-10-12)
    assert b.surface_temperature().max() <= b.material.T_feed + 1.0   # liquid left as it formed
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

    # the default cap, a quarter of the mean surface edge (about 1.0 mm of 4.0 mm): every hot corner needs more than that
    cap = skinbody.ABSORB_MAX * b.columns.edge_length.mean()
    assert np.all(raw[hot] > cap)
    depth = b._absorb_depth()
    assert np.all(depth[~hot] == 0.0)                                      # below T_rigid: nothing to absorb
    assert depth[hot] == pytest.approx(cap, rel=1e-12) and b.absorb_capped == hot.sum()

    # a cap too high to bind for the corners with a real gradient: their depth is the excess over the gradient
    monkeypatch.setattr(skinbody, "ABSORB_MAX", 100.0)
    cap = skinbody.ABSORB_MAX * b.columns.edge_length.mean()                # 0.4 m
    b.absorb_capped = 0
    depth = b._absorb_depth()
    below = hot & (raw < cap)
    capped = hot & ~below
    assert below.sum() > 100 and capped.sum() > 100 and (capped & ~usable).sum() > 100     # the last: gradient under the floor
    assert np.all(depth[~hot] == 0.0)
    assert depth[below] == pytest.approx(raw[below], rel=1e-9) and np.all(depth[below] > 0.0)
    assert depth[capped] == pytest.approx(cap, rel=1e-12)                  # a gradient under the floor absorbs up to the cap
    assert b.absorb_capped == capped.sum()


# -- what triggers the interface repeat (Task 7, measured in Task 10) --------------------------------------------------
def test_every_repeat_rule_keeps_the_books_and_records_what_triggered_it(coarse_sphere_mesh):
    """Each undamped repeat rule (measured in Task 10; the damped rule is the default) leaves the books exact, and the
    first pass's mismatch is recorded whatever the rule. 2 MW/m2 from 880 K on 4 mm elements: without a repeat, 3 MW/m2
    already runs away -- the 3D base overshoots T_feed and every sub-step removes and redraws a whole skin (measured
    while planning)."""
    results = {}
    for rule in ("column", "total", "never"):
        b = make(coarse_sphere_mesh, repeat_rule=rule)
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
        b = make(coarse_sphere_mesh, repeat_rule=rule)
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
    the parts give what the same shorter steps give when taken one by one with the same heating and the same share of
    the absorption cap (ABSORB_MAX is per macro step)."""
    b = make(coarse_sphere_mesh)
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
    monkeypatch.setattr(skinbody, "ABSORB_MAX", skinbody.ABSORB_MAX / n)  # each with a part's share of the cap
    c = make(coarse_sphere_mesh)
    for k in range(n):
        c.advance(0.5 * (k + 1) / n, 0.5 / n, loads)
    assert c.removed_mass == pytest.approx(b.removed_mass, rel=1e-6)
    assert c.energy() == pytest.approx(b.energy(), rel=1e-9)
    assert np.abs(c.moving.mesh.points[c.moving.corners] - b.moving.mesh.points[b.moving.corners]).max() < 1e-9
