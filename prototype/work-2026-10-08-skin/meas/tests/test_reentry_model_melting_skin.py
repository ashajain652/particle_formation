"""MeltingBody with skins (Step 3 sub-plan 18): identical to the element model until the first skin; skins at the
solidus drawing from their elements; exact books; deaths handing skins on; the interface repeat; the edge cases of the
plan's Review Focus."""
import numpy as np
import pytest

from reentry_model import body, heating, material, mesh, skin, thermal
from test_reentry_model_coupled import MASS_100MM, simulator

pytest.importorskip("cantera")

DEV = skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"], substep=0.05)
# The 3D step solved directly: the default AMG-preconditioned CG stops at a relative residual of 1e-10, which leaves
# 1e-8 to 1e-7 of a 0.5 s step's heat unbooked on a body held in the mushy range (measured at 880-900 K with and without
# skins) -- more than assert_books allows, and enough to hide an error in the skins' books. Solved directly, the same
# steps close to about 1e-13, so the energy check sees the skins' books.
EXACT = dict(linear_solver="direct")


def make_body(the_mesh, surface_model="skin", name="AA7075_scheil", solver="skfem", solver_options=None, **settings):
    mat = material.Material.from_drama_json(name)
    m = mesh.VolumeMesh(the_mesh.points, the_mesh.tets, dict(the_mesh.params), the_mesh.element_layer.copy())
    return body.MeltingBody(m, mat, thermal.thermal_solver(solver, **(solver_options or {})), MASS_100MM,
                            settings=body.MeltSettings(surface_model=surface_model, skin=DEV, **settings))


def uniform(b, q):
    return heating.HeatingResult(np.full(b.surface.n_patches, q), q, 0.0, 0.0, 0.0)


def start_at(b, T):
    """A body held at T throughout, its books reckoned from here."""
    b.solver.set_temperature(T)
    b.energy0 = b.energy()


def run(b, q, steps, t0=0.0):
    t = t0
    for _ in range(steps):
        t += 0.5
        b.advance(t, 0.5, uniform(b, q))
    return t


def assert_books(b):
    assert b.mass(0.0) == pytest.approx(b.mass0 - b.removed_mass, rel=1e-12)
    assert abs(b.energy_balance_residual()) < 1e-8


def test_settings_validate_the_surface_model():
    with pytest.raises(ValueError):
        body.MeltSettings(surface_model="magic")
    with pytest.raises(ValueError):
        body.MeltSettings(surface_model="skin", removal="instant")
    assert body.MeltSettings().surface_model == "elements"


def test_identical_to_the_element_model_until_the_first_skin(coarse_sphere_mesh):
    """Seeded alike before every step (the AMG set-up draws on numpy's generator, as the model's --seed handles)."""
    be, bs = make_body(coarse_sphere_mesh, "elements"), make_body(coarse_sphere_mesh, "skin")
    assert be.skins is None and bs.skins is not None and bs.skins.n == 0
    t, compared = 0.0, 0
    while t < 120.0:
        t += 0.5
        np.random.seed(int(2 * t))
        be.advance(t, 0.5, uniform(be, 4e5))
        np.random.seed(int(2 * t))
        bs.advance(t, 0.5, uniform(bs, 4e5))
        if bs.skins.n:
            break
        assert np.array_equal(be.solver.temperature(), bs.solver.temperature())
        compared += 1
    assert compared >= 3 and bs.skins_created > 0


def test_skins_draw_from_their_elements_and_the_books_close(coarse_sphere_mesh):
    b = make_body(coarse_sphere_mesh)
    start_at(b, 740.0)
    run(b, 1.5e6, 12)
    assert b.skins.n > 0 and b.skin_drawn_mass > 0.0 and (b.m_f + b.m_d).sum() > 0.0
    assert_books(b)
    assert (b.skins.thickness() >= 0.0).all()
    assert (b.patch_of_skin >= 0).all() and int((b.skin_of_patch >= 0).sum()) == b.skins.n
    assert b.solver.film_mass[np.unique(b.surface.faces[b.skin_of_patch >= 0])].sum() == pytest.approx(
        b._film_nodal()[np.unique(b.surface.faces[b.skin_of_patch >= 0])].sum())


def test_deaths_hand_skins_to_the_faces_they_expose(coarse_sphere_mesh):
    b = make_body(coarse_sphere_mesh)
    start_at(b, 860.0)
    n0, t = b.mesh.n_active, 0.0
    while b.mesh.n_active == n0 and t < 30.0:
        t = run(b, 3e6, 1, t)
    assert b.mesh.n_active < n0 and b.skin_handover_mass > 0.0
    assert_books(b)
    assert (b.skins.m.sum(axis=1) > 0.0).all() and (b.patch_of_skin >= 0).all()


def test_a_large_interface_mismatch_repeats_the_step(coarse_sphere_mesh, monkeypatch):
    monkeypatch.setattr(skin, "INTERFACE_TOLERANCE", 0.0)
    b = make_body(coarse_sphere_mesh)
    start_at(b, 760.0)
    run(b, 1e6, 4)
    assert b.skins.n > 0 and b.interface_repeats >= 3
    assert_books(b)


def test_liquid_on_a_cold_patch_moves_onto_its_new_skin_with_its_energy(coarse_sphere_mesh):
    b = make_body(coarse_sphere_mesh)
    b.m_f[5] = 1e-6
    b.solver.set_film_mass(b._film_nodal())
    E, M = b.energy() + b.pending_load.sum(), b.mass(0.0)
    assert b._create_skins(b.solver.temperature()) == 1 and b.skin_of_patch[5] >= 0
    assert b.energy() + b.pending_load.sum() == pytest.approx(E, rel=1e-12) and b.mass(0.0) == pytest.approx(M, rel=1e-14)
    assert b.solver.film_mass.sum() == 0.0                      # the liquid now rides the skin's top, not the 3D nodes


def test_a_skin_its_element_cannot_supply_is_made_thin_and_counted(coarse_sphere_mesh):
    b = make_body(coarse_sphere_mesh)
    e = b.surface.owner[0]
    b.phi[e] = body.PHI_DEATH + 0.01
    b.solver.set_fractions(b.phi)
    b.m_f[0] = 1e-9
    b.solver.set_film_mass(b._film_nodal())
    b.mass0, b.energy0 = b.mass(0.0), b.energy()
    b._create_skins(b.solver.temperature())
    r = b.skin_of_patch[0]
    assert r >= 0 and 0.0 < b.skins.thickness()[r] < 0.5 * DEV.thickness
    assert b.phi[e] >= body.PHI_DEATH - 1e-12
    b.advance(0.5, 0.5, uniform(b, 1e5))
    assert b.skin_short_events >= 1                             # it drew nothing: its element had nothing to give
    assert b.phi[e] >= body.PHI_DEATH - 1e-12 or not b.mesh.active[e]   # ... and died at PHI_DEATH if it reached it
    assert_books(b)


def test_coupled_flight_steps_with_skins_keep_the_books(coarse_sphere_mesh):
    """Physics-mode loads at 71 km, as test_film_spraying_death_and_balances, with skins: film forms and is sprayed."""
    b = make_body(coarse_sphere_mesh)
    sim = simulator(b, t_max=60.0)
    sim.advance(43.0)
    start_at(b, 880.0)
    model = heating.PhysicsHeating()
    for _ in range(8):
        sim.advance(0.5)
        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
        loads = model.evaluate(a, b.theta, b.surface_temperature(), 0.05, T_mean=b.mean_temperature())
        b.advance(sim.t, 0.5, loads, state=a)
    assert b.skins.n > 0 and b.sprayed_mass > 0.0
    assert_books(b)
    stats = b.melt_stats()
    assert np.isfinite(stats["nonrigid_depth_wet_median_um"]) and np.isfinite(stats["liquid_depth_mean_um"])


def test_radiated_power_counts_skin_tops_not_the_facets_under_them(coarse_sphere_mesh):
    """The history's radiated power with skins: bare facets plus skin tops, never the facets under skins, after the
    melt step's surface refresh (Task 6 review ruling)."""
    b = make_body(coarse_sphere_mesh)
    start_at(b, 740.0)
    run(b, 1.5e6, 4)
    assert b.skins.n > 0
    on = b.skin_of_patch >= 0
    Tf = b.surface.facet_mean(b.solver.temperature())
    sig = b.emissivity * thermal.SIGMA_SB
    expected = float((sig * Tf[~on] ** 4 * b.surface.areas[~on]).sum() + (sig * b.skins.T[:, 0] ** 4 * b.skins.area).sum())
    assert b.radiated_power() == pytest.approx(expected, rel=1e-12)


def test_an_element_its_skins_drained_dies(coarse_sphere_mesh):
    """Skins draw an element down to PHI_DEATH only to rounding; a few ulps above it the element must still die and
    hand its remainder to its skins (Task 7 review ruling)."""
    b = make_body(coarse_sphere_mesh)
    start_at(b, 860.0)
    b._create_skins(b.solver.temperature())
    e = b.surface.owner[0]
    b.phi[e] = body.PHI_DEATH * (1.0 + 1e-13)
    b.solver.set_fractions(b.phi)
    b.mass0, b.energy0 = b.mass(0.0), b.energy()
    b.advance(0.5, 0.5, uniform(b, 1.5e6))
    assert not b.mesh.active[e]
    assert_books(b)


FREEZE_T0 = 880.0     # K: below T_low (906 K) and near enough to it for the top node's energy to bind. From 760 K or 850 K
                      # all ten skin-masses freeze: liquid that far below the liquidus still holds its latent heat, and
                      # warming it to T_low takes more than freezing releases there (measured: 100 % at both, 68 % here)


def test_freeze_back_into_a_skin_is_limited_by_energy(coarse_sphere_mesh):
    """Liquid piled on a skin whose top is below the feed ramp freezes into it on the sub-steps only as far as the top
    node's energy allows (spec 5.4, decided 2026-10-08): some freezes, the rest stays liquid, and the books close."""
    b = make_body(coarse_sphere_mesh, solver_options=EXACT)
    start_at(b, FREEZE_T0)
    b._create_skins(b.solver.temperature())
    face = b.surface.face_ids[0]
    s = DEV.thickness
    L0 = 10.0 * b.material.rho * b.surface.areas[0] * s
    b.m_f[0] = L0
    b.solver.set_film_mass(b._film_nodal())
    b.mass0, b.energy0 = b.mass(0.0), b.energy()
    b.advance(0.5, 0.5, uniform(b, 1e5))
    p = b.patch_of_face[face]
    assert p >= 0 and b.skin_frozen_mass > 0.0
    assert b.m_f[p] + b.m_d[p] > 0.0                                   # not all of it
    assert b.skin_frozen_mass == pytest.approx(L0 - (b.m_f[p] + b.m_d[p]), rel=1e-9)
    assert b.frozen_mass == pytest.approx(b.skin_frozen_mass, rel=1e-14)  # one counter of re-solidified liquid
    # the liquid left waits at T_low to freeze in the next sub-steps: not film the model cannot freeze (review F2-1)
    assert b.film_frozen_fraction() == 0.0
    assert_books(b)


def thick_skin(the_mesh, full, T0=760.0):
    """A body at T0 (760 K unless given) with skins, the skin of the patch whose element has the most metal per unit of
    its area made 2.5 s thick by moving metal from that element into its bottom (at the element's own enthalpy, so
    nothing is booked); with `full`, the element is then made full (phi = 1). The books are reckoned from there."""
    b = make_body(the_mesh, solver_options=EXACT)
    start_at(b, T0)
    b._create_skins(b.solver.temperature())
    owner = b.surface.owner
    p = int(np.argmax(b.phi[owner] * b.element_mass[owner] / b.surface.areas))
    r, e = b.skin_of_patch[p], owner[p]
    extra = 1.5 * b.material.rho * b.skins.area[r] * DEV.thickness
    h_e = float(b.material.enthalpy(b.solver.temperature())[b.mesh.tets[e]].mean())
    b.skins.add_bottom(np.array([r]), np.array([extra]), np.array([extra * h_e]), np.array([0.0]))
    b.phi[e] = 1.0 if full else b.phi[e] - extra / b.element_mass[e]
    assert b.phi[e] > 2.0 * body.PHI_DEATH
    b.solver.set_fractions(b.phi)
    b.mass0, b.energy0 = b.mass(0.0), b.energy()
    return b, b.surface.face_ids[p], e


def test_a_skin_above_twice_its_thickness_returns_metal_to_its_element(coarse_sphere_mesh):
    """Spec 5.4, "Thickness bounds": a skin above 2s returns its bottom metal in excess of 2s to its element where the
    element has room; where the element is full, the skin keeps it and the event is counted."""
    s2 = skin.MAX_THICKNESS_FACTOR * DEV.thickness
    b, face, e = thick_skin(coarse_sphere_mesh, full=False)
    r = b.skins.rows_of(np.array([face]))[0]
    assert b.skins.thickness()[r] > 1.2 * s2
    phi0 = b.phi[e]
    b.advance(0.5, 0.5, uniform(b, 1e5))
    r = b.skins.rows_of(np.array([face]))[0]
    assert r >= 0 and b.skins.thickness()[r] == pytest.approx(s2, rel=1e-9)          # back at the bound, not below it
    assert b.skin_returned_mass > 0.0 and b.skin_full_events == 0
    assert b.phi[e] - phi0 == pytest.approx(b.skin_returned_mass / b.element_mass[e], rel=1e-12)
    assert_books(b)
    b, face, e = thick_skin(coarse_sphere_mesh, full=True)
    b.advance(0.5, 0.5, uniform(b, 1e5))
    r = b.skins.rows_of(np.array([face]))[0]
    assert r >= 0 and b.skins.thickness()[r] > 1.2 * s2                # it kept its metal
    assert b.skin_returned_mass == 0.0 and b.skin_full_events >= 1 and b.phi[e] == 1.0
    assert_books(b)


def test_a_thick_skin_that_melts_is_recorded(coarse_sphere_mesh):
    """Skins on elements with no room grow beyond 2s by freeze-back (accepted 2026-10-08 on condition that it is
    measured whether such skins later melt): a skin above 2s at the start of a step that hands liquid to the film or
    the deep account is counted that step, with its start-of-step thickness, and in the cumulative events."""
    s2 = skin.MAX_THICKNESS_FACTOR * DEV.thickness
    b, face, e = thick_skin(coarse_sphere_mesh, full=True)
    assert b.skin_thick_feed_events == 0 and b.last_thick_feeding[0] == 0 and np.isnan(b.last_thick_feeding[1])
    r = b.skins.rows_of(np.array([face]))[0]
    th0 = float(b.skins.thickness()[r])
    assert th0 > 1.2 * s2
    b.advance(0.5, 0.5, uniform(b, 1e5))                               # cold: the thick skin feeds nothing
    assert b.last_thick_feeding[0] == 0 and np.isnan(b.last_thick_feeding[1])
    assert b.skin_thick_feed_events == 0 and b.melt_stats()["skin_thick_feeding_count"] == 0.0
    assert np.isnan(b.melt_stats()["skin_thick_feeding_max_mm"])
    b, face, e = thick_skin(coarse_sphere_mesh, full=True, T0=900.0)
    r = b.skins.rows_of(np.array([face]))[0]
    th0 = float(b.skins.thickness()[r])
    assert th0 > 1.2 * s2
    b.advance(0.5, 0.5, uniform(b, 3e6))
    assert b.skin_thick_feed_events >= 1
    stats = b.melt_stats()
    assert stats["skin_thick_feeding_count"] >= 1.0 and stats["skin_thick_feed_events"] == float(b.skin_thick_feed_events)
    assert stats["skin_thick_feeding_max_mm"] == pytest.approx(th0 * 1e3, rel=1e-9)   # the start-of-step thickness
    assert b.last_thick_feeding[0] == int(stats["skin_thick_feeding_count"])
    assert_books(b)


def test_round_off_above_twice_the_thickness_is_not_an_excess(coarse_sphere_mesh):
    """A skin returned to 2s sits at it only to rounding; an excess below 1e-9 of the bound -- the margin of
    skin_thick_count -- is none: nothing is returned where the element has room, and no event is counted where it is
    full (review F2-2)."""
    for full in (False, True):
        b = make_body(coarse_sphere_mesh)
        start_at(b, 760.0)
        b._create_skins(b.solver.temperature())
        r, e = b.skin_of_patch[0], b.surface.owner[0]
        bound = b.material.rho * b.skins.area[r] * skin.MAX_THICKNESS_FACTOR * DEV.thickness
        T = b.solver.temperature()
        h_e = b.material.enthalpy(T)[b.mesh.tets].mean(axis=1)
        extra = bound * (1.0 + 1e-15) - b.skins.m[r].sum()
        b.skins.add_bottom(np.array([r]), np.array([extra]), np.array([extra * h_e[e]]), np.array([0.0]))
        assert b.skins.m[r].sum() > bound                                  # above the bound by round-off only
        b.phi[e] = 1.0 if full else 0.5
        b._return_from_skins(h_e)
        assert b.skin_full_events == 0 and b.skin_returned_mass == 0.0
