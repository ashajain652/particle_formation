"""MeltingBody with skins (Step 3 sub-plan 18): identical to the element model until the first skin; skins at the
solidus drawing from their elements; exact books; deaths handing skins on; the interface repeat; the edge cases of the
plan's Review Focus."""
import numpy as np
import pytest

from reentry_model import body, heating, material, mesh, skin, thermal
from test_reentry_model_coupled import MASS_100MM, simulator

pytest.importorskip("cantera")

DEV = skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"], substep=0.05)


def make_body(the_mesh, surface_model="skin", name="AA7075_scheil", solver="skfem", **settings):
    mat = material.Material.from_drama_json(name)
    m = mesh.VolumeMesh(the_mesh.points, the_mesh.tets, dict(the_mesh.params), the_mesh.element_layer.copy())
    return body.MeltingBody(m, mat, thermal.thermal_solver(solver), MASS_100MM,
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
