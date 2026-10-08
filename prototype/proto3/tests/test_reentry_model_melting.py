"""body.MeltingBody: feed, instant removal in the lumped limit, film and spraying, element death with hand-over,
energy and mass balances, demise, and the shape feedback -- Knudsen length, nose-cap radius and drag (sections 8-11,
18)."""
import math
import os

import numpy as np
import pytest

from reentry_model import body, heating, material, mesh, thermal
from reentry_model.body import PHI_DEATH as PHI_DEATH_FOR_TEST
from test_reentry_model_coupled import MASS_100MM, simulator

pytest.importorskip("cantera")


@pytest.fixture(scope="module")
def layered_mesh(tmp_path_factory):
    return mesh.sphere_mesh(0.05, 4e-3, 20e-3, str(tmp_path_factory.mktemp("meshes")), layers=2)


def melting_body(the_mesh, name="AA7075_range", k_scale=1.0, **settings):
    mat = material.Material.from_drama_json(name)
    if k_scale != 1.0:
        mat.k_table = mat.k_table * k_scale
    m = mesh.VolumeMesh(the_mesh.points, the_mesh.tets, dict(the_mesh.params), the_mesh.element_layer.copy())
    return body.MeltingBody(m, mat, thermal.thermal_solver("skfem"), MASS_100MM, settings=body.MeltSettings(**settings))


def test_settings_and_construction(layered_mesh):
    with pytest.raises(ValueError):
        body.MeltSettings(removal="magic")
    with pytest.raises(ValueError):
        body.MeltSettings(demise_fraction=1.5)
    with pytest.raises(ValueError):
        body.MeltingBody(layered_mesh, material.Material.from_drama_json("AA7075_nomelt"), thermal.thermal_solver("skfem"), MASS_100MM)
    b = melting_body(layered_mesh)
    assert b.mass(0.0) == pytest.approx(b.mass0) and b.mass0 == pytest.approx(MASS_100MM, rel=3e-3)
    assert b.reference_area() == pytest.approx(math.pi * 0.05 ** 2, rel=2e-3) and b.equivalent_radius() == pytest.approx(0.05, rel=1e-3)
    assert b.m_f.shape == (b.surface.n_patches,) and b.m_f.sum() == 0.0 and not b.demised() and b.melt_front_depth() == 0.0
    stats = b.melt_stats()
    assert stats["mass_kg"] == b.mass(0.0) and stats["n_active_elements"] == layered_mesh.n_elements and stats["sprayed_mass_kg"] == 0.0


def test_instant_removal_in_the_lumped_limit_follows_q_over_l(layered_mesh):
    """k x 1e4, uniform 3e5 W/m2 on the whole surface, AA7075: once the body sits on the 850 K plateau the mass leaves
    at (Q_conv - Q_rad) / L_f with the geometry intact (SESAM's lumped law), the energy balance exact."""
    b = melting_body(layered_mesh, "AA7075", k_scale=1e4, removal="instant", runoff=False, size_feedback="initial")
    q = np.full(b.surface.n_patches, 3e5)
    loads = heating.HeatingResult(q, 3e5, 0.0, 0.0, 0.0)
    t, dt = 0.0, 0.5
    while b.melt_onset is None:
        t += dt
        b.advance(t, dt, loads)
    m1 = b.mass(0.0)
    for _ in range(20):
        t += dt
        b.advance(t, dt, loads)
    Q_net = b.last.Q_conv - b.last.Q_rad
    assert (m1 - b.mass(0.0)) / (20 * dt) == pytest.approx(Q_net / 4e5, rel=2e-2)
    assert b.removed_mass == pytest.approx(b.mass0 - b.mass(0.0)) and b.m_f.sum() == 0.0 and b.sprayed_mass == 0.0
    assert abs(b.energy_balance_residual()) < 1e-8 and b.mesh.n_active == layered_mesh.n_elements    # nothing dead yet: phi shrinks uniformly
    assert 0.5 < b.phi.min() < b.phi.max() < 1.0 and 848.0 < b.mean_temperature() < 852.0
    assert b.reference_area() == pytest.approx(math.pi * 0.05 ** 2, rel=2e-3)
    assert b.melt_onset is not None and b.spray_onset is None


def test_film_spraying_death_and_balances(layered_mesh):
    """Physics-mode loads at 71 km on the real body for 8 s: film forms, is stripped, the outer layer dies on the
    windward face with the film handed over, mass and energy balances hold, the source table grows."""
    b = melting_body(layered_mesh)
    sim = simulator(b, t_max=60.0)
    sim.advance(43.0)                                                                 # to 71 km on the constant-mass trajectory
    b.solver.set_temperature(880.0)                                                   # a hot body: the surface melts at once
    b.energy0 = b.energy()                                                            # the balance is reckoned from here
    model = heating.PhysicsHeating()
    for _ in range(16):
        sim.advance(0.5)
        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
        loads = model.evaluate(a, b.theta, b.surface_temperature(), 0.05, T_mean=b.mean_temperature())
        b.advance(sim.t, 0.5, loads, state=a)
    stats = b.melt_stats()
    assert b.sprayed_mass > 0.0 and b.n_released > 0.0 and len(b.source_rows) > 0 and stats["n_released"] == b.n_released
    assert b.removed_mass == pytest.approx(b.sprayed_mass) and b.mass(0.0) == pytest.approx(b.mass0 - b.sprayed_mass, rel=1e-6)
    assert b.mesh.n_active < layered_mesh.n_elements and b.mesh.element_layer[~b.mesh.active].max() <= 1     # only layer elements have died
    assert abs(b.energy_balance_residual()) < 1e-7
    assert b.surface.n_patches == b.m_f.size == len(b.solver.areas) and b.theta.size == b.surface.n_patches
    assert stats["theta_cr_deg"] < 45.0 and stats["spraying_area_m2"] > 0.0 and 20.0 < stats["r_median_um"] < 2000.0
    assert stats["closure_fraction_girin"] + stats["closure_fraction_couette_slip"] + stats["closure_fraction_couette_fm"] == pytest.approx(1.0)
    assert stats["kn_body"] > 0.0 and stats["kn_local_stag"] < stats["kn_body"] and stats["re_shock"] > 0.0      # the wall gas is compressed
    assert stats["flow_branch"] in (0.0, 1.0, 2.0) and stats["p_w_stag_Pa"] > stats["kn_body"] * 0.0
    assert stats["film_thickness_max_mm"] >= stats["film_thickness_mean_mm"] >= 0.0
    h_out = stats["removed_enthalpy_J"] / b.sprayed_mass                                           # the droplets leave at the film's
    assert b.material.h_liquid < h_out <= b.material.enthalpy(b.surface_temperature().max())       # own temperature, superheat and all
    assert stats["film_T_max_K"] == pytest.approx(b.surface_temperature().max()) and stats["film_T_max_K"] > b.material.T_liquidus
    rows = np.array(b.source_rows)
    assert rows.shape[1] == 22 and np.all(rows[:, 14] > 0.0) and np.all(np.isfinite(rows[:, 12]))
    assert b.hist_n.sum() == pytest.approx(b.n_released) and b.hist_m.sum() == pytest.approx(b.sprayed_mass)
    assert not b.demised() and b.reference_area() < math.pi * 0.05 ** 2                            # the windward face has receded


def test_a_dying_element_leaves_its_last_solid_on_the_face_it_exposes(layered_mesh):
    """An element is removed while it still holds up to PHI_DEATH of its material (typically far less). That remainder
    has to end up where the surface receded to -- on the faces the element itself has just exposed -- exactly like the
    film it was already carrying, and it does, because the remainder is credited to the element's own patches first and
    is then handed over with them (spec amendment 20). Held below the feed ramp so that nothing else melts this step."""
    b = melting_body(layered_mesh, name="AA7075")
    owner = int(b.surface.owner[0])
    b.m_f[:] = 0.0
    b.solver.set_temperature(b.material.T_feed - 30.0)             # below the feed ramp: no other element melts
    b.phi[owner] = 0.5 * PHI_DEATH_FOR_TEST                        # below the death threshold, with solid left
    b.solver.set_fractions(b.phi)
    rest = float(b.phi[owner] * b.element_mass[owner])
    faces = b.mesh.element_faces([owner])[0]
    mass_before = b.mass(0.0)
    b.melt_step(0.0, 0.5, None)
    exposed = np.array([q for q in b.patch_of_face[faces] if q >= 0])
    assert not b.mesh.active[owner] and exposed.size > 0            # it died and uncovered at least one face
    assert b.m_f.sum() == pytest.approx(rest) and b.m_f[exposed].sum() == pytest.approx(rest)
    assert b.mass(0.0) + b.removed_mass == pytest.approx(mass_before, rel=1e-12)


def test_demise_and_consumption(layered_mesh):
    b = melting_body(layered_mesh, "AA7075", k_scale=1e4, removal="instant", runoff=False, demise_fraction=0.5, size_feedback="initial")
    loads = heating.HeatingResult(np.full(b.surface.n_patches, 2e6), 2e6, 0.0, 0.0, 0.0)
    t = 0.0
    while not b.demised():
        t += 0.5
        b.advance(t, 0.5, loads)
        assert t < 200.0
    assert b.mass(0.0) < 0.5 * b.mass0 and not b.consumed and b.mesh.n_active > 0


def test_size_feedback_nose_fit_and_knudsen_length(layered_mesh):
    """Intact sphere: the windward-cap fit returns R0 and the Knudsen length is D0; a flattened front (the cap patches
    pushed onto the plane x = 0.8 R0, then 0.6 R0) fits a larger radius, capped at NOSE_CAP_FACTOR x the transverse
    radius; with size_feedback 'initial' both stay at their initial values. The trajectory's Kn follows the body's length."""
    b = melting_body(layered_mesh)
    assert b.nose_radius() == pytest.approx(0.05, rel=5e-3) and b.reference_length() == pytest.approx(0.1, rel=1e-3)
    assert b.transverse_radius == pytest.approx(0.05, rel=2e-3) and np.linalg.norm(b.mass_centre) < 1e-3
    centre, r = body.fit_sphere(b.surface.centroids)
    assert np.linalg.norm(centre) < 1e-4 and r == pytest.approx(0.05, rel=5e-3)
    original = b.surface.centroids.copy()
    flat = original.copy()
    flat[original[:, 0] > 0.8 * 0.05, 0] = 0.8 * 0.05
    b.surface.centroids = flat
    b._fit_nose()
    assert b.fitted_nose_radius > 0.06                                             # flat centre + curved shoulders: ~84 mm
    assert b.nose_radius() == pytest.approx(min(b.fitted_nose_radius, body.NOSE_CAP_FACTOR * b.transverse_radius))
    flat = original.copy()
    flat[original[:, 0] > 0.6 * 0.05, 0] = 0.6 * 0.05                               # a wider flat face: the fit exceeds the cap
    b.surface.centroids = flat
    b._fit_nose()
    assert b.fitted_nose_radius > body.NOSE_CAP_FACTOR * b.transverse_radius and b.nose_radius() == pytest.approx(body.NOSE_CAP_FACTOR * b.transverse_radius)
    fixed = melting_body(layered_mesh, size_feedback="initial")
    fixed.surface.centroids = flat
    fixed._fit_nose()
    assert fixed.nose_radius() == 0.05 and fixed.reference_length() is None
    with pytest.raises(ValueError):
        body.MeltSettings(size_feedback="magic")
    # the trajectory's Knudsen number uses the body's length: halve the body's mass and compare
    sim = simulator(b)
    a0 = sim.aero_state(0.0, sim.y[:3], sim.y[3:])
    b.phi *= 0.125                                                                  # equivalent diameter halves
    a1 = sim.aero_state(0.0, sim.y[:3], sim.y[3:])
    assert a1.kn == pytest.approx(2.0 * a0.kn, rel=1e-6) and a0.kn == pytest.approx(0.0298, rel=0.02)
    assert melt_stats_has_radii(b)


def melt_stats_has_radii(b):
    s = b.melt_stats()
    return s["nose_radius_mm"] > 0.0 and s["transverse_radius_mm"] > 0.0 and s["fitted_nose_radius_mm"] > 0.0


def test_drag_shape_factor_between_the_two_atdb_endpoints(layered_mesh):
    """The shape factor is exactly 1 for the intact sphere (so the ATDB sphere table is reproduced bit for bit) and
    reaches the ATDB disc entry at the flat-face limit: HTG's continuum database is modified Newtonian, its
    sphere/disc ratio being 0.49897 at every Mach number, and the same integral over a meshed flat plate gives 2.00."""
    import json
    from reentry_model import aero, mesh as mesh_mod
    from helpers import REPO_ROOT
    b = melting_body(layered_mesh)
    assert b.drag_shape_factor() == pytest.approx(1.0, abs=1e-12)
    cd, area = b.newtonian_drag()
    assert cd == pytest.approx(0.92, rel=0.02) and area == pytest.approx(math.pi * 0.05 ** 2, rel=0.02)
    sphere = json.load(open(os.path.join(REPO_ROOT, "reentry_model", "data", "atdb_sphere.json")))
    disc = json.load(open(os.path.join(REPO_ROOT, "reentry_model", "data", "atdb_disc.json")))
    ratios = np.array(sphere["cd_continuum"]) / np.array(disc["cd_continuum"])
    assert np.allclose(ratios, 0.49897, atol=1e-5) and disc["mach"] == sphere["mach"]
    assert np.allclose(np.array(disc["cd_free_molecular"]) / np.array(sphere["cd_free_molecular"]), 1.03, atol=0.04)
    # a flat plate normal to the flow: every windward element is square-on, so the integral is C_p,max itself
    plate = mesh_mod.box_mesh(0.004, 0.06, 0.06, 0.002)
    flat = body.MeltingBody(plate, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"),
                            1.0, settings=body.MeltSettings())
    flat_cd = flat.newtonian_drag()[0]
    assert flat_cd == pytest.approx(body.CP_MAX_NEWTONIAN, rel=1e-6)
    assert flat_cd / b.newtonian_drag_sphere == pytest.approx(1.0 / 0.49897, rel=0.02)     # 2.00 vs ATDB's 2.0041
    fixed = melting_body(layered_mesh, size_feedback="initial")
    assert fixed.drag_shape_factor() == 1.0 and fixed.reference_area() == pytest.approx(math.pi * 0.05 ** 2, rel=2e-3)
    # the trajectory multiplies the continuum entry only
    tables, bridging = aero.SphereDragTables.from_json(), aero.SesamTable()
    assert aero.drag_coefficient(0.0, 10.0, tables, bridging, 2.0) == pytest.approx(2.0 * tables.cd_continuum(10.0))
    assert aero.drag_coefficient(1e6, 10.0, tables, bridging, 2.0) == pytest.approx(tables.cd_free_molecular(10.0))


def test_film_temperature_freeze_back_and_the_netted_transfer(layered_mesh):
    """2.2 MW/m2 for 6 s on a 100 mm body, then the heating off for 60 s, with no flow (so nothing runs off or
    sprays): the film takes the surface's temperature and its enthalpy is the mean of the nodes' -- not the enthalpy
    of the mean temperature, which differs by the whole latent heat on the ramp -- it freezes back onto its owner
    element once the surface falls through the feed ramp, and feed and freeze are netted, so neither runs while the
    other does. Mass and the coupled energy balance are exact throughout."""
    import types
    b = melting_body(layered_mesh, name="AA7075")
    b.solver.set_temperature(800.0)                                    # just below AA7075's melting point
    b.energy0 = b.energy()
    loads = lambda q: types.SimpleNamespace(q_conv=np.full(b.surface.n_patches, q))
    peak, tail = 0.0, 0.0
    for i in range(132):
        b.advance(i * 0.5, 0.5, loads(2.2e6 if i < 16 else 0.0))
        peak = max(peak, b.m_f.sum())
        if i >= 112:                                                   # the last 20 steps: nothing is being heated
            tail += b.last_melt["feed_mass"]
        assert abs(b.energy_balance_residual()) < 1e-6
        assert b.mass(0.0) + b.removed_mass == pytest.approx(b.mass0, rel=1e-12)
    assert tail / 20.0 < 1e-3 * peak   # netted: the feed decays away instead of cycling the film every step (run as
    #                                    two independent rates it sat at 3 % of the body per step with the film static)
    T, mat = b.solver.temperature(), b.material
    assert b.film_enthalpy() == pytest.approx(mat.enthalpy_liquid(T)[b.surface.faces].mean(axis=1))
    mix = mat.enthalpy(T)[b.surface.faces].mean(axis=1)                 # the film carries its latent heat wherever it
    assert (b.film_enthalpy() >= mix - 1e-6).all()                      # sits, so it is never below the mixture's
    assert (b.film_enthalpy() > mix + 0.2 * mat.latent_heat).any()      # and is well above it where the surface cooled
    assert b.film_energy() == pytest.approx((b.m_f * b.film_enthalpy()).sum())
    assert b.film_temperature() == pytest.approx(b.surface_temperature())
    assert peak > 0.0 and b.frozen_mass > 0.0 and b.phi.max() <= 1.0     # the film came back and phi never overflows
    assert b.solver.film_mass.sum() == pytest.approx(b.m_f.sum())        # the solver carries the same film
    assert b.melt_stats()["film_T_max_K"] == pytest.approx(b.surface_temperature().max())


def test_a_deferred_melt_load_never_moves_a_node_more_than_the_cap(layered_mesh):
    """A drained node has nothing to heat with the enthalpy its melt delivered: no node is moved more than
    LOAD_DT_MAX in one step, the remainder waits in `pending_load`, and the balance counts it either way."""
    import types
    b = melting_body(layered_mesh, name="AA7075")
    loads = lambda q: types.SimpleNamespace(q_conv=np.full(b.surface.n_patches, q))
    b.advance(0.0, 0.5, loads(0.0))
    b.pending_load[:] = 0.0
    b.pending_load[np.unique(b.surface.faces)] = 1.0e4                   # 10 kJ on every boundary node
    queued, T0 = b.pending_load.sum(), b.solver.temperature().copy()
    room = body.LOAD_DT_MAX * b.solver.nodal_capacity()
    b.advance(0.5, 0.5, loads(0.0))
    assert b.pending_load.sum() == pytest.approx(queued - np.minimum(1.0e4, room[np.unique(b.surface.faces)]).sum())
    assert (b.solver.temperature() - T0).max() < 1.05 * body.LOAD_DT_MAX
    assert b.pending_load.sum() + b.total_applied_load == pytest.approx(queued)     # nothing is lost on the way


def test_a_dying_element_hands_its_film_to_the_faces_it_exposes(layered_mesh):
    """When an element is removed the surface recedes into the film it carried, so that film belongs on the faces of
    that same element which have just become boundary -- its inward face and the walls of the pit it opens -- not on
    whichever surviving centroid happens to be nearest (which is a statement about the mesh's numbering, and which
    once put a third of a collapsing body's film onto one 0.6 mm2 facet)."""
    b = melting_body(layered_mesh, name="AA7075")
    owner = b.surface.owner[0]                                      # kill one patch owner, with film only on its patch
    b.m_f[:] = 0.0
    b.m_f[b.surface.owner == owner] = 1.0e-4
    carried, faces = float(b.m_f.sum()), b.mesh.element_faces([owner])[0]
    gone_normal = b.surface.normals[b.surface.owner == owner][0]
    b._kill(np.array([owner]))
    exposed = np.array([q for q in b.patch_of_face[faces] if q >= 0])
    assert exposed.size and b.m_f.sum() == pytest.approx(carried)    # nothing lost, and it all landed on that element's
    assert b.m_f[exposed].sum() == pytest.approx(carried)            # own newly exposed faces
    # shared by the area each exposed face presents to the vanished patch, i.e. its area projected on that patch's
    # outward normal: the floor the recession uncovered takes the film, the pit walls standing perpendicular take none
    cosine = np.clip(b.surface.normals[exposed] @ gone_normal, 0.0, None)
    w = b.surface.areas[exposed] * cosine
    assert b.m_f[exposed] == pytest.approx(carried * w / w.sum())
    floor = exposed[int(np.argmax(cosine))]                         # the face most nearly parallel to the old patch
    assert b.m_f[floor] > 0.5 * carried and cosine.max() > 0.9      # takes the bulk of it
