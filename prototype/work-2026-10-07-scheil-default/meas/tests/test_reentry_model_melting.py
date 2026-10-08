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
    # the histogram bins 1 um to 10 mm and drops the rest by design. With the rigid substrate on (amendment of 2026-10-06)
    # this body, 880 K throughout, is slurry to the core: off Girin's closure its film is held from the shear modes,
    # piles up, and the front-surface Rayleigh-Taylor mode sheds some of it above 10 mm (measured in two unseeded runs:
    # 5 and 6 droplets, 5.2 and 6.1 % of the mass, 10.1-10.9 mm), so the bins hold exactly what was released inside
    # their range and nothing is lost below it
    from reentry_model import spray
    r_rows, dn_rows, dm_rows = rows[:, 12], rows[:, 13], rows[:, 14]
    inside = (r_rows >= spray.R_MIN) & (r_rows <= spray.R_MAX)
    assert b.hist_n.sum() == pytest.approx(dn_rows[inside].sum()) and b.hist_m.sum() == pytest.approx(dm_rows[inside].sum())
    assert dn_rows.sum() == pytest.approx(b.n_released) and dm_rows.sum() == pytest.approx(b.sprayed_mass)
    assert (r_rows[~inside] > spray.R_MAX).all() and (r_rows <= 0.25 * 0.05).all()              # above the bins, below R/4
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
    # the molten cascade (amendment of 2026-10-03) stays off here: this device puts 2.2 MW/m2 on every facet, crater walls
    # and lee included, of a body 50 K below its melting point, and once the surface may recede through molten material
    # within the step the whole body melts in the 8 s of heating (measured: consumed at step 16, the books exact
    # throughout), leaving no film to freeze back; the cascade's own books are tested with the amendment's tests below
    b = melting_body(layered_mesh, name="AA7075", molten_cascade=False)
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


# ---------------------------------------------------------------------------------------------------------------
# Amendment of 2026-10-02: the liquid below the conjugate depth runs off, and only the skin sprays


def girin_state(b, t=50.0):
    """The 100 mm flight at 69.8 km (t = 50 s on the intact sphere's trajectory): the shock-layer branch, every windward
    patch under Girin's closure with a conjugate depth of 0.17-0.30 mm, and G positive from the nose to about 80 deg
    (measured 2026-10-02)."""
    sim = simulator(b, t_max=90.0)
    sim.advance(t)
    return sim, sim.aero_state(sim.t, sim.y[:3], sim.y[3:])


def molten_pool(b, T_pool=960.0, T_rest=850.0, depth=3e-3, theta_max=60.0):
    """A molten pool: the nodes within `depth` of the surface and `theta_max` of the flight direction (+x) above the
    liquidus, the rest of the body below it."""
    p = b.mesh.points
    r = np.linalg.norm(p, axis=1)
    near = (r > 0.05 - depth) & (p[:, 0] > r * math.cos(math.radians(theta_max)))
    b.solver.set_temperature(np.where(near, T_pool, T_rest))


def flow_and_delta_m(b, a):
    from reentry_model import spray
    flow = b.flow.evaluate(a, b.theta, b.nose_radius(), b.liquid.rho, b.surface_temperature())
    return flow, spray.melt_layer(flow, b.liquid)[0]


def deep_stage_args(b, delta_m):
    """What melt_step hands the deep stage: the shear-reach gate, the temperatures and the three enthalpies."""
    T, mat = b.solver.temperature(), b.material
    h_node = mat.enthalpy(T)
    return b._shear_reaches(delta_m), T, h_node, h_node[b.mesh.tets].mean(axis=1), b.film_enthalpy(mat.enthalpy_liquid(T))


def test_the_deep_runoff_moves_the_liquid_below_the_conjugate_depth_and_books_it_exactly(layered_mesh):
    """Fact 28(b) holds melt in its element unless the gas shear reaches it; fact 29 asked that the pressure-gradient
    and deceleration-driven flux integrate over the whole liquid depth. A pool 3 mm deep under the nose at 70 km: the
    deep stage takes liquid only from the buried elements the feed gate held back (never from an owner, never from an
    element the shear reaches, never below PHI_MIN), moves it outward as G pushes it there, and books mass and
    enthalpy exactly -- what the elements lost is on the patches, and every enthalpy difference is in the deferred
    loads. Only the skin is film afterwards: deep liquid became film only as far as the conjugate depth had room."""
    from reentry_model import surface_flow as sf
    b = melting_body(layered_mesh)
    sim, a = girin_state(b)
    molten_pool(b)
    flow, delta_m = flow_and_delta_m(b, a)
    assert flow.branch == sf.BRANCH_SHOCK_LAYER and np.isfinite(delta_m[b.windward]).all()
    reach, T, h_node, h_e, h_p = deep_stage_args(b, delta_m)
    phi0, mass0 = b.phi.copy(), b.mass(0.0)
    E0 = b.energy() + b.pending_load.sum()
    seen = b._deep_runoff(0.5, flow, delta_m, reach, T, h_node, h_e, h_p)
    b.solver.set_fractions(b.phi)                                   # as melt_step does at its end
    assert seen > 0.0 and b.deep_runoff_mass > 0.0                  # liquid lay below the conjugate depth, and it moved
    assert b.deep_runoff_mass == pytest.approx(((phi0 - b.phi) * b.element_mass).sum(), rel=1e-12)   # what left the elements
    assert b.deep_surfaced_mass == pytest.approx(b.m_f.sum(), rel=1e-12)                              # what became film
    assert b.mass(0.0) == pytest.approx(mass0, rel=1e-12)
    assert b.energy() + b.pending_load.sum() == pytest.approx(E0, rel=1e-12)
    lost = phi0 - b.phi
    owners = b.owner_area > 0.0
    assert (lost >= 0.0).all() and lost.sum() > 0.0 and not lost[owners].any() and not lost[reach > 0.0].any()
    assert (b.phi[lost > 0.0] >= body.PHI_MIN * (1.0 - 1e-12)).all()
    left = lost * b.element_mass                                    # where it came from and where it went: outward
    c = b.mesh.points[b.mesh.tets].mean(axis=1)
    theta_from = np.arccos(c[:, 0] / np.linalg.norm(c, axis=1))
    arrived = b.m_d + b.m_f
    assert (arrived * b.theta).sum() / arrived.sum() > (left * theta_from).sum() / left.sum()
    skin = b.liquid.rho * b.surface.areas * np.where(b.windward, delta_m, np.inf)
    deep = b.m_d > 0.0
    assert deep.any() and (b.m_f[deep] >= skin[deep] * (1.0 - 1e-12)).all() and (b.m_f <= skin * (1.0 + 1e-12)).all()
    assert not b.m_d[~b.windward].any()                             # no conjugate depth on the lee: no deep liquid there


def test_the_deep_liquid_is_never_sprayed_only_the_skin_is(layered_mesh):
    """A millimetre of deep liquid on every windward patch and no film, over a body with nothing molten beneath its
    surface: the deep stage tops the film up to the conjugate depth from the top of the deep liquid and leaves the
    rest deep; spraying then strips only the film and leaves the deep liquid exactly as it was."""
    b = melting_body(layered_mesh)
    sim, a = girin_state(b)
    b.solver.set_temperature(905.0)                                 # below T_feed: no element holds liquid
    flow, delta_m = flow_and_delta_m(b, a)
    liq = b.liquid
    b.m_d = np.where(b.windward, 1e-3 * liq.rho * b.surface.areas, 0.0)
    reach, T, h_node, h_e, h_p = deep_stage_args(b, delta_m)
    total = b.m_d.sum()
    seen = b._deep_runoff(0.5, flow, delta_m, reach, T, h_node, h_e, h_p)
    assert seen == pytest.approx(total) and b.m_f.sum() + b.m_d.sum() == pytest.approx(total, rel=1e-12)
    skin = liq.rho * b.surface.areas * np.where(b.windward, delta_m, np.inf)
    deep = b.m_d > 0.0
    assert deep.any() and b.m_f[deep] == pytest.approx(skin[deep], rel=1e-12)       # where any stays deep, the skin is full
    film0, deep0 = b.m_f.copy(), b.m_d.copy()
    assert b.deep_surfaced_mass == pytest.approx(film0.sum(), rel=1e-12) and b.deep_runoff_mass == 0.0   # none from elements
    assert 0.0 <= b.deep_blob_fraction() <= 1.0
    released = b._film_and_spray(sim.t, 0.5, a, h_p, flow, delta_m)
    assert 0.0 < released <= film0.sum() * (1.0 + 1e-12)            # the skin sprays ...
    assert np.array_equal(b.m_d, deep0)                             # ... and not a gram of the deep liquid
    assert b.m_f.sum() + released == pytest.approx(film0.sum(), rel=1e-12)


def test_nothing_changes_where_no_liquid_lies_below_the_conjugate_depth(layered_mesh):
    """The deep runoff adds a transfer; it changes nothing else. At 70 km, with Girin's closure on every windward patch
    and a surface just above the liquidus over a body that is not fully molten beneath its owner elements, there is
    no liquid below the conjugate depth: a whole melt step is bit-identical with the deep runoff on and off. (Where the
    closure is not Girin's there is no conjugate depth at all; the 50 mm flight, which never has it, is the
    flight-level check, measured 2026-10-02.)"""
    runs = []
    for deep in (True, False):
        b = melting_body(layered_mesh, deep_runoff=deep)
        sim, a = girin_state(b)
        boundary = np.zeros(b.mesh.n_nodes, dtype=bool)
        boundary[np.unique(b.surface.faces)] = True
        b.solver.set_temperature(np.where(boundary, 930.0, 880.0))  # owners molten, nothing beneath them
        b.melt_step(sim.t, 0.5, a)
        runs.append(b)
    on, off = runs
    assert on.last_melt["deep_liquid"] == 0.0 and on.deep_runoff_mass == 0.0 and on.sprayed_mass > 0.0
    for name in ("phi", "m_f", "m_d", "pending_load"):
        assert np.array_equal(getattr(on, name), getattr(off, name)), name
    assert on.sprayed_mass == off.sprayed_mass and on.removed_enthalpy == off.removed_enthalpy
    assert on.runoff_mass == off.runoff_mass and np.array_equal(on.last_delta_m, off.last_delta_m, equal_nan=True)


def test_full_steps_with_the_deep_runoff_keep_the_books_exact(layered_mesh):
    """Six coupled macro steps from 70 km with a molten pool under the nose and the physics loads: the deep runoff
    acts (liquid below the conjugate depth is found and moved), elements die with deep liquid on their patches, and
    the energy balance and the mass stay exact by the same measures as every other melting test."""
    b = melting_body(layered_mesh)
    sim, a = girin_state(b)
    molten_pool(b)
    b.energy0 = b.energy()
    model = heating.PhysicsHeating()
    seen = 0.0
    for _ in range(6):
        sim.advance(0.5)
        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
        loads = model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature())
        b.advance(sim.t, 0.5, loads, state=a)
        seen += b.last_melt["deep_liquid"]
        assert abs(b.energy_balance_residual()) < 1e-7
        assert b.mass(0.0) + b.removed_mass == pytest.approx(b.mass0, rel=1e-12)
        assert b.m_d.shape == (b.surface.n_patches,) and (b.m_d >= 0.0).all()
        assert b.solver.film_mass.sum() == pytest.approx(b.m_f.sum() + b.m_d.sum())   # the solver carries both
    assert seen > 0.0 and b.deep_runoff_mass > 0.0 and b.sprayed_mass > 0.0 and b.mesh.n_active < layered_mesh.n_elements
    stats = b.melt_stats()
    assert stats["deep_runoff_mass_kg"] == b.deep_runoff_mass and stats["deep_mass_kg"] == pytest.approx(b.m_d.sum())
    assert stats["deep_surfaced_mass_kg"] == b.deep_surfaced_mass and 0.0 <= stats["deep_blob_fraction"] <= 1.0
    assert stats["deep_liquid_kg"] == b.last_melt["deep_liquid"]


# ---------------------------------------------------------------------------------------------------------------
# Amendment of 2026-10-03: the molten cascade, the surface receding through molten material within the step


def fully_molten(b):
    """Per element: fully molten as `molten_depth` has it -- the mean nodal temperature at or above the top of the feed
    ramp -- which is the test the molten cascade uses."""
    return b.solver.temperature()[b.mesh.tets].mean(axis=1) >= b.material.T_feed


def test_the_molten_cascade_empties_a_fully_molten_column_within_one_step(layered_mesh):
    """Fact 50: the melt step runs after the conduction, so the material beneath the wall-owning element melts as well,
    and the element the owner's death exposes was buried a moment earlier -- the feed gate held its melt in place --
    so it waited, fully molten, for the next step's feed: the surface receded through molten material by one element
    per macro step. The molten cascade feeds such an element within the step, so that it dies in turn and exposes the
    next. A pool 4 mm deep under the nose, every node in it above the top of the feed ramp, over a body at 850 K: with
    the cascade one melt step leaves no newly exposed wall-owning element fully molten -- the surface has receded
    through the whole pool to material that is not -- every element the cascade fed was fully molten, and its whole
    mass went to the film; without it the elements the first deaths exposed are still there, fully molten. Mass and
    enthalpy books exact."""
    out = {}
    for cascade in (True, False):
        b = melting_body(layered_mesh, molten_cascade=cascade)
        molten_pool(b, depth=4e-3)
        full, start = fully_molten(b), b.owner_area > 0.0
        mass0, E0 = b.mass(0.0), b.energy() + b.pending_load.sum()
        b.melt_step(0.0, 0.5, None)
        assert b.mass(0.0) == pytest.approx(mass0, rel=1e-12)
        assert b.energy() + b.pending_load.sum() == pytest.approx(E0, rel=1e-12)
        out[cascade] = (b, full, start)
    (on, full, start), (off, _, _) = out[True], out[False]
    exposed = lambda b: (b.owner_area > 0.0) & b.mesh.active & ~start          # the wall this step's deaths exposed
    assert (exposed(off) & full).sum() > 0                              # without: the exposed layer waits, fully molten
    assert exposed(on).any() and not (exposed(on) & full).any()        # with: the surface receded through the pool
    lm = on.last_melt
    assert lm["cascade_passes"] >= 2 and not lm["cascade_capped"] and on.cascade_mass == lm["cascade_mass"] > 0.0
    fed = off.phi - on.phi > 0.0                                        # what the cascade took beyond the step without it
    assert fed.any() and full[fed].all() and not on.mesh.active[fed].any()     # fully molten elements only, now gone
    assert on.cascade_mass == pytest.approx(((off.phi - on.phi) * on.element_mass)[fed].sum(), rel=1e-12)
    assert on.m_f.sum() - off.m_f.sum() == pytest.approx(on.cascade_mass, rel=1e-9)   # ... all of it into the film
    assert off.cascade_mass == 0.0 and off.last_melt["cascade_passes"] == 0 and not off.last_melt["cascade_capped"]
    stats = on.melt_stats()
    assert stats["cascade_passes"] == lm["cascade_passes"] and stats["cascade_mass_kg"] == on.cascade_mass


def test_the_molten_cascade_changes_nothing_where_the_exposed_elements_are_not_fully_molten(layered_mesh):
    """The cascade adds a feed and changes nothing else. At 69.8 km under Girin's closure, with the surface nodes at
    960 K over a body at 850 K and a millimetre of deep liquid on every windward patch: the wall-owning elements (three
    surface nodes of four) are molten, so the feed takes three quarters of what they hold and they die; every element
    their deaths expose has at most two surface nodes, a mean temperature of at most 905 K, and so is not fully molten;
    the deep runoff moves the deep liquid and the film sprays -- and the whole melt step is bit-identical with the
    cascade on and off."""
    runs = []
    for cascade in (True, False):
        b = melting_body(layered_mesh, molten_cascade=cascade)
        sim, a = girin_state(b)
        boundary = np.zeros(b.mesh.n_nodes, dtype=bool)
        boundary[np.unique(b.surface.faces)] = True
        b.solver.set_temperature(np.where(boundary, 960.0, 850.0))
        owners = b.owner_area > 0.0
        b.phi[owners] = 0.15                                             # the feed takes 3/4 of it: the owners die
        b.solver.set_fractions(b.phi)
        b.m_d = np.where(b.windward, 1e-3 * b.liquid.rho * b.surface.areas, 0.0)
        b.melt_step(sim.t, 0.5, a)
        runs.append(b)
    on, off = runs
    assert on.last_melt["n_dead"] > 0 and on.sprayed_mass > 0.0 and on.deep_surfaced_mass > 0.0
    assert on.cascade_mass == 0.0 and on.last_melt["cascade_passes"] == 0 and not on.last_melt["cascade_capped"]
    for name in ("phi", "m_f", "m_d", "pending_load"):
        assert np.array_equal(getattr(on, name), getattr(off, name)), name
    assert np.array_equal(on.mesh.active, off.mesh.active) and on.sprayed_mass == off.sprayed_mass
    assert on.removed_enthalpy == off.removed_enthalpy and on.runoff_mass == off.runoff_mass
    assert on.deep_runoff_mass == off.deep_runoff_mass and on.deep_surfaced_mass == off.deep_surfaced_mass


def test_the_molten_cascade_feeds_only_the_wall_and_stops_at_material_that_is_not_fully_molten(layered_mesh, monkeypatch):
    """The cascade respects the feed gate of fact 28(b) by construction: it feeds an element only once a death has made
    it a wall owner -- in the sheared layer by definition, which is what the gate admits -- and it stops at the first
    exposed element that is not fully molten, so molten material separated from the wall by material that is not is
    never fed (the contiguity rule of `molten_depth`). A fully molten skin 2 mm deep, a layer at 880 K beneath it and a
    fully molten core below that, at 69.8 km under Girin's closure, with the deep runoff on and off: one melt step feeds
    the skin through the cascade -- every element it fed a wall owner at that moment and fully molten -- and leaves
    every element of the buried core exactly as it was."""
    for deep in (True, False):
        b = melting_body(layered_mesh, deep_runoff=deep)
        sim, a = girin_state(b)
        r = np.linalg.norm(b.mesh.points, axis=1)
        b.solver.set_temperature(np.where(r > 0.048, 960.0, np.where(r > 0.040, 880.0, 960.0)))
        full = fully_molten(b)
        core = (r[b.mesh.tets] < 0.040).all(axis=1)                     # every node in the buried molten core
        calls, feed = [], b._feed_exposed

        def spy(exposed, h_e, h_p, b=b, feed=feed, calls=calls):
            calls.append(((b.owner_area[exposed] > 0.0).all(), full[exposed].all(), exposed.size))
            return feed(exposed, h_e, h_p)

        monkeypatch.setattr(b, "_feed_exposed", spy)
        phi0 = b.phi.copy()
        b.melt_step(sim.t, 0.5, a)
        assert calls and all(owner and molten for owner, molten, _ in calls), deep
        assert sum(n for _, _, n in calls) > 0 and b.cascade_mass > 0.0
        assert core.sum() > 100 and full[core].all()
        assert np.array_equal(b.phi[core], phi0[core]) and b.mesh.active[core].all(), deep   # the buried pool is untouched
        assert not ((b.owner_area > 0.0) & b.mesh.active & full).any()  # and the skin above the 880 K layer is gone


def test_with_the_molten_cascade_the_deep_liquid_is_still_never_sprayed(layered_mesh, monkeypatch):
    """The cascade feeds the film, never the deep account, and it runs after the spray; so with both on, a step in which
    the deep runoff moves the liquid below the conjugate depth and the cascade carries the surface through a molten
    pool still sprays only the film. Checked stage by stage on a pool 4 mm deep at 69.8 km with a millimetre of deep
    liquid already on every windward patch: the spray stage leaves the deep account bit-identical, every element death
    moves deep liquid without changing its total, and each cascade feed lands in the film alone."""
    b = melting_body(layered_mesh)
    sim, a = girin_state(b)
    molten_pool(b, depth=4e-3)
    b.m_d = np.where(b.windward, 1e-3 * b.liquid.rho * b.surface.areas, 0.0)
    seen = {"spray": [], "kill": [], "feed": []}
    spray_stage, kill, feed = b._film_and_spray, b._kill, b._feed_exposed

    def spy_spray(*args, **kw):
        before = b.m_d.copy()
        out = spray_stage(*args, **kw)
        seen["spray"].append(np.array_equal(before, b.m_d))
        return out

    def spy_kill(dead):
        before = float(b.m_d.sum())
        kill(dead)
        seen["kill"].append((before, float(b.m_d.sum())))

    def spy_feed(exposed, h_e, h_p):
        deep0, film0 = b.m_d.copy(), float(b.m_f.sum())
        out = feed(exposed, h_e, h_p)
        seen["feed"].append((np.array_equal(deep0, b.m_d), float(b.m_f.sum()) - film0, out))
        return out

    monkeypatch.setattr(b, "_film_and_spray", spy_spray)
    monkeypatch.setattr(b, "_kill", spy_kill)
    monkeypatch.setattr(b, "_feed_exposed", spy_feed)
    b.melt_step(sim.t, 0.5, a)
    assert seen["spray"] == [True] and b.sprayed_mass > 0.0             # the spray never touched the deep account
    assert seen["kill"] and all(after == pytest.approx(before, rel=1e-12) for before, after in seen["kill"])
    assert seen["feed"] and all(same and gain == pytest.approx(out, rel=1e-9) for same, gain, out in seen["feed"])
    assert b.deep_runoff_mass > 0.0 and b.m_d.sum() > 0.0 and b.cascade_mass > 0.0


def test_full_steps_with_the_molten_cascade_keep_the_books_exact(layered_mesh):
    """Six coupled macro steps from 69.8 km with a pool 6 mm deep under the nose and the physics loads: the cascade acts
    (deaths expose fully molten elements and the surface recedes through them within the step), the deep runoff and
    the spray act with it, and the energy balance and the mass stay exact by the same measures as every other melting
    test; the solver carries the film and the deep liquid, and the history's two new columns report the cascade. The
    mass is exact to the film transport's round-off, which its direct solve conserves only to its conditioning (plan
    fact 48): this 6 mm pool moves up to a quarter of a kilogram of film per step through strongly coupled patch pairs,
    and with the runoff flux of 2026-10-05 the six steps drift by up to 4.4e-12 of the body (measured over ten states
    of numpy's generator; up to 2.8e-12 before that amendment), so the mass is held to 1e-11 here."""
    b = melting_body(layered_mesh)
    sim, a = girin_state(b)
    molten_pool(b, depth=6e-3)
    b.energy0 = b.energy()
    model = heating.PhysicsHeating()
    passes = []
    for _ in range(6):
        sim.advance(0.5)
        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
        loads = model.evaluate(a, b.theta, b.surface_temperature(), b.nose_radius(), T_mean=b.mean_temperature())
        b.advance(sim.t, 0.5, loads, state=a)
        passes.append(b.last_melt["cascade_passes"])
        assert abs(b.energy_balance_residual()) < 1e-7
        assert b.mass(0.0) + b.removed_mass == pytest.approx(b.mass0, rel=1e-11)
        assert b.solver.film_mass.sum() == pytest.approx(b.m_f.sum() + b.m_d.sum())
    assert b.cascade_mass > 0.0 and max(passes) >= 2 and b.cascade_capped_steps == 0
    assert b.sprayed_mass > 0.0 and b.deep_runoff_mass > 0.0
    stats = b.melt_stats()
    assert stats["cascade_mass_kg"] == b.cascade_mass and stats["cascade_passes"] == passes[-1]


# ---------------------------------------------------------------------------------------------------------------
# Amendment of 2026-10-06: the thin branch needs a rigid substrate


def linear_field(b, i, T_wall, gradient):
    """Nodal temperatures falling linearly with depth below patch i's plane: T_wall on it, `gradient` [K/m] inward. P1
    represents a linear field exactly, so the half-liquid point along the patch's inward normal is known to round-off."""
    return T_wall + gradient * ((b.mesh.points - b.surface.centroids[i]) @ b.surface.normals[i])


def test_the_nonrigid_depth_follows_the_temperature_field_to_the_half_liquid_point(layered_mesh):
    """From the patch centre straight inward, through the elements the line actually crosses, to where the field falls
    to T_rigid (829 K on AA7075_range): read inside the element, not counted in whole elements. Inside the second prism
    layer and several elements down alike; a wall below T_rigid has no slurry under it at all."""
    b = melting_body(layered_mesh)
    i = b.i_stag
    for depth in (0.4e-3, 3e-3):
        b.solver.set_temperature(linear_field(b, i, 900.0, (900.0 - b.material.T_rigid) / depth))
        assert b.nonrigid_depth()[i] == pytest.approx(depth, rel=1e-9)
    b.solver.set_temperature(b.material.T_rigid - 1.0)                       # coherent mush or solid everywhere
    assert not b.nonrigid_depth().any()


def test_the_nonrigid_depth_weighs_each_element_by_what_is_left_of_it(layered_mesh):
    """phi_e of an element's mass is still in it; the rest went to the film, which the layer counts already, so each
    element's stretch of the line counts phi_e of its length (plan fact 58's double count, removed here)."""
    b = melting_body(layered_mesh)
    i = b.i_stag
    b.solver.set_temperature(linear_field(b, i, 900.0, (900.0 - b.material.T_rigid) / 3e-3))
    full = b.nonrigid_depth()[i]
    b.phi[:] = 0.5
    assert b.nonrigid_depth()[i] == pytest.approx(0.5 * full, rel=1e-12)


def test_the_nonrigid_depth_stops_at_the_first_rigid_point(layered_mesh):
    """Only slurry continuous with the wall counts: the line stops where the field first falls to T_rigid, whatever
    lies deeper -- as `molten_depth` stops at the first element that is not fully molten."""
    b = melting_body(layered_mesh)
    i = b.i_stag
    T = linear_field(b, i, 900.0, (900.0 - b.material.T_rigid) / 0.4e-3)
    s = -(b.mesh.points - b.surface.centroids[i]) @ b.surface.normals[i]    # each node's depth below the patch's plane
    b.solver.set_temperature(np.where(s > 1e-3, 900.0, T))                 # hot again below the two prism layers
    assert b.nonrigid_depth()[i] == pytest.approx(0.4e-3, rel=1e-9)


def slurry_under_the_surface(b, T_wall=905.0, depth=1e-3):
    """The surface just below the feed ramp (nothing fully molten, so the liquid layer is the film alone) over material
    that stays above T_rigid for `depth` below it: slurry, radially."""
    d = 0.05 - np.linalg.norm(b.mesh.points, axis=1)
    b.solver.set_temperature(T_wall - (T_wall - b.material.T_rigid) * d / depth)


def test_a_film_on_slurry_deeper_than_the_conjugate_depth_is_thick_in_the_melt_step(layered_mesh):
    """At 70 km, under Girin's closure on every windward patch (conjugate depth 0.17-0.30 mm), a 20 um film over a
    millimetre of slurry: with the rigid substrate on, every wet windward patch is thick, and only because of the
    slurry; with it off the liquid layer is the film alone and every such patch is thin."""
    assert body.MeltSettings().rigid_substrate
    seen = {}
    for on in (True, False):
        b = melting_body(layered_mesh, rigid_substrate=on)
        sim, a = girin_state(b)
        slurry_under_the_surface(b)
        flow, delta_m = flow_and_delta_m(b, a)
        assert np.isfinite(delta_m[b.windward]).all() and delta_m[b.windward].max() < 0.5e-3
        b.m_f = np.where(b.windward, 2e-5 * b.liquid.rho * b.surface.areas, 0.0)
        h_p = deep_stage_args(b, delta_m)[4]
        b._film_and_spray(sim.t, 0.5, a, h_p, flow, delta_m)
        seen[on] = b
    on, off = seen[True], seen[False]
    assert on.last_melt["thick_branch_fraction"] == 1.0 and off.last_melt["thick_branch_fraction"] == 0.0
    assert on.last_melt["slurry_thick_fraction"] == 1.0 and on.last_melt["rigid_thin_fraction"] == 0.0
    assert 0.9 < on.last_melt["nonrigid_depth_mean_mm"] < 1.1                 # the millimetre of slurry, read back
    assert on.last_melt["slurry_held_mass_kg"] == 0.0                         # Girin's closure everywhere: nothing held
    assert np.isnan(off.last_melt["slurry_thick_fraction"]) and np.isnan(off.last_melt["nonrigid_depth_mean_mm"])
    from reentry_model import spray
    wet = on.windward & (on.last_b >= spray.B_MIN)
    assert not np.isin(on.last_spray.branch[wet], [spray.BRANCH_THIN, spray.BRANCH_RAREFIED]).any()
    assert on.last_nonrigid is not None and on.last_nonrigid.shape == (on.surface.n_patches,) and off.last_nonrigid is None
    for key in ("nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg"):
        assert key in on.melt_stats()


def test_off_girins_closure_a_film_on_slurry_does_not_spray_and_one_on_a_rigid_wall_does(layered_mesh):
    """At 30 s of the 100 mm flight no patch has Girin's closure, so there is no conjugate depth. A film over slurry is
    held from the shear modes, and the mass held is recorded; the same film over a rigid wall (surface below T_rigid)
    takes the thin mode as before, and so does the film over slurry with the rigid substrate off."""
    from reentry_model import spray
    out = {}
    for label, on, T_wall in (("slurry", True, 905.0), ("rigid", True, 820.0), ("off", False, 905.0)):
        b = melting_body(layered_mesh, rigid_substrate=on)
        sim, a = girin_state(b, t=30.0)
        slurry_under_the_surface(b, T_wall=T_wall)
        flow, delta_m = flow_and_delta_m(b, a)
        assert not np.isfinite(delta_m).any()
        b.m_f = np.where(b.windward, 2e-5 * b.liquid.rho * b.surface.areas, 0.0)
        film = b.m_f.copy()
        h_p = deep_stage_args(b, delta_m)[4]
        b._film_and_spray(sim.t, 0.5, a, h_p, flow, delta_m)
        out[label] = (b, film)
    b, film = out["slurry"]
    wet = b.windward & (b.last_b >= spray.B_MIN)
    assert wet.any() and not np.isin(b.last_spray.branch[wet], [spray.BRANCH_THIN, spray.BRANCH_RAREFIED]).any()
    assert b.last_melt["slurry_held_mass_kg"] == pytest.approx(b.m_f[wet].sum() + b.last_spray.dm[wet].sum(), rel=1e-12)
    for label in ("rigid", "off"):
        b, film = out[label]
        wet = b.windward & (b.last_b >= spray.B_MIN)
        assert np.isin(b.last_spray.branch[wet], [spray.BRANCH_THIN, spray.BRANCH_RAREFIED]).all()
        assert b.sprayed_mass > 0.0
    assert out["rigid"][0].last_melt["slurry_held_mass_kg"] == 0.0


# ---------------------------------------------------------------------------------------------------------------
# Amendment of 2026-10-07: freeze-back never leaves a negative film


def test_freeze_back_of_film_and_deep_liquid_leaves_neither_negative(layered_mesh):
    """Freeze-back draws a patch's deep liquid first and its film after, capped at what the two hold. Taking the cap as
    their sum and the film's part as that sum less the deep liquid leaves the film one rounding unit below zero whenever
    the deep liquid dwarfs the film -- for about half of all such pairs -- and the solver's film-mass check rejects the
    negative: it stopped the 0.00625 s run of the 2026-10-07 time-step series at its first fine step, a film of -1.3e-25
    kg on one patch. A patch whose film and deep liquid all freeze back must be left with exactly none of either, and its
    owner element must gain exactly their mass."""
    b = melting_body(layered_mesh)
    b.solver.set_temperature(800.0)                                 # below the feed ramp: every patch's liquid freezes
    i = b.i_stag
    owner = b.surface.owner[i]
    b.phi[owner] = 0.5                                              # room for it in the owner element
    b.m_f[:], b.m_d[:] = 0.0, 0.0
    b.m_f[i], b.m_d[i] = 2.3351899704880006e-14, 1.8789852661499484e-09   # (m_f + m_d) - m_d > m_f in floating point
    assert b.m_f[i] - ((b.m_f[i] + b.m_d[i]) - b.m_d[i]) < 0.0      # the pair the old arithmetic got wrong
    before = b.m_f[i] + b.m_d[i]
    solid = (1.0 - b.material.feed_fraction(b.film_temperature())) * (b.m_f + b.m_d)
    want = np.bincount(b.surface.owner, solid, b.mesh.n_elements)
    _, _, _, h_e, h_p = deep_stage_args(b, np.full(b.surface.n_patches, np.nan))
    frozen = b._freeze_back(want, solid, h_e, h_p)
    assert (b.m_f >= 0.0).all() and (b.m_d >= 0.0).all()
    assert b.m_f[i] == 0.0 and b.m_d[i] == 0.0
    assert frozen == pytest.approx(before, rel=1e-15)
    assert (b.phi[owner] - 0.5) * b.element_mass[owner] == pytest.approx(before, rel=1e-9)
