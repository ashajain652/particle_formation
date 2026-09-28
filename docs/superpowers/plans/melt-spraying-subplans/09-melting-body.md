# Sub-plan: Task 9 — The melting body

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 4349–5703). Read `00-shared-context.md` first — this is the largest and highest-risk task in the plan and needs the most context, not the least.


> **Amended 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md`
> (measured facts 38–45 in `00-shared-context.md`). The block below takes precedence over the extracted body.

**Changes to Task 9.**

1. **`PHI_DEATH = 0.05` becomes `PHI_DEATH = 0.50`** (fact 44). A patch owner dies at half consumed, not at 95 %
   consumed. Two measured reasons: the frozen-geometry conduction error is lowest there (+1.97 K typical against
   +2.97 K at 0.05, +9.87 K against +14.83 K at peak flux, a shallow optimum at a half), and fact 42's quality
   envelope wants it if Design A is ever adopted. The `PHI_MIN = 1e-3` floor for interior elements is unchanged, as
   is fact 4's rule that the remainder joins the film and that deaths cascade within the step.
2. **Check the film can absorb it.** The film now receives half an element at a time rather than a twentieth.
   `film_blob_fraction` and `film_frozen_fraction` must stay within the ranges fact 25 records, and the stranded-film
   mass below the 0.47 % of initial mass measured there. This is an acceptance test, not an assumption.
3. **The shape consumers read the derived surface.** `nose_radius()`, `newtonian_drag()`, `reference_area()`,
   `drag_shape_factor()`, `equivalent_radius()`, `transverse_radius` and `mass_centre` take their geometry from
   `mesh.surface(derived=True)` and its `smoothed_normals()`, not from the raw element faces. `theta` and `t_hat`
   likewise. The **mass** account is untouched: it remains `Σ φ_e ρ V_e` over the **original** element volumes, and
   `φ_e` remains authoritative (spec §5, Design B).
4. **New attribute** `last_recession_residual` — the median absolute per-patch mismatch from Task 16's least-squares
   solve — and a `melt_stats()` entry for it, so how well the geometry tracks the mass is recorded per step rather
   than assumed. Fact 43: 0.03 % for a smooth recession field, 22 % for a patchy one, and which of those the real
   flight produces is unmeasured.
**Also carried here: measured fact 37 (2026-09-27), which the committed master plan is missing.** It was
regenerated into `prototype/plan3/plan.md` but never copied over `2026-09-20-melt-spraying.md`, so the extracted
body below predates it. `MeltingBody` gains the attribute `last_face_ids` and the method `on_current_surface`; add
both to the `Produces:` list, and note that `last_spray`/`last_flow` are now read through the latter:

```python
        self.last_face_ids = None                                          # in __init__, beside last_flow/last_spray
        self.last_face_ids = self.surface.face_ids.copy()                  # at step (iii), before the deaths of step (v)

    def on_current_surface(self, values, fill):
        """A per-patch array from the last flow and spray evaluation (`last_flow`, `last_spray`), carried onto the
        current surface. Element deaths later in the same melt step rebuild the surface, so the arrays no longer line
        up with it by index; patches are matched by face id instead, and a patch that was not on the surface at the
        evaluation (a face the deaths exposed) gets `fill`. Matching by index instead left every frame written after
        a step with a death -- every melting frame -- at the defaults: 0 of 104 frames of the 50 mm physics flight
        carried a wall pressure (measured 2026-09-27)."""
        n = self.surface.n_patches
        if values is None or self.last_face_ids is None or not self.last_face_ids.size:
            return np.full(n, fill, dtype=float)
        index = np.full(len(self.patch_of_face), -1, dtype=np.int64)
        index[self.last_face_ids] = np.arange(self.last_face_ids.size)
        j = index[self.surface.face_ids]
        return np.where(j >= 0, np.asarray(values, dtype=float)[np.maximum(j, 0)], fill)

```

Without it, every frame written after a step with a death — which is every melting step — carries the defaults:
**0 of 104** frames of the 50 mm physics flight had a wall pressure. Output only; the history columns move by at
most 1.0e-8 relative, against 1.6e-8 between two runs of one build.


5. **New `melt_stats()` entry** `conduction_dT_err_K` = max over live patch owners of `q L (1 − φ_e) / k`, the
   surface temperature error Design B accepts by conducting through material that is no longer there. Spec §10.3
   makes this the measured trigger for reconsidering Design A.
---

**Depends on:** Tasks 1 through 7 — all of them. This is the integration point of the whole Step 3 plan.
**Produces, for later tasks:** the `MeltingBody` class: feed/freeze mass transfer between elements, the film's temperature and energy balance, spraying integration, the element-death cascade with hand-over of mass and heat to neighbors, nose-cap fitting, and deferred-load bookkeeping. Task 10's coupled loop is built directly on this class.
**Character:** physics and bookkeeping, the densest task in the plan — its code block runs to roughly 1,200 lines, larger than any other single task.
**Read before implementing:** the master plan's own text says, in its own words, "Read facts 25 and 26 before reading `melt_step`" — treat that as a direct instruction, not a suggestion. Facts 4, 5, 25, and 26 in the shared context describe two specific wrong implementations that were tried and rejected here during prototyping (debiting only the destination element; booking the film at mixture rather than liquid-only enthalpy), each of which silently reproduced a different bug. An agent working from this task's text alone, without that history, could plausibly re-derive and reintroduce either mistake while "simplifying" the bookkeeping.
**Recommendation carried over from the build-plan review:** do not split this task further across multiple sub-agents despite its size — the feed/freeze bookkeeping, film energy balance, and element-death hand-over are tightly coupled to each other, and a review checkpoint should happen after this task and before Task 10 begins.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan, and reproduce the "read facts 25/26 first" instruction and the two rejected-approach warnings directly in that plan's text so a future implementer doesn't have to go find them in the master document.

---

### Task 9: The melting body

**Files:**
- Modify (replace): `reentry_model/body.py`
- Modify: `reentry_model/trajectory.py` (three statements), `reentry_model/aero.py` (the Mach-1 step and the shape factor)
- Create: `reentry_model/data/atdb_disc.json` (extracted from DRAMA)
- Test: `tests/test_reentry_model_melting.py` (new); `tests/test_reentry_model_aero.py` (one assertion); `tests/test_reentry_model_fenicsx.py::test_melting_run_matches_the_skfem_backend` (from Task 3) now runs in `fenicsx_env`

**Interfaces:**
- Consumes: Tasks 1–7 (`mesh.deactivate/surface/active_nodes`, `Material.feed_fraction/liquid_fraction/enthalpy/h_liquid/liquid`, `thermal` `set_fractions/element_energies/step(nodal_load=)/pinned/temperature`, `surface_flow.SurfaceFlow`, `spray.SprayModel/melt_layer/source_rows/histogram/N_BINS`, `film.lubrication/Runoff`), `heating.HeatingResult`, `trajectory.AeroState`.
- Produces: `body.REMOVAL_NAMES = ("girin", "instant")`, `SIZE_FEEDBACK_NAMES = ("current", "initial")`, `NOSE_CAP_ANGLE = 30.0`, `NOSE_CAP_FACTOR = 1.67`, `CP_MAX_NEWTONIAN = 1.84`, `PHI_MIN = 1e-3`, `PHI_DEATH = 0.05`, `NEAREST_PATCHES = 4`, `LOAD_DT_MAX = 1000.0`, `FEED_DEPTH_NAMES = ("conjugate", "all")`; `fit_sphere(points) -> (centre, radius)`; `MeltSettings(removal, runoff, demise_fraction, particles, size_feedback="current", feed_depth="conjugate")`; `MeltingBody(mesh, material, solver, mass_kg, flow=None, spray_model=None, settings=None, T0, emissivity, T_ambient, v_hat)` with `.advance(t, dt, loads, state=None)`, `.melt_step(t, dt, state)`, `.mass(t)`, `.reference_area()`, `.reference_length()` (2 R_eq, or None with `initial`), `.nose_radius()` (the windward-cap fit, or R₀ with `initial`), `.newtonian_drag() -> (C_D, projected area)` (modified Newtonian over the windward convex hull), `.drag_shape_factor()` (that value over the meshed sphere's own, so exactly 1 while intact and 2.00 at the flat limit; 1 with `initial`), `.equivalent_radius()`, `.energy()` (FEM + film), `.film_temperature()` (the surface's own: the film is thermally thin), `.film_enthalpy(h_node=None)` (per patch, the mean of the *liquid* nodal enthalpies), `.film_energy()`, `.film_frozen_fraction()`, `.film_blob_fraction()` (the share of the film deeper than its patch is wide), `.molten_depth(Te=None, max_levels=8)` (the contiguous liquid depth inward from each patch, hopping element to element and stopping at the first one that is not molten), `.liquid_layer_depth()` (that depth plus the film, which is what the branch test and the Rayleigh-Taylor criterion see), `.region_extent(selected)` (each patch's own contiguous region as an equivalent diameter 2 sqrt(A/pi), which is the domain a wave has to fit inside -- the molten region, not the facet) and `.unstable_region_extent(unstable)` (its maximum, for the reported diagnostics; facts 30, 32 and 33), `.total_applied_load()`, `.total_dropped_load()` (the deferred-load accounts the energy balance closes against), `.mean_temperature()`, `.melt_front_depth()`, `.film_thickness_max/mean()`, `.demised()`, `.energy_balance_residual()`, `.melt_stats() -> dict` (the `coupled.MELT_COLUMNS` values plus `mass_kg`), attributes `phi, m_f, surface, theta, t_hat, windward, mass0, mass_centre, transverse_radius, fitted_nose_radius, cap_nose_radius, sprayed_mass, runoff_mass, removed_mass, removed_enthalpy, n_released, source_rows, hist_n, hist_m, melt_onset, spray_onset, consumed, last_flow, last_spray, last_b, last_melt, pending_load, liquid, flow, spray, runoff`; `Body.reference_area()`/`reference_length()`/`drag_shape_factor()` in the protocol (`ConstantBody`/`ThermalBody` return None; `ThermalBody.nose_radius()` returns the sphere's radius; `ThermalBody.advance` accepts `state=None`); `Simulator.aero_state` uses the body's reference area, Knudsen length and drag shape factor, and gives a consumed body no drag; `aero.drag_coefficient(kn, ma, tables, bridging, shape_factor=1.0)` scales the continuum entry only.

The film's temperature, its re-solidification and the way every transfer is booked are the subject of measured facts 25 and 26: the film has no energy equation of its own (its *mass* goes to the solver, which carries it on the boundary nodes with the liquid capacity), it holds `enthalpy_liquid` wherever it sits, the feed and the freeze are netted into one transfer per element, each transfer books only the enthalpy difference it carries and books it at the destination, and no deferred load may move a node more than `LOAD_DT_MAX` in one step. Read those two facts before reading `melt_step`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_melting.py`:

```python
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
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_melting.py -q`
Expected: AttributeError (`body.MeltSettings`).

- [ ] **Step 3: Replace `reentry_model/body.py`**

```python
"""The body the trajectory carries, as a stepper (spec section 4). ConstantBody keeps mass and one temperature
(Step 1 behaviour, `--thermal none`); ThermalBody wraps a mesh, a material and a conduction solver and advances the
temperature field with each macro step's per-patch convective loads, keeping the heat bookkeeping SESAM reports."""
import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Body(Protocol):
    def mass(self, t) -> float: ...
    def advance(self, t, dt, loads) -> None: ...          # loads: heating.HeatingResult applied over [t - dt, t]
    def reference_area(self): ...                        # m2 drag reference area, or None for the fixed pi D^2/4
    def reference_length(self): ...                      # m length of the body Knudsen number, or None for the fixed D
    def drag_shape_factor(self): ...                     # continuum C_D of the current shape / a sphere's (1 for a sphere)
    def surface_temperature(self) -> np.ndarray: ...     # K per patch
    def mean_temperature(self) -> float: ...             # K, energy-equivalent (spec 6.4)
    def energy(self) -> float: ...                       # J stored above the material's reference temperature
    def field(self): ...                                 # nodal temperatures, or None


def sphere_mass(diameter_m, density_kgm3):
    """Mass of a solid sphere [kg]."""
    return density_kgm3 * 4.0 / 3.0 * math.pi * (diameter_m / 2.0) ** 3


@dataclass
class ConstantBody:
    mass_kg: float
    temperature_K: float = 300.0

    def mass(self, t):
        return self.mass_kg

    def reference_area(self):
        return None

    def reference_length(self):
        return None

    def drag_shape_factor(self):
        return 1.0

    def advance(self, t, dt, loads):
        return None

    def surface_temperature(self):
        return np.array([self.temperature_K])

    def mean_temperature(self):
        return self.temperature_K

    def energy(self):
        return 0.0

    def field(self):
        return None


class ThermalBody:
    """Mesh + material + ThermalSolver. `v_hat` is the direction of motion in the body frame (fixed attitude: the
    stagnation patch is the one whose normal is along v_hat, spec 13.1)."""

    def __init__(self, mesh, material, solver, mass_kg, T0=300.0, emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        self.mesh, self.material, self.solver, self.mass_kg = mesh, material, solver, mass_kg
        self.emissivity = material.emissivity if emissivity is None else float(emissivity)
        self.T_ambient = float(T_ambient)
        self.surface = mesh.surface()
        self.theta = self.surface.angles_to(v_hat)
        self.i_stag, self.i_back = self.surface.patch_toward(v_hat), self.surface.patch_toward(-np.asarray(v_hat, dtype=float))
        self.i_centre = mesh.centre_node()
        self.volume = mesh.volume()
        self.solver.setup(mesh, material, self.emissivity)
        self.solver.set_temperature(T0)
        self.energy0 = self.solver.energy()
        self.integrated_heat = 0.0       # J, integral of Q_conv dt (SESAM's integrated_heat_J)
        self.absorbed_heat = 0.0         # J, integral of (Q_conv - Q_rad) dt
        self.radiated_heat = 0.0         # J
        self.iterations = []             # Newton iterations per step
        self.last = None                 # thermal.StepResult of the last step

    def mass(self, t):
        return self.mass_kg

    def reference_area(self):
        return None

    def reference_length(self):
        """Length scale of the body Knudsen number, or None for the fixed initial diameter."""
        return None

    def drag_shape_factor(self):
        """Continuum drag of the body's shape relative to a sphere's: 1 for the sphere of Steps 1-2."""
        return 1.0

    def nose_radius(self):
        """Stagnation-point radius of curvature [m] for the heating and the surface flow: the sphere's."""
        return float(self.mesh.params.get("radius_m", 0.05))

    def advance(self, t, dt, loads, state=None):
        res = self.solver.step(dt, loads.q_conv, self.T_ambient)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res

    def surface_temperature(self):
        return self.surface.facet_mean(self.solver.temperature())

    def mean_temperature(self):
        return float(self.material.temperature_from_enthalpy(self.energy() / (self.material.rho * self.volume)))

    def energy(self):
        return self.solver.energy()

    def field(self):
        return self.solver.temperature()

    def radiated_power(self):
        return self.solver.radiated_power(self.T_ambient)

    def energy_balance_residual(self):
        """(E - E0 - absorbed heat) / absorbed heat: zero for an exact discrete balance."""
        return (self.energy() - self.energy0 - self.absorbed_heat) / self.absorbed_heat if self.absorbed_heat else 0.0

    def surface_stats(self):
        T, Tf = self.solver.temperature(), self.surface_temperature()
        return {"surface_T_max_K": float(Tf.max()), "surface_T_min_K": float(Tf.min()),
                "surface_T_mean_K": float((Tf * self.surface.areas).sum() / self.surface.area),
                "T_stagnation_K": float(Tf[self.i_stag]), "T_back_K": float(Tf[self.i_back]), "T_centre_K": float(T[self.i_centre])}


SIZE_FEEDBACK_NAMES = ("current", "initial")
FEED_DEPTH_NAMES = ("conjugate", "all")     # how deep melt may leave an element for the film (Step 3, 2026-09-22)
CP_MAX_NEWTONIAN = 1.84          # only the ratio to the meshed sphere's own value is used, so this cancels
NOSE_CAP_ANGLE = 30.0            # deg: the windward cap fitted for the nose radius (depth (1 - cos 30 deg) R_t behind the front)
NOSE_CAP_FACTOR = 1.67           # a flat face of radius R_t heats like a sphere of 1.67 R_t (its stagnation velocity gradient is
                                 # ~0.6 x a sphere's of the same radius, Boison & Curtiss 1959): the cap on the fitted radius
PHI_MIN = 1.0e-3                 # element fraction kept by elements that do not own a patch (they cannot die: no cavities)
NEAREST_PATCHES = 4              # patches that receive an interior element's liquid (area-weighted)
LOAD_DT_MAX = 1000.0             # K, the most a deferred melt load may move one node in one macro step
PHI_DEATH = 0.05                 # a patch owner below this fraction dies (its remainder goes to the film): keeps the surface
                                 # nodes' thermal mass above 5 % of an element's, which the Newton iteration needs (measured
                                 # 2026-09-20: owners at 1e-3 made surface nodes swing by hundreds of K between iterates)
REMOVAL_NAMES = ("girin", "instant")


@dataclass
class MeltSettings:
    removal: str = "girin"           # girin: film + runoff + spraying; instant: liquid removed as it forms (verification device)
    runoff: bool = True
    demise_fraction: float = 0.01    # the run ends when the mass falls below this fraction of the initial mass
    particles: bool = True           # keep the source-table rows
    size_feedback: str = "current"    # current: Kn on the equivalent diameter of the remaining mass and the nose radius fitted
                                      # to the windward cap; initial: D0 and R0 throughout (SESAM's convention, for the devices)
    feed_depth: str = "conjugate"     # conjugate: only liquid the gas shear reaches leaves its element (fact 28); all: any depth

    def __post_init__(self):
        if self.removal not in REMOVAL_NAMES:
            raise ValueError("removal must be one of {}, got {!r}".format(REMOVAL_NAMES, self.removal))
        if self.size_feedback not in SIZE_FEEDBACK_NAMES:
            raise ValueError("size_feedback must be one of {}, got {!r}".format(SIZE_FEEDBACK_NAMES, self.size_feedback))
        if self.feed_depth not in FEED_DEPTH_NAMES:
            raise ValueError("feed_depth must be one of {}, got {!r}".format(FEED_DEPTH_NAMES, self.feed_depth))
        if not 0.0 < self.demise_fraction < 1.0:
            raise ValueError("demise_fraction must be within (0, 1)")


def fit_sphere(points):
    """Algebraic least-squares sphere through `points` (n, 3): (centre, radius)."""
    x = np.asarray(points, dtype=float)
    A = np.column_stack([2.0 * x, np.ones(len(x))])
    p = np.linalg.lstsq(A, (x * x).sum(axis=1), rcond=None)[0]
    c = p[:3]
    return c, float(math.sqrt(max(p[3] + c @ c, 0.0)))


class MeltingBody(ThermalBody):
    """ThermalBody with melting (spec Step 3 sections 6, 8-11): element fractions phi_e, the film per patch, the
    gas-side surface flow, spraying, element death and the mass/area/energy accounting.

    Per macro step (`advance`): the conduction step with the deferred melt loads of the previous step -> melt step:
    (i) every element's liquid inventory f_feed(T_e) phi_e rho V_e becomes film on its patches (owners) or the
    nearest patch (interior elements keep PHI_MIN so no cavity can open), netted against (iv) below; what leaves is
    the element's molten part, at the enthalpy that part carries, so the melt front moves at the energy-limited rate
    whatever the element size; (ii) surface flow, delta_m, lubrication, runoff transport; (iii) spraying and release;
    (iv) re-solidification: the fraction 1 - f_feed(T_patch) of each patch's film returns to its owner element
    (raising phi_e, capped at 1) once the patch falls back through the ramp -- the mirror of the feed rule, netted
    with it so that no element both melts and freezes in one step; (v) patch owners at phi <= PHI_DEATH die: the
    mesh's active set, the surface, the film (handed to the nearest surviving patch) and the solver's fractions are
    refreshed. With `removal = "instant"` the liquid leaves the body at h_liquid instead, and the difference to the
    element's own enthalpy is a nodal load on its nodes over the next step.

    The film has the surface's own temperature, because it is thermally thin: q b / k_l = 0.22 K across a 10 um film
    at 2 MW/m2 (22 K even at 1 mm) and b^2/alpha = 0.3 ms against the 0.5 s macro step. Rather than give it an energy
    equation, its *mass* is handed to the solver (`set_film_mass`), which weighs it with the liquid heat capacity, so
    the film rides the boundary nodes at the surface temperature by construction; it heats, cools and re-solidifies
    with the surface, and the droplets carry away whatever superheat the surface has. A film is liquid by
    construction, so it holds the *liquid* enthalpy h(T) + L_f (1 - f_l) wherever it sits (`film_enthalpy`): feeding
    it costs the latent heat its mass has not yet paid, which is what keeps a body that merely sits on the melting
    ramp from turning into film for free. `film_energy` sums that enthalpy over the patches -- identical to the nodal
    sum the solver's capacity matrix carries -- and every transfer -- feed, freeze-back, runoff between patches at
    different temperatures, the hand-over when a patch dies -- books the enthalpy difference it carries, at the
    patch it arrives on (`melt_step`), so the coupled balance stays exact (measured 1e-10). `removal = "instant"`
    removes the liquid inventory of every element as it forms, without film, runoff or spraying: the lumped-melting
    device that reproduces SESAM's Q/L_f law (SESAM hollows the sphere at fixed outer geometry, measured
    2026-09-20, so the geometry is kept until elements die).
    mass(t) = sum phi rho V + sum m_f; the drag reference area is the current surface's projection on the flight
    direction (pi R^2 while intact). Size feedback (decided 2026-09-21, replacing spec section 10's fixed R0): the body
    Knudsen number uses the equivalent diameter of the remaining mass, and the stagnation radius for the heating and
    the surface flow is fitted to the current windward cap -- a least-squares sphere through the patch centroids within
    (1 - cos 30 deg) R_t of the front-most point (R_t the transverse radius about the mass centre), bounded to
    [0.1 R_t, NOSE_CAP_FACTOR R_t] (the front erodes fastest and flattens, which lowers the stagnation heating as R^-1/2;
    a flat face heats like a sphere of 1.67 x its radius). The sphere drag tables are kept. `size_feedback = "initial"`
    keeps D0 and R0 (SESAM's convention, used by the verification devices)."""

    def __init__(self, mesh, material, solver, mass_kg, flow=None, spray_model=None, settings=None, T0=300.0,
                 emissivity=None, T_ambient=0.0, v_hat=(1.0, 0.0, 0.0)):
        if not material.melts or material.liquid is None:
            raise ValueError("MeltingBody needs a material with a latent heat and liquid properties (AA7075 or AA7075_range)")
        super().__init__(mesh, material, solver, mass_kg, T0, emissivity, T_ambient, v_hat)
        from . import film as film_mod, spray as spray_mod, surface_flow
        self.settings = settings or MeltSettings()
        self.flow = flow or surface_flow.SurfaceFlow()
        self.spray = spray_model or spray_mod.SprayModel(material.liquid)
        self._film_mod = film_mod
        self.v_hat = np.asarray(v_hat, dtype=float) / np.linalg.norm(v_hat)
        self.liquid = material.liquid
        self.element_mass = material.rho * mesh.element_volumes()          # kg at phi = 1
        self.phi = np.ones(mesh.n_elements)
        self.mass0 = float(self.element_mass.sum())
        self.m_f = np.zeros(self.surface.n_patches)
        self.pending_load = np.zeros(mesh.n_nodes)                          # J, deferred melt energy for the next step
        self.sprayed_mass = self.runoff_mass = self.removed_mass = self.removed_enthalpy = self.frozen_mass = 0.0
        self.n_released = 0.0
        self.source_rows = []
        self.hist_n, self.hist_m = np.zeros(spray_mod.N_BINS), np.zeros(spray_mod.N_BINS)
        self.last_flow = self.last_spray = None
        self.last_melt = {"n_dead": 0, "runoff_substeps": 0, "released_mass": 0.0, "n_released": 0.0, "feed_mass": 0.0}
        self.melt_onset = self.spray_onset = None
        self.consumed = False
        self._refresh_geometry()
        self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation

    # -- geometry -----------------------------------------------------------------------------------------------
    def _refresh_geometry(self):
        self.surface = self.mesh.surface()
        self.theta = self.surface.angles_to(self.v_hat)
        self.t_hat = self.surface.tangent_from(self.v_hat)
        self.windward = self.theta <= 0.5 * np.pi
        self.i_stag, self.i_back = self.surface.patch_toward(self.v_hat), self.surface.patch_toward(-self.v_hat)
        self.runoff = self._film_mod.Runoff(self.surface, self.mesh.points, self.windward) if self.settings.runoff else None
        owners = self.surface.owner
        self.owner_area = np.bincount(owners, self.surface.areas, self.mesh.n_elements)      # patch area per owner element
        from scipy.spatial import cKDTree
        self._patch_tree = cKDTree(self.surface.centroids)
        self.patch_of_face = np.full(len(self.mesh._face_nodes), -1, dtype=np.int64)     # face id -> patch index (-1: not a patch)
        self.patch_of_face[self.surface.face_ids] = np.arange(self.surface.n_patches)
        self._fit_nose()

    def _fit_nose(self):
        """Mass centre, transverse radius and the windward-cap radius of the current surface (module docstring)."""
        m, s, v = self.mesh, self.surface, self.v_hat
        w = self.phi * m.element_volumes()
        self.mass_centre = (w[:, None] * m.points[m.tets].mean(axis=1)).sum(axis=0) / max(w.sum(), 1e-300)
        X = m.points[np.unique(s.faces)] - self.mass_centre
        self.transverse_radius = float(np.sqrt(np.maximum(np.linalg.norm(X, axis=1) ** 2 - (X @ v) ** 2, 0.0)).max())
        x = (s.centroids - self.mass_centre) @ v
        band = x >= x.max() - (1.0 - math.cos(math.radians(NOSE_CAP_ANGLE))) * self.transverse_radius
        r0 = float(m.params.get("radius_m", 0.05))
        if band.sum() >= 4 and self.transverse_radius > 0.0:
            _, r = fit_sphere(s.centroids[band])
            self.fitted_nose_radius = r
            self.cap_nose_radius = float(min(max(r, 0.1 * self.transverse_radius), NOSE_CAP_FACTOR * self.transverse_radius))
        else:
            self.fitted_nose_radius = self.cap_nose_radius = r0

    def nose_radius(self):
        if self.settings.size_feedback == "initial":
            return float(self.mesh.params.get("radius_m", 0.05))
        return self.cap_nose_radius

    def newtonian_drag(self):
        """Modified-Newtonian drag coefficient of the current windward silhouette: sum C_p,max (n.v)^3 A / sum (n.v) A
        over the convex hull of the surface (0.92 for a sphere, 1.84 for a flat face, with C_p,max = 1.84).

        The hull, not the raw facets: the staircase left by element death scatters the normals and biases the integral
        low and non-monotonically (0.79-1.07 against the hull's smooth 0.92-1.73, measured 2026-09-21), and the flow
        sees the silhouette -- a shallow cavity recovers roughly the stagnation pressure at its mouth. It is therefore
        an upper bound where the face is cratered."""
        from scipy.spatial import ConvexHull
        points = self.mesh.points[np.unique(self.surface.faces)]
        try:
            hull = ConvexHull(points)
        except Exception:                                        # degenerate remnant: fall back to the sphere
            return CP_MAX_NEWTONIAN / 2.0, self.surface.projected_area(self.v_hat)
        tri = points[hull.simplices]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        area = 0.5 * np.linalg.norm(n, axis=1)
        n = n / (2.0 * np.maximum(area, 1e-300))[:, None]
        inward = np.einsum("ij,ij->i", n, tri.mean(axis=1) - points.mean(axis=0)) < 0.0
        n[inward] *= -1.0
        c = n @ self.v_hat
        w = c > 0.0
        projected = float((area[w] * c[w]).sum())
        if projected <= 0.0:
            return CP_MAX_NEWTONIAN / 2.0, self.surface.projected_area(self.v_hat)
        return float((CP_MAX_NEWTONIAN * c[w] ** 3 * area[w]).sum() / projected), projected

    def drag_shape_factor(self):
        """The continuum drag of the current shape relative to a sphere's, for aero.drag_coefficient: exactly 1 while
        the body is intact (it is normalised by the meshed sphere's own Newtonian value, so the mesh discretisation
        error cancels and the ATDB sphere table is reproduced bit for bit) and 2.00 at the flat-face limit, where the
        ATDB disc entry sits. `size_feedback = "initial"` pins it to 1 (SESAM's convention)."""
        if self.settings.size_feedback == "initial":
            return 1.0
        return self.newtonian_drag()[0] / self.newtonian_drag_sphere

    def reference_length(self):
        if self.settings.size_feedback == "initial":
            return None
        return 2.0 * self.equivalent_radius()

    def mass(self, t):
        return float((self.phi * self.element_mass).sum() + self.m_f.sum())

    def reference_area(self):
        """Drag reference area: the hull's projected area, paired with the C_D of drag_shape_factor (the raw surface's
        projection counts forward-facing patches inside a crater and is 2-4 % larger)."""
        if self.settings.size_feedback == "initial":
            return self.surface.projected_area(self.v_hat)
        return self.newtonian_drag()[1]

    def equivalent_radius(self):
        return (3.0 * self.mass(0.0) / (4.0 * math.pi * self.material.rho)) ** (1.0 / 3.0)

    def energy(self):
        return self.solver.energy() + self.film_energy()

    def mean_temperature(self):
        m = self.mass(0.0)
        return float(self.material.temperature_from_enthalpy(self.energy() / m)) if m > 0.0 else 0.0

    def melt_front_depth(self):
        """Deepest molten point [m]: the distance from the original surface of the innermost element with f_l > 0."""
        f_l = self.material.liquid_fraction(self.solver.temperature())[self.mesh.tets].max(axis=1)
        molten = self.mesh.active & (f_l > 0.0)
        if not molten.any():
            return 0.0
        r = np.linalg.norm(self.mesh.points[self.mesh.tets[molten]].mean(axis=1), axis=1)
        return float(self.mesh.params.get("radius_m", r.max()) - r.min())

    # -- the step --------------------------------------------------------------------------------------------------
    def advance(self, t, dt, loads, state=None):
        # A deferred melt load is energy the transferred mass delivered to the nodes it arrived at; a node that has
        # since melted away has nothing to heat with it, so no node is asked to move more than LOAD_DT_MAX in one
        # step and the remainder waits (it is applied as soon as the node has the capacity, dropped into Q_dropped
        # if the node dies first, and counted either way -- the balance sees `pending_load` and `Q_dropped` alike).
        room = LOAD_DT_MAX * self.solver.nodal_capacity()
        applied = np.clip(self.pending_load, -room, room)
        self.pending_load = self.pending_load - applied
        res = self.solver.step(dt, loads.q_conv, self.T_ambient, nodal_load=applied / dt)
        self.integrated_heat += res.Q_conv * dt
        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt
        self.radiated_heat += res.Q_rad * dt
        self.iterations.append(res.iterations)
        self.last = res
        self._applied_total = getattr(self, "_applied_total", 0.0) + res.Q_extra * dt
        self._dropped_total = getattr(self, "_dropped_total", 0.0) + res.Q_dropped * dt
        self.melt_step(t, dt, state)

    def melt_step(self, t, dt, state):
        """Feed, film transport, spraying, re-solidification and element death (class docstring).

        Energy is moved between three places -- the finite-element solid, the film, and the outside world -- and every
        transfer is booked by the same rule: mass that moves at one temperature books nothing, and mass that arrives
        somewhere colder or hotter than it left books the difference it carries, as a deferred load on the nodes it
        arrives at. So a transfer of m kilograms from enthalpy h_src to h_dst changes the accounted energy by
        m (h_dst - h_src) and applies m (h_src - h_dst) to the destination: the balance closes exactly, and nothing
        larger than the difference ever touches the temperature field. Booking the two halves separately -- the
        solid's m h_src on its element's nodes, the film's m h_dst on its patch's nodes -- closes the balance just as
        exactly and wrecks the field: the same mass is spread by 1/4 over four nodes on one side and by 1/3 over
        three on the other, which left +-m h/12 on every surface node and swung the body to 1618 K and -1202 K in two
        macro steps (measured 2026-09-22). The enthalpies are the ones the two energy functionals actually use: the
        element's is the mean of h(T_i) over its four nodes, the film's the mean of the *liquid* h over its patch's
        three, and what leaves an element is its molten part, at the enthalpy that part carries."""
        mat, liq, s = self.material, self.liquid, self.settings
        from . import spray as spray_mod
        T = self.solver.temperature()
        tets = self.mesh.tets
        # The surface flow is evaluated once per step, here rather than inside the film step, because the feed needs
        # Girin's conjugate melt-layer depth: the gas shear penetrates the liquid only that far, so only liquid within
        # that depth of the wall can be carried away, and material deeper than it keeps its melt until the surface has
        # receded to it. Feeding from any depth let the interior drain through an intact skin -- by 82 s of the 100 mm
        # flight the core was being consumed faster than the outermost shell (measured 2026-09-22).
        flow = delta_m = None
        if state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
            flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
            delta_m, _ = spray_mod.melt_layer(flow, liq)
        h_node = mat.enthalpy(T)
        h_e = h_node[tets].mean(axis=1)                          # as solver.energy() weighs the solid
        h_liq = mat.enthalpy_liquid(T)
        h_p = self.film_enthalpy(h_liq)                          # as film_energy() weighs the film: liquid, with L_f
        # (i) feed and (iv) re-solidification, netted. The feed hands the molten fraction of what each element still
        # holds to the film (f_feed of phi_e, so an element that is a quarter molten hands over a quarter of its
        # remainder); the mirror rule hands back the fraction of the film that has fallen below the ramp. What stops
        # the feed from eating a body that merely sits on the ramp is not the rule but the energy: film is liquid and
        # carries L_f wherever it sits (`Material.enthalpy_liquid`), so every kilogram fed debits the latent heat it
        # has not paid and the surface falls back to the solidus. Book the film at the mixture enthalpy instead and
        # melting is free on the ramp -- a 100 mm sphere turned entirely to film on a quarter of its latent heat
        # (measured 2026-09-22). The two directions are netted per element because they are one equilibrium seen
        # from opposite sides: run separately they cycled 3 % of the body's mass through the film every step with no
        # net effect, pinning the surface at T_feed and paying Newton iterations for it.
        fn = mat.feed_fraction(T)[tets]
        f = fn.mean(axis=1) * self.mesh.active
        if s.feed_depth == "conjugate" and s.removal != "instant":
            f = f * self._shear_reaches(delta_m)
        h_hot = (fn * h_node[tets]).mean(axis=1) * self.mesh.active     # what the molten part carries, per kg of element
        cap = np.where(self.owner_area > 0.0, self.phi, np.maximum(self.phi - PHI_MIN, 0.0))
        gross = np.minimum(f * self.phi, cap) * self.element_mass
        solid = (1.0 - mat.feed_fraction(self.film_temperature())) * self.m_f if s.removal != "instant" else np.zeros(0)
        want = np.bincount(self.surface.owner, solid, self.mesh.n_elements) if solid.size else np.zeros(self.mesh.n_elements)
        net = gross - want
        fed, want = np.maximum(net, 0.0), np.maximum(-net, 0.0)
        # the enthalpy the fed mass carries: the *molten* part of the element, not its mean. The mass that leaves sits
        # at the nodes that are above the ramp, so it carries h(T_i) there -- and what stays behind is the colder
        # rest, which is why the element must be debited the difference. Debit only the destination and the element
        # keeps a melt fraction it no longer has and feeds it again next step: the surface then melted three times
        # faster than the heat allowed and the body was gone in ten steps (measured 2026-09-22).
        with np.errstate(divide="ignore", invalid="ignore"):
            fed_h = np.where(f > 0.0, fed * h_hot / np.where(f > 0.0, f, 1.0), 0.0)
        self.phi = self.phi - fed / self.element_mass
        feed_mass = float(fed.sum())
        if feed_mass > 0.0 and self.melt_onset is None:
            self.melt_onset = t
        if s.removal == "instant":
            self.removed_mass += feed_mass
            self.removed_enthalpy += feed_mass * mat.h_liquid
            self._defer_to_elements(fed * (h_e - mat.h_liquid))  # it leaves the body at the liquidus, not at h_e
            released = feed_mass
        else:
            self._defer_to_elements(fed * h_e - fed_h)           # the element keeps only what the melt left behind
            delta, carried = self._add_to_film(fed, fed_h)
            self._defer_to_patches(carried - delta * h_p)        # ... and the melt arrives at its patch's temperature
            released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
        frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
        # (v) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed
        n_dead = 0
        while not self.consumed:
            dead = np.flatnonzero((self.owner_area > 0.0) & self.mesh.active & (self.phi <= PHI_DEATH))
            if dead.size == 0:
                break
            rest = self.phi[dead] * self.element_mass[dead]
            if s.removal == "instant":
                self.removed_mass += float(rest.sum())
                self.removed_enthalpy += float(rest.sum()) * mat.h_liquid
                self._defer_to_elements(rest * (h_e[dead] - mat.h_liquid), dead)
            else:
                extra, extra_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
                extra[dead], extra_h[dead] = rest, rest * h_e[dead]
                delta, carried = self._add_to_film(extra, extra_h)
                self._defer_to_patches(carried - delta * h_p)
            self.phi[dead] = 0.0
            before = float((self.m_f * h_p).sum())
            self._kill(dead)                                     # the film of a dead patch moves to another patch...
            if self.m_f.size:                                    # ... which is at its own temperature
                h_p = self.film_enthalpy(h_liq)
                self._spread_to_patches(before - float((self.m_f * h_p).sum()), self.m_f)
            n_dead += int(dead.size)
        self.solver.set_fractions(self.phi)
        self.solver.set_film_mass(self._film_nodal())
        self.frozen_mass += frozen
        self.last_melt.update({"n_dead": n_dead, "released_mass": released, "feed_mass": feed_mass, "frozen_mass": frozen})

    def _freeze_back(self, want, solid, h_e, h_p):
        """Give each element back the `want` kilograms of film that have fallen below the feed ramp, drawn from its
        own patches in proportion to what each wants to give (`solid`). An element can only recover film that is
        still there -- what ran off or sprayed away is gone -- and phi_e is capped at 1, because a patch may hold
        what its neighbours' runoff delivered and the mesh cannot grow a layer outside itself, so film with nowhere
        to go simply stays film (`film_frozen_fraction` records how much)."""
        if not self.m_f.size or not want.any():
            return 0.0
        owner = self.surface.owner
        have = np.bincount(owner, solid, self.mesh.n_elements)
        room = np.maximum(1.0 - self.phi, 0.0) * self.element_mass
        with np.errstate(divide="ignore", invalid="ignore"):
            share = np.where(have > 0.0, np.minimum(want, room) / np.where(have > 0.0, have, 1.0), 0.0)
        taken = np.minimum(solid * share[owner], self.m_f)
        if not taken.any():
            return 0.0
        per_element = np.bincount(owner, taken, self.mesh.n_elements)
        self.phi = np.clip(self.phi + per_element / self.element_mass, 0.0, 1.0)
        self.m_f = self.m_f - taken
        # the film arrives in the element it froze onto, which is colder than the surface it left
        self._defer_to_elements(np.bincount(owner, taken * h_p, self.mesh.n_elements) - per_element * h_e)
        return float(taken.sum())

    def _defer_to_elements(self, energy, elements=None):
        """Book a per-element energy [J] on the elements' nodes, to be applied over the next step."""
        e = np.asarray(energy, dtype=float)
        if not e.any():
            return
        nodes = self.mesh.tets if elements is None else self.mesh.tets[elements]
        np.add.at(self.pending_load, nodes.ravel(), np.repeat(e / 4.0, 4))

    def _spread_to_patches(self, energy, weights):
        """Book one energy [J] over the patches that received mass, in proportion to how much each received. Used
        where the transfer's per-edge detail is not carried back (runoff, and the hand-over when a patch dies): the
        total is exact and it lands on the arriving liquid, which is where the difference is released."""
        total = float(weights.sum())
        if total > 0.0 and energy != 0.0:
            self._defer_to_patches(energy * weights / total)

    def _defer_to_patches(self, energy):
        """Book a per-patch energy [J] on the patches' nodes, to be applied over the next step."""
        e = np.asarray(energy, dtype=float)
        if e.any():
            np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))

    def _film_nodal(self):
        """The film's mass per node [kg]: each patch's film shared over its three nodes. The film's thermal state
        lives on the nodes because that is where the solver's capacity and the solver's enthalpy live."""
        if not self.m_f.size:
            return np.zeros(self.mesh.n_nodes)
        return np.bincount(self.surface.faces.ravel(), np.repeat(self.m_f / 3.0, 3), self.mesh.n_nodes)



    def _add_to_film(self, fed, carried):
        """Distribute each element's freed liquid to the film: owners to their own patches by area, interior elements
        to the NEAREST_PATCHES nearest patches by area. `carried` is the enthalpy that liquid takes with it [J per
        element]; it rides the same weights, so the caller knows what arrived on each patch. Returns the per-patch
        increment [kg] and the per-patch enthalpy carried in [J]."""
        before, got = self.m_f.copy(), np.zeros_like(self.m_f)
        k = np.flatnonzero(fed > 0.0)
        if k.size == 0:
            return np.zeros_like(self.m_f), got
        owner = self.owner_area[k] > 0.0
        # owners: shared by patch area; interior elements: the nearest patch
        if owner.any():
            ko = k[owner]
            share, share_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
            share[ko], share_h[ko] = fed[ko] / self.owner_area[ko], carried[ko] / self.owner_area[ko]
            self.m_f += share[self.surface.owner] * self.surface.areas
            got += share_h[self.surface.owner] * self.surface.areas
        if (~owner).any():                                    # interior elements: the NEAREST_PATCHES nearest patches, by area
            ki = k[~owner]
            kk = min(NEAREST_PATCHES, self.surface.n_patches)
            _, near = self._patch_tree.query(self.mesh.points[self.mesh.tets[ki]].mean(axis=1), k=kk)
            near = near.reshape(len(ki), kk)
            w = self.surface.areas[near] / self.surface.areas[near].sum(axis=1, keepdims=True)
            np.add.at(self.m_f, near.ravel(), (fed[ki][:, None] * w).ravel())
            np.add.at(got, near.ravel(), (carried[ki][:, None] * w).ravel())
        return self.m_f - before, got

    def _film_and_spray(self, t, dt, state, h_p, flow=None, delta_m=None):
        liq, s, mat = self.liquid, self.settings, self.material
        from . import spray as spray_mod
        if state is None or self.m_f.sum() <= 0.0 or flow is None:
            self.last_flow = self.last_spray = None
            return 0.0
        areas = self.surface.areas
        # the depth of liquid under each patch: its film plus the contiguous molten material beneath it. The branch
        # tests below ask whether the gas shear reaches the bottom of the liquid, so they belong on this; the fluxes
        # and the release stay on the film, which is the mass that can actually move (decided 2026-09-22).
        molten = self.molten_depth()
        # (ii) lubrication and runoff
        n_sub, moved = 0, 0.0
        if self.runoff is not None:
            q_of_b = lambda bb: self._film_mod.lubrication(flow.tau, flow.G, bb, delta_m, liq.mu, b_layer=bb + molten)[1]
            before = self.m_f
            self.m_f, n_sub, moved = self.runoff.transport(self.m_f, q_of_b, self.t_hat, liq.rho, areas, dt)
            self.runoff_mass += moved                                              # mass that arrived on another patch
            # film that runs to a colder patch takes its enthalpy with it and arrives at that patch's temperature;
            # the difference is released where it lands (the per-edge detail is not carried back, so it is shared
            # over the patches that gained film, which is exact in total and second order in the attribution)
            d = self.m_f - before
            self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
        b = self.m_f / (liq.rho * areas)
        layer = b + molten
        v_s, q, _, thick = self._film_mod.lubrication(flow.tau, flow.G, b, delta_m, liq.mu, b_layer=layer)
        # the molten surface: the area over which a wave could form at all, either wetted by the film or molten in its
        # own right. Its contiguous extent is what every mode's wavelength is measured against (2026-09-24).
        wetted = (b >= spray_mod.B_MIN) | (molten > 0.0)
        extent = self.region_extent(wetted)
        # Girin & Kopyt's W sin(Theta): the deceleration normal to the film, which on this body is W cos(phi). It is the
        # whole deceleration at the stagnation point and vanishes at the equator, where W lies in the surface.
        w_n = flow.deceleration * np.maximum(np.cos(self.theta), 0.0)
        # (iii) spraying
        res = self.spray.evaluate(flow, state, b, delta_m, v_s, self.windward, dt, areas, self.m_f,
                                  self.transverse_radius, b_layer=layer, extent=extent, deceleration_n=w_n)
        released = float(res.dm.sum())
        if released > 0.0:
            if self.spray_onset is None:
                self.spray_onset = t
            if s.particles:
                self.source_rows += spray_mod.source_rows(t, state.h, state.V, self.theta, self.surface.centroids, self.t_hat, flow, res, b, liq, state)
            hn, hm = spray_mod.histogram(res.r, res.dn, res.dm)
            self.hist_n += hn
            self.hist_m += hm
            self.removed_enthalpy += float((res.dm * h_p).sum())   # the droplets keep the surface's superheat, and go
            self.m_f = self.m_f - res.dm
            self.sprayed_mass += released
            self.removed_mass += released
            self.n_released += float(res.dn.sum())
        self.m_f[self.m_f < 1e-30] = 0.0                       # no denormal films (they made 0/0 coefficients in the runoff)
        self.last_flow, self.last_spray, self.last_b = flow, res, b
        # The Rayleigh-Taylor criterion is reported, never applied, so report it usefully: a body-level "any patch"
        # boolean says nothing about how much melt is involved or whether the unstable wave even fits on the nose.
        # `rt_mass_fraction` is the share of the film on unstable patches; `rt_wavelength_over_nose` is the shortest
        # unstable wavelength against the nose diameter -- above 1 the mode cannot develop on the cap at all.
        rt_on, rt_lam, rt_tau = spray_mod.rayleigh_taylor(w_n, layer, liq)
        # the depth criterion above only says the pool is deep enough for the fastest mode to see it as deep
        # (W h^2 rho > 3 Sigma is h > lambda*/2pi); whether a wave fits across the pool is a second, lateral question,
        # answered against the molten region each patch belongs to rather than the body's largest one
        rt_bound, rt_bl, rt_bt = spray_mod.rayleigh_taylor_bounded(w_n, layer, extent, liq)
        rt_extent = float(extent[np.asarray(rt_on) & (self.m_f > 0.0)].max()) if (np.asarray(rt_on) & (self.m_f > 0.0)).any() else 0.0
        wet = self.m_f > 0.0
        rt_sel = np.asarray(rt_on) & wet
        rt_mass = float(self.m_f[rt_sel].sum() / self.m_f.sum()) if self.m_f.sum() > 0.0 else 0.0
        with np.errstate(invalid="ignore"):
            rt_fit = float(np.nanmin(rt_lam[rt_sel]) / (2.0 * self.nose_radius())) if rt_sel.any() else float("nan")
        r = res.r[res.dm > 0.0]
        self.last_melt.update({
            "runoff_substeps": n_sub, "n_released": float(res.dn.sum()), "regime_fractions": flow.closure_fractions(areas, self.windward),
            "kn_body": flow.kn_body, "re_shock": flow.re_shock, "flow_branch": float(flow.branch),
            "kn_local_stag": float(flow.kn_local[self.i_stag]), "p_w_stag": float(flow.p_w[self.i_stag]),
            "phi_sonic_deg": math.degrees(flow.phi_sonic) if np.isfinite(flow.phi_sonic) else float("nan"),
            "theta_cr_deg": float(np.degrees(self.theta[res.unstable].min())) if res.unstable.any() else float("nan"),
            "spraying_area_m2": float(areas[res.unstable].sum()), "rt_active": float(res.rt_active),
            "r_median_um": float(np.median(r) * 1e6) if r.size else float("nan"), "r_max_um": float(r.max() * 1e6) if r.size else float("nan"),
            "film_thickness_max_mm": float(b.max() * 1e3),
            "rt_mass_fraction": rt_mass, "rt_wavelength_over_nose": rt_fit,
            "rt_growth_ms": float(np.nanmedian(np.asarray(rt_tau)[rt_sel]) * 1e3) if rt_sel.any() else float("nan"),
            "rt_region_mm": rt_extent * 1e3,
            "rt_bounded_fraction": float(self.m_f[np.asarray(rt_bound) & (self.m_f > 0.0)].sum() / self.m_f.sum())
            if self.m_f.sum() > 0.0 else 0.0,
            "rt_bounded_growth_ms": float(np.nanmedian(np.asarray(rt_bt)[np.asarray(rt_bound) & (self.m_f > 0.0)]) * 1e3)
            if (np.asarray(rt_bound) & (self.m_f > 0.0)).any() else float("nan"),
            "spray_growth_ms": float(np.nanmedian(res.growth[res.dm > 0.0]) * 1e3)
            if res.growth is not None and (res.dm > 0.0).any() and np.isfinite(res.growth[res.dm > 0.0]).any() else float("nan"),
            "molten_depth_max_mm": float(molten.max() * 1e3),
            "molten_depth_mean_mm": float((self.m_f * molten).sum() / self.m_f.sum() * 1e3) if self.m_f.sum() > 0.0 else 0.0,
            "delta_m_mean_um": self._mean_conjugate_um(delta_m),
            "thick_branch_fraction": float(thick[self.windward & (self.m_f > 0.0)].mean())
            if (self.windward & (self.m_f > 0.0)).any() else float("nan")})
        return released

    def _kill(self, dead):
        old_surface, old_m_f = self.surface, self.m_f
        gone, new = self.mesh.deactivate(dead)
        if self.mesh.n_active == 0:                              # nothing left: the run ends (demise) without a surface
            self.consumed = True                                 # the film still on it leaves with the body
            self.removed_mass += float(self.m_f.sum())
            self.removed_enthalpy += float((self.m_f * self.film_enthalpy()).sum())
            self.m_f = np.zeros(0)
            return
        self._refresh_geometry()
        m_f = np.zeros(self.surface.n_patches)
        keep = self.patch_of_face[old_surface.face_ids]
        kept = keep >= 0
        m_f[keep[kept]] += old_m_f[kept]
        lost = ~kept & (old_m_f > 0.0)
        if lost.any():
            # The film on a vanished patch has not moved: the surface receded *into* it, so it now rests on the faces
            # of that same element which have just become boundary -- the face on its inward side, whose neighbour is
            # the next layer in, and the walls of the pit its lateral neighbours expose -- shared by their areas.
            # Handing it to the nearest surviving centroid instead sends it sideways to whatever triangle happens to
            # be closest, which put a third of a collapsing body's film onto one 0.6 mm2 face (measured 2026-09-22),
            # and is a rule about the mesh's numbering rather than about where the liquid is.
            idx = np.flatnonzero(lost)
            targets = self.patch_of_face[self.mesh.element_faces(old_surface.owner[idx])]   # (n, 4), -1: not a patch
            ok = targets >= 0
            safe = np.where(ok, targets, 0)
            # weight by each exposed face's area *projected on the vanished patch's own outward normal*, not by its raw
            # area: the face the recession uncovered underneath lies parallel to the patch that vanished (its new
            # outward normal points the same way), while the pit walls stand perpendicular to it. A prism element's
            # side faces can out-area its inward face -- 0.25 x 1.5 mm slivers against a 1 mm2 floor -- so raw areas
            # route most of the film sideways into walls no shear is pushing it towards, and concentrate it by the
            # element's aspect ratio (3 to 8 here). The projection is the footprint each face offers the liquid.
            cosine = np.clip(np.einsum("nij,nj->ni", self.surface.normals[safe], old_surface.normals[idx]), 0.0, None)
            area = np.where(ok, self.surface.areas[safe] * cosine, 0.0)
            total = area.sum(axis=1)
            flat = ok & (total[:, None] <= 0.0)                  # no face faces the right way: fall back to raw areas
            if flat.any():
                area = np.where(flat, self.surface.areas[safe], area)
                total = area.sum(axis=1)
            has = total > 0.0
            if has.any():
                share = old_m_f[idx[has]][:, None] * area[has] / total[has][:, None]
                np.add.at(m_f, targets[has][ok[has]], share[ok[has]])
            if (~has).any():                                     # the element exposed nothing: its death opened a
                orphan = idx[~has]                               # hole right through, so fall back to the nearest patch
                np.add.at(m_f, self._patch_tree.query(old_surface.centroids[orphan])[1], old_m_f[orphan])
        self.m_f = m_f

    # -- reporting -----------------------------------------------------------------------------------------------
    def film_temperature(self):
        """The film's temperature per patch: the surface's own (the film is thermally thin, class docstring)."""
        return self.surface_temperature()

    def film_enthalpy(self, h_node=None):
        """Specific enthalpy of each patch's film [J/kg]: the mean over the patch's three nodes of the *liquid*
        enthalpy h(T_i) + L_f (1 - f_l) -- the film is liquid by construction, so it carries its latent heat wherever
        it sits. Two things this must not be: the enthalpy of the mean temperature (the mean of the enthalpies and
        the enthalpy of the mean differ by the whole latent heat across the ramp, which left a -0.4 residual in the
        coupled balance), and the mixture enthalpy h(T) (which makes melting free on the ramp)."""
        if h_node is None:
            h_node = self.material.enthalpy_liquid(self.solver.temperature())
        return h_node[self.surface.faces].mean(axis=1) if self.m_f.size else np.zeros(0)

    def film_energy(self):
        """Enthalpy held by the film [J]. Identical to the nodal sum the solver carries (`set_film_mass` gives each
        node sum_p m_p/3), so the film's sensible heat is inside 1^T M dT and the balance closes on it exactly."""
        return float((self.m_f * self.film_enthalpy()).sum()) if self.m_f.size else 0.0

    def film_frozen_fraction(self):
        """Fraction of the film sitting on patches below the feed ramp: mass the enthalpy calls solid but the model
        still treats as liquid, because the mesh cannot grow a crust outside itself and its own element may be gone
        (section 9.5). It is the honest measure of that limitation, so it is recorded every step."""
        total = float(self.m_f.sum())
        if total <= 0.0:
            return 0.0
        return float(self.m_f[self.material.feed_fraction(self.film_temperature()) <= 0.0].sum() / total)

    @staticmethod
    def _mean_conjugate_um(delta_m):
        """Mean conjugate melt-layer depth [um] over the patches that have one; nan where the gate denies them all."""
        ok = np.isfinite(delta_m)
        return float(delta_m[ok].mean() * 1e6) if ok.any() else float("nan")

    def _shear_reaches(self, delta_m):
        """Per element, may its melt leave for the film this step? An element that owns part of the wall is in the
        sheared layer by definition and always may. A buried element may only if its centroid lies within Girin's
        conjugate melt-layer depth of the nearest patch, because that is how far the gas shear penetrates the liquid;
        deeper melt is stagnant and stays where it is until the surface reaches it. Where the wall Knudsen gate denies
        the Girin closure there is no conjugate depth, and then only wall-owning elements may feed, which is the
        conservative reading of the same statement."""
        owner = self.owner_area > 0.0
        if delta_m is None or not self.surface.n_patches:
            return owner.astype(float)
        dist, near = self._patch_tree.query(self.mesh.points[self.mesh.tets].mean(axis=1))
        d = delta_m[near]
        return (owner | (np.isfinite(d) & (dist <= d))).astype(float)

    def region_extent(self, selected):
        """Lateral extent [m] of *each* patch's own contiguous selected region, as that region's equivalent diameter
        2 sqrt(A/pi); zero where the patch is not selected.

        This is the domain a surface wave has to fit inside, and the region -- not the mesh facet -- is the physical
        boundary of it: melt is continuous across a facet edge, so a facet width is a bookkeeping artefact and using it
        would make the limit mesh-dependent (decided 2026-09-24). Contiguity is taken on the patch adjacency graph, so
        melt separated by solid does not add up to one wide pool. Each patch gets its own region rather than the
        body's largest, since a patch in an isolated puddle must not be credited with the width of the main pool."""
        n = self.surface.n_patches
        sel = np.asarray(selected, dtype=bool)
        out = np.zeros(n)
        if not n or not sel.any():
            return out
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        _, i, j = self.surface.edges()
        keep = sel[i] & sel[j]
        i, j = i[keep], j[keep]
        g = coo_matrix((np.ones(i.size * 2), (np.concatenate([i, j]), np.concatenate([j, i]))), shape=(n, n))
        count, label = connected_components(g, directed=False)
        area = np.bincount(label[sel], self.surface.areas[sel], count)
        out[sel] = 2.0 * np.sqrt(area[label[sel]] / np.pi)
        return out

    def unstable_region_extent(self, unstable):
        """Lateral extent [m] of the largest contiguous flagged region: the maximum of `region_extent`, kept for the
        reported diagnostics. A Rayleigh-Taylor wave has to fit inside the region where the melt is deep enough to
        support it, not on the nose as a whole: comparing its wavelength against the nose diameter asks whether the
        wave would fit on a surface most of which carries no pool (corrected 2026-09-23)."""
        e = self.region_extent(unstable)
        return float(e.max()) if e.size else 0.0

    def molten_depth(self, Te=None, max_levels=8):
        """Depth of *contiguous* fully molten material under each patch [m], measured inward from the patch.

        Marches inward through the mesh from the patch's owner element, hopping to the face-adjacent element furthest
        along the inward normal, and stops at the first element that is not fully molten (below `Material.T_feed`) or
        at the mesh boundary. Only contiguous liquid counts: molten material separated from the patch by solid is
        blocked from it, cannot be part of the layer the gas shear sees, and must not be credited to it -- accumulating
        every nearby molten element instead would count melt that cannot reach the wall (decided 2026-09-22).

        This is the thickness that belongs in the thick/thin instability test against Girin's conjugate melt-layer
        depth delta_m: the film mass on a patch is the mobile inventory, not the depth of liquid beneath the wall.
        Measured on the 100 mm physics flight, the two differ by an order of magnitude and invert the branch choice."""
        mat, mesh = self.material, self.mesh
        if not self.surface.n_patches or not mat.melts:
            return np.zeros(self.surface.n_patches)
        if Te is None:
            Te = self.solver.temperature()[mesh.tets].mean(axis=1)
        molten = mesh.active & (Te >= mat.T_feed)
        centroid = mesh.points[mesh.tets].mean(axis=1)
        inward = -self.surface.normals
        cur = self.surface.owner.copy()
        alive = molten[cur]
        depth = np.zeros(self.surface.n_patches)
        rows = np.arange(len(cur))
        for _ in range(max_levels):
            if not alive.any():
                break
            reach = np.einsum("ij,ij->i", centroid[cur] - self.surface.centroids, inward)      # projected depth so far
            depth = np.where(alive, np.maximum(depth, reach), depth)
            pair = mesh._face_elements[mesh.element_faces(cur)]                               # (n, 4, 2)
            nb = np.where(pair[:, :, 0] == cur[:, None], pair[:, :, 1], pair[:, :, 0])
            ok = nb >= 0
            step = np.einsum("nkj,nj->nk", centroid[np.where(ok, nb, 0)] - centroid[cur][:, None, :], inward)
            pick = np.where(ok, step, -np.inf).argmax(axis=1)
            nxt = nb[rows, pick]
            good = (nxt >= 0) & (step[rows, pick] > 0.0)
            alive = alive & good & molten[np.where(good, nxt, 0)]
            cur = np.where(good, nxt, cur)
        # the chain's last centroid sits half an element short of the far wall of that element: add that half
        return np.maximum(depth, 0.0) * 1.5

    def liquid_layer_depth(self, Te=None):
        """Depth of liquid at each patch [m]: its film plus the contiguous molten material beneath it. What the gas
        shear sees, and so what decides whether the film is thick or thin against the conjugate melt-layer depth."""
        if not self.m_f.size:
            return np.zeros(0)
        return self.m_f / (self.liquid.rho * self.surface.areas) + self.molten_depth(Te)

    def film_thickness_max(self):
        """Thickest film [m] over the patches where a thickness means anything: at least a tenth of the median patch
        area (slivers excluded) and no deeper than the patch is wide (b <= sqrt(A)). A film deeper than its patch is
        wide is not a film -- m_f/(rho_l A) is then a volume per area, not a depth, and the lubrication picture the
        number belongs to has already failed. Those patches are reported separately by `film_blob_fraction` rather
        than quietly averaged in: without the second test this diagnostic read 2431 mm on a 50 mm sphere in the last
        step of its flight, where thousands of dying patches handed their film to a few 0.6 mm2 survivors (measured
        2026-09-22; the mass was real, the depth was not)."""
        if not self.m_f.size:
            return 0.0
        a = self.surface.areas
        b = self.m_f / (self.liquid.rho * a)
        ok = (a >= 0.1 * np.median(a)) & (b <= np.sqrt(a))
        return float(b[ok].max()) if ok.any() else 0.0

    def film_blob_fraction(self):
        """Fraction of the film sitting deeper than its patch is wide (b > sqrt(A)) -- melt the lubrication film
        model cannot describe, because it no longer fits on the facet it is booked to. It is a geometry limit of the
        patch graph, not a mass error: the mass and the energy are still accounted exactly. Measured on the 50 mm
        physics flight it is zero until the last 5 % of the flight, where the body collapses from 3026 patches to 56."""
        total = float(self.m_f.sum())
        if total <= 0.0:
            return 0.0
        a = self.surface.areas
        return float(self.m_f[self.m_f / (self.liquid.rho * a) > np.sqrt(a)].sum() / total)

    def film_thickness_mean(self):
        """Film mass over the area it meaningfully wets [m]: patches at least B_MIN deep, the same floor below which
        spraying ignores a film. Counting every patch holding a nonzero number instead makes this diagnostic
        hypersensitive -- a patch with 1e-20 kg on it adds its whole area to the denominator, so a last-bit difference
        in the temperature field that nudges one element across the feed ramp moves the reported mean by 1-2 % while
        the film mass is bit-identical. That is how a reporting artefact was mistaken for a physics change
        (2026-09-23); the mass-based diagnostics (`film_blob_fraction`, `film_frozen_fraction`) never had the problem,
        because a patch holding nothing contributes nothing to a mass fraction."""
        from .spray import B_MIN
        if not self.m_f.size:                                    # no film account at all (e.g. --removal instant)
            return 0.0
        wet = self.m_f > self.liquid.rho * self.surface.areas * B_MIN
        return float(self.m_f[wet].sum() / (self.liquid.rho * self.surface.areas[wet].sum())) if wet.any() else 0.0

    def demised(self):
        """True once the body's material (film excluded) is below the demise fraction of the initial mass."""
        return self.consumed or float((self.phi * self.element_mass).sum()) < self.settings.demise_fraction * self.mass0

    def energy_balance_residual(self):
        """(E_body + E_film - E0 - absorbed heat + removed enthalpy + pending load + dropped load) / absorbed heat:
        zero for an exact discrete balance (the deferred loads cancel between the FEM and the removed material)."""
        if not self.absorbed_heat:
            return 0.0
        return (self.energy() - self.energy0 - self.absorbed_heat + self.removed_enthalpy + self.pending_load.sum()
                + self.total_dropped_load) / self.absorbed_heat

    @property
    def total_applied_load(self):
        return getattr(self, "_applied_total", 0.0)

    @property
    def total_dropped_load(self):
        return getattr(self, "_dropped_total", 0.0)

    def melt_stats(self):
        lm = self.last_melt
        fractions = lm.get("regime_fractions", [0.0, 0.0, 0.0])
        return {"mass_kg": self.mass(0.0), "film_mass_kg": float(self.m_f.sum()), "sprayed_mass_kg": self.sprayed_mass,
                "runoff_mass_kg": self.runoff_mass, "removed_mass_kg": self.removed_mass,
                "melt_front_depth_max_mm": self.melt_front_depth() * 1e3, "equivalent_radius_mm": self.equivalent_radius() * 1e3,
                "n_active_elements": float(self.mesh.n_active), "spraying_area_m2": lm.get("spraying_area_m2", 0.0),
                "theta_cr_deg": lm.get("theta_cr_deg", float("nan")), "n_released": self.n_released,
                "released_mass_kg": lm.get("released_mass", 0.0), "r_median_um": lm.get("r_median_um", float("nan")),
                "r_max_um": lm.get("r_max_um", float("nan")), "closure_fraction_girin": fractions[0],
                "closure_fraction_couette_slip": fractions[1], "closure_fraction_couette_fm": fractions[2], "rt_active": lm.get("rt_active", 0.0),
                "kn_body": lm.get("kn_body", float("nan")), "kn_local_stag": lm.get("kn_local_stag", float("nan")),
                "re_shock": lm.get("re_shock", float("nan")), "flow_branch": lm.get("flow_branch", float("nan")),
                "p_w_stag_Pa": lm.get("p_w_stag", float("nan")), "phi_sonic_deg": lm.get("phi_sonic_deg", float("nan")),
                "drag_shape_factor": self.drag_shape_factor(), "frozen_mass_kg": self.frozen_mass,
                "film_T_max_K": float(self.film_temperature().max()) if self.m_f.size else float("nan"),
                "film_T_mean_K": float((self.m_f * self.film_temperature()).sum() / self.m_f.sum()) if self.m_f.sum() > 0.0 else float("nan"),
                "film_frozen_fraction": self.film_frozen_fraction(), "unapplied_load_J": float(self.pending_load.sum()),
                "film_blob_fraction": self.film_blob_fraction(),
                "rt_mass_fraction": lm.get("rt_mass_fraction", float("nan")),
                "rt_growth_ms": lm.get("rt_growth_ms", float("nan")),
                "rt_region_mm": lm.get("rt_region_mm", float("nan")),
                "rt_bounded_fraction": lm.get("rt_bounded_fraction", float("nan")),
                "rt_bounded_growth_ms": lm.get("rt_bounded_growth_ms", float("nan")),
                "spray_growth_ms": lm.get("spray_growth_ms", float("nan")),
                "rt_wavelength_over_nose": lm.get("rt_wavelength_over_nose", float("nan")),
                "molten_depth_max_mm": lm.get("molten_depth_max_mm", float("nan")),
                "molten_depth_mean_mm": lm.get("molten_depth_mean_mm", float("nan")),
                "delta_m_mean_um": lm.get("delta_m_mean_um", float("nan")),
                "thick_branch_fraction": lm.get("thick_branch_fraction", float("nan")),
                "removed_enthalpy_J": self.removed_enthalpy,
                "film_thickness_max_mm": self.film_thickness_max() * 1e3, "film_thickness_mean_mm": self.film_thickness_mean() * 1e3,
                "nose_radius_mm": self.nose_radius() * 1e3, "transverse_radius_mm": self.transverse_radius * 1e3,
                "fitted_nose_radius_mm": self.fitted_nose_radius * 1e3,
                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active)}
```


- [ ] **Step 4: Extract the flat-disc endpoint of the drag family**

```bash
"$PY" - <<'EOF'
import scipy.io as sio, numpy as np, json
c = sio.netcdf_file("/Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/ATDB_CYLINDER.nc", "r", mmap=False)
ld = np.asarray(c.variables["log_LengthToDiameter"][:])
j = int(np.argmin(10.0 ** ld))                       # thinnest disc
ma = np.asarray(c.variables["MachNumber"][:])
get = lambda k: [round(float(v), 6) for v in np.asarray(c.variables[k][:])[0, :, j]]
doc = {
    "_provenance": ("Values of /Applications/DRAMA-4.1.4/TOOLS/SARA/REENTRY/data/ATDB_CYLINDER.nc (ESA DRAMA 4.1.4 aerothermal "
                    "database for the primitive CYLINDER, Hyperschall Technologie Goettingen GmbH) at angle of attack 0 (the flat "
                    "face normal to the flow) and the thinnest tabulated aspect ratio L/D = %.1e, read on 2026-09-22. This is the "
                    "flat-disc limit of the sphere-to-disc shape family: the face-on continuum C_D is independent of L/D over six "
                    "decades (1.824 at Ma 10 for every thickness), because hypersonic drag on a blunt body is pressure drag on the "
                    "frontal area while the side wall is parallel to the flow and the base sits in a near-vacuum wake. HTG's "
                    "continuum database is modified Newtonian: sphere/disc = 0.4990-0.4991 at every Mach number, i.e. exactly the "
                    "Newtonian ratio 1/2, so cd_continuum here IS C_p,max(Ma)." % (10.0 ** ld[j])),
    "mach": [float(v) for v in ma], "cd_free_molecular": get("FMF_CD"), "cd_continuum": get("CON_CD"),
    "heat_flux_factor_free_molecular": get("FMF_AvHeatFlux"), "heat_flux_factor_continuum": get("CON_AvHeatFlux"),
}
json.dump(doc, open("reentry_model/data/atdb_disc.json", "w"), indent=2)
print(doc["cd_continuum"])
EOF
```

Expected: `[1.801341, 1.824148, 1.828404, 1.829896, 1.830587, 1.830962]`, and the sphere table divided by it is 0.49897 at every Mach number.

In `reentry_model/aero.py` add, next to `DEFAULT_ATDB`:

```python
DISC_ATDB = os.path.join(DATA_DIR, "atdb_disc.json")     # the flat-disc limit (ATDB_CYLINDER at zero angle of attack)
```

and give `drag_coefficient` the shape factor:

```python
def drag_coefficient(kn, ma, tables, bridging, shape_factor=1.0):
    """SESAM's sphere C_D: continuum and free-molecular table values blended by f(Kn).

    `shape_factor` (Step 3) scales the continuum entry for a body that is no longer a sphere: it is the
    modified-Newtonian drag of the current windward shape divided by that of a sphere, so 1 for a sphere and 2.00 for
    a flat face (the ATDB disc entry is exactly 2 x the sphere entry at every Mach number -- HTG's continuum database
    is itself modified Newtonian). The free-molecular entry is left alone: face-on, a disc and a sphere differ by only
    3-4 % there (2.24 vs 2.15 at Ma 10), and a melting body is deep in the continuum by the time it flattens."""
    cd_c = tables.cd_continuum(ma) * shape_factor
    cd_fm = tables.cd_free_molecular(ma)
    f = bridging(kn) if kn > 0.0 else 0.0
    return cd_c + (cd_fm - cd_c) * f
```

- [ ] **Step 5: Edit `reentry_model/trajectory.py`**

In `Simulator.aero_state`, replace the three lines

```python
            kn = aero.knudsen(fs.rho, fs.m_bar, self.settings.diameter)
            cd = aero.drag_coefficient(kn, ma, self.tables, self.bridging)
            a_drag = -0.5 * fs.rho * V * v_rel * cd * self.area / self.body.mass(t)
```

with

```python
            kn = aero.knudsen(fs.rho, fs.m_bar, self.body.reference_length() or self.settings.diameter)   # a melting body's current size (Step 3)
            cd = aero.drag_coefficient(kn, ma, self.tables, self.bridging, self.body.drag_shape_factor())
            area = self.body.reference_area() or self.area          # a melting body's projected area (Step 3), else pi D^2/4
            m = self.body.mass(t)
            a_drag = -0.5 * fs.rho * V * v_rel * cd * area / m if m > 0.0 else np.zeros(3)      # a consumed body (Step 3) has no drag
```

- [ ] **Step 6: Smooth SESAM's Mach-1 drag step in `reentry_model/aero.py`**

A light melting remnant hovers at its terminal velocity near Ma 1, where the factor-2 step in `cd_continuum` stalls the adaptive integrator (1.3e5 RHS evaluations in one macro step, measured). After `DEFAULT_ATDB = os.path.join(DATA_DIR, "atdb_sphere.json")` add

```python
MACH_SWITCH_LO, MACH_SWITCH_HI = 0.98, 1.02      # SESAM halves C_D below Ma 1; the step is smoothed over this band (below)
```

replace `SphereDragTables.cd_continuum` with

```python
    def cd_continuum(self, ma):
        if ma < MACH_SWITCH_LO:
            return 0.5 * float(self.cd_c[0])
        if ma < MACH_SWITCH_HI:                                  # SESAM's factor-2 step at Ma 1, smoothed over +-2 % (module docstring)
            s = (ma - MACH_SWITCH_LO) / (MACH_SWITCH_HI - MACH_SWITCH_LO)
            return (0.5 + 0.5 * s * s * (3.0 - 2.0 * s)) * float(self.cd_c[0])
        if ma < self.mach[0]:
            return float(self.cd_c[0])
        return float(np.interp(ma, self.mach, self.cd_c))
```

and extend the class docstring's last sentence to: `simply clamped (Kn is negligible wherever Ma < 5). The factor-2 step is applied as a smooth (cubic) ramp over Ma 0.98-1.02: a discontinuous C_D stalls the adaptive integrator when a light body hovers at its transonic terminal velocity (a melting remnant, Step 3: 1.3e5 RHS evaluations in one macro step, measured 2026-09-21); the reference spheres cross Ma 1 in a fraction of a second, where the ramp changes nothing measurable (the drag test excludes |Ma - 1| <= 0.02 rows for SESAM's 3-decimal Mach column)."""`. In `tests/test_reentry_model_aero.py` replace the assertion `assert t.cd_continuum(3.0) == 0.898818 and t.cd_continuum(1.0) == 0.898818` with

```python
        assert t.cd_continuum(3.0) == 0.898818 and t.cd_continuum(1.02) == 0.898818
        assert t.cd_continuum(1.0) == pytest.approx(0.75 * 0.898818) and t.cd_continuum(0.98) == 0.5 * 0.898818   # the Ma-1 step smoothed over +-2 % (Step 3)
```

- [ ] **Step 7: Run the tests**

```bash
"$PY" -m pytest tests/test_reentry_model_melting.py tests/test_reentry_model_coupled.py tests/test_reentry_model_trajectory.py tests/test_reentry_model_aero.py -q
"$PY" -m pytest -m reference tests/test_reentry_model_reference.py -q                        # the Step 1 flights, ~2 min: unchanged by the ramp
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
```

Expected: 6 passed (melting), the Step 2 coupled/trajectory/aero tests unchanged, the eight Step 1 reference flights within their thresholds; in `fenicsx_env` 8 passed — the two backends give the same melting run (mass to 1e-6, temperatures to 0.5 K, the same dead elements).

- [ ] **Step 8: Commit**

```bash
git add reentry_model/body.py reentry_model/trajectory.py reentry_model/aero.py reentry_model/data/atdb_disc.json tests/test_reentry_model_melting.py tests/test_reentry_model_aero.py tests/test_reentry_model_data.py
git commit -m "Add the melting body: feed, film, spraying, element death, accounting and the projected area (Step 3 Task 9)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

