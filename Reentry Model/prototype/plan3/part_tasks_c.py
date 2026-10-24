# ---------------------------------------------------------------------------------------------------------- Task 9
emit(r'''### Task 9: The melting body

**Files:**
- Modify (replace): `reentry_model/body.py`
- Modify: `reentry_model/trajectory.py` (three statements), `reentry_model/aero.py` (the Mach-1 step and the shape factor)
- Create: `reentry_model/data/atdb_disc.json` (extracted from DRAMA)
- Test: `tests/test_reentry_model_melting.py` (new); `tests/test_reentry_model_aero.py` (one assertion); `tests/test_reentry_model_fenicsx.py::test_melting_run_matches_the_skfem_backend` (from Task 3) now runs in `fenicsx_env`

**Interfaces:**
- Consumes: Tasks 1–7 (`mesh.deactivate/surface/active_nodes`, `Material.feed_fraction/liquid_fraction/enthalpy/h_liquid/liquid`, `thermal` `set_fractions/element_energies/step(nodal_load=)/pinned/temperature`, `surface_flow.SurfaceFlow`, `spray.SprayModel/melt_layer/source_rows/histogram/N_BINS`, `film.lubrication/Runoff`), `heating.HeatingResult`, `trajectory.AeroState`.
- Produces: `body.REMOVAL_NAMES = ("girin", "instant")`, `SIZE_FEEDBACK_NAMES = ("current", "initial")`, `NOSE_CAP_ANGLE = 30.0`, `NOSE_CAP_FACTOR = 1.67`, `CP_MAX_NEWTONIAN = 1.84`, `PHI_MIN = 1e-3`, `PHI_DEATH = 0.05`, `NEAREST_PATCHES = 4`, `LOAD_DT_MAX = 1000.0`, `FEED_DEPTH_NAMES = ("conjugate", "all")`; `fit_sphere(points) -> (centre, radius)`; `MeltSettings(removal, runoff, demise_fraction, particles, size_feedback="current", feed_depth="conjugate")`; `MeltingBody(mesh, material, solver, mass_kg, flow=None, spray_model=None, settings=None, T0, emissivity, T_ambient, v_hat)` with `.advance(t, dt, loads, state=None)`, `.melt_step(t, dt, state)`, `.mass(t)`, `.reference_area()`, `.reference_length()` (2 R_eq, or None with `initial`), `.nose_radius()` (the windward-cap fit, or R₀ with `initial`), `.newtonian_drag() -> (C_D, projected area)` (modified Newtonian over the windward convex hull), `.drag_shape_factor()` (that value over the meshed sphere's own, so exactly 1 while intact and 2.00 at the flat limit; 1 with `initial`), `.equivalent_radius()`, `.energy()` (FEM + film), `.film_temperature()` (the surface's own: the film is thermally thin), `.film_enthalpy(h_node=None)` (per patch, the mean of the *liquid* nodal enthalpies), `.film_energy()`, `.film_frozen_fraction()`, `.film_blob_fraction()` (the share of the film deeper than its patch is wide), `.on_current_surface(values, fill)` (a per-patch array of the last flow and spray evaluation carried onto the current surface by face id, `fill` on faces the step's element deaths exposed; fact 37),`.molten_depth(Te=None, max_levels=8)` (the contiguous liquid depth inward from each patch, hopping element to element and stopping at the first one that is not molten), `.liquid_layer_depth()` (that depth plus the film, which is what the branch test and the Rayleigh-Taylor criterion see), `.region_extent(selected)` (each patch's own contiguous region as an equivalent diameter 2 sqrt(A/pi), which is the domain a wave has to fit inside -- the molten region, not the facet) and `.unstable_region_extent(unstable)` (its maximum, for the reported diagnostics; facts 30, 32 and 33), `.total_applied_load()`, `.total_dropped_load()` (the deferred-load accounts the energy balance closes against), `.mean_temperature()`, `.melt_front_depth()`, `.film_thickness_max/mean()`, `.demised()`, `.energy_balance_residual()`, `.melt_stats() -> dict` (the `coupled.MELT_COLUMNS` values plus `mass_kg`), attributes `phi, m_f, surface, theta, t_hat, windward, mass0, mass_centre, transverse_radius, fitted_nose_radius, cap_nose_radius, sprayed_mass, runoff_mass, removed_mass, removed_enthalpy, n_released, source_rows, hist_n, hist_m, melt_onset, spray_onset, consumed, last_flow, last_spray, last_b, last_face_ids, last_melt, pending_load, liquid, flow, spray, runoff`; `Body.reference_area()`/`reference_length()`/`drag_shape_factor()` in the protocol (`ConstantBody`/`ThermalBody` return None; `ThermalBody.nose_radius()` returns the sphere's radius; `ThermalBody.advance` accepts `state=None`); `Simulator.aero_state` uses the body's reference area, Knudsen length and drag shape factor, and gives a consumed body no drag; `aero.drag_coefficient(kn, ma, tables, bridging, shape_factor=1.0)` scales the continuum entry only.

The film's temperature, its re-solidification and the way every transfer is booked are the subject of measured facts 25 and 26: the film has no energy equation of its own (its *mass* goes to the solver, which carries it on the boundary nodes with the liquid capacity), it holds `enthalpy_liquid` wherever it sits, the feed and the freeze are netted into one transfer per element, each transfer books only the enthalpy difference it carries and books it at the destination, and no deferred load may move a node more than `LOAD_DT_MAX` in one step. Read those two facts before reading `melt_step`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_melting.py`:

''')
emit(file_block("tests/test_reentry_model_melting.py"))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_melting.py -q`
Expected: AttributeError (`body.MeltSettings`).

- [ ] **Step 3: Replace `reentry_model/body.py`**

''')
emit(file_block("reentry_model/body.py"))
emit(r'''
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
''')

# ---------------------------------------------------------------------------------------------------------- Task 10
emit(r'''### Task 10: The coupled loop — melt columns, demise, VTK melt fields, particle files

**Files:**
- Modify (replace): `reentry_model/coupled.py`
- Test: `tests/test_reentry_model_coupled.py` (append)

**Interfaces:**
- Consumes: `MeltingBody` (Task 9: `advance(state=)`, `demised()`, `melt_stats()`, `source_rows`, `m_f`, `liquid`, `last_spray/last_flow`, `on_current_surface(values, fill)`, `phi`, `material.liquid_fraction`), `spray.SOURCE_COLUMNS/BIN_EDGES/N_BINS/histogram` (Task 7).
- Produces: `coupled.MELT_COLUMNS` (23 names), `CoupledRun.melting` (property), `loads_at` passing `body.nose_radius()` to the heating, `melt_results(history)` (melt/spraying onsets, demise, masses, `n_released`, `r_median_um`, `n_dead_elements`, `n_source_rows`, `removed_enthalpy_J`, `melt_energy_balance_residual`), end reason `"demise"`; `write_vtk_frame` writes active cells only with `liquid_fraction`/`phi` and the surface fields `film_thickness, film_T, we_s, closure, kn_local, p_w, tau, r_droplet, release_rate` when melting, the flow and spray ones carried across the step's element deaths by face id (`MeltingBody.on_current_surface`, fact 37); `write_particles(run_dir, body, history, window=10.0) -> {"particles", "particles_summary", "size_distribution"}`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_coupled.py`:

''')
emit(block(tail_from("tests/test_reentry_model_coupled.py", "# ---------------------------------------------------------------------------------------------------------------\n# Step 3: the melting run")))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_coupled.py -q`
Expected: the two new tests fail (`advance()` gets no `state`, no `MELT_COLUMNS`, no `write_particles`).

- [ ] **Step 3: Replace `reentry_model/coupled.py`**

''')
emit(file_block("reentry_model/coupled.py"))
emit(r'''
- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_coupled.py tests/test_reentry_model_viz.py tests/test_reentry_model_cli.py -q`
Expected: all pass (Step 2's runs are unchanged: `melting` is False for a `ThermalBody`).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/coupled.py tests/test_reentry_model_coupled.py
git commit -m "Carry the melt step, the demise end and the particle writers through the coupled loop (Step 3 Task 10)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

