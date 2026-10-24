# ---------------------------------------------------------------------------------------------------------- Task 1
emit(r'''### Task 1: Prism-layered mesh, active set, box mesh

**Files:**
- Modify (replace): `reentry_model/mesh.py`
- Test: `tests/test_reentry_model_mesh.py` (append)

**Interfaces:**
- Consumes: gmsh generation and `VolumeMesh`/`SurfaceMesh` of Step 2 (kept: `generate_sphere_mesh`, `load_mesh`, `boundary_faces`, `mesh_file_name`, `centre_node`, `boundary_nodes`, `angles_to`, `facet_mean`, `patch_toward`).
- Produces: `sphere_mesh(radius, h_surface, h_core, mesh_dir, layers=0, layer_thickness=0.25e-3, growth=2.0)`; `layer_thicknesses(layers, layer_thickness, growth)`; `split_prisms(bottom, top)`; `add_prism_layers(inner, radius, layers, layer_thickness, growth)`; `box_mesh(lx, ly, lz, h)`; `VolumeMesh.active` (bool per element), `.element_layer` (−1 core, 0 outermost layer), `.n_active`, `.deactivate(elements) -> (gone_face_ids, new_face_ids)`, `.active_nodes()`, `.element_faces(elements) -> (n, 4) face ids` (the stable ids of an element's own four faces, which the death hand-over of Task 10 uses to find the patches a vanishing element exposes), `.volume()` (active), `.surface()` over the active set with `SurfaceMesh.owner` (element per patch) and `.face_ids` (stable ids into the face table `_face_nodes`); `SurfaceMesh.tangent_from(v_hat)`, `.projected_area(v_hat)`, `.edges() -> (edge node pairs, patch i, patch j)`; constants `DEFAULT_LAYERS = 4`, `DEFAULT_LAYER_THICKNESS = 0.25e-3`, `DEFAULT_LAYER_GROWTH = 2.0`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_mesh.py`:

''')
emit(block(tail_from("tests/test_reentry_model_mesh.py", "# ---------------------------------------------------------------------------------------------------------------\n# Step 3: prism layers")))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: the five new tests fail (`layer_thicknesses`, `split_prisms`, `box_mesh` missing; `sphere_mesh` has no `layers`).

- [ ] **Step 3: Replace `reentry_model/mesh.py`**

''')
emit(file_block("reentry_model/mesh.py"))
emit(r'''
- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_mesh.py -q`
Expected: 10 passed (the five Step 2 tests unchanged; the layered fixture builds a 4 mm / 20 mm inner sphere with two layers, ~5 k nodes). Then the unit tier: `"$PY" -m pytest -m "not drama and not reference" -q` — everything that used `mesh.surface()` still passes (the boundary now comes from the face table; the coarse-mesh fixtures are unchanged).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/mesh.py tests/test_reentry_model_mesh.py
git commit -m "Add prism layers, the active element set and the box mesh to the sphere mesh (Step 3 Task 1)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

# ---------------------------------------------------------------------------------------------------------- Task 2
emit(r'''### Task 2: Material with latent heat, melting ranges and liquid properties

**Files:**
- Modify (replace): `reentry_model/material.py`
- Create: `reentry_model/data/materials/AA7075.json`, `reentry_model/data/materials/AA7075_range.json`
- Test: `tests/test_reentry_model_material.py` (append)

**Interfaces:**
- Consumes: `reentry_model/data/materials/AA7075_nomelt.json` (Step 2).
- Produces: `Material` fields `latent_heat`, `T_solidus`, `T_liquidus`, `liquid: LiquidProperties(rho, mu, sigma)`; properties `melts`, `T_feed`, `h_liquid`; methods `liquid_fraction(T)`, `feed_fraction(T)`, `cp_eff(T)`, `enthalpy(T)` (exact, latent slope inside the range), `temperature_from_enthalpy(h)`, and the melt film's four: `enthalpy_liquid(T)` = h(T) + L_f (1 - f_l(T)) (what a kilogram of *liquid* holds at T -- the film carries its latent heat wherever it sits), `enthalpy_mixed(T, w)` and `cp_mixed(T, w)` for a node holding a fraction `w` of film and 1 - w of material, and `temperature_from_enthalpy_mixed(h, w)`, their exact inverse in T for every w (round-trip 4e-12 K, measured); `Material.from_drama_json(path=None)` accepting the names in `MATERIAL_NAMES` (`AA7075_nomelt`, `AA7075`, `AA7075_range`); constants `MELT_RAMP = 2.0`, `NO_MELT_ABOVE = 5000.0`. Single-temperature materials get a ±MELT_RAMP ramp; `feed_fraction` is a ±MELT_RAMP ramp ending at `T_feed` (the liquidus, +2 K for range materials).

- [ ] **Step 1: Write the two material files**

Generate them from the packaged no-melt file (the DRAMA tables verbatim, the 1e5 K holding row dropped, the liquid properties added):

```bash
"$PY" - <<'EOF'
import json
d = json.load(open("reentry_model/data/materials/AA7075_nomelt.json"))
base = {k: v for k, v in d.items() if k != "_comment"}
base["specificHeatCapacity"] = [r for r in d["specificHeatCapacity"] if r[0] <= 850.0]
base["heatConductivity"] = [r for r in d["heatConductivity"] if r[0] <= 850.0]
liquid = {"density": 2400.0, "viscosity": 1.3e-3, "surfaceTension": 0.86,
          "_sources": "pure aluminium near the liquidus: rho_l 2375-2400 kg/m3 (Smithells Metals Reference Book, 8th ed., Table 14.1), mu_l 1.2-1.4 mPa s at 933-1000 K (Assael et al. 2006, J. Phys. Chem. Ref. Data 35, 285), sigma 0.86-0.91 N/m (Smithells; ASM Handbook Vol. 2). Alloy corrections for AA7075 (5.6 % Zn, 2.5 % Mg, 1.6 % Cu) are within 10 %; used as assumptions (spec 2026-09-20 sections 2 and 16)."}
a = dict(base); a["name"] = "AA7075"; a["meltingTemperature"] = 850.0; a["liquid"] = liquid
a["_comment"] = "DRAMA 4.1.4 drama-AA7075 verbatim (TOOLS/material_database.xml: density, cp(T) and k(T) to 850 K, meltingHeat 400 kJ/kg, meltingTemperature 850 K, emissivity 0.4, catalycity 1, oxidation off) plus the liquid-phase properties DRAMA does not carry (`liquid`). Above 850 K the tables hold their last value (cp 1131.6 J/kgK, k 128.19 W/mK). Step 3 material (spec 2026-09-20 section 6): single melting temperature, numerical ramp of +-2 K in the model. DRAMA's database itself is untouched."
r = dict(a); r["name"] = "AA7075_range"; r["solidusTemperature"] = 750.0; r["liquidusTemperature"] = 908.0; r["meltingTemperature"] = 908.0
r["_comment"] = "AA7075 (see AA7075.json) with the alloy's melting range instead of DRAMA's single temperature: solidus 750 K, liquidus 908 K (ASM Handbook Vol. 2, Properties and Selection: Nonferrous Alloys, AA7075 477-635 C), latent heat 400 kJ/kg spread linearly across the range. meltingTemperature is set to the liquidus. Default material of --melt on (spec 2026-09-20 sections 6 and 14); melt and runoff start at the liquidus (section 8)."
for name, doc in (("AA7075", a), ("AA7075_range", r)):
    keys = ["_comment", "name", "materialType", "catalycity", "density", "meltingHeat", "meltingTemperature"] + (["solidusTemperature", "liquidusTemperature"] if "solidusTemperature" in doc else []) + ["specificHeatCapacity", "heatConductivity", "emissivity", "oxideActivationTemperature", "oxideEmissivity", "oxideHeatOfFormation", "oxideReactionProbability", "liquid"]
    json.dump({k: doc[k] for k in keys}, open("reentry_model/data/materials/%s.json" % name, "w"), indent=2)
print("written")
EOF
```

Expected: both files exist; `AA7075.json` has 29 c_p rows ending at 850 K, `meltingHeat` 400000, `meltingTemperature` 850; `AA7075_range.json` adds `solidusTemperature` 750 and `liquidusTemperature` 908 (`meltingTemperature` 908).

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_reentry_model_material.py`:

''')
emit(block(tail_from("tests/test_reentry_model_material.py", "# ---------------------------------------------------------------------------------------------------------------\n# Step 3: latent heat")))
emit(r'''
- [ ] **Step 3: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_material.py -q`
Expected: the five new tests fail (`from_drama_json("AA7075")` is a bad path; no `melts`, `feed_fraction`, `liquid`, `enthalpy_liquid`).

- [ ] **Step 4: Replace `reentry_model/material.py`**

''')
emit(file_block("reentry_model/material.py"))
emit(r'''
- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_material.py tests/test_reentry_model_thermal.py tests/test_reentry_model_cli.py -q`
Expected: all pass (the no-melt material behaves exactly as in Step 2: `latent_heat` 0, `melts` False, enthalpy unchanged).

- [ ] **Step 6: Commit**

```bash
git add reentry_model/material.py reentry_model/data/materials/AA7075.json reentry_model/data/materials/AA7075_range.json tests/test_reentry_model_material.py
git commit -m "Add the latent heat, melting ranges, feed fraction and liquid properties to the material (Step 3 Task 2)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

# ---------------------------------------------------------------------------------------------------------- Task 3
emit(r'''### Task 3: Thermal core — nodal enthalpy, element fractions, pinned nodes, nodal loads

**Files:**
- Modify (replace): `reentry_model/thermal/__init__.py`, `reentry_model/thermal/skfem_backend.py`, `reentry_model/thermal/fenicsx_backend.py`
- Modify: `reentry_model/cli.py` (the `--lumped-mass` flag becomes `--consistent-mass`, three lines)
- Test: `tests/test_reentry_model_thermal.py` (one test changed, five appended), `tests/test_reentry_model_fenicsx.py` (one line changed, two tests appended)

**Interfaces:**
- Consumes: `VolumeMesh.active`, `.surface()` (Task 1); `Material.enthalpy/enthalpy_liquid/enthalpy_mixed/cp/cp_eff/cp_mixed/temperature_from_enthalpy/temperature_from_enthalpy_mixed/melts/T_solidus/T_liquidus` (Task 2).
- Produces: `StepResult(T, Q_conv, Q_rad, iterations, Q_extra=0.0, Q_dropped=0.0)`; solver methods `set_fractions(phi)` (0 = dead; refreshes the boundary and the pinned nodes), `element_energies(T=None)` (φ_e ρ V_e mean_i h(T_i)), `step(dt, q_conv, T_amb, dirichlet=None, nodal_load=None)` (nodal_load in W, mesh node order); attributes `phi`, `pinned`, `faces`, `areas`, `last_damping`; `SkfemThermalSolver(lumped_mass=True)` default with `consistent_mass = not lumped_mass`; `element_matrices(points, tets) -> (vol, Ke, Mk)` with `Mk[e, k]` the nodal-coefficient mass matrices and `nodal_mass_weights()`; `mass_matrix(c_nodal, c_film=None)`; `operators(T, T_old=None) -> (K, M_tan[, E])`; and, for the melt film (Step 3, Task 9): `set_film_mass(mass)` (one non-negative value per node, kg -- it joins the nodes' capacity and keeps a node live even when its elements are gone), `nodal_capacity()` (J/K per node, material + film, which the body uses to bound the melt loads it defers) and `film_weight()` (the film's share of each node's mass, which the enthalpy-consistent Newton update inverts against). `facet_temperature(T=None)` (the mean of a facet's three nodal temperatures, defaulting to the solver's own stored field) is on both backends with the same signature -- they stand behind one protocol, so a call that works on either must work on both, and a test asserts it. The FEniCSx backend has the same public surface (`lumped_mass=False` raises `ValueError`).

The film's capacity is the **liquid** one -- tangent c_p(T) in M, the secant of `enthalpy_liquid` in E -- never the mixture's c_p,eff: the film has already paid its latent heat and must not pay it again on the ramp. The same weighting makes 1^T E the exact increment of (solid nodal enthalpy + film liquid enthalpy), which is what `MeltingBody.energy()` sums, so the coupled balance closes on it. CG gets one fresh AMG hierarchy and then a direct solve if it still stalls (melting drains interior elements to phi ~ 1e-3 with unscaled conduction; `direct_fallbacks` counts it).

- [ ] **Step 1: Update the conformance test and append the melting tests**

In `tests/test_reentry_model_thermal.py` replace the body of `test_operators_match_scikit_fem_assembly` — the lines

```python
def test_operators_match_scikit_fem_assembly(coarse_sphere_mesh):
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json())
    T = 300.0 + 400.0 * np.random.default_rng(1).random(coarse_sphere_mesh.n_nodes)
    K, M = s.operators(T)
    K_ref, M_ref = s.reference_operators(T)
    assert abs(K - K_ref).max() < 1e-10 * abs(K_ref).max() and abs(M - M_ref).max() < 1e-10 * abs(M_ref).max()
    assert abs(np.asarray(K.sum(axis=1))).max() < 1e-9 * abs(K).max()          # rows of K sum to zero
    assert M.sum() == pytest.approx((RHO * s.material.cp(T[s.tets].mean(axis=1)) * s.vol).sum(), rel=1e-12)
```

with

```python
def test_operators_match_scikit_fem_assembly(coarse_sphere_mesh):
    """K (element-mean k) against scikit-fem's assembly; the lumped capacity matrix is diag(sum_e V_e/4 rho c_p(T_i))
    (the nodal-enthalpy lumping, not the row sums of the consistent matrix), and the consistent option matches
    scikit-fem's consistent mass with the nodal (P1) c_p."""
    import scipy.sparse as sp
    s = solver(coarse_sphere_mesh, material.Material.from_drama_json())
    T = 300.0 + 400.0 * np.random.default_rng(1).random(coarse_sphere_mesh.n_nodes)
    K, M = s.operators(T)
    K_ref, M_ref = s.reference_operators(T)
    nodal = np.bincount(s.tets.ravel(), np.repeat(s.vol / 4.0, 4) * (RHO * s.material.cp(T))[s.tets].ravel(), minlength=s.points.shape[0])
    assert abs(K - K_ref).max() < 1e-10 * abs(K_ref).max() and abs(M - sp.diags(nodal)).max() < 1e-10 * abs(M_ref).max()
    _, M_c = solver(coarse_sphere_mesh, material.Material.from_drama_json(), lumped_mass=False).operators(T)
    assert abs(M_c - M_ref).max() < 1e-10 * abs(M_ref).max()
    assert abs(np.asarray(K.sum(axis=1))).max() < 1e-9 * abs(K).max()          # rows of K sum to zero
    assert M.sum() == pytest.approx((RHO * s.material.cp(T[s.tets]).mean(axis=1) * s.vol).sum(), rel=1e-12)   # nodal c_p, V_e/4 per node
```

(the remaining three lines of the test — the linear-field check on `K1` — stay). Then append to the file:

''')
emit(block(tail_from("tests/test_reentry_model_thermal.py", "# ---------------------------------------------------------------------------------------------------------------\n# Step 3: element fractions")))
emit(r'''
In `tests/test_reentry_model_fenicsx.py` change the constructor test's line `thermal.thermal_solver("fenicsx", lumped_mass=True)` to `thermal.thermal_solver("fenicsx", lumped_mass=False)           # the nodal-enthalpy (lumped) capacity only` and append:

''')
emit(block(tail_from("tests/test_reentry_model_fenicsx.py", "def test_melting_run_matches_the_skfem_backend")))
emit(r'''
(The appended FEniCSx test needs Tasks 5–9; it stays failing in `fenicsx_env` until Task 9 and is skipped in `drama_env`.)

- [ ] **Step 2: Run the thermal tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q`
Expected: the conformance test and the four new tests fail (no `set_fractions`, no `nodal_load`, `lumped_mass` default False).

- [ ] **Step 3: Replace `reentry_model/thermal/__init__.py`**

''')
emit(file_block("reentry_model/thermal/__init__.py"))
emit(r'''
- [ ] **Step 4: Replace `reentry_model/thermal/skfem_backend.py`**

''')
emit(file_block("reentry_model/thermal/skfem_backend.py"))
emit(r'''
- [ ] **Step 5: Replace `reentry_model/thermal/fenicsx_backend.py`**

''')
emit(file_block("reentry_model/thermal/fenicsx_backend.py"))
emit(r'''
- [ ] **Step 6: Rename the CLI flag**

In `reentry_model/cli.py`: replace `th.add_argument("--lumped-mass", action="store_true")` with

```python
    th.add_argument("--consistent-mass", action="store_true", help="consistent capacity matrix instead of the lumped nodal-enthalpy one (skfem only, analytic checks)")
```

replace `lumped_mass=args.lumped_mass` (in `build_thermal`) with `lumped_mass=not args.consistent_mass`, and `"lumped_mass": args.lumped_mass}` (the info dict) with `"lumped_mass": not args.consistent_mass}`.

- [ ] **Step 7: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py tests/test_reentry_model_coupled.py tests/test_reentry_model_cli.py tests/test_reentry_model_viz.py -q`
Expected: all pass — the Step 2 analytic cases (lumped limit, radiative cooling, Carslaw–Jaeger, balance) are unchanged by the lumping for constant c_p, and the melting cases converge in ≤ 8 iterations. Then in `fenicsx_env`:

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q --deselect tests/test_reentry_model_fenicsx.py::test_melting_run_matches_the_skfem_backend
```

Expected: 7 passed (the seven Step 2 conformance tests, including the skfem cross-check, with the nodal scheme).

- [ ] **Step 8: Commit**

```bash
git add reentry_model/thermal tests/test_reentry_model_thermal.py tests/test_reentry_model_fenicsx.py reentry_model/cli.py
git commit -m "Move the thermal core to the nodal enthalpy with element fractions, pinned nodes and nodal loads (Step 3 Task 3)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

# ---------------------------------------------------------------------------------------------------------- Task 4
emit(r'''### Task 4: Girin's dispersion relation and its cached table

**Files:**
- Create: `reentry_model/dispersion.py`, `reentry_model/data/girin_dispersion.json` (generated)
- Test: `tests/test_reentry_model_dispersion.py`

**Interfaces:**
- Produces: `cubic_coefficients(delta, we) -> (a2, a1, a0)`; `growth_rates(we, delta)`; `fastest_mode(we, delta=DELTA_GRID) -> (Delta_f, Im Omega_f, Re Omega_f)` (nan/0/0 when stable); `build_table()`, `write_table(path)`; `DispersionTable(path=TABLE_PATH)` callable on scalars/arrays returning `(Delta_f, Im Omega_f, Re Omega_f)` with attribute `we_onset`; constants `WE_CRITICAL_THEORY = 3.08`, `WE_CRITICAL_PRACTICAL = 4.62`, `TABLE_PATH`, and the grids the table is built on, `WE_GRID`, `DELTA_GRID` and `UNSTABLE_TOL`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_dispersion.py`:

''')
emit(file_block("tests/test_reentry_model_dispersion.py"))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_dispersion.py -q`
Expected: ImportError (`reentry_model.dispersion`).

- [ ] **Step 3: Create `reentry_model/dispersion.py`**

''')
emit(file_block("reentry_model/dispersion.py"))
emit(r'''
- [ ] **Step 4: Generate the packaged table and run the tests**

```bash
"$PY" -c "from reentry_model import dispersion; t = dispersion.write_table(); print(len(t['we']), t['delta_f'][-1], t['im_omega_f'][-1])"
"$PY" -m pytest tests/test_reentry_model_dispersion.py -q
```

Expected: `200 1.2256... 0.2469...` (2–3 s), then 3 passed.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/dispersion.py reentry_model/data/girin_dispersion.json tests/test_reentry_model_dispersion.py
git commit -m "Solve Girin's dispersion relation numerically and cache the fastest-mode table (Step 3 Task 4)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')
