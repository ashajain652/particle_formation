# Sub-plan: Task 10 — The coupled loop: melt columns, demise, VTK melt fields, particle files

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 5704–6029). Read `00-shared-context.md` first. **Amended 2026-10-02: the deep runoff and the per-patch conjugate depth** (section below). **Amended 2026-10-03: the molten cascade** (the section after it).


> **Amended 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md`
> (measured facts 38–45 in `00-shared-context.md`). The block below takes precedence over the extracted body.

**Changes to Task 10.**

1. **One new stage in the melt step**, after element death (step v) and before the history row: the recession solve
   and derived-surface update of Task 16, then the remesh trigger check of Task 17, then a remesh if triggered. The
   trajectory, heating and conduction stages are unchanged in order and content.
2. **A remesh rebuilds the world mid-run.** After Task 17 returns a new mesh, `CoupledRun` must re-setup the thermal
   solver, re-derive the surface helpers, and re-seat the film per patch, `φ_e` per element and `pending_load` per
   node from the transferred fields. This is the only place in the run where the solver's mesh changes, and both
   backends must survive it identically (fact 45; spec §10.1 item 6 makes the two-backend agreement across a remesh
   a non-negotiable acceptance test).
3. **New `MELT_COLUMNS`**: `surface_edge_max_mm`, `surface_edge_mean_mm`, `recession_residual_frac`,
   `phi_surface_min`, `conduction_dT_err_K`, `remesh_count`, `remesh_volume_residual_frac`,
   `remesh_energy_residual_frac`. Task 15's audit (`prototype/plan3/audit_docs.py`) checks the README's column list
   against `MELT_COLUMNS`, so both must change together.
5. **Also carried here: measured fact 37 (2026-09-27)**, missing from the committed master plan for the reason given
   in Task 9's amendment. `MELT_COLUMNS` is unchanged, but the VTK surface writer must read the flow and spray
   fields through `MeltingBody.on_current_surface`, which carries them across the step's element deaths by face id
   rather than by index, with the old defaults on faces the deaths exposed:

```python
        res, flow, carry = body.last_spray, body.last_flow, body.on_current_surface
        poly.cell_data["we_s"] = carry(res.we_s if res is not None else None, 0.0)
        poly.cell_data["closure"] = carry(flow.closure.astype(float) if flow is not None else None, 1.0)
        poly.cell_data["kn_local"] = carry(flow.kn_local if flow is not None else None, np.nan)
        poly.cell_data["p_w"] = carry(flow.p_w if flow is not None else None, 0.0)
        poly.cell_data["tau"] = carry(flow.tau if flow is not None else None, 0.0)
        poly.cell_data["r_droplet"] = carry(np.where(res.dm > 0.0, res.r, np.nan) if res is not None else None, np.nan)
        poly.cell_data["release_rate"] = carry(res.dm if res is not None else None, 0.0) / surface.areas
```

4. **The VTK series changes topology at a remesh.** `viz` reads only exported files (CLAUDE.md), so this is legal,
   but the writer must start a new series index rather than assume constant connectivity.
---

## Amendment of 2026-10-02 — the conjugate depth and the deep liquid in every surface frame, and five deep columns

> Part of the deep-runoff amendment (sub-plan 09's amendment of this date; facts 46–53 in `00-shared-context.md`).
> The code below is the tested change, as a diff against the prototype's `coupled.py` (the body below with the
> fact-37 lines of this sub-plan's amendment of 2026-09-27 applied).

**(A) The frames.** The large-fragment model will classify melted material by Asha's three-zone rule, for which it
needs Girin's conjugate depth patch by patch; until now it existed only as the flight-wide mean `delta_m_mean_um`.
`surface_<k>.vtp` gains two cell fields:

- **`delta_m`** [m]: the step's own conjugate depth per patch (`MeltingBody.last_delta_m`, from `spray.melt_layer`),
  carried across that step's element deaths by face id exactly as `p_w` and `tau` are (`on_current_surface`). It is
  NaN wherever the patch's closure is not Girin's — `melt_layer`'s own convention, because no conjugate depth exists
  there — and on faces the step's deaths exposed, and everywhere before the first evaluation. One difference from the
  other carried fields is deliberate: `closure`, `p_w` and `tau` come from the spray step and are written only when
  the step had film, while `delta_m` is written from every step that evaluated the surface flow.
- **`deep_thickness`** [m]: the deep liquid on the patch, m_d/(ρ_l A). It lives on the current surface (it is handed
  to the exposed faces at a death, like the film), so it is not carried. Read it as a mass per area, not a depth,
  where `deep_blob_fraction` says the liquid does not fit on its facet (fact 47).

**(B) The history and the results.** `MELT_COLUMNS` gains five columns (Task 15's audit compares the README's column
list with `MELT_COLUMNS`, so the README changes with them — sub-plan 15's amendment of this date):

- `deep_liquid_kg` — the liquid below the conjugate depth that the step's deep runoff saw: the held liquid of the
  contiguous molten chains under Girin patches plus the deep account (zone 2 of the rule, and zone 3 where it piles).
- `deep_mass_kg` — the deep account after the step: liquid the runoff has moved and that is not yet film.
- `deep_runoff_mass_kg` — cumulative mass the deep runoff took out of the elements. (A first version counted the mass
  that arrived on another patch in each of the four sub-steps, as `runoff_mass_kg` does for the film; because the
  deep account is re-transported every step, that counted the same liquid tens of times — 38 kg on a flight that
  sprays 1.1 kg — and was replaced.)
- `deep_surfaced_mass_kg` — cumulative deep liquid that became film from the top.
- `deep_blob_fraction` — the share of the deep account deeper than its patch is wide, the counterpart of
  `film_blob_fraction` (fact 27).

`melt_results` gains `deep_runoff_mass_kg`, `deep_surfaced_mass_kg` and `deep_mass_kg` (the final deep account).

Code (the tested change to `reentry_model/coupled.py`):

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -26,7 +26,8 @@
                 "kn_body", "kn_local_stag", "re_shock", "flow_branch", "p_w_stag_Pa", "phi_sonic_deg",
                 "drag_shape_factor", "frozen_mass_kg", "film_T_max_K", "film_T_mean_K", "film_frozen_fraction",
                 "unapplied_load_J", "film_blob_fraction", "rt_mass_fraction", "rt_wavelength_over_nose", "rt_growth_ms", "spray_growth_ms", "rt_region_mm", "rt_bounded_fraction", "rt_bounded_growth_ms", "molten_depth_max_mm", "molten_depth_mean_mm",
-                "delta_m_mean_um", "thick_branch_fraction", "n_dead_elements"]
+                "delta_m_mean_um", "thick_branch_fraction", "n_dead_elements", "deep_liquid_kg", "deep_mass_kg",
+                "deep_runoff_mass_kg", "deep_surfaced_mass_kg", "deep_blob_fraction"]
 PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'
 
 
@@ -106,6 +107,8 @@
                 "demise_altitude_km": float(c["altitude_km"][-1]) if history.end_reason == "demise" else None,
                 "initial_mass_kg": body.mass0, "final_mass_kg": float(c["mass_kg"][-1]), "sprayed_mass_kg": body.sprayed_mass,
                 "runoff_mass_kg": body.runoff_mass, "removed_mass_kg": body.removed_mass, "film_mass_kg": float(body.m_f.sum()),
+                "deep_runoff_mass_kg": body.deep_runoff_mass, "deep_surfaced_mass_kg": body.deep_surfaced_mass,
+                "deep_mass_kg": float(body.m_d.sum()),
                 "n_released": body.n_released, "time_of_peak_release_s": float(c["time_s"][i_peak]),
                 "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                 "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
@@ -149,9 +152,10 @@
 def write_vtk_frame(output_dir, k, body, loads):
     """field_<k>.vtu: nodal T on the active volume mesh (plus liquid fraction and element fractions when melting);
     surface_<k>.vtp: q_conv, q_rad, T per patch, and when melting the film thickness and temperature, We_s, closure,
-    local Knudsen number, wall pressure, shear, droplet radius and release rate -- the flow and spray fields of the
-    step's own evaluation, carried across that step's element deaths by face id (`MeltingBody.on_current_surface`),
-    with the defaults on faces the deaths exposed (PyVista/VTK XML)."""
+    local Knudsen number, wall pressure, shear, droplet radius, release rate and Girin's conjugate depth delta_m -- the
+    flow and spray fields of the step's own evaluation, carried across that step's element deaths by face id
+    (`MeltingBody.on_current_surface`), with the defaults on faces the deaths exposed (nan for delta_m, which is also
+    nan wherever the closure is not Girin's) -- and the thickness of the deep liquid beneath the film (PyVista/VTK XML)."""
     import pyvista as pv
     from .thermal import SIGMA_SB
     os.makedirs(output_dir, exist_ok=True)
@@ -183,6 +187,10 @@
         poly.cell_data["tau"] = carry(flow.tau if flow is not None else None, 0.0)
         poly.cell_data["r_droplet"] = carry(np.where(res.dm > 0.0, res.r, np.nan) if res is not None else None, np.nan)
         poly.cell_data["release_rate"] = carry(res.dm if res is not None else None, 0.0) / surface.areas
+        # Girin's conjugate depth per patch, nan where the patch's closure is not his (no conjugate depth exists there),
+        # and the deep liquid lying below it (amendment of 2026-10-02): the skin and the runoff layer beneath it
+        poly.cell_data["delta_m"] = carry(body.last_delta_m, np.nan)
+        poly.cell_data["deep_thickness"] = body.m_d / (body.liquid.rho * surface.areas)
     poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))
 
 
```

Tests (the tested change to `tests/test_reentry_model_coupled.py`):

```diff
--- a/tests/test_reentry_model_coupled.py
+++ b/tests/test_reentry_model_coupled.py
@@ -115,7 +115,8 @@
     grid = pv.read(os.path.join(run_dir, "vtk", "field_1.vtu"))
     assert "liquid_fraction" in grid.point_data and "phi" in grid.cell_data and grid.n_cells == int(c["n_active_elements"][20])
     poly = pv.read(os.path.join(run_dir, "vtk", "surface_1.vtp"))
-    for key in ("film_thickness", "we_s", "closure", "kn_local", "p_w", "tau", "r_droplet", "release_rate"):
+    for key in ("film_thickness", "we_s", "closure", "kn_local", "p_w", "tau", "r_droplet", "release_rate", "delta_m",
+                "deep_thickness"):
         assert key in poly.cell_data
     # the flow fields are the frame's own step's, carried across that step's element deaths: the wall pressure is positive
     # on the windward patches and nowhere above the stagnation value the history recorded from the same evaluation
@@ -146,6 +147,36 @@
         assert fh.readline().startswith("time_s,altitude_km,velocity_kms,released_mass_kg")
 
 
+def test_surface_frames_carry_the_conjugate_depth_per_patch(coarse_sphere_mesh, tmp_path):
+    """Girin's conjugate depth is the boundary between the skin the shear strips and the liquid that only runs off, and
+    the large-fragment model that reads these frames needs it patch by patch (amendment of 2026-10-02). Three seconds
+    from 69.8 km, where the windward face is under Girin's closure: the frame carries the step's own delta_m -- carried
+    across the step's element deaths like p_w and tau, nan wherever the closure is not Girin's, where no conjugate depth
+    exists -- and the thickness of the deep liquid on each patch. The body starts at 880 K so that it melts at once:
+    the closure field is the spray step's, written only where there is film, while delta_m is written from every step
+    that evaluates the flow."""
+    pytest.importorskip("cantera")
+    import pyvista as pv
+    from reentry_model import spray, surface_flow as sf
+    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
+    b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"), MASS_100MM, T0=880.0)
+    sim = simulator(b, t_max=53.0)
+    sim.advance(50.0)
+    out = os.path.join(str(tmp_path), "vtk")
+    hist = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), coupled.CoupledSettings(dt=0.5, frames_every=2, output_dir=out)).run()
+    assert hist.end_reason == "t_max" and hist.results["n_frames"] == 4               # steps 0, 2, 4 and 6
+    first, last = pv.read(os.path.join(out, "surface_0.vtp")), pv.read(os.path.join(out, "surface_3.vtp"))
+    assert np.isnan(np.asarray(first.cell_data["delta_m"])).all()        # before any step there is no evaluation yet
+    dm, closure = np.asarray(last.cell_data["delta_m"]), np.asarray(last.cell_data["closure"])
+    girin = closure == sf.CLOSURE_GIRIN
+    assert girin.sum() > 100 and (~girin).any()                           # windward patches with it, the lee without
+    assert np.isfinite(dm[girin]).all() and (dm[girin] > 1e-4).all() and (dm[girin] < 1e-3).all()   # 0.1-1 mm
+    assert np.isnan(dm[~girin]).all()
+    expected = b.on_current_surface(spray.melt_layer(b.last_flow, b.liquid)[0], np.nan)   # the step's own, carried
+    np.testing.assert_array_equal(dm, expected)
+    assert (np.asarray(last.cell_data["deep_thickness"]) >= 0.0).all()
+
+
 def test_demise_ends_the_run(coarse_sphere_mesh):
     """The lumped instant-removal device with a 60 % demise fraction: the loop stops with end_reason demise."""
     m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
```

**Measured (2026-10-02, throwaway copy of the prototype).** `test_coupled_melting_run_and_writers` (now also checking
that both fields are written) and the new `test_surface_frames_carry_the_conjugate_depth_per_patch` pass; both fail on
the unamended writer. On the frames of the 100 mm physics flight to 120 s (13 frames) `delta_m` is finite on exactly the
25 650 patch-frames with Girin's closure and NaN on all the others — no mismatch in either direction — and lies between
143 and 421 µm; it is NaN on every patch until the gate opens at 49.5 s, and on every patch of every frame of the 50 mm
flight, which never has Girin's closure. `deep_thickness` is non-zero on 600–930 patches from 60 to 110 s (139 at 120 s)
and reaches 0.69 m, a mass per area rather than a depth (fact 47) — at the default step; at 0.125 s the deep account
never exceeds 0.014 g (fact 50). With numpy's random seed fixed the history of the 50 mm flight is bit-identical with
and without the amendment, and so is the 100 mm flight's with `--deep-runoff off` (fact 49): the writer and the new
columns change no result.

## Amendment of 2026-10-03 — the molten cascade's two history columns and three results

> Part of the molten-cascade amendment (sub-plan 09's amendment of this date; facts 54–61 in `00-shared-context.md`).
> The code below is the tested change, as a diff against the prototype's `coupled.py` with the deep-runoff amendment
> of 2026-10-02 above applied.

The cascade adds no field to the frames — it moves liquid into the film, which the frames already carry — but its work
must be visible per step, because how far the surface recedes within a step, and whether the safety cap ever stops it,
is what the amendment changes. `MELT_COLUMNS` gains two columns (Task 15's audit compares the README's column list with
`MELT_COLUMNS`, so the README changes with them — sub-plan 15's amendment of this date):

- `cascade_passes` — the passes of the step's death loop in which the cascade fed something, i.e. how many layers of
  fully molten elements the surface receded through within the step beyond the first.
- `cascade_mass_kg` — cumulative mass the cascade fed to the film.

`melt_results` gains `cascade_mass_kg` (the total), `cascade_passes_max` (the most passes in any step) and
`cascade_capped_steps` (the steps in which `MAX_CASCADE_PASSES` stopped the cascade, leaving fully molten exposed
elements for the next step's feed).

Code (the tested change to `reentry_model/coupled.py`):

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -27,7 +27,7 @@
                 "drag_shape_factor", "frozen_mass_kg", "film_T_max_K", "film_T_mean_K", "film_frozen_fraction",
                 "unapplied_load_J", "film_blob_fraction", "rt_mass_fraction", "rt_wavelength_over_nose", "rt_growth_ms", "spray_growth_ms", "rt_region_mm", "rt_bounded_fraction", "rt_bounded_growth_ms", "molten_depth_max_mm", "molten_depth_mean_mm",
                 "delta_m_mean_um", "thick_branch_fraction", "n_dead_elements", "deep_liquid_kg", "deep_mass_kg",
-                "deep_runoff_mass_kg", "deep_surfaced_mass_kg", "deep_blob_fraction"]
+                "deep_runoff_mass_kg", "deep_surfaced_mass_kg", "deep_blob_fraction", "cascade_passes", "cascade_mass_kg"]
 PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'
 
 
@@ -109,6 +109,8 @@
                 "runoff_mass_kg": body.runoff_mass, "removed_mass_kg": body.removed_mass, "film_mass_kg": float(body.m_f.sum()),
                 "deep_runoff_mass_kg": body.deep_runoff_mass, "deep_surfaced_mass_kg": body.deep_surfaced_mass,
                 "deep_mass_kg": float(body.m_d.sum()),
+                "cascade_mass_kg": body.cascade_mass, "cascade_passes_max": int(c["cascade_passes"].max()),
+                "cascade_capped_steps": body.cascade_capped_steps,
                 "n_released": body.n_released, "time_of_peak_release_s": float(c["time_s"][i_peak]),
                 "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                 "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
```

Test (the tested change to `tests/test_reentry_model_coupled.py`):

```diff
--- a/tests/test_reentry_model_coupled.py
+++ b/tests/test_reentry_model_coupled.py
@@ -110,6 +110,9 @@
     assert c["mass_kg"][-1] == pytest.approx(c["mass_kg"][0] - c["sprayed_mass_kg"][-1], rel=1e-6)
     assert r["melt_onset_altitude_km"] is not None and r["spraying_onset_time_s"] >= r["melt_onset_time_s"] and r["demise_time_s"] is None
     assert r["sprayed_mass_kg"] == c["sprayed_mass_kg"][-1] and r["n_source_rows"] == len(b.source_rows) > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
+    # the molten cascade (amendment of 2026-10-03): its cumulative mass and its passes per step, in the history and results
+    assert r["cascade_mass_kg"] == c["cascade_mass_kg"][-1] and np.all(np.diff(c["cascade_mass_kg"]) >= 0.0)
+    assert r["cascade_passes_max"] == int(c["cascade_passes"].max()) and r["cascade_capped_steps"] == b.cascade_capped_steps
     assert np.all(c["convective_heat_W"] > 0.0) and c["convective_heat_W"][-1] == pytest.approx(b.last.Q_conv)
     import pyvista as pv
     grid = pv.read(os.path.join(run_dir, "vtk", "field_1.vtu"))
```

**Measured (2026-10-03, throwaway copy).** `test_coupled_melting_run_and_writers` passes with the two new assertions
(the result equals the last row of the cumulative column, which never decreases, and the most passes equals the
column's maximum); `MELT_COLUMNS` is a subset of the history columns as before. On the 100 mm flight to 120 s
at the default step (seeded, deep runoff off) `cascade_passes` is non-zero in 146 of
the 240 steps, from 33.5 s to 108 s, with a median of 3 (90th percentile 4, at most 6), `cascade_mass_kg` ends at
0.245 kg, `cascade_passes_max` is 6 and `cascade_capped_steps` 0 (157 g, at most 6 passes and 0 capped steps with the
deep runoff on). With `--molten-cascade off` the history is bit-identical with the deep-runoff amendment's in every
shared column (fact 55): the two columns are the only addition.

---

**Depends on:** Task 9 (`MeltingBody`).
**Produces, for later tasks:** melt-history columns, demise detection, VTK melt-field output, and particle output files, which Task 12's metrics and plots read from.
**Character:** plumbing/IO — straightforward glue code once Task 9 exists.
**Read before implementing:** none beyond the shared context.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan.

---

### Task 10: The coupled loop — melt columns, demise, VTK melt fields, particle files

**Files:**
- Modify (replace): `reentry_model/coupled.py`
- Test: `tests/test_reentry_model_coupled.py` (append)

**Interfaces:**
- Consumes: `MeltingBody` (Task 9: `advance(state=)`, `demised()`, `melt_stats()`, `source_rows`, `m_f`, `liquid`, `last_spray/last_flow`, `phi`, `material.liquid_fraction`), `spray.SOURCE_COLUMNS/BIN_EDGES/N_BINS/histogram` (Task 7).
- Produces: `coupled.MELT_COLUMNS` (23 names), `CoupledRun.melting` (property), `loads_at` passing `body.nose_radius()` to the heating, `melt_results(history)` (melt/spraying onsets, demise, masses, `n_released`, `r_median_um`, `n_dead_elements`, `n_source_rows`, `removed_enthalpy_J`, `melt_energy_balance_residual`), end reason `"demise"`; `write_vtk_frame` writes active cells only with `liquid_fraction`/`phi` and the surface fields `film_thickness, we_s, regime, tau, r_droplet, release_rate` when melting; `write_particles(run_dir, body, history, window=10.0) -> {"particles", "particles_summary", "size_distribution"}`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_coupled.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: the melting run, its columns, results, VTK fields and particle files

def test_coupled_melting_run_and_writers(coarse_sphere_mesh, tmp_path):
    """20 s of the 100 mm flight from 71 km with a warm body (700 K; physics heating, AA7075_range, Girin removal) on the
    coarse mesh: melt columns, results, demise bookkeeping, VTK melt fields and the three particle files."""
    pytest.importorskip("cantera")
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_range"), thermal.thermal_solver("skfem"), MASS_100MM, T0=700.0)
    sim = simulator(b, t_max=63.0)
    sim.advance(43.0)
    run_dir = str(tmp_path / "run")
    settings = coupled.CoupledSettings(dt=0.5, frames_every=20, output_dir=os.path.join(run_dir, "vtk"))
    run = coupled.CoupledRun(sim, b, heating.PhysicsHeating(), settings)
    assert run.melting
    hist = run.run()
    c, r = hist.columns, hist.results
    assert set(coupled.MELT_COLUMNS) <= set(c) and hist.end_reason == "t_max"
    assert c["mass_kg"][-1] < c["mass_kg"][0] and c["sprayed_mass_kg"][-1] > 0.0 and np.all(np.diff(c["sprayed_mass_kg"]) >= 0.0)
    assert c["removed_mass_kg"][-1] == pytest.approx(c["sprayed_mass_kg"][-1]) and c["n_dead_elements"][-1] > 0
    assert c["mass_kg"][-1] == pytest.approx(c["mass_kg"][0] - c["sprayed_mass_kg"][-1], rel=1e-6)
    assert r["melt_onset_altitude_km"] is not None and r["spraying_onset_time_s"] >= r["melt_onset_time_s"] and r["demise_time_s"] is None
    assert r["sprayed_mass_kg"] == c["sprayed_mass_kg"][-1] and r["n_source_rows"] == len(b.source_rows) > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
    assert np.all(c["convective_heat_W"] > 0.0) and c["convective_heat_W"][-1] == pytest.approx(b.last.Q_conv)
    import pyvista as pv
    grid = pv.read(os.path.join(run_dir, "vtk", "field_1.vtu"))
    assert "liquid_fraction" in grid.point_data and "phi" in grid.cell_data and grid.n_cells == int(c["n_active_elements"][20])
    poly = pv.read(os.path.join(run_dir, "vtk", "surface_1.vtp"))
    for key in ("film_thickness", "we_s", "closure", "kn_local", "p_w", "tau", "r_droplet", "release_rate"):
        assert key in poly.cell_data
    files = coupled.write_particles(run_dir, b, hist)
    for key in ("particles", "particles_summary", "size_distribution"):
        assert os.path.isfile(files[key])
    from reentry_model import spray
    npz = np.load(files["particles"])
    assert set(npz.files) == set(spray.SOURCE_COLUMNS) and npz["dm_kg"].sum() == pytest.approx(b.sprayed_mass, rel=1e-9)
    with open(files["size_distribution"]) as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "window_start_s,window_end_s,r_lo_m,r_hi_m,dn,dm_kg" and (len(lines) - 1) % spray.N_BINS == 0
    total = [l.split(",") for l in lines[1:]][-spray.N_BINS:]
    assert sum(float(row[5]) for row in total) == pytest.approx(b.sprayed_mass, rel=1e-6)
    with open(files["particles_summary"]) as fh:
        assert fh.readline().startswith("time_s,altitude_km,velocity_kms,released_mass_kg")


def test_demise_ends_the_run(coarse_sphere_mesh):
    """The lumped instant-removal device with a 60 % demise fraction: the loop stops with end_reason demise."""
    m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    b = body.MeltingBody(m, mat, thermal.thermal_solver("skfem"), MASS_100MM,
                         settings=body.MeltSettings(removal="instant", runoff=False, demise_fraction=0.6, size_feedback="initial"))
    hist = coupled.CoupledRun(simulator(b), b, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=0.5)).run()
    assert hist.end_reason == "demise" and hist.results["end_reason"] == "demise" and hist.results["demise_altitude_km"] < 71.0
    assert 0.55 * MASS_100MM < hist.columns["mass_kg"][-1] < 0.6 * MASS_100MM and hist.results["final_mass_kg"] == hist.columns["mass_kg"][-1]
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_coupled.py -q`
Expected: the two new tests fail (`advance()` gets no `state`, no `MELT_COLUMNS`, no `write_particles`).

- [ ] **Step 3: Replace `reentry_model/coupled.py`**

```python
"""Lockstep coupling of the trajectory stepper with the thermal body (spec section 7).

Per macro step of dt: Simulator.advance(dt) (DOP853, events truncate the last step) -> aero state at the end of the
step -> heating with that state and the wall temperatures at the start of the step -> ThermalBody.advance (radiation
implicit in the solver) -> history row, and every `frames_every` steps a VTK frame (nodal T on the volume mesh,
q_conv / q_rad / T per patch on the surface). First-order operator splitting; the dt-halving test bounds its error.
Mass is constant in Step 2; the loop already carries the body's mass and the mesh so Step 3 can change both.
`body` must be a ThermalBody: CoupledRun uses `theta`, `surface`, `radiated_power()`, `surface_stats()`,
`nose_radius()` and `integrated_heat`, beyond what the `Body` protocol declares. `ConstantBody` is for `Simulator.run()` only."""
import os
import time
from dataclasses import dataclass, field

import numpy as np

from . import trajectory as tj

THERMAL_COLUMNS = ["convective_heat_W", "rad_cooling_W", "integrated_heat_J", "absorbed_heat_J",
                   "surface_T_max_K", "surface_T_min_K", "surface_T_mean_K", "T_stagnation_K", "T_back_K",
                   "T_centre_K", "q_stag_Wm2", "heating_blend_f"]
MELT_COLUMNS = ["film_mass_kg", "sprayed_mass_kg", "runoff_mass_kg", "removed_mass_kg", "melt_front_depth_max_mm",
                "equivalent_radius_mm", "n_active_elements", "spraying_area_m2", "theta_cr_deg", "n_released",
                "released_mass_kg", "r_median_um", "r_max_um", "closure_fraction_girin", "closure_fraction_couette_slip",
                "closure_fraction_couette_fm", "rt_active", "removed_enthalpy_J", "film_thickness_max_mm",
                "film_thickness_mean_mm", "nose_radius_mm", "transverse_radius_mm", "fitted_nose_radius_mm",
                "kn_body", "kn_local_stag", "re_shock", "flow_branch", "p_w_stag_Pa", "phi_sonic_deg",
                "drag_shape_factor", "frozen_mass_kg", "film_T_max_K", "film_T_mean_K", "film_frozen_fraction",
                "unapplied_load_J", "film_blob_fraction", "rt_mass_fraction", "rt_wavelength_over_nose", "rt_growth_ms", "spray_growth_ms", "rt_region_mm", "rt_bounded_fraction", "rt_bounded_growth_ms", "molten_depth_max_mm", "molten_depth_mean_mm",
                "delta_m_mean_um", "thick_branch_fraction", "n_dead_elements"]
PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'


@dataclass
class CoupledSettings:
    dt: float = 0.5                  # s, macro step
    frames_every: int = 0            # VTK frame every n macro steps (0: none)
    output_dir: str = None           # directory of the VTK series (required when frames_every > 0)
    frames: list = field(default_factory=list)


class CoupledRun:
    def __init__(self, sim, body, heating_model, settings=None):
        """`body` must be a ThermalBody (uses `theta`, `surface`, `radiated_power()`, `surface_stats()` and
        `integrated_heat`, beyond the `Body` protocol); `ConstantBody` is for `Simulator.run()` only."""
        self.sim, self.body, self.heating = sim, body, heating_model
        self.settings = settings or CoupledSettings()
        if self.settings.frames_every and not self.settings.output_dir:
            raise ValueError("frames_every > 0 needs an output_dir")

    def loads_at(self, t, y):
        a = self.sim.aero_state(t, y[:3], y[3:])
        return a, self.heating.evaluate(a, self.body.theta, self.body.surface_temperature(), self.body.nose_radius(),
                                        T_mean=self.body.mean_temperature())

    def row(self, t, y, a, loads):
        body = self.body
        r = self.sim.sample_row(t, y, a)
        Q = body.last.Q_conv if body.last is not None else loads.total(body.surface.areas)     # the step's applied power (the surface may have changed since)
        r.update({"convective_heat_W": Q, "rad_cooling_W": -body.radiated_power(),
                  "integrated_heat_J": body.integrated_heat, "absorbed_heat_J": body.absorbed_heat,
                  "q_stag_Wm2": loads.q_stag, "heating_blend_f": loads.blend})
        r.update(body.surface_stats())
        if self.melting:
            r.update(body.melt_stats())
        return r

    @property
    def melting(self):
        return hasattr(self.body, "melt_stats")

    def run(self):
        started = time.perf_counter()
        sim, body, s = self.sim, self.body, self.settings
        a, loads = self.loads_at(sim.t, sim.y)
        rows, states, step = [self.row(sim.t, sim.y, a, loads)], [sim.y.copy()], 0
        self._frame(step, sim.t, loads)
        while sim.end_reason is None:
            dt = sim.advance(s.dt)
            a, loads = self.loads_at(sim.t, sim.y)
            body.advance(sim.t, dt, loads, state=a)
            step += 1
            if len(loads.q_conv) != body.surface.n_patches:       # elements died this step: loads on the new surface for the record
                a, loads = self.loads_at(sim.t, sim.y)
            rows.append(self.row(sim.t, sim.y, a, loads))
            states.append(sim.y.copy())
            self._frame(step, sim.t, loads)
            if self.melting and body.demised():
                sim.end_reason = "demise"
        self._collection()
        columns = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
        history = tj.History(columns, np.array(states), sim.end_reason)
        history.results = sim.results(history, sim.nfev, time.perf_counter() - started)
        history.results.update(self.thermal_results(history))
        if self.melting:
            history.results.update(self.melt_results(history))
        return history

    def melt_results(self, history):
        c, body = history.columns, self.body
        first = lambda mask: (float(c["time_s"][mask][0]), float(c["altitude_km"][mask][0])) if mask.any() else (None, None)
        t_melt, h_melt = first(c["removed_mass_kg"] + c["film_mass_kg"] > 0.0)
        t_spray, h_spray = first(c["sprayed_mass_kg"] > 0.0)
        i_peak = int(np.argmax(c["released_mass_kg"]))
        return {"melt_onset_time_s": t_melt, "melt_onset_altitude_km": h_melt, "spraying_onset_time_s": t_spray,
                "spraying_onset_altitude_km": h_spray, "demise_time_s": float(c["time_s"][-1]) if history.end_reason == "demise" else None,
                "demise_altitude_km": float(c["altitude_km"][-1]) if history.end_reason == "demise" else None,
                "initial_mass_kg": body.mass0, "final_mass_kg": float(c["mass_kg"][-1]), "sprayed_mass_kg": body.sprayed_mass,
                "runoff_mass_kg": body.runoff_mass, "removed_mass_kg": body.removed_mass, "film_mass_kg": float(body.m_f.sum()),
                "n_released": body.n_released, "time_of_peak_release_s": float(c["time_s"][i_peak]),
                "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
                "removed_enthalpy_J": body.removed_enthalpy, "frozen_mass_kg": body.frozen_mass,
                "melt_energy_balance_residual": body.energy_balance_residual()}

    def thermal_results(self, history):
        c, body = history.columns, self.body
        i_surf, i_mean, i_q = int(np.argmax(c["surface_T_max_K"])), int(np.argmax(c["temperature_K"])), int(np.argmax(c["convective_heat_W"]))
        return {
            "peak_surface_T_K": float(c["surface_T_max_K"][i_surf]), "time_of_peak_surface_T_s": float(c["time_s"][i_surf]),
            "altitude_of_peak_surface_T_km": float(c["altitude_km"][i_surf]),
            "peak_mean_T_K": float(c["temperature_K"][i_mean]), "time_of_peak_mean_T_s": float(c["time_s"][i_mean]),
            "peak_convective_heat_W": float(c["convective_heat_W"][i_q]), "time_of_peak_heating_s": float(c["time_s"][i_q]),
            "altitude_of_peak_heating_km": float(c["altitude_km"][i_q]),
            "integrated_heat_J": body.integrated_heat, "absorbed_heat_J": body.absorbed_heat, "radiated_heat_J": body.radiated_heat,
            "energy_balance_residual": body.energy_balance_residual(),
            "n_macro_steps": len(body.iterations), "mean_newton_iterations": float(np.mean(body.iterations)) if body.iterations else 0.0,
            "n_frames": len(self.settings.frames),
        }

    def _frame(self, step, t, loads):
        s = self.settings
        if not s.frames_every or step % s.frames_every:
            return
        k = len(s.frames)
        write_vtk_frame(s.output_dir, k, self.body, loads)
        s.frames.append((t, k))

    def _collection(self):
        s = self.settings
        if not s.frames:
            return
        for name in ("field", "surface"):
            entries = "".join('<DataSet timestep="{:.6g}" file="{}_{}.{}"/>\n'.format(t, name, k, "vtu" if name == "field" else "vtp")
                              for t, k in s.frames)
            with open(os.path.join(s.output_dir, name + ".pvd"), "w") as fh:
                fh.write(PVD_TEMPLATE.format(entries))


def write_vtk_frame(output_dir, k, body, loads):
    """field_<k>.vtu: nodal T on the active volume mesh (plus liquid fraction and element fractions when melting);
    surface_<k>.vtp: q_conv, q_rad, T per patch, and when melting the film thickness, We_s, regime, shear, droplet
    radius and release rate (PyVista/VTK XML)."""
    import pyvista as pv
    from .thermal import SIGMA_SB
    os.makedirs(output_dir, exist_ok=True)
    mesh, surface, T = body.mesh, body.surface, body.field()
    melting = hasattr(body, "melt_stats")
    tets = mesh.tets[mesh.active]
    cells = np.hstack([np.full((len(tets), 1), 4), tets]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(tets), pv.CellType.TETRA), mesh.points)
    grid.point_data["T"] = T
    if melting:
        grid.point_data["liquid_fraction"] = body.material.liquid_fraction(T)
        grid.cell_data["phi"] = body.phi[mesh.active]
    grid.save(os.path.join(output_dir, "field_{}.vtu".format(k)))
    faces = np.hstack([np.full((surface.n_patches, 1), 3), surface.faces]).ravel()
    poly = pv.PolyData(mesh.points, faces)
    Tf = surface.facet_mean(T)
    poly.point_data["T"] = T
    poly.cell_data["q_conv"] = loads.q_conv
    poly.cell_data["q_rad"] = body.emissivity * SIGMA_SB * (Tf ** 4 - body.T_ambient ** 4)
    poly.cell_data["T_patch"] = Tf
    if melting:
        poly.cell_data["film_thickness"] = body.m_f / (body.liquid.rho * surface.areas)
        poly.cell_data["film_T"] = body.film_temperature()
        n = surface.n_patches
        res, flow = body.last_spray, body.last_flow
        same = res is not None and res.r.size == n
        poly.cell_data["we_s"] = res.we_s if same else np.zeros(n)
        poly.cell_data["closure"] = flow.closure.astype(float) if same else np.ones(n)
        poly.cell_data["kn_local"] = flow.kn_local if same else np.full(n, np.nan)
        poly.cell_data["p_w"] = flow.p_w if same else np.zeros(n)
        poly.cell_data["tau"] = flow.tau if same else np.zeros(n)
        poly.cell_data["r_droplet"] = np.where(res.dm > 0.0, res.r, np.nan) if same else np.full(n, np.nan)
        poly.cell_data["release_rate"] = res.dm / surface.areas if same else np.zeros(n)
    poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))


def write_particles(run_dir, body, history, window=10.0):
    """particles.npz (the source table, SOURCE_COLUMNS as arrays), particles_summary.csv (per macro step) and
    size_distribution.csv (dn, dM per log bin per `window`-second window and in total). Returns the paths."""
    import csv
    from . import spray
    os.makedirs(run_dir, exist_ok=True)
    rows = np.array(body.source_rows, dtype=float) if body.source_rows else np.zeros((0, len(spray.SOURCE_COLUMNS)))
    npz = os.path.join(run_dir, "particles.npz")
    np.savez_compressed(npz, **{name: rows[:, i] for i, name in enumerate(spray.SOURCE_COLUMNS)})
    c = history.columns
    summary = os.path.join(run_dir, "particles_summary.csv")
    with open(summary, "w", newline="") as fh:
        w = csv.writer(fh)
        keys = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um", "theta_cr_deg",
                "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
        w.writerow(keys)
        for i in range(len(history)):
            w.writerow(["{:.9g}".format(c[k][i]) for k in keys])
    dist = os.path.join(run_dir, "size_distribution.csv")
    edges = spray.BIN_EDGES
    t = rows[:, 0] if len(rows) else np.zeros(0)
    t_end = float(c["time_s"][-1])
    starts = np.arange(0.0, t_end + 1e-9, window)
    with open(dist, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["window_start_s", "window_end_s", "r_lo_m", "r_hi_m", "dn", "dm_kg"])
        for t0 in list(starts) + [None]:
            if t0 is None:
                sel = np.ones(len(rows), dtype=bool)
                lo, hi = 0.0, t_end
            else:
                lo, hi = t0, t0 + window
                sel = (t >= lo) & (t < hi)
            n_hist, m_hist = spray.histogram(rows[sel, 12], rows[sel, 13], rows[sel, 14]) if sel.any() else (np.zeros(spray.N_BINS), np.zeros(spray.N_BINS))
            for j in range(spray.N_BINS):
                w.writerow(["{:.9g}".format(lo), "{:.9g}".format(hi), "{:.9g}".format(edges[j]), "{:.9g}".format(edges[j + 1]),
                            "{:.9g}".format(n_hist[j]), "{:.9g}".format(m_hist[j])])
    return {"particles": npz, "particles_summary": summary, "size_distribution": dist}
```


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

