# ---------------------------------------------------------------------------------------------------------- Task 12
emit(r'''### Task 12: Comparison metrics, plots and visualisation of the melting run

**Files:**
- Modify: `reentry_model/sesam_io.py` (two edits), `reentry_model/compare.py` (append + two constants)
- Modify (replace): `reentry_model/viz.py`
- Test: `tests/test_reentry_model_compare.py` (append), `tests/test_reentry_model_viz.py` (append)

**Interfaces:**
- Consumes: the melting reference files (Task 11), histories with `MELT_COLUMNS` (Task 10), VTK series with the melt fields (Task 10).
- Produces: `sesam_io.Reference.mass`, `.thickness`; `compare.MELT_PLOT_NAMES` (7, with `closures_time.png` carrying the closure area fractions and both Knudsen numbers on a twin log axis), `DEMISE_FRACTION = 0.01`, `has_melt(history, reference=None)`, `melt_metrics(history, reference) -> dict` (`mass.max_rel_m0`, `onset_altitude_diff_km`, `demise_time_diff_s`, `demise_time_rel`, ...), `plot_melt(history, outdir, title, reference=None, size_distribution_csv=None) -> paths`; `viz.render_frame(..., melting=False)`, `render_film_frame`, `render_section_frame(..., melting=False)`, `animate(..., melting=False)`, `animate_section(..., melting=False)`, `animate_film(run_dir, history, radius, fps, animation, stills)`, `still_marks(history, melting=False)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_compare.py`:

''')
emit(block(tail_from("tests/test_reentry_model_compare.py", "# ---------------------------------------------------------------------------------------------------------------\n# Step 3: the melting reference")))
emit(r'''
Append to `tests/test_reentry_model_viz.py`:

''')
emit(block(tail_from("tests/test_reentry_model_viz.py", "# ---------------------------------------------------------------------------------------------------------------\n# Step 3: overlays")))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py -q`
Expected: the five new tests fail (`Reference` has no `mass`, no `has_melt`, `render_frame` has no `melting`).

- [ ] **Step 3: Edit `reentry_model/sesam_io.py`**

In the `Reference` dataclass, after `integrated_heat: np.ndarray = None      # J, integral of the convective heat (SESAM's integrated_heat_J)` add:

```python
    mass: np.ndarray = None                 # kg (SESAM's mass_kg; decreases while the sphere melts)
    thickness: np.ndarray = None            # m (SESAM's thick_mm: the shell thickness of the melting sphere, Step 3 facts)
```

and in `load_reference`, after the `convective_heat=..., integrated_heat=optional("integrated_heat_J"),` line add:

```python
        mass=optional("mass_kg"), thickness=col("thick_mm", 1e-3) if "thick_mm" in rows[0] else None,
```

- [ ] **Step 4: Extend `reentry_model/compare.py`**

Replace the two constant lines `THERMAL_PLOT_NAMES = (...)` / `CONTINUUM_KN = 0.01` with:

```python
THERMAL_PLOT_NAMES = ("heating_time.png", "temperature_time.png", "integrated_heat.png")
MELT_PLOT_NAMES = ("mass_time.png", "mass_altitude.png", "mass_budget.png", "spraying_time.png", "closures_time.png",
                   "droplet_size_time.png", "size_distribution.png")
CONTINUUM_KN = 0.01
DEMISE_FRACTION = 0.01
```

and append to the file:

''')
emit(block(tail_from("reentry_model/compare.py", "# ---------------------------------------------------------------------------------------------------------------\n# Step 3: mass loss")))
emit(r'''
- [ ] **Step 5: Replace `reentry_model/viz.py`**

''')
emit(file_block("reentry_model/viz.py"))
emit(r'''
- [ ] **Step 6: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py tests/test_reentry_model_sesam_io.py -q`
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add reentry_model/sesam_io.py reentry_model/compare.py reentry_model/viz.py tests/test_reentry_model_compare.py tests/test_reentry_model_viz.py
git commit -m "Add the mass-loss comparison, the melt plots and the melting overlays of the videos (Step 3 Task 12)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')

# ---------------------------------------------------------------------------------------------------------- Task 13
emit(r'''### Task 13: Command line

**Files:**
- Modify (replace): `reentry_model/cli.py`
- Test: `tests/test_reentry_model_cli.py` (append)

**Interfaces:**
- Consumes: everything above.
- Produces: `run` flags `--melt off|on`, `--material` (names or path; default `AA7075_nomelt`, `AA7075_range` with `--melt on`), `--removal`, `--runoff`, `--rarefied-shear`, `--we-critical`, `--kr`, `--kt`, `--prism-layers` (default 4 with melting, 0 otherwise), `--layer-thickness` (mm), `--demise-fraction`, `--particles/--no-particles`, `--size-feedback current|initial` (default `current` with `girin`, `initial` with `instant`; it governs the Knudsen length, the nose-cap radius and the drag shape factor together), `--gamma-pm` (default 1.15), `--kn-body-shock` (default 0.01), `--k-scale`, `--consistent-mass`; run names end in `_melt-<removal>`; the JSON `settings` carry the melt settings and the liquid properties, `results` the melt results, `files` the particle files and melt plots (and `film`); `compare` handles melting histories; `model_run_name(..., heating_name=None, melt=None)`; `_fmt`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_cli.py`:

''')
emit(block(tail_from("tests/test_reentry_model_cli.py", "# ---------------------------------------------------------------------------------------------------------------\n# Step 3: melting flags")))
emit(r'''
- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: the new tests fail (`--melt` unknown).

- [ ] **Step 3: Replace `reentry_model/cli.py`**

''')
emit(file_block("reentry_model/cli.py"))
emit(r'''
- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: all pass (the bookkeeping device test runs the whole 100 mm flight on the coarse mesh in ~10 s).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/cli.py tests/test_reentry_model_cli.py
git commit -m "Add the melting flags, wiring and summaries to the command line (Step 3 Task 13)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---
''')
