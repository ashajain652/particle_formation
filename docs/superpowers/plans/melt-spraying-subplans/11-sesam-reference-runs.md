# Sub-plan: Task 11 — The two melting SESAM references

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 6030–6129). Read `00-shared-context.md` first — pay particular attention to the Global Constraints' rule that DRAMA's databases and the wrapper scripts (`sphere_reentry.py`, `sphere_sweep.py`) must never be modified.
**Depends on:** nothing from this Step 3 plan's code — this task generates data using the project's existing, unmodified wrapper scripts.
**Produces, for later tasks:** the two committed melting SESAM reference runs that Task 12 compares the coupled model's output against.
**Character:** data generation, not code — fundamentally different in kind from every other task in this plan.
**Read before implementing:** Measured facts 1 and 2 in the shared context: the two melting references hollow a sphere at fixed outer geometry, and SESAM's `--removal instant` device must remove liquid from every element as it forms, not just surface elements — get this detail right or the reference data itself will be wrong.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan. Since this task has no code dependency on Tasks 1–9, note in the plan that it can be executed at any time, in parallel with everything else, even though Task 12 won't consume its output until much later.

---

### Task 11: The two melting SESAM references

**Files:**
- Create: `data/reference_runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv/.json`, `data/reference_runs/sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind.csv/.json`
- Test: `tests/test_reentry_model_data.py` (append), `tests/test_reentry_model_aero.py` (one guard)

**Interfaces:**
- Consumes: the wrapper `sphere_reentry.py` (untouched) with DRAMA's `drama-AA7075` on the static US76 table, winds off — the break-off state of the existing references.
- Produces: the two reference runs (CSV + JSON) that Tasks 12–14 compare against: 100 mm melt onset 71.005 km at 43.55 s, mass 0 at 66.476 km / 67.315 s (70 rows); 50 mm onset 77.104 km at 178.55 s, mass 0 at 73.094 km / 191.55 s (194 rows).

- [ ] **Step 1: Generate the references with the wrapper (or reuse the ones generated on 2026-09-20)**

The runs already exist under `sphere_sweep_output/reference_AA7075/runs/` (generated 2026-09-20 with exactly the commands below; `sphere_sweep_output/` is git-ignored). If they are missing, generate them (each takes under a minute; the wrapper needs `drama_env`):

```bash
mkdir -p sphere_sweep_output/reference_AA7075/runs sphere_sweep_output/reference_AA7075/raw
"$PY" sphere_reentry.py --velocity 7.5 --altitude 77.500133 --temperature 300 --diameter 100 --flight-path-angle -0.959331 \
    --heading 347.168296 --lat 29.546067 --lon -82.134333 --epoch 2024-08-01T12:53:07 --atmosphere static --no-wind \
    --outdir sphere_sweep_output/reference_AA7075/runs --raw-dir sphere_sweep_output/reference_AA7075/raw --keep-raw
"$PY" sphere_reentry.py --velocity 7.5 --altitude 115 --temperature 300 --diameter 50 --flight-path-angle -0.959331 \
    --heading 347.168296 --lat 29.546067 --lon -82.134333 --epoch 2024-08-01T12:53:07 --atmosphere static --no-wind \
    --outdir sphere_sweep_output/reference_AA7075/runs --raw-dir sphere_sweep_output/reference_AA7075/raw --keep-raw
```

Expected console lines: `sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind: demised after 67.3 s (uncritical); Tmax 850.0 K, 23.8 s at melt, final mass 0 kg` and `sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind: demised after 191.6 s (uncritical); Tmax 850.0 K, 13.0 s at melt, final mass 0 kg`. Then copy the four files:

```bash
cp sphere_sweep_output/reference_AA7075/runs/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.* data/reference_runs/
cp sphere_sweep_output/reference_AA7075/runs/sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind.* data/reference_runs/
```

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_reentry_model_data.py`:

```python
MELTING_REFERENCE_NAMES = ["sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind",
                           "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind"]


def test_melting_reference_runs_are_committed_and_consistent():
    """The two Step 3 melting references: DRAMA's drama-AA7075 (850 K, 400 kJ/kg) on the US76 table, winds off; the
    sphere melts completely (mass 0 at the end, 'demised')."""
    for name, onset_km, end_km in zip(MELTING_REFERENCE_NAMES, (71.005, 77.104), (66.476, 73.094)):
        csv_path = os.path.join(REFERENCE_DIR, name + ".csv")
        doc = json.load(open(os.path.join(REFERENCE_DIR, name + ".json")))
        assert doc["run_name"] == name and doc["status"] == "ok" and doc["inputs"]["material"] == "drama-AA7075"
        assert doc["inputs"]["atmosphere"] == "static" and doc["inputs"]["use_wind"] is False
        assert doc["results"]["melting_temperature_K"] == 850.0 and abs(doc["results"]["altitude_first_melt_km"] - onset_km) < 1e-3
        assert abs(doc["results"]["altitude_last_melt_km"] - end_km) < 1e-3
        rows = list(csv.DictReader(open(csv_path)))
        assert float(rows[-1]["mass_kg"]) == 0.0 and float(rows[-1]["thick_mm"]) == 0.0 and float(rows[-1]["altitude_km"]) > 60.0
        assert abs(float(rows[0]["mass_kg"]) - doc["inputs"]["initial_mass_kg"]) < 1e-3


def test_disc_aerothermal_table():
    """data/atdb_disc.json: ATDB_CYLINDER at zero angle of attack and the thinnest aspect ratio -- the flat-face limit
    of the sphere-to-disc family (Step 3). HTG's continuum database is modified Newtonian, so the sphere entry is
    0.49897 of the disc entry at every Mach number, and the free-molecular entries differ by only 3 %."""
    sphere = json.load(open(os.path.join(REPO_ROOT, "reentry_model", "data", "atdb_sphere.json")))
    disc = json.load(open(os.path.join(REPO_ROOT, "reentry_model", "data", "atdb_disc.json")))
    assert disc["mach"] == sphere["mach"] == [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    assert "ATDB_CYLINDER" in disc["_provenance"] and "angle of attack 0" in disc["_provenance"]
    cd_c, cd_fm = disc["cd_continuum"], disc["cd_free_molecular"]
    assert cd_c[1] == pytest.approx(1.824148) and cd_fm[1] == pytest.approx(2.202072)
    assert all(abs(s / d - 0.49897) < 1e-5 for s, d in zip(sphere["cd_continuum"], cd_c))
    assert all(0.99 < d / s < 1.08 for s, d in zip(sphere["cd_free_molecular"], cd_fm))
    assert all(v == pytest.approx(0.329561) for v in disc["heat_flux_factor_continuum"])       # Mach-independent, like the sphere's
```


In `tests/test_reentry_model_aero.py::TestDragCoefficient::test_reproduces_the_reference_drag_column` (it globs every reference CSV; the melting references demise while still hypersonic, so their low-speed list is empty) replace the last assertion

```python
        assert np.abs(low).max() <= 1.5e-3, name       # the reference prints C_D with 3 decimals
```

with

```python
        if low:                                        # the melting references (Step 3) demise while still hypersonic
            assert np.abs(low).max() <= 1.5e-3, name   # the reference prints C_D with 3 decimals
```

- [ ] **Step 3: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_data.py tests/test_reentry_model_aero.py tests/test_reentry_model_sesam_io.py -q`
Expected: all pass (the melting references' hypersonic C_D rows reproduce SESAM's drag rule like the others; `test_reference_runs_are_committed_and_consistent` still covers only the no-melt names).

- [ ] **Step 4: Commit**

```bash
git add data/reference_runs tests/test_reentry_model_data.py tests/test_reentry_model_aero.py
git commit -m "Add the two melting SESAM references (drama-AA7075, US76, winds off) (Step 3 Task 11)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

