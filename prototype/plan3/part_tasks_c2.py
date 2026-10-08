# ---------------------------------------------------------------------------------------------------------- Task 11
emit(r'''### Task 11: The two melting SESAM references

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

''')
emit(block(tail_from("tests/test_reentry_model_data.py", "MELTING_REFERENCE_NAMES = [")))
emit(r'''
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
''')
