# Near-liquidus viscosity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the skins' default slurry viscosity (Li et al.'s averaged law with the 906.4-908 K "no data" bridge) with
the composite law of the design: Li et al.'s per-rate fits where they measured, the literature elsewhere, and each
cell's shear-rate history choosing its branch of the non-monotone flow curve. Keep the present law as `--rheology bridge`,
identical to the last bit, and probe the saved d050 state under both.

**Architecture:** The law and its continuation root are pure functions in `reentry_model/shared/rheology.py`. The skins
carry a per-cell `shear_rate` through every cell operation. Under `composite`, the momentum step and the exposed slurry's
mode find each cell's rate by continuing along its flow curve from its previous rate. Under `bridge`, today's code path
runs untouched. The command line selects the law (`--rheology`, `--structure-exponent`), and a probe script compares the
laws on the saved state's skins alone.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy, pytest (with pytest-xdist for the unit tier).

**Spec:** `Reentry Model/docs/superpowers/specs/2026-11-12-near-liquidus-viscosity-design.md`. Every executor reads
it; sections are cited below as "spec 3.4" and so on.

## Global Constraints

- Work in the worktree `.claude/worktrees/near-liquidus-viscosity` on the branch `near-liquidus-viscosity` (off
  `main` f1ab736). Commit locally. Never push, merge or commit to `main`. Never use a bare `git stash`.
- Run everything from `Reentry Model/` with `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`.
- Never edit `Reentry Model/prototype/` or anything in another worktree. The d050 state in the subplan20 worktree is
  read only.
- The present law gives today's numbers to the last bit. This holds under `--rheology bridge` and under the library
  default, `SkinSettings()`, which stays `"bridge"` (ruling R1 below). Existing tests keep passing unchanged, except
  CLI tests that Task 6 pins to `--rheology bridge`.
- `shared/rheology.py` imports nothing from the rest of `reentry_model`; `tests/test_reentry_model_shared.py` checks
  this. Its present names (`li_viscosity`, `viscosity`, `LI_A`, `LI_B`, `LI_N`, `LI_C`, `GAMMA_MAX`, `BRIDGE_LOW`,
  `BRIDGE_HIGH`, `RIGID_SOLID_FRACTION`) keep their meaning and numbers.
- The law's constants, verbatim from spec 3 and 6:
  - `LI_RATES` = (61.24, 122.47, 244.94, 367.41, 489.88) per second;
  - `LI_RATE_A` = (0.66711, 0.59623, 0.33326, 0.14323, 0.04346) Pa s;
  - `LI_RATE_B` = (4.78, 3.52, 3.45, 3.66, 3.23);
  - `LI_SOLID_MIN` = 0.1, `SUSPENSION_S` = 3.0, `SUSPENSION_PHI_M` = 0.55, `EXCESS_RATE_EXPONENT` = -0.85;
  - `STRUCTURE_EXPONENTS` = (2.0, 3.0), `STRUCTURE_EXPONENT` = 2.0, `GAMMA_FLOOR` = 1.0e-3 per second;
  - `LAWS` = ("composite", "bridge"), `DEFAULT_LAW` = "composite";
  - `SHEAR_TOL` stays 1e-3.
- Wherever the law is described, it states that two things are assumed and not supported by literature (Asha,
  2026-11-12). This covers the module docstring, `model_assumptions.md`, fact 196, and any analysis script's report.
  - Above 489.88 per second, the position of the flow curve's valley and the slope n_x are assumed. No option is added
    for them.
  - Below 10 % solid, the behaviour is assumed, and the uncertainties are those of spec 3.10.
- Bad arguments exit 2 with a message (`parser.error`).
- Outputs go under `Reentry Model/output/` (git-ignored).
- Commit messages are a scope prefix and a sentence on the behaviour that changed (`reentry_model: ...`, `docs: ...`,
  `analysis: ...`), ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

### Rulings made while planning (record any later ones in the ledger)

- **R1.** `SkinSettings.rheology` defaults to `"bridge"` in the library, and the run command defaults to
  `"composite"` (`rheology.DEFAULT_LAW`). This departs from spec 6's wording, which gives the dataclass `"composite"`.
  The reasons:
  - Every existing test that builds `SkinSettings()` keeps today's numbers.
  - A settings object unpickled from a state saved before this change, which has no `rheology` in its `__dict__`,
    reads the dataclass's class default `"bridge"`, as spec 6 requires.
  - With `"composite"` as the class default, such a state would silently read the new law.

  The spec was amended to match in the commit that added this plan.
- **R2.** The continuation root treats the law as monotone between sub-knots (section "Task 2"). Below about 1 % solid
  and 1 Pa, the falling pieces have shallow interior dips: at f_s = 0.006 the stress is 0.737 Pa at 244.94 per
  second, 0.730 at 281 per second and 0.757 at 367.41 per second. A sub-knot method may step over them. Cells in those
  dips carry stresses far below any gas load, so the error is immaterial. The oracle tests exclude stresses below
  2 Pa at f_s below 0.01.
- **R3.** Under `composite`, the sideways transport carries `m * shear_rate` in a second solve with the same matrix,
  so `m`, `H`, `S` and `U` are untouched to the last bit. Under `bridge` the rates stay zero and that solve does not
  run.

## Review Focus

1. A state saved before this change (`SkinSettings` without `rheology`, `Skins` without `shear_rate`) loads as the
   bridge law with zero rates, never silently as composite. Tested in Task 3.
2. `momentum_step` solves the live columns on `sk.rows(live)` copies; their converged rates must be written back, and
   the anchored columns' rates must be zero. Tested in Task 4.
3. Degenerate inputs to `continued_rate`: a stress of zero, a liquid cell (f_s = 0), a rigid or padded cell (f_s at
   least 0.5; padding sits at 300 K, so f_s = 1). Each returns a finite rate, never NaN. Tested in Task 2.
4. Mass-weighted carrying with zero-mass or padded cells (merge, transport, new top cells) leaves no NaN or inf in
   `shear_rate`. Tested in Task 3.
5. `--structure-exponent 3` together with the bridge law, whether given or carried by an old state on resume, exits 2.
   Tested in Task 6.

---

### Task 1: The composite law

**Files:**
- Modify: `reentry_model/shared/rheology.py` (constants, three functions, the module docstring)
- Test: `tests/test_reentry_model_rheology.py`

**Interfaces:**
- Produces:
  - `suspension_viscosity(f_s, mu_liquid) -> ndarray`, which is spec 3.2's eta_inf;
  - `composite_viscosity(f_s, gamma, mu_liquid, p=STRUCTURE_EXPONENT, n_x=None) -> ndarray`, in Pa s;
  - `flow_stress(f_s, gamma, mu_liquid, p=STRUCTURE_EXPONENT, n_x=None) -> ndarray`, which is composite_viscosity times
    gamma, in Pa.

  `n_x=None` means `EXCESS_RATE_EXPONENT`, read at call time, not bound when the function is defined. The analysis
  scripts can then change the module constant for their sensitivity runs, and the skins' calls follow; Task 7 relies
  on this. `p` is passed explicitly by every caller.

  All three broadcast their array-like arguments and return float arrays. `composite_viscosity` is `inf` at
  f_s >= `RIGID_SOLID_FRACTION` and `mu_liquid` at f_s <= 0. It evaluates the law at max(gamma, `GAMMA_FLOOR`)
  (spec 3.7). `flow_stress` multiplies by gamma itself, so below the floor the stress falls linearly to zero.

- [ ] **Step 1: Write the failing tests**

```python
import numpy as np, pytest
from reentry_model.shared import rheology as R
MU = 1.3e-3

def test_composite_reproduces_li_per_rate_fits():
    for g, a, b in zip(R.LI_RATES, R.LI_RATE_A, R.LI_RATE_B):
        for f in (0.1, 0.2, 0.3, 0.4, 0.4999):
            assert R.composite_viscosity(f, g, MU) == pytest.approx(a * np.exp(b * f), rel=1e-12)

@pytest.mark.parametrize("f, row", [
    (0.025, (6.053e5, 1707, 52.72, 23.62, 3.879, 1.294, 1.109)),
    (0.057, (3.147e6, 8870, 269.7, 118.4, 15.75, 2.312, 1.349)),
    (0.0975, (9.207e6, 2.595e4, 786.8, 344.2, 43.95, 4.636, 1.819)),
    (0.1, (9.685e6, 2.730e4, 827.7, 362.0, 46.18, 4.818, 1.854)),
    (0.2, (1.562e7, 4.404e4, 1335, 511.1, 63.78, 6.705, 2.616)),
    (0.45, (5.159e7, 1.454e5, 4410, 1211, 143.0, 18.79, 9.889))])
def test_composite_matches_spec_table(f, row):       # spec 3.8, p = 2
    rates = (1e-3, 1.0, 61.24, 244.94, 489.88, 1e4, 1e5)
    np.testing.assert_allclose(R.composite_viscosity(f, np.array(rates), MU) / MU, row, rtol=1e-3)

def test_composite_cube_row():
    rates = np.array((1e-3, 1.0, 61.24, 244.94, 489.88, 1e4, 1e5))
    np.testing.assert_allclose(R.composite_viscosity(0.025, rates, MU, p=3.0) / MU,
                               (1.513e5, 427.6, 13.99, 6.713, 1.779, 1.132, 1.086), rtol=1e-3)
```

Also write tests for the following:

- `test_composite_limits`:
  - eta equals `MU` exactly at f_s = 0 and at f_s = -0.01, at any rate;
  - eta is `inf` at f_s = 0.5 and 0.9;
  - eta at 1e9 per second over `suspension_viscosity(0.1, MU)` is within 1e-3 of 1;
  - `suspension_viscosity(np.array([0.02, 0.05, 0.1, 0.2, 0.5]), MU) / MU` is (1.0605, 1.1714, 1.3667, 1.9412,
    17.5) to 1e-3;
  - eta at 1e-4 per second equals eta at 1e-3 per second exactly;
  - `flow_stress(0.1, 1e-4, MU)` is a tenth of `flow_stress(0.1, 1e-3, MU)`.
- `test_composite_is_continuous`:
  - at f_s 0.1 - 1e-12 against 0.1, and at each rate edge 61.24 and 489.88 plus or minus 1e-9 relative, the values
    agree to 1e-8 relative;
  - for f_s in `np.linspace(0, 0.499, 200)` at rates 1, 300 and 1e4, eta has no relative jump above 1e-2 between
    neighbouring f_s;
  - eta is never below `suspension_viscosity`.
- `test_flow_curve_peaks_and_valleys` (spec 3.9; absolute tolerance 0.1 Pa): `flow_stress` at 244.94 and at 489.88
  per second is 115.3 and 29.4 Pa at f_s 0.1, 37.7 and 10.0 at 0.057, and 7.5 and 2.5 at 0.025; with p = 3 it is 21.6
  and 6.0 at 0.057. At f_s 0.49 the knot stresses at the five Li rates are (425.05, 409.75, 442.61, 316.26, 103.64)
  Pa, to 0.1.
- `test_n_x_keyword`: with `n_x=-1.0`, `flow_stress(0.1, g, MU) - suspension_viscosity(0.1, MU) * g` is the same at
  1e3 and at 1e4 per second, to 1e-9 relative. The excess stress is then constant.

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_rheology.py -q`
Expected: the new tests FAIL with `AttributeError` (no `composite_viscosity`); the existing ones pass.

- [ ] **Step 3: Implement**

Write the constants (Global Constraints) and the three functions in `reentry_model/shared/rheology.py` exactly as spec
3.2 to 3.7 define them. The construction:

- Li et al.'s window interpolates ln eta linearly in ln gamma between `LI_RATES`, with `np.interp` per solid fraction.
- Outside the rates, the excess is (L at the clipped rate minus eta_inf) times (gamma / clipped rate) to the power
  n_x.
- Below `LI_SOLID_MIN`, the excess is the excess at 0.1 times (f_s / 0.1) to the power p.

Rewrite the module docstring. It describes both laws:

- `composite`, region by region, each region with its source and status as in spec 4's table: the paper, the page
  and the table or figure, and whether the value was measured, fitted or assumed. Two regions say plainly that they
  are assumed, not supported by literature:
  - above 489.88 per second: the valley's position and the slope, with the bounds of spec 3.10 (n_x from -0.4 to
    -1.6);
  - below 10 % solid: with spec 3.10's nine uncertainties, numbered as there.

  The constants `EXCESS_RATE_EXPONENT` and `STRUCTURE_EXPONENT` carry a one-line comment saying they are assumed.
- `bridge`, the present text, kept.

It also names `DEFAULT_LAW` as the run command's default and `SkinSettings`' `"bridge"` as the library default
(ruling R1).

- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_rheology.py tests/test_reentry_model_shared.py -q`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/shared/rheology.py tests/test_reentry_model_rheology.py
git commit -m "reentry_model: the composite near-liquidus viscosity law (Li et al.'s per-rate fits where they measured, a Mori-Ototake suspension in the real liquid with a gamma^-0.85 structural excess outside their rates, the excess falling as (f_s/0.1)^p below 10 % solid) beside the present bridge law, with each region's source in the docstring" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 2: The continuation root, and the law's analysis script

**Files:**
- Modify: `reentry_model/shared/rheology.py`
- Create: `analysis/near_liquidus_viscosity.py`
- Test: `tests/test_reentry_model_rheology.py`

**Interfaces:**
- Consumes: `flow_stress` and `composite_viscosity` (Task 1).
- Produces: `continued_rate(f_s, tau, gamma_from, mu_liquid, p=STRUCTURE_EXPONENT, n_x=None) -> ndarray`, spec 5.2's
  root (`n_x=None` as in Task 1), broadcast over its arguments, in per second. It returns:
  - 0 where f_s >= 0.5 or tau <= 0;
  - tau / mu_liquid where f_s <= 0;
  - otherwise a rate g with `flow_stress(f_s, g, ...)` equal to tau within 1e-9 relative.

  It never returns NaN.

- [ ] **Step 1: Write the failing tests**

```python
def test_continued_rate_hysteresis():            # f_s = 0.1: peak 115.3 Pa at 244.94, valley 29.4 Pa at 489.88
    up = R.continued_rate(0.1, 100.0, 0.0, MU)
    assert up < 244.94 and R.flow_stress(0.1, up, MU) == pytest.approx(100.0, rel=1e-9)   # structured: 115.65 1/s
    broke = R.continued_rate(0.1, 120.0, 0.0, MU)
    assert broke > 489.88                                                                    # past the peak: 36 836 1/s
    assert R.continued_rate(0.1, 100.0, broke, MU) > 489.88                                  # stays broken down: 26 980
    assert R.continued_rate(0.1, 50.0, broke, MU) > 489.88                                   # above the valley: 5 230
    assert R.continued_rate(0.1, 20.0, broke, MU) < 244.94                                   # below the valley: 0.0219

def test_continued_rate_two_peaks():             # f_s = 0.49: peaks 425.05 (61.24) and 442.61 (244.94), dip 409.75
    assert R.continued_rate(0.49, 420.0, 0.0, MU) < 61.24                                   # 56.6 1/s
    assert 122.47 < R.continued_rate(0.49, 430.0, 0.0, MU) < 244.94                         # 188.9 1/s
    assert R.continued_rate(0.49, 450.0, 0.0, MU) > 489.88                                  # 15 490 1/s

def test_continued_rate_degenerate_inputs():
    out = R.continued_rate(np.array([0.0, 0.6, 1.0, 0.1, 0.1]), np.array([13.0, 50.0, 50.0, 0.0, -1.0]),
                           np.zeros(5), MU)
    np.testing.assert_array_equal(np.isfinite(out), True)
    assert out[0] == pytest.approx(13.0 / MU) and out[1] == 0.0 and out[2] == 0.0 and out[3] == 0.0 and out[4] == 0.0
```

Also write these tests:

- `test_continued_rate_vectorised_mixed_branches`: f_s, tau and gamma_from as (4, 3) arrays mixing structured and
  broken-down starts give the same values as the cell-by-cell scalar calls, exactly.
- `test_continued_rate_against_dense_oracle`. The oracle, written in the test:
  - evaluate `flow_stress` on 4000 log-spaced rates from 1e-3 to 1e9 per second;
  - walk from the start's grid position in the direction the stress demands;
  - bisect the first bracketing interval.

  For f_s in (0.005, 0.02, 0.05, 0.0975, 0.2, 0.35, 0.49), tau in (2, 10, 30, 61, 101, 170, 300, 450) Pa, and
  gamma_from in (0, 1e6), the result matches the oracle to 1e-6 relative. Per ruling R2, cases with tau < 2 Pa at
  f_s < 0.01 are left out.

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_rheology.py -q -k continued`
Expected: FAIL with `AttributeError`.

- [ ] **Step 3: Implement `continued_rate`**

The algorithm, which the signature does not determine:

1. Build the knot list: `GAMMA_FLOOR`, then 16 log-spaced sub-knots inside each interval between consecutive values
   of (`GAMMA_FLOOR`, 61.24, 122.47, 244.94, 367.41, 489.88), the five rates themselves, and a top of 1e12 per second.
2. Evaluate `flow_stress` at every knot for every cell, giving an array of shape (cells, knots).
3. Locate gamma_from's interval.
4. If the stress at gamma_from is below tau, take the first knot above it whose stress reaches tau. If it is above
   tau, take the first knot below it whose stress falls to tau. Below `GAMMA_FLOOR` the stress is linear in gamma, so
   that segment's root is tau / eta(floor).
5. Solve inside the bracketing interval by vectorised bisection in ln gamma until the bracket's ratio is below
   1 + 1e-12 (about 60 iterations). Then return the end whose stress is nearer tau, which meets the 1e-9 test.

Each interval is treated as monotone (ruling R2).

- [ ] **Step 4: Write `analysis/near_liquidus_viscosity.py`**

The script imports the package's functions and copies none of them. It writes into
`output/near_liquidus_viscosity/`:

- `law.png`: the draft figure's two panels, redrawn from `rheology`. The draft is
  `/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/reports/near_liquidus_viscosity_composite.py`,
  in the main checkout and not committed.
- `law_tables.txt`: spec 3.8's table and spec 3.9's peaks and valleys, printed from the functions.

It takes `--n-x` (default `EXCESS_RATE_EXPONENT`) for the sensitivity of spec 11. Run it once and check that the
printed table equals spec 3.8.

- [ ] **Step 5: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_rheology.py tests/test_reentry_model_shared.py -q` and
`"$PY" analysis/near_liquidus_viscosity.py`
Expected: all PASS; the script exits 0 and its table matches spec 3.8.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/shared/rheology.py tests/test_reentry_model_rheology.py analysis/near_liquidus_viscosity.py
git commit -m "reentry_model: the continuation root of the composite law (a cell's new shear rate is the first root of its flow curve in the direction the stress demands, so a structured cell breaks down only past its peak and a broken-down one recovers only below its valley), and analysis/near_liquidus_viscosity.py for the law's figure and tables" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 3: The skins carry each cell's shear rate, and the settings carry the law

**Files:**
- Modify: `reentry_model/skyscrapers/skin.py` (`SkinSettings`, `Skins` and its cell operations, `transport`,
  `cell_viscosity`)
- Test: `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: `rheology.composite_viscosity`, `rheology.LAWS`, `rheology.STRUCTURE_EXPONENTS` (Task 1).
- Produces:
  - `SkinSettings.rheology: str = "bridge"` and `SkinSettings.structure_exponent: float = 2.0`, both validated in
    `__post_init__` against `rheology.LAWS` and `rheology.STRUCTURE_EXPONENTS` (ValueError);
  - `Skins.shear_rate`, an (n, K) float array in per second, zero where a cell is not valid;
  - `cell_viscosity(sk, gamma)`, which under `"composite"` returns `rheology.composite_viscosity(f_s, gamma,
    mat.liquid.mu, sk.settings.structure_exponent)` and under `"bridge"` does exactly what it does today;
  - `Skins.__setstate__`, which gives a `shear_rate` of zeros of `m`'s shape to a state pickled without one.

- [ ] **Step 1: Write the failing tests**

- `test_settings_carry_the_law`:
  - `SkinSettings().rheology == "bridge"` and `SkinSettings().structure_exponent == 2.0`;
  - `SkinSettings(rheology="composite", structure_exponent=3.0)` is accepted;
  - `rheology="li"` and `structure_exponent=2.5` each raise ValueError;
  - a `SkinSettings` whose `__dict__` has `rheology` and `structure_exponent` deleted, pickled and unpickled, reads
    `"bridge"` and 2.0.
- `test_shear_rate_travels_with_the_cells`, on a composite-settings `Skins` filled with known masses and
  per-cell rates:
  - `_split` halves keep their parent's rate;
  - `_merge` gives the mass-weighted mean;
  - `take_top`, `take_top_mass` (a partly taken cell keeps its rate) and `take_bottom` keep the rates of the cells
    they leave;
  - `add_bottom` gives the bottom cell the mass-weighted mean with arriving metal at rate 0;
  - `state()`/`set_state()` round-trip `shear_rate`, and a 7-element state from before the change sets it to zero;
  - `rows(idx)` copies it.
- `test_shear_rate_padding_and_zero_mass`: after a merge that empties a column's trailing cells and after `_width`
  padding, `shear_rate` is finite everywhere and zero on every invalid cell.
- `test_transport_carries_shear_rate_mass_weighted`. Reuse an existing transport test's two-column setup with
  composite settings:
  - the body's sum of m times rate is conserved to round-off;
  - a cell that only receives from one donor ends at that donor's rate;
  - the same transport under composite settings leaves `m`, `H`, `T` and `U` bit-identical (`np.array_equal`) to the
    run under bridge settings: the second solve of ruling R3 touches nothing else.
- `test_old_skins_unpickle_with_zero_rates`: a `Skins` with `shear_rate` deleted from its `__dict__`, pickled and
  unpickled, has `shear_rate` equal to zeros of `m`'s shape.

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q -k "law or shear_rate or old_skins"`
Expected: FAIL (`SkinSettings` has no `rheology`; `Skins` has no `shear_rate`).

- [ ] **Step 3: Implement**

- Add the two `SkinSettings` fields and their validation.
- Add `shear_rate` to `Skins`. It is carried by every operation the way `S` is carried, intensive and mass-weighted:
  - `__init__` (zeros) and `_install`, which takes a `rate` argument with zeros when not given, as `U` defaults to
    rest;
  - `_width`, `_put_rows`, `_merge`, `_split`, `take_top`, `take_top_mass`, `add_bottom` and `take_bottom`
    (`exchange_bottom` goes through the last two);
  - `state()` (an eighth element) and `set_state()` (an absent eighth element gives zeros);
  - `rows()` and `__setstate__`.
- `landing_conduct`'s block needs no rate; it writes back only `H`, `T` and `S`.
- In `transport`, when `sk.settings.rheology == "composite"`, solve a second right-hand side, `cell(sk.m *
  sk.shear_rate)`, with the same matrix and `_solve` call (ruling R3). Install each cell's value divided by its new
  mass, including the landed top cells.
- In `cell_viscosity`, dispatch on `sk.settings.rheology`.

- [ ] **Step 4: Run the skin tests**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py tests/test_reentry_model_rheology.py -q -n 4`
Expected: all PASS, the existing ones unchanged.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/skyscrapers/skin.py tests/test_reentry_model_skin.py
git commit -m "reentry_model: the skins carry each cell's shear rate through every cell operation (split halves keep it, merges and the sideways transport mass-weight it, metal from the 3D model arrives at rest), and their settings carry the viscosity law, the bridge by default in the library so old states and existing tests keep the present law" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 4: The momentum step under the composite law

**Files:**
- Modify: `reentry_model/skyscrapers/skin.py` (`MomentumBook`, `momentum_step`, `_momentum_columns`, `FlowBook`,
  `flow_substeps`)
- Test: `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: `rheology.continued_rate` (Task 2); `Skins.shear_rate` and `SkinSettings.rheology` (Task 3).
- Produces:
  - `COMPOSITE_SHEAR_PASSES = 5`, the pass cap under composite, measured in Task 7;
  - `MomentumBook.unconverged: int = 0`, the flowing cells whose viscosity moved by more than `SHEAR_TOL` in the
    last pass, computed under both laws;
  - `FlowBook.momentum_passes: int = 0` (the most momentum passes in any flow sub-step) and `FlowBook.unconverged:
    int = 0` (summed over the flow sub-steps), carried by `FlowBook.add`.

- [ ] **Step 1: Write the failing tests**

All on one column of six cells, each 0.1 mm thick, with solid fractions (0.025, 0.057, 0.1, 0.1, 0.3, 0.6). Set the
temperatures by inverting the AA7075 Scheil curve. The last cell is rigid and anchored. Settings are composite with
p = 2. The gas shear is 100 Pa along a tangent and there is no push. Run `momentum_step` sub-steps of 1 ms until the
velocities change by less than 1e-9 relative, at most 2000.

- `test_composite_column_branches_by_history`:
  - starting from zero rates, cells 0 and 1 end above 489.88 per second (their peaks, 7.5 and 37.7 Pa, are below
    100 Pa);
  - cells 2 and 3 (peak 115.3 Pa) end below 244.94;
  - cell 4 (peak 229.8 Pa, valley 56.1 Pa) ends below 244.94;
  - every flowing cell's stored rate satisfies `flow_stress` = 100 Pa to 1e-3, the steady Couette stress;
  - repeated with cells 2 to 4 starting at a rate of 1e4, those three end above 489.88 per second;
  - the rigid cell's rate is 0.
- `test_composite_momentum_books_close`: on the same column, and on a two-column skin where one column is all rigid,
  the impulses (push, shear and wall) add up to the momentum change to round-off, as the existing momentum tests
  require. The rigid column's rates are 0, and the moving column's rates are written back from the `sk.rows(live)`
  copy, so they are non-zero.
- `test_breakdown_in_a_pass_stands`: one 1 ms sub-step on a cell at f_s 0.1 starting structured, with a gas shear of
  116 Pa. This is just above the 115.3 Pa peak, a load that breaks it down in the first pass. The stored rate after
  the sub-step is above 489.88 per second whatever `COMPOSITE_SHEAR_PASSES` is (monkeypatch it to 4 and to 5).
- `test_bridge_momentum_keeps_rates_zero`: under bridge settings, `momentum_step` on the composite test's column
  leaves `shear_rate` all zeros. The existing momentum tests, which assert today's values under the default (bridge)
  settings, pass unchanged; they are the regression of the bridge path.

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q -k "composite or breakdown or keeps_rates_zero"`
Expected: the composite tests FAIL (the rates stay 0).

- [ ] **Step 3: Implement the composite path in `_momentum_columns`**

The bridge path keeps today's lines unchanged. Under composite:

- Start `eta` from `cell_viscosity(sk, sk.shear_rate)`, with crust cells at `CRUST_VISCOSITY`. Do not start from rest.
- Each pass solves as today, then takes `face` (the mean stress over a flowing cell's two faces, as `rates()`
  computes it). It sets the cell's `gamma = rheology.continued_rate(f_s, face, gamma_prev_pass, mu_l, p)`, where
  `gamma_prev_pass` is the stored rate in the first pass and the last pass's rate after that (spec 5.3).
- It then sets `eta` to the law at that rate.
- Stop when the relative change is at most `SHEAR_TOL`, or after `COMPOSITE_SHEAR_PASSES` passes.
- Store the last rates in `sk.shear_rate` for flowing cells, and zero for every other cell.

In `momentum_step`, write `sub.shear_rate` back for the live columns and zero the others, as it does for `U`.

Compute `unconverged` from the last pass's change under both laws. Then carry `momentum_passes` and `unconverged`
into `FlowBook` from each flow sub-step's `MomentumBook` in `flow_substeps`.

- [ ] **Step 4: Run the skin tests**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q -n 4`
Expected: all PASS, the existing ones unchanged.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/skyscrapers/skin.py tests/test_reentry_model_skin.py
git commit -m "reentry_model: under the composite law the momentum step starts each cell from its stored shear rate and sets its new rate by continuing along its flow curve to the stress it carries, keeping a breakdown found in a pass, and records how many cells had not converged" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 5: The exposed slurry's mode under the composite law

**Files:**
- Modify: `reentry_model/skyscrapers/skinbody.py` (`slurry_modes`)
- Test: `tests/test_reentry_model_skinbody.py`

**Interfaces:**
- Consumes: `rheology.continued_rate`, `rheology.composite_viscosity` (Tasks 1-2); `Skins.shear_rate` (Task 3).
- Produces: no new names. Under composite, `slurry_modes`' returned `mu` and `gamma` per column come from the
  continuation.

- [ ] **Step 1: Write the failing test**

`test_slurry_modes_composite_uses_history`: build the setup of the existing test at
`tests/test_reentry_model_skinbody.py` near line 1115 with composite settings, and run it twice.

- With the exposed cells' stored rates at zero, every returned `gamma` satisfies `flow_stress(f_s, gamma) == mu *
  gamma` to 1e-9 relative, and `mu == composite_viscosity(f_s, gamma)`.
- With stored rates of 1e5, each `gamma` is at least the structured run's, and `mu` at most.
- No `gamma` is clipped at `GAMMA_MAX`: at least one exceeds 367 per second in the broken-down run.

The existing test, which uses the default (bridge) settings, stays unchanged.

- [ ] **Step 2: Run it to make sure it fails**

Run: `"$PY" -m pytest tests/test_reentry_model_skinbody.py -q -k slurry_modes_composite`
Expected: FAIL (the rate is held at 367 per second).

- [ ] **Step 3: Implement**

Under `sk.settings.rheology == "composite"`:

- Start the fixed point from `sk.shear_rate[rows, k]`.
- In each pass take the layer's stress as `mu * layer(mu)[2]` (the unclipped shear). Then set `gamma =
  continued_rate(f_s, stress, gamma_prev_pass, mu_l, p)` and `mu = composite_viscosity(f_s, gamma, mu_l, p)`.
- Use the same `SHEAR_PASSES` and `SHEAR_TOL`, and do not clip at `GAMMA_MAX`.

Under bridge, keep today's lines as they are. Leave the `T_slurry` gate unchanged (spec 5.4). Update the function's
docstring.

- [ ] **Step 4: Run the skinbody tests**

Run: `"$PY" -m pytest tests/test_reentry_model_skinbody.py -q -n 4 -m "not reference"`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/skyscrapers/skinbody.py tests/test_reentry_model_skinbody.py
git commit -m "reentry_model: under the composite law the exposed slurry's spraying mode starts from the cell's stored shear rate and continues along its flow curve, with no cap at 367 1/s" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 6: The command line, and the regression of the present law

**Files:**
- Modify: `reentry_model/cli.py` (run and resume options, `skin_settings`, `skin_name_suffix`, `run_name_of`, the run
  JSON's `skin` record, resume's `changed` handling and printout, the usage text)
- Test: `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_restart.py`

**Interfaces:**
- Consumes: `SkinSettings.rheology`/`structure_exponent` (Task 3); `rheology.LAWS`, `DEFAULT_LAW`,
  `STRUCTURE_EXPONENTS`.
- Produces:
  - `run --rheology {composite,bridge}` (default `rheology.DEFAULT_LAW`) and `--structure-exponent {2,3}` (default
    2), both passed through `skin_settings(args)`;
  - `skin_name_suffix(..., rheology="composite", structure_exponent=2.0)`, which adds `_rheology-bridge` and
    `_structexp3`;
  - `info["skin"]["rheology"]` and `info["skin"]["structure_exponent"]`, always recorded for skin runs;
  - `resume --rheology`/`--structure-exponent`, defaulting to the saved flight's, which reads `"bridge"` and 2 for an
    old state through ruling R1, applied to `run.body.settings` (the same object as `run.body.skins.settings`) and
    recorded in `changed`.

- [ ] **Step 1: Write the failing tests**

In `test_reentry_model_cli.py`:

- `test_rheology_options`:
  - the defaults parse to composite and 2;
  - `skin_settings` carries them;
  - `--rheology bridge --structure-exponent 3` exits 2;
  - `--structure-exponent 2.5` exits 2 (argparse choices);
  - with `--surface-model elements`, or without `--melt on`, both options are accepted and ignored (spec 6), and the
    run name carries no mark for them.
- `test_rheology_run_names`:
  - `skin_name_suffix` defaults give no mark;
  - `rheology="bridge"` adds `_rheology-bridge`;
  - `structure_exponent=3.0` adds `_structexp3`;
  - a default melting run's full name is unchanged from f1ab736's for the same arguments. Take that name from the
    existing name tests.
- `test_rheology_json`: a run's `info["skin"]` holds both keys, using the existing JSON test's harness.

In `test_reentry_model_restart.py`:

- `test_resume_rheology`: on the existing resume test's saved state, with its `SkinSettings` stripped of the two
  attributes to stand for an old state:
  - with no option, the resumed body's law is bridge;
  - `--rheology composite` sets it, records `changed["--rheology"]`, and the state's `Skins.shear_rate` is zeros;
  - `--structure-exponent 3` alone (the law staying bridge) exits 2.

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py -q -k rheology`
Expected: FAIL (unknown option `--rheology`).

- [ ] **Step 3: Implement the options, names, JSON and resume handling**

Follow the existing options' patterns. `--slurry-spray` and `--leeward-test` show the run, name, JSON, resume and
printout handling.

- [ ] **Step 4: Pin CLI tests that assert the present law's numbers**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py -q -n 4`

Every test that fails because the command line's default law changed gets `--rheology bridge` added to its
arguments. A test whose failure is anything other than a changed melting result is a bug to fix, not to pin. List
each pinned test in the commit message.

Expected: all PASS.

- [ ] **Step 5: Run the regression of the present law**

Run the short melting smoke flight, the CLAUDE.md command with `--t-max 0.125` and `--frames-every 0`, twice:

1. With `--slurry-spray on` added, in the old code: the `skyscrapers-fixes` worktree, which is at f1ab736. Use `cd`
   into its `Reentry Model/` and `--outdir <scratch>/old --name smoke_spray`.
2. In this branch, with the same arguments plus `--rheology bridge`, into `<scratch>/new`.

Repeat both without `--slurry-spray on` under `--name smoke`. Then run:

Run: `"$PY" analysis/skyscrapers_port/same_results.py compare <scratch>/old <scratch>/new --report`
Expected: every `.csv` identical. In the `.json`, only `settings`/`skin.rheology` and `skin.structure_exponent` (and
any key recording the command's arguments) may differ. Anything else is a defect.

Then run the smoke flight under the default (composite) into `<scratch>/composite` and confirm it exits 0. Its
differences are expected; record them in the ledger.

- [ ] **Step 6: Run the unit tier**

Run: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4`
Expected: all pass, 1 skipped (FEniCSx). Record the counts and the time.

- [ ] **Step 7: Commit**

```bash
git add reentry_model/cli.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py
git commit -m "reentry_model: run and resume take --rheology composite|bridge (composite by default) and --structure-exponent 2|3, recorded in the run JSON and marked in the run name when not the default; a resumed old state keeps the bridge law unless told otherwise" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 7: The probe of the d050 state (level 4)

**Files:**
- Create: `analysis/near_liquidus_probe.py`
- Output: `output/near_liquidus_viscosity/probe/` and `output/near_liquidus_viscosity/probe.md` (git-ignored)

**Interfaces:**
- Consumes:
  - `restart.load_state` (or the package's loader that `cmd_resume` uses);
  - `run.body.flow_forcing` or the body's saved forcing pair `run.body._forcing`;
  - `skin.Flow`, `skin.flow_substeps` with `release=None`, and `skin.momentum_step`;
  - `skinbody.girin_plan(..., slurry=True)`, following the body's own call near `skinbody.py:915`;
  - `FlowBook.momentum_passes`/`unconverged` (Task 4).

- [ ] **Step 1: Write the probe script (spec 8, level 4)**

`analysis/near_liquidus_probe.py --state <pkl> --case {bridge,composite_rest,composite_broken,composite_p3}
--flow-time 0.1 --out <dir>` does the following:

- Loads the state.
- Sets the skins' settings to the case's law. `composite_broken` sets every flowing cell's `shear_rate` to 1e6.
- Builds `skin.Flow` with the forcing the body holds for the saved step's end (`run.body._forcing[1]` if it exists,
  else `flow_forcing` at the run's aero state at 178.5 s; record which), the body's columns, top areas, `v_hat`,
  `rho_l` and `sigma`, `release=None` and adaptive cycles.
- Calls `flow_substeps` once per 10 ms heat sub-step until `--flow-time`, with no heat, no 3D conduction and no
  trajectory.
- At the start, at the end and every 20 ms, evaluates `girin_plan(..., slurry=True)` on the skins as they stand,
  releasing nothing. It records, on the windward columns (angle from the stagnation point below 60 degrees):
  - the liquid's and the slurry's Weber numbers against the critical value, and the release rates;
  - the top three cells' f_s, viscosity over mu_l and rate;
  - the stress each cell carries (face stress);
  - the share of slurry cells above 489.88 per second (broken down) and below 244.94 (structured);
  - the surface speeds;
  - the momentum passes and unconverged counts;
  - the flow sub-step counts and the wall time per flow sub-step.

  It writes these to `<out>/<case>.npz` and `<out>/<case>.json`.
- A `--single-column` mode picks the windward column with the median flowing depth among those holding slurry cells.
  It runs `momentum_step` alone on `sk.rows([i])` with that column's loads until steady (relative change below 1e-9,
  at most 2000 sub-steps of 1 ms), for bridge, composite from rest, composite broken down, and composite with n_x at
  -0.4 and -1.6, the literature bounds of spec 3.10 on the assumed slope. The n_x cases set
  `rheology.EXCESS_RATE_EXPONENT` in the script's own process, never in the package source; this works because of
  `n_x=None` (Task 1). It writes the profiles of velocity, viscosity, rate and stress.

- [ ] **Step 2: Measure the cost before running**

Run each case for 0.01 s of flow (`--flow-time 0.01`) and time it. If any case would take more than about 2 hours at
0.1 s, shorten `--flow-time` for all cases to fit, not below 0.03 s. Record the ruling in the ledger. Also check the
free memory before running cases in parallel.

- [ ] **Step 3: Run the four cases and the single column**

The state is
`/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/.claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/d050__handdown-off.state.pkl`.
Run the four cases as parallel processes into `output/near_liquidus_viscosity/probe/`, then
`--single-column`. Each must exit 0.

- [ ] **Step 4: Write `output/near_liquidus_viscosity/probe.md`**

The report gives, for each case, the quantities of Step 1 as medians with 10th to 90th percentiles over the windward
columns, and the single column's profiles. It then answers three questions with numbers:

1. Do the top cells break down under the stresses they actually carry?
2. Does the liquid or the slurry reach its critical Weber number under either law?
3. Did the momentum iteration converge, and if not, how many cells did not?

If more than 0.1 % of the flowing cells end a sub-step unconverged, rule on `COMPOSITE_SHEAR_PASSES` (raise it and
re-run the affected case), and record the ruling.

- [ ] **Step 5: Commit the script**

```bash
git add analysis/near_liquidus_probe.py
git commit -m "analysis: near_liquidus_probe.py runs the skins' flow alone from a saved state under each viscosity law, evaluating the liquid's and the slurry's spraying modes without releasing metal, and a single column to steady state" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 8: Documentation

**Files:**
- Modify:
  - `docs/model_assumptions.md` (a new section on the near-liquidus viscosity, shared physics);
  - `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` (the layered-skin design);
  - `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md` (fact 196);
  - `../README.md` and `CLAUDE.md`.

**Interfaces:**
- Consumes: Task 6's regression and unit-tier counts, and Task 7's `probe.md`.

- [ ] **Step 1: Write the documents (spec 9)**

- `model_assumptions.md`: a new section on the slurry viscosity. It gives:
  - the composite law region by region, with spec 4's table of sources and statuses;
  - spec 3.10 in full: the behaviour above 489.88 per second marked as assumed and not supported by literature, with
    its sensitivity table, and the nine uncertainties below 10 % solid;
  - the history rule;
  - the zero-rate floor marked as numerical;
  - the limitations of spec 10;
  - the probe's measured headline numbers, each with its measurement.
- The layered-skin design: a dated amendment (2026-11-12) at each passage that states the present law (the table near
  line 67, section 6.4 near line 343, the slurry mode near lines 611-617, line 958, and assumption 5 near line 1082).
  It says the run command's default is now the composite law, and the bridge is `--rheology bridge`.
- `00-shared-context.md` fact 196. It records:
  - the change, the commits and the options;
  - the regression (the present law identical on the smoke flights);
  - the unit tier's counts;
  - the convergence measured in the probe;
  - the probe's answers to Task 7's three questions;
  - that the behaviour above 489.88 per second and below 10 % solid is assumed (spec 3.10), and how far the single
    column's n_x bounds moved the result.
- `README.md` and `CLAUDE.md`: the two options and the default in the command surface.

- [ ] **Step 2: Check the documents against the measurements**

Every number written in Step 1 must appear in `probe.md`, the Task 6 logs or the spec. Run
`"$PY" -m pytest tests/test_reentry_model_rheology.py -q`, a sanity run after the documentation edits.
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add docs ../README.md CLAUDE.md
git commit -m "docs: the composite near-liquidus viscosity is the run command's default (model_assumptions section with each region's source, the layered-skin design amended, fact 196 with the regression, the unit tier and the d050 probe's results), the bridge kept as --rheology bridge" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
