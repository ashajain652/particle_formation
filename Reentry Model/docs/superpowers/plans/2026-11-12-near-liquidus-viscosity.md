# Near-liquidus viscosity Implementation Plan (revised)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the composite viscosity law the run command's default for melting skins: a suspension law on the grains'
envelope fraction below 10 % solid, and Li et al.'s 489.88-per-second fit from 10 % to half solid, both independent of
shear rate. Keep today's law as `--rheology bridge`, identical to the last bit, and probe the saved d050 state under both.

**Architecture:** The law is a pure function in `reentry_model/shared/rheology.py`. One helper in `skin.py` dispatches
on `SkinSettings.rheology`, and both `cell_viscosity` and `skinbody.slurry_modes` use it. Nothing else in the solvers
changes, because the law ignores the shear rate. The command line selects the law, and a probe script compares the laws
on the saved state's skins alone.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy, pytest (with pytest-xdist for the unit tier).

**Spec:** `Reentry Model/docs/superpowers/specs/2026-11-12-near-liquidus-viscosity-design.md` (revised, commit
93591c4). Every executor reads it; sections are cited as "spec 3.2" and so on. This plan replaces the withdrawn plan of
the first design (commits 3d0037d and d8ca92e).

## Global Constraints

- Work in the worktree `.claude/worktrees/near-liquidus-viscosity` on the branch `near-liquidus-viscosity` (off `main`
  f1ab736). Commit locally. Never push, merge or commit to `main`. Never use a bare `git stash`.
- Run everything from `Reentry Model/` with `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`.
- Never edit `Reentry Model/prototype/` or anything in another worktree. The d050 state in the subplan20 worktree is
  read only.
- The present law gives today's numbers to the last bit, under `--rheology bridge` and under the library default.
  `SkinSettings().rheology` is `"bridge"`, and the run command defaults to `"composite"` (spec 6). Existing tests keep
  passing unchanged, except the command-line tests that Task 3 pins to `--rheology bridge`.
- All the viscosity physics, including the choice between the laws by name, lives in
  `reentry_model/shared/rheology.py`. The level-set architecture will use the same law (Asha, 2026-11-12). The
  skyscrapers code only passes its settings' law name and the cell state, and holds no viscosity formula or law branch
  of its own.
- `shared/rheology.py` imports nothing from the rest of `reentry_model`; `tests/test_reentry_model_shared.py` checks
  this. The present names keep their meaning and numbers: `li_viscosity`, `viscosity`, `LI_A`, `LI_B`, `LI_N`,
  `LI_C`, `GAMMA_MAX`, `BRIDGE_LOW`, `BRIDGE_HIGH` and `RIGID_SOLID_FRACTION`.
- The constants, verbatim from spec 6:
  - `LI_BROKEN_RATE` = 489.88, `LI_BROKEN_A` = 0.04346, `LI_BROKEN_B` = 3.23;
  - `LI_SOLID_MIN` = 0.1, `SUSPENSION_S` = 3.0, `SUSPENSION_PHI_M` = 0.55;
  - `INTERNAL_SOLID_FRACTION` = 0.5;
  - `LAWS` = ("composite", "bridge"), `DEFAULT_LAW` = "composite".
- Wherever the law is described, its assumptions and uncertainties are stated as spec 3.7 gives them (Asha,
  2026-11-12). This covers the module docstring, `model_assumptions.md`, fact 196, and the analysis scripts' reports.
  The jump at 10 % solid is marked unresolved. No command-line option is added beyond `--rheology`. In particular,
  none is added for the internal solid fraction.
- `--slurry-spray` keeps its default (off), and `slurry_modes`' closure is unchanged (Asha, 2026-11-12).
- Bad arguments exit 2 with a message (`parser.error`). Outputs go under `Reentry Model/output/` (git-ignored).
- Commit messages are a scope prefix and a sentence on the behaviour that changed, ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. A state saved before this change has a `SkinSettings` without `rheology`. Unpickled, it must read `"bridge"` through
   the dataclass's class default, never `"composite"`. Tested in Task 2.
2. The law's boundaries: f_s exactly 0, just below and exactly at 0.1, and exactly at 0.5. Padded cells sit at 300 K,
   where f_s = 1, so they are infinite, and `cell_viscosity` masks them as now. Tested in Tasks 1 and 2.
3. `internal_fraction=None` reads `INTERNAL_SOLID_FRACTION` at call time, so the probe's brackets take effect when it
   sets the module constant. Tested in Task 1.
4. Under `bridge`, `slurry_modes` keeps its clip at `GAMMA_MAX` and today's numbers; under `composite` the rate it
   records is not clipped. Tested in Task 2.
5. Resuming an old state with no `--rheology` keeps the bridge law and records it in the JSON; with `--rheology
   composite` it switches. Tested in Task 3.

---

### Task 1: The composite law, its docstring and its analysis script

**Files:**
- Modify: `reentry_model/shared/rheology.py`
- Create: `analysis/near_liquidus_viscosity.py`
- Test: `tests/test_reentry_model_rheology.py`

**Interfaces:**
- Produces:
  - `suspension_viscosity(f_s, mu_liquid, internal_fraction=None) -> ndarray`, in Pa s. It gives spec 3.2's law, with
    phi = f_s / f_i: `mu_liquid` at f_s <= 0, and `inf` where phi >= `SUSPENSION_PHI_M`.
  - `composite_viscosity(f_s, mu_liquid, internal_fraction=None) -> ndarray`, in Pa s. It gives
    `suspension_viscosity` below `LI_SOLID_MIN`, `LI_BROKEN_A * exp(LI_BROKEN_B * f_s)` from `LI_SOLID_MIN` up to
    `RIGID_SOLID_FRACTION`, and `inf` at or above it.

  Both broadcast `f_s` and return float arrays. `internal_fraction=None` reads `INTERNAL_SOLID_FRACTION` at call time.
  An internal fraction at or below 0.1 / 0.55, or above 1, raises ValueError.
  - `law_viscosity(law, T, f_s, gamma, mu_liquid, f_s_bridge) -> ndarray`, in Pa s:
    - for `"bridge"`, the result of `viscosity(T, f_s, gamma, mu_liquid, f_s_bridge)`, unchanged;
    - for `"composite"`, `composite_viscosity(f_s, mu_liquid)` broadcast to `np.shape(T)`, with `gamma` ignored;
    - any other `law` raises ValueError.
  - `rate_cap(law) -> float`: the shear rate a law's caller holds a recorded rate at. It is `GAMMA_MAX` for
    `"bridge"` and `inf` for `"composite"`, and any other `law` raises ValueError. With it, `slurry_modes` needs no law
    branch.
  - `bridge_solid_fraction(material) -> float`: `1.0 - float(material.liquid_fraction(np.array([BRIDGE_LOW]))[0])`,
    the expression both skyscrapers call sites use today. It uses only the passed object's `liquid_fraction`, so the
    module imports nothing.

- [ ] **Step 1: Write the failing tests**

```python
import numpy as np, pytest
from reentry_model.shared import rheology as R
MU = 1.3e-3
FS = np.array([0.01, 0.025, 0.038, 0.057, 0.0975, 0.1 - 1e-9])

@pytest.mark.parametrize("fi, row", [(0.5, (1.062, 1.165, 1.265, 1.431, 1.906, 1.943)),
                                     (1.0, (1.031, 1.079, 1.122, 1.191, 1.356, 1.367)),
                                     (0.25, (1.129, 1.367, 1.630, 2.168, 5.022, 5.400))])
def test_below_ten_percent_matches_spec_table(fi, row):           # spec 3.6
    np.testing.assert_allclose(R.composite_viscosity(FS, MU, internal_fraction=fi) / MU, row, rtol=1e-3)

def test_li_broken_fit_from_ten_percent():
    f = np.array([0.1, 0.2, 0.3, 0.4, 0.4999])
    np.testing.assert_allclose(R.composite_viscosity(f, MU), 0.04346 * np.exp(3.23 * f), rtol=1e-12)
    np.testing.assert_allclose(R.composite_viscosity(np.array([0.1, 0.2, 0.3, 0.4, 0.45, 0.4999]), MU) / MU,
                               (46.18, 63.78, 88.10, 121.7, 143.0, 168.0), rtol=1e-3)

def test_limits_and_jump():
    assert R.composite_viscosity(0.0, MU) == MU and R.composite_viscosity(-0.01, MU) == MU
    assert np.isinf(R.composite_viscosity(0.5, MU)) and np.isinf(R.composite_viscosity(1.0, MU))
    below, at = R.composite_viscosity(0.1 - 1e-12, MU), R.composite_viscosity(0.1, MU)
    assert at / below == pytest.approx(46.18 / 1.943, rel=1e-3)                            # spec 3.5: a factor of 24
```

Also write these tests:

- `test_internal_fraction_validation`: `internal_fraction=0.18` and `1.01` each raise ValueError; `0.1 / 0.55 +
  1e-6` and `1.0` are accepted.
- `test_internal_fraction_read_at_call_time`: after `monkeypatch.setattr(R, "INTERNAL_SOLID_FRACTION", 0.25)`,
  `composite_viscosity(0.057, MU) / MU` is 2.168 to 1e-3.
- `test_suspension_viscosity_diverges_at_packing`: `suspension_viscosity(0.275, MU, 0.5)` is `inf` (phi = 0.55).
- `test_law_viscosity_dispatch`. On (3, 4) arrays of T between 895 and 910 K, with f_s from the AA7075 Scheil material
  and a gamma array:
  - `law_viscosity("bridge", ...)` is `np.array_equal` to `viscosity(...)`;
  - `law_viscosity("composite", ...)` equals `composite_viscosity(f_s, MU)` with shape (3, 4), and is unchanged when
    gamma is replaced by 1e4;
  - `law_viscosity("li", ...)` raises ValueError;
  - `rate_cap("bridge") == GAMMA_MAX`, `rate_cap("composite") == inf`, and `rate_cap("li")` raises ValueError.
- `test_bridge_solid_fraction`: for the AA7075 Scheil material it equals
  `1 - float(mat.liquid_fraction(np.array([906.4]))[0])` exactly (0.09754 to 1e-4).

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_rheology.py -q`
Expected: the new tests FAIL with `AttributeError`; the existing ones pass.

- [ ] **Step 3: Implement the constants and functions, and rewrite the module docstring**

The docstring describes both laws:

- `composite`, region by region. Each region gives its sources and statuses as in spec 4's table, then spec 3.7's
  assumptions and uncertainties in full, numbered as there. The jump at 10 % solid is marked unresolved.
- `bridge`, the present text, kept.

It names `DEFAULT_LAW` as the run command's default and `"bridge"` as `SkinSettings`' library default.
`INTERNAL_SOLID_FRACTION` carries a one-line comment saying it is assumed.

- [ ] **Step 4: Write `analysis/near_liquidus_viscosity.py`**

The script imports the package's functions and copies none of them. It writes into `output/near_liquidus_viscosity/`:

- `law.png`: viscosity over the liquid's against solid fraction from 0 to 0.5, on a log axis. It shows:
  - the law for f_i = 0.5, 1 and 0.25;
  - Li et al.'s five per-rate fits at their measured solid fractions, as points for reference;
  - the present law at 0 and 367 per second (`rheology.viscosity`, with the bridge's f_s_bridge from the AA7075 Scheil
    material);
  - shading for the region below 10 % solid, labelled "assumed: no aluminium data (spec 3.7)", and a mark at 10 %
    labelled "jump unresolved".
- `law_tables.txt`: spec 3.6's table, printed from the functions.

- [ ] **Step 5: Run the tests and the script**

Run: `"$PY" -m pytest tests/test_reentry_model_rheology.py tests/test_reentry_model_shared.py -q` and
`"$PY" analysis/near_liquidus_viscosity.py`
Expected: all PASS; the script exits 0 and its table matches spec 3.6.

- [ ] **Step 6: Commit**

```bash
git add reentry_model/shared/rheology.py tests/test_reentry_model_rheology.py analysis/near_liquidus_viscosity.py
git commit -m "reentry_model: the composite near-liquidus viscosity law beside the present bridge law -- a Mori-Ototake suspension on the grains' envelope fraction below 10 % solid (internal solid fraction 0.5, assumed) and Li et al.'s 489.88 1/s fit from 10 % to half solid, independent of shear rate, with each region's sources, assumptions and uncertainties in the docstring; and analysis/near_liquidus_viscosity.py for its figure and table" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 2: The skins follow the selected law

**Files:**
- Modify: `reentry_model/skyscrapers/skin.py` (`SkinSettings`, a new helper, `cell_viscosity`)
- Modify: `reentry_model/skyscrapers/skinbody.py` (`slurry_modes`)
- Test: `tests/test_reentry_model_skin.py`, `tests/test_reentry_model_skinbody.py`

**Interfaces:**
- Consumes: `rheology.law_viscosity`, `rheology.rate_cap`, `rheology.bridge_solid_fraction`, `rheology.LAWS`
  (Task 1).
- Produces:
  - `SkinSettings.rheology: str = "bridge"`, validated against `rheology.LAWS` in `__post_init__` (ValueError);
  - `skin.viscosity_of(sk, T, f_s, gamma) -> ndarray`, in Pa s: a thin wrapper that returns
    `rheology.law_viscosity(sk.settings.rheology, T, f_s, gamma, mat.liquid.mu, rheology.bridge_solid_fraction(mat))`,
    with `mat = sk.material`. It holds no formula and no law branch of its own (Global Constraints).

- [ ] **Step 1: Write the failing tests**

In `tests/test_reentry_model_skin.py`:

- `test_settings_carry_the_law`:
  - `SkinSettings().rheology == "bridge"`;
  - `SkinSettings(rheology="composite")` is accepted;
  - `SkinSettings(rheology="li")` raises ValueError;
  - a `SkinSettings` with `rheology` deleted from its `__dict__`, pickled and unpickled, reads `"bridge"`.
- `test_cell_viscosity_follows_composite`: on an existing test's melting skin with composite settings,
  `cell_viscosity(sk, 0.0)` equals `cell_viscosity(sk, 1e4)` exactly. On valid cells it equals
  `composite_viscosity(1 - sk.props.liquid_fraction(sk.T), mu_l)`, and it is `inf` beyond each column's cells.
- `test_composite_momentum_single_pass`: `momentum_step` on an existing test's mixed flowing skin under composite
  settings returns `book.passes == 1`. Its impulses add up to the momentum change to round-off, as the existing
  momentum tests require.

In `tests/test_reentry_model_skinbody.py`, beside the existing slurry-mode test near line 1115:

- `test_slurry_modes_composite`: the same setup with composite settings gives `s["mu"]` equal to
  `composite_viscosity(f_s, mu_l)` at the exposed cell's solid fraction to 1e-12. The recorded `s["gamma"]` is not
  clipped: choose loads so that at least one exceeds 367 per second. The existing bridge test stays unchanged.

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py -q -m "not reference" -k "law or composite"`
Expected: FAIL (`SkinSettings` has no `rheology`).

- [ ] **Step 3: Implement**

- Add the field and its validation.
- Add `viscosity_of`. Make `cell_viscosity` call it in place of its direct `rheology.viscosity` call, keeping its
  non-melting and validity masks.
- In `slurry_modes`, replace the `viscosity` lambda's direct call and its own `f_bridge` line with
  `skin_mod.viscosity_of(sk, T, f_s, gamma)`. In `layer()`, replace `np.minimum(shear, rheology.GAMMA_MAX)` with
  `np.minimum(shear, rheology.rate_cap(sk.settings.rheology))`. That is the same operation with the same constant under
  `"bridge"`, and no clip under `"composite"`, with no law branch in `skinbody.py`.
- Update both docstrings.

- [ ] **Step 4: Run the skin and skinbody tests**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py tests/test_reentry_model_rheology.py -q -n 4 -m "not reference"`
Expected: all PASS, the existing ones unchanged.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/skyscrapers/skin.py reentry_model/skyscrapers/skinbody.py tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py
git commit -m "reentry_model: the skins' settings carry the viscosity law (the bridge by default in the library, so old states and existing tests keep the present law), and the cells' viscosity and the exposed slurry's mode both follow it" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 3: The command line, and the regression of the present law

**Files:**
- Modify: `reentry_model/cli.py`:
  - the run option and `skin_settings`;
  - `skin_name_suffix` and `run_name_of`;
  - the run JSON's `skin` record;
  - the resume option, its `changed` handling, saved-arguments update, printout and `thermal_info`;
  - the usage text;
  - the `--slurry-spray` help text at line 473.
- Test: `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_restart.py`

**Interfaces:**
- Consumes: `SkinSettings.rheology` (Task 2); `rheology.LAWS`, `DEFAULT_LAW` and `INTERNAL_SOLID_FRACTION` (Task 1).
- Produces:
  - `run --rheology {composite,bridge}` (default `rheology.DEFAULT_LAW`), passed through `skin_settings(args)`;
  - `skin_name_suffix(..., rheology="composite")`, which adds `_rheology-bridge`;
  - `info["skin"]["rheology"]` and `info["skin"]["internal_solid_fraction"]` (`rheology.INTERNAL_SOLID_FRACTION`),
    always recorded for skin runs;
  - `resume --rheology {composite,bridge}` (default `None`, which keeps the saved flight's law and is `"bridge"` for an
    old state through the class default). It sets `run.body.settings.rheology`, which is the same object as
    `run.body.skins.settings`, and records `changed["--rheology"]`.

- [ ] **Step 1: Write the failing tests**

In `test_reentry_model_cli.py`:

- `test_rheology_option`:
  - the default parses to `"composite"`, and `skin_settings` carries it;
  - `--rheology li` exits 2;
  - with `--surface-model elements`, or without `--melt on`, the option is accepted and the run name has no mark.
- `test_rheology_run_names`:
  - `skin_name_suffix` defaults give no mark, and `rheology="bridge"` adds `_rheology-bridge`;
  - a default melting run's full name equals what the existing name tests assert for the same arguments.
- `test_rheology_json`: a skin run's `info["skin"]` holds `rheology` and `internal_solid_fraction == 0.5`, using the
  existing JSON test's harness.
- `test_slurry_spray_help_names_the_law`: the `--slurry-spray` help no longer contains "no data".

In `test_reentry_model_restart.py`:

- `test_resume_rheology`: on the existing resume test's saved state, with `rheology` deleted from the saved
  `SkinSettings.__dict__` to stand for an old state:
  - with no option, the resumed body's law is `"bridge"` and the JSON records it;
  - `--rheology composite` sets it, and records `changed["--rheology"] == "composite"`.

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py -q -k "rheology or slurry_spray_help"`
Expected: FAIL (unknown option `--rheology`).

- [ ] **Step 3: Implement**

Follow the existing patterns for options that resume can change: `--slurry-spray`, `--leeward-test`. The new
`--slurry-spray` help names the law in use. For example: "it sprays only within about a kelvin of the liquidus, at the
viscosity `--rheology` selects (the composite law's suspension part below 10 % solid is assumed, design 3.7)". Keep the
default off.

- [ ] **Step 4: Pin command-line tests that assert the present law's numbers**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py -q -n 4`

Every test that fails only because the command line's default law changed gets `--rheology bridge` added to its
arguments. Any other failure is a bug to fix, not to pin. List each pinned test in the commit message.

Expected: all PASS.

- [ ] **Step 5: Run the regression of the present law**

Run the short melting smoke flight from `CLAUDE.md` ("a short melting smoke flight (65 km, 0.125 s)") four times. Use
a scratch folder `S` under `output/near_liquidus_viscosity/regression/`.

1. In the old code, the `skyscrapers-fixes` worktree at f1ab736: `cd` into its `Reentry Model/`, and run with
   `--outdir S/old --name smoke`, then again with `--slurry-spray on --name smoke_spray`.
2. In this branch, run the same two with `--rheology bridge` added, into `S/new`.

Run: `"$PY" analysis/skyscrapers_port/same_results.py compare S/old S/new --report`
Expected: every `.csv` identical. In the `.json`, only `skin.rheology`, `skin.internal_solid_fraction` and the keys
recording the command's arguments may differ. Anything else is a defect.

Then run the smoke flight under the default (composite) into `S/composite` and confirm it exits 0. Record its
differences from `S/new` in the ledger as expected.

- [ ] **Step 6: Run the unit tier**

Run: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4`
Expected: all pass, 1 skipped (FEniCSx). Record the counts and the time.

- [ ] **Step 7: Commit**

```bash
git add reentry_model/cli.py tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py
git commit -m "reentry_model: run and resume take --rheology composite|bridge, composite by default for melting skins, recorded in the run JSON with the internal solid fraction and marked in the run name for the bridge; a resumed old state keeps the bridge law unless told otherwise; the --slurry-spray help names the law" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 4: The probe of the d050 state (spec 8, level 4)

**Files:**
- Create: `analysis/near_liquidus_probe.py`
- Output: `output/near_liquidus_viscosity/probe/` and `output/near_liquidus_viscosity/probe.md` (git-ignored)

**Interfaces:**
- Consumes:
  - `cli.load_state` (or the loader it wraps);
  - `run.sim.aero_state(t, y[:3], y[3:])` at the saved step, as `coupled.CoupledRun.loads_at` calls it;
  - `run.body.flow_forcing`;
  - `skin.Flow` and `skin.flow_substeps` with `release=None`;
  - `skin.momentum_step`, `Skins.depth_above`;
  - `skinbody.girin_plan(..., slurry=True)`, following the body's own call near `skinbody.py:915`, and
    `skinbody.slurry_modes`;
  - `SkinSettings.rheology` (Task 2) and `rheology.INTERNAL_SOLID_FRACTION` (Task 1).

- [ ] **Step 1: Write the probe script**

`analysis/near_liquidus_probe.py --state <pkl> --case {bridge,composite,composite_fi1,composite_fi025} --flow-time 0.1
--out <dir> [--single-column]` does the following:

- Loads the state. Sets `run.body.settings.rheology` to the case's law. For the two bracket cases, sets
  `rheology.INTERNAL_SOLID_FRACTION` to 1.0 or 0.25 in the script's own process.
- Recomputes the gas loads with the current code: `flow_forcing` at the saved step's aero state. It records the
  loads' difference from the body's saved forcing (`run.body._forcing`, if present): the median and largest change of
  `tau` and `dp_ds`.
- Builds `skin.Flow` from those loads, the body's columns, top areas, `v_hat`, `rho_l` and `sigma`, with
  `release=None` and adaptive cycles.
- Calls `flow_substeps` once per 10 ms heat sub-step until `--flow-time`. It runs no heat, no 3D conduction and no
  trajectory step.
- At the start, at the end and every 20 ms, on the windward columns (angle from the stagnation point below 60 degrees),
  evaluates and records the following, writing them to `<out>/<case>.npz` and `<out>/<case>.json`:
  - The liquid's and the exposed slurry's modes, under two sets of release rules: the saved ones (gates "before",
    slip band off), and gates `2026-11-10` with the slip band on. Use `run.body.spray.spray_gates` and
    `run.body.spray.slip_band`, restoring them after. For each: the Weber number against the critical value, which
    gate decided the mode, and the release rate.
  - The closure check (spec 8): the momentum step's top-cell speeds, and the slurry mode recomputed through
    `slurry_modes` with its flowing depth `F` ended at the 10 % solid depth. That depth is
    `sk.depth_above(T10)`, where `T10` is the temperature at 10 % solid on the material's Scheil curve (906.36 K for
    AA7075), found by inverting `liquid_fraction`.
  - The top three cells' f_s and viscosity over mu_l, and the stress each cell carries (the face stresses, as
    `_momentum_columns`' `rates()` computes them).
  - The flow sub-step counts and the wall time per flow sub-step.
- `--single-column` picks the windward column with the median flowing depth among those holding slurry cells. It runs
  `momentum_step` alone on `sk.rows([i])` with that column's loads until steady (relative change below 1e-9, at most
  2000 sub-steps of 1 ms), and writes the profiles of velocity, viscosity and stress.

- [ ] **Step 2: Measure the cost before running**

The state is
`/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/.claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/d050__handdown-off.state.pkl`.
Run the `bridge` and `composite` cases for 0.01 s of flow and time them. If a case would take more than about 2 hours
at 0.1 s, shorten `--flow-time` for all cases to fit, not below 0.03 s, and record the ruling in the ledger. Check
the free memory before running cases in parallel.

- [ ] **Step 3: Run the four cases and the single column**

Run the four cases as parallel processes, as memory allows, into `output/near_liquidus_viscosity/probe/`. Then run
`--single-column` for each case. Each must exit 0.

- [ ] **Step 4: Write `output/near_liquidus_viscosity/probe.md`**

The report states at its top that the law below 10 % solid and the internal fraction are assumed (spec 3.7). For each
case and each rule set it gives the quantities of Step 1, as medians with 10th to 90th percentiles over the windward
columns, and the single column's profiles. It then answers five questions with numbers:

1. How viscous are the top cells, and how fast does the surface move, under each law?
2. Under which law, internal fraction and rule set does the liquid or the slurry reach its spraying gate, and which gate
   decides?
3. How large is the closure artefact? Compare the mode's surface speed and release against the momentum step's
   top-cell speed and against the mode with its layer ended at 10 % solid.
4. How far did the recomputed loads differ from the saved forcing?
5. What did a flow sub-step cost under each law?

- [ ] **Step 5: Commit the script**

```bash
git add analysis/near_liquidus_probe.py
git commit -m "analysis: near_liquidus_probe.py runs the skins' flow alone from a saved state under each viscosity law and internal fraction, evaluating the liquid's and the slurry's spraying modes under two sets of release rules without releasing metal, the slurry mode's closure against the momentum step's speeds and against a layer ended at 10 % solid, and a single column to steady state" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 5: Documentation

**Files:**
- Modify:
  - `docs/model_assumptions.md`;
  - `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md`;
  - `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`;
  - `../README.md` and `CLAUDE.md`.

**Interfaces:**
- Consumes: Task 3's regression and unit-tier counts, and Task 4's `probe.md`.

- [ ] **Step 1: Write the documents (spec 9)**

- `model_assumptions.md`: a new section on the near-liquidus viscosity (shared physics). It gives:
  - the law region by region;
  - spec 4's table of sources and statuses;
  - spec 3.7 in full, with the jump at 10 % solid marked unresolved;
  - the probe's headline numbers, each with its measurement.
- The layered-skin design: a dated amendment (2026-11-12) at each passage that states the present law. These are the
  table near line 67, section 6.4 near line 343, the slurry mode near lines 611-617, line 958, and assumption 5 near
  line 1082 on the falling stress and the fixed point's cap. The amendment says the run command's default is now the
  composite law, which does not depend on the rate, and that the bridge is `--rheology bridge`.
- `00-shared-context.md` fact 196. It records:
  - the change, the commits and the option;
  - the regression (the present law identical on the smoke flights);
  - the unit tier's counts;
  - the probe's answers to Task 4's five questions;
  - that `--slurry-spray` stays off by default, pending the Girin session's new spraying law.
- `README.md` and `CLAUDE.md`: the option and the changed default in the command surface.

- [ ] **Step 2: Check the documents against the measurements**

Every number written in Step 1 must appear in `probe.md`, the Task 3 logs or the spec. Run
`"$PY" -m pytest tests/test_reentry_model_rheology.py -q`.
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add docs ../README.md CLAUDE.md
git commit -m "docs: the composite near-liquidus viscosity is the run command's default -- a model_assumptions section with each region's sources, assumptions and uncertainties, the layered-skin design amended, and fact 196 with the regression, the unit tier and the d050 probe's results; the bridge kept as --rheology bridge" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
