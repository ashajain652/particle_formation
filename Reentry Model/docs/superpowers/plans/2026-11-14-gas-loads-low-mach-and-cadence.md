# Gas Loads: Low-Mach Pressure, Solve Cadence and Reference Flights Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the interim gas loads' pressure step at 90 degrees below about Mach 4 by using Daniels &
Schwartzentruber's pressure at low Mach, solve the boundary layer at most once per `--dt` of flight instead of at every
continuum step, skip the solve when no part of the wall can be heated, and fly the default loads in the Step 2 reference
tests.

**Architecture:** All physics changes are in `reentry_model/shared/gas_loads.py` (the provider) and
`reentry_model/shared/heating.py` (`SimilarityHeating`); `AeroState` gains the flight time so the provider can hold a
solution; `cli.py` builds the provider with the hold; the reference test file gains the default-loads flights.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, Cantera, pytest (+ pytest-xdist).

**Spec:** Asha's decisions of 2026-11-14, verbatim: "Decision 1: let's use the Daniel's pressure law from their paper for
low mach speeds. Decision 2: I think the aerodynamics should only be calculated every 0.5 macrostep, not at every
sub-step. how will we know if the gas reaching the surface is the same temperature as the surface without doing the heat
calculation? Decision 3: add the test". These amend Part A of `docs/superpowers/specs/2026-11-12-shoulder-melt-loss-design.md`
(the pressure of A3/A-3 and the cost of A7). The known limitation they answer is `docs/model_assumptions.md` 14.1 and
14.9 (branch head c7a699e). Branch `interim-gas-loads`, base for this plan 40ec2f4.

Controller's reading of the decisions (ledger rulings of this plan):
- "Low Mach" is below Mach 3.5, just above the step's onset (Mach 3.39), with a linear blend in Mach to the current
  joined law at Mach 5. An abrupt switch would change the windward pressure in one step by up to 45 % between 60 and 90
  degrees (measured at 30 km: at Mach 3.39 Daniels & Schwartzentruber give 5.77 times the freestream pressure at 60
  degrees against the joined law's 4.02). Their pressure is anchored to the model's unreacted stagnation pressure:
  their stagnation coefficient differs from it by 0.3 to 2 %, and the similarity march takes the nose velocity gradient
  from the 4-degree bin, where the joined law's fall is only about 0.5 %, so an unanchored offset would distort the
  stagnation heat flux.
- "Every 0.5 macrostep" is once per `--dt` (default 0.5 s) of flight time: the provider's solution is held over the
  continuum steps (`--dt-continuum`, default 0.0125 s) in between. Flights whose step is `--dt` throughout (every
  non-melting flight, and melting flights before the continuum switch) are unchanged.
- The skip in Task 3 answers Asha's question: it needs no heat calculation (see Task 3). It runs last so that it can be
  dropped if Asha declines it.

## Global Constraints

- Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`; run tests from `Reentry Model/` with
  `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest ... -p no:cacheprovider` (add `-n 4` for whole modules).
- `--gas-loads classic` must still reproduce main (4e857e1) byte for byte: `analysis/skyscrapers_port/same_results.py
  compare` differing only in `settings.gas_loads` and the absolute `settings.fap_day` path, as in model_assumptions
  14.5 (Task 4).
- Shared physics stays in `reentry_model/shared/` and imports nothing from the rest of `reentry_model`;
  `tests/test_reentry_model_shared.py` must pass.
- At and above Mach 5 (the perfect-gas freestream Mach of `daniels_schwartzentruber.condition(...).mach_inf`) the
  provider's pressure, shear and heat flux must be bit-identical to c7a699e for the same call: the existing fixture tests
  in `tests/test_reentry_model_gas_loads.py` (Holman & Boyd states) pass unchanged.
- Pinned values: `DS_PRESSURE_MACH_LOW = 3.5`, `DS_PRESSURE_MACH_HIGH = 5.0`; hold tolerance `HOLD_TOL = 1e-6` (a
  fraction of the hold); everything pinned in the earlier plan (`2026-11-13-interim-gas-loads.md`, Global
  Constraints) still holds except the joined-pressure line, which now applies at and above Mach 5 only.
- The droplet stage is untouched; `PhysicsHeating` is not changed.
- Commit messages: a scope prefix and a sentence describing the behaviour that changed, ending with the line
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Do not merge or push.

## Review Focus

1. Subsonic states (Mach 0.5 and 0.8) on the new pressure: finite loads, a separation angle found, no exception, the
   stagnation bin at the unreacted stagnation pressure. Test in Task 1.
2. The blend's edges: Mach exactly 3.5 and 5.0 give weights 1 and 0, Mach 4.25 gives 0.5, and the pressure is
   continuous in Mach across both edges. Test in Task 1.
3. A melting flight saved in the middle of a hold window and resumed reproduces the unbroken flight bit for bit. Test in
   Task 2.
4. A state saved before Task 2 (a pickled `GasLoads` without the hold attributes) resumes and solves at every new key as
   before. Test in Task 2.
5. The skip at a state where the wall enthalpy equals the total enthalpy, and at a state with no hold active: zero
   heating, blend 1.0, no solve, exactly as the code before Task 3 reports for that state. Test in Task 3.

---

### Task 1: Daniels & Schwartzentruber's pressure below Mach 5

**Files:**
- Modify: `reentry_model/shared/gas_loads.py` (`joined_pressure` stays as it is; `MeridianLoads`; `GasLoads._solve`;
  module docstring)
- Test: `tests/test_reentry_model_gas_loads.py`
- Docs: `docs/model_assumptions.md` 14.1 and 14.9; the README's gas-loads line that points at the limitation; a dated
  amendment line in the design's A3 pressure paragraph (`docs/superpowers/specs/2026-11-12-shoulder-melt-loss-design.md`)

**Interfaces:**
- Produces: `DS_PRESSURE_MACH_LOW = 3.5`, `DS_PRESSURE_MACH_HIGH = 5.0`;
  `ds_pressure_weight(mach: float) -> float` (1 at and below the low value, 0 at and above the high value, linear in
  Mach between); `anchored_ds_pressure(cp_ds: np.ndarray, p_stag: float, p_inf: float, q_inf: float) -> np.ndarray`
  returning `p_inf + cp_ds * q_inf * (p_stag - p_inf) / (cp_ds[0] * q_inf)` (raises `ValueError` if `cp_ds[0] <= 0`);
  `MeridianLoads.ds_pressure_weight: float = 0.0` as the last field (with a default, so positional construction
  elsewhere is unchanged); `meet_deg` is nan when the weight is 1.
- In `_solve`: `w = ds_pressure_weight(cond.mach_inf)`; `p_w = joined` if `w == 0`, `anchored` if `w == 1`, else
  `w * anchored + (1 - w) * joined`, all before the existing `P_FLOOR` clamp. The zero-speed branch reports weight 0.

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_gas_loads.py`, using the module's `aero_state`
  helper at the 30 km freestream (density 0.0184 kg/m3, 226.5 K), radius 0.05 m, wall 900 K, speed
  `M * math.sqrt(1.4 * 287.05 * 226.5)`:
  - `test_ds_pressure_weight_edges`: weights 1.0 at Mach 0.5, 3.5; 0.5 at 4.25 (abs 1e-12); 0.0 at 5.0 and 8.0.
  - `test_low_mach_pressure_is_anchored_daniels_schwartzentruber`: at Mach 2 and 3, `ds_pressure_weight == 1.0`,
    `p_w` equals `anchored_ds_pressure` of `ds.coefficients` at the same bins (rtol 1e-12), `p_w[0]` equals the
    unreacted stagnation pressure `UnreactedAir().stagnation(rho, T, V).p` (rtol 1e-12), and `meet_deg` is nan.
  - `test_no_step_at_90_degrees_below_mach_4`: at Mach 2 and 3, for both `p_w` and `tau` on the bins from 80 to 94
    degrees, the ratio of each bin to the one before changes by less than 0.05 from one pair of bins to the next (the
    old law fails this: at Mach 2 its pressure ratio jumps from 1.00 to 0.78 at 90 degrees and its shear rises from
    8.95 to 56.5 Pa between 88 and 90 degrees; the new law's ratios change by about 0.01 to 0.02).
  - `test_blend_between_mach_3_5_and_5`: at Mach 4.25, `p_w` equals `0.5 * anchored + 0.5 * joined_pressure(...)[0]`
    (rtol 1e-12); at Mach 6, `p_w` equals `joined_pressure(...)[0]` exactly and the weight is 0.
  - `test_pressure_continuous_across_the_blend_edges`: `p_w` at Mach 3.5 and 3.5 + 1e-9, and at 5.0 and 5.0 - 1e-9,
    agree to rtol 1e-6.
  - `test_subsonic_states` (Review Focus 1): at Mach 0.5 and 0.8, all fields finite, `separation_deg` finite and below
    180, `p_w[0]` the stagnation pressure.
- [ ] **Step 2: Run them**: `"$PY" -m pytest tests/test_reentry_model_gas_loads.py -k "weight or anchored or step_at_90 or blend or continuous or subsonic" -q`.
  Expected: FAIL (`ds_pressure_weight` not defined).
- [ ] **Step 3: Implement** the Interfaces block in `reentry_model/shared/gas_loads.py`, and rewrite the module docstring's
  pressure paragraph to say which law applies at which Mach and why the Daniels & Schwartzentruber pressure is anchored.
- [ ] **Step 4: Run the whole module and the neighbours**: `"$PY" -m pytest tests/test_reentry_model_gas_loads.py tests/test_reentry_model_surface_flow.py tests/test_reentry_model_heating.py tests/test_reentry_model_similarity.py tests/test_reentry_model_shared.py -q -n 4`.
  Expected: all pass, the existing fixture tests unchanged (Global Constraints). If a test that pins a value at a state
  below Mach 5 fails, report its old and new value; do not loosen it without the controller's ruling.
- [ ] **Step 5: Document.** In `docs/model_assumptions.md` 14.1, replace the known-limitation bullet with the new rule
  (Daniels & Schwartzentruber below Mach 3.5, anchored; linear blend to Mach 5; the joined law above), with measured
  numbers from a short scratch script at the 30 km state: the Mach 2 and 3 shear from 80 to 94 degrees before and after
  (Step 1's states), the stagnation heat flux change at Mach 4 and 4.25 (expected about -3 % and -2 %), and the note that
  Daniels & Schwartzentruber validated their model against CFD at Mach 0.37, 2, 3, 4, 7 and 10. In 14.9 remove the
  limitation and add what remains (the blend is a choice, not a fitted law; the trajectories cross Mach 5 to 3.5 in about
  8 s). Update the README line and add the design amendment line (date 2026-11-14, Asha's decision 1).
- [ ] **Step 6: Commit** (`reentry_model: ...` for the code and tests, `docs: ...` for the documents).

### Task 2: One boundary-layer solve per `--dt` of flight

**Files:**
- Modify: `reentry_model/shared/trajectory.py` (`AeroState`, `Simulator.aero_state`)
- Modify: `reentry_model/shared/gas_loads.py` (`GasLoads`)
- Modify: `reentry_model/cli.py` (`build_thermal`: `GasLoads(gamma_pm=args.gamma_pm, hold_s=args.dt)`; a similarity
  run records `settings.gas_loads_hold_s`)
- Test: `tests/test_reentry_model_gas_loads.py`, `tests/test_reentry_model_restart.py` (rewrite
  `test_a_macro_step_of_a_skin_flight_solves_one_meridian_shared_by_the_heating_and_the_surface_flow`)
- Docs: `docs/model_assumptions.md` section 14 (where the provider is described: the solve cadence and what it costs in
  accuracy), 14.8 (cost)

**Interfaces:**
- Consumes: Task 1's `GasLoads._solve` unchanged in signature.
- Produces: `AeroState.t: float = None` as the last field (the flight time in s; `aero_state(t, r, v)` fills it);
  `GasLoads(gamma_pm=..., kn_switch=..., bin_deg=..., hold_s=None)`; class attributes `hold_s = None`, `_held = None`,
  `_held_t = None` (so a `GasLoads` pickled before this task has no hold); `GasLoads.held_loads(state) ->
  MeridianLoads | None`: the held solution when `hold_s` is set, `state.t` is not None, a solution is held and
  `0 <= state.t - _held_t < hold_s * (1 - HOLD_TOL)`, else None; `meridian(state, radius, T_wall)` returns
  `held_loads(state)` when that is not None, otherwise as today (key cache, then solve), and then, when `hold_s` is set
  and `state.t` is not None, holds the result at `state.t`. `__getstate__` keeps `_held` and `_held_t` (it still drops
  the key cache and the Cantera objects).
- While a solution is held the provider ignores changes in the freestream, nose radius and wall temperature; the
  heating's per-patch hot-wall factor and the free-molecular flux still use the current wall temperatures and freestream
  each step, as today.

- [ ] **Step 1: Write the failing tests.**
  - In `test_reentry_model_gas_loads.py`: `test_hold_returns_the_held_solution_within_the_window` (with `hold_s=0.5`:
    calls at t = 0, 0.2 and 0.4999 with different density, radius and wall temperature return the same object and
    `_solve` runs once, counted by monkeypatching `GasLoads._solve`; a call at t = 0.5 solves again);
    `test_no_hold_without_time_or_hold` (`hold_s=None`, or `state.t` None: behaves as the key cache does today, as in
    `test_cache_shared`); `test_pickle_keeps_the_held_solution` (after `pickle.loads(pickle.dumps(g))` within the
    window, `meridian` returns an equal solution without solving); `test_state_pickled_before_the_hold_solves_every_key`
    (Review Focus 4: delete `hold_s`, `_held`, `_held_t` from an instance's `__dict__`, round-trip it through pickle, and
    check that two calls at different keys solve twice).
  - In `test_reentry_model_restart.py`, rewrite the T7-4 test as
    `test_a_skin_flight_solves_one_meridian_per_dt_shared_by_the_heating_and_the_surface_flow`: the same flight with
    `--t-max 0.0625` gives 5 rows (t = 0, 0.025, 0.0375, 0.05, 0.0625) and exactly 3 solves (at t = 0, the heating's
    alone, then at 0.025 and 0.05); every call is at a solved key or within a held window.
  - In `test_reentry_model_restart.py`: `test_a_flight_resumed_inside_a_hold_window_matches_the_unbroken_flight`
    (Review Focus 3): that flight saved with `--save-state-at 0.025` (just after the solve at 0.025, so the resumed step
    to 0.0375 must use the held solution) and resumed to `--t-max 0.0625` gives rows identical to the unbroken flight's
    (use the module's existing resume-comparison helpers).
- [ ] **Step 2: Run them**; expected FAIL (no `hold_s`; 5 solves instead of 3).
- [ ] **Step 3: Implement** the Interfaces block, `HOLD_TOL = 1e-6` in `gas_loads.py`, and the `cli.py` changes. Check
  the `same_results.py` smoke cases and every other caller of `AeroState(...)` still construct it (the new field has a
  default).
- [ ] **Step 4: Run** `tests/test_reentry_model_gas_loads.py`, `tests/test_reentry_model_restart.py`,
  `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_trajectory.py` and `tests/test_reentry_model_shared.py`
  with `-n 4`. Expected: all pass. Tests of melting flights whose rows change after a continuum switch are expected to
  change only where they pin provider-dependent values after the switch; report each with old and new values.
- [ ] **Step 5: Document** the cadence in section 14: the solution is held for one `--dt` (0.5 s by default) and what
  that costs: on the Step 2 trajectories the stagnation heat flux changes between consecutive 0.5 s steps by at most
  7.0 % (100 mm) and 6.6 % (50 mm) while it is above a tenth of its peak, with medians of 0.5 % and 1.5 % (measured
  from the non-melting flights' histories), so a held solution lags by at most that much at the end of a window. Give
  the solve count per 0.5 s after the continuum switch (40 before, 1 now) in 14.8.
- [ ] **Step 6: Commit.**

### Task 3: Skip the solve when no part of the wall can be heated

**Files:**
- Modify: `reentry_model/shared/gas_loads.py` (`GasLoads.total_enthalpy`)
- Modify: `reentry_model/shared/heating.py` (`SimilarityHeating.evaluate`, class docstring)
- Test: `tests/test_reentry_model_heating.py`
- Docs: `docs/model_assumptions.md` 14.8 (cost) and the section on the heating

**Interfaces:**
- Consumes: Task 2's `GasLoads.held_loads(state)`.
- Produces: `GasLoads.total_enthalpy(state) -> float`: the unreacted stagnation enthalpy `air.stagnation(rho, T_inf,
  V).h` exactly as `_solve` computes `H_e` (no boundary-layer solve).
- `SimilarityHeating.evaluate`, after the zero-speed check: `L = gas_loads.held_loads(state)`; if `L` is None, compute
  `H_e = gas_loads.total_enthalpy(state)` and `h_w = float(gas_loads.wall_enthalpy(T_stag))` at the stagnation patch's
  wall temperature; if `H_e - h_w <= 0` return the zero result `HeatingResult(zeros, 0.0, nan, 0.0, 1.0)` without
  calling `meridian`; otherwise call `meridian` and continue exactly as today.

Why this needs no heat calculation (for the docstring and the docs): the total enthalpy of the gas reaching the body,
the freestream enthalpy plus half the speed squared, does not change across the shock, so it is known from the
trajectory alone; the enthalpy of gas at the wall's temperature is a table lookup. The boundary-layer solve only says
how much of that difference arrives and where. At the stagnation point the edge speed is zero, so a wall at or above
the total enthalpy receives no heat; elsewhere the gas recovers less than the total enthalpy (recovery factor about
the square root of the Prandtl number), so it receives none either. The model already returns zero heating in exactly
this case (`heating.py:194`), after paying for the solve; the skip makes the same decision before it. It neglects, as
the model already does, the convective cooling of a wall hotter than the gas.

- [ ] **Step 1: Write the failing tests** in `tests/test_reentry_model_heating.py`:
  `test_similarity_heating_skips_the_solve_when_the_wall_is_hotter_than_the_total_enthalpy` (a slow state, for example
  Mach 2 at 30 km, with a 1500 K stagnation wall: zero fluxes, blend 1.0, and `GasLoads._solve` never called, counted by
  monkeypatch); `test_skip_matches_the_solved_result` (for two states, one each side of `H_e = h_w`, the result's
  `q_conv`, `q_stag` and `blend` equal those of the same call with the skip disabled by monkeypatching
  `GasLoads.total_enthalpy` to return infinity, which sends it down the solving path);
  `test_skip_when_the_wall_enthalpy_equals_the_total_enthalpy` (Review Focus 5: `GasLoads.wall_enthalpy` monkeypatched
  to return exactly `total_enthalpy(state)`: zero heating, blend 1.0, no solve);
  `test_skip_is_not_taken_inside_a_hold_window` (within a held window the held solution is used and `total_enthalpy` is
  not called).
- [ ] **Step 2: Run them**; expected FAIL (`_solve` called).
- [ ] **Step 3: Implement** the Interfaces block.
- [ ] **Step 4: Run** `tests/test_reentry_model_heating.py`, `tests/test_reentry_model_gas_loads.py`,
  `tests/test_reentry_model_restart.py` and `tests/test_reentry_model_shared.py` with `-n 4`. Expected: all pass.
- [ ] **Step 5: Document** the skip and its reasoning (the paragraph above, in plain words) in the heating part of
  section 14; in 14.8 give the share of steps it skips on the Step 2 non-melting flights, measured from their histories
  as the rows with zero `q_stag_Wm2` and `heating_blend_f` 1.0: 398 of 734 steps from 168 s on the 100 mm flight and 553
  of 1119 from 283 s on the 50 mm flight. Note that these walls reach 1878 K because the reference flights do not melt;
  a melting wall stays near 900 K, so the skip starts much later there (only below about 1.2 km/s).
- [ ] **Step 6: Commit.**

### Task 4: The Step 2 reference flights on the default loads, and verification

**Files:**
- Modify: `tests/test_reentry_model_reference_thermal.py` (`coupled_run` gains `heating_name == "similarity"`;
  new test)
- Docs: `docs/model_assumptions.md` 14.7 and 14.8, `Reentry Model/CLAUDE.md` (the reference tier's time and the unit
  tier's count), the root `README.md` gas-loads paragraph

**Interfaces:**
- Consumes: Tasks 1-3 (`GasLoads(hold_s=...)`, `SimilarityHeating`).
- `coupled_run(ref, "similarity", ...)` builds `heating.SimilarityHeating(GasLoads(hold_s=dt))`.

- [ ] **Step 1: Add** `test_similarity_mode_is_reported` (marker `reference`, parametrized over `NAMES`): the flight
  ends on the ground, `abs(energy_balance_residual) < 1e-6`, `0.4 < ratio_end < 1.1`, and the metrics are recorded as
  `<key>__similarity.json` by `record`. Its docstring says the ratio is reported, not thresholded, and quotes the values
  measured in Step 2.
- [ ] **Step 2: Run** `"$PY" -m pytest tests/test_reentry_model_reference_thermal.py -m reference -k similarity -q`
  with `/usr/bin/time -p`; expected PASS. Record each flight's `ratio_end` and wall time (from `results.runtime_s`);
  before this plan they were 0.862 and 0.882, in 1068 s and 1684 s.
- [ ] **Step 3: Verify classic** with `analysis/skyscrapers_port/same_results.py` (run with `--append "--gas-loads
  classic"` on the cases `step1,skins,element,nomelt`, then `compare` against the main outputs already in
  `$SCRATCH/task8/sr_main`, where `$SCRATCH` is the session scratchpad the controller gives you): the only differences
  are `settings.gas_loads` and `settings.fap_day`, as in 14.5.
- [ ] **Step 4: Run the unit tier** (`"$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4`)
  and the droplet tests (`"$PY" -m pytest reentry_model/droplet_reentry/tests -q`); expected all pass.
- [ ] **Step 5: Document** the measured ratios and times in 14.7 and 14.8 (old and new, each measured once, and the
  machine's load), the reference tier's new time and the unit tier's new count in CLAUDE.md, and one line in the
  README.
- [ ] **Step 6: Commit.**
