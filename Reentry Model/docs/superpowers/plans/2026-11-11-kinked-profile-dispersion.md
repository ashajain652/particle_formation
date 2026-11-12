# Girin's gradient instability on the skins' own melt profile — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** On spraying skin runs, the model solves Girin's gradient mode on the piecewise-linear profile of the skins' own cell speeds, by default. A floor test places a wall at the shallowest floor-like kink, and the faster of that mode and Girin & Kopyt's gas-driven mode releases liquid and slurry together, down to a threshold temperature. `--dispersion table` switches back to Girin's theory as it stands today, bit for bit.

**Architecture:** All of the spraying's physics goes in the shared package `reentry_model/shared/`, written on plain arrays, so the level-set architecture can call it later (Asha, 2026-11-11, decision 8):
- `shared/dispersion.py` gains the broken-line solver and the floor test.
- A new `shared/melt_profile.py` holds the straight-segment fit, the layer's properties, the release-run rule, the mode inputs and the history summary.
- `shared/spray.py` gains the mode choice and the record layout.
- The setting lives on the shared spray model, as `spray_gates` and `slip_band` do.

The skin modules (`skin.py`, `skinbody.py`, `cli.py`) only read their cells into plain arrays, take cells off the top, compute their own shear rates and wire the option. The new path is one branch inside existing functions, keyed on an optional argument that defaults to today's behaviour, so `table` stays identical.

**Tech Stack:** Python 3, numpy, scipy (as already used), pytest with pytest-xdist. Environment: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`.

**Spec:** `Reentry Model/docs/superpowers/specs/2026-11-11-kinked-profile-dispersion-design.md`, with decisions 8 (shared placement) and 9 (default `profile`) added on 2026-11-11. Read it before any task. Section numbers below ("spec 6.4") are the spec's.

## Global Constraints

- Paths are relative to `Reentry Model/` unless they start with the repository top. The work goes on branch `kinked-dispersion`, in worktree `<repository top>/.claude/worktrees/kinked-dispersion`, cut from `main`. Commits stay local; nothing is pushed or merged without Asha.
- **Module paths.** The skyscrapers subpackage (`docs/superpowers/plans/2026-11-11-skyscrapers-subpackage.md`) landed on `main` on 2026-11-12 (655d1a2): `skin.py`, `skinbody.py`, `restart.py`, `coupled.py` and the other skyscrapers-only modules live in `reentry_model/skyscrapers/`, with alias files at the old paths, and `shared/` gained `lubrication`. The paths in this plan name the moved files. Never edit an alias file.
- Commit messages follow the repository's style: one or two full sentences saying what changed and why, then the line `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Shared placement (decision 8).** Every physics step lives in `reentry_model/shared/` and takes plain arrays:
  - the solver;
  - the floor test;
  - the fit;
  - the layer's properties above the wall;
  - the mode inputs and the mode choice;
  - the release-run rule and the choice of threshold temperature;
  - the record layout and the history summary.

  The skin modules only read cells into `melt_profile.LayerColumns`, take cells off the top (the cascade), compute the skins' own shear rates (`skin.cell_shear_rates`, which are the skins' flow numerics, not spray physics) and carry the option.
- The shared package imports nothing from the rest of `reentry_model` except the package-root constants. The new module `melt_profile` is added to `SHARED_MODULES` in `tests/test_reentry_model_shared.py`, which raises its file count from 23 to 24, and to the shared-module list in `CLAUDE.md`.
- **The default (decision 9).** On a spraying skin run (`--surface-model skin --melt on --skin-flow on`), the command line's default is `--dispersion profile`. `--dispersion table` is the switch back.
  - The setting is `SprayModel.dispersion`. Its class attribute and its constructor default are both `"table"`, so saved states and objects built directly in tests keep Girin's theory. The command line passes `"profile"` explicitly.
  - The element model's spray model is always `"table"`.
  - `--dispersion` given explicitly on any run that is not a spraying skin run is refused with exit code 2.
- With `--dispersion table`, a run is bit-identical to today's: the same history columns and values, the same 27 record columns (`spray.SOURCE_COLUMNS` unchanged) and the same JSON, except that a spraying skin run's JSON always records `dispersion`. The run name marks the switch back with the suffix `_dispersion-table`; `profile` is unmarked.
- These constants are labelled "assumed until stage 6" wherever they are defined: `KINK_DROP = 0.2`, `FLOOR_RE = 1.0`, `FIT_TOL = 0.03` (as a fraction of the top cell's speed) and `N_MAX = 6`.
- Unchanged: `K_R = 0.17`, `K_T = 1.1`, `WE_CRITICAL_THEORY = 3.08` (table path only), and the surface tension of 0.80 N/m (the liquid's, and `skinbody.SLURRY_SIGMA`). The solver always takes the surface tension as an argument, with no default.
- The plan is evaluated once per macro step, at its start, on windward wet columns only. The added cost must stay under 5 % of a macro step's wall time on the d050 state.
- Out of scope: the element model (`body.MeltingBody`) has no resolved profile and stays on the table, the level-set architecture is not wired, and the Rayleigh–Taylor branch and the rarefied branch's physics are unchanged.
- Unit tier: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4` (about 24 minutes). A single module: `"$PY" -m pytest tests/<file> -q`.
- Reports to Asha follow `~/.claude/CLAUDE.md`: full sentences, every number explained and terms spelled out. The order is what was asked, what was found, the evidence, what it means, and her decisions.

## Review Focus

These are the inputs most likely to break the profile path that the spec implies but does not test. Each one has a test in the task named.

1. **Plug columns.** The speed peaks below the surface, so the top segment's slope is negative, and the temperature need not fall with depth. The fit and the kink rule use slope magnitudes. If the top cell is colder than the threshold, the release run is empty: nothing leaves, nothing fails, and the books close. Tests are in Tasks 4, 6 and 7.
2. **Edge columns.** A column with a single flowing cell gets its surface speed from the line through that cell and the zero at the front. A kink at the last flowing cell has a rigid cell below it, so its viscosity is infinite and Re_f is 0. Both give finite outputs and raise no numpy warnings. Tests are in Tasks 2 and 4.
3. **Stable or straight layers, and empty steps.** A column whose layer above the wall is one straight segment, or whose waves are all stable, gets growth 0, a release depth of NaN and NaN record columns, with no divide warnings. A macro step with no resolved column writes NaN medians and a share of NaN. Tests are in Task 7.
4. **Resuming across a change of setting.** Switching `profile` to `table` cuts in-memory records to 27 columns, and switching the other way pads them with NaN. An old pickled `SourceSpool` with no `names` attribute reads as 27 columns. A state saved before the option existed resumes on `table`. Tests are in Tasks 5 and 8.
5. **Degenerate pencils.** Zero surface tension and infinite eigenvalues, a wall exactly at a breakpoint's depth, and very thin segments at large wavenumber. The solver returns finite eigenvalues only and raises no `LinAlgError`. Tests are in Task 1.

---

### Task 0: Worktree, branch and baseline

**Files:**
- Commit: this plan and the amended spec (decisions 8 and 9), on the branch.

- [ ] **Step 1:** Use superpowers:using-git-worktrees to create the worktree `.claude/worktrees/kinked-dispersion` on a new branch `kinked-dispersion` from `main`.
- [ ] **Step 2:** Copy in this plan and the amended spec if they are not there yet, then commit both.
- [ ] **Step 3: Record the baseline.** Run this before any code change, on mains power, with nothing else editing:
  `"$PY" analysis/skyscrapers_port/same_results.py run --code "$PWD" --out output/kinked_check/before --cases skins,element,resume --state "<repository top>/.claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_snapshots/d100__physics__skin.state.pkl"`
  Expect exit 0 in each case's log.
- [ ] **Step 4:** Run the unit tier. Expect 702 passed and 1 skipped, and record the count.

### Task 1: The broken-line solver (spec 6.2, stage 1; shared)

**Files:**
- Modify: `reentry_model/shared/dispersion.py`. Append the new code; today's table code stays untouched.
- Test: `tests/test_reentry_model_dispersion.py`

**Interfaces:**
- Produces, in `reentry_model.shared.dispersion`:
  - `DISPERSION_NAMES = ("table", "profile")` and `DEFAULT_DISPERSION = "profile"` (the command line's default on spraying skin runs, decision 9)
  - `@dataclass class KinkedProfiles` holds `count (n,) int`, `depth (n, B) m`, `u_top (n,) m/s`, `u_break (n, B) m/s`, `wall (n,) m` (inf means no wall), `u_wall (n,) m/s` and `rho (n, B+1) kg/m3`. Arrays are padded with NaN beyond each column's `count`. Segment j, with j at most `count`, runs from breakpoint j−1 to breakpoint j (the surface to breakpoint 1 for j = 0, and the last breakpoint to the wall for j = `count`).
    - `rows(idx) -> KinkedProfiles`
    - `slopes() -> (n, B+1)` gives dU/dy per segment with y pointing up, so a speed falling with depth has a positive slope.
    - The classmethod `one(depths, speeds, wall=np.inf, rho=1.0)` builds one profile, with `speeds = [u_top, *u_break, u_wall]`.
    - Validation raises ValueError if the depths do not increase, if `wall <= depths[-1]`, if `wall = inf` while `u_wall != u_break[-1]` (an unbounded bottom must be uniform), or if `count = 0` with `wall = inf`.
  - `broken_line_eigenvalues(profiles, k, sigma) -> complex (n, m, B+2)` gives every finite eigenvalue c [m/s], padded with NaN. `k` is in 1/m, (m,) or (n, m); `sigma` is in N/m, scalar or (n,).
  - `broken_line_c(profiles, k, sigma) -> complex (n, m)` gives the eigenvalue with the largest imaginary part.
  - `@dataclass class Waves` holds `k (n,) 1/m` (NaN where stable), `growth (n,) 1/s` (Im omega_f, 0 where stable), `c (n,) complex` and `release_depth (n,) m` (NaN where stable).
  - `fastest_waves(profiles, sigma) -> Waves`
  - `wave_structure(profiles, k, c, y) -> (phi complex (n, len(y)), displacement (n, len(y)))`. phi is normalised so phi(0) = 1, and displacement is |phi/(U − c)| over its surface value. y is depth in m.
  - `UNSTABLE_REL = 1e-6`: an eigenvalue is unstable where Im c exceeds this times the profile's largest |U|.
  - `K_SPAN = (0.02, 5.0)`, `K_PER_DECADE = 24`, `REFINE_ROUNDS = 3` and `REFINE_POINTS = 9`.

The physics is spec section 3. The method is the builder's (spec 8), within these requirements:
- Each segment's exponentials use a basis normalised inside that segment.
- Infinite eigenvalues are discarded.
- The search grid is log-spaced at `K_PER_DECADE` from `K_SPAN[0] / D` to `K_SPAN[1] / s_min`. D is the wall's depth, or the deepest breakpoint's without a wall; s_min is the thinnest segment.
- Refinement takes `REFINE_ROUNDS` rounds of `REFINE_POINTS` log-spaced points between the best point's two neighbours, batched over columns.
- One route that batches well: eliminate the c-free conditions (continuity of phi, and the wall), which leaves an (N+2)-square pencil. Solve it with `np.linalg.eig` where its B-matrix is well conditioned, and with `scipy.linalg.eig(A, B)` per column otherwise.
- The release depth is the first depth below the displacement's peak where the displacement falls to 1/e of its surface value. It is sampled finely in each segment and to 10/k below an unbounded bottom, interpolated linearly, and capped at the wall.

Test helpers, in the test file:

```python
def girin(we):            # Girin's units: delta_m = V_s = rho = 1, sigma = 1 / We_s
    return dp.KinkedProfiles.one([1.0], [1.0, 0.0, 0.0]), 1.0 / we

def slurry_over_wall(beta, eps, r=1.0, we=1e4):   # notes fulltext_girin_papers.md, lines 79-106
    S1 = 1.0 / (1.0 + eps * (beta - 1.0)); u_h = eps * (beta - 1.0) * S1
    return dp.KinkedProfiles.one([1.0], [1.0, u_h, 0.0], wall=beta, rho=[1.0, r]), 1.0 / we, S1, u_h

def ghosh(K2, d=0.4):     # Ghosh 2017 in units of V and H: kink at depth d, top slope K2 times the bottom's, sigma 0
    S1 = 1.0 / (d + (1.0 - d) / K2)
    return dp.KinkedProfiles.one([d], [1.0, (1.0 - d) * S1 / K2, 0.0], wall=1.0), 0.0
```

- [ ] **Step 1: Write the failing tests.**
  - `test_girin_profile_gives_the_roots_of_his_cubic`: for (Δ, We) in ((0.5, 5), (1.2, 20), (1.225, 1e4), (3.0, 50)), the finite eigenvalues times Δ, sorted, equal `np.roots([1, *dp.cubic_coefficients(Δ, We)])` sorted, within 1e-10 absolute. There are exactly 3.
  - `test_girin_profile_fastest_wave_matches_the_cached_table`: at We = 1e4, `abs(w.k - Δ_table) <= 0.002` (the table's grid spacing) and `Im_table - 1e-9 <= w.growth <= Im_table + 1e-4`, with `Δ_table, Im_table = dp.fastest_mode(1e4)[:2]`.
  - `test_girin_threshold_is_three`: at We = 3.00, `broken_line_c(...).imag <= 1e-10` for k in `np.logspace(-4, 0.7, 400)`. At We = 3.01, `broken_line_c(..., k=0.0033).imag > 1e-4`. Measured on the cubic, 2026-11-11: at 3.01 the unstable band is Δ below 0.0044, with Im Omega at most 3.2e-6. That lies below the fastest-wave grid, so this test calls the eigenvalue function directly.
  - `test_slurry_over_a_wall_gives_the_closed_form_roots`: for (β, ε, r) in ((1.5, 0, 1), (2, 1e-3, 1), (5, 1e-3, 1.1), (20, 1e-3, 0.95)) at Δ in (0.5, 1.2, 3.0), every finite eigenvalue c makes the closed form, multiplied by (1 − c)(u_h − c), vanish to 1e-10 relative to its largest term. The closed form is `(1 − c)(Δ tanh Δ + q) − [S1 + Δ²/(We (1 − c))](1 + (q/Δ) tanh Δ)` with `q = r Δ coth(Δ(β − 1)) + S1 (1 − r ε)/(u_h − c)`.
  - `test_slurry_over_a_wall_fastest_waves`: at We = 1e4, (k, growth) is (1.21, 0.196) at β = 1.5, (1.24, 0.234) at β = 2 and (1.225, 0.247) at β = 20. The tolerance is ±0.01 on k and ±0.002 on growth. These are the side session's values.
  - `test_ghosh_limits`: `ghosh(1.5)` peaks at growth 0.229 ± 0.002 at k 2.6 ± 0.15. `ghosh(K2)` for K2 = 1.1, 2.5 and 3.0 peaks at 0.109, 0.349 and 0.383 (±0.003). These come from the side session's calculation, which matched Ghosh's determinant; his Fig. 5a, read from the image, peaks just below 0.25. `ghosh(0.5)` and the straight film `KinkedProfiles.one([], [1.0, 0.0], wall=1.0)` have Im c at most 1e-10 for k in `np.logspace(-2, 2, 200)` and σ in (0, 0.01, 1).
  - `test_release_depth_on_girin_profile`: at We = 1e4 and the fastest wave, the displacement peaks at 1.6 ± 0.05 times its surface value, 0.67 ± 0.05 down. `w.release_depth` is 1.82 ± 0.02, and |phi| falls to 1/e at 1.29 ± 0.02.
  - `test_a_wall_closing_on_the_kink_kills_the_growth`: `slurry_over_wall(1.001, 1e-3)` grows at less than 1e-2 of `slurry_over_wall(20, 1e-3)`.
  - `test_deep_walls_and_short_waves_stay_accurate`: Girin's profile with a wall at H = 100/k, for k in (1, 10, 100), equals the unbounded result to 1e-8 relative. Segments 1e-3 and 10 thick at kH = 1000 give finite eigenvalues under `np.errstate(all="raise")`.
  - `test_degenerate_pencils_give_finite_eigenvalues_only` (Review Focus 5): with σ = 0, a wall 1e-12 below a breakpoint, and breakpoints 1e-9 apart, the output has no inf; NaN appears only as padding.
- [ ] **Step 2:** Run `"$PY" -m pytest tests/test_reentry_model_dispersion.py -q`. Expect FAIL with AttributeError.
- [ ] **Step 3:** Implement. Add a docstring giving the conditions, the sources and the note "derivation ours, checked against Girin Eq. 1 and Ghosh 2017".
- [ ] **Step 4:** Run the module and `tests/test_reentry_model_shared.py`. Expect PASS.
- [ ] **Step 5:** Commit.

### Task 2: The floor test and the wall (spec 6.3, stage 2; shared)

**Files:** Modify `reentry_model/shared/dispersion.py`. Test in `tests/test_reentry_model_dispersion.py`.

**Interfaces:**
- Consumes: `KinkedProfiles` (Task 1).
- Produces:
  - `KINK_DROP = 0.2`, `FLOOR_RE = 1.0` and `FASTEST_KD = 1.225`. The first two are labelled assumed until stage 6.
  - `kinks(profiles, kink_drop=KINK_DROP) -> bool (n, B)`. Breakpoint j is a kink where `|S_below| <= kink_drop * |S_above|` and `S_above != 0`.
  - `floor_reynolds(profiles, nu_below) -> (n, B)`: Re_f = |u_top − u_break_j| d_j / (FASTEST_KD ν_below_j), NaN off kinks and in padding. ν_below is in m²/s, (n, B); inf gives 0.
  - `@dataclass class Wall` holds `floor (n,) bool`, `index (n,) int` (−1 for the front), `depth (n,) m`, `speed (n,) m/s` and `re_f (n,)`. `re_f` is taken at the wall for floor-like columns, and at the first kink for fluid-like ones (NaN with no kink).
  - `place_wall(profiles, nu_below, floor_re=FLOOR_RE, kink_drop=KINK_DROP) -> Wall`. Kinks are tested from the top down, and the first with Re_f < `floor_re` is the wall; otherwise the wall is the profile's own (the rigid front).
  - `above_wall(profiles, wall) -> KinkedProfiles`. Floor-like columns get `count = index`, `wall = depth` and `u_wall = speed`; fluid-like ones are unchanged.

- [ ] **Step 1: Write the failing tests.**
  - `test_floor_reynolds_reproduces_the_crossover_speeds`: take ν_l = 1.3e-3/2400, a floor at rest and d in (5, 50, 150, 450) µm. The speed u_top at which Re_f = 1 is (66, 6.6, 2.2, 0.74) m/s at ν_below = 500 ν_l, and (530, 53, 18, 5.9) m/s at 4000 ν_l, each within 2 % (spec 3).
  - `test_a_five_fold_drop_is_a_kink_and_a_four_fold_is_not`: a slope ratio of 0.2 is a kink and 0.25 is not. A plug profile (negative top slope, a five-fold drop in magnitude below it) is a kink.
  - `test_the_wall_is_the_shallowest_floor_like_kink`: take two kinks with Re_f of (3, 0.5); the wall is the deep one. With both below 1 it is the shallow one. With neither, it is the front, `index == -1`, and `re_f` is the first kink's. The velocity jump is measured against the moving floor.
  - `test_above_wall_cuts_the_profile_at_the_wall` (Review Focus 2): the cut gives count = `index`, with the wall at the kink's depth and speed. A kink with ν_below = inf gives Re_f = 0 and is floor-like, under `np.errstate(all="raise")`.
- [ ] **Step 2:** Run the module. Expect FAIL. Implement. Run the module and `tests/test_reentry_model_shared.py`. Expect PASS.
- [ ] **Step 3:** Commit.

### Task 3: The gas-shear check (Asha's decision 7; gates Task 4)

This is the carry-over item "MELT SURFACE SPEED vs GAS WALL SHEAR" (`memory/skyscrapers-port-to-main.md`; the spec's header). First check whether another session has completed it (the memory index, `research_notes/`, recent commits on main). If it has, record its result and go on to Task 4.

**Files:**
- Modify: `reentry_model/skyscrapers/skin.py` (a new function after `interface_drag`)
- Create: `analysis/kinked_dispersion/__init__.py` and `analysis/kinked_dispersion/shear_check.py`
- Test: `tests/test_reentry_model_skin.py`

**Interfaces:**
- Produces: `skin.cell_shear_rates(sk, tau) -> (n, K) 1/s`. Each flowing cell's shear rate comes from the present velocities by the rule `_momentum_columns.rates` uses: the mean over its two faces of the stress there, over its viscosity.
  - The top face carries |tau|.
  - An interior face carries `series_drag` times the speed difference between the centres.
  - The bottom face carries the base drag over the bottom cell's half-thickness times its speed.
  - The viscosities come from a fixed point from rest, in at most `SHEAR_PASSES` passes, stopping at `SHEAR_TOL`.
  - The rate is 0 in cells that do not flow. `tau` is in Pa, (n,).

  `_momentum_columns` itself is not edited. This function is the skins' own flow numerics, so it stays in `skin.py`.

- [ ] **Step 1: Write the failing tests.**
  - `test_cell_shear_rates_on_a_couette_column`: an all-`LIQUID` column holds the steady Couette profile for τ = 50 Pa. Every rate equals τ/μ_l within 1e-10 relative.
  - `test_a_column_of_liquid_over_slurry_settles_to_the_speed_its_viscosities_give`: take `stacked(col, [[LIQUID] * 6 + [SLURRY] * 10 + [RIGID] * 4])` with 50 Pa of shear along x and no push. Step `momentum_step` until no speed changes by more than 1e-10 relative. The top cell's speed then equals τ (A_top/A) times the sum of dz_0/(2η_0) and dz_i/η_i over the flowing cells i ≥ 1, with η at `cell_shear_rates`, within 1e-6.
- [ ] **Step 2:** Run. Expect FAIL. Implement. Run again and expect PASS.
- [ ] **Step 3: The probe.** `shear_check.py --state <d050>` loads the state (`restart.load`) and selects the spec's 622 columns: angle below 60 degrees, and sheared depth over flowing depth between 0.12 and 0.35. For each column it predicts the steady top-cell speed from the face balance of `momentum_step`. The stress at a face is the gas shear times A_top/A, plus the push per unit volume along the top cell's direction times the face's depth. The push comes from `skin.push` with `body._forcing[1]` (falling back to `body.last_flow.G`), and the viscosities from `cell_shear_rates`. It prints the 10th, 50th and 90th percentiles of the actual over the predicted top speed, with the gas shear alone and with the push.
  Run: `"$PY" -m analysis.kinked_dispersion.shear_check --state "<repository top>/.claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/d050__handdown-off.state.pkl"`.
- [ ] **Step 4: The gate.** If the median ratio with the push lies within [0.8, 1.25], go on to Task 4. Otherwise **stop** and report the percentiles and the decomposition to Asha (spec decision 7). The fix is a separate decision.
- [ ] **Step 5:** Commit. Record the percentiles in the message.

### Task 4: The straight-segment fit, shared, and reading the cells (spec 6.4, stage 3)

**Files:**
- Create: `reentry_model/shared/melt_profile.py`. Add it to `SHARED_MODULES` in `tests/test_reentry_model_shared.py`, change the file count from 23 to 24, and add it to CLAUDE.md's shared-module list.
- Modify: `reentry_model/skyscrapers/skinbody.py` (a new `layer_columns`, after `sheared_depth`)
- Create: `analysis/kinked_dispersion/extract_fixture.py` and `data/reference_values/d050_layer_columns.npz`
- Test: `tests/test_reentry_model_melt_profile.py` (new; plain arrays only, no skins) and `tests/test_reentry_model_skinbody.py` (`layer_columns`)

**Interfaces:**
- Consumes: `dispersion.KinkedProfiles` (Task 1) and `skin.cell_shear_rates` (Task 3).
- Produces, in `reentry_model.shared.melt_profile`:
  - `FIT_TOL = 0.03` and `N_MAX = 6`, labelled assumed until stage 6.
  - `@dataclass class LayerColumns` describes n columns of layered melt, top first, all (n, K):
    - `z` (cell-centre depth, m), `dz` (m), `u` (each cell's speed along the top cell's direction, m/s), `T` (K), `m` (kg), `rho` (kg/m3);
    - `mu` (Pa s at the cell's shear rate, inf where it does not flow);
    - `flowing` (bool, the contiguous flowing cells from the top), `liquid` (bool, the contiguous fully liquid cells from the top);
    - `front` (n,): the depth of the rigid front, in m.

    It has `rows(idx)`. Its constructor accepts `**np.load(...)`, so the fixture loads directly.
  - `reducible(cols) -> (n,) bool`: the top cell flows and moves, and the front lies below its centre.
  - `@dataclass class Fit` holds `profiles` (a `KinkedProfiles` with B = `N_MAX`, its wall the front at speed 0) and `cell (n, N_MAX) int` (the cell at each breakpoint, −1 beyond).
  - `fit_profiles(cols, tol=FIT_TOL, n_max=N_MAX) -> Fit` (on the reducible rows; others get count 0 and NaN).
  - `nu_below(cols, cell) -> (n, N_MAX) m²/s`: the kinematic viscosity μ/ρ of the cell just below each breakpoint, inf where that cell is rigid, NaN in padding.
  - `couette_viscosity(cols, depth) -> (nu_eff (n,), rho (n,))`, over the cells whose centres lie above `depth`. μ_eff is their Σdz over Σ(dz/μ), ρ is their mass-weighted density, and ν_eff = μ_eff/ρ.
- Produces, in `skinbody`: `layer_columns(sk, tau) -> melt_profile.LayerColumns` for every column. It projects each cell's velocity on the top cell's direction. Its viscosities are `skin.cell_viscosity(sk, skin.cell_shear_rates(sk, tau))`, its densities `sk.density(sk.T)`, and its front `depth_above(T_rigid)`.

The fit (spec 6.4):
- The surface speed is extrapolated linearly from the top two cell centres; with one flowing cell, it comes from the line through that centre and zero at the front.
- The fit is continuous and piecewise linear, through the surface speed, each breakpoint cell's speed, and zero at the front.
- Breakpoints are added greedily at the cell farthest from the fit until every flowing cell lies within `FIT_TOL` times |the top cell's speed| or `N_MAX` is reached. A pruning pass then drops any breakpoint whose removal keeps the tolerance. **[plan call: greedy insertion plus pruning stands in for "the fewest".]**
- A segment's density is the mass-weighted density of the cells whose centres lie in its closed depth interval.

- [ ] **Step 1: Write the failing shared tests** (synthetic arrays: K = 40 cells with dz = 10 µm and z = (i + ½) dz, 30 of them flowing, and the front at 0.3 mm):
  - `test_girin_linear_profile_gives_one_breakpoint_at_its_depth`: u = 3 clip(1 − z/D, 0) with D = 0.2 mm gives count 1, a breakpoint within one cell of D, and `u_top == pytest.approx(3.0, rel=1e-9)`.
  - `test_a_layer_moving_as_one_piece_has_no_kink`: uniform u = 0.5 gives one breakpoint, at the last flowing cell's centre, and `dispersion.kinks` is False there. **[plan call: the spec's "gives none" is read as "no kink"; the zero at the front forces the one breakpoint.]**
  - `test_a_plug_column_and_a_one_cell_column_fit_without_warnings` (Review Focus 1 and 2), under `np.errstate(all="raise")`.
  - `test_segment_densities_are_mass_weighted` and `test_couette_viscosity_of_two_layers` (exact values for a two-layer column).
  - `test_nu_below_reads_the_cell_under_each_breakpoint`: a rigid cell gives inf.
  - `test_every_fit_holds_its_tolerance_or_its_cap_on_the_d050_fixture`: load `LayerColumns(**np.load(fixture))`. Every fit is within tolerance or has `N_MAX` breakpoints, and the depths increase. The test prints, as a measurement (spec stage 3), two medians against the depth of the 906.4 K isotherm: the first kink's depth, and the wall's from `dispersion.place_wall`. The builder check measured 0.94 for the first kink.
- [ ] **Step 2: Write the failing skin test** `test_layer_columns_reads_the_cells_as_the_skins_see_them`. Velocities rotated off the x axis project onto the top cell's direction. `mu` equals `cell_viscosity(sk, cell_shear_rates(sk, tau))`. The front is `depth_above(T_rigid)`, and `liquid` is the contiguous liquid run.
- [ ] **Step 3:** Run. Expect FAIL. Implement `melt_profile` (it imports only shared modules and numpy) and `layer_columns`. Run the new module, the skinbody module and `tests/test_reentry_model_shared.py`. Expect PASS.
- [ ] **Step 4: The fixture.** `extract_fixture.py --state <d050> --out data/reference_values/d050_layer_columns.npz` saves `layer_columns` for 36 columns: 24 from the spec's thick set, 6 of the 14 plug columns and 6 other wet windward columns, plus their `theta` and `tau`. It should be about 0.2 MB. Run it, re-run the fixture test, and commit. The message gives the two measured medians.

### Task 5: The mode choice, the setting and the record layout (spec 6.1, 6.5, 6.6, 6.8; shared)

**Files:** Modify `reentry_model/shared/spray.py`. Test in `tests/test_reentry_model_spray.py`.

**Interfaces:**
- Produces, in `reentry_model.shared.spray`:
  - `SprayModel.dispersion = "table"` is a class attribute, and the constructor gains `dispersion="table"`, validated against `dispersion.DISPERSION_NAMES`. Callers read it; `evaluate` does not.
  - `BRANCH_PROFILE = 5`.
  - `@dataclass class ProfileModes` holds, all (n,):
    - `resolved` (bool);
    - `growth` (1/s, 0 where stable or straight);
    - `k` (1/m);
    - `wall_depth` (m);
    - `nu_eff` (m²/s);
    - `we_s` (the recorded Weber number, spec 6.7);
    - `cap_gradient` and `cap_driven` (m): the run's depth at the step's start for each mode.
  - `viscous_film(lam, tau, h, liquid, nu=None)`: with `nu` given, it replaces `liquid.mu / liquid.rho`; without it the result is bit-identical.
  - `SprayResult.cap_depth = None` is a new field: the crest cap's depth per column, set on the profile path.
  - `SprayModel.evaluate(..., profile=None, rt_extent=None)`. With `profile=None`, nothing changes.
  - `PROFILE_COLUMNS = ("depth_to_wall_m", "floor_like_wall", "floor_re", "kink_depth_m", "kink_T_K", "release_depth_m", "growth_rate_per_s", "girin_delta_m_m", "girin_v_s_ms")`
    (The first two were renamed from their `wall_`-prefixed names in the final review's fix round, 2026-11-12, because `wall_` marks the wall-clock timers that the comparisons of runs skip; `ProfileRecord`'s field followed.)
  - `record_columns(width) -> list` gives `SOURCE_COLUMNS` for 27 and `SOURCE_COLUMNS + PROFILE_COLUMNS` for 36; any other width raises ValueError.
  - `source_block(..., extra=None)` appends `extra` (n, 9) to each released row.
  - `SourceSpool(directory, names=None)` stores `self.names`, defaulting to `SOURCE_COLUMNS`. A class attribute `names = tuple(SOURCE_COLUMNS)` serves spools pickled before this change. `append` checks the width against `names`, and `columns()` keys by `names`. `source_table` and `source_columns` take the width from the rows or the spool.
  - `fit_width(rows, width) -> list` pads each block of a list with NaN columns, or cuts it, to `width`. Other forms come back unchanged.

The profile block inside `evaluate` applies only to `profile.resolved` columns. All other columns take today's code, with the resolved columns removed from its thick, regime-2, thin, rarefied and held masks.
- **Gate (decision 6):** a resolved column may shear-spray only if it is windward, `wall_depth >= B_MIN` and `delta_m` is finite, whatever `spray_gates` says.
- **Gradient mode:** λ = 2π/k, r = k_r λ, a growth time of k_t/growth, and a rate of `liq.rho * pi * r**2 / (lam * t)`. It is active where growth > 2 ν_eff k² and `wave_fits(λ, extent)`. Its branch is `BRANCH_PROFILE`, and `we_s` is `profile.we_s`. No critical Weber number is applied.
- **Gas-driven mode:** today's thin branch (or rarefied where `flow.branch` is free-molecular), with b = `wall_depth` in `min(b, λ*/8)`, and `viscous_film(λ*, τ*, wall_depth, liq, nu=nu_eff)` with h_floor = `wall_depth`.
- **Choice:** the shorter growth time wins (k_t/growth against τ*), the quantity the Rayleigh–Taylor comparison already uses. **[plan call.]** Regime 2 is not evaluated on resolved columns, because the thick/thin gate is not used there; it is unreachable for aluminium (measured 2026-09-25). **[plan call.]**
- The Rayleigh–Taylor block runs as today over all columns, reading `rt_extent` where given.
- **Cap:** h = `cap_gradient` where branch 5, `cap_driven` where a resolved column took branch 1 or 2, and the caller's `cap_depth` otherwise. It is stored in `cap_depth`.

- [ ] **Step 1: Write the failing tests.**
  - `test_profile_none_is_todays_evaluate`: with `profile=None` and `rt_extent=None`, every SprayResult field is `np.array_equal` to the existing inputs' results, using the FakeFlow inputs of `test_melt_layer_and_branches`.
  - `test_the_gradient_mode_on_the_profile`: with growth 2000 1/s, k = 1/(50 µm) and ν_eff = 5e-6, the rate equals `rho * pi * (0.17 * lam)**2 / (lam * 1.1 / 2000)`, the branch is 5 and `we_s == profile.we_s`. Growth equal to `2 * nu_eff * k**2` makes the mode inactive.
  - `test_the_gas_driven_mode_reads_the_wall_and_nu_eff`: the rate uses `min(wall_depth, lam/8)` and `gamma_tau` equals `viscous_film(..., nu=nu_eff)`. It wins where τ* < k_t/growth.
  - `test_no_delta_m_no_shear_mode_on_resolved_columns`: under `spray_gates="before"`, the branch is −1 unless the Rayleigh–Taylor mode takes the column.
  - `test_the_rayleigh_taylor_mode_reads_rt_extent` and `test_the_crest_cap_follows_the_chosen_mode`.
  - `test_the_dispersion_setting_defaults_to_table_and_is_validated`: the class default is `"table"`, an unknown name raises ValueError, and a model unpickled without the attribute reads `"table"`.
  - `test_records_append_the_profile_columns_and_spools_keep_their_names` (Review Focus 4): check `record_columns(27)` and `record_columns(36)`. A spool round-trips with the profile names. A spool unpickled without `names` reads 27 columns. A 30-column block is refused. `fit_width` pads and cuts bit for bit.
- [ ] **Step 2:** Run. Expect FAIL. Implement. Run the spray module, `tests/test_reentry_model_melting.py` and `tests/test_reentry_model_shared.py`. Expect PASS.
- [ ] **Step 3:** Commit.

### Task 6: The release run: the rule shared, the cascade in the skins (spec 6.6, decision 5)

**Files:**
- Modify: `reentry_model/shared/melt_profile.py` and `reentry_model/skyscrapers/skin.py` (`Cascade`, `cascade_mass`, `release_top` and `flow_substeps`, plus a new `mode_run`)
- Test: `tests/test_reentry_model_melt_profile.py` and `tests/test_reentry_model_skin.py`

**Interfaces:**
- Produces, shared: `melt_profile.threshold_run(T, flowing, liquid, T_run) -> bool (n, K)`. Where `T_run` [K, (n,)] is finite, it gives the contiguous flowing cells from the top with T ≥ T_run. Where it is NaN, it gives the `liquid` run (today's).
- Produces, in skin.py:
  - `mode_run(sk, T_run)` applies `threshold_run` to the skins' present cells.
  - `Cascade.T_run: np.ndarray = None` is a new last field, so old pickles read None.
  - `cascade_mass(sk, rate, rate_slurry, A_top, dt, slurry=False, T_slurry=None, T_run=None)`. With `T_run`, the mode's cells are `mode_run` at `rate`. Beneath them, with `slurry`, flowing cells at least `T_slurry` warm go at `rate_slurry`, and the first colder cell ends the column. Without `slurry`, the run ends at the mode's last cell. With `T_run=None` the result is bit-identical.
  - `release_top(sk, dm, slurry=False, run=None)`. With `run` given, it admits the flowing cells and splits the book by `run` (phase 0 is the mode's cells, phase 1 the rest). With None it is bit-identical.
  - `flow_substeps`: for a `Cascade` with `T_run`, compute `run = mode_run(sk, given.T_run)` in every flow sub-step and pass it to both functions.

- [ ] **Step 1: Write the failing tests.**
  - Shared: `test_threshold_run_takes_the_warm_cells_from_the_top_and_nan_means_the_liquid_run`, on arrays.
  - Skin: `test_the_profile_run_takes_cells_down_to_its_threshold`. A column at [907.5, 907.3, 907.0, 906.5, 905.0, 902.0, 880.0] K with T_run 907.0 takes cells 0–2 without slurry. With slurry and T_slurry 905.0, cells 3–4 go at the slurry's rate and cell 5 ends the column. Phase 0 equals the mass taken from cells 0–2.
  - Skin: `test_a_nan_threshold_is_todays_liquid_run` (bit for bit) and `test_a_plug_top_colder_than_its_threshold_releases_nothing` (Review Focus 1).
  - Skin: `test_the_profile_release_keeps_the_books`. The FlowBook's `released` mass, enthalpy and momentum equal the column's changes to round-off.
- [ ] **Step 2:** Run. Expect FAIL. Implement. Run the melt_profile and skin modules. Expect PASS.
- [ ] **Step 3:** Commit.

### Task 7: The mode inputs and summary, shared, and the skins' wiring (spec 6.4–6.8)

**Files:**
- Modify: `reentry_model/shared/melt_profile.py` and `reentry_model/skyscrapers/skinbody.py` (`girin_plan`, `SprayPlan` and `SkinBody`)
- Test: `tests/test_reentry_model_melt_profile.py` and `tests/test_reentry_model_skinbody.py`

**Interfaces:**
- Consumes: Tasks 1, 2 and 4–6.
- Produces, shared, in `melt_profile`:
  - `@dataclass class ProfileRecord` holds (n,) arrays: `wall_depth`, `floor_like_wall` (1, 0, or NaN where unresolved), `floor_re`, `kink_depth`, `kink_T` (the first kink's, NaN with none), `release_depth`, `growth`, `T_gradient` and `T_driven`.
  - `profile_modes(cols, windward, sigma) -> (spray.ProfileModes, ProfileRecord)`.
    - It resolves the windward `reducible` columns.
    - It calls `fit_profiles`, `place_wall` with `nu_below`, and `above_wall`.
    - It runs `fastest_waves` only where the layer above the wall has a breakpoint; a straight layer gets growth 0.
    - It computes `couette_viscosity` above the wall and the recorded Weber number: the top segment's density times (u_top − u at the first breakpoint or the wall)² times that depth, over σ.
    - `T_gradient` is the temperature of the cell containing the release depth, or `T_driven` where the release depth reaches the wall. `T_driven` is the temperature of the deepest cell whose centre lies strictly above the wall.
    - The caps are the depths of `threshold_run` at those thresholds.
  - `select_run(branch, modes, record) -> T_run (n,)`: `T_gradient` where branch 5, `T_driven` where a resolved column took branch 1 or 2, and NaN elsewhere.
  - `run_properties(cols, run) -> (rho (n,), mu (n,))`, mass-weighted over the run's cells. **[plan call: the run at the step's start stands in for "the released cells".]**
  - `profile_rows(record, delta_m, v_s) -> (n, 9)`, in `PROFILE_COLUMNS` order.
  - `DISPERSION_COLUMNS = ("floor_like_share", "fluid_like_columns", "depth_to_wall_median_mm", "floor_re_median", "we_s_profile_median", "sprayed_gradient_kg", "sprayed_gas_driven_kg")` and `summary(modes, record, branch, mass) -> dict` over the resolved columns. A step with none gives a share of NaN, a count of 0 and NaN medians.
- Produces, in `skinbody` (wiring only):
  - `girin_plan` reads `model.dispersion`. On `"profile"`, it builds `layer_columns(sk, |flow.tau|)` and calls `profile_modes`.
    - The wave-fit extent is `region_extent(edges, windward & (D >= B_MIN), A_top)`, where D is the wall depth on resolved columns and F elsewhere. `rt_extent` is today's liquid extent.
    - It calls `evaluate(..., profile=modes, rt_extent=...)`, then sets `T_run = select_run(...)`, `run = skin.mode_run(sk, T_run)` and `n_run = run.sum(axis=1)`.
    - It passes `n_run` to `slurry_modes` as its `n_liq`, so the slurry's mode reads the first flowing cell below the run.
    - The liquid mode's `props` become `LiquidProperties(rho_run, mu_run, sigma)` on resolved columns and the liquid's elsewhere.
  - `SprayPlan` gains `T_run = None`, `resolved = None` and `profile = None` (the (n, 9) rows). `__post_init__` passes `T_run` to the `Cascade`. In `source_block`, the liquid mode's rows get `extra = profile`; the slurry's get the same array with `release_depth_m` and `growth_rate_per_s` set to NaN. **[plan call.]**
  - `SkinBody.record_names()` gives `record_columns(27)` on `table` and `record_columns(36)` on `profile`, as `self.spray.dispersion` says. The spool is created with those names.
    - `__init__` raises ValueError for a spray model on `"profile"` with the flow off.
    - `advance` counts droplets with `p.count(p.props.rho)`, which is identical on `table` because `props` is the liquid there.
    - `dispersion_stats()` returns {} on `table` and `melt_profile.summary(...)` on `profile`. It goes last in `skin_stats`, after `spread_stats`.

- [ ] **Step 1: Write the failing shared tests** on synthetic `LayerColumns`:
  - `test_a_floor_like_column_gets_its_wall_at_the_kink_and_no_gradient_growth`: a bridge-like top over slurry with a five-fold drop and Re_f below 1 gives a straight layer above the wall and growth 0. This is the spec's anti-goal: no gradient release from a straight film on a floor. `T_driven` is the temperature of the cell above the kink.
  - `test_a_fluid_like_curved_column_gets_its_release_depth_and_threshold`
  - `test_select_run_and_run_properties`
  - `test_summary_of_an_empty_step_is_nan` (Review Focus 3)
- [ ] **Step 2: Write the failing skin tests:**
  - `test_on_profile_a_floor_like_column_sprays_by_the_gas_driven_mode_at_its_kink` (branch 1, using `girin_flow` with the Girin closure)
  - `test_on_profile_a_column_whose_top_does_not_move_falls_back_to_today`: its SprayResult entries equal those of `table`.
  - `test_profile_records_have_36_columns_and_slurry_rows_carry_no_release_depth`
  - `test_a_profile_skin_body_sprays_and_keeps_its_books`: on `coarse_sphere_mesh`, built like `test_a_skin_body_sprays_by_girin_...` with `SprayModel(..., dispersion="profile")`. Mass and energy close, `DISPERSION_COLUMNS` come last in the history row, and every record has 36 columns.
  - `test_a_table_skin_body_has_no_dispersion_columns`
- [ ] **Step 3:** Run. Expect FAIL. Implement. Run the melt_profile, skin, skinbody, spray, coupled and shared modules. Expect PASS.
- [ ] **Step 4:** Commit.

### Task 8: Command line (default `profile`), existing tests, documentation, same results and cost

**Files:**
- Modify: `reentry_model/cli.py`, `analysis/skyscrapers_port/same_results.py` (and its `test_port_tools.py` and README), `README.md` (the spray section), `docs/model_assumptions.md` and `CLAUDE.md` (the options list)
- Create: `analysis/kinked_dispersion/plan_cost.py`
- Test: `tests/test_reentry_model_cli.py` and `tests/test_reentry_model_restart.py`

**Interfaces:**
- `run --dispersion table|profile`, with argparse default None, sits after `--slip-band` in the melt group. Add `'--dispersion'` after `'--slip-band'` in `BASE_RUN_OPTION_STRINGS`, which pins the list and its order.
  - `dispersion_of(args)` returns the given value, else `dispersion.DEFAULT_DISPERSION` (Task 1) on a spraying skin run, else `"table"`.
  - A given value on a run that is not a spraying skin run is refused (exit 2) with "--dispersion needs the skins' spraying: --surface-model skin --melt on --skin-flow on".
  - `spray_model_of` passes `dispersion_of(args)` for the skins and `"table"` for the element model.
  - `skin_name_suffix(..., dispersion="profile")` appends `_dispersion-table` when `table`.
  - `info["skin"]["dispersion"] = dispersion_of(args)` is always written on spraying skin runs ("always: the default changed").
- `resume --dispersion table|profile`, default None (the saved model's; `table` for a state saved before the option existed), is refused (exit 2) unless the saved body is a `SkinBody` with a spray model.
  - Applying it sets `run.body.spray.dispersion`.
  - It pops `DISPERSION_COLUMNS` from `state["row"]` and updates the row with `dispersion_stats()`.
  - In-memory blocks go through `spray.fit_width(rows, len(run.body.record_names()))`, and the new spool takes `run.body.record_names()`.
  - `thermal_info["dispersion"]` is always written for skins with a spray model. `saved.dispersion` is set in the arguments loop, and the start-up print line gains "dispersion {}".
- `same_results.py run --append "ARGS"` appends the given arguments to every listed case.

- [ ] **Step 1: Write the failing tests.**
  - cli: `test_dispersion_is_profile_by_default_on_spraying_skin_runs_and_reaches_the_spray_model`, modelled on `test_the_spray_gates_are_2026_11_10_by_default_...`. It checks skins with the flow on, the element model on `"table"`, and the explicit `table` reaching the model.
  - cli: `test_dispersion_is_refused_off_the_skins_spraying` (element model, heat-only skins, no melting: exit 2).
  - cli: `test_the_switch_back_is_named_and_every_skin_run_records_its_dispersion`.
  - restart: `test_resuming_with_dispersion_table_cuts_its_columns_and_a_state_without_it_resumes_on_table(flown, tmp_path)` (Review Focus 4). `flown` now flies on `profile`. Resuming it with `--dispersion table` writes 27 npz keys and no `DISPERSION_COLUMNS`. A state whose spray model lost its `dispersion` attribute (simulated as the spray-gates test does) resumes on `table`. Resuming with nothing changed repeats the flight to the last bit, which the existing test already covers.
  - restart: add `--dispersion profile` on a saved flight without spraying to the cases of `test_the_resume_command_refuses_what_it_cannot_resume`.
  - port tools: `test_append_adds_arguments_to_every_listed_case`.
- [ ] **Step 2:** Run. Expect FAIL. Implement.
- [ ] **Step 3: Existing tests under the new default.** Run `tests/test_reentry_model_cli.py`, `tests/test_reentry_model_restart.py`, `tests/test_reentry_model_coupled.py` and `tests/test_reentry_model_thermal_verification.py`. A test that fails because it pinned Girin's theory's spraying (record counts, masses, 27 record columns, settings keys) gets `--dispersion table` added to its run, since its subject is that theory. A test about the flight machinery keeps the default, with only its expected record width or settings keys updated. List every test changed, and how, in the commit message.
- [ ] **Step 4: Documentation.** In `docs/model_assumptions.md`, state that:
  - the default is the profile path (Asha, 2026-11-11);
  - the floor test is a scaling argument, assumed until stage 6;
  - the broken-line derivation is ours, checked against Girin Eq. 1 and Ghosh 2017;
  - the kink rests on the bridge's 906.4 K end point;
  - the threshold is 3.00 on `profile` and 3.08 on `table`;
  - the constants are assumed.

  The README's spray section describes the option, the switch back and the records. CLAUDE.md lists `--dispersion` and the new shared module.
- [ ] **Step 5: Same results.** Run `same_results.py run --cases skins --append "--dispersion table"` and `--cases element,resume` (no flag) into `output/kinked_check/after`, with the same `--state` as Task 0. Then run `same_results.py compare output/kinked_check/before output/kinked_check/after --report`.
  Expect exactly one difference per skin case, in its JSON's new `dispersion` key, and none anywhere else: every history column, every record, and the `element` case's whole output. If anything else differs, stop and report; never edit a baseline.
- [ ] **Step 6: Cost.** `plan_cost.py --state <d050>` calls `girin_plan` five times with the spray model on `"table"` and five times on `"profile"`. It uses the body's own inputs: `body.skins`, `body.last_flow`, `body.spray_plan.state`, `body.spray_plan.theta`, `body.spray_plan.area`, `body.columns.edges`, `body.transverse_radius()` and `nose_radius=body.nose_radius()`. It reports the median added time.
  The macro step's wall time is the median increase of `wall_skins_s + wall_3d_s + wall_move_s + wall_rebuild_s` between consecutive rows of `d050__handdown-off.rows.csv` from 177.5 to 179.5 s, per macro step. Measured 2026-11-11: about 3,100 s per 0.5 s of flight. Pass if the added time is under 5 % of the step's wall time, and record both numbers.
- [ ] **Step 7:** Run the unit tier. Expect the Task 0 count plus the new tests, all passing.
- [ ] **Step 8:** Commit. Then run the whole-branch review (requesting-code-review) before Task 9.

### Task 9: Comparison flights (stage 5; Mac Studio by runbook)

**Files:** Create `docs/superpowers/plans/<date written>-kinked-dispersion-flights-macstudio.md`, with the section layout of `2026-11-10-skin-fixes-test-flights-macstudio.md`.

- [ ] **Step 1:** Write the runbook for four flights (each over two hours of work in total; Asha, 2026-10-31):
  - **(a) and (b):** resume `d050__handdown-off.state.pkl` (t = 178.5 s, copied to the Mac Studio by hand as the d100 state was) to `--t-max 179.0`, with `--spray-gates 2026-11-10 --slip-band on` and `--dispersion table` or `profile`. That is about an hour each at the measured rate.
  - **(c) and (d):** a fresh 100 mm skin flight with today's options plus `--spray-gates 2026-11-10 --slip-band on`, with each setting, to `--t-max 30` (first liquid at 24.5 s).
  - The runbook covers pushing the branch, the checkout, a timing check and bringing the results back. Commit it, and ask Asha before pushing the branch for the Mac Studio.
- [ ] **Step 2:** Once the results are back, write a report to Asha in `reports/` following her rules. It compares the two settings on the floor-like share, where the walls land against the 906.4 K isotherm, the mass sprayed by each mode (gradient, gas-driven, Rayleigh–Taylor and slurry) and the droplet sizes. The default is already `profile` (decision 9). The report says what the change does and flags anything that argues against the default.

### Task 10: The viscous (Orr–Sommerfeld) check (stage 6; analysis only)

**Files:** Create `analysis/kinked_dispersion/orr_sommerfeld.py`, `analysis/kinked_dispersion/test_orr_sommerfeld.py`, `analysis/kinked_dispersion/floor_sweep.py` and a report in `reports/`.

**Interfaces:** `growth_rate(k, y, U, mu, rho, sigma, top="free") -> complex c`. This is the most unstable Orr–Sommerfeld eigenvalue of a column.
- y is the cell faces (m, downward), U the base speed at the cell centres, and mu and rho are per cell.
- The bottom is a no-slip wall. The top is a free surface with tension (`"free"`) or a moving wall (`"wall"`, for the Couette validations).
- Velocity and both stresses are continuous at the faces. The method is the builder's.

- [ ] **Step 1:** If Albert & Charru (2000) can be obtained, read it first and note anything that changes the validations.
- [ ] **Step 2: Validation tests.**
  - Girin's profile at a Reynolds number V_s δ/ν of 1e4 and of 1e5 is within 5 % of the broken-line growth, with a smaller error at 1e5.
  - Two deep fluids in uniform shear (equal densities, no tension) recover Hooper & Boyd's short-wave growth σ = a2 m(1 − m)²/[2(1 + m)] × a2/(k² ν2) within 5 % at α = k(ν2/a2)^(1/2) of at least 5, where m = μ1/μ2 (notes `fulltext_viscous_papers_read.md`, line 37).
  - The sign of the long-wave two-layer Couette mode follows Charru & Hinch: stable where the thinner layer is the less viscous, and unstable otherwise.
- [ ] **Step 3:** `floor_sweep.py` runs a few dozen sampled flight columns, from d050 and d100 at several states, including the Task 4 fixture, through `melt_profile.LayerColumns`. It also runs synthetic two-layer sweeps across Re_f = 1. It writes the viscous growth against the floor Reynolds number and against the slope ratio at the kink.
- [ ] **Step 4:** Write a report to Asha in `reports/`, following her rules. It says where the viscous growth vanishes against Re_f, recommends `FLOOR_RE` and `KINK_DROP` with their evidence, and says what changes if they replace 1.0 and 0.2. Changing them is her decision.
- [ ] **Step 5:** Commit.

---

## Calls made in this plan (beyond the spec's own)

1. The fit inserts breakpoints greedily and then prunes them, instead of an exact search for the fewest (Task 4).
2. A layer moving as one piece has one breakpoint but no kink, because the fit ends at zero at the rigid front. The spec says "none" (Task 4).
3. "Grows faster" compares the growth times the Rayleigh–Taylor rule already uses (Task 5).
4. Regime 2 is not evaluated on resolved columns (Task 5).
5. The liquid mode's droplet properties are the run's at the step's start (Task 7).
6. The slurry's record rows carry the column's wall and kink values, but NaN release depth and growth (Task 7).
7. The gas-shear gate passes at a median ratio of actual to predicted top speed between 0.8 and 1.25 (Task 3).
8. The setting's constructor and class defaults stay `"table"`, so saved states and directly built objects keep Girin's theory. The command line sets `"profile"` on spraying skin runs, which is what makes it the default for every run (Tasks 5 and 8).
9. The skins' per-cell shear rate stays in `skin.py` as flow numerics, not spray physics. The level-set architecture supplies its own (Task 3).
