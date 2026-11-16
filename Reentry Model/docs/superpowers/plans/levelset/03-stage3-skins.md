# Level-set stage 3 — skins on markers, both modes, flow solve off — Implementation Plan

> **Superseded in part (2026-11-15).** The design's amendment of 2026-11-13 (§21, one skin depth everywhere, the gas
> surface φ evolved and the skin base β derived) retired the two skin modes this plan builds. Its Tasks 1 to 9 and 9b
> stand as built; **Tasks 10 to 12 are superseded** by `03b-stage3-one-skin-depth.md`: Task 10 (the two modes) is
> retired by that plan's Tasks 5 and 6, and Tasks 11 (the verification cases and cost) and 12 (the documents and the
> ledger) are redone as its Tasks 7, 11 and 12. What stage 3 is as built is recorded in `docs/model_assumptions.md`
> section 16 and in the design's §21.8 and §21.12. This plan is kept as the record and is not edited further.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the level-set body its skins: layered columns of cells carried by markers on the gas surface, in the
two modes of design §20.6 (full-depth over thin melt, top-zone over thick melt), heated by the shared gas laws, coupled
to stage 2's cut-cell heat solve, flowing and running off between markers, and spraying through the shared spray
model — with the volume flow solve still off (thick melt held), so that stage 3's verification (§20.11 item 3) passes.

**Architecture:** Three layers. (1) `reentry_model/shared/` holds every physical law the level set uses; Task 1 moves
the few laws still inside skyscrapers modules there, behaviour-preserving. (2) `reentry_model/levelset/` holds the level
set's own numerics: `skin.py` (the cell columns: storage, grading, curved cross-section, conduction, layered flow),
`markers.py` (markers on φ: create, snap, carry, area, split, merge, trap, neighbours), `runoff.py` (exchange between
markers' skins), `loads.py` (per-marker gas side from the shared laws), `spraying.py` (the per-marker spray plan
through the shared `SprayModel`), `modes.py` (the thin test per marker and the mode changes) and `skinbody.py`
(`LevelSetSkinBody`, one macro step: loads, trial pass, cut-cell solve, real pass, base motion, surfaces, markers,
modes). (3) Verification scripts in `analysis/levelset/` and tests in `tests/levelset/`. No run driver yet: the
command-line `--architecture levelset` still exits 2 (the driver comes with the flights, stage 5).

**Tech Stack:** Python 3.12 (`drama_env`), numpy, scipy, pyamg, scikit-fem (via stage 2); pytest with pytest-xdist.
No new dependency.

**Spec:** `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`: §9 (markers, §9.1–9.2,
§9.4–9.6), §20.3 (β, φ, ψ, crusts), §20.5 (the thin test), §20.6 (skins in two modes), §20.8 (shadow test), §20.10 (the
parity target), §20.11 item 3 (this stage's verification), §13 item 3 (its tests), §12.2 (history columns), §4.2 (the
shared package; skins never shared). Parity ledger: `docs/superpowers/plans/levelset/parity-ledger.md`. Reference for
the skyscrapers' skin physics (what parity means, read on 2026-11-11 from `main` ae20dd3, with file and line numbers):
the summary handed to each implementer as `stage3-skyscrapers-reference.md` in the plan's SDD workspace (Task 0).

## The parity target (pinned here, design §20.10 and §20.14 item 1)

Pinned at **`main` ae20dd3** (2026-11-11): the skyscrapers build through sub-plan 22 Task 5k (Tasks 5g–5j included),
plus the stage 0 + S refactor. Not in the target: branch `we-critical-308` (unmerged; ledger row 8), the kinked-profile
dispersion design of f6699e3 (not built), the refrozen-metal-in-skins draft (ledger row 7). The skin rules the level
set reproduces (by amendment, never by import): cells of 10 µm (production; 20 µm dev) over a 0.4 mm fine zone, each
deeper cell 10 % thicker; a velocity per cell; implicit conduction with radiation from the top cell; the layered column
flow (gas shear on the top cell; pressure gradient, deceleration, thickness and capillary pushes through the depth;
viscosity per cell at its own shear rate; anchored rigid cells at rest; crusts floating as one piece); runoff by
implicit upwind exchange with the arrival rule and landing conduction; the release cascade (liquid first, slurry with
`--slurry-spray on`); Girin's thick and thin branches, the Rayleigh–Taylor droplet, the gates of `--spray-gates
2026-11-10`, the slip band (`--slip-band`, default off); the windward loads from a fitted sphere; `--leeward centroid`
as the one compatibility switch. Rules the level set replaces (ledger rows 2, 5, 7): the hand-down and part count of
the columns (§20.6's hand-down below the stub instead), the leeward centroid test (the shadow test by default).

## Global Constraints

- Paths are relative to `Reentry Model/` in the worktree `.claude/worktrees/levelset-stage3` (branch
  `levelset-stage3`, from `main` ae20dd3). Interpreter `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`, run
  from `Reentry Model/`. Paths contain spaces: quote them.
- **Shared physics only (Asha, 2026-11-11: "make sure the code for level set is using the physics in the shared
  folder").** Every physical law, material property, gas-side quantity and physical constant the level set uses comes
  from `reentry_model.shared` (material and its `T_rigid`, `T_feed`, `T_REF`, densities, liquid properties; `rheology`;
  `surface_flow`; `heating`; `gas`; `spray` and `dispersion`; `droplets`; `thermal.SIGMA_SB`; and what Task 1 moves).
  `reentry_model/levelset/` defines no physical constant of its own; its numbers are numerical settings (cell sizes,
  tolerances, sub-step lengths, the thin test's coefficients) and each is named in a module-level constant or a
  settings field with its source. Task 1's test `test_levelset_imports_only_shared_physics` enforces the import side.
- **Skins are never shared** (design §4.2): nothing in `reentry_model/levelset/` imports or copies `skin.py`,
  `skinbody.py`, `columns.py`, `film.py`, `spreading.py`, `movemesh.py`, `mesh.py`, `restart.py`, `coupled.py`, `cli.py`,
  `body.py` or `thermal/*_backend.py`. Reading them to match the physics is expected; copying code is not.
- **Changes to `shared/` are behaviour-preserving** and are accepted only if
  `analysis/skyscrapers_port/same_results.py compare` against the refactor's baseline says `identical` (Task 1) and the
  shared-package tests pass.
- Do not touch `DRAMA Sphere Experiments/`, `Spheral/`, `reentry_model/droplet_reentry/`, `reentry_model/wake_column/`,
  the `subplan20` worktree, `prototype/`, or other sessions' files. No prototype folder: stage 3 is built in the package
  directly (Asha's rule of 2026-11-11 for work after the refactor; it supersedes design §14's "prototyped in
  `prototype/levelset/`" for this stage).
- SI units throughout; level sets negative inside metal; normals point out of the metal; enthalpy above the shared
  `T_REF` (293 K) per kilogram, as the shared `Material`.
- Tests in `tests/levelset/`, one file per module (`test_levelset_<module>.py`). The fast ones (all of a task's unless
  stated) run in the unit tier; whole-plate and box runs are marked `reference`. Use the test-only materials of stage 2
  (`ConstMaterial`, `LinearMeltMaterial`, moved in Task 2 to `tests/levelset/materials.py`) where a test needs simple
  properties, and the shared `Material.from_drama_json()` (AA7075 Scheil) where parity needs the real alloy.
- Run tests only on mains power (`pmset -g batt | head -1` says `AC Power`). A run longer than about two hours goes to
  the Mac Studio through a runbook (Asha, 2026-10-31).
- Commit after each task (named paths only, never `git add -A`). Message: one sentence on the physics or behaviour
  that changed, a blank line, `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Do not push.

## Decisions for Asha (asked with this plan; the plan follows the recommendation until she answers)

- **D1 — ledger row 1, the "piecewise shoulder fit".** It does not exist in the skyscrapers code, designs or plans (a
  search of `reentry_model/` and `analysis/` on 2026-11-11 finds none). Recommendation: drop it; the level set fits one
  sphere to φ's windward half (her rule of 2026-11-11) with the 1.67 flat-face cap, as the skyscrapers fit one sphere
  to their windward cap. Cost if wrong: one function in `loads.py` later.
- **D2 — ledger row 4, the slip-band depth.** The momentum-balance depth is already shared code
  (`shared/spray.momentum_layer` and `spray_layer`, fed by `surface_flow.shear_moment`); it reads only the gas side,
  never a skin layer. What the level set supplies is only the film depth b and the regime depth, which §20.6 already
  fixes (the skin's flowing depth; on a top-zone skin, φ to ψ). Recommendation: close row 4 as **shared**, no
  amendment needed.
- **D3 — ledger row 6, deposit spreading.** The skyscrapers' spreading moves only handed-down metal where a column base
  would fold, and is off by default (it needs `--hand-down on`). The level set's volume mesh does not move, so nothing
  folds. Recommendation: **not applicable**; revisit if a level-set surface shows a reason for it.

- **D4 — where the release cascade lives** (passed from Asha's forked session, 2026-11-11). The cascade (liquid cells
  first from the top, then the flowing slurry with `--slurry-spray on`, never rigid metal, each cell at its own enthalpy
  and velocity) could move into `shared/` as a function on cell arrays, or stay in each architecture's skin code with a
  test that compares the two releases on the same column. Recommendation (the forked session's, which I share): keep
  it per architecture — it acts on each skin's own cells, which are never shared — and pin it with
  `test_the_cascade_matches_the_skyscrapers_release_on_one_column` in Task 8 (the same column of cells released by both
  over one macro step: the same masses, enthalpies and momenta within 1e-12; the test may import skyscrapers' `skin.py`
  because tests are not the level-set package). Until Asha answers, Task 8 builds it per architecture with that test.

## Rulings made while planning (recorded; each with its cost if wrong)

1. **Markers sit on φ** (the gas surface), one per φ-patch at creation, and carry their skin's base at x − H n.
   Loads, spray and curvature read φ; the base defines β. Cost: re-deriving β if markers belong on β instead.
2. **Marker area and the faces between markers** come from a Voronoi construction in each marker's tangent plane over
   its graph neighbours (Task 5), symmetrised so a face is shared by both markers, and scaled so the markers' areas
   sum to φ's patch area. A skin's area therefore follows the surface; a growing area thins the skin and the skin then
   draws metal at its base (§9.2), with mass and enthalpy exact. Cost: a different tessellation in one function.
3. **Flow off at the thin–thick boundary.** With the volume flow held (this stage), a full-depth skin exchanges with a
   neighbouring top-zone skin only through the top 0.4 mm, layer by layer; its deeper cells do not run into the held
   volume. §20.6's flux boundary into the volume flow comes with the flow solve (stage 4). Cost: none for stage 3's
   tests; the ring-forming flux into pools is stage 4's.
4. **Momentum given to the held volume** (thin → thick, a top-zone base) is booked as `momentum_to_volume` and reported,
   as the skyscrapers book the normal part to the wall. Cost: none (reported, not lost silently).
5. **The q6p box** (§13 item 3: 8 × 10 × 10 mm, 6 MW/m², 2 s) is run on the level set, which does not fold at edges;
   its recession at the heated face's centre is judged against the cold plate's one-dimensional reference within 5 %
   (the face centre is 4 mm or more from any edge), and its books exactly. The skyscrapers dropped the box for the
   plate because their columns folded at its edges, so there is no skyscrapers box to compare with. Cost: one test.
6. **No run driver in stage 3.** `LevelSetSkinBody` is driven by tests and analysis scripts at prescribed gas states
   and decelerations; the trajectory coupling and `--architecture levelset` come with the flights. Cost: none now.
7. **Parts within a macro step** follow the skyscrapers' rule as physics of the skins' time stepping: no part longer
   than 0.0625 s while any top cell is at or above `T_rigid`, and a part whose base moves exceed a quarter of the local
   marker spacing is rerun shorter (at most six levels). Cost: one constant.

## Review Focus

- **A marker whose skin has nothing flowing and nothing heating it** must stay bit-for-bit unchanged through a macro
  step (no drift from regrid, snap or area) — pinned by `test_rigid_motion_leaves_skins_unchanged` (Task 5).
- **A skin thinner than one fine cell** (just created, or nearly all sprayed) must keep closed books and positive
  cell masses — pinned by `test_a_skin_reduced_to_one_cell_keeps_its_books` (Task 3).
- **A marker trapped between two surfaces** (a fold) must return its whole skin to the volume with its enthalpy —
  pinned by `test_a_trapped_marker_returns_its_skin_exactly` (Task 5).
- **A mode change on the step a marker is created** — pinned by `test_a_new_marker_may_change_mode_only_after_the_hold`
  (Task 10).
- **The shared spray model fed through the level set** must give exactly what it gives fed directly — pinned by
  `test_marker_release_equals_the_shared_model` (Task 8).

---

## Task 0 (controller only): the workspace

- [ ] Copy `/private/tmp/claude-501/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/eb7e9b74-487b-431a-863c-b8a3214cb080/scratchpad/skyscrapers_skin_reference.md`
  into the SDD workspace as `stage3-skyscrapers-reference.md`; every implementer brief points at it. Record
  `same_results.py`'s baseline directory:
  `../.claude/worktrees/shared-refactor/Reentry Model/output/refactor_check/before` (five cases, 26 files).

### Task 1: Shared physics only

**Files:**
- Create: `reentry_model/shared/lubrication.py` (moved `lubrication` from `film.py`), `reentry_model/shared/nose.py`
  (moved `fit_sphere`, `NOSE_CAP_ANGLE`, `NOSE_CAP_FACTOR` from `body.py`; plus the extracted rules below)
- Modify: `reentry_model/film.py`, `reentry_model/body.py`, `reentry_model/skinbody.py` (import from shared, re-export
  the moved names); `reentry_model/shared/spray.py` (gains `SLURRY_SIGMA`, `region_extent`, `slurry_properties`);
  `reentry_model/levelset/surfaces.py`, `reentry_model/levelset/cutthermal.py`
- Test: `tests/test_reentry_model_shared.py` (append), `tests/levelset/test_levelset_imports.py` (new)

**Interfaces:**
- Produces: `shared.lubrication.lubrication(tau, G, b, delta_m, mu_l, b_layer=None) -> (V_s, q, shear_rate, thick)`;
  `shared.nose.fit_sphere(points) -> (centre, radius)`, `NOSE_CAP_ANGLE = 30.0`, `NOSE_CAP_FACTOR = 1.67`,
  `shared.nose.nose_radius(fitted, bound, concave, n_points, fallback, size_feedback) -> float` (the selection now
  inline in `SkinBody.nose_fit`: `min(fitted, bound)` if convex, `bound` if concave or coplanar, `fallback` with fewer
  than 4 points, `fallback` throughout under `initial`); `shared.nose.LEEWARD_TESTS = ("facing", "centroid")`,
  `shared.nose.ahead_of(points, centre, v_hat) -> bool array`, `shared.nose.load_angles(theta, windward) -> theta'`
  (θ where windward, else max(θ, π/2 + 1e-6)); `shared.spray.SLURRY_SIGMA = 0.80`,
  `shared.spray.region_extent(edges, selected, areas)`, `shared.spray.slurry_properties(material, T, gamma, rho)
  -> LiquidProperties` (the viscosity and density the skins build in `slurry_modes`); `levelset.surfaces.
  molten_rigid_surface(mesh, phi, T, T_rigid, core_point=None, band=None)` with `T_rigid` now required (callers pass
  `material.T_rigid`); `cutthermal.SIGMA_SB` imported from `shared.thermal`.

- [ ] **Step 1: Write the failing tests.** In `test_reentry_model_shared.py`: `test_moved_film_and_nose_names_are_the_shared_objects`
  (`reentry_model.film.lubrication is reentry_model.shared.lubrication.lubrication`, the same for `body.fit_sphere`,
  `body.NOSE_CAP_ANGLE`, `body.NOSE_CAP_FACTOR`, `skinbody.SLURRY_SIGMA`, `skinbody.region_extent`,
  `skinbody.LEEWARD_TESTS`); `test_shared_nose_rules_match_the_old_inline_rules` (on a hand-built set of points and
  angles: `nose_radius` and `load_angles` equal what `SkinBody.nose_fit` and `SkinBody.load_angles` gave before the
  change — record the old values first, from the unmodified tree, into the test). In `test_levelset_imports.py`:
  `test_levelset_imports_only_shared_physics` (ast scan of every `reentry_model/levelset/*.py`: every `reentry_model`
  import is `reentry_model.shared[.x]`, `reentry_model.levelset[.x]`, a relative import inside `levelset`, or the
  package root constants; and none of the skins' modules named in Global Constraints) and
  `test_levelset_defines_no_physical_constants` (no module-level name in `levelset/` equals the value of a shared
  physical constant: `SIGMA_SB`, 895.1, 910.0, 293.0).
- [ ] **Step 2: Run them; they fail** (names missing; `T_RIGID` and `SIGMA_SB` defined in `levelset`).
- [ ] **Step 3: Move and extract.** `git mv`-free moves (the functions are cut and pasted unchanged, the old modules
  import and re-export them). Extract the three inline rules of `skinbody.py` (`nose_fit`'s selection, `_leeward`'s
  plane test, `load_angles`, the slurry property build in `slurry_modes`) into the shared functions above and make
  `skinbody.py` call them, changing no arithmetic. In `levelset/`: delete `T_RIGID`, make `T_rigid` required and update
  its callers and tests (the tests pass the value they mean, from their test material); import `SIGMA_SB` from
  `reentry_model.shared.thermal`. Update the README and `docs/model_assumptions.md` §10 where they name `T_RIGID`.
- [ ] **Step 4: Prove no behaviour changed.** `"$PY" analysis/skyscrapers_port/same_results.py run --code . --out
  output/stage3_check/task1` for its five cases (about 5 minutes), then `compare` against the baseline of Task 0 →
  `identical: 26 files compared`, exit 0. Run `tests/test_reentry_model_shared.py`, `tests/test_reentry_model_cli.py`,
  `tests/test_reentry_model_spray.py`, `tests/test_reentry_model_film.py`, `tests/test_reentry_model_skinbody.py -m "not
  reference"` and `tests/levelset` → all pass. Any difference: stop and report; never "fix" a result.
- [ ] **Step 5: Commit.**

### Task 2: The skin's cells — storage, grading, curved cross-section, regrid

**Files:** Create `reentry_model/levelset/skin.py`; `tests/levelset/materials.py` (the stage 2 test materials, moved
out of `test_levelset_cutthermal.py`); test `tests/levelset/test_levelset_skin.py`.

**Interfaces:**
- Produces: `SkinSettings(cell=10e-6, fine_zone=0.4e-3, growth=1.10, stub=0.8e-3, substep=0.010, swell=True)` with
  `SKIN_PRESETS = {"production": 10e-6, "dev": 20e-6}`; `Skins(material, area, kappa, settings)` holding n columns in
  padded arrays `m (n,K)`, `H (n,K)` (J above `T_REF`), `U (n,K,3)`, `count (n,)`, plus `area (n,)` (top area A₀),
  `kappa (n,2)` (φ's principal curvatures, positive convex); methods `T()`, `density(T)` (the shared mixture law
  1/((1−f_l)/ρ_s + f_l/ρ_l), or ρ_s with `swell=False`), `dz()` (each cell's thickness from its mass and the
  cross-section at its depth), `depth_above(T_level)`, `thickness()`, `targets()`, `regrid()`, `add_bottom(i, m, H)`,
  `take_top(i, ...)`, `totals() -> (mass, enthalpy, momentum)`.
- The cross-section at depth z below the top is A(z) = A₀ (1 − κ₁ z)(1 − κ₂ z) (§20.6). A cell's thickness solves
  ∫ A(z) dz = m/ρ over its span, top to bottom (closed form: the cubic in z, solved by Newton from the flat value).
  Target thickness `cell + (growth − 1) · max(z_top − fine_zone, 0)`; `regrid` merges a cell below 0.5 of its target
  into its neighbour and splits one above 1.5, exactly in mass, enthalpy and momentum, at most 64 passes.

- [ ] **Step 1: Write the failing tests:** `test_cells_grade_ten_percent_below_the_fine_zone` (a 3 mm solid column of 10
  µm cells regrids to 40 cells in the top 0.4 mm and each deeper target 10 % thicker within 1e-12 relative);
  `test_regrid_conserves_mass_enthalpy_and_momentum` (random column, 1e-14 relative); `test_the_curved_cross_section_holds_a_spherical_shell`
  (κ₁ = κ₂ = 1/R, R = 25 mm: a column of mass m has the thickness t with A₀ ∫₀ᵗ (1 − z/R)² dz = m/ρ, analytic
  shell volume (4π/3)(R³ − (R − t)³) per 4πR² of top area, within 1e-12 relative; and a concave κ = −1/R grows with
  depth); `test_a_flat_skin_has_the_flat_thickness` (κ = 0: t = m/(ρ A₀) exactly); `test_density_is_the_shared_mixture_law`
  (`density(T)` equals the formula on the shared Material's `liquid_fraction` at five temperatures).
- [ ] **Step 2: Run; they fail.** **Step 3: Implement.** **Step 4: Run; they pass.** **Step 5: Commit.**

### Task 3: The skin's heat — conduction, the base law, removal, draw and hand-down

**Files:** Modify `reentry_model/levelset/skin.py`; test `tests/levelset/test_levelset_skin.py`.

**Interfaces:**
- Produces: `Skins.conduct(dt, q_top, T_amb, emissivity, T_b, trial=False) -> HeatBook` (backward Euler, Newton in T
  with a backtracking line search; the top cell gets q_top·A_top and loses εσ(T₀⁴ − T_amb⁴)·A_top with the shared
  `SIGMA_SB`; harmonic conductances between cells with the cross-section at each interface; the base held at T_b
  through the bottom half-cell; with `trial=True` it carries S = dT/dT_b and the book's `dE_base_dTb`). `HeatBook`
  fields per column: `E_top, E_rad, E_base, E_side` (J), `dE_base_dTb` (J/K, trial only); E_base is the residual of the
  books so energy closes exactly. `Skins.base_law(dt) -> (a, b)` per column (W, W/K; heat into the volume, b ≤ 0 clipped
  for diagonal dominance), as the skyscrapers' trial pass. `Skins.remove_feed()` (the contiguous top cells at or above
  the shared `T_feed`, the flow-off removal of the plates) `-> (mass, enthalpy)`. `Skins.draw(T_b, g_base, m_min)` (a
  column below `m_min` takes the shortfall at its bottom at the profile enthalpy of T(s) = T_b − g_base·s capped at
  `T_rigid`, 8-point Gauss–Legendre) `-> (mass, enthalpy)` per column; `Skins.anchored()` (the unbroken run of rigid
  cells from the bottom) and `Skins.crust()` (rigid, not anchored); `Skins.hand_down(stub) -> (mass, enthalpy, volume)`
  per column: the anchored cells below the deepest `stub` of the anchored run leave at their own enthalpy (§20.6).
  `minimum_mass = ρ_s A_base · stub` (as flown by the skyscrapers: no δ_m term).

- [ ] **Step 1: Write the failing tests:** `test_steady_ablation_recedes_at_the_energy_limited_speed_with_the_exponential_profile`
  (one flat column, q = 1e8 W/m², T₀ = 300 K, a sharp-melting test material, 20 µm cells, 0.5 ms sub-steps, 5 mm all
  fine, `remove_feed` each sub-step: recession speed within 3 % of q/(ρ(h(T_feed) − h(T₀))), the profile within 2 % of
  (T_feed − T₀) of T₀ + (T_feed − T₀) exp(−v z/α)); `test_conduction_books_close` (radiation and a base law on, 1e-12
  relative); `test_the_trial_sensitivity_is_the_derivative_of_the_base_heat` (central difference of E_base in T_b,
  1e-6 relative); `test_draw_and_hand_down_are_exact` (mass and enthalpy leaving equal what the volume is told, 1e-14);
  `test_a_crust_is_never_handed_down` (rigid–molten–rigid column: only the bottom run beyond the stub leaves);
  `test_a_skin_reduced_to_one_cell_keeps_its_books`.
- [ ] **Step 2–4: fail, implement, pass.** **Step 5: Commit.**

### Task 4: The layered column flow

**Files:** Modify `reentry_model/levelset/skin.py`; test `tests/levelset/test_levelset_skin.py`.

**Interfaces:**
- Produces: `Skins.momentum_step(dt, tau, push, g_n)` per component, implicit:
  m_i (u_i' − u_i)/dt = A [D_(i−1)(u'_(i−1) − u'_i) + D_i(u'_(i+1) − u'_i)] + push·V_i + τ A_top [top cell], with
  D = 1/(dz_i/(2η_i) + dz_j/(2η_j)) at each interface (a flowing cell on rigid metal drags over its own half-thickness),
  η per cell from `shared.rheology.viscosity` at the cell's own shear rate by fixed point (at most 5 passes, tolerance
  1e-3), anchored cells at rest, a crust one unknown with 1e9 Pa s at its faces; returns the momentum given to the
  wall. `Skins.flowing()` (solid fraction below 0.5, i.e. T above the shared `T_rigid`), `Skins.liquid()` (T at or
  above the shared liquidus). `push` (N/m³ per column, tangent) is built by Task 6 from the shared surface flow.
- [ ] **Step 1: Write the failing tests** (constant-viscosity test material unless stated): `test_linear_couette`
  (shear τ on a liquid layer on rigid metal: steady profile linear, top speed τ b/μ within 1e-6); `test_half_parabola_second_order`
  (a body force alone: the steady profile's error falls by four when cells halve); `test_liquid_over_slurry_kinks`
  (shared AA7075 rheology: the profile's slope jumps at the liquid–slurry interface by the viscosity ratio, 1e-3);
  `test_stokes_first_problem` (an impulsively sheared deep layer against u = (2τ/μ)√(νt/π)·… the erfc solution,
  2 % at three times); `test_a_floating_crust_moves_as_one_piece` (its cells' speeds equal, momentum exact 1e-12).
- [ ] **Step 2–4.** **Step 5: Commit.**

### Task 5: Markers

**Files:** Create `reentry_model/levelset/markers.py`; test `tests/levelset/test_levelset_markers.py`.

**Interfaces:**
- Consumes: `cutcell.patches`, `levelset.redistance`, `transfer.locate/interpolate`, `criteria.principal_curvatures`.
- Produces: `Markers` with `x (n,3)`, `normal (n,3)`, `area (n,)`, `kappa (n,2)`, `skin_index (n,)` (row in `Skins`),
  `patch (n,)`, counters `splits, merges, trapped`. Functions: `Markers.create(mesh, phi, select)` (one marker per
  selected φ-patch at its centroid); `snap(mesh, phi)` (x ← x − φ ∇φ/|∇φ|, iterated once, §9.1); `carry(velocity, dt)`;
  `tessellate(mesh, phi)` (Ruling 2: neighbours from the patch graph within 1.5 target spacings; Voronoi in each
  tangent plane; returns `faces` as `(i, j, length, unit_normal_ij)` and sets `area`, scaled to sum to φ's patch area);
  `split(target_area)` (area above twice the target: two markers half a spacing apart along the longest Voronoi
  direction, each with the same profile and half the mass and area); `merge(target_spacing)` (closer than a third of the
  spacing: columns mixed layer by layer by mass and enthalpy, the hand-over rule); `trapped(mesh, phi, cell_size)`
  (no zero contour within half a cell) → indices; `assign(mesh, phi_patches)` (the patch whose cut cell contains each
  marker).
- [ ] **Step 1: Write the failing tests:** `test_rigid_motion_leaves_skins_unchanged` (a sphere's markers carried by a
  rigid translation and then a rotation, with φ moved the same way: every skin array bit-for-bit unchanged, areas within
  1e-10 relative); `test_uniform_stretching_thins_skins_and_draws_the_implied_metal` (a sphere growing from 25 to 26 mm
  radius with markers carried radially: each skin's top area grows by (26/25)², and after one draw at m_min each column
  holds its minimum and the drawn mass equals ρ times the implied volume within 1e-12); `test_split_and_merge_conserve_mass_and_enthalpy`
  (1e-14, counters right); `test_a_trapped_marker_returns_its_skin_exactly` (two spheres overlapping: markers in the
  overlap trapped, their mass and enthalpy returned as a deferred load, 1e-14); `test_areas_sum_to_the_surface_area`
  (sphere, 1e-12 of the patch area sum).
- [ ] **Step 2–4.** **Step 5: Commit.**

### Task 6: Runoff between markers, the pushes and landing

**Files:** Create `reentry_model/levelset/runoff.py`; test `tests/levelset/test_levelset_runoff.py`.

**Interfaces:**
- Consumes: `Skins` (Tasks 2–4), `Markers.tessellate` faces.
- Produces: `pushes(markers, skins, flow, decel_vec, sigma, phi_curvature) -> push (n,3)` (per unit volume, tangent: the
  pressure term −dp_w/ds from the shared `SurfaceFlowResult.dp_ds` along the meridian away from the stagnation point;
  deceleration ρ_l (g − (g·n)n); thickness ρ_l (g·n) ∇_s H; capillary −σ ∇_s κ with κ = ∇·(∇φ/|∇φ|) on φ (§20.6),
  only where the top flows); `thickness_step(...)` (the thickness and capillary terms implicit in the column heights, a
  two-point push per face, Newton on the piecewise-linear flux, at most 8 linearisations, tolerance 1e-3);
  `transport(skins, faces, dt) -> TransportBook` (linearly implicit upwind exchange (I + dt C) m' = m with
  c = (u_k · f)⁺ / A_donor across each face; the arrival rule: a flowing donor's span of fractions of its flowing depth
  joins the receiver's flowing layer between the same fractions; on a rigid top everything lands as a new top cell,
  crust metal always lands on top; anchored metal neither leaves nor receives; landed metal flows on at its arrival
  velocity; arriving momentum projected onto the receiver's tangent plane, the normal part booked to the wall);
  `landing_conduct(skins, landed, dt)` (the top 0.4 mm, at least two cells, adiabatic top and bottom);
  `flow_substeps(...)` (the count from the last sub-step's largest Courant number (target 6) and Rayleigh–Taylor
  growth (target 0.02), relaxed at most twofold, clipped to 1–40; 10 when nothing flowed before).
- [ ] **Step 1: Write the failing tests:** `test_runoff_onto_cold_metal_freezes_where_it_lands_and_the_books_close`
  (a flat strip of markers 1 mm apart, 0.2 mm of liquid at 950 K on the first 10, 700 K rigid beyond, shear 20 Pa, one
  10 ms heat sub-step of ten 1 ms flow sub-steps, shared AA7075: with landing conduction only the edge marker stays
  non-rigid and the deposit centroid is under 0.4 of the no-landing case's; mass and enthalpy 1e-13, momentum 1e-12 of
  the impulses); `test_landed_metal_flows_on_at_its_arrival_velocity` (per-marker ratio 5/6 within 1 %, the
  skyscrapers' case); `test_nothing_passes_through_anchored_metal`; `test_a_pushed_face_runs_downstream`
  (a flat patch at the skyscrapers' `PUSHES["flight"]`: dp/ds −3.424e5 Pa/m, shear 80.29 Pa, deceleration 29.77 m/s²,
  at 45°: books as above, downstream gain above 10 times upstream, face flux within 10° of the push).
- [ ] **Step 2–4.** **Step 5: Commit.**

### Task 7: Gas loads per marker

**Files:** Create `reentry_model/levelset/loads.py`; test `tests/levelset/test_levelset_loads.py`.

**Interfaces:**
- Consumes: shared `heating` (`PhysicsHeating.evaluate`, `SesamEquivalentHeating.evaluate`), `surface_flow.SurfaceFlow
  .evaluate(state, theta, radius, rho_liquid, T_wall)`, `nose.fit_sphere`, `nose.nose_radius`, `nose.load_angles`,
  `nose.ahead_of`; `shadow.shadowed`.
- Produces: `marker_angles(markers, v_hat) -> theta`; `leeward(mesh, phi, P, markers, v_hat, test="shadow"|"centroid",
  centroid=None) -> bool` (shadow: the stage 1 ray test on each marker's patch; centroid: behind the plane through the
  mass centroid, through `shared.nose`); `fit_nose(markers, centroid, v_hat, size_feedback="current",
  fallback) -> NoseFit(fitted, bound, used, n_points)` (ledger row 1: the windward half of φ — points ahead of the plane
  through the mass centroid normal to the flight direction — fitted by `fit_sphere`, unweighted; bound = 1.67 times
  the largest transverse distance of those points; `used` by `shared.nose.nose_radius`); `GasLoads(heating, surface_flow)
  .evaluate(state, markers, T_top, theta, nose) -> (q_conv per marker, SurfaceFlowResult)`.
- [ ] **Step 1: Write the failing tests:** `test_a_sphere_fits_its_own_radius` (25 mm sphere: `fitted` within 1e-6,
  `used` = fitted); `test_a_flat_face_is_capped_at_1_67_transverse_radii`; `test_loads_are_the_shared_laws` (on a
  sphere at a fixed state, `q_conv` and `tau` equal the shared functions called directly with the same angles, radius
  and wall temperatures, bit-for-bit); `test_leeward_centroid_matches_shared_rule`; `test_on_a_sphere_shadow_and_facing_agree`.
- [ ] **Step 2–4.** **Step 5: Commit.**

### Task 8: Spraying per marker, and the droplet records

**Files:** Create `reentry_model/levelset/spraying.py`; test `tests/levelset/test_levelset_spraying.py`.

**Interfaces:**
- Consumes: shared `spray.SprayModel.evaluate`, `spray_layer`, `region_extent`, `slurry_properties`, `SLURRY_SIGMA`,
  `source_block`, `histogram`, `droplets.SKIN_SUMMARY_COLUMNS`; `shared.lubrication.lubrication` for v_s.
- Produces: `spray_plan(skins, markers, flow, state, model, windward, w_n, dt, radius, nose_radius, slurry=False,
  b_mode="skin"|"phi_psi", depth_phi_psi=None) -> SprayPlan` (once a macro step, as the skyscrapers' `girin_plan`:
  b and the regime depth = the skin's flowing depth `depth_above(T_rigid)`, or φ to ψ on a top-zone skin (§20.6);
  δ_m from `spray_layer` (slip band through it); v_s from `lubrication(tau, G, F, δ_m, μ_l, b_layer=F)[0]`; cap depth
  the liquid depth, or the flowing depth with slurry; the RT pool the liquid depth; extent from `region_extent` on the
  markers' face graph); `SprayPlan.release(skins, dt)` (the cascade: liquid cells from the top first, the slurry beneath
  at its own rate with `slurry`, never rigid metal; each cell leaves at its own enthalpy and velocity) `-> ReleaseBook`;
  `SprayPlan.records(t, state, ...)` → rows by the shared `source_block`; `step_summary() -> dict` with every
  `SKIN_SUMMARY_COLUMNS` key.
- **No renewal limit** (ruling of 2026-11-11, on a question passed from Asha's forked session). Design §9.5 caps
  liquid deeper than δ_m below the top at ρ_l δ_m A per renewal time δ_m²/ν_l; §20.6 says spraying acts "as in §9.5".
  The skyscrapers code at the parity commit ae20dd3 applies no renewal limit (searched: `skin.py`, `skinbody.py`,
  `shared/spray.py`, `film.py`); the limit belongs to the 2026-10-07 two-zone skin design, where it stood in for the
  momentum diffusion of an unresolved deep layer. The layered skins resolve the velocity through the depth cell by cell,
  so the release is capped only as the skyscrapers cap it (the cap depth and the cascade). §9.5's renewal clause is
  therefore not built; Task 12 records this as an amendment of §9.5 in the design.
- [ ] **Step 1: Write the failing tests:** `test_marker_release_equals_the_shared_model` (for a set of markers with
  given depths and a fixed gas state, the plan's `mdot`, `r`, `branch` equal `SprayModel.evaluate` called directly with
  the same arrays, bit-for-bit); `test_the_cascade_takes_liquid_first_and_never_rigid_metal`;
  `test_the_cascade_matches_the_skyscrapers_release_on_one_column` (D4); `test_release_is_not_capped_by_a_renewal_time`
  (a 2 mm liquid column on rigid metal whose shared rate over one 0.5 s macro step releases more than ρ_l δ_m A per
  δ_m²/ν_l releases the shared rate's mass, within 1e-12); `test_release_books_close`
  (1e-14); `test_droplet_records_go_through_the_shared_writer` (a fake history of three steps and the plan's rows:
  `shared.droplets.write_particles` writes `particles.npz`, `particles_summary.csv`, `size_distribution.csv`, and the
  released mass in them equals the books'); `test_girin_cases_through_the_level_set_path` (the 1994 Table 2 case's
  λ* within 2 % and τ* within 1 %, fed through `spray_plan` on a single flat marker with the case's depth).
- [ ] **Step 2–4.** **Step 5: Commit.**

### Task 9: The skinned body — one macro step coupled to the cut-cell solve

**Files:** Create `reentry_model/levelset/skinbody.py` (`LevelSetSkinBody`); test `tests/levelset/test_levelset_skinbody.py`.

**Interfaces:**
- Consumes: everything above; `CutThermal` (stage 2), `surfaces.outer_surface`, `surfaces.molten_rigid_surface`,
  `criteria.make_target`, `boxmesh.build`, `levelset.correct_volumes/pieces`.
- Produces: `LevelSetSkinBody(material, base, beta0, T0, settings, loads=None, spray_model=None, leeward="shadow",
  emissivity, fine_margin=4e-3, fine_size=0.7e-3, cell_floor=0.25e-3, remove="spray"|"feed")`;
  `advance(dt, state=None, decel_vec=None, q_prescribed=None) -> StepBook` (one macro step, design §10.1 with §20:
  loads once per marker; Girin's plan once; parts per Ruling 7; in each part: the trial pass gives each marker's (a, b),
  each β-patch's law is q0 = Σ a / A_patch, q1 = Σ b / A_patch over its markers (Ruling: per area, the skyscrapers'
  per-column law converted), `CutThermal.step` with the pending load, the real pass with T_b read at each marker's base
  from the 3D field, mismatch E_base − (a + b T_b)dt booked as next part's pending, repeat by the damped Aitken rule
  (ω₀ 0.5, ω in [0.02, 1], stop at 1e-3 relative or 0.05 K or 12 passes); then draws, hand-downs and trapped returns
  move β (β − δ extended along the normal, reset, `correct_volumes` per piece to the exact volume target) with
  `CutThermal.move`, its enthalpy difference from what the skins took booked as pending; φ = `outer_surface(β, H)`;
  ψ = `molten_rigid_surface(φ, T, material.T_rigid)`; markers snapped, tessellated, split, merged; the mesh rebuilt
  on the fine region (ψ > −d) and curvature with `CutThermal.remesh`). `StepBook`: mass and enthalpy in, out (sprayed,
  removed, radiated), stored (volume + skins), the residuals, `momentum_to_volume`, `momentum_to_wall`, the history
  columns of §12.2 (`n_markers`, `marker_splits`, `marker_merges`, `markers_trapped`, `interface_mismatch_J`,
  `volume_correction_max`, `n_pieces`, `frontal_area_m2`) and `SKIN_SUMMARY_COLUMNS`.
- [ ] **Step 1: Write the failing tests:** `test_a_macro_step_closes_mass_and_enthalpy` (a 10 mm cube of the AA7075
  Scheil material at 850 K heated 2 MW/m² on one face for 0.5 s, dev cells: mass 1e-12, enthalpy 1e-8 relative);
  `test_a_cold_body_is_unchanged_by_the_skins` (no heat: T, β, the skins unchanged); `test_pending_load_reaches_the_volume`;
  `test_markers_follow_the_recession` (after a step every marker within 1e-9 m of φ = 0).
- [ ] **Step 2–4.** **Step 5: Commit.**

### Task 10: Two modes — the thin test per marker and exact mode changes

**Files:** Create `reentry_model/levelset/modes.py`; modify `levelset/skinbody.py`; test `tests/levelset/test_levelset_modes.py`.

**Interfaces:**
- Consumes: `thin.thin_candidate`, `thin.ModeState` (now keyed by marker), `thin.patch_measures` at the markers.
- Produces: `ModeSettings(c_kappa=0.1, c_slope=0.2, hold=3, min_markers=7)`; `to_thick(skins, i) -> Handover` (the
  full-depth skin's cells below its fine zone go to the volume: β steps out by their volume; enthalpy exact through the
  pending load; tangential momentum to `momentum_to_volume`); `to_thin(skins, i, depth) -> Handover` (the volume melt
  beneath the marker, down to the stub below ψ, into the skin, filled from the 3D field along the normal; β steps in);
  top-zone skins keep 0.4 mm (draw and return at their base), their base velocity is zero (flow off, Ruling 3), and
  their b for spraying is φ to ψ (D2).
- [ ] **Step 1: Write the failing tests:** `test_mode_changes_conserve_mass_enthalpy_and_tangential_momentum` (both
  directions on a column of melt over rigid metal: mass and enthalpy 1e-14, tangential momentum exact, the normal part
  booked); `test_a_new_marker_may_change_mode_only_after_the_hold`; `test_a_small_group_does_not_change_mode`
  (six markers wanting thick stay thin); `test_thin_and_thick_neighbours_exchange_only_their_top_zone` (Ruling 3).
- [ ] **Step 2–4.** **Step 5: Commit.**

### Task 11: Verification cases and cost

**Files:** Create `analysis/levelset/stage3_cases.py` (subcommands `plate-cold`, `plate-hot`, `box`, `shoulder`,
`cost`), test `tests/levelset/test_levelset_stage3_cases.py` (marked `reference`).

- [ ] **Step 1: The plates at 0.5 s macro steps (Decision A of stage 2).** The level-2 plate's centre column (8 mm deep,
  insulated walls on grid lines, shared AA7075 Scheil, no radiation, `swell=False`, flow off, `remove="feed"`, skins
  removing liquid continuously within their sub-steps), dev and production cells: cold (700 K, 6 MW/m²) recession
  within 5 % of `REFERENCE_1D = {1.0: 0.789, 1.5: 2.712, 2.0: 5.219}` mm and of `REFERENCE_1D_LIQUIDUS = {1.0: 0.837,
  1.5: 2.772, 2.0: 5.298}` mm (the numbers of `analysis/skin_plate.py` 86–95, copied as data with that reference;
  the script itself imports skins and is not imported); hot (880 K, 2 MW/m²) recession within 5 % at 1.5 and 2.0 s
  (1.0 s reported) and non-rigid depth within 10 % of `data/reference_values/plate_hot_1d.json`; mass 1e-12, energy
  1e-8. Tests `test_cold_plate_at_half_second_steps`, `test_hot_plate_at_half_second_steps`.
- [ ] **Step 2: The q6p box** (Ruling 5): 8 × 10 × 10 mm, 700 K, 6 MW/m² on the 10 × 10 face, 2 s, dev cells, flow
  off with the skins' own flow on: recession at the face centre within 5 % of the cold reference at 1.0, 1.5, 2.0 s,
  books as above, no marker trapped, edges reported. Test `test_q6p_box`.
- [ ] **Step 3: Runoff over a shoulder.** A 25 mm sphere segment at the pushed plate's loads (`PUSHES["flight"]` at the
  45° station) with liquid fed on the windward cap and cold metal behind: the deposit freezes on the cold metal, books
  close, and the deposit's position is reported against the skyscrapers' pushed-plate result. Test `test_shoulder_runoff_freezes_a_deposit`.
- [ ] **Step 4: Cost.** On the 50 mm sphere at a melting state (dev and production cells): time and memory per macro
  step split into loads, skins (trial, real, repeats), cut-cell solve, surfaces (φ, ψ), markers, remesh; the cost of
  cells at the 0.175 mm floor on a rim of radius 0.5 mm (design §20.13); whether ψ, the shadow test and the remesh can
  run every second step (the change in a cold-plate recession when they do). Written to `output/levelset/stage3_cost_*.json`.
- [ ] **Step 5: Run** (mains power; if the estimate exceeds two hours, write a Mac Studio runbook instead and stop for
  Asha). **Commit** the script, the test and its recorded numbers.

### Task 12: Documentation and the ledger

**Files:** Modify the repository `README.md` (the level-set section: stage 3), `docs/model_assumptions.md` (new §12,
stage 3: each rule with its test or marked assumed, the rulings, the measured cost), the design (§20.11 item 3: built,
the measured results; §20.13: the cost; §14's prototype line amended), `docs/superpowers/plans/levelset/parity-ledger.md`
(rows pinned at ae20dd3; rows 1, 4, 6 by Asha's answers to D1–D3; rows 3 and 9 closed with their tests),
`Reentry Model/CLAUDE.md` (the level set is built through stage 3; the shared-physics rule).
- [ ] **Step 1: Write the documents**, every number with the test or run that measured it, once-only measurements said so.
- [ ] **Step 2: Run** `tests/levelset`, `tests/test_reentry_model_shared.py`, `tests/test_reentry_model_cli.py` → all pass.
- [ ] **Step 3: Commit.**
