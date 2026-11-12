# Girin's gradient instability solved on the skins' own melt profile — design

**Date:** 2026-11-11. **Status:** draft for Asha's review; nothing built.
**Lands after:** the combined refactor (stage 0 + S, `docs/superpowers/plans/levelset/00s-combined-refactor-shared-package.md`)
is on `main`. Built on a short branch off `main` (Asha, 2026-11-11: no prototype folder).
**Prerequisite, gating stage 3 only (decision 7):** the single-column check of the gas wall shear against the skins'
surface speeds (a runny layer over slurry), which is on the post-refactor carry-over list ("MELT SURFACE SPEED vs GAS
WALL SHEAR"). Stages 1 and 2 read no skin speeds and start at once; stage 3's committed fixture records the speeds, so
the check comes before it, and if it finds the speeds wrong they are fixed first. A rough read of the d050 state
(2026-11-11, one state) already points to consistency: the gas shear on the 622 thick columns is 61–142 Pa (median 101),
and that shear pushed through the cells' own viscosities at rest, as a constant-stress layer, gives top speeds of
0.07–0.17 m/s — a lower bound, since viscosity at rest is the highest — against the skins' 0.18–0.49 m/s, the rest
plausibly from the pressure gradient (driving term median about 28 kPa/m). The side session's mismatch (0.4–58 m/s
predicted) assumed the top layer was liquid at 1.3 mPa s; it is the bridge, at 5 to 600 times that.

**Sources.** The side session's report `reports/Kinked profile film dispersion 2.md` (repository root) and its notes in
`research_notes/Kinked profile film dispersion/`, chiefly `fulltext_girin_papers.md` (the derivation), `fulltext_viscous_papers_read.md`
(Ghosh 2017 check, viscous modes) and `fulltext_usha_hooper1985.md` (the floor test). Its check scripts were in a session
scratchpad (`d2.py`, `ext.py`, `ghosh_check.py`, `floor.py`) and are not in the repository; the builder re-derives and
re-tests everything here rather than relying on them. Asha's profile plot: d050, hand-down off, t = 178.5 s.

**Asha's decisions (2026-11-11), each the recommended option:**

1. The kink and the floor are found **from the speed profile**: straight segments fitted to each column's cell speeds,
   breakpoints where the shear rate drops sharply, the floor test run at each breakpoint from the top down, and the
   shallowest floor-like breakpoint becomes the wall (otherwise the rigid front).
2. Where the skins resolve the flow, the gradient mode is solved with **the skins' own speeds and depths**. Girin's
   theoretical delta_m and V_s stay as a recorded check and as the fallback where a column has no resolved flow.
3. In fluid-like columns the slurry leaves **with the wave, in one event**: everything above the depth where the unstable
   wave's amplitude falls to 1/e of its surface value, capped at the flowing depth.
4. An offline **viscous (Orr–Sommerfeld) check** that calibrates the floor threshold is in scope, as the last stage,
   analysis only.

After the builder check (same day):

5. The release is **heat-limited**: cells at least as warm as the cell at the release depth (or just above the wall) at
   the macro step's start leave at the mode's rate; colder cells wait for the next macro step (section 6.6).
6. The profile path **keeps today's regime limit**: its shear modes act only where delta_m is finite — Girin's continuum
   closure, or the slip band with `--slip-band on` (section 6.5). Lifting it is a later, separate decision.
7. The formal **gas-shear check** stays, before stage 3; stages 1 and 2 start at once.

Items marked **[my call]** were settled by me while writing and are listed again in section 11 for Asha to check.

---

## 1. Why

The spray model's gradient mode today solves Girin's (2017, A&A 606, A63) Eq. 1 — a single linear melt boundary layer
of depth delta_m over stagnant, infinitely deep liquid, with a free surface on top — using Girin's *theoretical* layer
depth delta_m and surface speed V_s (his Eq. 2 conjugate-layer closure through `film.lubrication`, or the slip band's
momentum balance). The layered skins compute a different profile. In the "thick" windward columns of the 50 mm flight
(sheared depth about 0.15 to 0.30 of the flowing depth b) almost all the shear sits in a thin runny layer on top, below
which slurry hundreds to thousands of times more viscous shears gently down to the rigid front. Further columns have a
speed peak below the surface ("plug" columns).

The side session showed four things (its report, sections 1–4):

- Girin's Eq. 1 is exactly an inviscid **broken-line** model (straight velocity segments, all the physics in matching
  conditions at the corners and at the free surface, gas density neglected). Its exact threshold is a surface Weber
  number of 3.00, not 3.08; the fastest wave at large Weber number is wavenumber times layer depth 1.226 with growth
  rate 0.247 surface speeds per layer depth (Girin: 1.225 and 0.24).
- The method extends to slurry and a wall, and the extension agrees with Ghosh et al. (2017, Meccanica), who solved the
  same geometry independently, to between 1e-16 and 1e-12.
- Inviscidly, the slurry's slow shear barely matters, but a wall close under the runny layer suppresses the gradient
  mode completely: a single straight film on a wall is stable at every wavelength (Ghosh, Eq. 20).
- A slurry 500 to 4000 times more viscous than the liquid behaves as a **floor** (a wall) on the wave's time scale when
  its viscous diffusion outruns the wave's advection. That holds throughout the slip regime by a factor of at least 2.5
  and fails only for films of 150 micrometres or more under hard deceleration. This is a scaling argument, not a solved
  viscous problem.

Asha's context: the flights seen so far are in the slip regime, before the hard deceleration that drives the shear up;
the melt does spray later, so the kink has to be treated correctly for that phase. The goal is therefore to solve the
gradient mode on the profile the skins actually compute, treat the slurry as a floor where it is one, and let it move
with the wave where it is not.

## 2. What the flights show, and what that changes (measured 2026-11-11, one state)

Probe of the saved state `d050__handdown-off.state.pkl` at t = 178.5 s (subplan20 worktree,
`prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/`), on the 622 windward
columns (angle from the stagnation point below 60 degrees) whose sheared depth is 0.12 to 0.35 of the flowing depth.
This is one state of one flight; the numbers are medians with 10th–90th percentiles.

- **The runny layer is not liquid.** Every one of the 622 tops is below the 908 K liquidus (top cell 907.1–907.6 K), so
  the depth to the liquidus is zero in all of them. The steep layer is the rheology's **bridge**, the 906.4–908 K band
  in which `rheology.viscosity` interpolates log-linearly, without data, from Li et al.'s slurry law at 10 % solid down to
  the liquid. Its depth (to the 906.4 K isotherm) is 55 micrometres (24–67); the flowing depth b is 830 micrometres
  (770–1980). Across the bridge the viscosity rises from about 5 to about 600 times the liquid's (Asha's plot), so the
  runny layer's own speed profile is strongly curved, not straight.
- **The slurry moves faster than surface shear alone would make it.** The speed at the bridge's base is 0.21 of the top
  speed (0.11–0.68). A constant-stress (Couette) layer with a thousand-fold viscosity jump would give about 0.01, so
  something else drives the slurry — presumably the pressure gradient and the deceleration term in the skins' momentum
  step. The broken-line profile takes the speed at the kink from the skins and needs no model for this.
- **Top speeds** are 0.23 m/s (0.18–0.49).
- **The kink rule of section 6 lands at the bridge's base** (builder check, same state): a fit by section 6.4's rules
  needs 3 to 5 breakpoints per column; 94 % of the columns have a breakpoint with at least a five-fold drop in shear
  rate, at a median 0.94 times the depth of the 906.4 K isotherm. Of the 957 wet windward columns, 14 are plug columns
  and none has a top cell at rest.

Consequences for the design, beyond the side session's proposal:

- A wall "at the liquidus" would put b at zero in exactly these columns. Asha's decision 1 replaces it with breakpoints
  from the speed profile.
- **Two existing blockers.** The cascade's liquid run (`skin.cascade_mass`) and the wave-fit region
  (`skinbody.girin_plan`, `region_extent(edges, L > 0, ...)`) are both defined by cells at or above the liquidus. In these
  columns both are empty, so today the liquid's mode could release nothing there even if it were unstable. Under the new
  setting both are redefined by the wall (section 6.6).
- The bridge's runny layer is viscous, not inviscid: with its top viscosity of about 5 times the liquid's, the top
  20 micrometres have a Reynolds number (velocity jump times depth over kinematic viscosity) of order 1. The broken-line
  solver is inviscid, so the viscous gate of section 6.5 and the viscous check of stage 6 matter here.
- The kink's depth, and so every input to the floor test, rests on the bridge's unmeasured 906.4 K end point. This
  design does not change the bridge; it records the kink's depth and temperature so the dependence stays visible.

## 3. The physics the design rests on

The derivations below are the side session's and are not published; each is checked against Girin's printed equation
or against Ghosh (2017). The builder re-derives them, and the tests of stage 1 catch sign errors (a sign error in the
kinematic condition was caught exactly this way in the side session).

**Broken-line model.** Coordinates: y upward, free surface at y = 0, wall at y = −H. The base speed U(y) is continuous
and piecewise linear, with breakpoints at depths d_1 < d_2 < ... < d_N (y_j = −d_j), segment slopes S_0 (top segment)
to S_N (bottom segment), segment densities rho_0 to rho_N, and U(−H) = U_w, the wall's own speed (zero at the rigid
front; the floor's speed when the wall is a floor). A disturbance stream function phi(y) exp[ik(x − ct)], with c the
complex wave speed, obeys the Rayleigh equation; because U'' = 0 inside every segment, phi'' = k² phi there, so each
segment's solution is a pair of exponentials. Conditions:

- wall: phi(−H) = 0 (no flow through it);
- at each breakpoint y_j: phi continuous, and pressure continuous,
  rho_above [(U_j − c) phi'_above − S_above phi_j] = rho_below [(U_j − c) phi'_below − S_below phi_j],
  which for equal densities is Girin's (U − c)[phi'] = [U'] phi;
- free surface, with the surface displacement eta: kinematic, (U_0 − c) eta = −phi(0); normal stress with surface
  tension sigma and the gas pressure neglected, rho_0 [−(U_0 − c) phi'(0) + S_0 phi(0)] = sigma k² eta.

With two coefficients per segment plus eta there are 2N + 3 unknowns and 2N + 3 conditions, each linear in c: a
generalised linear eigenvalue problem (A − c B) x = 0 with N + 2 finite eigenvalues. Girin's profile (one breakpoint,
no wall: the bottom segment keeps only its decaying exponential) gives 3, his cubic. The growth rate of a wave is
k Im(c); the fastest wave maximises it over k. The eigenvector gives phi(y), the wave's depth structure.

**The floor test.** A layer of depth d above a breakpoint, with a velocity jump DeltaU across it relative to the
material below, drives a wave of wavenumber about k = 1.225/d (Girin's fastest). The material below acts as a no-slip
floor on that wave's time scale when its viscous diffusion outruns the wave's advection,

  floor Reynolds number Re_f = DeltaU / (nu_below k) = DeltaU d / (1.225 nu_below)   is small,

with nu_below the kinematic viscosity (viscosity over density) just below the breakpoint. This is the stricter of the
side session's two versions (the other, comparing diffusion with the growth rate, is about 5 times weaker). Surface
speeds at which Re_f = 1, for slurry 500 / 4000 times the liquid's viscosity (`fulltext_usha_hooper1985.md`, table in
section 3): film 5 micrometres, 66 / 530 m/s; 50 micrometres, 6.6 / 53 m/s; 150 micrometres, 2.2 / 18 m/s;
450 micrometres, 0.74 / 5.9 m/s. It assumes a Newtonian slurry; a shear-thinning or yielding slurry would weaken the
floor. Stage 6 calibrates the threshold.

## 4. Goals, non-goals, anti-goals

**Goals**

1. A shared, inviscid broken-line dispersion solver for any piecewise-linear profile with a free surface and a wall,
   which reproduces Girin's Eq. 1 exactly in his limit.
2. The floor test, and the wall it places, per column.
3. The skins' cell speeds reduced to that profile, and the gradient mode solved on it with the skins' speeds.
4. A choice between the gradient mode and Girin & Kopyt's (1994) gas-driven mode by which grows faster.
5. In fluid-like columns, one release taking liquid and slurry together down to the wave's depth.
6. Records and history columns that make the wall, the kink and the floor test visible.
7. The offline viscous check that calibrates the floor threshold (stage 6, analysis only).

**Non-goals**

- The element model (`body.MeltingBody`) is unchanged: it has no resolved profile and keeps the table. **[my call]**
- The level-set architecture: the solver is shared physics so it can be called there later, but wiring it in is not
  here.
- Deceleration inside the solver. The front-surface Rayleigh–Taylor branch stays a separate branch, unchanged (it still
  reads the depth to the liquidus — see section 10). Adding the term later is one term in the surface condition (as
  Ghosh's gravity term). **[my call]**
- The gas side (density ratio neglected, as Girin), surface tension at kinks (unknown; zero), the oxide skin,
  non-Newtonian slurry, nonlinear crest stripping (Ishii and Grolmes type), and any change to the rheology bridge.
- The exposed slurry's own mode (`skinbody.slurry_modes`, `--slurry-spray`) keeps its physics; on the `profile` path
  only the cell it reads changes (section 6.6). **[my call]**

**Anti-goals** (failures even if it runs)

- Any change of result under the existing setting. With the new setting off, runs are bit-identical
  (`analysis/skyscrapers_port/same_results.py` cases).
- The bridge's construction silently deciding where spray comes from: the kink's depth, temperature and the floor
  Reynolds number are recorded every macro step.
- Girin's large-Weber growth rate applied without surface tension: the solver always includes it.
- A gradient-mode release from a straight film resting on a floor (Ghosh: stable at every wavelength).

## 5. Constraints

- Shared modules (`reentry_model/shared/`) import no skyscrapers module (refactor rule, pinned by its independence
  test). The solver and the floor test take plain arrays; turning skin cells into a profile lives in `skinbody.py`.
- Evaluated once per macro step at its start, like the rest of the spray plan (spec 8), on windward wet columns only.
- Cost: measured on the d050 state; the added time per macro step should stay under about 5 % of a macro step's wall
  time. **[my call]** Columns whose layer above the wall is one straight segment skip the eigen-solve (stable by Ghosh).
- A run saved before the setting existed resumes with today's behaviour (the class-default pattern of `spray_gates`).
- The run name encodes the setting when it is not the default (repository convention).
- `--dispersion profile` requires `--skin-flow on`, and the CLI refuses it otherwise with exit code 2, as it does for
  `--slurry-spray on`. **[my call]**

## 6. Design

### 6.1 The setting

A new run and resume option `--dispersion table|profile`, default `table` (today's behaviour, bit for bit). `profile`
switches the skins' liquid mode to sections 6.2–6.7. The default stays `table` until Asha has seen the comparison
flights of stage 5 and decides. **[my call]**

### 6.2 The solver (`reentry_model/shared/dispersion.py`)

Inputs: breakpoint depths, speeds at the surface, at each breakpoint and at the wall, densities per segment, surface
tension, and a wavenumber (or a range). Outputs per profile: the fastest wave's wavenumber k_f, growth rate
Im(omega_f) = k_f Im(c), phase speed Re(c), the eigenfunction phi(y) normalised at the surface, and the release depth.
The release depth is read from the wave's vertical **displacement**, |phi(y) / (U(y) − c)|, not from phi (which is the
vertical velocity): it is the depth, below the displacement's peak, where the displacement first falls to 1/e of its
surface value, capped at the wall. On Girin's profile at his fastest wave the displacement peaks at 1.6 times its
surface value 0.67 layer depths down and falls to 1/e at 1.82 layer depths, against 1.29 for phi (reviewer's
calculation, 2026-11-11). **[my call]**

Requirements on the numerics (the method is the builder's):

- Each segment's exponentials are written in a basis normalised inside that segment, so that no matrix entry grows
  with k times the segment's depth; correct at kH up to at least 100.
- Infinite eigenvalues of the generalised problem are discarded; among the finite ones the most unstable is taken.
- The fastest wave is found on a log-spaced grid that brackets every segment's scale (at least k d_1 from 0.02 to 5,
  and down to kH = 0.02), then refined. Columns are solved in batches; a batched standard eigenproblem after eliminating
  the conditions free of c is one way.
- The existing table code stays as it is for the `table` path and for the element model.

### 6.3 The floor test and the wall (shared, beside the solver)

Given the profile, each breakpoint's kinematic viscosity just below it, and a threshold:

- A breakpoint is a **kink** (a wall candidate) when the shear rate below it is at most `KINK_DROP` times the shear
  rate above it (default 0.2: a five-fold drop). Other breakpoints only describe curvature. **[my call]**
- Kinks are tested from the top down. At kink j: d = its depth, DeltaU = |U(0) − U(y_j)| (the jump relative to the
  candidate floor, which moves), nu_below = the kinematic viscosity of the cell immediately below it (the least viscous
  part of the would-be floor), Re_f = DeltaU d / (1.225 nu_below). **[my call: relative jump, cell just below]**
- The first kink with Re_f < `FLOOR_RE` is the wall (**floor-like** column); the profile below it is dropped and the
  wall moves at U(y_j). If none qualifies, the wall is the rigid front at the flowing depth (**fluid-like** column).
- `FLOOR_RE` defaults to 1.0, the crossover of section 3, labelled assumed until stage 6. **[my call]**

### 6.4 From the skins to a profile (`skinbody.py`)

- Speeds: each flowing cell's velocity component along the top cell's direction, at its centre; the surface speed
  extrapolated linearly from the top two cell centres; zero at the rigid front (`depth_above(T_rigid)`). **[my call]**
- A continuous piecewise-linear fit through the surface speed and the zero at the rigid front, breakpoints at cell
  centres, the fewest that keep every cell within `FIT_TOL` (default 3 % of the top speed) of the fit, at most `N_MAX`
  (default 6) breakpoints. **[my call]**
- Per segment: the mass-weighted density of its cells. Per kink: the cell viscosity just below it, from
  `skin.cell_viscosity` at the cell's shear rate as the skins use it.
- **Fallback (decision 2):** a wet windward column whose top cell does not move (no resolved shear) takes today's
  theoretical path for that column.

### 6.5 Choosing the mode

Per wet windward column, after the wall is placed:

- **Regime limit (decision 6):** both shear modes below act only where `spray.spray_layer` gives a finite delta_m —
  Girin's continuum closure, or the slip band's momentum balance with `--slip-band on` — exactly where they may act
  today, and the continuum gate of `--spray-gates 2026-11-10` applies as today. delta_m decides only *whether* a column
  may spray; the profile decides *how*.
- **Gradient mode** on the profile above the wall (skipped, and stable, when that profile is one straight segment).
  It counts as unstable only if its growth exceeds the viscous damping of its own wave,
  Im(omega_f) > 2 nu_eff k_f², with nu_eff the layer's Couette viscosity (depth over the sum of each cell's thickness
  divided by its viscosity, over density). No critical Weber number is applied; the solver's own threshold replaces it.
  **[my call]**
- **Gas-driven mode:** Girin & Kopyt's thin branch as today (wavelength from the edge state, growth time from their
  Eq. 12, the continuum gate of `spray_gates`), with the film depth b = the depth to the wall and the viscous-film check
  (`spray.viscous_film`) on that depth with nu_eff in place of the liquid's viscosity. **[my call: nu_eff]**
- The mode that grows faster is taken. In a floor-like column with a straight layer that is always the gas-driven mode
  (the side session's proposal (b)); in a curved or fluid-like column either can win.
- The rarefied branch and the Rayleigh–Taylor branch are unchanged.

### 6.6 Release

- **Gradient mode:** wavelength lambda = 2 pi / k_f, droplet radius r = k_r lambda, growth time k_t / Im(omega_f), rate
  per area Girin's torus form rho pi r² / (lambda t) — today's formulas with delta_m / (V_s Im Omega) replaced by
  1 / Im(omega_f), so the same k_r = 0.17 and k_t = 1.1.
- **Release run (decisions 3 and 5, and the blockers of section 2):** at the macro step's start the plan fixes a
  depth — for the gradient mode the release depth of section 6.2, capped at the wall; for the gas-driven mode the wall —
  and a **threshold temperature**, that of the cell at that depth (the last cell above the wall for the gas-driven
  mode). During the step, cells from the top down leave at the mode's rate while they are at least as warm as the
  threshold, liquid and slurry together, whatever their temperature relative to the liquidus; the first colder cell
  ends the column's release until the next macro step's evaluation. This is the slurry cascade's existing rule
  (`skin.Cascade.T_slurry`) applied to the mode's run: the release is limited by heating, not by the number of macro
  steps. (Why it matters: the torus rate amounts to the surface receding at about 0.1 times the surface speed, about
  29 micrometres per millisecond at 0.27 m/s, so any run defined by depth alone is used up within milliseconds; a fixed
  budget per step would make the sprayed total scale with the step count, and a window following the top would strip
  all non-rigid metal in hundredths of a second — reviewer's estimates, 2026-11-11.) The one-crest cap
  ((4/3) pi r³ at most lambda² h) uses h = the run's depth at the step's start.
- Below the run, exposed slurry cascades by its own mode with `--slurry-spray on`. On the `profile` path
  `slurry_modes` reads the first flowing cell **below the run** (its layer from that cell to the rigid front, its
  temperature as the cascade's threshold), not the first cell below the liquid run: in the d050 thick columns there is
  no liquid run, so today's choice would read the top cell, inside the run, and set a threshold (about 907.4 K) that
  every cell below the run is colder than. **[my call]**
- The wave-fit region (`wave_fits`) is the contiguous region of wet windward columns with a layer above the wall at least
  `spray.B_MIN` deep, not the liquid-pool region. **[my call]**
- Droplet properties for the records: the released cells' mass-weighted density and viscosity (at their shear rate);
  surface tension 0.80 N/m, which is both the liquid's and `SLURRY_SIGMA`. **[my call]**

### 6.7 What happens to delta_m, b, the thick/thin gate and the Weber threshold

- **b** is no longer always the flowing depth: it is the depth to the wall — the kink's depth in floor-like columns, the
  flowing depth in fluid-like ones.
- **The thick/thin gate** (layer against delta_m) is not used on the `profile` path: the wall enters the solver
  explicitly, so its stabilising effect comes out continuously. This reverses, on that path, the 2026-10-06 rule that a
  film on slurry deeper than delta_m takes the thick branch.
- **delta_m** (Girin's Eq. 2, or the slip band's momentum balance) stops being a depth or speed scale where the skins
  resolve the profile. It still decides whether a column may shear-spray at all (decision 6), is recorded beside the
  resolved depths as a check on the skins' momentum, and is the input wherever the fallback applies. The flow-weighted `sheared_depth` is not its replacement, because it
  averages across the kink.
- **`--we-critical`** (default 3.08) applies only to the `table` path. The record still carries a surface Weber number
  for comparison: density times the square of the top segment's velocity jump times its depth, over surface tension.

### 6.8 Records and history

- A new branch code for the gradient mode solved on the profile (5), so its rows are told apart from today's thick
  branch (0). The gas-driven mode keeps code 1. **[my call]**
- On the `profile` path only, the source table gains columns appended at the end (existing columns keep their
  meaning): depth to the wall, wall kind (floor or front), floor Reynolds number at the wall, depth and temperature of
  the first kink, release depth, the solved growth rate, and Girin's theoretical delta_m and V_s.
- History, on the `profile` path only, per macro step: the share of wet windward columns that are floor-like, the count
  of fluid-like ones, medians of the depth to the wall, the floor Reynolds number and the recorded Weber number, and the
  mass sprayed by each mode.
- `table` runs, and the element model, write exactly today's columns and settings keys, so `same_results.py` passes
  unchanged (it compares headers and keys exactly) and `spray.SOURCE_COLUMNS` stays the element model's. The spool,
  the coupled writer and every reader handle both layouts, as the spreading history columns already exist only with a
  spread limit. A state saved before the setting existed has 27-column records, which are padded on resume when it is
  switched to `profile` (as when the size bins grew from 40 to 50). **[my call]**

## 7. Stages and verification

1. **Solver** (shared). Tests: Girin's cubic roots reproduced to 1e-10 at several (wavenumber, Weber number) pairs; the
   cached table's fastest wave (1.226, 0.247 at Weber number 1e4) recovered within its grid spacing; the threshold
   for Girin's profile between 3.00 and 3.01; the side session's closed form for slurry over a wall
   ((1 − c)(Delta tanh Delta + q) − [S_1 + Delta²/(We_s (1 − c))](1 + (q/Delta) tanh Delta) = 0, its notation)
   reproduced to 1e-10 for several wall depths, shear ratios and density ratios; Ghosh's limits (a single straight film
   on a wall stable at all wavenumbers; more viscous fluid on top stable; peak growth 0.229 at slope ratio 1.5 and
   interface at 0.4 of the depth with no surface tension); a wall approaching the kink drives growth to zero; correct at
   kH = 100.
2. **Floor test** (shared). Tests: section 3's crossover speeds reproduced; a kink below a five-fold shear drop is a
   candidate and a gentler one is not; the wall is the shallowest floor-like kink.
3. **Profile from the skins** (`skinbody.py`). A committed fixture of a few dozen columns extracted from the d050 state
   at 178.5 s (the state itself is git-ignored and on this laptop only). Tests: Girin's linear profile on synthetic
   cells gives one breakpoint at its depth; a skin moving as one piece gives none; on the fixture, report where the wall
   lands relative to the 906.4 K depth (a measurement, not a pass/fail).
4. **Integration** under `--dispersion profile`: mode choice, release run, cascade change, records, history, CLI,
   restart default. Acceptance: `table` bit-identical on the three `same_results.py` cases; the unit tier passes;
   the per-step cost measured against the 5 % budget.
5. **Comparison flights**, both settings, from the d050 hand-down-off state at 178.5 s over a short window, and an early
   100 mm flight. Both settings run with `--spray-gates 2026-11-10` and `--slip-band on`, so the comparison isolates
   the dispersion change. The d050 state was saved with `spray_gates` "before" and the slip band off, and at 178.5 s
   the whole body is in the slip band (body Knudsen number between 0.01 and 0.1, every column on the Couette closure),
   so with the slip band off neither setting can spray there (decision 6); both are switched on resume. **[my call]** anything over about two hours runs on the Mac Studio through a runbook (Asha, 2026-10-31). Reported:
   floor-like share, where walls land, mass sprayed by mode, droplet sizes. Asha then decides the default.
6. **Viscous check** (analysis only, `analysis/kinked_dispersion/`). An Orr–Sommerfeld solver for the whole column:
   viscosity varying cell by cell, free surface with tension, no-slip wall at the rigid front, base speeds from the
   skins. Validated against Girin's inviscid growth at high Reynolds number, the Hooper and Boyd (1983) short-wave
   asymptote, and the sign of Charru and Hinch's (2000) thin-layer mode. Run on a few dozen sampled flight columns (d050
   and d100, several times) and on synthetic sweeps across Re_f = 1. Output: where viscous growth vanishes against the
   floor Reynolds number, and recommended `FLOOR_RE` and `KINK_DROP`. Albert and Charru (2000) is read first if it can
   be obtained.

Documentation with stage 4: `docs/model_assumptions.md` (the floor test as assumed until stage 6; the derivation ours;
the kink resting on the bridge) and the README's spray section.

## 8. Left to the builder

The eigen-solver's formulation and batching; the wavenumber grid and refinement; the fit algorithm within section 6.4's
requirements; module and function names; how the cascade carries the new run (section 6.6) — provided the books stay
closed; the fixture's column selection.

## 9. Open items and known limits

- The gas-shear check (gates stage 3) — the floor test and every Weber number depend on the skin speeds; the rough d050
  read in the header points to consistency but is one state.
- The regime limit (decision 6) means the profile path cannot spray where Girin's closure gives no delta_m; whether to
  lift it is a later decision.
- The kink's depth rests on the bridge's 906.4 K end point, which has no data behind it.
- The floor test is a scaling argument until stage 6, and assumes a Newtonian slurry.
- Albert and Charru (2000) and Aslanov and Girin (1981) are unread.
- Density differences are carried per segment; surface tension at kinks is unknown and set to zero.
- The threshold is 3.00 on the `profile` path and 3.08 on the `table` path.
- The runny layer's Reynolds number is of order 1 in the d050 thick columns, where the inviscid solver is least
  trustworthy; stage 6 is the check.

## 10. Related observation, outside this design

The front-surface Rayleigh–Taylor branch reads the liquid pool as the depth to the liquidus, which is zero in the d050
thick columns, so it cannot release there either. Changing it is outside this design (section 4); it is listed so that
Asha can decide whether it needs its own change.

## 11. Calls I made (for Asha to check)

1. The element model keeps the table; the Rayleigh–Taylor branch is unchanged; deceleration stays out of the solver;
   the exposed slurry's own mode keeps its physics but on the `profile` path reads the first flowing cell below the run.
2. The floor test uses the velocity jump relative to the moving floor and the viscosity of the cell just below the kink.
3. Defaults `KINK_DROP` 0.2, `FLOOR_RE` 1.0, `FIT_TOL` 3 % of the top speed, `N_MAX` 6, all labelled assumed until stage 6.
4. The gradient mode's gate is growth above the viscous damping 2 nu_eff k² rather than a critical Weber number;
   nu_eff, the layer's Couette viscosity, also replaces the liquid's viscosity in the gas-driven mode's viscous check.
5. The release run and the wave-fit region are defined by the wall, not by the liquidus; droplet properties are the
   released cells' mass-weighted ones.
6. `--dispersion table|profile`, default `table` until stage 5; refused without skin flow; old saved runs resume on
   `table`.
7. Records: new branch code 5; appended source-table and history columns on the `profile` path only, so `table` runs
   and the element model write exactly today's files; cost budget 5 % of a macro step.
8. The release depth is read from the wave's vertical displacement, not its stream function (about 1.4 times deeper
   on Girin's profile).
9. The stage 5 comparison runs both settings with `--spray-gates 2026-11-10` and `--slip-band on`.
