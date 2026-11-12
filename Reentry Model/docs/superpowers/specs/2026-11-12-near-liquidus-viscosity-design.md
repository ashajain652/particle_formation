# Near-liquidus viscosity: Li et al.'s law where it was measured, the literature elsewhere (design, 2026-11-12)

Status: decided with Asha on 2026-11-12, to be planned and built on the branch `near-liquidus-viscosity` (off `main`
f1ab736, which includes the skyscrapers-fixes work). Nothing here has been built. The literature behind it is the report `reports/Near liquidus aluminium alloy
viscosity.md` and the full-text notes in `research_notes/Near liquidus aluminium alloy viscosity/` (both at the repository
root, not committed), and the draft figure `reports/near_liquidus_viscosity_composite.png` with its script.

## 1. Why

The skins' flow and the exposed slurry's spraying mode take every cell's viscosity from `reentry_model/shared/rheology.py`.
Today that is Li et al.'s (2014) averaged law, eta = [0.871 - 0.00849 gamma^0.74924] exp(3.7311 f_s) Pa s, where
gamma is the shear rate (per second) and f_s the solid fraction. The shear rate is held at no more than 367 per second.
Between 906.4 K, where the model's Scheil curve stands at 10 % solid, and the 908 K liquidus, a log-linear "bridge" in
temperature joins that law to the liquid's 1.3 mPa s. The bridge was labelled "no data" when it was written. It sets
the top cells of a melting skin (907.1 to 907.6 K, 2.5 to 5.7 % solid) at 3.7 to 48 times the liquid, depending on the
shear rate. The no-spray result of the 50 mm and 100 mm skin flights rests on those cells being that viscous.

The literature check of 2026-11-12 found three things. First, no measurement exists for an aluminium alloy between the
liquid and 10 % solid. Second, Li et al.'s averaged law misrepresents their own data above 367 per second: it reaches
zero at 483 per second, although they measured 0.043 Pa s at 490 per second. Third, everything Li et al. measured lies
inside 10 to 50 % solid and 61 to 490 per second, while the gas drives the top cells at thousands to hundreds of
thousands per second.

Asha's decision (2026-11-12): use Li et al.'s law only where it is valid, and construct the law elsewhere from the other
literature. She accepted all five recommendations put to her the same day:

1. Inside Li et al.'s measured window, use their per-rate fits (their Table 3), not the averaged law.
2. Below 10 % solid, the structural part of the viscosity falls as the square of the solid fraction by default (Pan et
   al.'s form, the stiffer choice), with the cube (the Loue-Sigworth form) as an option.
3. The fluid under the high-shear limit is the real liquid (1.3 mPa s), not the 13 to 18 mPa s that rheometer fits
   imply.
4. The present law stays selectable by name for comparison runs; the new law is the default, and a run name marks the
   setting when it is not the default.
5. The solver keeps history: each cell's previous shear rate is carried, so a structured cell breaks down only when the
   stress on it passes the peak of its flow curve, and then stays broken down.

Asha also asked (2026-11-12) that the sources justifying each choice be documented. Section 4 does this. The same
attributions go into the module's docstring and into `docs/model_assumptions.md` (section 9 of this document).

## 2. Goals, non-goals and anti-goals

Goals:

- A new default viscosity law, `composite`. It is defined for every solid fraction from 0 to just under 0.5 and every
  shear rate from 0 upward. It reproduces Li et al.'s per-rate fits at every point inside their measured window, equals
  the liquid's viscosity at zero solid, and is continuous in both solid fraction and shear rate.
- The cells' shear-rate history, as in section 5: carried by every cell operation, saved with the flight, and used to
  choose which branch of the flow curve a cell sits on.
- The present law kept as `bridge`. A flight run with `--rheology bridge` gives the same numbers as today to the last
  bit.
- The sources of every choice written where a reader of the code or of the thesis will find them.
- A cheap first comparison of the two laws: a probe of the skins alone on the laptop (section 8, level 4). Resumed
  flights come afterwards, only if Asha decides (level 5).

Non-goals:

- No thixotropic kinetics: no rebuilding of structure while a cell rests, and no transient stress overshoot after a
  sudden rise in shear rate. Section 10 records both as limitations.
- No change to the liquid's viscosity value (1.3 mPa s from the material file). The following are also unchanged: the
  Scheil curve, the rigid threshold of half solid, the slurry's surface tension (`SLURRY_SIGMA`), and the element model
  (`--surface-model elements` does not use this module).
- No wiring into the level-set architecture. The law sits in the shared package, ready for it.
- No production flights. The Mac Studio flights that follow, if any, are Asha's decision and run from a runbook.
- No edits to other branches. The `kinked-dispersion` branch's analysis scripts read `rheology.BRIDGE_LOW`, so that
  constant and the present law's other names stay.

Anti-goals, which would make the build a failure even if it ran:

- The results of `--rheology bridge` changing at all.
- The law below 10 % solid reading as measured. Each region's source and its status (measured, fitted or assumed) must
  be visible in the docstring, in `model_assumptions.md` and in the run's JSON.
- A numerical floor or clamp acting as physics without being stated. Here that means the zero-rate floor of section
  3.7.

## 3. The law (decided)

### 3.1 Variables

The law takes four inputs:

- the cell's solid fraction f_s, from its temperature through the material's Scheil curve, as now;
- the cell's shear rate gamma, in per second;
- the liquid's viscosity mu_l, which is the material's `liquid.mu`, 1.3 mPa s for AA7075;
- the structure exponent p (section 3.5).

Temperature enters only through f_s. Li et al. report their data by solid fraction from their own DSC scan, which puts
10 % solid at 630.3 C (903.45 K). The model's Scheil curve puts it at 906.4 K. Matching by solid fraction rather than by
temperature is a builder's call recorded in section 11.

The law has two parts, eta = eta_inf(f_s) + X(f_s, gamma). The first, eta_inf, is the viscosity of loose solid grains
in the liquid, the high-shear limit. The second, X, is the structural excess above it.

### 3.2 The suspension part (the high-shear limit)

    eta_inf(f_s) = mu_l [1 + S / (1/f_s - 1/phi_m)],   S = 3,  phi_m = 0.55,   and eta_inf(0) = mu_l.

This is the Mori-Ototake law. S = 3 is its value for spheres, and phi_m is the maximum packing fraction. It gives 1.06,
1.17, 1.37 and 1.94 times the liquid at 2, 5, 10 and 20 % solid, and 17.5 times at 50 % solid.

### 3.3 Inside Li et al.'s measured window (0.1 at most f_s, below 0.5; 61.24 to 489.88 per second)

Li et al. fitted eta = A_i exp(B_i f_s) Pa s at each measured rate gamma_i. Their Table 3 gives:

| gamma_i (per second) | 61.24 | 122.47 | 244.94 | 367.41 | 489.88 |
|---|---|---|---|---|---|
| A_i (Pa s) | 0.66711 | 0.59623 | 0.33326 | 0.14323 | 0.04346 |
| B_i | 4.78 | 3.52 | 3.45 | 3.66 | 3.23 |

Between the measured rates, ln eta is interpolated linearly in ln gamma:

    L(f_s, gamma) = exp( interp(ln gamma; ln gamma_i, ln A_i + B_i f_s) ).

Inside the window the law is exactly L, so X = L - eta_inf. Li et al.'s fit is reproduced at all 25 of their
(solid fraction, rate) points, to round-off.

### 3.4 Outside Li et al.'s rates, at 10 % solid or more

The excess at the nearer edge of the window continues as a power of the shear rate:

    X(f_s, gamma) = [L(f_s, gamma_c) - eta_inf(f_s)] (gamma / gamma_c)^n_x,
    gamma_c = min(max(gamma, 61.24), 489.88),   n_x = -0.85.

With n_x = -0.85, the excess stress X gamma grows only as gamma^0.15, nearly a constant stress. At high rates the
viscosity therefore falls toward eta_inf, and at low rates it rises as gamma^-0.85 without a plateau.

### 3.5 Below 10 % solid

The structural excess is Li et al.'s at 10 % solid, scaled to vanish at the liquid:

    X(f_s, gamma) = X(0.1, gamma) (f_s / 0.1)^p,   p = 2 by default, 3 as the option.

The law is continuous at 10 % solid by construction. At zero solid it is the liquid.

### 3.6 Rigid

At half solid or more the viscosity is infinite, as now (`RIGID_SOLID_FRACTION` = 0.5). The jump at half solid is
unchanged from the present law.

### 3.7 The zero-rate floor

The law is evaluated at gamma* = max(gamma, GAMMA_FLOOR), with GAMMA_FLOOR = 1e-3 per second. Below the floor the
viscosity is constant (Newtonian creep). The floor exists for two reasons:

- Without it, X diverges at zero rate. Every test the skin code makes with `np.isfinite(cell_viscosity(sk, 0.0))`
  (eight call sites in `skin.py`, one in `skinbody.py`) would then call every slurry cell rigid.
- Two flowing cells with infinite viscosity would lose their drag on each other (`series_drag`).

At the floor the viscosity is about 1.3e4 Pa s at 10 % solid and 6.7e4 Pa s at 45 % solid. That is far below the
1e9 Pa s crust coupling whose round-off sub-plan 21 Task 6b measured. The floor rate moves a 1 mm layer by about
1 micrometre per second, or 0.2 mm over a three-minute flight. Its value is a builder's call (section 11).

### 3.8 Reference values for the tests

The table gives eta / mu_l (p = 2) at the solid fractions of the top cells and of Li et al.'s window, across rates from
the floor to 1e5 per second. It was computed from sections 3.2 to 3.7 by an independent script (scratch, 2026-11-12). A
test should match it to a relative 1e-3, the precision of four printed figures.

| f_s | 0.001 per s (floor) | 1 | 61.24 | 244.94 | 489.88 | 1e4 | 1e5 |
|---|---|---|---|---|---|---|---|
| 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 0.025 | 6.053e5 | 1707 | 52.72 | 23.62 | 3.879 | 1.294 | 1.109 |
| 0.057 | 3.147e6 | 8870 | 269.7 | 118.4 | 15.75 | 2.312 | 1.349 |
| 0.0975 | 9.207e6 | 2.595e4 | 786.8 | 344.2 | 43.95 | 4.636 | 1.819 |
| 0.1 | 9.685e6 | 2.730e4 | 827.7 | 362.0 | 46.18 | 4.818 | 1.854 |
| 0.2 | 1.562e7 | 4.404e4 | 1335 | 511.1 | 63.78 | 6.705 | 2.616 |
| 0.45 | 5.159e7 | 1.454e5 | 4410 | 1211 | 143.0 | 18.79 | 9.889 |
| 0.025, p = 3 | 1.513e5 | 427.6 | 13.99 | 6.713 | 1.779 | 1.132 | 1.086 |

Two conclusions follow from the table. Above about 1e4 per second, metal at up to 10 % solid is within a factor of 5 of
the liquid. Below Li et al.'s lowest rate, the same metal is hundreds to millions of times the liquid. Which of these a
cell sees depends on the stress on it and on its history (section 5).

### 3.9 The shape of the flow curve

The stress the law gives, tau(gamma) = eta gamma, is not monotone in the shear rate:

- It rises from zero to a peak at Li et al.'s 244.94 per second.
- It falls to a valley at 489.88 per second.
- It rises again without limit.

The peak and valley stresses are:

| f_s | 0.025 | 0.038 | 0.057 | 0.0975 | 0.1 | 0.2 | 0.3 | 0.4 | 0.45 |
|---|---|---|---|---|---|---|---|---|---|
| Peak stress (Pa), p = 2 | 7.5 | 16.9 | 37.7 | 109.6 | 115.3 | 162.7 | 229.8 | 324.5 | 385.5 |
| Valley stress (Pa), p = 2 | 2.5 | 4.8 | 10.0 | 28.0 | 29.4 | 40.6 | 56.1 | 77.5 | 91.1 |
| Peak / valley (Pa), p = 3 | 2.1 / 1.1 | 6.7 / 2.3 | 21.6 / 6.0 | 106.9 / 27.3 | | | | | |

Close to half solid (0.49), the curve has a second, smaller peak at 61.24 per second (425 Pa) and a dip at 122.47 per
second (410 Pa). This happens because Li et al.'s B_i is largest at their lowest rate. The branch rule of section 5 must
therefore handle any number of rising and falling pieces, not just one peak.

Between a valley and its peak, a given stress has more than one root (a slow, structured state and a fast, broken-down
one), and the middle root is unstable. The present solver already meets the falling piece of Li et al.'s averaged law
(229 to 367 per second): its fixed point drifts off the unstable root and stops at its five-pass cap. The fact for
sub-plan 21 Task 6d in `00-shared-context.md` records this. The new law's falling piece is wider and steeper, which is
why section 5 changes the solver.

## 4. Sources for each choice

The PDFs are in `Research/Literature Review/Near Liquidus Viscosity/` (above the repository). The page-by-page
extraction, with every value read from a figure and its uncertainty, is in the research notes at the repository root:
`fulltext_li_2014.md`, `fulltext_moon.md`, `fulltext_joly_koke_pan.md` and `fulltext_fan_wang_assael_atkinson.md`.
"My calculation" in those notes marks a number derived from a paper rather than stated in it. None of these papers
measured an Al-Zn-Mg-Cu (7xxx) alloy except Li et al. (rheology) and Atkinson et al. (microstructure only).

| Choice | What it takes | Source, with location | Status and caveat |
|---|---|---|---|
| Use Li et al.'s data inside 10-50 % solid and 61.24-489.88 per second | The measured window: f_s 0.1, 0.2, 0.3, 0.4, 0.5 at 630.3 to 620.7 C; five rates | Li et al. 2014, Table 2 (pp. 80-81) and the method section (pp. 79-81); 10 % solid is a measured point (Figs. 3, 5, 8, 9) | Measured on AA7075. The slurry was solidified while sheared (50 C superheat, cooled at 5 C/min), not remelted. No calibration with a known liquid; the outer cup rotates, so there are no Taylor vortices |
| Per-rate fits, not the averaged law | A_i, B_i at each rate | Li et al. 2014, Table 3 (pp. 82-83); the averaged law is their Eq. (5) on the same pages | Fitted to the data. Eq. (5) reaches zero at 483 per second although A(489.88) = 0.04346 Pa s was measured; at 61 per second it reads 1.4 to 1.9 times low at 30-50 % solid |
| Log-log interpolation between Li et al.'s rates | Straight lines in ln eta against ln gamma through the five fits | None: Li et al. give no form between rates (their Fig. 7 fits a power law to A(gamma) for the averaged law) | Assumed. It is the least assumption that passes through every fit |
| Suspension part eta_inf: Mori-Ototake, S = 3, phi_m = 0.55 | The high-shear limit and the dilute limit | Moon 1990 (MIT PhD thesis), pp. 101-102, Eq. 51 and Fig. 5.14 (p. 191); the law's form in Table 2.1 (p. 209); phi_m = 0.55 from the infiltrated compact (p. 102); Moon, Cornie and Flemings 1991, MSE A144, 253-265 | Measured on rigid SiC particles in liquid Al-6.5Si, not on solidifying metal. High-shear limits after pre-shear at 900 per second lay on the sphere line (S = 3); after 180 per second, on the angular-particle line (S = 11.8). Moon's ratios are taken against his rig's 0.020 Pa s "liquid", 15-16 times the true value (Table 5.1, p. 216), so they hold only if that error is proportional |
| S = 3 (spheres) rather than 11.8 (angular) | Rounded grains | Atkinson, Burke and Vaneetveld 2008, MSE A490, 266-276: remelted extruded 7075 is fully spheroidal by about 5 % liquid (abstract; pp. 269-274), grains 50-100 micrometres (Fig. 6, p. 272; Fig. 11, p. 275) | Inferred. Atkinson et al. held samples for 10 min; a re-entering sphere heats far faster |
| The liquid under the high-shear limit is the real liquid | mu_l = 1.3 mPa s inside eta_inf, not 13-18 mPa s | Koke and Modigell 2003, JNNFM 112, 141-160, Eq. 16 (p. 158): eta tends to eta_l at zero solid and at infinite rate. Fan 2002, Int. Mater. Rev. 47, 49-85, pp. 52 and 54: the asymptotic viscosity depends only on solid fraction and is a suspension law in the liquid's viscosity (Chen and Fan's and Brown's models) | Asha's decision 3. The rheometer-implied 13.4 mPa s matrix (Fan's Fig. 8, p. 54, by my calculation) is about 5 times the Sn-Pb liquid, and Moon's rig read the liquid at 15-16 times true; both are taken as instrument floors |
| The liquid's value, 1.3 mPa s | Unchanged material value | Assael et al. 2006, J. Phys. Chem. Ref. Data 35, 285-300, Eq. 19 (p. 295): 1.42 mPa s at 908 K, extrapolated 25 K below its range, uncertainty 14 % | Pure aluminium. 1.3 mPa s is 8 % below, inside the uncertainty. Wang and Overfelt 2002 (Int. J. Thermophys. 23, 1063-1076) show alloying moves it by up to 25 %; no 7xxx data exist |
| The excess falls as gamma^-0.85 outside Li et al.'s rates | n_x = -0.85: the excess stress grows as gamma^0.15 | Joly and Mehrabian 1976, J. Mater. Sci. 11, 1393-1418, Table IV (p. 1415): n = -0.82 at 10-40 per second, 45 % solid (Sn-15Pb). Fan 2002, Fig. 4 (p. 52), Turng and Wang's Sn-15Pb: -0.84 to -0.94 between 105 and 300 per second (my calculation). Moon 1990, Eq. 47 (pp. 83, 89): n' = -0.85 for Al-6.5Si at 40 % solid; Table 5.2 (p. 217): -0.68 to -0.85 for SiC suspensions | Measured on other alloys at 10-1400 per second; beyond that it is extrapolation. Above 300 per second, Turng and Wang's excess falls faster than this (about 12 Pa at 105 per second to 2 Pa at 1900 per second over their asymptote, at 17 % solid), so -0.85 errs on the viscous side at 1e3-1e5 per second |
| No zero-shear plateau below 61 per second | X keeps rising as the rate falls | Joly and Mehrabian 1976, Table IV and section 1.8 of the notes (near-constant 58-74 Pa at 10-40 per second); Fan 2002, Fig. 1 (p. 51): flat stress at about 35 % solid below 50 per second; Moon 1990, p. 99 and paper p. 258: his low-shear Newtonian region was assumed, never measured | Measured as a yield-like stress on other alloys; no data below 10 per second |
| Below 10 % solid, scale the excess by (f_s/0.1)^p | p = 2 default, 3 option | Pan, Apelian and Alexandrou 2004, Metall. Mater. Trans. B 35, 1187-1202: Eq. 9 (p. 1192) makes the yield stress vanish as f_s squared; the Loue-Sigworth law, Eq. 11 (pp. 1197-1198), as f_s cubed; both give about 19-24 Pa at 10 % solid for Al-Si, near the bottom of Li et al.'s 25-100 Pa. Koke and Modigell 2003, p. 158: at zero solid the alloy is the Newtonian liquid | Assumed. No source measured below 31 % solid (Pan) or 10 % solid (any). The exponent is a property of the fitted forms, not of data. p = 2 is Asha's decision 2 |
| Rigid at half solid | Infinite viscosity | Unchanged from the layered-skin design (spec 6.4); Li et al. saw torque fluctuations from particle bridging above 50 % solid (pp. 82-83) | As before |
| The flow curve's fall, and history decides the branch | A stress peak near 245 per second, a valley at 490; a cell keeps its branch until the stress forces a change | Li et al. 2014, Table 3 (the falling stress) and Fig. 11 (finer, rounder grains at higher rates: each experiment had its own structure). Koke and Modigell 2003, Fig. 15 (p. 155): static yield stress about 33 Pa going up, dynamic about 16 Pa coming down (hysteresis). Joly and Mehrabian 1976, pp. 1405-1407: each isothermal point is its own structure, and changes after a rate step take minutes. Moon 1990, pp. 75-76 and 84, Fig. 4.7: structure curves | The fall is in Li et al.'s data. Choosing the branch by continuation (section 5.2) is this design's way of carrying out Asha's decision 5 |
| Zero-rate floor of 1e-3 per second | Newtonian creep below it | Koke and Modigell 2003, pp. 155-156: "a minimal creeping behaviour" below the yield stress. Pan et al. 2004, Eq. 12 (p. 1198): a regularisation (Papanastasiou) used as a numerical convenience | Numerical. The value is a builder's call |
| Solid fraction, not temperature, carries Li et al.'s data | Their 10 % solid lands at 906.4 K | Li et al. 2014, Table 2 (pp. 80-81): 10 % at 903.45 K by their DSC. Fan 2002, Table 1 (p. 68): 7075's Scheil liquidus is 634.9 C, matching the model's 908 K | Builder's call. Mapping by temperature would move Li et al.'s window 3 K colder on the model's curve |

## 5. History: the cells' shear rate and the solver

### 5.1 The state

Each skin cell carries its last shear rate, `Skins.shear_rate`, an (n, K) array in per second, the same shape as `T`. It
is part of the cells' state, like `U` and `S`:

- Split halves keep it.
- Merged cells take the mass-weighted mean.
- The sideways transport carries it with the metal, mass-weighted, as it carries `m S`.
- Metal drawn from the 3D model arrives at rest (zero).
- A rigid cell, including a floating crust's cells, holds zero. This means metal that freezes and remelts starts
  structured.
- A liquid cell holds its own rate (its stress over mu_l). Metal that cools from the liquid while sheared therefore
  starts on the broken-down branch. That agrees with Li et al.'s finer, rounder grains under faster shear.

The array is part of `state()`, `set_state()`, `rows()`, the flow's write-back of a block of columns, every top and
bottom operation, and the restart file. A state saved before this change, such as the d050 state at 178.5 s, has no
such array. On load it gets zeros: every cell starts structured, and the first momentum sub-step breaks down those
whose stress passes their peak.

### 5.2 The continuation root

Given a cell's solid fraction, a stress tau and the rate it starts from, gamma_0, the cell's new rate is found by
moving along its own flow curve in the direction the stress demands, to the first root:

- If tau(gamma_0) is below tau, the rate rises from gamma_0 to the first rate where the curve reaches tau.
- If tau(gamma_0) is above tau, the rate falls from gamma_0 to the first rate where the curve comes down to tau.

Because tau(gamma) goes to zero as the rate goes to zero and grows without limit at high rate, the root always exists.

What follows:

- A structured cell whose stress passes its peak jumps to the broken-down branch.
- A broken-down cell stays there until its stress falls below the valley.
- The unstable middle root is never chosen.

A test must cover the curve with two peaks near half solid (section 3.9). The vectorised method (bracketing on the
law's pieces, then Brent or Newton in ln gamma) is the builder's choice.

### 5.3 The momentum step (`skin._momentum_columns`) under `composite`

1. The cells' viscosities start from their stored rates (`shear_rate`), not from rest.
2. Each pass solves the velocities with those viscosities. It then takes every flowing cell's stress, the mean over
   its two faces as `rates()` computes it now. It sets the cell's rate to the continuation root of that stress from the
   cell's rate in the previous pass (the stored rate, in the first pass), and its viscosity to the law at that rate.
3. Passes end when no flowing cell's viscosity moves by more than `SHEAR_TOL`.
4. The converged rates are stored in `shear_rate` at the end of the sub-step.

Continuing from the previous pass's rate makes the passes keep the same hysteresis as the sub-steps. A cell that
breaks down in one pass stays broken down for the rest of the sub-step unless its stress falls below the valley. A
structured cell stays structured unless its stress passes the peak.

Continuing from the start-of-sub-step rate in every pass would let a cell near its threshold flip between passes. In a
flow sub-step (0.25 to 10 ms) shorter than a broken-down layer's viscous time (about 0.02 s for 0.1 mm), the stress a
cell carries depends on its own viscosity. So a cell loaded just above its peak would:

1. carry a stress above its peak with its structured viscosity, and break down;
2. carry a stress below its peak once broken down, and return to the structured branch.

The pass cap's parity, not the physics, would then decide whether it breaks down. The builder check of 2026-11-12
found this.

In a column driven at its top by the gas and through its depth by the push, the stresses are set mostly by the forces,
so this iteration on stress should converge in a few passes. Where it does not, the cap ends it, as now. The book must
record the passes and the number of flowing cells that had not converged, so a test and the flights can measure them.
The pass cap may be raised for `composite` if the measurement shows it is needed. `SHEAR_TOL` stays at 1e-3.

### 5.4 The exposed slurry's mode (`skinbody.slurry_modes`) under `composite`

The fixed point there starts from the exposed cell's stored rate instead of from rest. It uses the same continuation
root: the layer's stress from `lubrication()` (mu times its shear rate) gives the rate, and the law at that rate gives mu.
The rate is not held at `GAMMA_MAX`.

The cascade's gate `T_slurry` is unchanged. It assumes a colder slurry cell is more viscous, which history can break
(a colder cell may sit on the broken-down branch). The gate then only delays such a cell's release to the next macro
step's evaluation. It never releases metal it should not.

### 5.5 Under `bridge`

Everything is today's code path, unchanged:

- the present `rheology.viscosity` with the bridge and `GAMMA_MAX`;
- the fixed point from rest with its five passes;
- `slurry_modes` as now.

The new array may still be carried, but nothing under `bridge` reads it, and the numbers must be identical to today's
to the last bit.

## 6. Interfaces and options

`reentry_model/shared/rheology.py`:

- The present names stay with their present meaning and numbers: `li_viscosity`, `viscosity`, `LI_A`, `LI_B`, `LI_N`,
  `LI_C`, `GAMMA_MAX`, `BRIDGE_LOW`, `BRIDGE_HIGH` and `RIGID_SOLID_FRACTION`. The `kinked-dispersion` branch and the
  tests use them.
- New constants: `LI_RATES`, `LI_RATE_A` and `LI_RATE_B` (section 3.3), `LI_SOLID_MIN` = 0.1, `SUSPENSION_S` = 3.0,
  `SUSPENSION_PHI_M` = 0.55, `EXCESS_RATE_EXPONENT` = -0.85, `STRUCTURE_EXPONENTS` = (2.0, 3.0) with
  `STRUCTURE_EXPONENT` = 2.0 the default, `GAMMA_FLOOR` = 1e-3, `LAWS` = ("composite", "bridge") and `DEFAULT_LAW` =
  "composite".
- New functions:
  - `composite_viscosity(f_s, gamma, mu_liquid, p=STRUCTURE_EXPONENT)` returns the viscosity in Pa s, vectorised,
    infinite at half solid or more.
  - `flow_stress(f_s, gamma, mu_liquid, p)` is the viscosity times the rate.
  - `continued_rate(f_s, tau, gamma_from, mu_liquid, p)` is section 5.2's root, vectorised.
- The module docstring states each region with its source and status, as in section 4.
- The module still imports nothing from the rest of `reentry_model` (the shared-package rule that
  `tests/test_reentry_model_shared.py` checks).

`reentry_model/skyscrapers/skin.py`:

- `SkinSettings` gains `rheology` = "composite" and `structure_exponent` = 2.0, validated against `LAWS` and
  `STRUCTURE_EXPONENTS`.
- `cell_viscosity(sk, gamma)` dispatches on `sk.settings.rheology`.
- `Skins` gains `shear_rate` (section 5.1).
- A `Skins` or `SkinSettings` loaded from a state saved before this change reads as `bridge` for the law and zeros for
  the rates. How it does this is the builder's choice; one option is a default through `getattr` or a `__setstate__`.

Command line (`cli.py`):

- `run --rheology composite|bridge` (default composite) and `--structure-exponent 2|3` (default 2) apply to melting
  skin runs.
  - Giving `--structure-exponent` with `--rheology bridge` exits 2 with a message.
  - Giving either with `--surface-model elements` or without `--melt on` is accepted and ignored, as the other skin
    options are.
- `resume` takes both. Its default is the saved flight's setting, and `bridge` for a state saved before the option
  existed. A change is recorded in `changed`, as the other resume options are.
- Run names: `skin_name_suffix` adds `_rheology-bridge` for the present law and `_structexp3` for p = 3. The defaults
  are unmarked, as Task 5g left its changed default. A default melting run made after this change therefore has the
  same name as one made before it, under a different law. The run's JSON records `rheology` and `structure_exponent`
  under `skin` always, as it does `part_count_from`, so the two can be told apart.
  - `README.md`, `CLAUDE.md` and the run command's usage text list the two options.

Under `composite`, any code that needs a cell's viscosity at its working point reads the cell's stored rate. It never
iterates from rest, because iterating from rest would put every cell in the two-valued stress band on its structured
branch. For example, at 9.75 % solid under 61 Pa that gives 913 times the liquid instead of 4.55 times (section 7),
against Asha's decision 5. This matters for the `kinked-dispersion` branch (pushed at 45ea185, not merged), which edits
the same functions:

- It adds `skin.cell_shear_rates`, which iterates every cell's viscosity from rest to feed the profile-based
  dispersion.
- It changes `girin_plan`, the call to `slurry_modes`, `cascade_mass` and `release_top`.
- It records the depth of the 906.4 K isotherm.

Whichever branch merges second adapts:

- If `kinked-dispersion` lands first, this build makes `cell_shear_rates` read the stored rates under `composite`,
  leaving its numbers under `bridge` unchanged.
- If this build lands first, the dispersion session makes the same change when it merges.
- The record of the 906.4 K depth stays as that branch wrote it. What replaces it as the kink's depth under the new law
  is the dispersion session's decision.

An analysis script, `Reentry Model/analysis/near_liquidus_viscosity.py`, replaces the uncommitted draft at the
repository root. It plots the law against Li et al.'s points and the present law, prints section 3.8's table and
section 3.9's peaks and valleys, and calls the package's functions rather than copying them.

## 7. What the change should do (single-cell numbers, to be measured in the flights)

At a given stress, the law has a structured root (reached from rest) and a broken-down root (reached from a high
rate). These are equal where the stress is above the peak or below the valley. The table gives eta / mu_l at each,
structured / broken down, for p = 2.

| Stress on the cell | 2.5 % solid | 3.8 % | 5.7 % | 9.75 % | 20 % |
|---|---|---|---|---|---|
| 10 Pa | 1.45 / 1.45 | 118 / 2.46 | 19 900 / 19 900 (0.4 per second) | 9.2e6 (at the floor) | 1.6e7 (at the floor) |
| 61 Pa | 1.14 / 1.14 | 1.29 / 1.29 | 1.65 / 1.65 | 913 / 4.55 | 30 700 / 14.0 |
| 101 Pa (the median gas shear) | 1.12 / 1.12 | 1.23 / 1.23 | 1.46 / 1.46 | 545 / 2.69 | 1 780 / 5.49 |
| 170 Pa | 1.10 / 1.10 | 1.19 / 1.19 | 1.35 / 1.35 | 2.03 / 2.03 | 3.49 / 3.49 |

The present law gives the 2.5 % and 5.7 % solid cells 3.7 to 5.6 and 18.5 to 48 times the liquid. With p = 3, the
5.7 % cell at 10 Pa is 476 times structured and 3.5 times broken down.

Read with the history rule, the table says:

- The top cells (2.5 to 5.7 % solid) break down under any gas shear above about 38 Pa and then run within a factor of 1.7
  of the liquid.
- Metal near 10 % solid stays structured, at hundreds of times the liquid, unless the stress on it has passed about
  110 Pa, or unless it cooled from the liquid while sheared.

A caveat: the stress a cell carries is the gas shear plus the push (pressure gradient, deceleration, thickness and
capillary terms) on the metal above it, not the gas shear alone. The flights measure what the cells actually carry.
Faster top layers also raise the adaptive count of flow sub-steps (their Courant number), so a skin flight may cost
more per macro step.

## 8. Tests and acceptance

Level 1, the law (`tests/test_reentry_model_rheology.py`):

- Li et al.'s 25 fitted values reproduced to 1e-12.
- Section 3.8's table to 1e-3.
- eta equals mu_l at zero solid, and is infinite at 0.5 and above.
- Continuity at 10 % solid and at both rate edges of Li et al.'s window.
- eta tends to eta_inf as the rate grows, and the floor applies below 1e-3 per second.
- Section 3.9's peaks and valleys.
- The present law's existing tests unchanged and passing.

Level 2, the continuation root:

- A stress ramped up from rest past the peak jumps to the broken-down branch.
- Ramped back down, the rate stays on that branch until the stress passes the valley, then falls.
- The middle root is never returned.
- The two-peak curve near half solid is handled.
- The function is vectorised over cells with mixed starting branches.

Level 3, the skins:

- `shear_rate` carried exactly by every cell operation:
  - split, merge (mass-weighted), the top and bottom operations and the sideways transport (mass-weighted, books
    exact);
  - `state()`/`set_state()`, `rows()` and the flow's block write-back;
  - save and resume, and an old state loading with zeros.
- On a synthetic column under a fixed gas shear and push, under `composite`:
  - the momentum books close to round-off as now;
  - the iteration's passes and unconverged count are recorded;
  - a structured cell under a stress just below its peak stays structured, and just above, breaks down.
- Under `bridge`: the existing skin and skinbody tests pass unchanged, and the port's comparison cases give identical
  outputs (`analysis/skyscrapers_port/same_results.py`, as the skyscrapers-fixes branch checked).
- The command-line options, the run-name suffixes, the JSON record and resume's defaults.
- The unit tier passes (`pytest tests -m "not drama and not reference" -n 4`).

Level 4, a probe of the skins alone on the laptop (Asha, 2026-11-12, option A; no flight):

- Why not a resumed flight. The state below was flown with slurry spraying off, and the saved flight's rows show no
  liquid on any column at 178.5 and 179.0 s (the deepest liquid layer is 0 mm; 0.006 mm at 179.5 s). The top cells, at
  2.5 to 5.7 % solid, count as slurry, so with spraying off neither law would release anything. Each 0.5 s macro step
  from this state also took 63 to 76 minutes of skin time on the laptop under the present law, 93 % of it in the flow
  sub-steps. That is about 2.5 hours for two steps, and more under the new law (not measured).
- The state: `d050__handdown-off.state.pkl` at 178.5 s, in the subplan20 worktree under
  `prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/`, loaded with the package.
- The probe:
  - Run only the skins' flow sub-steps (momentum, sideways transport, thickness), with the gas loads (wall shear,
    pressure gradient, deceleration) held at the values the skin body computes from the saved state for its next macro
    step. Run no heat sub-step, no 3D conduction and no trajectory step.
  - Run about 0.1 s of flow, with the flow sub-steps chosen adaptively as in a flight.
  - Run four cases:
    - `bridge`;
    - `composite` with p = 2, the cells starting structured (zero rates, as a resumed flight would start);
    - `composite` with p = 2, every cell starting broken down (a very high rate), which brackets the history from the
      other side;
    - `composite` with p = 3, starting structured.
  - At the start, at the end and at a few times between, evaluate the liquid's and the exposed slurry's spraying modes
    from the skins as they stand (`slurry_modes` and the liquid's mode, as with `--slurry-spray on`), without releasing
    any metal.
- The script is `Reentry Model/analysis/near_liquidus_probe.py` (committed). Its report is
  `output/near_liquidus_viscosity/probe.md` (git-ignored). For each case it gives:
  - the top cells' viscosity ratios;
  - the stresses the cells carry (gas shear plus push);
  - the share of slurry cells on each branch;
  - the surface speeds;
  - the liquid's and the slurry's Weber numbers against the critical value, and the release rates the modes would
    give;
  - the passes and unconverged counts;
  - the flow sub-step counts and the wall time per flow sub-step.
- What it cannot show: feedback. No metal leaves, and the heating does not respond, so it describes the state at
  178.5 s, not the rest of the flight.
- A single-column script runs one representative windward column's profile to a steady state under both laws. It also
  runs n_x at -0.7 and -1.0 as a sensitivity check.

Level 5, resumed flights (after level 4, only if Asha decides it): resume the same state at 178.5 s under both laws
with `--slurry-spray on`, to 179.5 s, on the Mac Studio by a runbook. This captures the feedback that level 4 leaves
out. The estimated cost is about 2.5 hours under `bridge` and 6 to 9 hours under `composite`, not measured.

## 9. Documentation

- `docs/model_assumptions.md`: a new section on the near-liquidus viscosity (shared physics). It gives the law region by
  region, with section 4's table of sources and statuses, and the limitations of section 10.
- The layered-skin design (`2026-10-11-layered-skin-moving-mesh-design.md`): a dated amendment where it states the
  present law: the table at line 67, section 6.4 near line 343, the slurry mode near lines 611-617, line 958, and
  assumption 5 near line 1082. It says the default is now this design's law and the bridge is kept as `--rheology
  bridge`.
- `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`: a new fact recording the change, the measured
  convergence and the level-4 results. The next free number was 196 on 2026-11-12.
- `README.md` and `CLAUDE.md`: the options and the changed default.
- The `kinked-dispersion` session must be told. It places its kink at the depth of the bridge's 906.4 K, which the new
  law no longer has; the natural replacement is the depth where the stress stops breaking the structure.

## 10. Limitations that stay

- No data between the liquid and 10 % solid. The exponent p, and where the structured slurry gives way to loose grains,
  are inferences.
- Li et al.'s slurry was solidified while sheared. A re-entering sphere remelts wrought metal, whose grains start
  unsheared and probably agglomerated (Fan 2002, pp. 59-67: breakdown takes seconds, rebuilding hundreds to thousands of
  seconds). Atkinson et al. show the same class of rounded grains in liquid, but not whether the structure at 10 % solid
  is as strong.
- Li et al.'s rig was not calibrated with a known liquid, so the lowest values (0.05-0.2 Pa s at 490 per second) may
  include an instrument floor.
- The law is a steady-state law. Koke and Modigell (2003, Fig. 18, p. 157) found the stress just after a jump in rate
  rising as the rate to the power 2.07, so the first moments after a sudden load are stiffer than any law here.
- The history never rebuilds structure. A broken-down cell stays broken down at any stress above its valley. Joly and
  Mehrabian (1976, Figs. 19-20, p. 1408) found rest-built strength of tens to hundreds of pascals within seconds to
  minutes.
- Above 300 per second, the excess may fall faster than gamma^-0.85 (section 4).

## 11. Left to the builder, and the calls made in this design

Left to the builder:

- the root-finding method and any relaxation in the momentum iteration;
- the pass cap under `composite` (measured, section 5.3);
- how old states get their defaults;
- the internal structure of the tests;
- the single-column script;
- how the level-4 probe obtains the saved state's gas loads and drives the skins' flow alone, and its exact length
  (about 0.1 s).

Calls made in this design that Asha did not decide (she should check these):

- Solid fraction, not temperature, carries Li et al.'s data (section 3.1).
- GAMMA_FLOOR = 1e-3 per second (section 3.7).
- n_x = -0.85 is a constant, not a command-line option. Sensitivity at -0.7 and -1.0 can be run in the analysis script
  and the single-column test, not in flights.
- `--structure-exponent` takes only 2 or 3, the two sourced forms.
- Merged cells take the mass-weighted mean rate. Rigid cells hold zero. Old states start structured (section 5.1).
- Under `composite`, the cascade's `T_slurry` gate is unchanged (section 5.4).
- The default's run name is unmarked; `_rheology-bridge` and `_structexp3` mark the others (section 6).
- `resume` defaults to `bridge` for states saved before the option existed (section 6).
- Inside a sub-step, the passes continue from the previous pass's rate, so a breakdown found in one pass stands
  (section 5.3).
- Under `composite`, code needing a cell's working viscosity reads its stored rate; with `kinked-dispersion`,
  whichever branch merges second adapts (section 6).
