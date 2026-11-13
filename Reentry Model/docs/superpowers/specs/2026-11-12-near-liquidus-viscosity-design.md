# Near-liquidus viscosity: a suspension law below 10 % solid, Li et al.'s most broken-down 7075 above it (design, 2026-11-12, revised)

Status: decided with Asha on 2026-11-12. This document was revised the same day after two literature reviews (section
1.2) and replaces the version committed at f143e19, whose plan (3d0037d, d8ca92e) is withdrawn. It is to be planned
again and built on the branch `near-liquidus-viscosity` (off `main` f1ab736). Nothing here had been built when it was
written.

Built 2026-11-12 on that branch: commits 0f7cfff, cd56947, 136e9d8, 3e09ec6, 32c37f3, b67866e, 2d45858 and the final
review fix commit; the shared-context fact is fact 196 (docs/superpowers/plans/melt-spraying-subplans/
00-shared-context.md); the probe report is output/near_liquidus_viscosity/probe.md (git-ignored).

Merged into main 2026-11-13 (Asha, 2026-11-12: "integrate the branch, and merge branch onto origin/main. Default should be this new
law, not the old one. Let's wait to fly new long tests on the Mac Studio once more changes land on main from other
agents."). Two amendments came with the merge. First, `SkinSettings.rheology` defaults to "composite" in the library as
well as on the command line (section 6 and the first call of section 11 are amended); a state saved before the option
existed still resumes on "bridge", the law it flew, through `SkinSettings.__setstate__`. Second, main's kinked-profile
dispersion (merged there the same day) reads the cells' viscosity through `skin.cell_viscosity`, so its profile path
follows the selected law too. Level 5 waits until more changes land on main. The model_assumptions section is section 13,
since main's profile path took section 12.

The literature behind it lives in `reports/` and `research_notes/` at the repository root (not committed):

- the first check, `Near liquidus aluminium alloy viscosity.md`, and its full-text notes (`fulltext_*.md`);
- the review of stress-controlled branching, `Semisolid flow curve branching.md`;
- the review of the band below 10 % solid, `Aluminium viscosity below ten percent solid.md`.

## 1. Why

### 1.1 The present law

The skins' flow and the exposed slurry's spraying mode take every cell's viscosity from `reentry_model/shared/rheology.py`.
Today that is Li et al.'s (2014) averaged law, eta = [0.871 - 0.00849 gamma^0.74924] exp(3.7311 f_s) Pa s, where gamma
is the shear rate (per second) and f_s the solid fraction. The shear rate is held at no more than 367 per second.
Between 906.4 K, where the model's Scheil curve stands at 10 % solid, and the 908 K liquidus, a log-linear "bridge" in
temperature joins that law to the liquid's 1.3 mPa s. The bridge was labelled "no data" when it was written. It sets
the top cells of a melting skin (2.5 to 5.7 % solid) at 3.7 to 48 times the liquid.

That viscosity bears on spraying only when the slurry may spray. With `--slurry-spray off`, the default, only fully
liquid cells (at or above the 908 K liquidus) can be released, so the new law changes how the skin flows but not by
itself whether it sprays (the builder check of 2026-11-12). The saved d050 flight ran with it off and had no liquid at
178.5-179.5 s.

Asha decided (2026-11-12) that the default stays off in this build. The new spraying law the Girin session is building
will set a new default spraying pattern in which the slurry can spray, and it will take its viscosities from this law.

Asha decided (2026-11-12) to use Li et al.'s law only where it is valid and to build the rest from the literature. She
also asked:

- that every choice be tied to its sources (section 4);
- that behaviour not supported by literature be marked as assumed, with the uncertainties below 10 % solid stated, in
  the documents and the code (section 3.7);
- that no command-line flag be added for those assumptions.

### 1.2 What the two reviews changed

The first version of this design used Li et al.'s per-rate fits inside their measured window, with a non-monotone stress
(a peak near 245 per second and a valley at 490). Below 10 % solid it scaled Li's 10 % excess down as the square of the
solid fraction, and it kept each cell's shear-rate history to choose between a structured branch and a broken-down one.
Two reviews then found the following.

- **Stress never falls with rate within one structure.** No source shows a semi-solid metal whose stress falls as the
  rate rises within one microstructure. Every falling "steady" curve is a set of steady states of different
  structures. Li et al.'s own method confirms this: each rate was a separate slurry solidified while stirred at that
  rate, and their Fig. 11 shows a different grain structure at each.
- **Rate-jump tests rise.** Tests that read the stress before the structure can change show it rising with rate.
- **A falling stretch is unstable and is never occupied.** Under imposed stress, soft materials jump between a
  near-stopped branch and a fast one.
- **At these rates breakdown is effectively instant and rebuilding frozen.** Structure breaks in milliseconds and
  rebuilds over hundreds of seconds, so a layer the gas has sheared is broken down.
- **Below 10 % solid there is no aluminium measurement.** The analogue systems and the nearest aluminium fits behave as
  a dilute suspension when the crystals are compact. Branched crystals behave like a more crowded suspension because
  they trap liquid. Li et al.'s 40 to 750 times the liquid at 10 % solid is about a hundred times Hirai et al.'s
  aluminium-copper at the same fraction, very likely because the two had different structures.

Asha's decisions after the reviews (2026-11-12):

1. **No history.** The two-branch switch and the per-cell shear-rate history are dropped.
2. **Below 10 % solid, a suspension law in the grains' envelope fraction, independent of shear rate.** The internal
   solid fraction is 0.5 by default, as a constant in the code marked as assumed, not a run option. The test brackets it
   with 1 and 0.25.
3. **At 10 % solid and above, Li et al.'s most broken-down structure.** This is their fit at 489.88 per second, used
   at every shear rate. The jump this leaves at 10 % solid is marked as unresolved.

The earlier decisions stand:

- the present law stays selectable as `--rheology bridge`, and a run name marks it when it is not the default;
- the sources are documented;
- the assumptions are marked.

## 2. Goals, non-goals and anti-goals

Goals:

- A new default viscosity law, `composite`, defined for every solid fraction from 0 to just under 0.5. It is
  independent of the shear rate, equals the liquid at zero solid, and equals Li et al.'s 489.88-per-second fit at 10
  to 50 % solid.
- The present law kept as `bridge`. A flight run with `--rheology bridge` gives the same numbers as today to the last
  bit.
- The sources of every choice, and the status of every assumption, written in this document, the module docstring,
  `model_assumptions.md` and the shared-context fact.
- A cheap first comparison: a probe of the saved d050 state's skins on the laptop (section 8, level 4). Resumed flights
  come later, only if Asha decides (level 5).

Non-goals:

- **No history and no structure parameter.** A never-sheared, structured deposit is a separate question for later,
  only if re-frozen deposits turn out to matter. That deposit holds firm until it yields at kilopascals, and the
  branching review's "option two" would represent it with a transported structure parameter.
- **No rate dependence anywhere in the new law.**
- **No change to the momentum solver.** With a rate-independent law, its shear-rate fixed point ends at its first pass.
- **These stay unchanged:** the liquid's viscosity (1.3 mPa s from the material file), the Scheil curve, the rigid
  threshold of half solid, `SLURRY_SIGMA`, and the element model.
- **No wiring into the level-set architecture.** The law sits in the shared package, ready for it.
- **No change to `--slurry-spray`'s default (off)** (Asha, 2026-11-12). The Girin session's new spraying law will set
  the default spraying pattern.
- **No change to the exposed slurry mode's closure** (Asha, 2026-11-12). `skinbody.slurry_modes` gives the whole flowing
  depth, down to the half-solid front, the viscosity of the first non-liquid cell. Under the new law that cell is about
  1.3 times the liquid, while the cells from 10 % solid down are 46 to 168 times, so the mode's surface speed and release
  rate can be far above what the momentum step resolves cell by cell. The probe measures this artefact (section 8).
- **No production flights.**
- **No edits to other branches.** The `kinked-dispersion` branch reads `rheology.BRIDGE_LOW`, so the present law's
  names stay.

Anti-goals:

- The results of `--rheology bridge` changing at all.
- Any part of the law reading as measured when it is assumed.

## 3. The law (decided)

### 3.1 Variables

The law takes the cell's solid fraction f_s, from its temperature through the material's Scheil curve as now, and the
liquid's viscosity mu_l (the material's `liquid.mu`, 1.3 mPa s for AA7075). It does not take the shear rate. Li et al.'s
data are placed by solid fraction, not temperature. Their DSC puts 10 % solid at 903.45 K, and the model's Scheil curve
puts it at 906.36 K.

### 3.2 Below 10 % solid: the grains as a suspension, in their envelope fraction

    eta = mu_l [1 + S / (1/phi - 1/phi_m)],   phi = f_s / f_i,   S = 3,  phi_m = 0.55,  f_i = 0.5,

and eta = mu_l at f_s = 0.

This is the Mori-Ototake law (Moon's form, with S = 3 for spheres) applied to the envelope fraction phi. Each grain is
enclosed in an envelope through its tips. The internal solid fraction f_i is the share of that envelope that is solid,
and the liquid inside moves with the grain, so a suspension law sees phi = f_s / f_i. With f_i = 1 (compact grains) the
law is Moon's on the true solid fraction. f_i must exceed 0.1 / 0.55, about 0.18, so that phi stays below phi_m below
10 % solid.

### 3.3 At 10 % solid and above: Li et al.'s most broken-down structure

    eta = 0.04346 exp(3.23 f_s) Pa s,   0.1 <= f_s < 0.5,  at every shear rate.

This is Li et al.'s Table 3 fit at their highest rate, 489.88 per second, where their slurry had its finest, roundest
grains. Taking it at every rate assumes two things:

- the structure in flight is the most broken-down one;
- a fixed structure's viscosity does not change with rate.

Section 3.7 marks both as assumed. The law gives 46.2 times the liquid at 10 % solid and 168 times at 50 %.

### 3.4 Rigid

At half solid or more the viscosity is infinite, as now (`RIGID_SOLID_FRACTION` = 0.5).

### 3.5 The jump at 10 % solid (unresolved)

At 10 % solid the law jumps from the suspension value to Li et al.'s:

- from 1.94 to 46.2 times the liquid with f_i = 0.5, a factor of 24;
- from 1.37 to 46.2 times with f_i = 1, a factor of 34;
- from 5.40 to 46.2 times with f_i = 0.25, a factor of 8.6.

No data say how the two join. Li et al.'s 7075 at 10 % solid is about a hundred times Hirai et al.'s aluminium-copper,
and the reviews attribute that to different structures. The jump is left as it is and marked unresolved.

On the model's Scheil curve the jump sits at 906.36 K. That is 0.04 K below the lower end of the present law's bridge,
906.4 K, where the curve stands at 9.75 % solid.

### 3.6 Reference values for the tests

The table gives eta / mu_l, computed from sections 3.2 to 3.5 by a scratch script on 2026-11-12. A test should match it
to a relative 1e-3.

| f_s | 0.01 | 0.025 | 0.038 | 0.057 | 0.0975 | just below 0.1 | 0.1 | 0.2 | 0.3 | 0.4 | 0.45 | 0.4999 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| f_i = 0.5 (the law) | 1.062 | 1.165 | 1.265 | 1.431 | 1.906 | 1.943 | 46.18 | 63.78 | 88.10 | 121.7 | 143.0 | 168.0 |
| f_i = 1 (bracket) | 1.031 | 1.079 | 1.122 | 1.191 | 1.356 | 1.367 | 46.18 | (as above) | | | | |
| f_i = 0.25 (bracket) | 1.129 | 1.367 | 1.630 | 2.168 | 5.022 | 5.400 | 46.18 | (as above) | | | | |

In the table, the top cells stay within a factor of 2.2 of the liquid for any internal fraction from 0.25 to 1. Metal
at 10 % solid and above is 46 to 168 times the liquid.

### 3.7 What is assumed, and what is uncertain

Asha asked on 2026-11-12 that these be stated wherever the law is described: the module docstring,
`model_assumptions.md` and the shared-context fact. "Abstract level" below means the review read only the abstract or a
search snippet, as the reviews themselves state for most semi-solid sources.

**At 10 % solid and above.**

1. *Assumed: the structure in flight is Li et al.'s most broken-down one.* This rests on measured breakdown and
   rebuild times in other alloys:
   - breakdown in milliseconds: A357 fell to 6 Pa s in about 4 ms at 800 per second (Yurko and Flemings, cited in
     Kirkwood and Ward 2004);
   - rebuilding over about 500 seconds: A356's fitted rebuild rate is 0.00191 per second (Florides et al. 2024).

   Both are extrapolated to 7075 and to the gas's rates. Metal that has never been sheared (remelting from below, or a
   re-frozen deposit) would be stiffer, and the law does not represent it (section 2).
2. *Assumed: no rate dependence at a fixed structure.* Structural models are Newtonian at fixed structure (Kirkwood and
   Ward 2004, citing Chen and Fan). Rate-jump tests show the stress rising in proportion to the rate or faster
   (Atkinson and Favier 2016, abstract level; Koke and Modigell 2003).
3. *Upper bound above 489.88 per second.* In steady tin-lead data the viscosity was still falling at 2000 per second
   (Turng and Wang, in Fan 2002, Fig. 4). A 7075 layer sheared faster than 490 per second may break down further than
   this law allows.
4. *Li et al.'s lowest values are uncertain.* At 10 % solid and 490 per second their plotted points scatter from 0.05
   to 0.2 Pa s, and their rig was not calibrated with a known liquid.
5. *Different processing.* Li et al.'s slurry was solidified while stirred; the sphere remelts wrought metal.

**Below 10 % solid (no aluminium measurement exists).**

1. *Nothing has been measured.* No rheological measurement of any aluminium alloy between the liquid and 10 % solid
   was found. The nearest are fitted curves, and analogue systems (section 4).
2. *The internal solid fraction is unmeasured for 7075.* 0.5 is the solidification models' default. The only measured
   values, 0.1 to 0.55, are for a transparent analogue, succinonitrile-acetone (Badillo, Ceynar and Beckermann). The
   review infers about 1 for remelting wrought grains and 0.2 to 0.3 for re-frozen deposit. The bracket rows of
   section 3.6 show the effect.
3. *No rate dependence is assumed.* Hirai et al.'s aluminium-10 % copper fits are strongly rate dependent for branched
   crystals: 5.1 times the liquid at 75 per second against 1.22 times at 200 per second, at 5 % solid. These are
   values inferred from fitted curves, not data points. At the gas's rates, of order 10^3 per second and above, the
   sheared values are the relevant ones.
4. *Coherency may come earlier in coarse deposits.* Coherency, where grains first form a connected network, ranges from
   below 10 % solid in poorly grain-refined aluminium alloys to above 50 % in fine-grained ones (Arnberg, Chai and
   Bäckerud 1993; Chai et al. 1995). For 7075 it was measured at 55 to 85 % by a thermal method that reads high (Ahmad,
   Naher and Brabazon 2013). A coarse re-frozen deposit could be coherent below 10 % solid, which this law cannot
   represent.
5. *The jump at 10 % solid is unresolved* (section 3.5).
6. *The layer's position is uncertain by more than its width.* Measured 7075 liquidus temperatures span 903 to 913 K,
   about six times the 1.6 K the model's Scheil curve gives the 0-10 % solid layer.
7. *The liquid itself is uncertain.* Its viscosity is known to 14 % for pure aluminium (Assael et al. 2006), alloying
   moves it by up to 25 % (Wang and Overfelt 2002), and no 7xxx liquid has been measured.
8. *S = 3 comes from rigid particles.* It was measured on silicon carbide in liquid aluminium. Moon's ratios are taken
   against his rig's 0.020 Pa s "liquid", 15-16 times the true value, so they hold only if that error is proportional.
9. *The law is a steady-state law.* The stress just after a sudden jump in rate rises as the rate to the power 2.07
   (Koke and Modigell 2003, Fig. 18), which is not included.

## 4. Sources for each choice

The PDFs of the first check are in `Research/Literature Review/Near Liquidus Viscosity/` (above the repository), with
page-level notes in `research_notes/Near liquidus aluminium alloy viscosity/`. The sources of the two reviews were mostly
read at abstract or snippet level, as their reports state. Their notes are in `research_notes/Semisolid flow curve
branching/` and `research_notes/Aluminium viscosity below ten percent solid/`, with links in the reports.

| Choice | What it takes | Source | Status |
|---|---|---|---|
| Li et al.'s fit at 489.88 per second for 10-50 % solid | A = 0.04346 Pa s, B = 3.23 | Li et al. 2014, China Foundry 11, 79-84: Table 3 (pp. 82-83); method pp. 79-81; Fig. 11 for the structures | Fitted to 7075 data at 10, 20, 30, 40 and 50 % solid; read in full |
| The structure in flight is the most broken-down | That fit at every rate | Breakdown: Yurko and Flemings (A357), via Kirkwood and Ward 2004; Liu et al. 2003 (Sn-15Pb, fast stage under 1 s, depending only on the final rate). Rebuild: Florides et al. 2024 (A356, about 520 s) | Assumed for 7075; measured on other alloys, abstract level except Kirkwood and Ward |
| No rate dependence at fixed structure | eta independent of gamma | Kirkwood and Ward 2004 (structural model Newtonian at fixed lambda); Atkinson and Favier 2016 (rate jumps: Newtonian or thickening); Koke and Modigell 2003 | Assumed; abstract level except Koke and Modigell |
| No history, no falling branch | The switch dropped | Divoux et al. 2016 (falling branch "not found in experiment"); Yerushalmi, Katz and Shinnar 1970 (instability, via Divoux); Coussot et al. 2002 and Møller, Mewis and Bonn 2006 (stress-controlled jumps and history in soft matter); Li et al.'s method (one slurry per rate) | Review consensus in soft matter; no semi-solid metal measurement of a falling within-structure curve found |
| Suspension law below 10 % solid | Mori-Ototake, S = 3, phi_m = 0.55 | Moon 1990 (MIT thesis), pp. 101-102, Eq. 51, Fig. 5.14 (p. 191); Moon, Cornie and Flemings 1991, MSE A144, 253-265 | Measured on SiC in liquid aluminium; read in full |
| Analogues behave as dilute suspensions | The suspension form below 10 % | Ishibashi and Sato (basalt to 31 % crystals); AZ91D Couette (J. Wuhan Univ. Tech. 2015); ice slurries (Int. J. Air-Cond. Refrig. 2023) | Measured on analogue systems, abstract level |
| Envelope fraction phi = f_s / f_i | Morphology through the internal solid fraction | Beckermann and Wang 1997 (JOM); Wu and Ludwig 2009 (Acta Mater. 57, 5621-5631) | The solidification models' construction |
| f_i = 0.5 | Default internal solid fraction | Zheng et al. 2018 (Leoben: 0.5 held until the true fraction overtakes it); measured 0.1-0.55 in succinonitrile-acetone (Badillo, Ceynar and Beckermann) | Assumed for 7075 (Asha's decision 2); brackets 1 and 0.25 in the test |
| The nearest aluminium values below 10 % | Comparison, not used in the law | Hirai et al. 1993, ISIJ Int. 33 (Al-10Cu: 1.08, 1.22, 1.49 times at 2, 5, 10 % solid and 200 per second; 2.5, 5.1, 10.3 times at 75 per second) | Inferred from fitted curves |
| Coherency | The limit of the suspension picture | Arnberg, Chai and Bäckerud 1993; Chai et al. 1995, MMTA 26, 965-970; Malekan and Shabestari 2009 (A319, 16-21 %); Ahmad, Naher and Brabazon 2013 (7075, 55-85 %, thermal) | Measured; abstract level |
| The liquid, 1.3 mPa s | Unchanged | Assael et al. 2006, J. Phys. Chem. Ref. Data 35, 285-300, Eq. 19 (p. 295): 1.42 mPa s at 908 K, 14 % | Pure aluminium, extrapolated 25 K; read in full |
| Rigid at half solid | Infinite | Unchanged from the layered-skin design (spec 6.4) | As before |
| Solid fraction carries Li et al.'s data | Their 10 % at the model's 906.36 K | Li et al. 2014, Table 2 (pp. 80-81) | Builder's call |

## 5. No history, and what that leaves of the solver (decision 1)

No cell carries a shear rate or a structure parameter. The reasons, from the branching review:

- no within-structure falling stress has been measured;
- a falling stretch is unstable and cannot be occupied;
- breakdown at these rates takes milliseconds and rebuilding hundreds of seconds.

The skins' momentum step is unchanged. Its shear-rate fixed point evaluates `cell_viscosity` at each pass's rates, and
under `composite` the viscosity does not depend on the rate. The first pass therefore moves no viscosity, and the fixed
point ends there. The exposed slurry's mode (`skinbody.slurry_modes`) likewise takes the law's viscosity at the cell's
solid fraction. Its fixed point ends at once, and under `composite` the rate it records is the layer's own, not clipped
at `GAMMA_MAX`.

Under `bridge`, everything is today's code path, unchanged.

## 6. Interfaces and options

`reentry_model/shared/rheology.py`:

- The present names keep their meaning and numbers: `li_viscosity`, `viscosity`, `LI_A`, `LI_B`, `LI_N`, `LI_C`,
  `GAMMA_MAX`, `BRIDGE_LOW`, `BRIDGE_HIGH`, `RIGID_SOLID_FRACTION`. The `kinked-dispersion` branch and the tests use them.
- New constants:
  - `LI_BROKEN_RATE` = 489.88, `LI_BROKEN_A` = 0.04346 and `LI_BROKEN_B` = 3.23;
  - `LI_SOLID_MIN` = 0.1, `SUSPENSION_S` = 3.0 and `SUSPENSION_PHI_M` = 0.55;
  - `INTERNAL_SOLID_FRACTION` = 0.5, with a comment saying it is assumed;
  - `LAWS` = ("composite", "bridge") and `DEFAULT_LAW` = "composite".
- New functions, both vectorised:
  - `suspension_viscosity(f_s, mu_liquid, internal_fraction=None)` gives section 3.2's law for any f_s from 0 to the
    envelope's packing limit.
  - `composite_viscosity(f_s, mu_liquid, internal_fraction=None)` gives sections 3.2 to 3.4. It is infinite at half
    solid or more.

  `internal_fraction=None` means `INTERNAL_SOLID_FRACTION`, read at call time, so the probe can set the brackets in its
  own process. A value at or below 0.1 / 0.55 raises ValueError.
- The choice between the laws also lives here, so the level-set architecture can use the same law through the same
  call (Asha, 2026-11-12):
  - `law_viscosity(law, T, f_s, gamma, mu_liquid, f_s_bridge)`: for `"bridge"`, exactly today's `viscosity(T, f_s,
    gamma, mu_liquid, f_s_bridge)`; for `"composite"`, `composite_viscosity(f_s, mu_liquid)` broadcast to `T`'s shape,
    ignoring `gamma`; any other name raises ValueError.
  - `bridge_solid_fraction(material)`: the material's solid fraction at `BRIDGE_LOW`, as both skyscrapers call sites
    compute it today. It reads only `material.liquid_fraction`, so the module still imports nothing.
- The module docstring gives each region with its sources and statuses (section 4) and section 3.7 in full.
- The module imports nothing from the rest of `reentry_model`.

`reentry_model/skyscrapers/skin.py` and `skinbody.py`:

- `SkinSettings` gains `rheology`. Its dataclass default is "bridge", so that existing tests and states saved before
  this change keep the present law: an unpickled old settings object has no `rheology` and reads the class default. The
  run command's default is "composite" (`DEFAULT_LAW`).
  *Amendment at the merge into main (2026-11-13; Asha, 2026-11-12):* the dataclass default is "composite" (`DEFAULT_LAW`) too; an
  unpickled old settings object with no `rheology` is given "bridge" by `SkinSettings.__setstate__`, and the existing
  tests whose values were set on the bridge law name it.
- A thin skins wrapper passes `sk.settings.rheology`, the cell state and the material's values to
  `rheology.law_viscosity`, and holds no physics of its own. Both `cell_viscosity` and `slurry_modes` use it, so the
  slurry mode follows the law too. `slurry_modes` today calls `rheology.viscosity` directly.

Command line (`cli.py`):

- `run --rheology composite|bridge` (default composite) applies to melting skin runs. With `--surface-model elements`
  or without `--melt on` it is accepted and ignored.
- `resume --rheology` defaults to the saved flight's law, which is `bridge` for a state saved before the option
  existed. A change is recorded in `changed`.
- Run names: `skin_name_suffix` adds `_rheology-bridge` for the present law; the default is unmarked. The run's JSON
  always records `rheology` and `internal_solid_fraction` under `skin`, so runs made before and after the change can be
  told apart.
- `README.md`, `CLAUDE.md` and the usage text list the option.

An analysis script, `Reentry Model/analysis/near_liquidus_viscosity.py`, draws the law against Li et al.'s points and
the present law, prints section 3.6's table, and marks the assumed regions. It calls the package's functions and copies
none of them.

The `kinked-dispersion` branch needs no change for this law. Its iteration of each cell's viscosity from rest gives
the same values, since the law ignores the rate. Its kink at the 906.4 K isotherm lies 0.04 K above the new law's jump at
10 % solid (906.36 K) on the model's Scheil curve. That session should still be told.

## 7. What the change should do (to be measured in the probe)

The table gives each cell's viscosity as a multiple of the liquid's, at any shear rate under the new law. The present
law depends on the rate, so its column shows the range over the rates the solver uses, at most 367 per second.

| Cell | New law, f_i = 0.5 | f_i = 1 | f_i = 0.25 | Present law |
|---|---|---|---|---|
| 2.5 % solid (907.6 K) | 1.17 | 1.08 | 1.37 | 3.7-5.6 |
| 5.7 % solid (907.1 K) | 1.43 | 1.19 | 2.17 | 18.5-48 |
| 9.75 % solid | 1.91 | 1.36 | 5.02 | 180-964 |
| 10 % solid and above | 46-168 | (same) | (same) | 181 and more |

Under the median gas stress of 101 Pa, a cell at 10 % solid now shears at about 1700 per second.

The new law makes the top cells near-liquid for any internal fraction tested, and the slurry from 10 % solid upward
4 to 26 times runnier than the present law: 46 against 181-973 times at 10 % solid, and 168 against 806-4328 times at
50 % (the present law's `li_viscosity` at 367 and at 0 per second). Faster layers raise the adaptive count of flow sub-steps, so a skin
flight may cost more per macro step. The stress a cell carries is the gas shear plus the push on the metal above it.

## 8. Tests and acceptance

Level 1, the law (`tests/test_reentry_model_rheology.py`):

- section 3.6's table to 1e-3;
- Li et al.'s 489.88 fit reproduced exactly at 10, 20, 30, 40 and 49.99 % solid;
- the liquid at zero solid and below;
- infinite at half solid and above;
- the jump at 10 % solid as section 3.5 gives it;
- the same values at any shear rate through the skins' helper;
- an internal fraction at or below 0.1 / 0.55 rejected;
- the present law's existing tests unchanged.

Level 2, the skins:

- `cell_viscosity` and `slurry_modes` follow the law under `composite`;
- under `composite`, the momentum step's fixed point ends after its first pass, and its books close to round-off;
- under `bridge`, the existing skin and skinbody tests pass unchanged.

Level 3, the command line:

- the option, the run-name mark, the JSON record, and resume's defaults, including a state saved before the change
  loading as `bridge`;
- the regression: the short melting smoke flight, with and without `--slurry-spray on`, run on the old code (the
  `skyscrapers-fixes` worktree, f1ab736) and on the new code with `--rheology bridge`. The results are identical except
  for the two new JSON keys.
- the unit tier passes.

Level 4, a probe of the skins alone on the laptop (Asha, 2026-11-12):

- **The state:** `d050__handdown-off.state.pkl` at 178.5 s, in the subplan20 worktree under
  `prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/`.
- **Why not a resumed flight:** the state was flown with slurry spraying off and shows no liquid at 178.5-179.5 s, and
  a macro step there took 63 to 76 minutes under the present law.
- **The probe:**
  - Run only the skins' flow sub-steps for about 0.1 s, with the gas loads held at the saved step's. Run no heat
    sub-step, no 3D conduction and no trajectory step.
  - Run four cases: `bridge`; `composite` with f_i = 0.5; with f_i = 1; with f_i = 0.25.
  - At several times, evaluate the liquid's and the exposed slurry's spraying modes without releasing metal, under
    two sets of release rules:
    - the saved flight's own: release gates "before", slip band off;
    - the rules of the recent d050 work: gates 2026-11-10, slip band on, as in the kinked-dispersion resumes and the
      Mac Studio flight `d050__fixes__slipband-on`.

    The builder check of 2026-11-12 explains why both are needed. At 178.5 s the body's Knudsen number is 0.056, inside
    the slip band (0.01 to 0.1), and every patch is on the Couette closure. Under the saved rules, Girin's film depth is
    therefore undefined everywhere, both modes fall on the thin branch, and the critical Weber number is not the gate
    there. The probe records which gate decided each mode.
  - The closure artefact is measured (Asha's decision, 2026-11-12). Beside each slurry mode as the code computes it,
    the probe reports two things:
    - the momentum step's own top-cell speeds;
    - the slurry mode recomputed with its layer ended at the 10 % solid depth (the 906.36 K isotherm, by
      `Skins.depth_above`). This is done in the probe, without changing the model.
  - The gas loads, for the flow and for the modes, come from the surface flow recomputed with the current code at the
    saved step's aero state. The slip band needs the shear integral, which the saved result lacks. The probe reports
    the recomputed loads' difference from the saved forcing.
- **The report**, `output/near_liquidus_viscosity/probe.md`, gives for each case:
  - the top cells' viscosity ratios;
  - the stresses the cells carry;
  - the surface speeds;
  - the Weber numbers against the critical value, and the release rates the modes would give;
  - the flow sub-step counts and the wall time.

  A single column, the windward column with the median flowing depth, is also run to a steady state under the four
  cases.

Level 5, resumed flights (only if Asha decides after level 4): resume the same state under both laws with
`--slurry-spray on`, to 179.5 s, on the Mac Studio by a runbook. *Asha, 2026-11-12:* wait until more changes from other
agents land on main before flying long tests on the Mac Studio.

## 9. Documentation

- `docs/model_assumptions.md`: a new section on the near-liquidus viscosity (shared physics). It gives the law region by
  region, section 4's table, section 3.7 in full, and the jump at 10 % solid marked unresolved.
- The layered-skin design (`2026-10-11-layered-skin-moving-mesh-design.md`): a dated amendment where it states the
  present law. This covers the table near line 67, section 6.4 near line 343, the slurry mode near lines 611-617, line
  958, and assumption 5 near line 1082, which describes the falling stress and the fixed point's cap.
- `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`: a new fact (196 was free on 2026-11-12). It
  records the change, the regression, the unit tier and the probe's results.
- `README.md` and `CLAUDE.md`: the option and the changed default.
- The `kinked-dispersion` session is told (section 6).

## 10. Limitations

- Everything in section 3.7.
- The law describes metal that the gas has sheared. A never-sheared structured deposit, which would hold until it yields
  at kilopascals, is not represented.
- The exposed slurry's mode gives its whole flowing depth the top slurry cell's viscosity (section 2). Under the new law
  this can overstate the slurry's speed and release. The probe measures by how much.

## 11. Left to the builder, and the calls made in this design

Left to the builder:

- the helper's name and where it sits;
- the internal structure of the tests.

Calls in this design that Asha did not decide (she should check these):

- **The dataclass default.** `SkinSettings.rheology` defaults to "bridge" in the library, and the run command defaults
  to "composite" (section 6). *Decided by Asha at the merge, 2026-11-12:* "composite" in the library too.
- **The variable.** Solid fraction, not temperature, places Li et al.'s data (section 3.1).
- **The validation limit.** An internal fraction at or below 0.1 / 0.55 is rejected (section 6).
- **The default run name.** It is unmarked, and the JSON records the law and the internal fraction.
- **Resume.** It defaults to `bridge` for states saved before the option existed.
- **The slurry mode's recorded rate.** Under `composite` it is not clipped at 367 per second.
- **The probe's release rules and gas loads.** The probe evaluates the spraying modes under the saved rules and under
  the 2026-11-10 gates with the slip band on, both on gas loads recomputed with the current code (section 8).
- **The `--slurry-spray` help text.** Today it says the slurry sprays "on a viscosity law with no data there"; it is
  updated to name the law in use.
