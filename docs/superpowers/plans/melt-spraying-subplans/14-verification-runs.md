# Sub-plan: Task 14 — Verification and sensitivity drivers, the reference-tier test, the runs

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 7194–7488). Read `00-shared-context.md` first. **Amended 2026-10-02: the deep runoff and the per-patch conjugate depth** (section below). **Amended 2026-10-03: the molten cascade** (the section after it). **Amended 2026-10-05: the seed of numpy's generator** (the section after that). **Amended 2026-10-05 (runoff flux): the verification and sensitivity runs, and the switched-step study** (the section after the seed's). **Amended 2026-10-07: the drivers on Scheil's curve with the rigid substrate, and the switched-step series repeated** (the section after the runoff flux's). **Amended 2026-10-07 (continuum step): the drivers switch by default, and the harness retires** (the section after that).


> **Amended 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md`
> (measured facts 38–45 in `00-shared-context.md`). The block below takes precedence over the extracted body.

**Changes to Task 14.** Every number in the verification table was measured with 4 prism layers and
`PHI_DEATH = 0.05`; both change, so the table must be re-run rather than edited.

1. **Re-measure the whole table** with the new defaults. Fact 12's layer sensitivity (sprayed mass within 0.2 %,
   demise altitude within 0.3 % across layers 2/4/6) predicts the move is small, but predicting is not measuring.
2. **The shape A/B**, which is the measurement that says what this whole design bought: the 50 mm physics flight
   with `--derived-surface off` against `on`, comparing sprayed mass, demise time and median droplet radius.
3. **The `PHI_DEATH` A/B**: 0.05 against 0.50 on the same flight, with `film_blob_fraction`,
   `film_frozen_fraction` and the stranded-film mass reported against fact 25's ranges.
4. **The remesh conservation record**: for every remesh in a flight, the pre-correction volume and energy residuals.
   Spec §10.2 — these replace the unqualified 10⁻¹⁰ balance claim of facts 28 and 31 with "exact after an explicit
   correction, pre-correction residual X".
5. **The first number to obtain, before anything else** (fact 43): how smooth the real per-patch recession field is.
   The least-squares fit is excellent at 0.03 % for a smooth field and mediocre at 22 % for a patchy one, and which
   regime the flight sits in is currently unknown.
6. **The melting references KEEP the band — re-measure them and update the README** (decided 2026-09-28).
   `tests/test_reentry_model_reference_melt.py:24` calls `mesh.sphere_mesh(ref.diameter / 2.0)` with no band, so the
   **bookkeeping device** takes `DEFAULT_BAND` and meshes the way production melting runs mesh. That is deliberate
   and it is the decision: a verification device should exercise the configuration actually run, and the band exists
   for melting runs. **Do not pin it with `band=0.0`.** Its Step 2 sibling, `test_reentry_model_reference_thermal.py`,
   *is* pinned, because a no-melt run has no business acquiring a melting-run mesh — the asymmetry is intended.

   Consequences this task must carry out, not merely note:

   * **Re-measure the six committed bookkeeping figures** and update them wherever they appear — measured fact 12 in
     `00-shared-context.md`, the README verification table, and the thresholds in
     `tests/test_reentry_model_reference*.py`. The current values were measured on the old mesh and are superseded:
     100 mm mass within 0.98 % of m₀, onset +0.10 km, 1 %-mass time −1.3 %; 50 mm 1.29 %, +0.17 km, −0.2 %; against
     thresholds 2 % / 0.5 km / 2 %. Fact 12's "no later amendment can move them" no longer holds for this device and
     that sentence must be amended too, with this decision as the reason.
   * **Expect the numbers to move very little, and treat it as a finding if they move a lot.** The device multiplies
     the conductivity by 10⁴ (`mat.k_table * 1e4`), so the body is isothermal and the mesh has almost nothing left to
     influence. If the six figures shift materially, that says the bookkeeping device is more mesh-sensitive than its
     design implies — which is worth knowing before the thesis leans on it.
   * **Record the runtime cost.** The 100 mm reference mesh goes from **18 896 nodes / 87 632 tets** to
     **33 355 / 177 363**, a factor of 2.02 in elements (measured 2026-09-27). The element-proportional part of each
     macro step — principally the conduction solve — roughly doubles for the two melting flights. Time the melting
     reference tier before and after and write both figures into `CLAUDE.md`, whose "~15 min (Step 1) + ~35 min
     (Step 2)" line does not yet cover the melting references at all.

7. **The better Design A trigger** (spec §10.3): the **melt-rate** error under a surface pinned by the latent heat.
   While melting, a temperature error becomes a mass-flux error rather than a temperature error, and mass flux is
   what feeds the spraying. The 1-D probe behind fact 44 cannot measure it; this task can.
---

## Amendment of 2026-10-02 — the deep runoff in the verification and sensitivity runs

> Part of the deep-runoff amendment (sub-plan 09's amendment of this date; facts 46–53 in `00-shared-context.md`).

1. **The verification devices are unaffected, by construction.** `--removal instant` evaluates no surface flow and so
   has no conjugate depth, and the bookkeeping device also runs `--runoff off`, so the deep stage never runs in either.
   Their six figures and thresholds stand as this sub-plan's amendment of 2026-09-27 left them.
2. **The 50 mm physics flight is unaffected, and that is measured, not assumed.** It never has Girin's closure (the
   wall-Knudsen gate denies it throughout, fact 32), so it has no conjugate depth and no liquid below one; with numpy's
   random seed fixed in both runs the amended and unamended prototypes give it bit for bit (fact 49).
3. **The 100 mm physics flight moves, and its row must be re-measured with the deep runoff on** (fact 49): at the
   default step it sprays 12 % more and lands with 22.4 % instead of 30.8 % of its initial mass — but fact 50 shows that
   this is a time-step artefact, so item 5 comes first. Re-measure the deep pile-up of fact 47 after Task 16 as well:
   the derived surface takes θ and the film tangents from smoothed normals, which is expected to remove most of the
   element-death craters in which the deep liquid collects at the default step, and the size of that change is
   unknown.
4. **A sensitivity row, `nodeep` (`--deep-runoff off`)**, in `analysis/melt_sensitivity.py`, which bounds the new
   mechanism the way the `gammapm14` and `knbody003` rows bound the amendments of 2026-09-22:

```diff
--- a/analysis/melt_sensitivity.py
+++ b/analysis/melt_sensitivity.py
@@ -2,11 +2,12 @@
 """Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).
 
     "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
-        [--variants base,layers2,layers6,dt025,bridged,norunoff,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]
+        [--variants base,layers2,layers6,dt025,bridged,norunoff,nodeep,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]
 
 Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
 six layers from 0.125 mm, 15.9 mm in all -- eight layers of 0.25 mm with growth 2 would exceed the radius; `gammapm14` and
-`knbody003` bound the two modelling choices of the 2026-09-22 amendment, the Prandtl-Meyer gamma and the body gate); the table
+`knbody003` bound the two modelling choices of the 2026-09-22 amendment, the Prandtl-Meyer gamma and the body gate; `nodeep`
+turns off the deep runoff of the 2026-10-02 amendment, the liquid below the conjugate depth then staying in its elements); the table
 lists sprayed mass, median droplet radius, melt-onset, spraying-onset and demise altitudes and the runtime, with the
 change relative to `base`. Every run is a subprocess of `--python` (default: this interpreter); the `fenicsx` variant
 needs the fenicsx_env interpreter (run it separately with --variants fenicsx --python <fenicsx_env python>, with CC
@@ -24,7 +25,8 @@
 CASES = {"d100": ["--diameter", "100", "--altitude", "77.500133"], "d050": ["--diameter", "50", "--altitude", "115"]}
 VARIANTS = {
     "base": [], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
-    "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "we308": ["--we-critical", "3.08"],
+    "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "nodeep": ["--deep-runoff", "off"],
+    "we308": ["--we-critical", "3.08"],
     "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
     "AA7075": ["--material", "AA7075"], "sizeinitial": ["--size-feedback", "initial"], "gammapm14": ["--gamma-pm", "1.4"],
     "knbody003": ["--kn-body-shock", "0.003"], "fenicsx": ["--thermal-solver", "fenicsx"],
```

   Checked: `argv_for("d050", "nodeep", ...)` parses to `deep_runoff == "off"` with the run name `d050__nodeep`
   (17 variants). The runs themselves are Task 14's to make.
5. **A time-step study, before any droplet-population or deep-runoff result is quoted.** Fact 50 has it for the
   100 mm physics flight to 120 s at 0.5, 0.25 and 0.125 s, with `--deep-runoff on` and `off`, all runs seeded: the
   liquid the unamended model holds below δ_m vanishes as the step shrinks (a median 4.75, 0.23 and 0.004 g) and the
   deep runoff's effect with it (+7.1 %, −0.3 % and −0.16 % in sprayed mass), but the molten layer the branch test
   reads shrinks with the step as well, so the unamended model's droplet population has not converged even at
   0.125 s (droplet count 1.72e7, 2.82e7 and 5.85e7; median radius by number 180, 105 and 76 µm). Repeat the study
   on the prototype as Task 14 runs it (sub-plan 02's materials, Task 16's derived surface), for the whole flight and
   for the 50 mm sphere, and again after fact 53 (a)'s option (3) if it is taken; run the `dt025` sensitivity row
   with both settings; and report the held liquid (`deep_liquid_kg` with the flag on, or the contiguous molten layer
   `molten_depth_*` with it off), the deep runoff's effect and the unamended model's own change at each step. The
   cost is proportional to the number of steps: 10, 20 and 40 minutes for 120 s of the 100 mm flight at the three
   steps.
6. **Compare runs as seeded pairs.** Two runs of one build differ from the first melting step on, because pyamg draws
   random starting vectors from numpy's global generator (fact 52); the 50 mm flight's collapse can amplify that to
   per-cent level in the late-flight columns. Any before/after comparison this task reports should fix
   `np.random.seed` at the start of both runs (fact 52 has the measurement and a one-line fix), or else quote the
   run-to-run spread beside the difference.

## Amendment of 2026-10-03 — the molten cascade in the verification and sensitivity runs

> Part of the molten-cascade amendment (sub-plan 09's amendment of this date; facts 54–61 in `00-shared-context.md`).

1. **The verification devices are unaffected, by construction.** `--removal instant` feeds every element as it melts,
   so no fully molten element ever waits to be exposed and the cascade is skipped there; the six bookkeeping figures
   and their thresholds stand as this sub-plan's amendment of 2026-09-27 left them.
2. **The 50 mm physics flight changes, and must be re-measured with the cascade on** (fact 60): at the default step it
   demises 8.0 s earlier and 2.5 km higher (199.0 s and 70.7 km against 207.0 s and 68.2 km), sprays 1.0 % more, and
   makes 32 % fewer droplets of a 17 % larger median radius by number. Its row in the verification table and the
   `nocascade` sensitivity row must be re-measured, and the flight run at a smaller step, with the cascade on, before
   its droplet population is quoted.
3. **The 100 mm physics flight changes** (fact 57): to 120 s at the default step it sprays 3.6 % more (1.055
   against 1.018 kg, deep runoff off) and is 8 % lighter at 120 s (0.417 against 0.454 kg). The whole flight to the
   ground with the cascade was **not** run (fact 60: the machine was on battery); run it before quoting fact 24's
   headline (how much of the 100 mm sphere reaches the ground) again, with the deep runoff on and off.
4. **A sensitivity row, `nocascade` (`--molten-cascade off`)**, in `analysis/melt_sensitivity.py`, beside `nodeep`:

```diff
--- a/analysis/melt_sensitivity.py
+++ b/analysis/melt_sensitivity.py
@@ -2,12 +2,14 @@
 """Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).
 
     "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
-        [--variants base,layers2,layers6,dt025,bridged,norunoff,nodeep,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]
+        [--variants base,layers2,layers6,dt025,bridged,norunoff,nodeep,nocascade,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]
 
 Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
 six layers from 0.125 mm, 15.9 mm in all -- eight layers of 0.25 mm with growth 2 would exceed the radius; `gammapm14` and
 `knbody003` bound the two modelling choices of the 2026-09-22 amendment, the Prandtl-Meyer gamma and the body gate; `nodeep`
-turns off the deep runoff of the 2026-10-02 amendment, the liquid below the conjugate depth then staying in its elements); the table
+turns off the deep runoff of the 2026-10-02 amendment, the liquid below the conjugate depth then staying in its elements;
+`nocascade` turns off the molten cascade of the 2026-10-03 amendment, the surface then receding through molten
+material by one element per macro step); the table
 lists sprayed mass, median droplet radius, melt-onset, spraying-onset and demise altitudes and the runtime, with the
 change relative to `base`. Every run is a subprocess of `--python` (default: this interpreter); the `fenicsx` variant
 needs the fenicsx_env interpreter (run it separately with --variants fenicsx --python <fenicsx_env python>, with CC
@@ -26,7 +28,7 @@
 VARIANTS = {
     "base": [], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
     "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "nodeep": ["--deep-runoff", "off"],
-    "we308": ["--we-critical", "3.08"],
+    "nocascade": ["--molten-cascade", "off"], "we308": ["--we-critical", "3.08"],
     "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
     "AA7075": ["--material", "AA7075"], "sizeinitial": ["--size-feedback", "initial"], "gammapm14": ["--gamma-pm", "1.4"],
     "knbody003": ["--kn-body-shock", "0.003"], "fenicsx": ["--thermal-solver", "fenicsx"],
```

   Checked: `argv_for("d050", "nocascade", ...)` parses to `molten_cascade == "off"` with the run name
   `d050__nocascade` (18 variants). The runs themselves are Task 14's to make.
5. **The time-step study of item 5 of the 2026-10-02 amendment, begun again with the cascade** (facts 57–58): the 100 mm
   physics flight to 120 s at 0.5, 0.25 and 0.125 s, with the cascade on and off and the deep runoff on and off, every
   run seeded. Done so far: the cascade at 0.5 s with the deep runoff on and off, against the seeded runs of 2026-10-02
   at all three steps; **not** done, because the machine ran on battery (fact 60): the cascade at 0.25 and 0.125 s.
   Measured: the cascade removes the backlog (liquid held below the conjugate depth a median 1.23 g instead of 4.75 g at
   0.5 s) but not the droplet population's step dependence (18.1 million droplets at 0.5 s with the cascade, 58.5
   million at 0.125 s without it), because the branch test's layer counts one molten wall-owning element as 0.67–1.16 mm
   of liquid on this mesh, above the conjugate depth on every patch (fact 58). Make the two missing pairs first
   (`--dt 0.25` and `--dt 0.125`, `--deep-runoff on` and `off`, 120 s, seeded; about 30 and 45 CPU-minutes each on mains
   power), then the whole 100 mm flight and the 50 mm flight with the cascade. Repeat it on the prototype as Task 14
   runs it (sub-plan 02's materials, Task 16's derived surface and fact 44's `PHI_DEATH`), because fact 44's
   `PHI_DEATH = 0.50` changes how long a partly molten wall-owning element survives, which is half of fact 58's
   mechanism, and the derived surface changes which facets are wall patches at all; and repeat it once more after fact
   61 (a)'s decision on the branch test's liquid depth, which is now expected to decide whether the droplet population
   converges.
6. **Compare runs as seeded pairs** (item 6 of the 2026-10-02 amendment) still applies; every number of this amendment
   was measured that way, with numpy's generator seeded in the measurement harness and not in the model.

## Amendment of 2026-10-05 — the seed in the verification and sensitivity runs, and the scatter to quote results against

> Part of the seeding amendment (facts 62–68 in `00-shared-context.md`; sub-plan 13's amendment of this date holds the
> design and the CLI code).

1. **Seeded pairs are now automatic.** Item 6 of the 2026-10-02 amendment and item 6 of the 2026-10-03 amendment
   ("compare runs as seeded pairs") are met by the model itself: every run seeds numpy's generator at the start of
   `cli.cmd_run`, with 12345 by default — the seed the measurement harness of those amendments set — so two runs of one
   build agree bit for bit, and the seeded runs of facts 49–61 can be reproduced by the command line alone, with the
   flags that select the model each was made with (fact 64: verified on three of them). No measurement harness needs to
   seed numpy any more. What a before/after comparison must still do is quote its
   difference against the scatter between seeds (fact 65, item 5 below), not against zero: two runs that differ only in
   the seed land that far apart, so a smaller difference between two settings is not a result.
2. **The drivers are seeded by construction.** `analysis/melt_verification.py` and Step 2's
   `analysis/reentry_model_thermal_verification.py` call `cli.main` in-process and `analysis/melt_sensitivity.py` runs
   the CLI as subprocesses, so every run they make starts from the default seed, whatever ran before it in the same
   process; the verification and sensitivity tables are reproducible to the last digit. Their summary tables list
   results by case and mode or variant, not the configuration, so they gain no seed column.
3. **The reference-tier tests are not seeded, by decision.** `tests/test_reentry_model_reference_melt.py` and Step 2's
   `tests/test_reentry_model_reference_thermal.py` construct their runs directly rather than through `cmd_run`, so they
   draw from whatever state the generator is in. Measured on the bookkeeping device through the CLI (the 100 mm
   reference's initial state, SESAM-equivalent heating, `AA7075`, instant removal, `--runoff off`, `--k-scale 1e4`, no
   prism layers, 132 macro steps to demise): seed 12345 against seed 1 moves no history column by more than 2.4e-9 of
   its value, the body mass by at most 1.4e-11 of the initial mass, the melt onset (71.10 km at 43.0 s) and the demise
   (66.0 s) not at all, and the integrated heat by 1.3e-12 — nine to ten orders of magnitude below the device's
   thresholds of 2 % of the mass, 0.5 km and 2 % of the 1 %-mass time, and far below the digits the README prints. The
   Step 2 runs are smaller still (a 5 s run on the coarse mesh moves by 2e-16 to 6e-16, fact 62). Seeding them would
   make the metrics files they write repeat to the last bit, at the cost of a test-file change that cannot be run in the
   prototype until Task 11 commits the melting references; fact 68 (a) leaves it to Asha.
4. **A sensitivity row, `seed1` (`--seed 1`)**, in `analysis/melt_sensitivity.py`, beside `base`: the same flight with
   another seed, so that the table carries its own floor and every other row's change can be read against it.

```diff
--- a/analysis/melt_sensitivity.py
+++ b/analysis/melt_sensitivity.py
@@ -2,14 +2,16 @@
 """Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).
 
     "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
-        [--variants base,layers2,layers6,dt025,bridged,norunoff,nodeep,nocascade,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]
+        [--variants base,seed1,layers2,layers6,dt025,bridged,norunoff,nodeep,nocascade,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]
 
 Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
 six layers from 0.125 mm, 15.9 mm in all -- eight layers of 0.25 mm with growth 2 would exceed the radius; `gammapm14` and
 `knbody003` bound the two modelling choices of the 2026-09-22 amendment, the Prandtl-Meyer gamma and the body gate; `nodeep`
 turns off the deep runoff of the 2026-10-02 amendment, the liquid below the conjugate depth then staying in its elements;
 `nocascade` turns off the molten cascade of the 2026-10-03 amendment, the surface then receding through molten
-material by one element per macro step); the table
+material by one element per macro step; `seed1` reruns `base` with another seed of numpy's generator -- every run is
+seeded, `base` with the default 12345 (amendment of 2026-10-05) -- so its change relative to `base` is the run-to-run
+scatter, the floor against which every other row is read); the table
 lists sprayed mass, median droplet radius, melt-onset, spraying-onset and demise altitudes and the runtime, with the
 change relative to `base`. Every run is a subprocess of `--python` (default: this interpreter); the `fenicsx` variant
 needs the fenicsx_env interpreter (run it separately with --variants fenicsx --python <fenicsx_env python>, with CC
@@ -26,7 +28,7 @@
 
 CASES = {"d100": ["--diameter", "100", "--altitude", "77.500133"], "d050": ["--diameter", "50", "--altitude", "115"]}
 VARIANTS = {
-    "base": [], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
+    "base": [], "seed1": ["--seed", "1"], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
     "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "nodeep": ["--deep-runoff", "off"],
     "nocascade": ["--molten-cascade", "off"], "we308": ["--we-critical", "3.08"],
     "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
```

   Checked: `argv_for("d050", "seed1", ...)` parses to `seed == 1` with the run name `d050__seed1`, and `base` to the
   default 12345 (19 variants). The runs themselves are Task 14's to make.
5. **The scatter to quote results against** (fact 65). Measured on the 100 mm physics flight to 120 s at the default
   step, with every setting at its default and the seeds 12345, 1, 2 and 3: the range across the four runs, as a share
   of their mean, is 0.13 % in sprayed mass (1.0600 to 1.0614 kg), 0.33 % in the mass at 120 s (0.4107 to 0.4120 kg),
   2.0 % in droplet count (1.76e7 to 1.80e7), 0.29 % in the median radius by number (179.6 to 180.1 µm) and 1.4 % in
   re-solidified mass (6.38 to 6.46 g), with 6.7 % in the front-surface Rayleigh–Taylor release; the 50 mm whole flight,
   with one pair of seeds, moves by at most 2.5e-6 of any history value. Fact 52's figures (0.15 % in sprayed mass,
   0.33 % in the mass at 120 s, 1.2 % in droplet count and 14 % in re-solidified mass, from two unseeded runs of the
   model before the molten cascade) are superseded by these as the floor. Re-measure it on the prototype as Task 14 runs
   it (sub-plan 02's materials, Task 16's derived surface, fact 44's `PHI_DEATH`), for the 50 mm flight as well, and at
   whatever step the results are quoted at, since the amplification that makes it comes from the melt step's thresholds
   and will move with them; with four seeds it is a range, not a distribution (fact 68 (b)).
6. **The `fenicsx` row is seed-independent.** The FEniCSx backend draws nothing from numpy's generator and two of its
   runs with different seeds are identical (fact 66), so its difference from `base` contains `base`'s own scatter; read
   it against the `seed1` row, not against zero.

## Amendment of 2026-10-05 (runoff flux) — the verification and sensitivity runs, and the switched-step study

> Part of the runoff-flux amendment (sub-plan 06's amendment of this date holds the change; facts 69–77 in
> `00-shared-context.md`). No driver or sensitivity code changes: the change has no flag, so there is no new row.

1. **The bookkeeping devices are unaffected, by construction.** `--removal instant` evaluates no surface flow and so
   moves no film, and the bookkeeping device also runs `--runoff off`; their six figures and thresholds stand as this
   sub-plan's amendment of 2026-09-27 left them. The **resolved mode** (SESAM-equivalent heating with Girin's removal,
   D₀ and R₀) does transport a film, and where the flow branch is the shock layer it has a conjugate depth, so its two
   rows can move: re-measure them (not run here; the comparison needs the melting references of Task 11).
2. **The 50 mm physics flight is unaffected, measured rather than assumed.** It never has Girin's closure, so no patch
   is thick and the change never acts: the whole flight is bit-identical with and without it (fact 73).
3. **The 100 mm physics flight moves from 49.5 s, its first step under Girin's closure** (fact 73): at the default step,
   to 120 s, its masses by about the scatter between seeds of fact 65 (sprayed −0.14 % against 0.13 %, mass at 120 s
   +0.37 % against 0.33 %), its droplet count and median radius by number within it, and its front-surface
   Rayleigh–Taylor release (−16 %) and re-solidified mass (+3.7 %) outside it. Its row in the verification table and
   every sensitivity row of the 100 mm case are to be re-measured with the change.
4. **The switched-step study** (facts 69, 74 and 75), the time-step study of items 5 of the amendments of 2026-10-02 and
   2026-10-03 continued in the continuum regime, where the step matters. It was made with a measurement harness outside
   the package (`dtswitch.py`, which shortens the macro step once the body Knudsen number falls below the model's own
   continuum boundary of 0.01, and `compare.py`, which reports only measures that do not depend on the step); neither
   belongs in `analysis/`, because the switch is not a model option (fact 77 (b) asks whether it should become one).
   Measured with the change: the masses settle within their scatter from 0.05 s down (sprayed mass 1.050, 1.052, 1.032
   and 1.035 kg at 0.05, 0.025, 0.0125 and 0.00625 s, against 1.0585 kg at the default step); runoff stays minor at
   every step (at most 2.8–7.2 % of the film formed inside a latitude cap leaves it, and the equatorial ring sprays
   7.9–9.2 % of the mass); and the droplet population does not converge even at 0.00625 s (87.7, 106.6, 126.5 and 141.2
   million droplets over the continuum window; the thick branch's share of the sprayed mass 47.0, 35.3, 25.1 and
   20.5 %), because the step decides which branch takes the mass while each branch's droplets stay the same (fact 75).
   The 0.0125 s flight that crashed now runs to 120 s with its emptying rates bounded (fact 74); the runs, the probe and
   their cost are in fact 76. Repeat the series once fact 77 (a)'s branch-test decision is taken.
5. **What remains of the earlier time-step items.** The uniform-step runs that items 5 of 2026-10-02 and 2026-10-03
   asked for (the cascade at 0.25 and 0.125 s throughout, the whole 100 mm flight and the 50 mm flight with the cascade,
   at a smaller step) were not made; the switched-step series covers the continuum window of the 100 mm flight, where
   the step changes the answer, and keeps the default step before it, where the flight is in the merged branch and only
   the thin-film mode acts (the 87 g sprayed before 49.5 s, fact 69, whose own step dependence was not measured). Repeat
   the series on the prototype as Task 14 runs it (sub-plan 02's materials, Task 16's derived surface, fact 44's
   `PHI_DEATH`) once the time-step decision of fact 77 is taken.
6. **Quote every difference against the scatter at its own step.** Fact 65's scatter was measured at the default step;
   at fine steps it is not the same, so one more seed was run at 0.05 s and at 0.0125 s (seeds 12345 and 1, fact 75). At
   0.05 s the two runs differ by 0.008 % in sprayed mass, 0.72 % in droplet count, 0.37 % and 0.25 % in the median
   radius by number and by mass, and 17 % in the front-surface Rayleigh–Taylor release; at 0.0125 s by 1.45 % in sprayed
   mass and 3.5 % in the mass at 120 s — the shape feedback amplifies the round-off more at the finer step — but by only
   0.43 % in droplet count, 0.96 % and 0.47 % in the medians and 0.16 points in thick share. Run a seed pair at whatever
   step a result is quoted at; at a fine step the masses need it more than the droplet population does.

## Amendment of 2026-10-07 — the drivers on Scheil's curve with the rigid substrate, and the switched-step series repeated

> Asha's decisions of 2026-10-07 (facts 88–96 in `00-shared-context.md`; the rigid substrate is the amendment of
> 2026-10-06, facts 78–87). The code below is the tested code, as a diff against the copy of the rigid-substrate
> amendment. These changes stack on this sub-plan's amendments above, which still apply except where an item below
> replaces their wording.

1. **The drivers run the new default.** `melt_verification.py`'s `physics` mode and `melt_sensitivity.py`'s `base` run
   name `--material AA7075_scheil` (they named `AA7075_range`), so every physics row and every sensitivity row is now on
   Scheil's curve with the rigid substrate on. The bookkeeping device (`AA7075`, `--removal instant`, `--runoff off`) is
   untouched: with instant removal there is no film, so the rule never acts, and its rows and thresholds stand. The
   resolved mode (`AA7075`, Girin removal, D₀ and R₀) transports a film and, under Girin's closure, can now take the
   thick branch over slurry — on `AA7075` the slurry band is the 2 K above its 850 K half-liquid point, so the effect is
   expected to be small, which was not measured — and its two rows are to be re-measured with the melting references of
   Task 11.
2. **Two new sensitivity variants.** `norigid` (`--rigid-substrate off`: the regime test reads the fully molten depth,
   as before 2026-10-06) and `range` (`--material AA7075_range`: the linear law, `base`'s material before 2026-10-07).
   Between them they carry the two changes the table's earlier rows were measured without, so each can be read against
   `base` and `seed1` at the default step.
3. **What moves, measured rather than assumed** (facts 82–84 and 89–95). At the default step on the 100 mm flight to 120
   s the rule raises the sprayed mass by 1.2 % on the linear law and 0.6 % on Scheil's, lowers the droplet count by 32 %
   and 26 %, and raises the median radius by number by 4 % and 3 %; on the 50 mm flight, which never has Girin's
   closure, two-thirds of the sprayed mass leaves by the front-surface Rayleigh–Taylor mode and the median radius by
   mass is 7.4 times larger. Every physics and sensitivity row is therefore to be re-run; none was here.
4. **The switched-step series, repeated with the rule on** (facts 91–95), with the 2026-10-05 harness unchanged
   (`dtswitch.py` and `compare.py`, kept with the throwaway copy; fact 77 (b) still asks whether the switch should
   become a model option). It replaces item 4 of the 2026-10-05 amendment's measured values as the series for the
   default model; the 2026-10-05 series stays the record of the linear law with the rule off.
5. **Quote every difference against the scatter at its own step** (item 6 of 2026-10-05, unchanged): this series ran
   seed 1 at 0.05 and 0.0125 s again, and fact 92 gives the spreads.

In `analysis/melt_verification.py`:

```diff
--- a/analysis/melt_verification.py
+++ b/analysis/melt_verification.py
@@ -11,7 +11,8 @@
 time within 2 %); `resolved` -- the same heating and material with the real conductivity, film + runoff + Girin
 spraying on the default layered mesh, D0/R0 kept as SESAM keeps them (reported: the surface melts before the
 interior is hot, so the mass leaves earlier and, per unit heat, the interior's sensible heating delays the end);
-`physics` -- physics-mode heating, AA7075_range, Girin removal with the size feedback (the model proper; the SESAM
+`physics` -- physics-mode heating, AA7075_scheil (the melting default since 2026-10-07), Girin removal with the size
+feedback (the model proper; the SESAM
 overlay is context, not a target). Every run writes its
 overlay + residual plots and metrics JSON through the CLI; cases left out are read back from existing JSONs so
 summary.md / summary.json cover everything available."""
@@ -32,7 +33,7 @@
 MODES = {
     "bookkeeping": ["--heating", "sesam", "--material", "AA7075", "--removal", "instant", "--runoff", "off", "--k-scale", "1e4", "--prism-layers", "0"],
     "resolved": ["--heating", "sesam", "--material", "AA7075", "--removal", "girin", "--size-feedback", "initial"],
-    "physics": ["--heating", "physics", "--material", "AA7075_range", "--removal", "girin"],
+    "physics": ["--heating", "physics", "--material", "AA7075_scheil", "--removal", "girin"],
 }
 THRESHOLDS = {"mass_rel_m0": 0.02, "onset_km": 0.5, "demise_time_rel": 0.02}     # bookkeeping mode only
 COLUMNS = ["case", "mode", "max |dm| (of m0)", "melt onset [km] (model / SESAM)", "1 %-mass time [s] (model / SESAM)",
```

In `analysis/melt_sensitivity.py`:

```diff
--- a/analysis/melt_sensitivity.py
+++ b/analysis/melt_sensitivity.py
@@ -2,14 +2,16 @@
 """Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).
 
     "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
-        [--variants base,seed1,layers2,layers6,dt025,bridged,norunoff,nodeep,nocascade,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]
+        [--variants base,seed1,layers2,layers6,dt025,bridged,norunoff,nodeep,nocascade,norigid,we308,kr-30,kr+30,kt-30,kt+30,AA7075,range,sizeinitial,gammapm14,knbody003,fenicsx]
 
 Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
 six layers from 0.125 mm, 15.9 mm in all -- eight layers of 0.25 mm with growth 2 would exceed the radius; `gammapm14` and
 `knbody003` bound the two modelling choices of the 2026-09-22 amendment, the Prandtl-Meyer gamma and the body gate; `nodeep`
 turns off the deep runoff of the 2026-10-02 amendment, the liquid below the conjugate depth then staying in its elements;
 `nocascade` turns off the molten cascade of the 2026-10-03 amendment, the surface then receding through molten
-material by one element per macro step; `seed1` reruns `base` with another seed of numpy's generator -- every run is
+material by one element per macro step; `norigid` turns off the rigid substrate of the 2026-10-06 amendment, the
+regime test then reading the fully molten depth instead of the depth down to the half-liquid point; `range` runs the
+linear melting range instead of Scheil's curve, the melting default and `base`'s material since 2026-10-07; `seed1` reruns `base` with another seed of numpy's generator -- every run is
 seeded, `base` with the default 12345 (amendment of 2026-10-05) -- so its change relative to `base` is the run-to-run
 scatter, the floor against which every other row is read); the table
 lists sprayed mass, median droplet radius, melt-onset, spraying-onset and demise altitudes and the runtime, with the
@@ -30,9 +32,9 @@
 VARIANTS = {
     "base": [], "seed1": ["--seed", "1"], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
     "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "nodeep": ["--deep-runoff", "off"],
-    "nocascade": ["--molten-cascade", "off"], "we308": ["--we-critical", "3.08"],
+    "nocascade": ["--molten-cascade", "off"], "norigid": ["--rigid-substrate", "off"], "we308": ["--we-critical", "3.08"],
     "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
-    "AA7075": ["--material", "AA7075"], "sizeinitial": ["--size-feedback", "initial"], "gammapm14": ["--gamma-pm", "1.4"],
+    "AA7075": ["--material", "AA7075"], "range": ["--material", "AA7075_range"], "sizeinitial": ["--size-feedback", "initial"], "gammapm14": ["--gamma-pm", "1.4"],
     "knbody003": ["--kn-body-shock", "0.003"], "fenicsx": ["--thermal-solver", "fenicsx"],
 }
 KEYS = ["sprayed_mass_kg", "r_median_um", "melt_onset_altitude_km", "spraying_onset_altitude_km", "demise_altitude_km", "runtime_s"]
@@ -41,7 +43,7 @@
 
 def argv_for(key, variant, outdir):
     return ["run"] + CASES[key] + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
-                                   "--melt", "on", "--heating", "physics", "--material", "AA7075_range", "--no-particles",
+                                   "--melt", "on", "--heating", "physics", "--material", "AA7075_scheil", "--no-particles",
                                    "--outdir", outdir, "--name", "{}__{}".format(key, variant), "--quiet"] + VARIANTS[variant]
 
 
```

## Amendment of 2026-10-07 (continuum step) — the drivers switch by default, and the harness retires

> Asha's decision of 2026-10-07 on fact 96 (a) (facts 97–100 in `00-shared-context.md`). The code below is the tested
> code, as a diff against the copy of this date's earlier amendment. These items stack on this sub-plan's amendments
> above.

1. **Every physics and sensitivity run now switches.** `melt_verification.py`'s `physics` mode and every
   `melt_sensitivity.py` variant run with `--removal girin`, so from the first step under the continuum boundary they
   take 0.0125 s steps (sub-plan 13's default). The bookkeeping and resolved modes run `--removal instant` and Girin
   removal with `--size-feedback initial` respectively: the first keeps the default step throughout (the flag's default
   is off for instant removal), the second switches like the physics mode. The 50 mm flight never crosses the continuum
   boundary, so it is unchanged bit for bit (fact 99); the 100 mm flight switches at 49.5 s and to 120 s costs about
   twelve times as much (fact 99 gives the measured times).
2. **A `dtcoff` variant** (`--dt-continuum off`): the default 0.5 s step throughout, as every row before this amendment
   was run, so the table carries the step's effect against `base` like every other switch's. `dt025` now halves the step
   before the switch only.
3. **The switched-step harness retires.** `dtswitch.py` and `compare.py` (2026-10-05) wrapped the model from outside the
   package to switch the step; `--dt-continuum <s>` does the same from the model's own command line, and reproduces the
   harness's 0.0125 s and 0.5 s runs of the 2026-10-07 series bit for bit (fact 99). A step series is now a set of
   ordinary runs with `--dt-continuum 0.05`, `0.025`, `0.0125`, `0.00625` and `off`; the harness's runoff-by-latitude
   instrumentation stays in the throwaway copy, as it never belonged in `analysis/`.

In `analysis/melt_sensitivity.py`:

```diff
--- a/analysis/melt_sensitivity.py
+++ b/analysis/melt_sensitivity.py
@@ -2,18 +2,22 @@
 """Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).
 
     "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
-        [--variants base,seed1,layers2,layers6,dt025,bridged,norunoff,nodeep,nocascade,norigid,we308,kr-30,kr+30,kt-30,kt+30,AA7075,range,sizeinitial,gammapm14,knbody003,fenicsx]
+        [--variants base,seed1,layers2,layers6,dt025,dtcoff,bridged,norunoff,nodeep,nocascade,norigid,we308,kr-30,kr+30,
+                    kt-30,kt+30,AA7075,range,sizeinitial,gammapm14,knbody003,fenicsx]
 
 Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
 six layers from 0.125 mm, 15.9 mm in all -- eight layers of 0.25 mm with growth 2 would exceed the radius; `gammapm14` and
 `knbody003` bound the two modelling choices of the 2026-09-22 amendment, the Prandtl-Meyer gamma and the body gate; `nodeep`
 turns off the deep runoff of the 2026-10-02 amendment, the liquid below the conjugate depth then staying in its elements;
 `nocascade` turns off the molten cascade of the 2026-10-03 amendment, the surface then receding through molten
-material by one element per macro step; `norigid` turns off the rigid substrate of the 2026-10-06 amendment, the
-regime test then reading the fully molten depth instead of the depth down to the half-liquid point; `range` runs the
-linear melting range instead of Scheil's curve, the melting default and `base`'s material since 2026-10-07; `seed1` reruns `base` with another seed of numpy's generator -- every run is
-seeded, `base` with the default 12345 (amendment of 2026-10-05) -- so its change relative to `base` is the run-to-run
-scatter, the floor against which every other row is read); the table
+material by one element per macro step; `norigid` turns off the rigid substrate of the 2026-10-06 amendment, the regime
+test then reading the fully molten depth instead of the depth down to the half-liquid point; `range` runs the linear
+melting range instead of Scheil's curve, the melting default and `base`'s material since 2026-10-07; `dtcoff` keeps the
+default 0.5 s step through the continuum regime instead of the 0.0125 s continuum step every melting run with Girin
+removal now switches to there (amendment of 2026-10-07; `dt025` halves the step before the switch); `seed1` reruns
+`base` with another seed of numpy's generator -- every run is seeded, `base` with the default 12345 (amendment of
+2026-10-05) -- so its change relative to `base` is the run-to-run scatter, the floor against which every other row is
+read); the table
 lists sprayed mass, median droplet radius, melt-onset, spraying-onset and demise altitudes and the runtime, with the
 change relative to `base`. Every run is a subprocess of `--python` (default: this interpreter); the `fenicsx` variant
 needs the fenicsx_env interpreter (run it separately with --variants fenicsx --python <fenicsx_env python>, with CC
@@ -30,7 +34,7 @@
 
 CASES = {"d100": ["--diameter", "100", "--altitude", "77.500133"], "d050": ["--diameter", "50", "--altitude", "115"]}
 VARIANTS = {
-    "base": [], "seed1": ["--seed", "1"], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
+    "base": [], "seed1": ["--seed", "1"], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"], "dtcoff": ["--dt-continuum", "off"],
     "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "nodeep": ["--deep-runoff", "off"],
     "nocascade": ["--molten-cascade", "off"], "norigid": ["--rigid-substrate", "off"], "we308": ["--we-critical", "3.08"],
     "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
```

---

**Depends on:** Task 13 (CLI) and Task 12 (metrics).
**Produces:** the reference-tier test, the verification/sensitivity driver scripts, and the runs that produce the plan's headline numbers — including the finding that a 100 mm sphere, with size feedback and the amended surface flow, stops demising and has 29.9% of its mass reach the ground.
**Character:** testing/verification/analysis — also the task with the largest real-world runtime cost in the plan.
**Read before implementing:** Measured fact 23 in the shared context: one macro time-step costs about 2.0 seconds of wall-clock time, and a full flight whose remnant survives to the ground can take up to roughly 2,357 seconds (about 40 minutes) — far longer than the spec's original under-6-minute target, which the master plan flags as needing to be restated rather than treated as already fixed. Budget real wall-clock time for whichever agent or CI system executes this task. Measured fact 24 records the 29.9%-mass-surviving finding itself, and is explicit that it is an artifact of the model's fixed body attitude (no tumbling simulated yet) — carry that caveat into every result this task reports; do not state it as a general re-entry finding.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan, and resolve the runtime-target restatement (see the shared context's note on this) before finalizing the plan's acceptance criteria.

---

### Task 14: Verification and sensitivity drivers, the reference-tier test, the runs

**Files:**
- Create: `analysis/melt_verification.py`, `analysis/melt_sensitivity.py`, `tests/test_reentry_model_reference_melt.py`
- Outputs (git-ignored): `reentry_model_output/verification_melt/{summary.md, summary.json, <case>__<mode>.*}`, `.../sensitivity/{sensitivity.md, sensitivity.json, ...}`, `.../girin/` (Task 8), `.../reference_tests/`

**Interfaces:**
- Consumes: the CLI (Task 13), `compare.melt_metrics` (Task 12), the melting references (Task 11), `sesam_io.REFERENCE_DIR`.
- Produces: the verification table (bookkeeping thresholded, resolved and physics reported), the sensitivity table, the reference-tier test with the spec §13.1 thresholds (mass 2 % of m₀, onset 0.5 km, 1 %-mass time 2 %).

- [ ] **Step 1: Create `tests/test_reentry_model_reference_melt.py`**

```python
"""Melting model vs the two melting US76 SESAM references (marker: reference, ~2 min): the bookkeeping device
(SESAM-equivalent heating, AA7075, instant removal, k x 1e4, D0/R0 kept, the Step 2 default mesh) against the acceptance
thresholds of spec section 13.1 -- mass within 2 % of the initial at every reference time, melt-onset altitude within
0.5 km, the 1 %-mass time within 2 % (measured 2026-09-21: 0.98 % / +0.10 km / -1.3 % for 100 mm, 1.29 % / +0.17 km /
-0.2 % for 50 mm). The resolved and physics-mode runs are analysis/melt_verification.py's business (reported). Metrics
are written to reentry_model_output/verification_melt/reference_tests/."""
import json
import os

import pytest

from reentry_model import aero, atmosphere, body, compare, coupled, heating, material, mesh, sesam_io, thermal
from reentry_model import trajectory as tj

NAMES = {"d100": "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind",
         "d050": "sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind"}
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reentry_model_output", "verification_melt", "reference_tests")
THRESHOLDS = {"mass_rel_m0": 0.02, "onset_km": 0.5, "demise_time_rel": 0.02}


def bookkeeping_run(ref):
    initial = tj.InitialState(ref.initial.velocity, ref.initial.altitude, ref.initial.flight_path, ref.initial.heading,
                              ref.initial.lat, ref.initial.lon, ref.initial.epoch)
    the_mesh = mesh.sphere_mesh(ref.diameter / 2.0)
    mat = material.Material.from_drama_json("AA7075")
    mat.k_table = mat.k_table * 1e4
    the_body = body.MeltingBody(the_mesh, mat, thermal.thermal_solver("skfem"), body.sphere_mass(ref.diameter, ref.material_density),
                                settings=body.MeltSettings(removal="instant", runoff=False, size_feedback="initial"))
    sim = tj.Simulator(initial, the_body, atmosphere.US76TableAtmosphere(), aero.SphereDragTables.from_json(), aero.SesamTable(),
                       tj.Settings(diameter=ref.diameter))
    return coupled.CoupledRun(sim, the_body, heating.SesamEquivalentHeating(), coupled.CoupledSettings(dt=0.5)).run()


@pytest.mark.reference
@pytest.mark.parametrize("key", ["d100", "d050"])
def test_bookkeeping_device_follows_sesam(key):
    ref = sesam_io.load_reference(os.path.join(sesam_io.REFERENCE_DIR, NAMES[key] + ".csv"))
    hist = bookkeeping_run(ref)
    mm = compare.melt_metrics(hist, ref)
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, key + "__bookkeeping.json"), "w") as fh:
        json.dump({"case": ref.name, "results": hist.results, "melt_metrics": mm}, fh, indent=2, default=str)
    assert hist.end_reason == "demise" and abs(hist.results["melt_energy_balance_residual"]) < 1e-6
    assert mm["mass"]["max_rel_m0"] <= THRESHOLDS["mass_rel_m0"]
    assert abs(mm["onset_altitude_diff_km"]) <= THRESHOLDS["onset_km"]
    assert abs(mm["demise_time_rel"]) <= THRESHOLDS["demise_time_rel"]
```


- [ ] **Step 2: Create `analysis/melt_verification.py`**

```python
#!/usr/bin/env python3
"""Run the melting model against the two melting US76 SESAM references and tabulate the mass-loss errors (Step 3
verification, spec section 13.1).

    "$PY" analysis/melt_verification.py [--outdir reentry_model_output/verification_melt] [--cases d100,d050]
        [--modes bookkeeping,resolved,physics] [--animate]

Modes: `bookkeeping` -- SESAM-equivalent heating, AA7075 (DRAMA's single melting temperature), --removal instant,
--runoff off, --k-scale 1e4 (near-isothermal body): the thresholded check of the melting bookkeeping against SESAM's
lumped Q/L_f law (mass within 2 % of the initial at every reference time, onset altitude within 0.5 km, 1 %-mass
time within 2 %); `resolved` -- the same heating and material with the real conductivity, film + runoff + Girin
spraying on the default layered mesh, D0/R0 kept as SESAM keeps them (reported: the surface melts before the
interior is hot, so the mass leaves earlier and, per unit heat, the interior's sensible heating delays the end);
`physics` -- physics-mode heating, AA7075_range, Girin removal with the size feedback (the model proper; the SESAM
overlay is context, not a target). Every run writes its
overlay + residual plots and metrics JSON through the CLI; cases left out are read back from existing JSONs so
summary.md / summary.json cover everything available."""
import argparse
import json
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from reentry_model import cli, sesam_io  # noqa: E402

CASES = {
    "d100": ("sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind", ["--diameter", "100", "--altitude", "77.500133"]),
    "d050": ("sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind", ["--diameter", "50", "--altitude", "115"]),
}
MODES = {
    "bookkeeping": ["--heating", "sesam", "--material", "AA7075", "--removal", "instant", "--runoff", "off", "--k-scale", "1e4", "--prism-layers", "0"],
    "resolved": ["--heating", "sesam", "--material", "AA7075", "--removal", "girin", "--size-feedback", "initial"],
    "physics": ["--heating", "physics", "--material", "AA7075_range", "--removal", "girin"],
}
THRESHOLDS = {"mass_rel_m0": 0.02, "onset_km": 0.5, "demise_time_rel": 0.02}     # bookkeeping mode only
COLUMNS = ["case", "mode", "max |dm| (of m0)", "melt onset [km] (model / SESAM)", "1 %-mass time [s] (model / SESAM)",
           "sprayed / film left [kg]", "droplets (median r)", "runtime", "verdict"]


def run_case(key, mode, outdir, animate):
    name, args = CASES[key]
    ref = os.path.join(sesam_io.REFERENCE_DIR, name + ".csv")
    label = "{}__{}".format(key, mode)
    argv = ["run"] + args + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                             "--melt", "on"] + MODES[mode] + ["--reference", ref, "--outdir", outdir, "--name", label, "--quiet"]
    if animate:
        argv.append("--animate")
    rc = cli.main(argv)
    if rc != 0:
        raise SystemExit("run {} failed with exit code {}".format(label, rc))
    return label


def row_from_doc(key, mode, doc):
    mm, res = doc["comparison"]["melt_metrics"], doc["results"]
    fmt = lambda x, f="{:.1f}": "n/a" if x is None else f.format(x)
    verdict = ""
    if mode == "bookkeeping":
        ok = (mm["mass"]["max_rel_m0"] <= THRESHOLDS["mass_rel_m0"] and mm["onset_altitude_diff_km"] is not None
              and abs(mm["onset_altitude_diff_km"]) <= THRESHOLDS["onset_km"] and mm["demise_time_rel"] is not None
              and abs(mm["demise_time_rel"]) <= THRESHOLDS["demise_time_rel"])
        verdict = "pass" if ok else "FAIL"
    return {"case": key, "mode": mode, "max |dm| (of m0)": "{:.2%}".format(mm["mass"]["max_rel_m0"]),
            "melt onset [km] (model / SESAM)": "{} / {}".format(fmt(mm["onset_altitude_model_km"], "{:.2f}"), fmt(mm["onset_altitude_reference_km"], "{:.2f}")),
            "1 %-mass time [s] (model / SESAM)": "{} / {}{}".format(fmt(mm["demise_time_model_s"]), fmt(mm["demise_time_reference_s"]),
                                                              "" if mm["demise_time_rel"] is None else " ({:+.1%})".format(mm["demise_time_rel"])),
            "sprayed / film left [kg]": "{:.4f} / {:.4f}".format(res["sprayed_mass_kg"], res["film_mass_kg"]),
            "droplets (median r)": "{:.3g} ({} um)".format(res["n_released"], fmt(res["r_median_um"])),
            "runtime": "{:.0f} s, {} steps, {:.1f} it/step".format(res["runtime_s"], res["n_macro_steps"], res["mean_newton_iterations"]),
            "verdict": verdict}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--modes", default=",".join(MODES))
    p.add_argument("--animate", action="store_true")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    for key in args.cases.split(","):
        for mode in args.modes.split(","):
            run_case(key, mode, args.outdir, args.animate)
    rows = []
    for key in CASES:
        for mode in MODES:
            path = os.path.join(args.outdir, "{}__{}.json".format(key, mode))
            if os.path.isfile(path):
                doc = json.load(open(path))
                if doc.get("comparison") and "melt_metrics" in doc["comparison"]:
                    rows.append(row_from_doc(key, mode, doc))
    lines = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
    lines += ["| " + " | ".join(str(r[c]) for c in COLUMNS) + " |" for r in rows]
    md = "\n".join(lines) + "\n"
    open(os.path.join(args.outdir, "summary.md"), "w").write(md)
    json.dump(rows, open(os.path.join(args.outdir, "summary.json"), "w"), indent=2)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```


- [ ] **Step 3: Create `analysis/melt_sensitivity.py`**

```python
#!/usr/bin/env python3
"""Sensitivity and convergence table of the melting model (spec Step 3 section 13.5).

    "$PY" analysis/melt_sensitivity.py [--outdir reentry_model_output/verification_melt/sensitivity] [--cases d100,d050]
        [--variants base,layers2,layers6,dt025,bridged,norunoff,we308,kr-30,kr+30,kt-30,kt+30,AA7075,sizeinitial,gammapm14,knbody003,fenicsx]

Each variant is one physics-mode melting flight (US76, winds off) differing from `base` in one setting (`layers6`:
six layers from 0.125 mm, 15.9 mm in all -- eight layers of 0.25 mm with growth 2 would exceed the radius; `gammapm14` and
`knbody003` bound the two modelling choices of the 2026-09-22 amendment, the Prandtl-Meyer gamma and the body gate); the table
lists sprayed mass, median droplet radius, melt-onset, spraying-onset and demise altitudes and the runtime, with the
change relative to `base`. Every run is a subprocess of `--python` (default: this interpreter); the `fenicsx` variant
needs the fenicsx_env interpreter (run it separately with --variants fenicsx --python <fenicsx_env python>, with CC
exported). Variants left out are read back from existing JSONs."""
import argparse
import json
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)


CASES = {"d100": ["--diameter", "100", "--altitude", "77.500133"], "d050": ["--diameter", "50", "--altitude", "115"]}
VARIANTS = {
    "base": [], "layers2": ["--prism-layers", "2"], "layers6": ["--prism-layers", "6", "--layer-thickness", "0.125"], "dt025": ["--dt", "0.25"],
    "bridged": ["--rarefied-shear", "bridged"], "norunoff": ["--runoff", "off"], "we308": ["--we-critical", "3.08"],
    "kr-30": ["--kr", "0.119"], "kr+30": ["--kr", "0.221"], "kt-30": ["--kt", "0.77"], "kt+30": ["--kt", "1.43"],
    "AA7075": ["--material", "AA7075"], "sizeinitial": ["--size-feedback", "initial"], "gammapm14": ["--gamma-pm", "1.4"],
    "knbody003": ["--kn-body-shock", "0.003"], "fenicsx": ["--thermal-solver", "fenicsx"],
}
KEYS = ["sprayed_mass_kg", "r_median_um", "melt_onset_altitude_km", "spraying_onset_altitude_km", "demise_altitude_km", "runtime_s"]
COLUMNS = ["case", "variant", "sprayed [kg]", "median r [um]", "melt onset [km]", "spraying onset [km]", "demise [km]", "runtime [s]"]


def argv_for(key, variant, outdir):
    return ["run"] + CASES[key] + ["--velocity", "7.5", "--flight-path-angle", "-0.959331", "--atmosphere", "us76", "--thermal", "fem",
                                   "--melt", "on", "--heating", "physics", "--material", "AA7075_range", "--no-particles",
                                   "--outdir", outdir, "--name", "{}__{}".format(key, variant), "--quiet"] + VARIANTS[variant]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(REPO_ROOT, "reentry_model_output", "verification_melt", "sensitivity"))
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--variants", default=",".join(v for v in VARIANTS if v != "fenicsx"))
    p.add_argument("--python", default=sys.executable, help="interpreter for the runs (default: this one; the fenicsx variant needs fenicsx_env's)")
    args = p.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    for key in args.cases.split(","):
        for variant in args.variants.split(","):
            # every run in a fresh interpreter: one flight is ~250 s and a 46 k-node solver; a long in-process series was
            # seen to crawl (measured 2026-09-21)
            env = dict(os.environ, FI_PROVIDER="tcp")
            subprocess.run([args.python, "-m", "reentry_model"] + argv_for(key, variant, args.outdir), cwd=REPO_ROOT, check=True, env=env)
    rows = []
    for key in CASES:
        base = None
        for variant in VARIANTS:
            path = os.path.join(args.outdir, "{}__{}.json".format(key, variant))
            if not os.path.isfile(path):
                continue
            res = json.load(open(path))["results"]
            if variant == "base":
                base = res
            fmt = lambda k, f: "n/a" if res.get(k) is None else (f.format(res[k]) + ("" if base is None or variant == "base" or base.get(k) is None or not base[k] else " ({:+.1%})".format(res[k] / base[k] - 1.0)))
            rows.append({"case": key, "variant": variant, "sprayed [kg]": fmt("sprayed_mass_kg", "{:.4f}"), "median r [um]": fmt("r_median_um", "{:.1f}"),
                         "melt onset [km]": fmt("melt_onset_altitude_km", "{:.2f}"), "spraying onset [km]": fmt("spraying_onset_altitude_km", "{:.2f}"),
                         "demise [km]": fmt("demise_altitude_km", "{:.2f}"), "runtime [s]": fmt("runtime_s", "{:.0f}")})
    lines = ["| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)] + ["| " + " | ".join(r[c] for c in COLUMNS) + " |" for r in rows]
    md = "\n".join(lines) + "\n"
    open(os.path.join(args.outdir, "sensitivity.md"), "w").write(md)
    json.dump(rows, open(os.path.join(args.outdir, "sensitivity.json"), "w"), indent=2)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```


- [ ] **Step 4: Run the reference-tier test and the verification driver**

```bash
"$PY" -m pytest -m reference tests/test_reentry_model_reference_melt.py -q          # ~1 min: both bookkeeping runs
"$PY" analysis/melt_verification.py --animate                                         # ~6 min: 2 spheres x 3 modes, with the videos
```

Expected: 2 passed (measured 0.98 % / +0.10 km / −1.3 % and 1.29 % / +0.17 km / −0.2 %); the driver prints the six-row table with `pass` in both bookkeeping rows and the resolved and physics rows of the README's verification table (Task 15), and writes `summary.md`. The whole reference tier (`"$PY" -m pytest -m reference -q`) is 15 passed in ≈ 24 min. Check that `reentry_model_output/verification_melt/d100__physics/vtk/` holds `animation.mp4`, `film.mp4`, `section.mp4` and the stills (`melt_onset`, `spraying_onset`, `peak_release`, `film_*`, `section_*`), and that `d100__resolved/mass_time.png` shows the SESAM overlay with its residual panel.

- [ ] **Step 5: Run the sensitivity study**

```bash
"$PY" analysis/melt_sensitivity.py                                                    # ~60 min: 13 variants x 2 spheres
"$PY" analysis/melt_sensitivity.py --variants fenicsx --python "$FX"                  # the backend variant from fenicsx_env (CC and FI_PROVIDER exported)
```

(export `CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang` before the second command; the script sets `FI_PROVIDER=tcp` itself.) Expected: `sensitivity.md` with every variant's sprayed mass, median radius, onsets and demise altitude and the change relative to `base`; the `fenicsx` row equal to `base` to the printed digits; `layers2`/`layers6`/`dt025` within a few percent of `base` in demise altitude and sprayed mass; `we308`, `kt±30` nearly identical (spraying is melt-limited), `kr-30` moving the median radius; `AA7075` earlier onset (850 K vs 908 K liquidus) and a leeward remnant that reaches the ground (demise n/a); `sizeinitial` ending higher (the nose stays at R₀, the body keeps the sphere's drag). `gammapm14` (`--gamma-pm 1.4`) and `knbody003` (`--kn-body-shock 0.003`) bound the two new modelling choices: the first halves the wall pressure beyond the sonic point, the second shrinks the Girin-certified band. Record the table in Task 15.

- [ ] **Step 6: Run the FEniCSx tests once more and the whole unit tier**

```bash
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q
"$PY" -m pytest -m "not drama and not reference" -q
```

Expected: 8 passed; the unit tier green (~6 min).

- [ ] **Step 7: Commit**

```bash
git add analysis/melt_verification.py analysis/melt_sensitivity.py tests/test_reentry_model_reference_melt.py
git commit -m "Add the melt verification, sensitivity drivers and the reference-tier bookkeeping test (Step 3 Task 14)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

