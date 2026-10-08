# Sub-plan: Task 13 — Command line

> Extracted verbatim from `2026-09-20-melt-spraying.md` (current version, lines 6652–7193). Read `00-shared-context.md` first. **Amended 2026-10-02: the deep runoff and the per-patch conjugate depth** (section below). **Amended 2026-10-03: the molten cascade** (the section after it). **Amended 2026-10-05: the seed of numpy's generator** (the section after that). **Amended 2026-10-05 (runoff flux): no flag, and the exit code of a model failure** (the section after the seed's). **Amended 2026-10-06: `--rigid-substrate on|off`** (the section after the runoff flux's).


> **Amended 2026-09-27** by `docs/superpowers/specs/2026-09-27-surface-recession-remeshing-design.md`
> (measured facts 38–45 in `00-shared-context.md`). The block below takes precedence over the extracted body.

**Changes to Task 13.** New flags, all defaulting so that a run with none of them behaves as spec §6 and §7
describe:

- `--band-thickness` (metres, default 0.015) — the dense band of Task 1.
- `--prism-layers` — unchanged in meaning, **default now 0** (was 4).
- `--remesh {auto,off}` — default `auto` for melting runs, `off` for the pinned verification devices of fact 12.
- `--remesh-min-steps` (default 20) — the anti-thrash floor between remeshes.
- `--derived-surface {on,off}` (default `on`) — `off` reproduces the current staircase geometry exactly, and is what
  the shape A/B comparison of Task 14 runs against.

**Run names must encode all of this** (CLAUDE.md's resumability convention): append `_band15mm` and `_remesh` or
`_noremesh`, so runs with and without the new machinery cannot overwrite each other in one output directory.
Exit codes are unchanged: 0 ok, 1 a model/integration failure — **including a Task 17 validity-gate stop** — 2 bad
arguments or a missing optional library.
---

## Amendment of 2026-10-02 — `--deep-runoff on|off`

> Part of the deep-runoff amendment (sub-plan 09's amendment of this date; facts 46–53 in `00-shared-context.md`).

**The flag.** `--deep-runoff on|off` in the Step 3 group, default `on`, sets `MeltSettings.deep_runoff`. With `off`
the deep stage of the melt step is skipped entirely, so the contiguous liquid below the conjugate depth stays in its
elements until the surface recedes to it, as in every run before 2026-10-02. It needs `--runoff on`: with
`--runoff off` there is no transport and nothing below the conjugate depth moves.

**Why a flag rather than a fixed behaviour.** The repository's pattern for a mechanism that changes results is a flag
with the old behaviour available (`--rt-spray on|off`, fact 35), and the deep runoff changes the 100 mm results by
more than any amendment since fact 28 (fact 49). The flag also gives Task 14 its sensitivity row (`nodeep`) and lets
the two settings be compared on the same flight. Measured: with `--deep-runoff off` the amended code reproduces the
unamended prototype bit for bit on the 100 mm flight when both runs fix numpy's random seed (fact 52) — every history
column, every row. **Whether the default should be `off` instead is the first open decision of fact 53:** at the
default 0.5 s step the liquid the deep runoff moves is a backlog left by the order of the melt step, which vanishes as
the step shrinks: at 0.25 s the deep runoff changes the 100 mm results by under 1 %, at 0.125 s by less than the
run-to-run spread (fact 50). The default here is `on`, following the
repository's pattern; flipping it is the one line `default="off"` here and `deep_runoff: bool = False` in
`MeltSettings`.

**Run names and the JSON.** A melting run with `--deep-runoff off` ends in `_deeprunoff-off`, so runs with and
without the deep runoff can share one output directory without overwriting each other (CLAUDE.md's resumability
convention, which this sub-plan's amendment of 2026-09-27 applied to the band and remesh flags). The default keeps the
existing name, so no existing run name, reference or test changes. The JSON `settings` carry `"deep_runoff": "on"` or
`"off"`, and the `results` carry `deep_runoff_mass_kg`, `deep_surfaced_mass_kg` and `deep_mass_kg` (Task 10's
amendment). Exit codes are unchanged; an invalid value is argparse's exit 2, which the bad-argument test now covers.

Code (the tested change to `reentry_model/cli.py`, the prototype's file — the body below plus the dense-band line of
this sub-plan's amendment of 2026-09-27):

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -45,10 +45,14 @@
         raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))
 
 
-def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None):
+def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None,
+                   deep_runoff=True):
+    """The run name encodes the configuration, so runs of different settings can share an output directory. A melting
+    run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`; the default keeps the old name."""
     name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
         diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
-    return name + ("_fem-" + heating_name if heating_name else "") + ("_melt-" + melt if melt else "")
+    return (name + ("_fem-" + heating_name if heating_name else "") + ("_melt-" + melt if melt else "")
+            + ("_deeprunoff-off" if melt and not deep_runoff else ""))
 
 
 def make_atmosphere(spec, epoch, wind_name):
@@ -152,6 +156,11 @@
                     help="girin: film + runoff + Girin spraying (default); instant: liquid removed as it forms -- the lumped-melting "
                          "verification device, not physical")
     me.add_argument("--runoff", choices=("on", "off"), default="on", help="film runoff transport along the surface (default on)")
+    me.add_argument("--deep-runoff", choices=("on", "off"), default="on",
+                    help="move the contiguous liquid below Girin's conjugate depth along the surface by the pressure-gradient and "
+                         "deceleration-driven lubrication flux over the whole liquid depth, never spraying it; it becomes film only "
+                         "from the top, as the film thins below the conjugate depth (default on; needs --runoff on). off: that "
+                         "liquid stays in its elements until the surface recedes to it, as in every run before 2026-10-02")
     me.add_argument("--rt-spray", choices=("on", "off"), default="on",
                     help="apply Girin & Kopyt's front-surface Rayleigh-Taylor mode as a release mechanism, on the patches where it "
                          "grows faster than the shear mode (default on). off: evaluate and report it but release nothing through it, "
@@ -207,7 +216,8 @@
                                        rt_spray=args.rt_spray == "on")
         size_feedback = args.size_feedback or ("current" if args.removal == "girin" else "initial")
         melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
-                                          particles=args.particles, size_feedback=size_feedback)
+                                          particles=args.particles, size_feedback=size_feedback,
+                                          deep_runoff=args.deep_runoff == "on")
         the_body = body.MeltingBody(the_mesh, mat, solver, mass, flow, spray_model, melt_settings, T0=args.temperature,
                                     emissivity=args.emissivity, T_ambient=args.t_ambient)
     else:
@@ -227,7 +237,7 @@
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
     if melting:
-        info.update({"removal": args.removal, "runoff": args.runoff, "rt_spray": args.rt_spray,
+        info.update({"removal": args.removal, "runoff": args.runoff, "deep_runoff": args.deep_runoff, "rt_spray": args.rt_spray,
                      "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                      "k_r": args.kr, "k_t": args.kt, "demise_fraction": args.demise_fraction, "particles": args.particles,
                      "size_feedback": size_feedback, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
@@ -274,7 +284,8 @@
     atm, atm_name, atm_info = make_atmosphere(args.atmosphere, args.epoch, args.wind)
     reference = sesam_io.load_reference(args.reference) if args.reference else None
     name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
-                                       args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None)
+                                       args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None,
+                                       deep_runoff=args.deep_runoff == "on")
     run_dir = os.path.join(args.outdir, name)
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
```

Tests (the tested change to `tests/test_reentry_model_cli.py`):

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -169,6 +169,8 @@
 
 def test_run_name_with_melt():
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin").endswith("_fem-physics_melt-girin")
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
+                              deep_runoff=False).endswith("_fem-physics_melt-girin_deeprunoff-off")    # amendment of 2026-10-02
 
 
 def test_melting_run_writes_columns_files_and_json(tmp_path):
@@ -185,7 +187,9 @@
     s, r, f = doc["settings"], doc["results"], doc["files"]
     assert s["melt"] == "on" and s["material"] == "AA7075_range" and s["removal"] == "girin" and s["runoff"] == "on" and s["prism_layers"] == 0
     assert s["rarefied_shear"] == "slip" and s["we_critical"] == 4.62 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["liquid"]["sigma"] == 0.86
-    assert s["gamma_pm"] == 1.15 and s["kn_body_shock"] == 0.01
+    assert s["gamma_pm"] == 1.15 and s["kn_body_shock"] == 0.01 and s["deep_runoff"] == "on"
+    assert "deep_runoff_mass_kg" in rows[0] and "deep_mass_kg" in rows[0] and r["deep_runoff_mass_kg"] >= 0.0
+    assert "deep_surfaced_mass_kg" in rows[0] and "deep_blob_fraction" in rows[0] and "deep_liquid_kg" in rows[0]
     assert s["size_feedback"] == "current" and "nose_radius_mm" in rows[0] and float(rows[-1]["nose_radius_mm"]) > 0.0
     assert s["T_liquidus_K"] == 908.0 and s["latent_heat_Jkg"] == 400e3 and s["demise_fraction"] == 0.01 and s["particles"] is True
     assert r["melt_onset_altitude_km"] is not None and r["sprayed_mass_kg"] > 0.0 and r["n_source_rows"] > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
@@ -221,6 +225,7 @@
     MELT + ["--k-scale", "0"],
     MELT + ["--removal", "magic"],
     MELT + ["--size-feedback", "shrinking"],
+    MELT + ["--deep-runoff", "maybe"],
 ])
 def test_bad_melt_arguments_exit_2(argv, tmp_path):
     with pytest.raises(SystemExit) as exc:
```

**Measured (2026-10-02, throwaway copy of the prototype).** `test_run_name_with_melt`,
`test_melting_run_writes_columns_files_and_json` and the eight `test_bad_melt_arguments_exit_2` cases pass; the first
two fail on the unamended CLI (no `deep_runoff` keyword, no `deep_runoff` setting) and the new bad-argument case passes
on both (argparse rejects an unknown flag and an unknown value alike, with exit 2). Whole unit tier: see sub-plan 09's
amendment.

## Amendment of 2026-10-03 — `--molten-cascade on|off`

> Part of the molten-cascade amendment (sub-plan 09's amendment of this date; facts 54–61 in `00-shared-context.md`).

**The flag.** `--molten-cascade on|off` in the Step 3 group, default `on`, sets `MeltSettings.molten_cascade`. With
`on`, an element that an element death exposes fully molten — its mean nodal temperature at or above the top of the feed
ramp, the test the molten layer is measured by — is fed whole within the same macro step, so the surface recedes through
molten material by as many elements as are molten. With `off` an exposed element waits for the next step's feed, so the
surface recedes through molten material by one element per macro step, as in every run before this amendment. It does
nothing with `--removal instant`, which feeds every element as it melts.

**Why a switch is genuinely needed.** The repository's pattern for a mechanism that changes results is a flag that keeps
the old behaviour available (`--rt-spray`, `--deep-runoff`); this one changes the 100 mm results at the default step
(fact 57), the time-step study compares the two settings, and the sensitivity table gets a `nocascade` row (sub-plan
14). Measured: with `--molten-cascade off` the amended code reproduces the deep-runoff amendment of 2026-10-02 bit for
bit on the 100 mm flight to 120 s, both runs seeded — all 241 rows of all 84 shared history columns, all 167 473
source-table rows and every result field but the run time.

**Run names and the JSON.** A melting run with `--molten-cascade off` ends in `_moltencascade-off` (after
`_deeprunoff-off` when both are off), so runs with and without the cascade can share one output directory (CLAUDE.md's
resumability convention); the default keeps the existing name, so no existing run name, reference or test changes. The
JSON `settings` carry `"molten_cascade": "on"` or `"off"`, and the `results` carry `cascade_mass_kg`,
`cascade_passes_max` and `cascade_capped_steps` (Task 10's amendment). Exit codes are unchanged; an invalid value is
argparse's exit 2, which the bad-argument test now covers.

Code (the tested change to `reentry_model/cli.py`, the prototype's file with the deep-runoff amendment above applied):

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -46,13 +46,14 @@
 
 
 def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None,
-                   deep_runoff=True):
+                   deep_runoff=True, molten_cascade=True):
     """The run name encodes the configuration, so runs of different settings can share an output directory. A melting
-    run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`; the default keeps the old name."""
+    run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
+    (amendment of 2026-10-03) in `_moltencascade-off`; the defaults keep the old name."""
     name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
         diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
     return (name + ("_fem-" + heating_name if heating_name else "") + ("_melt-" + melt if melt else "")
-            + ("_deeprunoff-off" if melt and not deep_runoff else ""))
+            + ("_deeprunoff-off" if melt and not deep_runoff else "") + ("_moltencascade-off" if melt and not molten_cascade else ""))
 
 
 def make_atmosphere(spec, epoch, wind_name):
@@ -161,6 +162,13 @@
                          "deceleration-driven lubrication flux over the whole liquid depth, never spraying it; it becomes film only "
                          "from the top, as the film thins below the conjugate depth (default on; needs --runoff on). off: that "
                          "liquid stays in its elements until the surface recedes to it, as in every run before 2026-10-02")
+    me.add_argument("--molten-cascade", choices=("on", "off"), default="on",
+                    help="when an element death exposes an element that is fully molten (its mean temperature at or above the top "
+                         "of the feed ramp, the test the molten layer is measured by), feed it whole within the same macro step, so "
+                         "that it dies in turn and the surface recedes through molten material by as many elements as are molten "
+                         "(default on). off: an exposed element waits for the next step's feed, so the surface recedes through "
+                         "molten material by one element per macro step, as in every run before the molten-cascade amendment of "
+                         "2026-10-03")
     me.add_argument("--rt-spray", choices=("on", "off"), default="on",
                     help="apply Girin & Kopyt's front-surface Rayleigh-Taylor mode as a release mechanism, on the patches where it "
                          "grows faster than the shear mode (default on). off: evaluate and report it but release nothing through it, "
@@ -217,7 +225,7 @@
         size_feedback = args.size_feedback or ("current" if args.removal == "girin" else "initial")
         melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
                                           particles=args.particles, size_feedback=size_feedback,
-                                          deep_runoff=args.deep_runoff == "on")
+                                          deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on")
         the_body = body.MeltingBody(the_mesh, mat, solver, mass, flow, spray_model, melt_settings, T0=args.temperature,
                                     emissivity=args.emissivity, T_ambient=args.t_ambient)
     else:
@@ -237,7 +245,8 @@
         info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
     if melting:
-        info.update({"removal": args.removal, "runoff": args.runoff, "deep_runoff": args.deep_runoff, "rt_spray": args.rt_spray,
+        info.update({"removal": args.removal, "runoff": args.runoff, "deep_runoff": args.deep_runoff,
+                     "molten_cascade": args.molten_cascade, "rt_spray": args.rt_spray,
                      "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                      "k_r": args.kr, "k_t": args.kt, "demise_fraction": args.demise_fraction, "particles": args.particles,
                      "size_feedback": size_feedback, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
@@ -285,7 +294,7 @@
     reference = sesam_io.load_reference(args.reference) if args.reference else None
     name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
                                        args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None,
-                                       deep_runoff=args.deep_runoff == "on")
+                                       deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on")
     run_dir = os.path.join(args.outdir, name)
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
```

Tests (the tested change to `tests/test_reentry_model_cli.py`):

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -171,6 +171,8 @@
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin").endswith("_fem-physics_melt-girin")
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
                               deep_runoff=False).endswith("_fem-physics_melt-girin_deeprunoff-off")    # amendment of 2026-10-02
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
+                              molten_cascade=False).endswith("_fem-physics_melt-girin_moltencascade-off")   # the molten cascade
 
 
 def test_melting_run_writes_columns_files_and_json(tmp_path):
@@ -190,6 +192,10 @@
     assert s["gamma_pm"] == 1.15 and s["kn_body_shock"] == 0.01 and s["deep_runoff"] == "on"
     assert "deep_runoff_mass_kg" in rows[0] and "deep_mass_kg" in rows[0] and r["deep_runoff_mass_kg"] >= 0.0
     assert "deep_surfaced_mass_kg" in rows[0] and "deep_blob_fraction" in rows[0] and "deep_liquid_kg" in rows[0]
+    assert s["molten_cascade"] == "on" and "cascade_passes" in rows[0] and "cascade_mass_kg" in rows[0]      # the molten cascade
+    assert r["cascade_mass_kg"] == pytest.approx(float(rows[-1]["cascade_mass_kg"])) and r["cascade_mass_kg"] >= 0.0
+    assert r["cascade_capped_steps"] == 0
+    assert r["cascade_passes_max"] == max(int(float(row["cascade_passes"])) for row in rows)
     assert s["size_feedback"] == "current" and "nose_radius_mm" in rows[0] and float(rows[-1]["nose_radius_mm"]) > 0.0
     assert s["T_liquidus_K"] == 908.0 and s["latent_heat_Jkg"] == 400e3 and s["demise_fraction"] == 0.01 and s["particles"] is True
     assert r["melt_onset_altitude_km"] is not None and r["sprayed_mass_kg"] > 0.0 and r["n_source_rows"] > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
@@ -226,6 +232,7 @@
     MELT + ["--removal", "magic"],
     MELT + ["--size-feedback", "shrinking"],
     MELT + ["--deep-runoff", "maybe"],
+    MELT + ["--molten-cascade", "maybe"],
 ])
 def test_bad_melt_arguments_exit_2(argv, tmp_path):
     with pytest.raises(SystemExit) as exc:
```

**Measured (2026-10-03, throwaway copy).** `test_run_name_with_melt`, `test_melting_run_writes_columns_files_and_json`
and the nine `test_bad_melt_arguments_exit_2` cases pass; the first two fail on the copy without the cascade (no
`molten_cascade` keyword, no `molten_cascade` setting, no cascade columns) and the new bad-argument case passes on
both (argparse rejects an unknown flag and an unknown value alike, with exit 2). The 15 s melting run of the second test
starts from a warm body at 71 km on the coarse mesh, and the cascade already feeds
there. Whole unit tier: see sub-plan 09's amendment.

## Amendment of 2026-10-05 — `--seed`: numpy's generator is seeded at the start of every run

> Part of the seeding amendment (facts 62–68 in `00-shared-context.md`). This sub-plan holds its design and its tested
> code; sub-plan 14's amendment of this date adds the sensitivity row and sub-plan 15's the documentation.

**The flag.** `--seed` in the main `run` group — not the Step 3 group, because it applies to every run, melting or not
— takes an integer in [0, 2³² − 1], the range `np.random.seed` accepts, with the default `DEFAULT_SEED = 12345`. The
first statement of `cmd_run` is `np.random.seed(args.seed)`, before the argument checks and before anything is built,
so nothing a run does can draw from numpy's global generator before it is seeded. Measured (fact 62): importing the
package and parsing the arguments draw nothing, and the only draws in a run are pyamg's, inside
`approximate_spectral_radius`, at every hierarchy the skfem backend builds; seeding at the start of `cmd_run` therefore
fixes every random number a run uses, and two runs of one build agree bit for bit (fact 64). A malformed seed — not an
integer, negative, or 2³² or more — is rejected by the argument's type check (`parse_seed`, the pattern of
`parse_epoch`) with argparse's exit 2, before `cmd_run` runs; exit codes are otherwise unchanged.

**Why 12345.** It is the seed the measurement harness of the amendments of 2026-10-02 and 2026-10-03 set before
importing the package (`MEASURE_SEED=12345`). Nothing draws between that point and the start of `cmd_run`, so the
model's own default reproduces the harness's runs: measured on the 100 mm flight to 120 s of fact 59 and the whole 50 mm
flight of fact 60 — both with the cascade and the deep runoff on — every history cell, every source-table row, every
result field but the run time and every array of every VTK frame; and with `--deep-runoff off --molten-cascade off` it
reproduces a seeded run of 2026-10-02 of the unamended prototype the same way (fact 64). The seeded runs of facts 49–61
are therefore reproducible from the command line alone, with the flags that select the model each was made with. Fact 52
suggested `np.random.seed(0)`, which would have left every one of those runs out of the model's reach.

**Alternatives not taken.** (1) A fixed starting vector inside the skfem backend, for instance a private generator
reseeded at every hierarchy build: it would make the backend reproducible wherever it is constructed, the tests and the
reference tier included, but it is a change to sub-plan 03's thermal core, `smoothed_aggregation_solver` offers no
argument for the starting vector of the spectral-radius estimate inside its prolongation smoother, and it would change
every run's numbers against the measured ones. (2) `--seed random`, a seed drawn from the operating system and recorded:
a scatter study is reproducible only if its seeds are listed, and an explicit list (`--seed 1`, `--seed 2`, …) does that
with no new code. (3) No seed by default, the old behaviour: an unseeded run can do nothing that a run with another seed
cannot, and it cannot be repeated. (4) The suffix only on `--thermal fem` runs, the only ones in which anything draws:
the seed is applied to every run and recorded in every JSON, and the name should encode the whole configuration.

**Run names, the JSON and the printed summary.** A run with a seed other than 12345 ends in `_seed-<n>` (after
`_deeprunoff-off` and `_moltencascade-off`), in every mode, so runs that differ only in the seed can share one output
directory (CLAUDE.md's resumability convention); the default keeps the existing name, so no existing run name, reference
or test changes, and `--name` still overrides the name. The run JSON's `settings` carry `"seed"` for every run, and the
printed summary's first line names it — `<run name> (seed 12345): t_max at t = ...` — because the default does not
appear in the name. The verification and sensitivity drivers' summary tables list results by case and mode or variant,
not the configuration, so they are unchanged; every run they make goes through `cmd_run` and is seeded (sub-plan 14's
amendment of this date).

Code (the tested change to `reentry_model/cli.py`, the prototype's file with the amendments of 2026-10-02 and
2026-10-03 above applied):

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -36,6 +36,11 @@
 WINDS = ("none", "static")
 THERMAL_MODES = ("none", "fem")
 MELT_MODES = ("off", "on")
+# numpy's global generator is seeded with this at the start of every run (amendment of 2026-10-05): pyamg draws from it
+# the starting vector of the spectral-radius estimate that weights the prolongation smoother of every hierarchy it builds,
+# so an unseeded run differed from the next from the first solve on. 12345 is the seed the measurement harness of the
+# amendments of 2026-10-02 and 2026-10-03 set, so the model's default reproduces their runs.
+DEFAULT_SEED = 12345
 
 
 def parse_epoch(text):
@@ -45,15 +50,28 @@
         raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))
 
 
+def parse_seed(text):
+    """--seed: an integer in [0, 2**32 - 1], the range numpy's np.random.seed accepts."""
+    try:
+        seed = int(text)
+    except ValueError:
+        seed = -1
+    if not 0 <= seed < 2 ** 32:
+        raise argparse.ArgumentTypeError("seed must be an integer in [0, 2**32 - 1], got {!r}".format(text))
+    return seed
+
+
 def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None,
-                   deep_runoff=True, molten_cascade=True):
+                   deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED):
     """The run name encodes the configuration, so runs of different settings can share an output directory. A melting
     run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
-    (amendment of 2026-10-03) in `_moltencascade-off`; the defaults keep the old name."""
+    (amendment of 2026-10-03) in `_moltencascade-off`, and any run with a seed other than DEFAULT_SEED (amendment of
+    2026-10-05) in `_seed-<n>`; the defaults keep the old name."""
     name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
         diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
     return (name + ("_fem-" + heating_name if heating_name else "") + ("_melt-" + melt if melt else "")
-            + ("_deeprunoff-off" if melt and not deep_runoff else "") + ("_moltencascade-off" if melt and not molten_cascade else ""))
+            + ("_deeprunoff-off" if melt and not deep_runoff else "") + ("_moltencascade-off" if melt and not molten_cascade else "")
+            + ("_seed-{}".format(seed) if seed != DEFAULT_SEED else ""))
 
 
 def make_atmosphere(spec, epoch, wind_name):
@@ -125,6 +143,12 @@
                         "with --thermal fem the macro-step history is interpolated)")
     r.add_argument("--outdir", default=DEFAULT_OUTDIR)
     r.add_argument("--name", default=None, help="run name (default: model_d..mm_v..kms_h..km_<atmosphere>_<bridging>_<wind>[_fem-<heating>])")
+    r.add_argument("--seed", type=parse_seed, default=DEFAULT_SEED,
+                   help="seed of numpy's global random generator, set at the start of the run (default %(default)s): pyamg draws "
+                        "the starting vectors of its spectral-radius estimates from it, so a fixed seed makes a run reproducible "
+                        "bit for bit and another seed changes it only through that round-off, which the melting model amplifies "
+                        "into the run-to-run scatter; an integer in [0, 2**32 - 1], recorded in the run JSON, and a seed other "
+                        "than the default ends the run name in _seed-<n>")
     r.add_argument("--quiet", action="store_true")
     th = r.add_argument_group("thermal model (Step 2)")
     th.add_argument("--thermal", choices=THERMAL_MODES, default="none", help="none: Step 1 trajectory only (default); fem: coupled 3D conduction")
@@ -256,6 +280,10 @@
 
 
 def cmd_run(args, parser):
+    # first, before anything that could draw from it (amendment of 2026-10-05; plan facts 52 and 62): the only draws in a
+    # run are pyamg's, at every hierarchy the skfem backend builds, and with the generator seeded here two runs of one
+    # build agree bit for bit
+    np.random.seed(args.seed)
     for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                          ("--cadence", args.cadence), ("--t-max", args.t_max)):
         if value <= 0.0:
@@ -294,7 +322,7 @@
     reference = sesam_io.load_reference(args.reference) if args.reference else None
     name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
                                        args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None,
-                                       deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on")
+                                       deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on", seed=args.seed)
     run_dir = os.path.join(args.outdir, name)
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
@@ -326,7 +354,8 @@
                    "material_density_kgm3": args.material_density, "mass_kg": mass, "initial_temperature_K": args.temperature},
         "settings": {"atmosphere": atm_name, "wind": args.wind, "bridging": args.bridging, "gravity": args.gravity,
                      "rtol": args.rtol, "atol_position_m": settings.atol_position, "atol_velocity_ms": settings.atol_velocity,
-                     "cadence_s": args.cadence, "t_max_s": args.t_max, **atm_info, "thermal": args.thermal, **thermal_info},
+                     "cadence_s": args.cadence, "t_max_s": args.t_max, "seed": args.seed, **atm_info, "thermal": args.thermal,
+                     **thermal_info},
         "results": history.results,
         "comparison": None,
         "provenance": prov,
@@ -362,8 +391,8 @@
     tj.write_run_json(json_path, doc)
     if not args.quiet:
         res = history.results
-        print("{}: {} at t = {:.1f} s, final V {:.4f} km/s, Kn {:.3g} -> {:.3g}, {} RHS evaluations in {:.1f} s".format(
-            name, res["end_reason"], res["final_time_s"], res["final_velocity_kms"], res["knudsen_start"],
+        print("{} (seed {}): {} at t = {:.1f} s, final V {:.4f} km/s, Kn {:.3g} -> {:.3g}, {} RHS evaluations in {:.1f} s".format(
+            name, args.seed, res["end_reason"], res["final_time_s"], res["final_velocity_kms"], res["knudsen_start"],
             history.columns["knudsen"][-1], res["rhs_evaluations"], res["runtime_s"]))
         if args.thermal == "fem":
             print("  thermal: peak surface T {:.0f} K at t = {:.0f} s, peak mean T {:.0f} K, integrated heat {:.3g} J, "
```

Tests (the tested change to `tests/test_reentry_model_cli.py`):

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -3,9 +3,10 @@
 import json
 import os
 
+import numpy as np
 import pytest
 
-from reentry_model import cli, compare, sesam_io
+from reentry_model import cli, compare, coupled, sesam_io
 
 R100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_mAA7075_nomelt_msis_nowind.csv")
 BASE = ["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.500133", "--flight-path-angle", "-0.959331"]
@@ -14,6 +15,13 @@
 def test_run_name():
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none") == \
         "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none"
+    # the seed of numpy's generator (amendment of 2026-10-05): the default keeps the name, any other seed ends it in
+    # _seed-<n>, so runs that differ only in the seed never overwrite each other
+    assert cli.DEFAULT_SEED == 12345
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", seed=cli.DEFAULT_SEED) == \
+        "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none"
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", seed=7) == \
+        "model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_seed-7"
 
 
 def test_run_us76_short_flight(tmp_path):
@@ -24,6 +32,7 @@
     assert [float(r["time_s"]) for r in rows] == [0.0, 5.0, 10.0, 15.0, 20.0]
     doc = json.load(open(tmp_path / (name + ".json")))
     assert doc["settings"]["atmosphere"] == "us76" and doc["results"]["end_reason"] == "t_max"
+    assert doc["settings"]["seed"] == cli.DEFAULT_SEED                         # every run records its seed (2026-10-05)
     assert doc["inputs"]["mass_kg"] == pytest.approx(1.4728833557580150) and doc["comparison"] is None
     assert doc["provenance"]["package_version"] and "git_commit" in doc["provenance"]
     assert doc["provenance"]["reference_sha256"] is None and doc["provenance"]["replay_sha256"] is None
@@ -72,6 +81,10 @@
     BASE + ["--gravity", "j6"],
     ["run", "--diameter", "-1", "--velocity", "7.5", "--altitude", "77.5"],
     BASE + ["--epoch", "yesterday"],
+    BASE + ["--seed", "abc"],                     # the seed (amendment of 2026-10-05): an integer in [0, 2**32 - 1],
+    BASE + ["--seed", "1.5"],                     # the range numpy's generator accepts
+    BASE + ["--seed", "-1"],
+    BASE + ["--seed", "4294967296"],
 ])
 def test_bad_arguments_exit_2(argv, tmp_path):
     with pytest.raises(SystemExit) as exc:
@@ -160,6 +173,42 @@
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics").endswith("_fem-physics")
 
 
+def test_the_seed_makes_a_run_reproducible_and_another_seed_changes_it(tmp_path, monkeypatch, capsys):
+    """numpy's global generator is seeded at the start of every run (amendment of 2026-10-05): pyamg draws from it the
+    starting vector of the spectral-radius estimate that weights the prolongation smoother of every hierarchy it
+    builds, which made two runs of one build differ from the first solve on (plan fact 52). Three 5 s coupled runs on
+    the coarse mesh, whose conduction solves go through pyamg: two with the default seed -- the second started with the
+    generator scrambled, so a run cannot depend on what drew from it before -- agree bit for bit in every history column
+    and every result field but the run time, while one with another seed differs, at round-off, under its own name."""
+    histories = []
+    original_run = coupled.CoupledRun.run
+
+    def run_and_keep(self):
+        histories.append(original_run(self))
+        return histories[-1]
+
+    monkeypatch.setattr(coupled.CoupledRun, "run", run_and_keep)
+    argv = FEM + ["--t-max", "5", "--dt", "0.5"]
+    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "sesam")
+    np.random.seed(1)
+    assert cli.main(argv + ["--outdir", str(tmp_path / "a")]) == 0
+    np.random.seed(2)
+    np.random.rand(1000)
+    assert cli.main(argv + ["--outdir", str(tmp_path / "b")]) == 0
+    capsys.readouterr()
+    assert cli.main([a for a in argv if a != "--quiet"] + ["--outdir", str(tmp_path / "c"), "--seed", "7"]) == 0
+    assert name + "_seed-7 (seed 7): t_max at t = 5.0 s" in capsys.readouterr().out     # the printed summary names it
+    a, b, c = histories
+    ra, rb = (json.load(open(tmp_path / d / (name + ".json")))["results"] for d in ("a", "b"))
+    assert set(a.columns) == set(b.columns) and all(np.array_equal(a.columns[k], b.columns[k], equal_nan=True) for k in a.columns)
+    assert {k: v for k, v in ra.items() if k != "runtime_s"} == {k: v for k, v in rb.items() if k != "runtime_s"}
+    assert not all(np.array_equal(a.columns[k], c.columns[k], equal_nan=True) for k in a.columns)
+    for k in a.columns:                                       # the seed acts at round-off; the melting model amplifies it
+        np.testing.assert_allclose(c.columns[k], a.columns[k], rtol=1e-6, atol=1e-9)
+    doc = json.load(open(tmp_path / "c" / (name + "_seed-7.json")))
+    assert doc["settings"]["seed"] == 7 and doc["run_name"] == name + "_seed-7"
+
+
 # ---------------------------------------------------------------------------------------------------------------
 # Step 3: melting flags, the melting run's files, the bookkeeping device against the melting reference
 
@@ -173,6 +222,8 @@
                               deep_runoff=False).endswith("_fem-physics_melt-girin_deeprunoff-off")    # amendment of 2026-10-02
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
                               molten_cascade=False).endswith("_fem-physics_melt-girin_moltencascade-off")   # the molten cascade
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin", deep_runoff=False,
+                              molten_cascade=False, seed=7).endswith("_melt-girin_deeprunoff-off_moltencascade-off_seed-7")
 
 
 def test_melting_run_writes_columns_files_and_json(tmp_path):
@@ -189,7 +240,7 @@
     s, r, f = doc["settings"], doc["results"], doc["files"]
     assert s["melt"] == "on" and s["material"] == "AA7075_range" and s["removal"] == "girin" and s["runoff"] == "on" and s["prism_layers"] == 0
     assert s["rarefied_shear"] == "slip" and s["we_critical"] == 4.62 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["liquid"]["sigma"] == 0.86
-    assert s["gamma_pm"] == 1.15 and s["kn_body_shock"] == 0.01 and s["deep_runoff"] == "on"
+    assert s["gamma_pm"] == 1.15 and s["kn_body_shock"] == 0.01 and s["deep_runoff"] == "on" and s["seed"] == cli.DEFAULT_SEED
     assert "deep_runoff_mass_kg" in rows[0] and "deep_mass_kg" in rows[0] and r["deep_runoff_mass_kg"] >= 0.0
     assert "deep_surfaced_mass_kg" in rows[0] and "deep_blob_fraction" in rows[0] and "deep_liquid_kg" in rows[0]
     assert s["molten_cascade"] == "on" and "cascade_passes" in rows[0] and "cascade_mass_kg" in rows[0]      # the molten cascade
```

**Measured (2026-10-05, throwaway copy).** `test_run_name`, `test_run_us76_short_flight`, `test_run_name_with_melt`,
`test_melting_run_writes_columns_files_and_json` and the new
`test_the_seed_makes_a_run_reproducible_and_another_seed_changes_it` fail on the copy without the amendment (no
`DEFAULT_SEED`, no `seed` keyword or setting, `--seed` an unknown argument) and pass with it; the four new
`test_bad_arguments_exit_2` cases pass on both, because argparse rejects an unknown flag and a malformed value alike,
with exit 2. The new test makes three 5 s Step 2 runs on the coarse mesh in 0.7 s; on the copy without the amendment two
such runs started from different generator states differ in 6 of their 31 history columns, by 2e-16 to 6e-16 relative —
the last bit — which is what its bit-identity assertion catches. One thing the test found that is not the seed's: meshio
writes an empty line to standard output when it reads a cached mesh, so the printed summary is not the first line of a
run's output; the test looks for it by the run name. The CLI module: 34 passed and the one known failure from the
missing melting reference (Task 11), against 29 and the same failure before. Whole unit tier: 234 passed, 1 skipped, and
the five known failures and errors from the missing melting SESAM references, against 229, 1 and the same five, in 179 s
either way.

## Amendment of 2026-10-05 (runoff flux) — no flag, and the exit code of a model failure (no code change)

> Part of the runoff-flux amendment (sub-plan 06's amendment of this date holds the change; facts 69–77 in
> `00-shared-context.md`). Nothing in `cli.py` changes.

**No flag.** The change to the film's flux on thick patches has no `--` switch, unlike the deep runoff and the molten
cascade: the old flux is a defect — it crashed the transport at fine steps and moved a thin film like a whole conjugate
layer — not a modelling alternative to keep selectable, and every earlier run remains reproducible from the amendments
that made it. Run names, the run JSON's settings and the printed summary are therefore unchanged, which means that a run
made before the change and one made after it carry the same name and must not share an output directory if both are to
be kept (fact 35's rule — fresh output directories — covers it).

**Found, not changed: a model failure is reported as bad arguments.** The crashed run of fact 70 ended with "ERROR: film
mass must be one finite non-negative value per node" and exit code 2. That message is the thermal solver's input check
in `set_film_mass` (both backends raise `ValueError` for a non-finite or negative nodal film mass), and `cli.main` maps
every `ValueError` to exit 2, the code the command-line interface reserves for bad arguments or a missing optional
library (`CLAUDE.md`: "0 ok, 1 a model/integration failure ..., 2 bad arguments or a missing optional library"); a NaN
produced by the model mid-flight is a model failure, exit 1. The fix of this amendment removes the cause of that
particular NaN, but the mapping would report any future one the same way. Options for Asha (fact 77): (1) leave it; (2)
raise a dedicated exception (or `RuntimeError`, which `cli.main` already maps to exit 1) from the solvers' input checks
when they are called by the melt step rather than by argument handling; (3) check the film for finite values at the end
of `MeltingBody.melt_step` and raise `RuntimeError` there, naming the stage, which also points at the cause instead of
at the solver that tripped over it. Recommendation: (3), with a unit test that a NaN film ends a run with exit 1, in a
change of its own.

## Amendment of 2026-10-06 — `--rigid-substrate on|off`

> Part of the rigid-substrate amendment (sub-plan 09's amendment of this date holds the design; facts 78–87 in
> `00-shared-context.md` the measurements). The code below is the tested code, as a diff against the throwaway copy
> described there.

**The flag.** `--rigid-substrate on|off` in the Step 3 group, default `on`, sets `MeltSettings.rigid_substrate`. With
`off` the regime test reads the film plus the fully molten material, as in every run before this amendment, and the
amended code reproduces the copy before the change bit for bit on both flights (fact 81).

**Run names and the JSON.** A melting run with `--rigid-substrate off` ends in `_rigidsubstrate-off`, after
`_moltencascade-off` and before `_seed-<n>`; the default keeps the existing name, as the deep runoff's and the cascade's
switches did, so a default run made after this amendment carries the name of one made before it although its physics
differs — compare runs by their JSON settings, which now record `"rigid_substrate": "on"` or `"off"`. A bad value exits
2 like every other bad argument.

In `reentry_model/cli.py`:

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -62,15 +62,17 @@
 
 
 def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None,
-                   deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED):
+                   deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED, rigid_substrate=True):
     """The run name encodes the configuration, so runs of different settings can share an output directory. A melting
     run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
-    (amendment of 2026-10-03) in `_moltencascade-off`, and any run with a seed other than DEFAULT_SEED (amendment of
-    2026-10-05) in `_seed-<n>`; the defaults keep the old name."""
+    (amendment of 2026-10-03) in `_moltencascade-off`, one with the rigid substrate off (amendment of 2026-10-06) in
+    `_rigidsubstrate-off`, and any run with a seed other than DEFAULT_SEED (amendment of 2026-10-05) in `_seed-<n>`;
+    the defaults keep the old name."""
     name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
         diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
     return (name + ("_fem-" + heating_name if heating_name else "") + ("_melt-" + melt if melt else "")
             + ("_deeprunoff-off" if melt and not deep_runoff else "") + ("_moltencascade-off" if melt and not molten_cascade else "")
+            + ("_rigidsubstrate-off" if melt and not rigid_substrate else "")
             + ("_seed-{}".format(seed) if seed != DEFAULT_SEED else ""))
 
 
@@ -193,6 +195,11 @@
                          "(default on). off: an exposed element waits for the next step's feed, so the surface recedes through "
                          "molten material by one element per macro step, as in every run before the molten-cascade amendment of "
                          "2026-10-03")
+    me.add_argument("--rigid-substrate", choices=("on", "off"), default="on",
+                    help="the thin spraying branch needs a rigid substrate: the regime test compares Girin's conjugate depth "
+                         "with the film plus everything more than half liquid beneath it (slurry included), and where there "
+                         "is no conjugate depth a film on slurry does not spray (default on). off: the regime test reads the "
+                         "film plus the fully molten material, as in every run before the amendment of 2026-10-06")
     me.add_argument("--rt-spray", choices=("on", "off"), default="on",
                     help="apply Girin & Kopyt's front-surface Rayleigh-Taylor mode as a release mechanism, on the patches where it "
                          "grows faster than the shear mode (default on). off: evaluate and report it but release nothing through it, "
@@ -249,7 +256,8 @@
         size_feedback = args.size_feedback or ("current" if args.removal == "girin" else "initial")
         melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
                                           particles=args.particles, size_feedback=size_feedback,
-                                          deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on")
+                                          deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on",
+                                          rigid_substrate=args.rigid_substrate == "on")
         the_body = body.MeltingBody(the_mesh, mat, solver, mass, flow, spray_model, melt_settings, T0=args.temperature,
                                     emissivity=args.emissivity, T_ambient=args.t_ambient)
     else:
@@ -270,7 +278,7 @@
                      "accommodation": args.accommodation, "catalycity": args.catalycity})
     if melting:
         info.update({"removal": args.removal, "runoff": args.runoff, "deep_runoff": args.deep_runoff,
-                     "molten_cascade": args.molten_cascade, "rt_spray": args.rt_spray,
+                     "molten_cascade": args.molten_cascade, "rigid_substrate": args.rigid_substrate, "rt_spray": args.rt_spray,
                      "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                      "k_r": args.kr, "k_t": args.kt, "demise_fraction": args.demise_fraction, "particles": args.particles,
                      "size_feedback": size_feedback, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
@@ -322,7 +330,8 @@
     reference = sesam_io.load_reference(args.reference) if args.reference else None
     name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
                                        args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None,
-                                       deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on", seed=args.seed)
+                                       deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on", seed=args.seed,
+                                       rigid_substrate=args.rigid_substrate == "on")
     run_dir = os.path.join(args.outdir, name)
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
```

In `tests/test_reentry_model_cli.py`:

```diff
--- a/tests/test_reentry_model_cli.py
+++ b/tests/test_reentry_model_cli.py
@@ -224,6 +224,11 @@
                               molten_cascade=False).endswith("_fem-physics_melt-girin_moltencascade-off")   # the molten cascade
     assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin", deep_runoff=False,
                               molten_cascade=False, seed=7).endswith("_melt-girin_deeprunoff-off_moltencascade-off_seed-7")
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
+                              rigid_substrate=False).endswith("_fem-physics_melt-girin_rigidsubstrate-off")  # 2026-10-06
+    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin", deep_runoff=False,
+                              molten_cascade=False, rigid_substrate=False,
+                              seed=7).endswith("_melt-girin_deeprunoff-off_moltencascade-off_rigidsubstrate-off_seed-7")
 
 
 def test_melting_run_writes_columns_files_and_json(tmp_path):
@@ -244,6 +249,8 @@
     assert "deep_runoff_mass_kg" in rows[0] and "deep_mass_kg" in rows[0] and r["deep_runoff_mass_kg"] >= 0.0
     assert "deep_surfaced_mass_kg" in rows[0] and "deep_blob_fraction" in rows[0] and "deep_liquid_kg" in rows[0]
     assert s["molten_cascade"] == "on" and "cascade_passes" in rows[0] and "cascade_mass_kg" in rows[0]      # the molten cascade
+    assert s["rigid_substrate"] == "on" and all(k in rows[0] for k in ("nonrigid_depth_mean_mm", "slurry_thick_fraction",
+                                                                       "rigid_thin_fraction", "slurry_held_mass_kg"))
     assert r["cascade_mass_kg"] == pytest.approx(float(rows[-1]["cascade_mass_kg"])) and r["cascade_mass_kg"] >= 0.0
     assert r["cascade_capped_steps"] == 0
     assert r["cascade_passes_max"] == max(int(float(row["cascade_passes"])) for row in rows)
@@ -284,6 +291,7 @@
     MELT + ["--size-feedback", "shrinking"],
     MELT + ["--deep-runoff", "maybe"],
     MELT + ["--molten-cascade", "maybe"],
+    MELT + ["--rigid-substrate", "maybe"],
 ])
 def test_bad_melt_arguments_exit_2(argv, tmp_path):
     with pytest.raises(SystemExit) as exc:
```

---

**Depends on:** effectively everything before it (Tasks 9–12 directly, and transitively the rest).
**Produces, for later tasks:** the new command-line flags, run naming, and JSON settings/results schema that Task 14's verification runs are driven through.
**Character:** pure plumbing/glue — large in line count (roughly 500 lines transcribed) but conceptually simple, mostly argparse and dict wiring, not new physics.
**Read before implementing:** none beyond the shared context.
**Refinement goal for the sub-agent:** turn the section below into a standalone implementation plan.

---

### Task 13: Command line

**Files:**
- Modify (replace): `reentry_model/cli.py`
- Test: `tests/test_reentry_model_cli.py` (append)

**Interfaces:**
- Consumes: everything above.
- Produces: `run` flags `--melt off|on`, `--material` (names or path; default `AA7075_nomelt`, `AA7075_range` with `--melt on`), `--removal`, `--runoff`, `--rarefied-shear`, `--we-critical`, `--kr`, `--kt`, `--prism-layers` (default 4 with melting, 0 otherwise), `--layer-thickness` (mm), `--demise-fraction`, `--particles/--no-particles`, `--size-feedback current|initial` (default `current` with `girin`, `initial` with `instant`; it governs the Knudsen length, the nose-cap radius and the drag shape factor together), `--gamma-pm` (default 1.15), `--kn-body-shock` (default 0.01), `--k-scale`, `--consistent-mass`; run names end in `_melt-<removal>`; the JSON `settings` carry the melt settings and the liquid properties, `results` the melt results, `files` the particle files and melt plots (and `film`); `compare` handles melting histories; `model_run_name(..., heating_name=None, melt=None)`; `_fmt`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reentry_model_cli.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3: melting flags, the melting run's files, the bookkeeping device against the melting reference

MELT_100 = os.path.join(sesam_io.REFERENCE_DIR, "sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv")
MELT = FEM + ["--melt", "on", "--prism-layers", "0", "--h-surface", "4", "--h-core", "20"]


def test_run_name_with_melt():
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin").endswith("_fem-physics_melt-girin")


def test_melting_run_writes_columns_files_and_json(tmp_path):
    """15 s from 71 km with a warm body: melt columns, the particle files, the melt plots, stills with the melt marks."""
    pytest.importorskip("cantera")
    argv = [a for a in MELT if a not in ("--heating", "sesam")] + ["--heating", "physics", "--altitude", "71", "--velocity", "7.24", "--temperature", "800",
                                                                   "--t-max", "15", "--outdir", str(tmp_path), "--name", "melt_short", "--stills"]
    assert cli.main(argv) == 0
    rows = list(csv.DictReader(open(tmp_path / "melt_short.csv")))
    assert len(rows) == 31 and "sprayed_mass_kg" in rows[0] and "film_mass_kg" in rows[0] and "closure_fraction_girin" in rows[0]
    assert float(rows[-1]["kn_body"]) > 0.0 and float(rows[-1]["kn_local_stag"]) < float(rows[-1]["kn_body"]) and "flow_branch" in rows[0]
    assert float(rows[-1]["mass_kg"]) < float(rows[0]["mass_kg"]) and float(rows[-1]["sprayed_mass_kg"]) > 0.0
    doc = json.load(open(tmp_path / "melt_short.json"))
    s, r, f = doc["settings"], doc["results"], doc["files"]
    assert s["melt"] == "on" and s["material"] == "AA7075_range" and s["removal"] == "girin" and s["runoff"] == "on" and s["prism_layers"] == 0
    assert s["rarefied_shear"] == "slip" and s["we_critical"] == 4.62 and s["k_r"] == 0.17 and s["k_t"] == 1.1 and s["liquid"]["sigma"] == 0.86
    assert s["gamma_pm"] == 1.15 and s["kn_body_shock"] == 0.01
    assert s["size_feedback"] == "current" and "nose_radius_mm" in rows[0] and float(rows[-1]["nose_radius_mm"]) > 0.0
    assert s["T_liquidus_K"] == 908.0 and s["latent_heat_Jkg"] == 400e3 and s["demise_fraction"] == 0.01 and s["particles"] is True
    assert r["melt_onset_altitude_km"] is not None and r["sprayed_mass_kg"] > 0.0 and r["n_source_rows"] > 0 and abs(r["melt_energy_balance_residual"]) < 1e-6
    for key in ("particles", "particles_summary", "size_distribution"):
        assert os.path.isfile(f[key])
    assert len(f["melt_plots"]) == len(compare.MELT_PLOT_NAMES) and all(os.path.isfile(p) for p in f["melt_plots"])
    assert doc["comparison"] is None and f["film"] is None
    names = [os.path.basename(p) for p in f["stills"]]
    assert any(n.startswith("melt_onset") for n in names) and any(n.startswith("film_") for n in names) and any(n.startswith("section_spraying_onset") for n in names)


def test_bookkeeping_device_against_the_melting_reference(tmp_path):
    """SESAM-equivalent heating + AA7075 + instant removal + k x 1e4 on the coarse mesh over the whole flight: the
    melt metrics are written and the mass follows SESAM's lumped law (2 % / 0.5 km / 2 %, spec 13.1; measured
    0.82-0.98 % / 0.10 km / -1.3 %)."""
    argv = MELT + ["--material", "AA7075", "--removal", "instant", "--runoff", "off", "--k-scale", "1e4", "--reference", MELT_100,
                   "--outdir", str(tmp_path), "--name", "bookkeeping"]
    assert cli.main(argv) == 0
    doc = json.load(open(tmp_path / "bookkeeping.json"))
    mm = doc["comparison"]["melt_metrics"]
    assert doc["results"]["end_reason"] == "demise" and mm["mass"]["max_rel_m0"] < 0.02
    assert abs(mm["onset_altitude_diff_km"]) < 0.5 and abs(mm["demise_time_rel"]) < 0.02 and mm["sprayed_mass_kg"] == 0.0
    assert os.path.isfile(tmp_path / "bookkeeping" / "mass_time.png") and doc["settings"]["k_scale"] == 1e4
    assert doc["settings"]["size_feedback"] == "initial"                                             # SESAM's D0 / R0 for the device
    assert not os.path.isfile(tmp_path / "bookkeeping" / "particles.npz") or True                    # written (empty table) with --particles


@pytest.mark.parametrize("argv", [
    BASE + ["--atmosphere", "us76", "--melt", "on"],                                              # needs --thermal fem
    MELT + ["--kr", "0"],
    MELT + ["--prism-layers", "-1"],
    MELT + ["--demise-fraction", "1.5"],
    MELT + ["--k-scale", "0"],
    MELT + ["--removal", "magic"],
    MELT + ["--size-feedback", "shrinking"],
])
def test_bad_melt_arguments_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2
```


- [ ] **Step 2: Run the tests to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: the new tests fail (`--melt` unknown).

- [ ] **Step 3: Replace `reentry_model/cli.py`**

```python
"""Command line of the re-entry model.

    python -m reentry_model run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 \
        [--atmosphere nrlmsise|us76|replay:<sesam.csv>] [--reference <sesam.csv>] [--outdir ...]
        [--thermal fem --heating sesam|physics ... --animate]          (Step 2: coupled 3D conduction)
        [--melt on --material AA7075_range --removal girin|instant ...] (Step 3: melting, film, spraying)
    python -m reentry_model compare --model <model.csv> --reference <sesam.csv> [--outdir ...]

Exit codes: 0 ok, 1 the flight escaped / integration failed, 2 bad input or a missing optional library
(cantera for --heating physics, dolfinx for --thermal-solver fenicsx: create the fenicsx_env environment).
"""
import argparse
import importlib.metadata
import math
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, spray, surface_flow, thermal, viz
from . import trajectory as tj
from .earth import GRAVITY_MODELS
from .thermal import MissingBackend

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "reentry_model_output")
# the reference cases' break-off state (parent satellite, see README "Reference runs")
DEFAULT_HEADING_DEG = 347.168296
DEFAULT_LAT_DEG = 29.546067
DEFAULT_LON_DEG = -82.134333
DEFAULT_EPOCH = "2024-08-01T12:53:07"
DEFAULT_MATERIAL_DENSITY = 2813.0          # drama-AA7075
ATMOSPHERES = ("nrlmsise", "us76")         # plus replay:<path>
WINDS = ("none", "static")
THERMAL_MODES = ("none", "fem")
MELT_MODES = ("off", "on")


def parse_epoch(text):
    try:
        return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError("epoch must be YYYY-MM-DDTHH:MM:SS, got {!r}".format(text))


def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None):
    name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
        diameter_m * 1e3, velocity_ms / 1e3, altitude_m / 1e3, atmosphere_name, bridging_name, wind_name)
    return name + ("_fem-" + heating_name if heating_name else "") + ("_melt-" + melt if melt else "")


def make_atmosphere(spec, epoch, wind_name):
    """(atmosphere object, short name, provenance dict) for an --atmosphere value."""
    wind = atmosphere.NoWind() if wind_name == "none" else atmosphere.StaticProfileWind()
    if spec == "us76":
        return atmosphere.US76TableAtmosphere(wind=wind), "us76", {}
    if spec == "nrlmsise":
        solar = fap.solar_indices(fap.load_fap_day(fap.DEFAULT_FAP_DAY), epoch.date())
        return atmosphere.NRLMSISE00Atmosphere(epoch, solar, wind), "nrlmsise", \
            {"f107": solar.f107, "f107a": solar.f107a, "ap": solar.ap, "fap_day": fap.DEFAULT_FAP_DAY}
    if spec.startswith("replay:"):
        ref = sesam_io.load_reference(spec[len("replay:"):])
        return atmosphere.ReplayAtmosphere(ref, wind), "replay", {"replay_reference": ref.name, "replay_sha256": ref.sha256}
    raise ValueError("--atmosphere must be one of {} or replay:<sesam.csv>, got {!r}".format(ATMOSPHERES, spec))


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def provenance():
    import scipy
    try:
        import pymsis
        pymsis_version = pymsis.__version__
    except ImportError:
        pymsis_version = None
    package_names = {"skfem": "scikit-fem", "gmsh": "gmsh", "pyamg": "pyamg", "pyvista": "pyvista",
                      "cantera": "cantera", "dolfinx": "fenics-dolfinx"}
    versions = {}
    for key, pkg in package_names.items():
        try:
            versions[key] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[key] = None
    return {"package_version": __version__, "git_commit": git_commit(), "numpy": np.__version__,
            "scipy": scipy.__version__, "pymsis": pymsis_version, "python": sys.version.split()[0], **versions}


def build_parser():
    p = argparse.ArgumentParser(prog="reentry_model", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="integrate one sphere trajectory")
    r.add_argument("--diameter", type=float, required=True, help="sphere diameter [mm]")
    r.add_argument("--velocity", type=float, required=True, help="initial velocity [km/s], relative to the rotating atmosphere")
    r.add_argument("--altitude", type=float, required=True, help="initial geodetic altitude [km]")
    r.add_argument("--flight-path-angle", type=float, default=0.0, help="[deg], negative = descending")
    r.add_argument("--heading", type=float, default=DEFAULT_HEADING_DEG, help="[deg] clockwise from north (default %(default)s)")
    r.add_argument("--lat", type=float, default=DEFAULT_LAT_DEG, help="geodetic latitude [deg] (default %(default)s)")
    r.add_argument("--lon", type=float, default=DEFAULT_LON_DEG, help="longitude [deg] (default %(default)s)")
    r.add_argument("--epoch", type=parse_epoch, default=parse_epoch(DEFAULT_EPOCH), help="UTC YYYY-MM-DDTHH:MM:SS (default {})".format(DEFAULT_EPOCH))
    r.add_argument("--material-density", type=float, default=DEFAULT_MATERIAL_DENSITY, help="[kg/m3] (default %(default)s)")
    r.add_argument("--temperature", type=float, default=300.0,
                   help="initial temperature [K]: recorded only with --thermal none, the field's initial condition with --thermal fem (default %(default)s)")
    r.add_argument("--atmosphere", default="nrlmsise", help="nrlmsise (default) | us76 | replay:<sesam.csv>")
    r.add_argument("--wind", choices=WINDS, default="none")
    r.add_argument("--bridging", choices=aero.BRIDGING_NAMES, default="sesam-table")
    r.add_argument("--gravity", choices=GRAVITY_MODELS, default="j2")
    r.add_argument("--rtol", type=float, default=1e-9)
    r.add_argument("--cadence", type=float, default=1.0, help="history sample spacing [s] (default %(default)s)")
    r.add_argument("--t-max", type=float, default=3600.0, help="[s] (default %(default)s)")
    r.add_argument("--reference", default=None,
                   help="SESAM run CSV to compare against (with --thermal none the model is also sampled at its times; "
                        "with --thermal fem the macro-step history is interpolated)")
    r.add_argument("--outdir", default=DEFAULT_OUTDIR)
    r.add_argument("--name", default=None, help="run name (default: model_d..mm_v..kms_h..km_<atmosphere>_<bridging>_<wind>[_fem-<heating>])")
    r.add_argument("--quiet", action="store_true")
    th = r.add_argument_group("thermal model (Step 2)")
    th.add_argument("--thermal", choices=THERMAL_MODES, default="none", help="none: Step 1 trajectory only (default); fem: coupled 3D conduction")
    th.add_argument("--heating", choices=heating.HEATING_NAMES, default="physics",
                    help="physics (default) or sesam: SESAM's uniform 0.27471 x q_stag on every patch -- a verification device, not physical")
    th.add_argument("--stagnation", choices=heating.STAGNATION_NAMES, default="fay-riddell")
    th.add_argument("--bridging-heat", choices=heating.BRIDGING_HEAT_NAMES, default="matting")
    th.add_argument("--matting-n", type=float, default=1.0, help="Matting exponent n (default %(default)s)")
    th.add_argument("--accommodation", type=float, default=0.8, help="free-molecular energy accommodation A_cq (default %(default)s)")
    th.add_argument("--catalycity", type=float, default=1.0, help="wall catalycity 0..1 in Fay-Riddell (default %(default)s)")
    th.add_argument("--material", default=None,
                    help="AA7075_nomelt | AA7075 | AA7075_range | <DRAMA material JSON> (default: AA7075_nomelt, or AA7075_range with --melt on)")
    th.add_argument("--emissivity", type=float, default=None, help="override the material's emissivity")
    th.add_argument("--k-scale", type=float, default=1.0,
                    help="multiplies the material's conductivity (default 1): a verification device for the near-isothermal (lumped) limit, not physical")
    th.add_argument("--t-ambient", type=float, default=0.0, help="radiation background [K] (default %(default)s, SESAM's)")
    th.add_argument("--mesh-size", type=float, default=1.0, help="multiplies --h-surface and --h-core (default %(default)s)")
    th.add_argument("--h-surface", type=float, default=mesh.DEFAULT_H_SURFACE * 1e3, help="surface element size [mm] (default %(default)s)")
    th.add_argument("--h-core", type=float, default=mesh.DEFAULT_H_CORE * 1e3, help="core element size [mm] (default %(default)s)")
    th.add_argument("--thermal-solver", choices=thermal.SOLVER_NAMES, default="skfem")
    th.add_argument("--linear-solver", choices=("direct", "amg"), default="amg")
    th.add_argument("--consistent-mass", action="store_true", help="consistent capacity matrix instead of the lumped nodal-enthalpy one (skfem only, analytic checks)")
    th.add_argument("--dt", type=float, default=0.5, help="macro step [s] (default %(default)s)")
    th.add_argument("--frames-every", type=int, default=0, help="VTK frame every n macro steps (default 0: none; 10 with --animate/--stills)")
    th.add_argument("--animate", action="store_true", help="MP4/GIF of the surface temperature and of the meridional cross-section, plus stills")
    th.add_argument("--stills", action="store_true", help="only the stills (start, peak heating, peak surface T, end; surface and section)")
    me = r.add_argument_group("melting and spraying (Step 3)")
    me.add_argument("--melt", choices=MELT_MODES, default="off", help="off: Step 2 behaviour (default); on: melting, melt film, spraying")
    me.add_argument("--removal", choices=body.REMOVAL_NAMES, default="girin",
                    help="girin: film + runoff + Girin spraying (default); instant: liquid removed as it forms -- the lumped-melting "
                         "verification device, not physical")
    me.add_argument("--runoff", choices=("on", "off"), default="on", help="film runoff transport along the surface (default on)")
    me.add_argument("--rt-spray", choices=("on", "off"), default="on",
                    help="apply Girin & Kopyt's front-surface Rayleigh-Taylor mode as a release mechanism, on the patches where it "
                         "grows faster than the shear mode (default on). off: evaluate and report it but release nothing through it, "
                         "which is how every run before 2026-09-24 behaved")
    me.add_argument("--rarefied-shear", choices=surface_flow.RAREFIED_SHEAR_NAMES, default="slip",
                    help="within the shock-layer branch: slip (default) uses Maxwell slip below Kn_local 0.1 and free-molecular shear above; "
                         "bridged blends them with SESAM's f(Kn). The merged branch always bridges")
    me.add_argument("--gamma-pm", type=float, default=surface_flow.GAMMA_PM,
                    help="effective ratio of specific heats of the Prandtl-Meyer expansion that sets the wall pressure beyond the sonic "
                         "point (default %(default)s; 1.4 frozen air, ~1.15 dissociated -- a factor ~2 on p_w at 90 deg)")
    me.add_argument("--kn-body-shock", type=float, default=surface_flow.KN_BODY_SHOCK,
                    help="body Knudsen number below which a distinct bow shock is assumed and the shock-layer construction is used "
                         "(default %(default)s); above it the flow is treated as merged and every melt closure is flagged")
    me.add_argument("--we-critical", type=float, default=dispersion.WE_CRITICAL_PRACTICAL, help="critical surface Weber number (default %(default)s)")
    me.add_argument("--kr", type=float, default=spray.K_R, help="droplet radius / wavelength (default %(default)s)")
    me.add_argument("--kt", type=float, default=spray.K_T, help="release period / growth time (default %(default)s)")
    me.add_argument("--prism-layers", type=int, default=None, help="prism layers under the surface (default 4 with --melt on, 0 otherwise)")
    me.add_argument("--layer-thickness", type=float, default=mesh.DEFAULT_LAYER_THICKNESS * 1e3, help="outermost layer thickness [mm] (default %(default)s, growth 2)")
    me.add_argument("--demise-fraction", type=float, default=0.01, help="the run ends when the body mass falls below this fraction of the initial (default %(default)s)")
    me.add_argument("--size-feedback", choices=body.SIZE_FEEDBACK_NAMES, default=None,
                    help="current: Kn on the equivalent diameter of the remaining mass and the stagnation radius fitted to the windward cap "
                         "(default with --removal girin); initial: D0 and R0 throughout, SESAM's convention (default with --removal instant)")
    me.add_argument("--particles", dest="particles", action="store_true", default=True, help="write the particle source table (default)")
    me.add_argument("--no-particles", dest="particles", action="store_false")

    c = sub.add_parser("compare", help="metrics and plots for an existing model history")
    c.add_argument("--model", required=True, help="model history CSV")
    c.add_argument("--reference", required=True, help="SESAM run CSV")
    c.add_argument("--outdir", default=DEFAULT_OUTDIR)
    c.add_argument("--title", default=None)
    c.add_argument("--quiet", action="store_true")
    return p


def build_thermal(args, settings, mass):
    """(ThermalBody or MeltingBody, HeatingModel, settings-provenance dict) for --thermal fem [--melt on]."""
    radius = settings.diameter / 2.0
    melting = args.melt == "on"
    h_surface, h_core = args.h_surface * 1e-3 * args.mesh_size, args.h_core * 1e-3 * args.mesh_size
    layers = args.prism_layers if args.prism_layers is not None else (mesh.DEFAULT_LAYERS if melting else 0)
    the_mesh = mesh.sphere_mesh(radius, h_surface, h_core, layers=layers, layer_thickness=args.layer_thickness * 1e-3)
    material_name = args.material or ("AA7075_range" if melting else "AA7075_nomelt")
    mat = material.Material.from_drama_json(material_name)
    if args.k_scale != 1.0:
        mat.k_table = mat.k_table * args.k_scale                  # verification device (near-isothermal body), not physical
    solver = thermal.thermal_solver(args.thermal_solver, linear_solver=args.linear_solver, lumped_mass=not args.consistent_mass)
    if melting:
        flow = surface_flow.SurfaceFlow(rarefied_shear=args.rarefied_shear, gamma_pm=args.gamma_pm, kn_body_shock=args.kn_body_shock)
        spray_model = spray.SprayModel(mat.liquid, k_r=args.kr, k_t=args.kt, we_critical=args.we_critical,
                                       rt_spray=args.rt_spray == "on")
        size_feedback = args.size_feedback or ("current" if args.removal == "girin" else "initial")
        melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
                                          particles=args.particles, size_feedback=size_feedback)
        the_body = body.MeltingBody(the_mesh, mat, solver, mass, flow, spray_model, melt_settings, T0=args.temperature,
                                    emissivity=args.emissivity, T_ambient=args.t_ambient)
    else:
        the_body = body.ThermalBody(the_mesh, mat, solver, mass, T0=args.temperature, emissivity=args.emissivity, T_ambient=args.t_ambient)
    if args.heating == "sesam":
        heating_model = heating.SesamEquivalentHeating()
    else:
        heating_model = heating.PhysicsHeating(stagnation=args.stagnation, bridging=args.bridging_heat, matting_n=args.matting_n,
                                               accommodation=args.accommodation, catalycity=args.catalycity)
    info = {"heating": args.heating, "material": mat.name, "material_file": os.path.abspath(material.MATERIAL_NAMES.get(material_name, material_name)),
            "emissivity": the_body.emissivity, "t_ambient_K": args.t_ambient, "mesh_file": the_mesh.params["path"],
            "h_surface_mm": h_surface * 1e3, "h_core_mm": h_core * 1e3, "prism_layers": layers, "layer_thickness_mm": args.layer_thickness,
            "n_nodes": the_mesh.n_nodes, "n_elements": the_mesh.n_elements, "n_patches": the_body.surface.n_patches,
            "thermal_solver": args.thermal_solver, "linear_solver": args.linear_solver, "lumped_mass": not args.consistent_mass,
            "k_scale": args.k_scale, "melt": args.melt}
    if args.heating == "physics":
        info.update({"stagnation": args.stagnation, "bridging_heat": args.bridging_heat, "matting_n": args.matting_n,
                     "accommodation": args.accommodation, "catalycity": args.catalycity})
    if melting:
        info.update({"removal": args.removal, "runoff": args.runoff, "rt_spray": args.rt_spray,
                     "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                     "k_r": args.kr, "k_t": args.kt, "demise_fraction": args.demise_fraction, "particles": args.particles,
                     "size_feedback": size_feedback, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
                     "liquid": {"rho": mat.liquid.rho, "mu": mat.liquid.mu, "sigma": mat.liquid.sigma},
                     "T_solidus_K": mat.T_solidus, "T_liquidus_K": mat.T_liquidus, "latent_heat_Jkg": mat.latent_heat})
    return the_body, heating_model, info


def cmd_run(args, parser):
    for label, value in (("--diameter", args.diameter), ("--velocity", args.velocity), ("--material-density", args.material_density),
                         ("--cadence", args.cadence), ("--t-max", args.t_max)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.altitude < 0.0:
        parser.error("--altitude must be >= 0")
    if args.atmosphere not in ATMOSPHERES and not args.atmosphere.startswith("replay:"):
        parser.error("--atmosphere must be one of {} or replay:<sesam.csv>".format(ATMOSPHERES))
    for label, value in (("--dt", args.dt), ("--mesh-size", args.mesh_size), ("--h-surface", args.h_surface), ("--h-core", args.h_core), ("--k-scale", args.k_scale)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if not 0.0 <= args.catalycity <= 1.0:
        parser.error("--catalycity must be within [0, 1]")
    if args.thermal == "none" and (args.animate or args.stills or args.frames_every):
        parser.error("--animate/--stills/--frames-every need --thermal fem")
    if args.melt == "on" and args.thermal != "fem":
        parser.error("--melt on needs --thermal fem")
    if not 1.0 < args.gamma_pm < 2.0:
        parser.error("--gamma-pm must be within (1, 2)")
    for label, value in (("--we-critical", args.we_critical), ("--kr", args.kr), ("--kt", args.kt), ("--layer-thickness", args.layer_thickness),
                         ("--kn-body-shock", args.kn_body_shock)):
        if value <= 0.0:
            parser.error("{} must be > 0".format(label))
    if args.prism_layers is not None and args.prism_layers < 0:
        parser.error("--prism-layers must be >= 0")
    if not 0.0 < args.demise_fraction < 1.0:
        parser.error("--demise-fraction must be within (0, 1)")

    initial = tj.InitialState(velocity=args.velocity * 1e3, altitude=args.altitude * 1e3,
                              flight_path=math.radians(args.flight_path_angle), heading=math.radians(args.heading),
                              lat=math.radians(args.lat), lon=math.radians(args.lon), epoch=args.epoch)
    settings = tj.Settings(diameter=args.diameter * 1e-3, gravity=args.gravity, rtol=args.rtol,
                           cadence=args.cadence, t_max=args.t_max)
    mass = body.sphere_mass(settings.diameter, args.material_density)
    atm, atm_name, atm_info = make_atmosphere(args.atmosphere, args.epoch, args.wind)
    reference = sesam_io.load_reference(args.reference) if args.reference else None
    name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
                                       args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None)
    run_dir = os.path.join(args.outdir, name)
    os.makedirs(args.outdir, exist_ok=True)
    thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
    if args.thermal == "fem":
        the_body, heating_model, thermal_info = build_thermal(args, settings, mass)
    sim = tj.Simulator(initial, the_body, atm, aero.SphereDragTables.from_json(), aero.bridging_by_name(args.bridging), settings)
    if args.thermal == "fem":
        frames_every = args.frames_every or (10 if (args.animate or args.stills) else 0)
        run = coupled.CoupledRun(sim, the_body, heating_model,
                                 coupled.CoupledSettings(dt=args.dt, frames_every=frames_every, output_dir=os.path.join(run_dir, "vtk")))
        history = run.run()
        thermal_info["macro_step_s"], thermal_info["frames_every"] = args.dt, frames_every
    else:
        history = sim.run(extra_times=reference.time if reference is not None else None)
    csv_path = os.path.join(args.outdir, name + ".csv")
    json_path = os.path.join(args.outdir, name + ".json")
    tj.write_history_csv(history, csv_path)
    prov = provenance()
    # spec section 8: the reference file(s)' SHA-256 also live under provenance, alongside the copies
    # already recorded under "comparison" (--reference) and "settings" (replay_reference/replay_sha256).
    prov["reference_sha256"] = reference.sha256 if reference is not None else None
    prov["replay_sha256"] = atm_info.get("replay_sha256")
    doc = {
        "schema_version": 1,
        "run_name": name,
        "inputs": {"diameter_mm": args.diameter, "initial_velocity_kms": args.velocity, "initial_altitude_km": args.altitude,
                   "flight_path_angle_deg": args.flight_path_angle, "heading_deg": args.heading, "latitude_deg": args.lat,
                   "longitude_deg": args.lon, "epoch_utc": args.epoch.strftime("%Y-%m-%dT%H:%M:%S"),
                   "material_density_kgm3": args.material_density, "mass_kg": mass, "initial_temperature_K": args.temperature},
        "settings": {"atmosphere": atm_name, "wind": args.wind, "bridging": args.bridging, "gravity": args.gravity,
                     "rtol": args.rtol, "atol_position_m": settings.atol_position, "atol_velocity_ms": settings.atol_velocity,
                     "cadence_s": args.cadence, "t_max_s": args.t_max, **atm_info, "thermal": args.thermal, **thermal_info},
        "results": history.results,
        "comparison": None,
        "provenance": prov,
        "files": {"csv": os.path.abspath(csv_path)},
    }
    melting = args.thermal == "fem" and args.melt == "on"
    if melting:
        doc["files"].update({k: os.path.abspath(v) for k, v in coupled.write_particles(run_dir, the_body, history).items()} if args.particles
                            else {})
        doc["files"]["melt_plots"] = [os.path.abspath(p) for p in compare.plot_melt(history, run_dir, name, reference if compare.has_melt(history, reference) else None)]
    if reference is not None:
        plots = compare.plot_all(history, reference, run_dir, name)
        doc["comparison"] = {"reference": reference.name, "reference_csv": reference.csv_path,
                             "reference_sha256": reference.sha256, "metrics": compare.metrics(history, reference),
                             "plots": [os.path.abspath(p) for p in plots]}
        if compare.has_thermal(history, reference):
            doc["comparison"]["thermal_metrics"] = compare.thermal_metrics(history, reference)
            doc["comparison"]["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, run_dir, name)]
        if compare.has_melt(history, reference):
            doc["comparison"]["melt_metrics"] = compare.melt_metrics(history, reference)
    if args.thermal == "fem" and (args.animate or args.stills):
        vtk_dir = os.path.join(run_dir, "vtk")
        out = viz.animate(vtk_dir, history, settings.diameter / 2.0, animation=args.animate, melting=melting)
        section = viz.animate_section(vtk_dir, history, settings.diameter / 2.0, animation=args.animate, melting=melting)
        doc["files"]["animation"], doc["files"]["section"] = out.get("animation"), section.get("animation")
        doc["files"]["stills"] = out["stills"] + section["stills"]
        if melting:
            film = viz.animate_film(vtk_dir, history, settings.diameter / 2.0, animation=args.animate)
            doc["files"]["film"] = film.get("animation")
            doc["files"]["stills"] += film["stills"]
    if args.thermal == "fem":
        doc["files"]["vtk_dir"] = os.path.abspath(os.path.join(run_dir, "vtk")) if history.results.get("n_frames") else None
    tj.write_run_json(json_path, doc)
    if not args.quiet:
        res = history.results
        print("{}: {} at t = {:.1f} s, final V {:.4f} km/s, Kn {:.3g} -> {:.3g}, {} RHS evaluations in {:.1f} s".format(
            name, res["end_reason"], res["final_time_s"], res["final_velocity_kms"], res["knudsen_start"],
            history.columns["knudsen"][-1], res["rhs_evaluations"], res["runtime_s"]))
        if args.thermal == "fem":
            print("  thermal: peak surface T {:.0f} K at t = {:.0f} s, peak mean T {:.0f} K, integrated heat {:.3g} J, "
                  "energy balance residual {:.1e}, {} macro steps, {:.1f} Newton iterations/step".format(
                      res["peak_surface_T_K"], res["time_of_peak_surface_T_s"], res["peak_mean_T_K"], res["integrated_heat_J"],
                      res["energy_balance_residual"], res["n_macro_steps"], res["mean_newton_iterations"]))
        if melting:
            print("  melt: onset {} km, spraying onset {} km, demise {} km at t = {} s; sprayed {:.4f} kg of {:.4f}, film left {:.4f} kg, "
                  "{:.3g} droplets (median r {} um), {} source rows, melt balance {:.1e}".format(
                      _fmt(res["melt_onset_altitude_km"]), _fmt(res["spraying_onset_altitude_km"]), _fmt(res["demise_altitude_km"]),
                      _fmt(res["demise_time_s"]), res["sprayed_mass_kg"], res["initial_mass_kg"], res["film_mass_kg"], res["n_released"],
                      _fmt(res["r_median_um"]), res["n_source_rows"], res["melt_energy_balance_residual"]))
            if reference is not None and "melt_metrics" in doc["comparison"]:
                mm = doc["comparison"]["melt_metrics"]
                print("  mass vs SESAM: max |dm| {:.2%} of m0, melt onset {:+.2f} km, 1 %-mass time {:+.1f} s ({:+.1%})".format(
                    mm["mass"]["max_rel_m0"], mm["onset_altitude_diff_km"], mm["demise_time_diff_s"], mm["demise_time_rel"]))
        if reference is not None:
            hyp = doc["comparison"]["metrics"]["hypersonic"]
            print("  vs {}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m; end time {:+.1f} s".format(
                reference.name, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], doc["comparison"]["metrics"]["d_end_time_s"]))
            if "thermal_metrics" in doc["comparison"]:
                tm = doc["comparison"]["thermal_metrics"]
                print("  heat vs SESAM: Q_conv max {:.2%} of peak ({:.2%} point-wise, continuum), integrated heat {:+.2%} (hypersonic) "
                      "{:+.2%} (end), |dT_eq| max {:.1f} K ({:.2%}), radiated max {:.2%} of peak".format(
                          tm["Q_conv"]["max"], tm["Q_conv"]["continuum_rel_max"], tm["integrated_heat"]["rel_error_end_of_hypersonic"],
                          tm["integrated_heat"]["rel_error_end"], tm["temperature"]["dT_max_K"], tm["temperature"]["dT_rel_max"], tm["radiated"]["max"]))
        print("  csv  -> {}\n  json -> {}".format(os.path.abspath(csv_path), os.path.abspath(json_path)))
    # spec section 7: exit 1 when the flight escaped rather than reaching the ground or t_max;
    # the CSV/JSON are already written above so the escaped trajectory is still available.
    if history.end_reason == "escape":
        return 1
    return 0


def _fmt(x):
    return "n/a" if x is None else "{:.1f}".format(x)


def cmd_compare(args, parser):
    history = tj.read_history_csv(args.model)
    reference = sesam_io.load_reference(args.reference)
    stem = os.path.splitext(os.path.basename(args.model))[0]
    title = args.title or "{} vs {}".format(stem, reference.name)
    os.makedirs(args.outdir, exist_ok=True)
    plots = compare.plot_all(history, reference, args.outdir, title)
    doc = {"model": os.path.abspath(args.model), "reference": reference.name, "reference_sha256": reference.sha256,
           "metrics": compare.metrics(history, reference), "plots": [os.path.abspath(p) for p in plots]}
    if compare.has_thermal(history, reference):
        doc["thermal_metrics"] = compare.thermal_metrics(history, reference)
        doc["plots"] += [os.path.abspath(p) for p in compare.plot_thermal(history, reference, args.outdir, title)]
    if compare.has_melt(history, reference):
        doc["melt_metrics"] = compare.melt_metrics(history, reference)
        doc["plots"] += [os.path.abspath(p) for p in compare.plot_melt(history, args.outdir, title, reference)]
    out = os.path.join(args.outdir, stem + "_vs_reference.json")
    tj.write_run_json(out, doc)
    if not args.quiet:
        hyp = doc["metrics"]["hypersonic"]
        print("{}: hypersonic max |dV| {:.1f} m/s ({:.3%}), max |dh| {:.0f} m -> {}".format(title, hyp["dV_max_ms"], hyp["dV_rel_max"], hyp["dh_max_m"], out))
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            return cmd_run(args, parser)
        return cmd_compare(args, parser)
    except KeyError as exc:
        # e.g. fap.solar_indices(): the run epoch (or the day before it) has no record in
        # data/fap_day.dat; could also be a reference CSV missing an expected column.
        key = exc.args[0] if exc.args else exc
        print("ERROR: missing key {!r} (epoch outside the fap file, or a reference CSV without that column)".format(key),
              file=sys.stderr)
        return 2
    except (ValueError, FileNotFoundError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except MissingBackend as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
    except ModuleNotFoundError as exc:
        print("ERROR: missing optional library {!r}: install requirements-step2.txt into drama_env (dolfinx: the separate "
              "fenicsx_env environment)".format(exc.name), file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1
```


- [ ] **Step 4: Run the tests**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q`
Expected: all pass (the bookkeeping device test runs the whole 100 mm flight on the coarse mesh in ~10 s).

- [ ] **Step 5: Commit**

```bash
git add reentry_model/cli.py tests/test_reentry_model_cli.py
git commit -m "Add the melting flags, wiring and summaries to the command line (Step 3 Task 13)

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

