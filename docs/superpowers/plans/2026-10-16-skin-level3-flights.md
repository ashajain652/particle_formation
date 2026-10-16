# Level 3 of the layered skins: three Step 2 flights on the Mac Studio — runbook

> **For the agent running this.** You are on Asha's Mac Studio. Your job is to run three long, already-written
> verification flights, watch them, gather their results, and report. **You write no model code.** If anything in
> this runbook does not work as written, stop and tell Asha what happened rather than changing code or settings to make
> it work. Read the whole runbook once before running anything.

**Asked by:** Asha, 2026-10-16 ("go with option B. I want to run the other three flights on my macstudio with another
agent there").

**Design of record:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` §11.1 (level 3) and §11.2
(build stage 2 ends with level 3); sub-plan 20, Task 10 Step 4, in
`docs/superpowers/plans/melt-spraying-subplans/20-layered-skin-heat.md`.

## 1. What this is and why it matters

The repository models what happens to the solid aluminium (AA7075) fragments of a re-entering satellite. Step 2 of the
model is a 3D finite-element heat model of a sphere, verified against ESA's SESAM code on two reference flights that do
not melt: a **100 mm** sphere entering at 77.5 km and a **50 mm** sphere entering at 115 km, both at 7.5 km/s, US76
atmosphere, no wind. Each flight is run with two heating models:

- **`sesam`** — "SESAM-equivalent" heating: SESAM's surface-averaged heat input spread evenly over the surface. This is a
  verification device, not physics; it lets the model's temperature be compared with SESAM's single lumped temperature.
  These runs are **judged** against Step 2's thresholds.
- **`physics`** — the model proper (stagnation-point heating with a distribution around the sphere). These runs are
  **reported**, with only a loose sanity bound on the total heat.

The new part of the model, **layered skins** (sub-plan 20), puts a column of very thin cells (10 micrometres at the
top) on every corner of the surface mesh, above a 3D mesh that can move as the surface recedes. Before flow in the
molten layer is built on top of the skins (sub-plan 21), the spec requires **level 3**: the same Step 2 flights with
`--surface-model skin`, against the same thresholds, and a measured cost for a whole flight. That is what you run. The
flights do not melt (the material is `AA7075_nomelt`), so the skins only conduct heat; nothing is removed.

Four skin flights make up level 3. Asha's laptop is running one of them (**100 mm, `sesam`**, started 2026-10-16 at
18:07 EDT). You run the other three:

| Label | Sphere | Heating | Judged? | Macro steps (0.5 s each) | Expected wall time (estimate) |
|---|---|---|---|---|---|
| `d100__physics__skin` | 100 mm from 77.5 km | physics | reported | 732 (ground at about 366 s) | 12–18 h |
| `d050__sesam__skin` | 50 mm from 115 km | sesam | **judged** | about 1 118 (ground at about 559 s) | 5–8 h |
| `d050__physics__skin` | 50 mm from 115 km | physics | reported | about 1 118 | 5–8 h |

The wall times are extrapolations, not measurements: on the laptop (Apple M3) a 0.5 s step of the 50 mm skin flight
cost 15 s at the start of the flight, and a step of the 100 mm sphere below melting cost 60–90 s when the skins were
timed on 2026-10-12. Each flight uses about one CPU core (measured: 56 s of CPU in 66 s of wall time), so the three
can run side by side. Peak memory measured: 1.1 GB for the 50 mm smoke check; the 100 mm flight held 2.0 GB a few
minutes in. You also run the four **element-model** flights (the Step 2 model without skins), which take minutes, so
that the skin results can be put beside the element results from the same machine.

## 2. Rules

- Do not edit anything under `prototype/`, `reentry_model/`, `analysis/`, `tests/` or `docs/`, except creating the
  results directory of Step 9. The code is verified and its plan is generated from it byte for byte; a hand edit breaks
  that. If you think code is wrong, report it.
- Do not change any flight setting (macro step, mesh, cells, heating) and do not re-run a flight with different
  settings to make it pass. A `FAIL` verdict is a result to report, not an error to fix.
- Do not touch DRAMA, pyDRAMA, `sphere_reentry.py` or `sphere_sweep.py`; nothing here needs them.
- Git: never use bare `git stash`. Push only the results branch named in Step 9, never `main` or
  `step3-rigid-substrate-scheil`. Do not merge anything.
- Quote every path: the repository's path on Asha's machines contains spaces (for example `University Documents /MIT`).
- When in doubt, stop and ask Asha. Report in her style: full sentences, every number with what was measured, its units
  and why it matters, terms spelled out (her `~/.claude/CLAUDE.md` gives the full rules).

## 3. Get the code

- [ ] **Step 1.** In the Mac Studio's clone of the repository (`git@github.com:ashajain652/particle_formation.git`):

```bash
cd "<the clone's root>"
git status --short          # must show no modified tracked files; if it does, stop and ask Asha
git fetch origin
git switch step3-rigid-substrate-scheil
git pull --ff-only
```

Then set the two paths every later step uses, and confirm the level-3 option is present:

```bash
REPO="$(git rev-parse --show-toplevel)"
CODE="$REPO/prototype/work-2026-10-12-skin-heat-plan/code"
grep -c "summary-only" "$CODE/analysis/reentry_model_thermal_verification.py"     # expect a number above 0
```

`$CODE` is a complete, runnable copy of the model package with the skins: the tracked copy from which sub-plan 20 is
generated, byte-identical to the laptop's build repository (`prototype/work-2026-10-12-skin-heat/code` at f6e12fb,
which is git-ignored and so not in your clone). Every command below runs from `$CODE`, never from `$REPO`.

## 4. The Python environment

- [ ] **Step 2.** The laptop uses the conda environment `drama_env` at `/Users/ashajain/miniforge3/envs/drama_env`.
Find the Mac Studio's:

```bash
conda env list
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python      # adjust if the Mac Studio's path differs
"$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera, matplotlib, meshio, pytest; print(sys.version.split()[0], numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__, matplotlib.__version__)"
```

On the laptop this prints Python 3.12.14, numpy 2.5.3, scipy 1.18.1, scikit-fem 12.0.2, gmsh 4.15.2, pyamg 5.3.0,
cantera 3.2.0, matplotlib 3.11.1. Small differences in numpy, scipy or matplotlib versions are acceptable; a different
**gmsh** or **scikit-fem** version is worth noting in the report, because the meshes and the assembly come from them.
If the environment is missing or an import fails, create or repair it (this is the only environment work allowed):

```bash
conda create -n drama_env python=3.12 numpy scipy matplotlib pytest -c conda-forge
"$PY" -m pip install -r "$REPO/requirements-step2.txt"      # scikit-fem, gmsh, meshio, pyamg, pyvista, imageio(-ffmpeg), cantera
```

pyDRAMA, pymsis and FEniCSx are not needed.

## 5. Check the machine before spending a day on it

- [ ] **Step 3: the unit tier.**

```bash
cd "$CODE"
"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py \
    --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py
```

Expected: **403 passed, 1 skipped** (the skip is the FEniCSx module, because `dolfinx` is not in `drama_env`). The
three ignored modules need the DRAMA wrapper scripts at the repository root and are out of scope. Any failure: stop and
report it with its output.

- [ ] **Step 4: the smoke checks.** Two flights cut off at 2 s of flight time, one with skins and one without, into a
scratch directory. They also build and cache the 50 mm meshes, so that the two 50 mm skin flights of Step 5 do not both
try to build the same mesh file at once.

```bash
SMOKE="$CODE/reentry_model_output/level3_smoke"
/usr/bin/time -l "$PY" analysis/reentry_model_thermal_verification.py --outdir "$SMOKE" --cases d050 --modes sesam \
    --surface-model skin --t-max 2
"$PY" analysis/reentry_model_thermal_verification.py --outdir "$SMOKE" --cases d050 --modes sesam --t-max 2
"$PY" -c "import json,sys; r=json.load(open(sys.argv[1]))['results']; print(r['end_reason'], r['n_macro_steps'], r['skin_columns'], r['energy_balance_residual'], r['runtime_s'])" "$SMOKE/d050__sesam__skin__t2.json"
```

Expected, as measured on the laptop: the skin run ends by `t_max` after **4** macro steps on **2 367** skin columns, with
an energy-balance residual of about 3e-12 (it must be below 1e-9) and a run time of about a minute; `/usr/bin/time -l`
reports a "maximum resident set size" of about 1.1e9 bytes. Each run writes `<label>.json`, `<label>.csv`,
`<label>.progress` and a plot directory. A different number of columns means a different mesh: stop and report.

## 6. Run the three skin flights

- [ ] **Step 5: launch.** First check the machine has room: `sysctl -n hw.perflevel0.physicalcpu hw.memsize` (the three
flights need three free cores and about 10 GB of memory between them; if there are fewer than four performance cores or
less than 16 GB, run the 50 mm flights one after the other instead, and tell Asha). Nothing else heavy should be
running. Then:

```bash
cd "$CODE"
OUT="$CODE/reentry_model_output/verification_thermal"
LOGS="$CODE/reentry_model_output/level3_logs"
mkdir -p "$OUT" "$LOGS"
launch() {   # case, mode
    nohup caffeinate -i /usr/bin/time -l "$PY" analysis/reentry_model_thermal_verification.py --outdir "$OUT" \
        --cases "$1" --modes "$2" --surface-model skin > "$LOGS/$1_$2.log" 2>&1 &
    echo "$1 $2 pid $!  started $(date '+%F %T')" | tee -a "$LOGS/launches.txt"
}
launch d100 physics
launch d050 sesam
launch d050 physics
```

`caffeinate -i` keeps the Mac from idle sleep for as long as each flight runs; `/usr/bin/time -l` writes the flight's
wall time and peak memory at the end of its log. The flights must outlive the shell command that started them: if your
harness ends background jobs when a command returns, start each `launch` line with its own background-run facility
(in Claude Code, `run_in_background`) instead. Either way, a minute later confirm with
`ps -Ao pid,etime,pcpu,command | grep reentry_model_thermal_verification | grep -v grep` that three flights are
running, each near 100 % of a core. The 100 mm flight builds its own mesh (no other process here needs
that one); the 50 mm mesh was cached by Step 4. There is no checkpoint: a flight interrupted by a reboot or a crash
starts again from the beginning.

- [ ] **Step 6: the element-model flights**, while the skin flights run (they take one more core for a few minutes):

```bash
cd "$CODE"
"$PY" analysis/reentry_model_thermal_verification.py --outdir "$OUT" > "$LOGS/elements.log" 2>&1
cat "$OUT/summary.md"
```

This runs all four element flights one after another (`d100__sesam`, `d100__physics`, `d050__sesam`,
`d050__physics`) and writes `summary.md`. They should reproduce the README's Step 2 verification table (README.md,
"Verification" under Step 2), measured on the laptop:

| case | heating | Q_conv max (of peak) | Q_conv point-wise (continuum) | integrated heat (end of hypersonic / end) | max ΔT_eq | radiated (of peak) | peak surface T |
|---|---|---|---|---|---|---|---|
| d100 | sesam | 0.46 % | 2.24 % | +0.31 % / +0.30 % | 24.5 K (1.20 %) | 6.65 % | 2089 K at 139 s |
| d050 | sesam | 0.37 % | 2.36 % | +0.03 % / +0.02 % | 28.0 K (1.22 %) | 3.40 % | 2442 K at 245 s |
| d100 | physics | — | — | ratio 0.744 | — | — | 2319 K at 118 s |
| d050 | physics | — | — | ratio 0.770 | — | — | 2493 K at 239 s |

Agreement to the printed digits, or within one unit of the last digit, means the Mac Studio's environment reproduces
Step 2. A larger difference does not stop the skin flights, but report it: it would mean the two machines' results
cannot be compared digit for digit.

## 7. Watch them

- [ ] **Step 7.** Every flight appends a line to `$OUT/<label>.progress` every 20 macro steps (10 s of flight time):

```
step    120  t    60.00 s  h   61.234 km  V  7.4800 km/s  T   512.3 K  T_surface_max   640.1 K  wall      2400 s
```

`T` is the body's mean temperature and `T_surface_max` the hottest surface point. Check about once an hour:

```bash
for f in "$OUT"/*__skin.progress; do echo "$(basename "$f"): $(tail -n 1 "$f")"; done
ps -Ao pid,etime,pcpu,rss,command | grep reentry_model_thermal_verification | grep -v grep
```

From the last line, estimate each flight's finish: (steps still to go) × (wall seconds per step over the last few
lines), with 732 steps for the 100 mm sphere and about 1 118 for the 50 mm. The step cost is expected to rise through
peak heating (around 120–140 s of flight for the 100 mm sphere, 240–250 s for the 50 mm) and fall after it. Each
process should sit near 100 % of one core. Give Asha a one-line estimate after the first hour.

If a flight's process has gone, read the end of its log (`tail -n 60 "$LOGS/<case>_<mode>.log"`):

- **"failed with exit code 1"**, an `ERROR:` line or a Python traceback: a model failure (for example a Newton solve
  that did not converge, or a refused mesh). Do not retry and do not change settings. Report the error text, the last
  progress line (how far the flight got) and the log tail.
- **"failed with exit code 2"**: a bad argument or a missing library; fix the environment only (Step 2), then
  re-launch that flight.
- **no error, but the process vanished** (a reboot, or the machine slept): re-launch that flight from the beginning, and
  report that it was restarted and why.

A progress file that has not grown for more than two hours while its process is still at about 100 % of a core is not a
failure; note it with the time, and keep watching.

## 8. When all three have finished

- [ ] **Step 8: the summary.**

```bash
cd "$CODE"
"$PY" analysis/reentry_model_thermal_verification.py --outdir "$OUT" --surface-model skin --summary-only
cat "$OUT/summary_skin.md"
for f in "$LOGS"/d*.log; do echo "== $f"; grep -E "real|maximum resident set size|exit code|ERROR" "$f"; done
```

`summary_skin.md` has a row per flight, each skin flight below the element flight of the same case and heating:

- **verdict** — `pass` or `FAIL (<the metrics outside their threshold>)` for `sesam` heating, against Step 2's
  thresholds: convective power within 3 % of SESAM's peak at every point of the hypersonic phase (speed above 1 km/s)
  and within 3 % point by point where the flow is continuum; heat taken in by the end of the hypersonic phase within
  3 %; the body's energy-equivalent temperature within 2 % of SESAM's lumped temperature over the whole flight;
  radiated power within 8 % of SESAM's peak; the flight reaching the ground; the energy books closed to 1e-6. For
  `physics` heating it is `reported` if the heat taken in by the end is between 0.4 and 1.1 times SESAM's, which is the
  only check the Step 2 test makes of that model.
- **skins** — the number of skin columns, mesh rebuilds (expected 0: nothing recedes), interface repeats (extra passes
  of the skin–3D coupling), base faults (expected 0), the wall time spent in the skins and in the 3D solve, and the
  wall time a macro step.

## 9. Bring the results back

- [ ] **Step 9.** Copy the small result files into a tracked directory, record the machine, and commit on a new
branch:

```bash
cd "$REPO"
R="$REPO/prototype/work-2026-10-12-skin-heat-plan/level3/macstudio"
mkdir -p "$R"
cp "$OUT"/d*.json "$OUT"/d*.csv "$OUT"/d*.progress "$OUT"/summary*.md "$OUT"/summary*.json "$R"/
cp "$LOGS"/*.log "$LOGS"/launches.txt "$R"/
{ sysctl -n machdep.cpu.brand_string hw.perflevel0.physicalcpu hw.ncpu hw.memsize; sw_vers;
  "$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera; print(sys.version); print(numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__)";
  git rev-parse HEAD; } > "$R/environment.txt"
du -sh "$R"            # expect a few MB; PNG plots and VTK files are not copied (and are git-ignored)
git switch -c skin-level3-macstudio
git add "$R"
git commit -m "Level 3 of the layered skins on the Mac Studio: the 100 mm physics and both 50 mm flights with a skin on every corner, and the four element flights beside them"
git push -u origin skin-level3-macstudio
```

Push only this branch. Asha (or the laptop's agent) gathers these JSONs with the laptop's `d100__sesam__skin.json` into
one directory and runs `--summary-only` there for the combined table.

## 10. Report to Asha

- [ ] **Step 10.** Write the report in Asha's style (full sentences; every number with what it measures, its units, and
why it matters; terms spelled out; no symbols in prose). In this order:

1. Whether the three skin flights and the four element flights ran to the ground, and whether the results are pushed
   (branch and commit).
2. For `d050__sesam__skin`: the verdict, and each of the five measured values next to its threshold and next to the
   element flight's value from the same machine.
3. For the two `physics` flights: the ratio of heat taken in to SESAM's at the end, the peak surface temperature and
   when it occurred, each beside the element flight's.
4. The cost of each skin flight: total wall time, wall time per macro step (and how it changed through the flight, from
   the progress files), the split between the skins and the 3D solve, peak memory, and the element flight's wall time
   for comparison. Say that the three ran side by side on one machine, so their wall times are not those of a flight
   run alone.
5. Anything unexpected: a rebuild, base faults, a stalled progress file, a restart, an environment difference from the
   laptop (Step 2), or element results that differ from the README's table (Step 6).

Do not interpret the results beyond this: whether the skins pass level 3 is decided with the laptop's
`d100__sesam__skin` flight alongside, by Asha.
