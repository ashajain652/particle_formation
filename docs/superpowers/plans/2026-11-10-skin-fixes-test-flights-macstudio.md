# Sub-plan 22 Tasks 5g-5j: the two 100 mm test flights on the Mac Studio — runbook

> **For the agent running this.** You are on Asha's Mac Studio. Your job is to run two already-written melting
> flights of the 100 mm sphere, each expected to take one to three days, watch them, gather their results, push them on
> a results branch, and report. **You write no model code.** If anything in this runbook does not work as written, stop
> and tell Asha what happened rather than changing code or settings to make it work. Read the whole runbook once before
> running anything.

**Asked by:** Asha, 2026-11-10: re-fly the 100 mm sphere from its state saved at 22.0 s to 45 s with all of this week's
fixes, once with the slip band off and once with it on; run the flights on the Mac Studio (Asha's rule of 2026-10-31:
flights expected to take more than about two hours go there, side by side, one core each); push after the final review.

**Design of record:** `docs/superpowers/specs/2026-11-10-fitted-nose-radius-design.md` (all four fixes, with Asha's
decisions and their amendments) and the plan `docs/superpowers/plans/2026-11-10-fitted-nose-radius-build.md`, section
"After the final review", which names these two flights.

## 1. What this is and why it matters

The repository models what happens to the solid aluminium (AA7075) fragments of a re-entering satellite. Its newest
surface model, **layered skins** (also called the skyscrapers), puts a column of thin cells on every corner of the
surface mesh of a 3D heat model; molten metal flows in the columns and sprays off as droplets by Girin's model. On
2026-11-08 to 2026-11-10 four comparison flights on Asha's laptop showed problems that sub-plan 22 Tasks 5g to 5j fixed:

- **5g, the fitted nose radius.** The windward loads (stagnation heating, the surface flow's pressure gradient and shear)
  used the sphere's original radius whatever the shape. They now use a sphere fitted every macro step to the tops of the
  skin columns near the nose, bounded on the flat side by the flat-face equivalent (1.67 times the windward transverse
  radius, Boison and Curtiss 1959). Option `--size-feedback current`.
- **5h, new defaults and the Rayleigh-Taylor droplet.** For melting skins the re-frozen metal is no longer handed down
  into the core (`--hand-down off`) and the number of parts a macro step is split into is set by melting alone
  (`--part-count-from melting`); the Rayleigh-Taylor droplet is no longer capped at a quarter of the body's radius.
- **5i, the Knudsen-number gates and the slip band.** Girin and Kopyt's thin branch now sprays only in the continuum
  regime (`--spray-gates 2026-11-10`). With `--slip-band on`, both of Girin's branches may also spray in the slip band,
  where the body's Knudsen number (the ratio of the gas's mean free path to the body's size) is between 0.01 and 0.1 and
  the patch's own wall Knudsen number is below 0.1, with the melt layer's depth from a momentum balance.
- **5j, the leeward test by the centroid.** A surface corner counts as windward only if it lies ahead of the plane
  through the body's mass centroid (rigid and molten metal together) at right angles to the flight direction and also
  faces the flow; behind that plane there is no shear, no pressure gradient, no Girin spraying and no Rayleigh-Taylor
  spraying (`--leeward-test centroid`).

Both flights resume the 100 mm sphere's state saved at 22.0 s of flight (the step before its first liquid), on the
production mesh (9 417 skin columns), at 7.5 km/s from 77.5 km, flight-path angle -0.959331 degrees, US76 atmosphere,
physics-mode heating, and fly to 45 s. Over that window the body's Knudsen number falls from about 0.019 to about 0.011,
all inside the slip band, so **with the slip band off the flights can spray only by the Rayleigh-Taylor mode**, and with
it on Girin's shear spraying may also release droplets. That difference is the point of the pair.

| Label | Slip band | Ends | Saves its state |
|---|---|---|---|
| `d100__fixes__slipband-off` | off | at 45 s of flight | at 32.5 s |
| `d100__fixes__slipband-on` | on | at 45 s of flight | at 32.5 s |

**What the laptop measured for the same window** (the comparison flight `d100__handdown-off`, 2026-11-08 to 2026-11-10,
older code with the hand-down off, the original radius and the old spraying rules): 1 757 macro steps, about 36 hours
of computing (about 34 hours in the skins and 2.4 in the 3D solve, the log's running totals less the first 22 s) while sharing the laptop's four performance cores with
three other flights, peak memory 4.8 GB, no droplets released, and its frames (one per macro step) took about 1 GB. The
new flights should spray, so their droplet records will add disk and some time. **Expect one to three days each.**
Saving the state at 32.5 s writes about 1.4 GB per flight.

## 2. Rules

- Do not edit anything under `prototype/`, `reentry_model/`, `analysis/`, `tests/` or `docs/`, except creating the
  results directory of Step 9. The code is verified and its plan is generated from it byte for byte; a hand edit breaks
  that. If you think code is wrong, report it.
- Do not change any flight setting and do not re-run a flight with different settings. A flight that stops is a result
  to report, not an error to fix (Asha, 2026-10-31: a flight that stops on a fold of the surface keeps and judges what it
  made up to then).
- Do not change system settings (sleep, updates, energy). Asha keeps the machine from restarting while the flights run.
- Do not touch DRAMA, pyDRAMA, `sphere_reentry.py` or `sphere_sweep.py`; nothing here needs them.
- Run nothing else heavy on the machine while the flights run.
- Git: never use bare `git stash`. Push only the results branch named in Step 9, never `main` or
  `step3-rigid-substrate-scheil`. Do not merge anything. Never commit `*.state.pkl`, `particles.npz`, `source_rows/` or
  `vtk/` (all too large).
- Quote every path: the repository's path on Asha's machines contains spaces (for example `University Documents /MIT`).
- When in doubt, stop and ask Asha. Report in Asha's style: full sentences, every number with what was measured, its units
  and why it matters, terms spelled out, words rather than symbols in prose (Asha's `~/.claude/CLAUDE.md` gives the full
  rules). Use they/them for Asha.

## 3. Get the code

- [ ] **Step 1.** In the Mac Studio's clone of the repository (`git@github.com:ashajain652/particle_formation.git`):

```bash
cd "<the clone's root>"
git status --short          # must show no modified tracked files; if it does, stop and ask Asha
git fetch origin
git switch step3-rigid-substrate-scheil
git pull --ff-only
REPO="$(git rev-parse --show-toplevel)"
P="$REPO/prototype/work-2026-10-28-skin-spraying-plan"
CODE="$P/code"
FLIGHTS="$P/flights"
ls "$FLIGHTS/fly.py" "$FLIGHTS/launch.sh"          # both must exist
git merge-base --is-ancestor a772c86 HEAD && echo "code is at or after the reviewed commit"
git log --oneline -1                                # note this commit for the report
```

The `merge-base` line must print its message: `a772c86` is the commit that passed the final review of Tasks
5g-5j. `$CODE` is a complete, runnable copy of the model package. Every command below runs from `$CODE` unless it says
otherwise.

## 4. The saved state the flights start from

- [ ] **Step 2.** The flights resume from a file that is not in git (1.4 GB, an output). **Asha copies it** from the
laptop, where it is
`prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_snapshots/d100__physics__skin.state.pkl`
in the `subplan20` worktree, to the same path under the Mac Studio's clone. If it is not there, ask Asha for it and
wait. Then check it:

```bash
S100="$CODE/reentry_model_output/task5_snapshots/d100__physics__skin.state.pkl"
ls -l "$S100"               # 1380690172 bytes
shasum -a 256 "$S100"       # 512d217c40ada2a901f35bb305b5f242cae5da4ce7070a7c18f464f26a6eaa5f
```

Both must match exactly. If either differs, the copy is damaged or is the wrong file: stop and tell Asha. The file was
written on 2026-11-06 by the snapshot flight at commit 055c313, at 22.0 s of flight, macro step 44.

## 5. The Python environment

- [ ] **Step 3.** Find the Mac Studio's `drama_env`:

```bash
conda env list
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python      # adjust if the Mac Studio's path differs
export PY
"$PY" -m pip install pytest-xdist==3.8.0                    # the parallel test runner; a no-op if already there
"$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera, matplotlib, meshio, pytest, xdist; print(sys.version.split()[0], numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__, matplotlib.__version__, pytest.__version__, xdist.__version__)"
```

On the laptop this prints `3.12.14 2.5.3 1.18.1 12.0.2 4.15.2 5.3.0 3.2.0 3.11.1 9.1.1 3.8.0` (Python, numpy, scipy,
scikit-fem, gmsh, pyamg, cantera, matplotlib, pytest, pytest-xdist). Installing pytest-xdist is the only environment
change asked for; if an import fails, repair the environment as the level-3 runbook's Step 2 does
(`docs/superpowers/plans/2026-10-16-skin-level3-flights.md`) and say so in the report. A different **gmsh**,
**scikit-fem** or **numpy** version is worth noting in the report.

## 6. Check the machine before spending days on it

- [ ] **Step 4: the unit tier**, across six processes:

```bash
cd "$CODE"
PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest -m "not drama and not reference" -q -p no:cacheprovider -n 6 \
    --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py
```

Expected: **702 passed, 1 skipped** (the skip is the FEniCSx module: `dolfinx` is not in `drama_env`). On the laptop this took
about 30 minutes across four processes (48 when long runs shared the laptop). Any failure: stop and report it with its output.

- [ ] **Step 5: the smoke check.** The first 0.6 s of the slip-band-on flight (to 22.6 s), through the flight launcher,
from the saved state, with no frames:

```bash
SMOKE="$CODE/reentry_model_output/task5_fixes_smoke"
FIXES=(--hand-down off --part-count-from melting --size-feedback current --spray-gates 2026-11-10 \
       --leeward-test centroid --sideways-heat height)
"$FLIGHTS/launch.sh" "$SMOKE" "$SMOKE" smoke_slipband_on resume "$S100" "${FIXES[@]}" --slip-band on \
    --t-max 22.6 --frames-every 0
tail -n 30 "$SMOKE/smoke_slipband_on.log"
"$PY" -c "import json,sys; r=json.load(open(sys.argv[1]))['results']; print(r['end_reason'], r['n_macro_steps'], r['skin_columns'], r['released_mass_kg'], r['energy_balance_residual'], r['momentum_residual'])" "$SMOKE/smoke_slipband_on.json"
"$PY" -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1]))); r=rows[-1]; print(len(rows), r['time_s'], r['nose_radius_mm'], r['fitted_nose_radius_mm'], r['nose_bound_mm'], r['centroid_offset_mm'], r['leeward_corners'], r['surface_T_max_K'])" "$SMOKE/smoke_slipband_on.rows.csv"
```

Expected, as measured on the laptop at `a772c86` on 2026-11-10: the log ends `exit 0`; the first
line prints `t_max 54 9417 0.0 -1.1271905903431824e-12 2.1145612562630774e-17` (end reason, macro steps counted from
the start of the flight, skin columns, mass released in kilograms, energy-balance residual, momentum residual); the
second prints `2 22.6 49.759876263273746 49.759876263273746 83.74355777889993 -0.026606509784501024 28.0
899.1908971436608` (history rows, time in seconds, the radius the loads used, the fitted radius and the flat-face bound
in millimetres, the centroid's offset in millimetres, the corners counted leeward, the hottest surface point in
kelvin). The log's second line, beginning `resuming`, lists `part count from melting, hand-down off, sideways heat
height, size feedback current, leeward test centroid, spray gates 2026-11-10, slip band on`. On the laptop it took
about 4 minutes (251 s, peak memory 4.0 GB); the `wall time ... in the skins` in the log is the running total since
the start of the flight, carried in the saved state, not this run's. A different number of columns or of macro steps, a released mass or a fitted
radius that differs in its first three digits, or an energy residual above 1e-9: stop and report. A difference only in
the fourth digit or beyond: note it in the report and go on.

## 7. Run the two flights

- [ ] **Step 6: room.** `sysctl -n hw.perflevel0.physicalcpu hw.memsize` and `df -h "$CODE"`. The two flights need two
free cores, about 12 GB of memory between them, and at least 40 GB of free disk (frames about 1 GB, the saved state 1.4
GB, and the droplet records, whose size is not known in advance, per flight). If the machine has less, stop and tell
Asha.

- [ ] **Step 7: launch both, side by side.**

```bash
OUT="$CODE/reentry_model_output/task5_fixes"
LOGS="$CODE/reentry_model_output/task5_fixes_logs"
mkdir -p "$OUT" "$LOGS"
FIXES=(--hand-down off --part-count-from melting --size-feedback current --spray-gates 2026-11-10 \
       --leeward-test centroid --sideways-heat height)
launch() {   # label, then the resume arguments
    nohup "$FLIGHTS/launch.sh" "$OUT" "$LOGS" "$@" > /dev/null 2>&1 < /dev/null &
    echo "$1 pid $! started $(date '+%F %T')" | tee -a "$LOGS/launches.txt"
}
launch d100__fixes__slipband-off resume "$S100" "${FIXES[@]}" --slip-band off --t-max 45 --save-state-at 32.5 --frames-every 1
launch d100__fixes__slipband-on  resume "$S100" "${FIXES[@]}" --slip-band on  --t-max 45 --save-state-at 32.5 --frames-every 1
```

The flights must outlive the shell command that started them: if your harness ends background jobs when a command
returns, start each `launch` line with its own background-run facility (in Claude Code, `run_in_background`) instead.
A minute later confirm with `ps -Ao pid,etime,pcpu,rss,command | grep "flights/fly.py" | grep -v grep` that two flights
are running near 100 % of a core each. The second line of each log begins `resuming` and lists the options read back:
check that it says `hand-down off`, `part count from melting`, `sideways heat height`, `save state at 32.5 s`, `frames
every 1` and the flight's own slip-band setting, and report the line.

## 8. Watch them

- [ ] **Step 8.** Each flight appends a line to `$OUT/<label>.progress` every 20 macro steps (time, altitude, speed,
mass, mean temperature, hottest surface point, mass released so far, mean flowing depth, wall seconds). Check every few
hours:

```bash
for f in "$OUT"/*.progress; do echo "$(basename "$f"): $(tail -n 1 "$f")"; done
ps -Ao pid,etime,pcpu,rss,command | grep "flights/fly.py" | grep -v grep
df -h "$CODE" | tail -n 1
```

Estimate each flight's finish from the flight time still to go (to 45 s) and the wall seconds per flight second over the
last few lines, and give Asha a one-line estimate after the first several hours and then once a day. A process near 100
% of a core whose progress file has not grown for an hour or two is not a failure (20 macro steps of 12.5 ms can take
that long); note it and keep watching. Resident memory (the `rss` column, in kilobytes) above 15 GB for one flight, or
free disk below 15 GB: tell Asha.

When a flight's process has gone, read the end of its log (`tail -n 60 "$LOGS/<label>.log"`). The last lines are
`exit <code>` and the end time.

- **`exit 0`**: the flight ended normally (at 45 s); its `<label>.csv`, `<label>.json` and `<label>/` (the particle
  files and the frames) are in `$OUT`.
- **any other exit code, after the flight had started** (its `rows.csv` has more than a few rows): a model failure, for
  example a fold of the surface that the mesh rebuild refuses (the error reads "a corner move of up to ... degrades the
  mesh even after a rebuild"). **Do not retry and do not change settings.** Keep everything. Report the error text, the
  time of the last row in `<label>.rows.csv`, its largest `surface_fold_deg`, and the size of `<label>/source_rows/`.
- **exit 2 with no `rows.csv`**: a bad argument or a missing library; fix the environment only (Step 3), then launch that
  flight again.
- **no `exit` line, and the process has gone** (a reboot, or the machine slept): if `$OUT/<label>.state.pkl` exists (the
  flight passed 32.5 s), resume from it, which continues the flight to the last bit with the options it saved:
  `launch <label>__from32.5 resume "$OUT/<label>.state.pkl" --t-max 45 --frames-every 1`; otherwise launch the flight
  again from Step 7. Either way report that it was restarted, why, and from where.

## 9. Bring the results back

- [ ] **Step 9.** When both flights have ended, commit what they made (if one ends a day or more before the other,
commit it then and again at the end):

```bash
cd "$REPO"
R="$P/flights/macstudio-fixes"
mkdir -p "$R"
for f in "$OUT"/*.csv "$OUT"/*.json "$OUT"/*.progress; do cp "$f" "$R"/; done
for d in "$OUT"/*/; do l="$(basename "$d")"; mkdir -p "$R/$l"; cp "$d"/*.csv "$R/$l"/ 2>/dev/null; done
cp "$LOGS"/*.log "$LOGS"/launches.txt "$R"/
cp "$SMOKE"/smoke_slipband_on.log "$SMOKE"/smoke_slipband_on.json "$SMOKE"/smoke_slipband_on.rows.csv "$R"/
{ for f in "$OUT"/*/particles.npz "$OUT"/*.state.pkl; do [ -f "$f" ] && { ls -l "$f"; shasum -a 256 "$f"; }; done;
  du -sh "$OUT"/*/source_rows "$OUT"/*/vtk 2>/dev/null; } > "$R/large_files.txt"
{ sysctl -n machdep.cpu.brand_string hw.perflevel0.physicalcpu hw.ncpu hw.memsize; sw_vers;
  "$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera; print(sys.version); print(numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__)";
  git rev-parse HEAD; } > "$R/environment.txt"
du -sh "$R"            # expect tens of MB at most; the large files stay on the Mac Studio
git switch -c skin-fixes-flights-macstudio 2>/dev/null || git switch skin-fixes-flights-macstudio
git add "$R"
git commit -m "Sub-plan 22 Tasks 5g-5j test flights on the Mac Studio: <which flights have ended, and how each ended>"
git push -u origin skin-fixes-flights-macstudio
```

Push only this branch. The particle files, saved states, droplet records and frames stay on the Mac Studio;
`large_files.txt` records where they are and their checksums, and Asha fetches what they need.

## 10. Report to Asha

- [ ] **Step 10.** After each commit, write the report in Asha's style (full sentences; every number with what it
measures, its units and why it matters; terms spelled out; no symbols in prose). In this order:

1. Which flights ended and how (at 45 s or a failure), and whether the results are pushed (branch and commit).
2. For each flight, from its JSON's `results` and its history: the mass released (`released_mass_kg`), the number of
   droplets (`n_released`) and their median and largest radius (`r_median_um`, `r_max_um`, in micrometres); when the
   first droplet left (the first row with `n_released` above zero) and how the release was split over time.
3. For each flight, the shape and the loads: the radius the loads used (`nose_radius_mm`), the fitted radius
   (`fitted_nose_radius_mm`) and the flat-face bound (`nose_bound_mm`) at 22.5 s, 30 s, 35 s, 40 s and the end, and how
   many rows had a concave cap (`nose_cap_concave` equal to 1); the centroid's offset from the original centre
   (`centroid_offset_mm`) and the number of corners counted leeward (`leeward_corners`) at the same times.
4. For each flight, the mesh and the books: the largest fold (`surface_fold_max_deg`), the lowest triangle quality
   (`surface_quality_min`), the number of rebuilds, the largest `step_parts`, the largest held re-frozen mass
   (`held_mass_max_kg`) and depth (`held_depth_max_mm`), the energy residual, the momentum residual and the jump faults.
5. Whether each state was saved at 32.5 s (`results.state_saved_at_s` or `results.state_save_error`) and its size.
6. The cost: each flight's wall time and peak memory (the end of its log), how the wall time per flight second changed
   through the flight, and that the two ran side by side.
7. Anything unexpected: a restart, a stalled progress file, an environment difference from the laptop (Step 3), or a
   smoke check that differed from Step 5's expected values.

Do not interpret the results beyond this: the comparison with the laptop's flights and what the results mean for the
model are Asha's, with the laptop's agent.
