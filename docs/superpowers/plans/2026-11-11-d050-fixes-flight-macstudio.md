# Sub-plan 22 Tasks 5g-5j: the 50 mm test flight on the Mac Studio — runbook

> **For the agent running this.** You are on Asha's Mac Studio. Your job is to run one already-written melting flight
> of the 50 mm sphere, expected to take one to several days, watch it, gather its results, push them, and report.
> **You write no model code.** If anything here does not work as written, stop and tell Asha rather than changing code
> or settings. Read the whole runbook once before running anything. This runbook follows the 100 mm runbook
> (`docs/superpowers/plans/2026-11-10-skin-fixes-test-flights-macstudio.md`, on `step3-rigid-substrate-scheil`); its
> sections 1 and 2 (what the fixes are, and the rules) apply here unchanged. Read them there.

**Asked by:** Asha, 2026-11-11, after the two 100 mm flights ended normally: fly the 50 mm sphere on the Mac Studio with
the re-frozen metal kept in the skins (no hand-down), the slip band on and the part count set by melting alone, with
the other fixes of Tasks 5g-5j as in the 100 mm flights.

## 1. The flight

One flight, `d050__fixes__slipband-on`. It resumes the 50 mm sphere's state saved at **170.5 s** of flight (macro step
341) by the same snapshot flight (2026-11-06, commit 055c313) that saved the 100 mm state, on the production mesh (2 367
skin columns). The original flight: 50 mm AA7075 sphere at 7.5 km/s from 115 km, flight-path angle -0.959331 degrees,
US76 atmosphere, physics-mode heating. DRAMA's melting reference for this sphere has it melt between 77.1 and 73.1 km.

Options (the 100 mm flights' set, slip band on):
`--hand-down off --part-count-from melting --size-feedback current --spray-gates 2026-11-10 --leeward-test centroid
--sideways-heat height --slip-band on`, flown to the end of the flight (`--t-max 3600`: the flight ends on its own, at
demise or another end condition, long before), its state saved at **180 s** (`--save-state-at 180`), a frame every macro
step's frame interval (`--frames-every 1`).

**Cost, from the laptop:** the 50 mm comparison flight `d050__handdown-off` (older code, hand-down off, 2026-11-08 to
2026-11-09, sharing four cores with three other flights) took more than a day to go from 170.5 s to 187.0 s and was
getting slower (its last 0.25 s of flight took about 5 hours, with columns 20-60 mm deep). Expect **one to several
days** alone on a core. The smoke check used 0.54 GB; the saved state is about 0.3 GB.

## 2. Get the code

- [ ] **Step 1.** In the Mac Studio's clone:

```bash
cd "<the clone's root>"
git status --short          # no modified tracked files; if any, stop and ask Asha
git fetch origin
git switch skin-fixes-flights-macstudio
git pull --ff-only
REPO="$(git rev-parse --show-toplevel)"
P="$REPO/prototype/work-2026-10-28-skin-spraying-plan"
CODE="$P/code"
FLIGHTS="$P/flights"
git merge-base --is-ancestor a772c86 HEAD && echo "code is at or after the reviewed commit"
git log --oneline -1
```

This branch holds the same code as the 100 mm flights (`5b68cbe`) plus their results and this runbook.

## 3. The saved state

- [ ] **Step 2.** **Asha copies** from the laptop
`prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_snapshots/d050__physics__skin.state.pkl`
(in the `subplan20` worktree) to the same path under the Mac Studio's clone. If it is not there, ask Asha and wait.

```bash
S050="$CODE/reentry_model_output/task5_snapshots/d050__physics__skin.state.pkl"
ls -l "$S050"               # 290992357 bytes
shasum -a 256 "$S050"       # 3ac90ebe9a914bdd85544f214ee2ad800ad3ff9040fff5e282deb32d9d6409e0
```

Both must match exactly; otherwise stop and tell Asha.

## 4. Environment and checks

- [ ] **Step 3.** `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python; export PY` (as in the 100 mm runbook's Step 3;
the environment was checked on 2026-11-10 and needs nothing new). The unit tier passed on this code on 2026-11-10 and is
not re-run.

- [ ] **Step 4: the smoke check.** The first 0.6 s (to 171.1 s), no frames:

```bash
SMOKE="$CODE/reentry_model_output/d050_fixes_smoke"
FIXES=(--hand-down off --part-count-from melting --size-feedback current --spray-gates 2026-11-10 \
       --leeward-test centroid --sideways-heat height)
"$FLIGHTS/launch.sh" "$SMOKE" "$SMOKE" smoke_d050_slipband_on resume "$S050" "${FIXES[@]}" --slip-band on \
    --t-max 171.1 --frames-every 0
tail -n 30 "$SMOKE/smoke_d050_slipband_on.log"
"$PY" -c "import json,sys; r=json.load(open(sys.argv[1]))['results']; print(r['end_reason'], r['n_macro_steps'], r['skin_columns'], r['released_mass_kg'], r['energy_balance_residual'], r['momentum_residual'])" "$SMOKE/smoke_d050_slipband_on.json"
"$PY" -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1]))); r=rows[-1]; print(len(rows), r['time_s'], r['nose_radius_mm'], r['fitted_nose_radius_mm'], r['nose_bound_mm'], r['centroid_offset_mm'], r['leeward_corners'], r['surface_T_max_K'])" "$SMOKE/smoke_d050_slipband_on.rows.csv"
```

Expected, as measured on the laptop on 2026-11-11 (its code adds only Task 5k's diagnostic jump-fault change, which
does not touch these numbers): the log ends `exit 0`; the first line prints
`t_max 351 2367 0.0 -2.1686276564124965e-13 4.984814506545659e-17`; the second prints
`2 171.1 24.706279029862895 24.706279029862895 42.61478147927497 8.257488795135937e-05 7.0 898.1296711632019`. The log's
second line, beginning `resuming`, lists `part count from melting, hand-down off, sideways heat height, size feedback
current, leeward test centroid, spray gates 2026-11-10, slip band on`. On the laptop it took 87 s and 0.54 GB. A
different number of columns or macro steps, a fitted radius differing in its first three digits, or an energy residual
above 1e-9: stop and report. A difference only in the fourth digit or beyond: note it and go on.

## 5. Run the flight

- [ ] **Step 5: room.** One free core, about 8 GB of memory, at least 20 GB of free disk (`df -h "$CODE"`).

- [ ] **Step 6: launch.**

```bash
OUT="$CODE/reentry_model_output/task5_fixes_d050"
LOGS="$CODE/reentry_model_output/task5_fixes_d050_logs"
mkdir -p "$OUT" "$LOGS"
FIXES=(--hand-down off --part-count-from melting --size-feedback current --spray-gates 2026-11-10 \
       --leeward-test centroid --sideways-heat height)
nohup "$FLIGHTS/launch.sh" "$OUT" "$LOGS" d050__fixes__slipband-on resume "$S050" "${FIXES[@]}" --slip-band on \
    --t-max 3600 --save-state-at 180 --frames-every 1 > /dev/null 2>&1 < /dev/null &
echo "d050__fixes__slipband-on pid $! started $(date '+%F %T')" | tee -a "$LOGS/launches.txt"
```

The flight must outlive the command that started it (in Claude Code, use `run_in_background`). A minute later check
with `ps -Ao pid,etime,pcpu,rss,command | grep "flights/fly.py" | grep -v grep` that it runs near 100 % of a core, and
that the log's second line (`resuming ...`) says `hand-down off`, `part count from melting`, `slip band on`,
`save state at 180 s` and `frames every 1`; report the line.

## 6. Watch it

- [ ] **Step 7.** As the 100 mm runbook's section 8, with `$OUT`, `$LOGS` and the label above: check
`tail -n 1 "$OUT/d050__fixes__slipband-on.progress"`, the process and the disk every few hours; give Asha a one-line
finish estimate after the first several hours and then once a day (the flight's end is not known in advance: estimate
from the altitude falling towards about 73 km and the mass falling towards zero). The same rules on exits apply: `exit 0`
is a normal end; any other exit after the flight started is a result (do not retry, do not change settings; report the
error, the last row's time, its largest `surface_fold_deg` and the size of `source_rows/`); if the process vanished
without an `exit` line and `$OUT/d050__fixes__slipband-on.state.pkl` exists, resume it with
`resume "$OUT/d050__fixes__slipband-on.state.pkl" --t-max 3600 --frames-every 1` under the label
`d050__fixes__slipband-on__from180`, otherwise relaunch from Step 6, and report it.

## 7. Bring the results back

- [ ] **Step 8.** When the flight has ended (and, if it runs more than two days, once a day while it runs with what it
has so far, the progress and rows files):

```bash
cd "$REPO"
R="$P/flights/macstudio-fixes-d050"
mkdir -p "$R"
for f in "$OUT"/*.csv "$OUT"/*.json "$OUT"/*.progress; do [ -f "$f" ] && cp "$f" "$R"/; done
for d in "$OUT"/*/; do l="$(basename "$d")"; mkdir -p "$R/$l"; cp "$d"/*.csv "$R/$l"/ 2>/dev/null; done
cp "$LOGS"/*.log "$LOGS"/launches.txt "$R"/
cp "$SMOKE"/smoke_d050_slipband_on.log "$SMOKE"/smoke_d050_slipband_on.json "$SMOKE"/smoke_d050_slipband_on.rows.csv "$R"/
{ for f in "$OUT"/*/particles.npz "$OUT"/*.state.pkl; do [ -f "$f" ] && { ls -l "$f"; shasum -a 256 "$f"; }; done;
  du -sh "$OUT"/*/source_rows "$OUT"/*/vtk 2>/dev/null; } > "$R/large_files.txt"
git rev-parse HEAD > "$R/commit.txt"
git add "$R"
git commit -m "Sub-plan 22 Tasks 5g-5j 50 mm test flight on the Mac Studio: <how far it got and how it ended>"
git push origin skin-fixes-flights-macstudio
```

Push only this branch. The frames (`vtk/`, about 1 GB or more) stay on the Mac Studio; Asha copies them by hand.

## 8. Report to Asha

- [ ] **Step 9.** As the 100 mm runbook's Step 10 (same order and the same quantities: how it ended; droplets released,
their number, median and largest radius and when the first left; the radius the loads used, the fitted radius and the
flat-face bound at 171 s, 175 s, 180 s, 185 s and the end; folds, mesh quality, rebuilds, parts, held mass and depth,
energy and momentum residuals, jump faults; whether the state was saved at 180 s; wall time, peak memory and how the
wall time per flight second changed), in Asha's style. Do not interpret beyond that.
