# Kinked-profile dispersion Task 9: the four comparison flights on the Mac Studio — runbook

> **For the agent running this.** You are on Asha's Mac Studio. Your job is to run four already-written melting
> flights, two resumed from a saved state and two fresh, watch them, gather their results, push them on a results
> branch, and report. **You write no model code.** If anything in this runbook does not work as written, stop and tell
> Asha what happened rather than changing code or settings to make it work. Read the whole runbook once before running
> anything.

**Asked by:** Asha, 2026-11-11 (the kinked-profile dispersion's decisions) and 2026-10-31 (flights expected to take more
than about two hours go to the Mac Studio, side by side, one core each): fly the layered skins' spraying on both of its
theories, `table` (Girin's theory as it stood before) and `profile` (Girin's gradient mode solved on each column's own
speed profile), so that the two can be compared.

**Design of record:** `docs/superpowers/specs/2026-11-11-kinked-profile-dispersion-design.md` and its build plan; the
assumptions are `docs/model_assumptions.md` section 12.

## 1. What this is and why it matters

The repository models what happens to the solid aluminium (AA7075) fragments of a re-entering satellite. Its surface
model, **layered skins**, puts a column of thin cells on every corner of the surface mesh of a 3D heat model; molten metal
flows in the columns and sprays off as droplets by Girin's model. The branch `kinked-dispersion` adds a choice of the
theory the spraying solves, `--dispersion table|profile`:

- **`table`**: Girin's theory as it stood before, with his theoretical layer depth and surface speed (and
  `--we-critical`). The flight's results are bit for bit those of the code before the branch, apart from the JSON's
  `dispersion` key.
- **`profile`**, the default on a spraying skin run since Asha's decision of 2026-11-11: Girin's gradient mode solved on
  each windward column's own speed profile, fitted as straight segments, above a **wall** (the depth below which the melt
  does not flow like a layer; a floor test places it at the shallowest floor-like kink in the profile, else at the rigid
  front), and Girin and Kopyt's gas-driven mode on that wall, whichever is faster releasing the cells at least as warm as
  its threshold temperature.

The flights compare the two on the floor-like share (the fraction of resolved columns whose wall is a floor-like kink),
where the walls land against the 906.4 K isotherm (the temperature of the bridge between slurry and rigid metal), the mass
each mode sprayed (gradient, gas-driven, Rayleigh-Taylor and slurry) and the droplet sizes. The options below are
`--spray-gates 2026-11-10` and `--slip-band on` throughout, so that the only thing that differs inside a pair is
`--dispersion`.

| Label | Starts from | Dispersion | Ends | Wall time (an estimate, see below) |
|---|---|---|---|---|
| `d050__kd__table` (a) | the 50 mm sphere's state at 178.5 s | `table` | 179.0 s | about 1.3 h |
| `d050__kd__profile` (b) | the same state | `profile` | 179.0 s | about 1.3 h |
| `d100__kd__table` (c) | fresh, 100 mm, from 77.5 km | `table` | 30 s (state saved at 25 s) | about 4 to 6 h |
| `d100__kd__profile` (d) | fresh, 100 mm, from 77.5 km | `profile` | 30 s (state saved at 25 s) | about 4 to 6 h |

**Flights (a) and (b) are one macro step long.** The saved state is at 178.5 s, the macro step is 0.5 s, and the flight
ends at 179.0 s, so each writes one history row. That is what Asha asked for; it is enough to compare the first spray plan
of the two theories on the same saved body, not a spraying history. The d050 body released nothing to 179.5 s under the old
gates (the laptop's rows), so (a) and (b) may spray nothing; they compare the first spray plan on the same body, and (c) and
(d) carry the spraying comparison.

**Where the times come from (all of them estimates; none was measured on this branch's code).** The laptop's flight `d050__handdown-off` (older code, four flights sharing four cores)
spent 3 120 s of wall time in the macro step ending at 178.5 s, 4 643 s in the one ending at 179.0 s and 3 840 s in the one
ending at 179.5 s (the increase of `wall_skins_s + wall_3d_s + wall_move_s + wall_rebuild_s` in its `.rows.csv`); the
profile path adds about 0.8 s per macro step on top (measured by `analysis/kinked_dispersion/plan_cost.py`, Step 4).
The 100 mm snapshot flight `d100__physics__skin` (the same laptop, four flights sharing it) reached 24.15 s in 5 682 s and
32.9 s in 21 621 s of wall time; (c) and (d) stop at 30 s, so expect roughly 4 to 6 hours, which is an estimate from
that older code and not a measurement of this branch. **First liquid is at 24.5 s**; from there the macro step is
0.0125 s.

**Which can run at the same time.** All four. They share nothing: each has its own name, output folder and spool, and
each uses one core and 4 to 5 GB of memory (the laptop measured 4.8 GB for the 100 mm flight). Run all four side by side.
The only precaution is the mesh cache: the 100 mm skin mesh is cached under `output/meshes/` by geometry, so launch (c)
first and launch (d) once (c) has built it (Step 6).

## 2. Rules

- Do not edit anything under `prototype/`, `reentry_model/`, `analysis/`, `tests/` or `docs/`, except creating the results
  folder of Step 8. If you think code is wrong, report it.
- Do not change any flight setting and do not re-run a flight with different settings. A flight that stops is a result to
  report, not an error to fix.
- Do not change system settings (sleep, updates, energy). Asha keeps the machine from restarting while the flights run.
- Do not touch DRAMA, pyDRAMA, `sphere_reentry.py` or `sphere_sweep.py`; nothing here needs them.
- Run nothing else heavy on the machine while the flights run.
- Git: never use bare `git stash`. Push only the results branch named in Step 8, never `main`. Do not merge anything. Never
  commit `*.state.pkl`, `particles.npz`, `source_rows/` or `vtk/` (too large).
- Quote every path: the repository's path on Asha's machines contains spaces (for example `University Documents /MIT`).
- When in doubt, stop and ask Asha. Report in Asha's style: full sentences, every number with what was measured, its units
  and why it matters, terms spelled out, words rather than symbols in prose (Asha's `~/.claude/CLAUDE.md` gives the full
  rules). Use they/them for Asha.

## 3. Get the code

The branch `kinked-dispersion` must be **pushed first**, by the laptop's agent after Asha agrees; it was not pushed when
this runbook was written. If `git fetch` below does not show `origin/kinked-dispersion`, stop and tell Asha.

- [ ] **Step 1.** **Assume every fenced block below runs in a fresh shell:** no variable, array or function survives from
one block to the next. This step therefore writes an environment file, `$HOME/kinked-dispersion-flights/env.sh`
(outside the repository), holding every variable and the `launch` function, and **every later block begins by sourcing
it**. The blocks run in zsh (macOS's shell) and are also valid in bash. In the Mac Studio's clone of the repository
(`git@github.com:ashajain652/particle_formation.git`):

```bash
cd "<the clone's root>"
git status --short          # must show no modified tracked files; if it does, stop and ask Asha
git fetch origin
git switch kinked-dispersion
git pull --ff-only
test "$(git rev-parse HEAD)" = "$(git rev-parse origin/kinked-dispersion)" && echo "at the tip of the pushed branch"
REPO="$(git rev-parse --show-toplevel)"
CODE="$REPO/Reentry Model"
ls "$CODE/reentry_model/cli.py" "$CODE/analysis/kinked_dispersion/plan_cost.py" "$CODE/analysis/kinked_dispersion/fly.py"   # all three must exist
git rev-parse HEAD                                  # note this commit for the report

KD="$HOME/kinked-dispersion-flights"
mkdir -p "$KD/states" "$KD/flights" "$KD/logs"
ENVF="$KD/env.sh"
printf 'export REPO=%q\n' "$REPO" > "$ENVF"
cat >> "$ENVF" <<'EOF'
export CODE="$REPO/Reentry Model"
export KD="$HOME/kinked-dispersion-flights"
export STATES="$KD/states" OUT="$KD/flights" LOGS="$KD/logs"
export S050="$STATES/d050__handdown-off.state.pkl"
export PY="${PY:-/Users/ashajain/miniforge3/envs/drama_env/bin/python}"      # Step 3 says where to change it
FRESH=(--diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76
       --thermal fem --heating physics --melt on --surface-model skin
       --hand-down off --part-count-from melting --size-feedback current --leeward-test centroid --sideways-heat height
       --spray-gates 2026-11-10 --slip-band on --t-max 30 --frames-every 0 --save-state-at 25.0)
RESUME=(--spray-gates 2026-11-10 --slip-band on --t-max 179.0 --frames-every 0)
launch() {   # launch <label> <run arguments...>   or   launch <label> resume <state> <options...>
    nohup zsh "$KD/fly.sh" "$@" > /dev/null 2>&1 < /dev/null &
    echo "$1 pid $! started $(date '+%F %T')" | tee -a "$LOGS/launches.txt"
}
EOF
cat > "$KD/fly.sh" <<'EOF'
#!/bin/zsh
# fly.sh <label> <run arguments... | resume <state> options...>: one flight through fly.py, from $CODE, no bytecode written.
# It needs the exported variables of env.sh (the launch function exports them to it).
label="$1"; shift
log="$LOGS/$label.log"
echo "start $(date '+%F %T')" > "$log"
cd "$CODE" && PYTHONDONTWRITEBYTECODE=1 caffeinate -i /usr/bin/time -l "$PY" -B -m analysis.kinked_dispersion.fly "$OUT" "$label" "$@" >> "$log" 2>&1
echo "exit $?" >> "$log"
echo "end $(date '+%F %T')" >> "$log"
EOF
cat "$ENVF" | head -3                                # the first line names the clone's root
```

The tip check must print its message: the expected commit is the tip of the pushed `kinked-dispersion`, not a fixed
number; the commit it prints is the one to put in the report. `$CODE` is the `Reentry Model/` folder, the package's root.
(Unlike the sub-plan 22 flights, these use the package itself, not a planning copy's `flights/launch.sh`: the flight
launcher is `analysis/kinked_dispersion/fly.py`, the sub-plan 22 `fly.py` ported to the package, wrapped by `fly.sh`,
which logs the start, the end and the exit code, `/usr/bin/time -l`'s wall time and peak memory, and runs under
`caffeinate -i` so that the machine does not sleep.)

## 4. The saved state flights (a) and (b) start from

- [ ] **Step 2.** The state is not in git (270 MB, an output). **Asha copies it by hand** from the laptop, where it is, in
the `subplan20` worktree,

`prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/d050__handdown-off.state.pkl`

to this folder on the Mac Studio (create it):

`$HOME/kinked-dispersion-flights/states/d050__handdown-off.state.pkl`

**Asha also copies `d050__handdown-off.rows.csv`** (0.8 MB, from the same laptop folder) **to the same folder**; Step 4's
timing check reads it. The other files of that stem are **not** needed: a resume loads the whole run from the state
alone and writes its droplet records afresh (the `.progress` file is only a log, and the directory `d050__handdown-off/`
holds the old flight's droplet records, which a resume replaces). If the files are not there, ask Asha for them and
wait. Then check the state:

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
ls -l "$S050"               # 270157216 bytes
shasum -a 256 "$S050"       # 89aa20b54c1d20daa41984c68b59257abc994070c2e4a49fbdfb44724f180842
```

Both must match exactly. If either differs, the copy is damaged or is the wrong file: stop and tell Asha. The file was
written on 2026-11-09 by the laptop's flight `d050__handdown-off` at 178.5 s of flight, macro step 357 (50 mm sphere,
7.5 km/s from 115 km, hand-down off, sideways heat height, saved by code from before the dispersion option). A resume
keeps every setting the state saved except the ones the command line gives, so (a) and (b) change exactly the spray gates,
the slip band, the dispersion, the end time and the frames.

Fresh flights (c) and (d) need no file.

## 5. The Python environment, and the timing check

- [ ] **Step 3.** Find the Mac Studio's `drama_env`:

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
conda env list
echo "$PY"                  # env.sh's default is the laptop's path; if the Mac Studio's drama_env is elsewhere, edit the
                            # PY line of "$HOME/kinked-dispersion-flights/env.sh" and source it again
"$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera, matplotlib, meshio; print(sys.version.split()[0], numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__, matplotlib.__version__)"
```

On the laptop this prints `3.12.14 2.5.3 1.18.1 12.0.2 4.15.2 5.3.0 3.2.0 3.11.1` (Python, numpy, scipy, scikit-fem, gmsh,
pyamg, cantera, matplotlib). If an import fails, repair the environment as the level-3 runbook's Step 2 does
(`docs/superpowers/plans/2026-10-16-skin-level3-flights.md`) and say so in the report. A different **gmsh**,
**scikit-fem** or **numpy** version is worth noting in the report.

- [ ] **Step 4: the timing check** (a few minutes; it flies nothing): how long the spray plan takes on each theory on the
saved 50 mm body, against a macro step's wall time.

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
cd "$CODE"
PYTHONDONTWRITEBYTECODE=1 "$PY" -m analysis.kinked_dispersion.plan_cost --state "$S050"
```

Expected on the laptop on 2026-11-11: it prints `girin_plan on table` with a median of about 0.025 s, `girin_plan on
profile` with a median of about 0.815 s, an added time of about 0.79 s, a macro step wall time (median over the steps
from 177.5 to 179.5 s) of 3 480 s and `0.023 %` of it, against a budget of 5 %, then `PASS`. The Mac Studio's times may
differ; the check passes if the last line says `PASS`. A `FAIL`, or a traceback: stop and report it with its output. Note
the Mac Studio's two medians and its macro step wall time in the report.

(The package's unit tests are not part of this runbook: the laptop ran them on the branch. If Asha asks for them, the
command is in `CLAUDE.md`, "Tests".)

## 6. Run the four flights

- [ ] **Step 5: room.** `sysctl -n hw.perflevel0.physicalcpu hw.memsize` and `df -h "$HOME"` (the flights write under
`$HOME/kinked-dispersion-flights`). The four flights need four free cores, about 20 GB of memory between them (the laptop
measured 4.8 GB peak for the 100 mm flight, so this is a floor) and **at least 60 GB of free disk**. That figure is an
estimate, built from these parts: no frames are written; each fresh flight saves a state at 25 s, about 1.4 GB (an
estimate from the laptop's earlier 100 mm flights; `cli.py` says about 1.3 GB); and the droplet records go to a spool,
`<label>/source_rows`, as they are made. `cli.py` (lines 66 to 67) gives the 100 mm flight 35 to 70 million records at the
continuum step, 7.6 to 15.2 GB in memory, without naming the flight's length; **I assumed it means a flight to about 45 s,
liquid from 24.5 s**, and scaled by the liquid window: 5.5 s (24.5 to 30 s) against about 20.5 s (24.5 to 45 s) is about
0.27, which gives about 9 to 19 million records and about 2 to 4 GB of spool for each of (c) and (d). Nothing was measured
on this branch, so the 60 GB uses the unscaled 15.2 GB for each of the two spools (30 GB), the two states (3 GB), the
histories, particle files and logs (a few GB) and the two 50 mm resumes (small), with margin. If the machine has less, stop
and tell Asha.

- [ ] **Step 6: launch.** All four flights go through `analysis/kinked_dispersion/fly.py`, which runs the command line
in its own process and keeps each history row in `<label>.rows.csv` and a progress line every 20 rows in
`<label>.progress` as the flight runs. The common options of the fresh flights are the sub-plan 22 flights' (the 100 mm
sphere, 7.5 km/s from 77.500133 km, flight path angle -0.959331 degrees, US76 atmosphere, physics-mode heating) with the
fixes of 2026-11-10; the spray options are the same inside a pair. The fresh flights save their whole state at the end of
the first macro step at or after 25 s (`--save-state-at 25.0`, just after the first liquid at 24.5 s) as
`$OUT/<label>.state.pkl`, so that a failure after that can resume from there instead of flying again.

The flights must outlive the shell command that started them: if your harness ends background jobs when a command
returns, run each block below with its own background-run facility (in Claude Code, `run_in_background`) instead. Each
block is self-contained. Launch (c) first, alone, so that no two processes build one mesh file at once:

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
launch d100__kd__table "${FRESH[@]}" --dispersion table
```

Then wait until (c) has built the mesh and made its first history row (its `.progress` file then has a line), for at most
30 minutes, and stop if the flight has already ended:

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
i=0
while [ ! -s "$OUT/d100__kd__table.progress" ] && ! grep -q '^exit ' "$LOGS/d100__kd__table.log" && [ "$i" -lt 180 ]; do
    sleep 10; i=$((i + 1))
done
if [ -s "$OUT/d100__kd__table.progress" ]; then echo "mesh built, (c) is making rows"
else echo "STOP: (c) made no row (ended or 30 minutes passed); read its log"; tail -n 30 "$LOGS/d100__kd__table.log"; fi
```

If it says STOP, do not launch the others: read the log and tell Asha. Otherwise launch (d):

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
launch d100__kd__profile "${FRESH[@]}" --dispersion profile
```

The two 50 mm resumes load their mesh from the state and need no wait:

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
launch d050__kd__table   resume "$S050" "${RESUME[@]}" --dispersion table
launch d050__kd__profile resume "$S050" "${RESUME[@]}" --dispersion profile
```

A minute later confirm with `ps -Ao pid,etime,pcpu,rss,command | grep "analysis.kinked_dispersion.fly" | grep -v grep`
that four flights are running near 100 % of a core each (the `caffeinate` and `time` wrappers also show; count the
`python` lines). The first line of each resume log after `start` begins `resuming` and lists the options read back: check
that it says `spray gates 2026-11-10`, `slip band on`, `frames every 0`, the flight's own `dispersion` setting and `to 179
s`, and report the line for (a) and (b). For (c) and (d) the log's first lines are the `run` summary only at the end, so
check the JSON's `settings` and `thermal` entries then (Step 7); meanwhile `tail -n 2 "$OUT"/*.progress` shows that each
flight is making rows.

## 7. Watch them

- [ ] **Step 7.** Each flight appends every history row to `$OUT/<label>.rows.csv` as it is made, and a line to
`$OUT/<label>.progress` every 20 rows (time, altitude, speed, mass, mean temperature, hottest surface point, mass
released so far, mean flowing depth, and on `profile` flights the floor-like share, then wall seconds). Check every few
hours:

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
for f in "$OUT"/*.progress; do echo "$(basename "$f"): $(tail -n 1 "$f")"; done
ps -Ao pid,etime,pcpu,rss,command | grep "analysis.kinked_dispersion.fly" | grep -v grep
du -sh "$OUT"/*/source_rows 2>/dev/null          # the droplet spool, which grows once liquid sprays (from 24.5 s)
df -h "$HOME" | tail -n 1
```

Estimate each flight's finish from the flight time still to go (to 30 s, or to 179.0 s for (a) and (b)) and the wall
seconds per flight second over the last few lines, and tell Asha in one line after the first hour, then every few hours.
A process near 100 % of a core whose progress file has not grown for an hour or two is not a failure (20 macro steps of
12.5 ms can take that long); note it and keep watching. Resident memory (the `rss` column, in kilobytes) above 15 GB for
one flight, or free disk below 15 GB: tell Asha.

When a flight's process has gone, read the end of its log (`tail -n 60 "$LOGS/<label>.log"`). The last lines are
`exit <code>` and the end time, and `/usr/bin/time -l`'s wall time and `maximum resident set size` just above them.

- **`exit 0`**: the flight ended normally. Its `<label>.csv` (the history), `<label>.json` (settings, results and
  provenance) and `<label>/` (`particles.npz`, `particles_summary.csv`, `size_distribution.csv`) are in `$OUT`. Check
  in the JSON that `settings.spray_gates` is `2026-11-10`, `settings.slip_band` is `on` and `settings.dispersion` is the
  flight's own (`table` or `profile`).
- **any other exit code, after the flight had started** (its `rows.csv` has more than a few rows): a model failure, for
  example a fold of the surface that the mesh rebuild refuses (the error reads "a corner move of up to ... degrades the
  mesh even after a rebuild"). **Do not retry and do not change settings.** Keep everything: the history to its last
  completed macro step stays in `<label>.rows.csv`, and the droplet records made so far stay in `<label>/source_rows/`
  (the droplet spool survives an early stop; it is removed only once the particle files are written). Report the error
  text, the time of the last row in `<label>.rows.csv`, and the size of `<label>/source_rows/`.
- **exit 2 with no files**: a bad argument (fly.py passes argparse's own exit code through) or a missing library; fix the
  environment only (Step 3), then launch that flight again.
- **no `exit` line, and the process has gone** (a reboot, or the machine slept): for (a) and (b), which save no state,
  launch the flight again from Step 6. For (c) and (d), if `$OUT/<label>.state.pkl` exists (the flight passed 25 s), it
  resumes from there, with the options it saved, through fly.py's resume form (a new label, so that nothing is
  overwritten):
  `source "$HOME/kinked-dispersion-flights/env.sh"; launch <label>__from25 resume "$OUT/<label>.state.pkl" --t-max 30 --frames-every 0`;
  otherwise launch the flight again from Step 6. Either way report that it was restarted, why, and from where. A flight
  that stopped on a model failure is not restarted (above). **A `__from25` resume starts a new droplet spool**, so the
  droplet records from 24.5 s to 25 s stay in the first run's `<label>/source_rows/` and the resumed run's
  `particles.npz` holds only the records after 25 s: bring back both (the first run's spool is kept only if that run did
  not finish; tell Asha, who will combine them), and say so in the report.

## 8. Bring the results back

- [ ] **Step 8.** When all four flights have ended, commit what they made (if some end hours before the others, commit
them then and again at the end). The small files go in git; the large ones stay on the Mac Studio:

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
setopt nullglob 2>/dev/null || shopt -s nullglob      # zsh: setopt; bash: shopt. A glob that matches nothing then expands to nothing
cd "$REPO"
R="$CODE/reports/kinked_dispersion_flights_macstudio"
mkdir -p "$R"
for f in "$OUT"/*.csv "$OUT"/*.json "$OUT"/*.progress; do cp "$f" "$R"/; done
for d in "$OUT"/*/; do l="$(basename "$d")"; mkdir -p "$R/$l"; cp "$d"/particles_summary.csv "$d"/size_distribution.csv "$R/$l"/ 2>/dev/null; done
cp "$LOGS"/*.log "$LOGS"/launches.txt "$R"/
{ for f in "$OUT"/*/particles.npz "$OUT"/*.state.pkl; do [ -f "$f" ] && { ls -l "$f"; shasum -a 256 "$f"; }; done;
  du -sh "$OUT"/*/source_rows 2>/dev/null; } > "$R/large_files.txt"
{ sysctl -n machdep.cpu.brand_string hw.perflevel0.physicalcpu hw.ncpu hw.memsize; sw_vers;
  "$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera; print(sys.version); print(numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__)";
  git rev-parse HEAD; } > "$R/environment.txt"
du -sh "$R"            # expect a few MB at most
```

What comes back, per flight, and why:

- `<label>.rows.csv`, the history as the flight made it, and `<label>.progress` (these also exist for a flight that
  stopped early), and `<label>.csv`, the macro-step history written at the end. On `profile` flights it carries the history columns of
  `melt_profile.DISPERSION_COLUMNS`: `floor_like_share` (the share of resolved columns whose wall is a floor-like kink),
  `fluid_like_columns` (the count of columns that are fluid-like), `depth_to_wall_median_mm` (the median depth to the
  wall in millimetres), `floor_re_median` (the median floor Reynolds number), `we_s_profile_median` (the median surface
  Weber number the profile gave) and `sprayed_gradient_kg` and `sprayed_gas_driven_kg` (the mass sprayed by the gradient
  mode and by the gas-driven mode, in kilograms). Also `released_mass_kg`, `n_released`, `r_median_um` and `r_max_um`.
  On `table` flights those seven columns are absent.
- `<label>.json`: the run's settings, results (mass released, droplets, residuals) and provenance, with `settings.dispersion`.
- `<label>.log`: the end of the flight, its wall time and peak memory.
- `<label>/particles_summary.csv` and `size_distribution.csv`: droplets per macro step and per size bin.
- `<label>/particles.npz`, **the droplet records, one row per released patch** (27 columns on `table`; 36 on `profile`,
  the 27 followed by `wall_depth_m`, `wall_floor`, `floor_re`, `kink_depth_m`, `kink_T_K`, `release_depth_m`,
  `growth_rate_per_s`, `girin_delta_m_m`, `girin_v_s_ms`). The mass by mode and the droplet sizes come from it: `branch`
  (0 Girin thick, 1 thin, 2 rarefied, 3 Rayleigh-Taylor, 4 regime 2, 5 the gradient mode on the profile), `dm_kg` (mass
  released), `r_m` (droplet radius) and `solid_fraction` (above zero for slurry). **It is not committed**: it can be tens of
  MB. Check its size in `large_files.txt`; if it is under 50 MB per flight, copy the `particles.npz` files into the
  results folder too (`cp "$OUT/<label>/particles.npz" "$R/<label>/"`) and commit them, otherwise tell Asha where they are
  and how big they are, and Asha fetches them (Asha's report needs them).

```bash
source "$HOME/kinked-dispersion-flights/env.sh"
cd "$REPO"
R="$CODE/reports/kinked_dispersion_flights_macstudio"
git switch -c kinked-dispersion-flights-macstudio 2>/dev/null || git switch kinked-dispersion-flights-macstudio
git add "$R"
git commit -m "Kinked-profile dispersion Task 9 comparison flights on the Mac Studio: <which flights have ended, and how each ended>"
git push -u origin kinked-dispersion-flights-macstudio
```

Push only this branch. Particle files over the size above, droplet spools and the saved states at 25 s of (c) and (d)
(about 1.4 GB each, an estimate; never committed) stay on the Mac Studio; `large_files.txt` lists them with their checksums.

## 9. Report to Asha

- [ ] **Step 9.** After each commit, write the report in Asha's style (full sentences; every number with what it
measures, its units and why it matters; terms spelled out; no symbols in prose). In this order:

1. Which flights ended and how (at their end time or a failure), whether any was resumed from its 25 s state (and then
   that its first 0.5 s of droplet records are in a separate spool), and whether the results are pushed (branch and commit).
2. For each flight, from its JSON's `results` and its history: the mass released (`released_mass_kg`), the number of
   droplets (`n_released`), their median and largest radius (`r_median_um`, `r_max_um`, in micrometres) and, for (c) and
   (d), when the first droplet left (the first row with `n_released` above zero).
3. For the `profile` flights, the `DISPERSION_COLUMNS` at the first spraying step and at the end: the floor-like share,
   the number of fluid-like columns, the median depth to the wall, the median floor Reynolds number, the median surface
   Weber number, and the mass sprayed by the gradient and the gas-driven modes.
4. The cost: each flight's wall time and peak memory (the end of its log), the timing check's two medians (Step 4), and
   that the four ran side by side.
5. Anything unexpected: a restart, a failure, an environment difference from the laptop (Step 3), or a timing check that
   differed from Step 4's expected values.

Do not interpret the results beyond this: the comparison of the two theories (the floor-like share, where the walls land
against the 906.4 K isotherm, the mass by mode and the droplet sizes) is Asha's, with the laptop's agent.
