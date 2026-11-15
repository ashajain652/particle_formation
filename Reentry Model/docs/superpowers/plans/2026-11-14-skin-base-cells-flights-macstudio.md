# The skins' base cells: the two confirmation flights on the Mac Studio — runbook

> **For the agent running this.** You are on Asha's Mac Studio. Your job is to run two already-written melting flights,
> both resumed from one saved state, side by side; watch them; gather their results; push them on a results branch; and
> report. **You write no model code.** If anything in this runbook does not work as written, stop and tell Asha what
> happened rather than changing code or settings to make it work. Read the whole runbook once before running anything.

**Asked by:** Asha, 2026-11-14 ("Keep 0.2 mm; write the Mac Studio runbook for the flight"), and 2026-10-31 (flights
expected to take more than about two hours go to the Mac Studio, side by side, one core each).

**Design of record:** there is no separate design. The change, its single-column test and what is assumed are in
`docs/model_assumptions.md` section 16 (section 14 when this runbook was written) and fact 198 of `docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`;
open item 7 of `docs/superpowers/plans/2026-11-12-skyscrapers-open-items.md` is the item these flights finish.

## 1. What this is and why it matters

The repository models what happens to the solid aluminium (AA7075) fragments of a re-entering satellite. Its surface
model, **layered skins**, puts a column of thin cells on every corner of the surface mesh of a 3D heat model; the cells
hold the metal above the **rigid front** (the depth where the alloy is half solid, 895.1 K), and metal crosses between the
skins and the 3D model there. A skin's cells used to be sized from its top only: 10 micrometres in the top 0.4 mm, each
cell below 10 % thicker than the one above, so a 28 mm skin ends in a bottom cell 2 to 4 mm thick. The branch
`skin-base-cells` also grades the cells up from the skin's base, from a **base cell** of 0.2 mm (`--skin-base-cell`,
0.2 mm by default; `0` is the old grading), so a thick skin is fine at both ends.

Why: a laptop re-run of the 50 mm sphere's flight from 178.5 s to 186.5 s (older code, 11.7 hours) counted 30,848
**jump faults** under the old rule and 25,620 under Task 5k's rule. A jump fault is a column whose skin's bottom cell lies
more than 5 K outside the temperatures of the metal that crossed its base (`skinbody.jump_faults`); every one of these was a
bottom cell about 10 K too warm, in skins 24 to 32 mm thick. A one-dimensional test of the same kind of skin
(`analysis/skin_base_cells/`) then showed that the old grading puts the rigid front 17 to 64 % too far, and that the
0.2 mm base cell brings it within 1.2 %. Nothing has yet been flown with the base cell. These two flights are the first:
they show whether the faults go and what the finer cells cost on a real body.

| Label | Starts from | Base cell | Ends | Wall time (an estimate, see below) |
|---|---|---|---|---|
| `d050__base-off` (a) | the 50 mm sphere's state at 178.5 s | none (`--skin-base-cell 0`, the old grading) | 186.5 s | about 10 to 15 h |
| `d050__base-0.2mm` (b) | the same state | 0.2 mm (`--skin-base-cell 0.2`) | 186.5 s | about 10 to 20 h |

Both flights change, against the saved flight, exactly these: the base cell, the end time (186.5 s), the geometry stop
(off, so that both fly the whole span the laptop's re-run flew; `--geometry-stop` did not exist when the state was
saved), the frames (none), and a state saved once at 182.5 s (`--save-state-at 182.5`, so that a reboot does not cost the
first half). They also give `--hand-down off --sideways-heat height`, as the laptop's re-run did; the state was saved with
the same. Every other setting is the saved flight's own, as `resume` keeps it.

**What this branch does not contain.** It is `main` (commit 4e857e1) plus the base cells. The skins' windward loads from
a fitted sphere (the branch `fitted-sphere-loads`, `--windward-loads`) and the interim gas loads (the branch
`interim-gas-loads`) are not on `main` and not on this branch; these flights use the windward loads the saved flight used.

**Where the times come from (estimates; nothing was measured on this branch's code on the 50 mm body).** The laptop's
re-run of the same span, on older code, took 42,176 s (11.7 h), sharing the laptop with a test run for its first 2 hours.
On 2026-11-14 another session resumed the same state on code built on `main` for four macro steps of 0.5 s, 178.5 to
180.5 s: 8,148 s and 12,519 s on the laptop for its two settings, the first macro step the slowest (4,267 and 5,009 s).
The base cell gives a 28 mm skin about a quarter more cells (120 to 126 against 100 to 103 on the one-dimensional test),
and the skins' solve loops over the most cells of any column, so (b) may take up to about a quarter longer than (a).
That is a guess, not a measurement; measuring it is one of these flights' purposes.

**They run at the same time.** They share nothing but the state file, which each only reads: each has its own label,
files and log, and uses one core and a few gigabytes of memory (the laptop's re-run peaked at 2.6 GB).

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

The branch `skin-base-cells` must be **pushed first**, by the laptop's agent after Asha agrees; it was not pushed when
this runbook was written. If `git fetch` below does not show `origin/skin-base-cells`, stop and tell Asha.

- [ ] **Step 1.** **Assume every fenced block below runs in a fresh shell:** no variable, array or function survives from
one block to the next. This step therefore writes an environment file, `$HOME/skin-base-cells-flights/env.sh` (outside
the repository), holding every variable and the `launch` function, and **every later block begins by sourcing it**. The
blocks run in zsh (macOS's shell) and are also valid in bash. In the Mac Studio's clone of the repository
(`git@github.com:ashajain652/particle_formation.git`):

```bash
cd "<the clone's root>"
git status --short          # must show no modified tracked files; if it does, stop and ask Asha
git fetch origin
git switch skin-base-cells
git pull --ff-only
test "$(git rev-parse HEAD)" = "$(git rev-parse origin/skin-base-cells)" && echo "at the tip of the pushed branch"
REPO="$(git rev-parse --show-toplevel)"
CODE="$REPO/Reentry Model"
ls "$CODE/reentry_model/cli.py" "$CODE/analysis/flights/launch.sh" "$CODE/analysis/skin_base_cells/column.py"   # all three must exist
grep -c "skin-base-cell" "$CODE/reentry_model/cli.py"            # must print a number above 0
git rev-parse HEAD                                  # note this commit for the report

SB="$HOME/skin-base-cells-flights"
mkdir -p "$SB/states" "$SB/flights" "$SB/logs"
ENVF="$SB/env.sh"
printf 'export REPO=%q\n' "$REPO" > "$ENVF"
cat >> "$ENVF" <<'EOF'
export CODE="$REPO/Reentry Model"
export SB="$HOME/skin-base-cells-flights"
export STATES="$SB/states" OUT="$SB/flights" LOGS="$SB/logs"
export S050="$STATES/d050__handdown-off.state.pkl"
export PY="${PY:-/Users/ashajain/miniforge3/envs/drama_env/bin/python}"      # Step 3 says where to change it
RESUME=(--t-max 186.5 --hand-down off --sideways-heat height --geometry-stop off --frames-every 0 --save-state-at 182.5)
launch() {   # launch <label> resume <state> <options...>
    nohup zsh "$CODE/analysis/flights/launch.sh" "$OUT" "$LOGS" "$@" > /dev/null 2>&1 < /dev/null &
    echo "$1 pid $! started $(date '+%F %T')" | tee -a "$LOGS/launches.txt"
}
EOF
cat "$ENVF" | head -3                                # the first line names the clone's root
```

The tip check must print its message: the expected commit is the tip of the pushed `skin-base-cells`, not a fixed number;
the commit it prints is the one to put in the report. `$CODE` is the `Reentry Model/` folder, the package's root. The
flights go through the package's own launcher, `analysis/flights/launch.sh` (its `README.md` describes it), which runs
`analysis/flights/fly.py` under `caffeinate -i` (no idle sleep) and `/usr/bin/time -l` (wall time and peak memory at the
end of the log), keeps every history row in `<label>.rows.csv` as it is made, and logs the start, the end, the exit code
and the power source in `$LOGS/<label>.log`. It reads `$PY`, which `env.sh` exports.

## 4. The saved state both flights start from

- [ ] **Step 2.** The state is not in git (270 MB, an output). If the kinked-profile dispersion flights' copy of
2026-11-12 is still on this machine, copy it from there:

```bash
source "$HOME/skin-base-cells-flights/env.sh"
OLD="$HOME/kinked-dispersion-flights/states/d050__handdown-off.state.pkl"
[ -f "$OLD" ] && cp -p "$OLD" "$S050" && echo "copied from the kinked-dispersion flights' folder"
ls -l "$S050"
```

If there is no such file, **Asha copies it by hand** from the laptop, where it is, in the `subplan20` worktree,

`prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/d050__handdown-off.state.pkl`

to `$HOME/skin-base-cells-flights/states/d050__handdown-off.state.pkl` on the Mac Studio. Ask Asha for it and wait. No
other file of that stem is needed. Then check the state:

```bash
source "$HOME/skin-base-cells-flights/env.sh"
ls -l "$S050"               # 270157216 bytes
shasum -a 256 "$S050"       # 89aa20b54c1d20daa41984c68b59257abc994070c2e4a49fbdfb44724f180842
```

Both must match exactly (checked on the laptop on 2026-11-14). If either differs, the copy is damaged or is the wrong
file: stop and tell Asha. The file was written on 2026-11-09 by the laptop's flight `d050__handdown-off` at 178.5 s of
flight, macro step 357 (50 mm sphere, 7.5 km/s from 115 km, hand-down off, sideways heat height), by code from before the
base cell existed, so it holds the old grading until a resume's `--skin-base-cell` changes it.

## 5. The Python environment

- [ ] **Step 3.** Find the Mac Studio's `drama_env`:

```bash
source "$HOME/skin-base-cells-flights/env.sh"
conda env list
echo "$PY"                  # env.sh's default is the laptop's path; if the Mac Studio's drama_env is elsewhere, edit the
                            # PY line of "$HOME/skin-base-cells-flights/env.sh" and source it again
"$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera, matplotlib, meshio; print(sys.version.split()[0], numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__, matplotlib.__version__)"
```

On the laptop this prints `3.12.14 2.5.3 1.18.1 12.0.2 4.15.2 5.3.0 3.2.0 3.11.1` (Python, numpy, scipy, scikit-fem, gmsh,
pyamg, cantera, matplotlib). If an import fails, repair the environment as the level-3 runbook's Step 2 does
(`docs/superpowers/plans/2026-10-16-skin-level3-flights.md`) and say so in the report. A different **gmsh**, **scikit-fem**
or **numpy** version is worth noting in the report.

- [ ] **Step 4: the column check** (a few seconds; it flies nothing). The single-column test's own unit test, which runs
the branch's skins on a 28 mm pool against a one-dimensional reference:

```bash
source "$HOME/skin-base-cells-flights/env.sh"
cd "$CODE"
PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests/test_skin_base_cells_column.py -q -p no:cacheprovider
```

Expected: `2 passed` (in about 4 s on the laptop). Anything else: stop and report it with its output. (The package's
whole unit tier is not part of this runbook: the laptop ran it on the branch, 1,004 passed and 1 skipped, with the one
failure it showed fixed and re-run. If Asha asks for it, the command is in `CLAUDE.md`, "Tests".)

## 6. Run the two flights

- [ ] **Step 5: room.** `sysctl -n hw.perflevel0.physicalcpu hw.memsize` and `df -h "$HOME"`. The two flights need two
free cores, about 10 GB of memory between them (the laptop's re-run peaked at 2.6 GB; this allows for the finer cells
and margin) and **at least 10 GB of free disk** (each saves its state twice, at 182.5 s and at its end, about 0.3 GB each,
an estimate from the saved state's size; the histories and logs are a few MB; the laptop's re-run released almost
nothing, so the droplet records stay small). If the machine has less, stop and tell Asha.

- [ ] **Step 6: launch.** The flights must outlive the shell command that started them: if your harness ends background
jobs when a command returns, run the block with its own background-run facility (in Claude Code, `run_in_background`)
instead. Both load their mesh from the state, so they can start together:

```bash
source "$HOME/skin-base-cells-flights/env.sh"
launch d050__base-off   resume "$S050" "${RESUME[@]}" --skin-base-cell 0
launch d050__base-0.2mm resume "$S050" "${RESUME[@]}" --skin-base-cell 0.2
```

A minute later confirm with `ps -Ao pid,etime,pcpu,rss,command | grep "analysis/flights/fly.py" | grep -v grep` that two
flights are running near 100 % of a core each (the `caffeinate` and `time` wrappers also show; count the `python` lines).
Within a few minutes the log of each (`$LOGS/<label>.log`) has, after its `start` line, a line beginning `resuming` that
lists the options read back. Check that it says `to 186.5 s`, `hand-down off`, `sideways heat height`, `geometry stop
off`, `frames every 0`, `save state at 182.5 s`, and `skin base cell off` for (a) and `skin base cell 0.2 mm` for (b), and
put both lines in the report. If either line is missing after 15 minutes, or differs, stop both flights (`kill <pid>` of
each `python` process, which makes them save and write what they have) and tell Asha.

## 7. Watch them

- [ ] **Step 7.** Each flight appends every history row to `$OUT/<label>.rows.csv` as it is made, and a line to
`$OUT/<label>.progress` every 20 rows (about every 0.25 s of flight). Check every few hours:

```bash
source "$HOME/skin-base-cells-flights/env.sh"
for f in "$OUT"/*.progress; do echo "$(basename "$f"): $(tail -n 1 "$f")"; done
ps -Ao pid,etime,pcpu,rss,command | grep "analysis/flights/fly.py" | grep -v grep
df -h "$HOME" | tail -n 1
```

The first history row can take more than an hour (the laptop's first macro step took 4,267 to 5,009 s on code built on
`main`), so an empty progress file in the first two hours is not a failure. Estimate each flight's finish from the flight
time still to go (to 186.5 s) and the wall seconds per flight second over the last few progress lines, and tell Asha in
one line after the first three hours, then every few hours. A process near 100 % of a core whose progress file has not
grown for two hours is not a failure; note it and keep watching. Resident memory (the `rss` column, in kilobytes) above
10 GB for one flight, or free disk below 5 GB: tell Asha.

When a flight's process has gone, read the end of its log (`tail -n 60 "$LOGS/<label>.log"`). The last lines are `exit
<code>` and the end time, and `/usr/bin/time -l`'s wall time (`real`) and `maximum resident set size` just above them.

- **`exit 0`**: the flight ended normally at 186.5 s. Its `<label>.csv` (the history), `<label>.json` (settings, results
  and provenance) and `<label>/` (the droplet files) are in `$OUT`. Check in the JSON that `settings.skin_base_cell_mm`
  is `0.0` for (a) and `0.2` for (b), and `settings.geometry_stop` is `off`.
- **exit 1, after the flight had started** (its `rows.csv` has more than a few rows): a model failure, for example a fold
  of the surface that the mesh rebuild refuses (the error reads "a corner move of up to ... degrades the mesh even after a
  rebuild"). **Do not retry and do not change settings.** The flight has written its `<label>.csv` and `<label>.json`
  from its completed steps (end reason `failed`, the error and its time in `results.failure`) and saved its last
  completed step as `<label>.end.state.pkl`. Report the error text and the time of the last row in `<label>.rows.csv`,
  and let the other flight run on.
- **exit 2 with no files**: a bad argument or a missing library; fix the environment only (Step 3), then launch that
  flight again from Step 6.
- **no `exit` line, and the process has gone** (a reboot, or the machine slept): if `$OUT/<label>.state.pkl` exists (the
  flight passed 182.5 s), resume it from there, with the options it saved, under a new label so that nothing is
  overwritten:
  `source "$HOME/skin-base-cells-flights/env.sh"; launch <label>__from182.5 resume "$OUT/<label>.state.pkl" --t-max 186.5 --frames-every 0`;
  otherwise launch the flight again from Step 6. Either way report that it was restarted, why, and from where. A
  `__from182.5` resume's history starts at 182.5 s, so bring back both labels' files. A flight that stopped on a model
  failure is not restarted (above).

## 8. Bring the results back

- [ ] **Step 8.** When both flights have ended, commit what they made (if one ends hours before the other, commit it then
and again at the end). The small files go in git; the large ones stay on the Mac Studio:

```bash
source "$HOME/skin-base-cells-flights/env.sh"
setopt nullglob 2>/dev/null || shopt -s nullglob      # zsh: setopt; bash: shopt. A glob that matches nothing then expands to nothing
cd "$REPO"
R="$CODE/reports/skin_base_cells_flights_macstudio"
mkdir -p "$R"
for f in "$OUT"/*.csv "$OUT"/*.json "$OUT"/*.progress; do cp "$f" "$R"/; done
for d in "$OUT"/*/; do l="$(basename "$d")"; mkdir -p "$R/$l"; cp "$d"/particles_summary.csv "$d"/size_distribution.csv "$R/$l"/ 2>/dev/null; done
cp "$LOGS"/*.log "$LOGS"/launches.txt "$R"/
{ for f in "$OUT"/*/particles.npz "$OUT"/*.state.pkl; do [ -f "$f" ] && { ls -l "$f"; shasum -a 256 "$f"; }; done; } > "$R/large_files.txt"
{ sysctl -n machdep.cpu.brand_string hw.perflevel0.physicalcpu hw.ncpu hw.memsize; sw_vers;
  "$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera; print(sys.version); print(numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__)";
  git rev-parse HEAD; } > "$R/environment.txt"
du -sh "$R"            # expect a few MB at most
git switch -c skin-base-cells-flights-macstudio 2>/dev/null || git switch skin-base-cells-flights-macstudio
git add "$R"
git commit -m "Skins' base cells: the two confirmation flights on the Mac Studio: <which flights have ended, and how each ended>"
git push -u origin skin-base-cells-flights-macstudio
```

Push only this branch. The `particles.npz` files, if under 50 MB each, may be copied into the results folder too
(`cp "$OUT/<label>/particles.npz" "$R/<label>/"`); the saved states (`*.state.pkl`) are never committed and
`large_files.txt` lists them with their checksums.

What comes back, per flight: `<label>.rows.csv` (the history as the flight made it, also for a flight that stopped
early) and `<label>.progress`; `<label>.csv`, the history written at the end, whose first row is the saved step's at
178.5 s; `<label>.json`; the log, with the wall time and peak memory; and the droplet summaries.

## 9. Report to Asha

- [ ] **Step 9.** After each commit, write the report in Asha's style (full sentences; every number with what it measures,
its units and why it matters; terms spelled out; no symbols in prose). In this order:

1. Which flights ended and how (at 186.5 s or a failure, with the error and its time), whether either was restarted (and
   from where), the commit flown, and the results branch and commit.
2. The two `resuming` lines of Step 6.
3. For each flight, from its `<label>.csv` (first row: the saved step at 178.5 s; last row: the end): the increase of
   `jump_faults` (the faults counted under Task 5k's rule) and of `jump_faults_bottom_cell` (the old rule, the bottom
   cell alone) from the first row to the last, and the flight time of the first row at which each rises; `skin_cells`
   (the cells of all skins) and `skin_thickness_max_mm` (the thickest skin) at the first and last rows and their
   largest values; the increase of `drawn_mass_kg` and of `absorbed_mass_kg` (the metal taken into the skins from the 3D
   model, in kilograms); the increase of `rebuilds` and of `interface_repeats`; the largest `surface_fold_deg`; and the
   increase of `released_mass_kg`.
4. The cost: each flight's wall time and peak memory (the end of its log), and the increase from the first row to the
   last of `wall_skins_s`, `wall_3d_s`, `wall_move_s` and `wall_rebuild_s` (the wall seconds spent in the skins, the 3D
   solve, the mesh moves and the rebuilds), and that the two ran side by side.
5. Anything unexpected: a restart, a failure, an environment difference from the laptop (Step 3), or a Step 4 result other
   than `2 passed`.

For orientation only, not to judge against: the laptop's re-run of this span on older code (with the old grading)
counted 30,848 faults under the old rule and 25,620 under Task 5k's. Do not interpret the results beyond the report
above: the comparison of the two flights with each other and with that re-run is Asha's, with the laptop's agent.
