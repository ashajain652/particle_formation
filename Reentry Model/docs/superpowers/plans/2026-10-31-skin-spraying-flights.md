# Sub-plan 22 Task 5: the spraying skin flights on the Mac Studio — runbook

> **For the agent running this.** You are on Asha's Mac Studio. Your job is to run four already-written melting
> flights, two of which take days to weeks, watch them, gather their results, and report. **You write no model code.**
> If anything in this runbook does not work as written, stop and tell Asha what happened rather than changing code or
> settings to make it work. Read the whole runbook once before running anything.

**Asked by:** Asha, 2026-10-31: flights expected to take more than about two hours run on the Mac Studio, side by side,
one core each; after the cost probe, "Run as they are" (the base flights unchanged, at the cost measured).

**Design of record:** `docs/superpowers/specs/2026-10-11-layered-skin-moving-mesh-design.md` §11 (the flights) and §10.2
(their outputs); sub-plan 22 Task 5 in `docs/superpowers/plans/melt-spraying-subplans/22-skin-spraying-flight.md`
(generated from `prototype/work-2026-10-28-skin-spraying-plan/part_22.md`).

## 1. What this is and why it matters

The repository models what happens to the solid aluminium (AA7075) fragments of a re-entering satellite. Step 3 of the
model melts them: the surface melts, the molten metal flows over the body and sprays off as droplets, whose sizes and
numbers are the thesis's result. The new surface model, **layered skins**, puts a column of very thin cells (10
micrometres at the top) on every corner of the surface mesh of a 3D heat model whose mesh moves as the surface recedes;
in each column the molten and part-molten metal (the slurry) flows, and droplets leave from the top cells by Girin's
spraying model. Sub-plans 20-22 built it. These are its first whole melting flights.

You run four flights of the two reference cases, both at 7.5 km/s, flight-path angle -0.959331 degrees, US76
atmosphere, physics-mode heating, the AA7075 alloy on Scheil's melting curve (the default with `--melt on`):

| Label | Sphere | Surface model | Ends | Macro steps (estimated) | Wall time (estimated) |
|---|---|---|---|---|---|
| `d100__physics__skin` | 100 mm from 77.5 km | skins | at 120 s of flight, saving its state at 88 s | about 44 of 0.5 s, then about 7 840 of 12.5 ms from the first liquid (about 22 s) | **one to three weeks** |
| `d050__physics__skin` | 50 mm from 115 km | skins | demise or the ground | about 340 of 0.5 s, then 12.5 ms steps from the first liquid (about 170 s); the element model demises at about 200 s | half a day to two days |
| `d100__physics__elements` | 100 mm from 77.5 km | elements | at 120 s | 0.5 s, then 12.5 ms below a Knudsen number of 0.01 (about 49.5 s) | a few hours |
| `d050__physics__elements` | 50 mm from 115 km | elements | demise or the ground | 0.5 s throughout | under an hour |

The two element-model flights are the comparison: the model the skins replace, on the same machine and code.

**The wall times are extrapolations.** The cost probe of 2026-10-31 (one second of the 100 mm sphere's melting flight
on the production mesh, on Asha's laptop, alone) took 74 minutes. A 12.5 ms step cost 15 seconds at its start and 89
seconds at its end, rising with the depth of the slowly creeping slurry, which deepens as the sphere's interior warms
towards the rigid point (895 K). The 100 mm flight has about 7 840 such steps; at 89 seconds each that is about eight
days alone, more if the slurry keeps deepening, and the level-3 flights ran about 1.5 times slower a step side by side
on this machine than alone on the laptop. **Expect the 100 mm skin flight's steps to get slower through the flight**, so
an early estimate of its finish will be too early. Memory measured on the probe: 5.4 GB at most; disk: the droplet
records are written to disk as they come, about 1.5 GB for a 100 mm flight measured at the probe's rate (up to 15 GB
estimated in Task 4), and copied into the final file at the end, so allow twice that. There is **no checkpoint**: a
flight interrupted by a reboot or a crash starts again from the beginning.

Each flight runs through `prototype/work-2026-10-28-skin-spraying-plan/flights/fly.py`, which runs the command line
unchanged and also writes every history row to `<label>.rows.csv` as it is made and a progress line every 20 macro
steps to `<label>.progress` (checked on the laptop: the smoke check of Step 4 through it and through the command line
alone give the same history, results and droplet files to the last bit). The model itself writes its history only
when a flight ends, so a flight that stops early
keeps its history to its last step in `rows.csv` and its droplet records in `<label>/source_rows/`. **Never delete
either.**

## 2. Rules

- Do not edit anything under `prototype/`, `reentry_model/`, `analysis/`, `tests/` or `docs/`, except creating the
  results directory of Step 9. The code is verified and its plan is generated from it byte for byte; a hand edit breaks
  that. If you think code is wrong, report it.
- Do not change any flight setting (macro step, mesh, cells, heating, seed) and do not re-run a flight with different
  settings. A flight that stops is a result to report, not an error to fix (Asha, 2026-10-31: a flight that stops on a
  fold of the surface keeps and judges what it made up to then).
- Do not change system settings (sleep, updates, energy). Asha keeps the machine from restarting while the flights run.
- Do not touch DRAMA, pyDRAMA, `sphere_reentry.py` or `sphere_sweep.py`; nothing here needs them.
- Run nothing else heavy on the machine while the flights run.
- Git: never use bare `git stash`. Push only the results branch named in Step 9, never `main` or
  `step3-rigid-substrate-scheil`. Do not merge anything. Never commit `particles.npz` or `source_rows/` (too large).
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
REPO="$(git rev-parse --show-toplevel)"
P="$REPO/prototype/work-2026-10-28-skin-spraying-plan"
CODE="$P/code"
FLIGHTS="$P/flights"
ls "$FLIGHTS/fly.py" "$FLIGHTS/launch.sh"          # both must exist
git log --oneline -1                                # note this commit for the report
```

`$CODE` is a complete, runnable copy of the model package: the tracked copy from which sub-plan 22 is generated,
byte-identical to the laptop's build repository at `dfe3b66` (sub-plan 22 Task 5a, which adds `--save-state-at`; a
flight not given it is the same to the last bit as at `caad816`, Task 4b, from which three of the flights may already
have been launched). Every command below runs from `$CODE` unless it says otherwise.

## 4. The Python environment

- [ ] **Step 2.** Find the Mac Studio's `drama_env` (it was used for level 3 on 2026-10-16):

```bash
conda env list
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python      # adjust if the Mac Studio's path differs
export PY
"$PY" -m pip install pytest-xdist==3.8.0                    # the parallel test runner; brings execnet
"$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera, matplotlib, meshio, pytest, xdist; print(sys.version.split()[0], numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__, matplotlib.__version__, pytest.__version__, xdist.__version__)"
```

On the laptop this prints `3.12.14 2.5.3 1.18.1 12.0.2 4.15.2 5.3.0 3.2.0 3.11.1 9.1.1 3.8.0` (Python, numpy, scipy,
scikit-fem, gmsh, pyamg, cantera, matplotlib, pytest, pytest-xdist). Installing pytest-xdist is the only environment
change asked for; if an import fails, repair the environment as the level-3 runbook's Step 2 does
(`docs/superpowers/plans/2026-10-16-skin-level3-flights.md`) and say so in the report. A different **gmsh** or
**scikit-fem** version is worth noting in the report, because the meshes and the assembly come from them.

## 5. Check the machine before spending weeks on it

- [ ] **Step 3: the unit tier**, across six processes:

```bash
cd "$CODE"
PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest -m "not drama and not reference" -q -p no:cacheprovider -n 6 \
    --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py
```

Expected: **532 passed, 1 skipped** (the skip is the FEniCSx module: `dolfinx` is not in `drama_env`; a serial run
also reports 18 deselected, which the parallel runner does not print).
On the laptop this took 19.5 minutes serially and 7 minutes across four processes. Any failure: stop and report it with
its output.

- [ ] **Step 4: the smoke check.** The first 0.1 s of the cost probe through the flight launcher: the 100 mm sphere
from 65 km at 7.24 km/s with its body at 893 K, on the production mesh, eight macro steps. It also builds and caches the
100 mm skin mesh.

```bash
SMOKE="$CODE/reentry_model_output/task5_smoke"
"$FLIGHTS/launch.sh" "$SMOKE" "$SMOKE" probe_t01 --diameter 100 --velocity 7.24 --altitude 65 \
    --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --heating physics --melt on --temperature 893 \
    --dt 0.025 --t-max 0.1 --frames-every 0
tail -n 30 "$SMOKE/probe_t01.log"
"$PY" -c "import json,sys; r=json.load(open(sys.argv[1]))['results']; print(r['end_reason'], r['n_macro_steps'], r['skin_columns'], r['sprayed_mass_kg'], r['energy_balance_residual'])" "$SMOKE/probe_t01.json"
```

Expected, as measured on the laptop on 2026-10-31: the log ends `exit 0`; the run ends by `t_max` after **7** macro
steps on **9 417** skin columns, having sprayed **0.0013978** kg (about 1.4 grams), with an energy-balance residual of
about 1e-10 (it must be below 1e-9); `/usr/bin/time -l` reports a maximum resident set size of about 4.2e9 bytes. On
the laptop, with its mesh already built, it took 2.5 minutes; building the mesh adds a few. `$SMOKE` then holds
`probe_t01.csv`, `.json`, `.rows.csv` (8 rows), `.progress` (one line) and `probe_t01/` with `particles.npz`,
`particles_summary.csv` and `size_distribution.csv`. A different number of columns means a different mesh, and a
sprayed mass that differs in its first three digits means the machines would not fly the same flights: either way,
stop and report. A difference only in the fourth digit or beyond: note it in the report and go on.

## 6. Run the four flights

- [ ] **Step 5: room.** `sysctl -n hw.perflevel0.physicalcpu hw.memsize` and `df -h "$CODE"`. The four flights need
four free cores, about 25 GB of memory between them, and at least 40 GB of free disk. If the machine has less, stop
and tell Asha.

> **The 100 mm skin flight saves its state (Asha, 2026-11-01).** It was held until it could; it now can (sub-plan 22
> Task 5a). It launches last, in Step 7b, with `--save-state-at 88`: at the end of its first macro step at or after 88 s
> of flight -- the start of the window around peak heating (93.5 s) from which the convergence repeats will resume --
> it writes its whole state, about 1.4 GB, to `$OUT/d100__physics__skin.state.pkl` and flies on unchanged. **Keep that
> file and never commit it.** If the save fails the flight goes on, and its JSON's `results.state_save_error` says why:
> report it. If you launched the other three flights from an earlier commit, `git pull --ff-only` now, which is safe
> while they run (each loaded its code when it started), and go to Step 7b.

- [ ] **Step 6: launch the 50 mm skin flight.**

```bash
OUT="$CODE/reentry_model_output/task5_flights"
LOGS="$CODE/reentry_model_output/task5_logs"
mkdir -p "$OUT" "$LOGS"
COMMON=(--velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --heating physics --melt on)
launch() {   # label, then the run arguments
    nohup "$FLIGHTS/launch.sh" "$OUT" "$LOGS" "$@" > /dev/null 2>&1 &
    echo "$1 pid $! started $(date '+%F %T')" | tee -a "$LOGS/launches.txt"
}
launch d050__physics__skin --diameter 50 --altitude 115 "${COMMON[@]}"
```

The flights must outlive the shell command that started them: if your harness ends background jobs when a command
returns, start each `launch` line with its own background-run facility (in Claude Code, `run_in_background`) instead.
A minute later confirm with `ps -Ao pid,etime,pcpu,rss,command | grep "flights/fly.py" | grep -v grep` that the flight
is running near 100 % of a core.

- [ ] **Step 7: launch the two element flights** once the 50 mm skin flight has a first line in its progress file
(`ls "$OUT"/*.progress`; its mesh is then built, so no two processes build one mesh file at once):

```bash
launch d100__physics__elements --diameter 100 --altitude 77.500133 "${COMMON[@]}" --surface-model elements --t-max 120
launch d050__physics__elements --diameter 50 --altitude 115 "${COMMON[@]}" --surface-model elements
```

- [ ] **Step 7b: launch the 100 mm skin flight, saving its state at 88 s**, once the 100 mm element flight has a first
line in its progress file. Check `git log --oneline -1` shows this runbook's commit or later and
`"$PY" -B -m reentry_model run --help | grep -c save-state-at` prints 2 (run from `$CODE`). If the shell of Steps 6-7 is gone, set `OUT`,
`LOGS`, `COMMON` and `launch` again as in Step 6 first.

```bash
launch d100__physics__skin --diameter 100 --altitude 77.500133 "${COMMON[@]}" --t-max 120 --save-state-at 88
```

The 100 mm skin mesh was cached by the smoke check of Step 4. The save comes days into the flight; its progress lines
pass 88 s some way before it.

## 7. Watch them

- [ ] **Step 8.** Every flight appends a line to `$OUT/<label>.progress` every 20 macro steps:

```
step    840  t   32.500 s  h  74.012 km  V 7.4950 km/s  m  1.4721 kg  T  520.3 K  T_surface_max  912.0 K  released  0.0123 kg  flowing depth  0.250 mm  wall     12345 s
```

`T` is the body's mean temperature, `T_surface_max` the hottest surface point, `released` the mass sprayed so far and
`flowing depth` the mean depth of the moving metal (skins only; a field the history lacks is left out of the line).
Check every few hours:

```bash
for f in "$OUT"/*.progress; do echo "$(basename "$f"): $(tail -n 1 "$f")"; done
ps -Ao pid,etime,pcpu,rss,command | grep "flights/fly.py" | grep -v grep
df -h "$CODE" | tail -n 1
```

Estimate each flight's finish from the flight time still to go and the wall seconds per flight second over the last few
lines (120 s of flight for the 100 mm flights; the 50 mm flights end when the body has gone). Give Asha a one-line
estimate after the first day, and again each day for the 100 mm skin flight, saying that its steps are expected to keep
getting slower. A process near 100 % of a core whose progress file has not grown for several hours is not a failure at
12.5 ms steps late in the flight (20 steps can take an hour); note it with the time and keep watching. Resident memory
(the `rss` column, in kilobytes) above 20 GB for one flight, or free disk below 15 GB: tell Asha.

When a flight's process has gone, read the end of its log (`tail -n 60 "$LOGS/<label>.log"`). The last lines are
`exit <code>` and the end time.

- **`exit 0`**: the flight ended normally (by `t_max`, demise or the ground); its `<label>.csv`, `<label>.json` and
  `<label>/` (the particle files) are in `$OUT`.
- **any other exit code, after the flight had started** (its `rows.csv` has more than a few rows): a model failure,
  for example a fold of the surface that the mesh rebuild refuses, or a solve that did not converge. **Do not retry and
  do not change settings.** Keep everything. Report the error text (the `ERROR:` line or the traceback), the time and
  altitude of the last row in `<label>.rows.csv`, its largest `surface_fold_deg` (the fold angle of the surface), and the
  size of `<label>/source_rows/` (`du -sh`).
- **exit 2 before the flight started** (no `rows.csv`, or only its first row): a bad argument or a missing library;
  fix the environment only (Step 2), then launch that flight again.
- **no `exit` line, and the process has gone** (a reboot, or the machine slept): launch that flight again from the
  beginning, and report that it was restarted and why.

## 8. Bring the results back

- [ ] **Step 9.** When the three shorter flights have ended (expected within two days), commit what they made; when
the 100 mm skin flight ends, commit again. Each time:

```bash
cd "$REPO"
R="$P/flights/macstudio"
mkdir -p "$R"
for f in "$OUT"/*.csv "$OUT"/*.json "$OUT"/*.progress; do cp "$f" "$R"/; done
for d in "$OUT"/*/; do l="$(basename "$d")"; mkdir -p "$R/$l"; cp "$d"/*.csv "$R/$l"/ 2>/dev/null; done
cp "$LOGS"/*.log "$LOGS"/launches.txt "$R"/
{ for f in "$OUT"/*/particles.npz "$OUT"/*.state.pkl; do [ -f "$f" ] && { ls -l "$f"; shasum -a 256 "$f"; }; done; du -sh "$OUT"/*/source_rows 2>/dev/null; } > "$R/large_files.txt"
{ sysctl -n machdep.cpu.brand_string hw.perflevel0.physicalcpu hw.ncpu hw.memsize; sw_vers;
  "$PY" -c "import sys, numpy, scipy, skfem, gmsh, pyamg, cantera; print(sys.version); print(numpy.__version__, scipy.__version__, skfem.__version__, gmsh.__version__, pyamg.__version__, cantera.__version__)";
  git rev-parse HEAD; } > "$R/environment.txt"
du -sh "$R"            # expect tens of MB at most; particles.npz and source_rows/ stay on the Mac Studio
git switch -c skin-spraying-flights-macstudio 2>/dev/null || git switch skin-spraying-flights-macstudio
git add "$R"
git commit -m "Sub-plan 22 Task 5 on the Mac Studio: <which flights have ended, and how each ended>"
git push -u origin skin-spraying-flights-macstudio
```

Push only this branch. `particles.npz` (the droplet records, about 1.5 GB for a 100 mm flight) stays on the Mac Studio;
`large_files.txt` records where it is and its checksum, and Asha fetches it.

- [ ] **Step 10: beside SESAM.** For each flight that ended normally, the comparison with SESAM's melting reference
(cheap; if it fails, note the error and go on, as the flights' results do not depend on it):

```bash
cd "$CODE"
REF="$CODE/data/reference_runs"
for l in d100__physics__skin d100__physics__elements; do
    "$PY" -B -m reentry_model compare --model "$OUT/$l.csv" --outdir "$OUT/compare_$l" --quiet \
        --reference "$REF/sphere_d100.00mm_T0300.0K_v07.50000kms_h077.500km_nowind.csv"; done
for l in d050__physics__skin d050__physics__elements; do
    "$PY" -B -m reentry_model compare --model "$OUT/$l.csv" --outdir "$OUT/compare_$l" --quiet \
        --reference "$REF/sphere_d050.00mm_T0300.0K_v07.50000kms_h115.000km_nowind.csv"; done
```

Copy each `compare_<label>/*_vs_reference.json` into `$R` with the next commit.

## 9. Report to Asha

- [ ] **Step 11.** After each commit, write the report in Asha's style (full sentences; every number with what it
measures, its units and why it matters; terms spelled out; no symbols in prose). In this order:

1. Which flights ended and how (by `t_max`, demise, the ground, or a failure), and whether the results are pushed
   (branch and commit).
2. For each flight that ended: the mass sprayed and the mass left, from its JSON's `results`; the number of droplets
   and their median radius by number and by mass; and the same for the element flight of the same sphere beside it.
3. For a flight that stopped early: the error, when and where it stopped (time and altitude), its largest fold angle,
   and what its `rows.csv` and `source_rows/` hold.
4. For the 100 mm skin flight: whether its state was saved (`results.state_saved_at_s`, the flight time of the save,
   or `results.state_save_error`), and the state file's size and checksum.
5. The cost: each flight's wall time and peak memory (the end of its log), how the wall time per flight second changed
   through the flight (from the progress file), and that the flights ran side by side.
6. Anything unexpected: a restart, a stalled progress file, an environment difference from the laptop (Step 2), or a
   smoke check that differed from Step 4's expected values.

Do not interpret the results beyond this: the convergence checks, the comparison with the element model and the
acceptance of the flights are Asha's, with the laptop's agent.
