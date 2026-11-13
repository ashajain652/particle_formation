# The flight launcher

`fly.py` and `launch.sh` start a long flight of the model (a sub-plan 22 Task 5 flight, which takes days) so that it keeps
its history as it goes. They run the package in this repository's `Reentry Model/` folder. They replace, for flights from
`main`, the launcher in `prototype/work-2026-10-28-skin-spraying-plan/flights/`, which runs the retired planning copy of
the code in that folder's `code/` and stays there unedited as a record (ported 2026-11-12).

`fly.py` runs the command line (`reentry_model run`, or `reentry_model resume`) in its own process with a hook on
`reentry_model.skyscrapers.coupled.CoupledRun.row` that only reads each history row: it appends every row to
`<outdir>/<name>.rows.csv` as it is made (every digit) and a progress line every 20 rows to `<outdir>/<name>.progress`.
The exit code is the command line's, or 1 after the traceback on any other exception.
`launch.sh` runs `fly.py` from `Reentry Model/` with no bytecode written, under `caffeinate -i` (no idle sleep) and
`/usr/bin/time -l` (wall time and peak memory at the end of the log), and logs the start, the end, the exit code and
the power source (`pmset -g batt`) in `<logdir>/<label>.log`. `$PY` chooses the interpreter (default: the laptop's
`drama_env`).

The model writes its history and results when a flight ends, however it ends (since 2026-11-12). A flight that fails
(any error raised in a macro step: a fold of the surface, a failed solve, a moving-mesh error) or is stopped by hand
(SIGTERM sent to its Python process, as `kill <pid>` sends it, or Control-C) leaves `<outdir>/<name>.end.state.pkl`, its
state at the end of its last completed macro step, which `resume` continues by repeating the failing step, and
`<name>.csv`, `<name>.json` and the droplet files under `<name>/`, written from the completed steps with end reason
`failed` or `stopped` (the JSON's `results.failure` gives the error and the time of that step). It exits with code 1 for
a failure (after the traceback), 143 for SIGTERM and 130 for Control-C. The rows file can hold rows that `<name>.csv`
does not: the hook keeps a step's row as soon as it is made, before that step's frame, its demise and geometry tests,
its `--save-state-at` save and the copy of its state, so a failure or a stop in any of those leaves that step's row in
the rows file while `<name>.csv` ends at the step before. A process killed without warning (`kill -9`, a crash of Python
itself, a power cut) saves nothing: then the rows file holds the rows made so far and `<outdir>/<name>/source_rows` the
droplet records.

The two command forms, from `Reentry Model/` (relative folders are taken from where the command is run):

```bash
analysis/flights/launch.sh <outdir> <logdir> <label> [reentry_model run arguments...]      # a new flight
analysis/flights/launch.sh <outdir> <logdir> <label> resume <state> [resume options...]    # a saved flight continued
"$PY" -B analysis/flights/fly.py <outdir> <name> [run arguments...]                        # the same without the wrapper
"$PY" -B analysis/flights/fly.py <outdir> <name> resume <state> [resume options...]
```

`--outdir` and `--name` are added by `fly.py` (the label is the run's name); do not give them. Put the outputs under
`Reentry Model/output/` (for example `output/flights/<date>-<what>/`), which is git-ignored. A resumed flight's rows file
holds the rows made after the saved step; the saved step's own row is the first row of its history. Flights longer than
about two hours run on the Mac Studio by a runbook (`docs/superpowers/plans/2026-11-10-skin-fixes-test-flights-macstudio.md`
is the last one, written for the prototype's launcher). The tests are `tests/test_reentry_model_flight_launcher.py`.
