# The flight launcher

`fly.py` and `launch.sh` start a long flight of the model (a sub-plan 22 Task 5 flight, which takes days) so that it keeps
its history as it goes. They run the package in this repository's `Reentry Model/` folder. They replace, for flights from
`main`, the launcher in `prototype/work-2026-10-28-skin-spraying-plan/flights/`, which runs the retired planning copy of
the code in that folder's `code/` and stays there unedited as a record (ported 2026-11-12).

`fly.py` runs the command line (`reentry_model run`, or `reentry_model resume`) in its own process with a hook on
`reentry_model.skyscrapers.coupled.CoupledRun.row` that only reads each history row: it appends every row to
`<outdir>/<name>.rows.csv` as it is made (every digit) and a progress line every 20 rows to `<outdir>/<name>.progress`.
The model itself writes its history and results only when a flight ends, so a flight that stops early (a fold of the
surface, a failed solve, an interrupt) keeps its rows to its last completed macro step; its droplet records are in
`<outdir>/<name>/source_rows`. The exit code is the command line's, or 1 after the traceback on any other exception.
`launch.sh` runs `fly.py` from `Reentry Model/` with no bytecode written, under `caffeinate -i` (no idle sleep) and
`/usr/bin/time -l` (wall time and peak memory at the end of the log), and logs the start, the end, the exit code and
the power source (`pmset -g batt`) in `<logdir>/<label>.log`. `$PY` chooses the interpreter (default: the laptop's
`drama_env`).

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
