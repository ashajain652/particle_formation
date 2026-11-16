# A geometry stop for skin flights, and a saved state at every end of a flight: design

**Status: drafted 2026-11-12 from Asha's decisions of the same day, for her approval.** Open item 4 of
`docs/superpowers/plans/2026-11-12-skyscrapers-open-items.md`, with the saved state Asha added to it.

## Why

A skin flight runs on, with nothing to stop it, after its columns no longer describe the metal. The 50 mm comparison
flights of 2026-11-09 to 2026-11-10 left valid geometry at about 185.5-186 s and ran on for more than a day of laptop
time; the Mac Studio flight `d050__fixes__slipband-on` (2026-11-11) did the same: by 185.5 s its column tops enclosed
92 cm3 against 65 cm3 of metal, by 187.5 s 133 cm3, and a liquid layer appeared on its rear surface, 30 mm behind a core
below the liquidus, which nothing in the model can have heated there (open item 3, report
`output/open_items/rear_heating-report.md` in the skyscrapers-fixes worktree). Results after that point cannot be used,
and each such flight costs days of a machine.

Separately, a flight that fails (exit 1: the mesh-quality stop, the part-count runaway, a Newton non-convergence, any
unexpected error) or is stopped by hand leaves no saved state and no results files today: only the launcher's
`<name>.rows.csv` survives. A state exists only where `--save-state-at` asked for one. Asha, 2026-11-12: "the code
should always produce a saved state when it terminates."

## Goals

1. A skin flight ends itself, with a message, at whichever comes first of two tests (Asha, 2026-11-12: both tests):
   the core's front has passed the sphere's original centre, or the column tops enclose more than 1.25 times the
   volume of the metal.
2. Every flight of `run` and `resume` saves a state when it ends, however it ends, except where the process is killed
   without warning (below).
3. A flight that fails saves the state at the end of its last completed macro step (Asha, 2026-11-12), so that a resume
   from it repeats the failing step exactly.

## Non-goals

- Fixing the geometry itself (columns that fan out, folding, overhangs): the stop only ends the flight before it.
- The element model's geometry: it has no columns and gets no geometry stop. It does get the saved state at every end.
- A periodic checkpoint on disk. A forced kill (`kill -9`), a crash of the Python interpreter itself, or a power cut
  gives a process no chance to save, so those three ends save nothing; only a checkpoint written every so often would
  cover them, and it is not part of this design.

## Approaches considered

- **A body-level end condition checked after each macro step, as demise is (chosen).** `coupled.CoupledRun.run`
  already ends a flight normally when `body.demised()` is true, setting `sim.end_reason = "demise"` and writing the
  results. The geometry stop sits beside it: a normal end, exit 0, end reason `geometry`.
- **Raising an error (exit 1), as the mesh-quality stop does.** Rejected: it marks a valid history as a failure, and the
  Mac Studio runbooks treat every exit other than 0 as an error to report rather than an end.
- **A watcher outside the model that kills the flight on reading its rows.** Rejected: the rows carry neither quantity
  today, the tests could not exercise it, and the kill would lose the state.

For the saved state at a failure, three were considered: the last completed step (chosen by Asha), the state as it stood
at the failure (a step half done, which cannot be resumed), and both. The last completed step suffices for diagnosis
because a resume repeats the flight to the last bit (`tests/test_reentry_model_restart.py`), so resuming it re-runs the
failing step, which can then be instrumented.

## Design

### The two tests (skin flights)

Both are measured after every macro step, on the surface of that step's end, and reported in every history row.

- **The core's front**, `core_front_mm`: the largest distance of the 3D core's surface corners (the skins' base
  corners, `columns.points`) ahead of the sphere's original centre along the flight direction, in millimetres. The mesh
  is generated around the origin (`mesh.py`, `addSphere(0, 0, 0, radius)`), so the original centre is the origin. The
  test fires when the front is behind the centre, `core_front_mm < 0`. On `d050__fixes__slipband-on` that is about
  188.5 s (read from its cross-section figure).
- **The volume the tops enclose**, `tops_volume_ratio`: the volume enclosed by the top surface (`SkinBody.surface`, the
  gas boundary, by the divergence theorem) over the volume of the metal the model holds, the 3D core's active volume
  (`mesh.volume()`) plus the skins' cell volumes (each cell's mass over its own density, so liquid counts at the liquid
  density). A geometry that describes its metal gives a ratio of about 1 whether or not it has melted; the excess is
  space between the columns that holds no metal. The test fires when the ratio exceeds `GEOMETRY_VOLUME_RATIO = 1.25`.
  The cross-section figure's version of the ratio, which divided every cell's mass by the solid density, was 1.03 at
  179 s, 1.17 at 183 s, 1.42 at 185.5 s and 2.0 at 187.5 s on that flight, so the stop would fire at about 184 s; the
  model's own volumes count liquid at its lower density, so the real stop comes a little later. Not yet measured on a
  flight.

When either fires, `sim.end_reason = "geometry"`, the flight writes its results as at any normal end, exits 0, and prints
one line naming the test, the time and the value, for example `geometry stop at 184.03 s: the column tops enclose 1.26
times the metal's volume (limit 1.25)` or `geometry stop at 188.51 s: the core's front is 0.4 mm behind the original
centre`. The results JSON records `geometry_stop: {"test": "volume" | "centre", "time_s": ..., "value": ...}`.

### Saving a state at every end

- **A normal end** (geometry, demise, ground, `t_max`, escape): after the last step the run saves its whole state, as
  `_save` does today, to `<outdir>/<name>.end.state.pkl`. The name differs from `--save-state-at`'s
  `<name>.state.pkl`, so the mid-flight state the runbooks resume from is never overwritten.
- **A failure** (any exception raised while advancing a macro step): at the end of every completed step the run keeps an
  in-memory copy of its state, pickled as `restart.save` pickles it (`restart._Pickler`, protocol 5). On a failure it
  writes that copy to `<name>.end.state.pkl`, writes the results files from the completed steps' rows with end reason
  `failed` and the error message, prints the error as today, and exits 1. A new copy replaces the old only once it is
  complete, so a failure during the copy leaves the previous step's.
- **A stop by hand** (`SIGTERM`, as `kill` sends, or `SIGINT`, Control-C): caught; treated as a failure with end reason
  `stopped`, writing the last completed step's state and the results files, then exiting with 143 or 130 as the shell
  would report.
- **Every flight**: every `run` and `resume`, both surface models and the non-melting runs, with no option to turn it
  off (Asha: always). `settings.save_meta` (the command line's arguments, which a resume needs) is therefore always set.
- **A save that fails** is recorded (`state_save_error` in the JSON) and printed, and never fails the flight, as today.
- **Cost**, measured 2026-11-12 on the 50 mm state saved at 170.5 s (291 MB): 0.1 s to pickle into memory, 0.2 s to
  disk. The 50 mm macro steps measured on the laptop took about 9 s each (the smoke check: 10 steps in 87 s) to 122 s
  (one 0.5 s step from 170.5 s, 2026-11-12), and late steps of the 2026-11-09 flights minutes each, so the copy adds
  about 1 % of the time at most. It holds a second copy of the state in memory: about 0.3 GB on the 50 mm sphere, about 1.4 GB on the 100 mm
  sphere (the size of its saved states; its pickling time was not measured).

### Options and names

- `--geometry-stop on|off` on `run` and `resume`, skin flights only (refused, exit 2, with `--surface-model elements`),
  default `on` for both, including a resume of a state saved before this change. This departs from the rule that an old
  state resumes as it was flown; the stop changes nothing until it fires, and the flights it is for are resumed from old
  states. A flight with it off is named `_geomstop-off`.
- The threshold 1.25 is a module constant, not an option.

## The controller's calls (2026-11-12, for Asha to check)

- The measure of the metal's volume (the model's own volumes, not mass over the solid density). The threshold 1.25 is
  Asha's (2026-11-12, "1.25 is fine", after the comparison: 1.15 would stop the slip-band-on flight at about 182.5 s,
  1.25 at about 184 s, 1.5 at about 186 s, read between the figure's frames).
- The geometry stop is a normal end (exit 0); a failure and a stop by hand write the results files too (they write none
  today), with end reasons `failed` and `stopped`.
- One end-state name, `<name>.end.state.pkl`, for every kind of end; the JSON says which kind.
- The stop is on by default on resume, even for old states; the saved state at every end has no off switch.
- A signed volume is used for the tops' enclosed volume; where the top surface folds over itself the folded part counts
  with its sign, so the ratio understates the gaps there. Accepted, since the stop is meant to fire before folding.

## Left to the builder

- Where in `CoupledRun.run` the copy is made and the failure caught, so that the rows of the completed steps and the
  history are built as today and existing exit codes and messages are kept.
- The signal handlers' installation and removal, so that tests and library callers are not affected.
- Tests: each test firing on a small skin flight (a prescribed geometry or a lowered threshold), the end state of each
  kind of end resuming to the last bit, a failure injected into a step writing the previous step's state, and the
  message and JSON fields.
- `README.md`, `docs/model_assumptions.md` and `CLAUDE.md`, where they describe how a flight ends and its exit codes.
- Other branches touching `cli.py` and `coupled.py` (`kinked-dispersion`, `near-liquidus-viscosity`,
  `levelset-stage3`): whichever merges second adapts.
