# Handoff note for the refactor agent: the skyscrapers code is on main (2026-11-11)

This note is for the agent who will carry out the combined refactor of 2026-10-12 on the code that has just been
ported onto `main`. You have not seen the work that led here, so every term is explained where it first appears.
Read it first, then `Reentry Model/CLAUDE.md`, which was updated by the port and describes the layout.

## What the code is

"The skyscrapers code" is the Python package `reentry_model` (in `Reentry Model/reentry_model/`) as extended by
the Step 3 work: a model of a melting re-entering sphere of aluminium alloy AA7075 (a sphere that breaks off a
satellite, heats up, melts and sprays droplets into the wake). Its main part is the layered-skin moving-mesh model,
in which the sphere's outer metal is a set of thin columns ("skins", stacked layers) whose mesh moves with the
surface as it melts and flows. The package also keeps the older element model (a fixed mesh with a surface model
`--surface-model elements`). The Step 3 designs, plans and prototype folders came along in the same port, under
`Reentry Model/docs/superpowers/` and `Reentry Model/prototype/`.

## What landed, and where

The work is on the branch `skyscrapers-port`, which is to become `main` (it was not pushed, and nothing here merged
it). Before the port, `main` stood at commit `8ed75cf`; the port brought in the two merged lines, counted
and verified on 2026-11-11: 80 prototype-line commits and 96 rewritten build commits. The port's own first-parent
commits (below) are listed by `git log --first-parent 8ed75cf..main`. The first-parent commits, in order:

- `595928b`: the port's own tools, in `Reentry Model/analysis/skyscrapers_port/` (the path map, the checks, the
  three-way resolver; later `same_results.py`). Read its `README.md`.
- `592a77d` (merge 1): the prototype line, meaning the branch `step3-rigid-substrate-scheil` (tag
  `skyscrapers-prototype-end`, commit `16675d7`), merged into main's layout. It brought the designs, plans, specs and
  prototype folders under `Reentry Model/`, and the melting reference runs of the SESAM program (ESA's DRAMA
  re-entry module, whose output the model is verified against) into `Reentry Model/data/reference_runs/`.
- `d716a24` (merge 2): the "build repository", the separate git repository (branch `skin-spraying`, 96 commits)
  where the skyscrapers code was actually written, with every commit rewritten into main's paths. Each rewritten
  commit's message names its original in a `Build-repository-commit:` trailer, so `git log --grep` and `git blame`
  reach the original history. The original, unrewritten history is kept under the tag `skyscrapers-build-end`
  (commit `75519a5`). The merge added 47 files and changed 26; 34 files conflicted, 30 were resolved by the
  three-way resolver and 4 by hand (`tests/conftest.py` takes main's; `mesh.py` and
  `analysis/reentry_model_thermal_verification.py` take the build's code with main's output-path lines;
  `tests/test_reentry_model_data.py` keeps both sides' tests).
- `b92c4e9`: data paths in main's layout (the build's tests read `PROJECT_ROOT` and `data/`, outputs go to
  `OUTPUT_DIR`, the melting references come from `sesam_io.MELTING_REFERENCE_DIR`). These are the only code edits
  besides main's own changes reapplied; the port changed no physics, default, option or number.
- `419c1fd` and `2164b67`: documentation (`Reentry Model/CLAUDE.md`, the root `README.md`, a README in each of six
  retired prototype folders).
- `1e87c82`, `68d7ccb`, `33196a6`, `58d768a`: `same_results.py`, the tool that proves "same results" (next section).

Design of the port: `Reentry Model/docs/superpowers/specs/2026-11-11-skyscrapers-port-to-main-design.md`; its plan:
`Reentry Model/docs/superpowers/plans/2026-11-11-skyscrapers-port-to-main.md`. The path checks passed on the merged
tree: of the build's 145 files, 118 are the same in the port, 11 are the same with main's own change reapplied, 1
(the build's `Reentry Model/.gitignore`) was removed, and 15 are fix-ups (the four hand resolutions and the twelve
path repointings, with `data.py` in both lists). Only `.gitignore` and the root `README.md` differ from `main` outside
`Reentry Model/`. `earth.py`, `heating.py` and `atmosphere.py` are unchanged from the old main.

## What you are to do

Carry out the combined refactor decided on 2026-10-12. It is described in
`Reentry Model/docs/superpowers/plans/shared/README.md` (the table of six "shared" sub-plans S1 to S6, now
superseded by one refactor) and in section 14, item "0+S", of the level-set design
`Reentry Model/docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`: the Step 1 and 2
modules and Step 3's architecture-independent physics move into `reentry_model/shared/` with aliases at the old
paths; the spray physics is split from the skyscrapers bookkeeping; and `run --architecture skyscrapers|levelset` is
added. The skins are never shared. Do it on its own short branch off `main`. The plan for it is not written; the
design says to write it once the port has landed, in `docs/superpowers/plans/levelset/`. Moves of this size are
where silent changes of results happen, so use the check below after each stage.

## How to show the refactor changes nothing

The tool is `Reentry Model/analysis/skyscrapers_port/same_results.py` (documented in the README beside it). Its
`run` command runs named cases with `python -B -m reentry_model ...` inside a code folder and writes the outputs to a
folder; its `compare` command compares two such folders: every csv cell as text, every json field, every npz array
element, and any other file byte for byte, skipping only the following: the files ending `.log`, `.pkl`, `.progress`, `.vtu`, `.vtp` and `.pvd`; the json keys
`files` and `provenance`; every csv column or json key named `wall_*` or `runtime_s` (the run's elapsed time); and the
mesh and material file paths in the settings, which are compared by file name and SHA-256 checksum instead. Because
`.pkl` files are skipped, a change to a saved state is not caught by `compare`; this matters for defect 1 below, where
state saving is what is broken, so check any fix to it separately. `compare` exits 0 only on an identical result, and exits 2 on a missing or empty folder. `run`
refuses a case folder that already holds results, so give each run a new folder. The cases:

- `element`: a 0.125 s flight of a 100 mm sphere, melting on, with `--surface-model elements`.
- `skins`: the same flight with the layered skins (`--skin-preset dev --fine-size 2`).
- `resume`: the 100 mm sphere's saved state at 22.0 s, run from there to an end time of 22.6 s, about 0.6 s of flight (about 9 417 skin columns) with the
  options of the 2026-11-10 flight runbook. It needs `--state <file>`. The state is
  `prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_snapshots/d100__physics__skin.state.pkl`
  in the worktree `.claude/worktrees/subplan20` on Asha's laptop. It is 1 380 690 172 bytes and git-ignored, so it is
  not in the repository; read it only, never overwrite it.
- `nomelt`: a 20 s non-melting flight (100 mm, 7.5 km/s, 77.5 km). See defect 3 below: it fails by design of the
  default, so it is not used for "same results". (Since 2026-11-12 it runs the element model by default, the same
  flight as `nomelt_elements`.)
- `nomelt_elements`: the same with `--surface-model elements`; this one runs.

The Step 3 reference-marked tests (`tests/test_reentry_model_reference_melt.py` and two in `tests/test_reentry_model_skinbody.py`) have been collected but not run on the port.

The reference outputs are on Asha's laptop and git-ignored: the port's are in `Reentry Model/output/port_check/port/`
of the worktree `.claude/worktrees/skyscrapers-port`, and the prototype's in
`prototype/work-2026-10-28-skin-spraying/measure/port/prototype/` of the worktree `.claude/worktrees/subplan20`.
The baseline outputs in `Reentry Model/output/port_check/port` (worktree `skyscrapers-port`) name that worktree's mesh-cache and material files, whose checksums `compare` checks, and the resume case's state file lives in the `subplan20` worktree; so neither worktree may be removed until the refactor has been verified against the baseline.
The commands, from `Reentry Model/` of your worktree, with the interpreter of the conda environment `drama_env`:

```
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
"$PY" analysis/skyscrapers_port/same_results.py run --code "<your worktree>/Reentry Model" --out <a new folder> \
    --cases element,skins,resume --state <the 22.0 s state file named above>
"$PY" analysis/skyscrapers_port/same_results.py compare <that new folder> "<skyscrapers-port worktree>/Reentry Model/output/port_check/port"
```

`run` refuses (exit 2, nothing run, nothing deleted) a folder that already holds results, so use a new folder each time.
In words: run the three melting cases (`element`, `skins`, `resume`) on your refactored code into a new folder, then
`compare` it against the port's folder; the expected output is `identical: 22 files compared`, exit 0. If the port's
folder is missing, rebuild it first by running the same three cases on the port's own head (before your first
change), and check them against the prototype's if you can. Run it on mains power (`pmset -g batt | head -1` should
say `AC Power`); the three runs took 4.6 s, 21.1 s and 289.8 s of wall time on the port.

## What was measured on the port (2026-11-11)

These are the results of the port's verification (the report is in `.superpowers/sdd/2026-11-11-skyscrapers-port-to-main/task-7-report.md`
of the `subplan20` worktree, which is git-ignored). They are the baseline against which your refactor is judged.

1. Bit for bit against the prototype: the three melting cases give the identical 22 compared files (every csv cell,
   json field and npz array, plus the element case's pictures), against the prototype's code. The only difference
   found at first was `runtime_s`, an elapsed time, which the tool now ignores. This was one run of each case.
2. The port is no slower: wall times of 4.6 s (element), 23.4 s against 21.1 s (skins) and 303.1 s against 289.8 s
   (resume), prototype against port.
3. The unit tier (`PYTHONDONTWRITEBYTECODE=1 python -m pytest tests -m "not drama and not reference" -n 4 -p no:cacheprovider`,
   from `Reentry Model/`, the Python of the conda environment `drama_env`): 702 passed, 1 skipped (the FEniCSx module,
   because `dolfinx` is not installed) in about 24 minutes (1 440.83 s) on four processes. The prototype's
   count was 709; the difference of 7 is the tests of `test_fixtures.py`, which moved to `DRAMA Sphere Experiments/`.
4. Main's pinned reference tests of Step 2 (`tests/test_reentry_model_reference.py` and
   `tests/test_reentry_model_reference_thermal.py`, marker `reference`): 13 passed in 539.37 s.
5. The droplet and wake-column tests (`-m "not reference"`): 183 passed in 200.02 s. The Spheral tests: 356 passed,
   95 skipped (docker is not installed, so the container tests skip), 1 failed (defect 4).
6. Without melting, against the old main at `8ed75cf`: the trajectory columns (altitude, velocity, mass, drag,
   Knudsen number) are identical; 13 thermal columns differ slightly. The surface temperature runs up to 3.92 K (1.05
   percent) different at its maximum, the stagnation-point temperature 4.47 K (1.20 percent), the total heat absorbed
   by the sphere 0.031 percent, and the mean body temperature 0.040 K. This is the expected effect of the prototype's
   rewritten thermal solver, which came with the port, not of the port. It is one case, and with the element surface
   model, because main's command line has no `--surface-model`.

## Known defects, left unfixed

None of these is the port's doing, and none was fixed, so that the port changes no behaviour. Decide with Asha before
fixing any of them, because a fix changes behaviour and so changes what "same results" means. (Defects 1 to 3 were
fixed on 2026-11-12 with Asha's approval, on branch `skyscrapers-fixes`; each is marked below.)

1. **State saving and resuming with the element model.** A melting element-model flight cannot save its state: it
   fails with "cannot pickle module". A non-melting element-model resume raises `AttributeError` on `source_rows`.
   Both were found in sub-plan 22, Task 5k, while testing a 40-to-50-bin padding of the saved state for a non-skin
   body, and they are in the prototype's code as well. Skin flights save and resume correctly (the `resume` case).
   **Fixed on 2026-11-12** (Asha's approval; branch `skyscrapers-fixes`): `skyscrapers/body.py` imports the film module at
   module level instead of holding it on the body, and `cmd_resume` reads `source_rows` only where the body has them.
   Neither touches a calculation; `tests/test_reentry_model_restart.py` pins both (a melting and a non-melting element
   flight each resume to the last bit).
2. **Mesh cache race.** `reentry_model/mesh.py` caches each generated mesh in `Reentry Model/output/meshes/`,
   writing a temporary file with a fixed name (`.tmp.msh`) and then moving it into place with `os.replace`. Two
   processes meshing the same sphere into an empty cache can collide on that file. It broke one test in the first
   parallel unit-tier run (701 passed, 1 failed; the second run, with the cache full, gave 702 passed, 1 skipped). The
   old main has the same code, so the race is older than the port. Any parallel run on an empty cache can hit it.
   Warm the cache by running the unit tier once serially, or rerun after a failure; a fix would write the temporary
   file with a unique name. **Fixed on 2026-11-12** (Asha's approval; branch `skyscrapers-fixes`): each writer writes
   its own temporary file, named with its process id and a random suffix and still ending in `.msh`, and removes it if
   the write fails (`mesh.temporary_mesh_path`, pinned in `tests/test_reentry_model_mesh.py`).
3. **The default surface model without melting.** With `--thermal fem`, `--surface-model` defaults to `skin` even
   with `--melt off`; this became the prototype's default in sub-plan 22. The non-melting case (100 mm, 7.5 km/s,
   77.5 km altitude) then stops after about 50 s of wall time (51.6 s on the port, 48.6 s on the prototype) with
   "a corner move of up to 0.0008 m degrades the mesh even after a rebuild" (`movemesh.py`, line 473). Non-melting
   runs need `--surface-model elements`. Asha has not yet decided whether the default should change; it waits for her.
   **Decided and fixed on 2026-11-12** (Asha's decision; branch `skyscrapers-fixes`; fact 195 of the shared context).
   The cause was not melting being off but the mesh: under the 4 mm surface the default 0.7 mm fine band leaves flat
   elements at the surface, so the skins' first move (their 0.8 mm minimum taken at time zero) is refused, melting or
   not (`output/open_items/skins_nomelt-report.md`, git-ignored, of this worktree). Two changes: a thermal run that does
   not melt now runs the element model unless `--surface-model skin` is given (the skins stay the default with
   `--melt on`), so this case runs (40 macro steps, about 1 s of test time) and `nomelt` is the same run as
   `nomelt_elements`; and a skin run whose fine band leaves elements below the moving mesh's radius-ratio floor exits 2
   at once with the remedy (7.4 s here instead of exit 1 after about 50 s). No melting flight's numbers or name change.
4. **A Spheral fixture test.** `tests/test_spheral_frag_fixtures.py::test_regenerated_sphere_run_equals_committed_bitwise`
   fails because a regenerated csv ends its lines with LF where the committed fixture has CRLF (first difference at
   byte 211). It fails the same way on the old main.
5. **Two sets of facts numbered 101 to 107.** `Reentry Model/docs/superpowers/plans/melt-spraying-subplans/00-shared-context.md`
   holds main's frame-export amendment and the prototype's moving-mesh amendment, each with facts numbered 101 to 107.
   They were left unrenumbered so that references to either keep working. Cite them by their heading.

## Later work that is not yet done

- **The flight launcher.** The scripts that start long skin flights, `launch.sh` and `fly.py` in
  `Reentry Model/prototype/work-2026-10-28-skin-spraying-plan/flights/`, still run the retired planning copy of the
  code (`.../work-2026-10-28-skin-spraying-plan/code/`), not the package. That folder is kept as a record and must not
  be edited. A launcher that runs the package is later work. Flights longer than about two hours run on the Mac Studio
  by a runbook (`docs/superpowers/plans/2026-11-10-skin-fixes-test-flights-macstudio.md`). **Done on 2026-11-12**
  (branch `skyscrapers-fixes`): `Reentry Model/analysis/flights/` holds the same `fly.py` and `launch.sh`, running the
  package (see its README).
- **Branches made before the port.** These merge later, each on its own, and will meet the moved prototype paths
  (prototype folders now live under `Reentry Model/prototype/`) and the rewritten history: `levelset-two-interface`,
  `levelset-stage0`, `wake-residence`, `segmented-trail`, `janaf-vapour-pressure`, `we-critical-308` and
  `skin-fixes-flights-macstudio` (which exists only as the remote branch
  `origin/skin-fixes-flights-macstudio`; its results were pushed on 2026-11-11 at `39e8ded`). Do not merge them as part of the refactor. `levelset-stage0` is superseded by the
  refactor itself (design section 14). Resolve conflicts one branch at a time with the path map in
  `analysis/skyscrapers_port/path_map.py`. Merging such a branch can leave files in the old flat layout without git raising a conflict (dry runs on 2026-11-11: `levelset-two-interface` would leave 62 files at the repository root, under `prototype/levelset`, `reentry_model/levelset`, `tests/levelset` and `analysis/levelset`, and bring back a root `docs/model_assumptions.md` through a modify/delete conflict, with 6 conflicts in all; `origin/skin-fixes-flights-macstudio` would put 20 result files under a root `prototype/work-2026-10-28-skin-spraying-plan/flights/macstudio-fixes/` with 1 conflict). So after such a merge, before committing, run `"$PY" "Reentry Model/analysis/skyscrapers_port/checks.py" outside "$(git write-tree)" --base HEAD --allow .gitignore --allow README.md` from the repository top, and move each path it reports under `Reentry Model/` (the place `path_map.prototype_line_path` gives). (`--base HEAD`, added 2026-11-12, compares against the commit being merged into; without it the check compares against `MAIN0`, `main` before the port, and so also reports every later change of `main` outside `Reentry Model/`.)
- **Do not touch** `DRAMA Sphere Experiments/`, `Spheral/`, `reentry_model/droplet_reentry/` and
  `reentry_model/wake_column/` without Asha's say: they were checked unchanged against the old main.

## What to say to Asha, and when

Report to Asha with full sentences and every number explained. Stop and ask before changing any result, default or
option, since the proof of a refactor is that nothing changes.
