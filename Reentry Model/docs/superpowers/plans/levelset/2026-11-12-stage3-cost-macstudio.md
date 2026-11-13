# Level-set stage 3 Task 11, Step 4: the 50 mm sphere's cost on the Mac Studio — runbook

> **For the agent running this.** You are on Asha's Mac Studio. Your job is to run two already-written cost
> measurements of the level-set model's skinned body on the 50 mm sphere, each expected to take one to three hours, gather
> their results, push them on a results branch, and report. **You write no model code.** If anything in this runbook
> does not work as written, stop and tell Asha what happened rather than changing code or settings to make it work.
> Read the whole runbook once before running anything.

**Asked by:** the plan of level-set stage 3 (`docs/superpowers/plans/levelset/03-stage3-skins.md`, Task 11, Step 4):
"On the 50 mm sphere at a melting state (dev and production cells): time and memory per macro step split into loads,
skins (trial, real, repeats), cut-cell solve, surfaces (phi, psi), markers, remesh". Asha's rule of 2026-10-31: work
expected to take more than about two hours on the laptop goes to the Mac Studio. Task 11 estimated these two runs at
two to six hours together on the laptop (section 1), so they come here.

## 1. What this is and why it matters

The repository models the solid AA7075 fragments of a re-entering satellite as they heat, melt and spray droplets. Its
second surface model, the **level set** architecture, carries the body's outer surface as the zero of a level set on an
adaptive box mesh (phi, the gas surface; beta, the base of the surface layer; psi, where the metal turns rigid), with
**skins** (columns of thin cells holding the molten layer) standing on **markers** that sit on phi about 2 mm apart.
One macro step of 0.5 s couples the skins to a cut-cell heat solve of the volume beneath (design
`docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`, section 20).

Stage 3 of the level set was built on the branch `levelset-stage3`. Its Task 11 measures what a macro step costs on a
sphere the size of the thesis's smaller fragment, the 50 mm sphere, at a melting state, with the 20 micrometre
development skin cells (`dev`) and the 10 micrometre production cells (`production`). On the laptop a 0.5 s melting step
of a 10 mm cube took about 20 minutes with about 190 markers; the 50 mm sphere carries 2 404 markers on 1 075 620 cells, so
the laptop could not run it within its two-hour budget. What these runs measure (the split of each step's time between the gas
loads, the skins' passes, the volume's heat solve, the surfaces, the markers, the re-partition, the remesh, the base
moves and the mode changes, and the peak memory) decides where stage 4 must make the model faster before the flights.

| Label | Skin cells | Steps | Writes |
|---|---|---|---|
| `dev` | 20 micrometres | two of 0.5 s from 880 K | `Reentry Model/output/levelset/stage3_cost_sphere_dev.json` |
| `production` | 10 micrometres | two of 0.5 s from 880 K | `Reentry Model/output/levelset/stage3_cost_sphere_production.json` |

**What the laptop measured** (Task 11's probe, once, on mains power, while the laptop was shared with other sessions'
test runs, so the times are high by an unknown factor of up to about two): building the `dev` body took 1 197 s (20
minutes) and a peak of 2.4 GB, and gave 2 404 markers and 1 075 620 cells (118 767 active nodes); one 0.0625 s step from
880 K (one part, before any metal reached the rigid temperature) took 197 s, of which the remesh took about 82 s, psi
24 s, phi 14 s, the runoff 29 s, the cut-cell solves 14 s, the markers 11 s and the shadow test 6 s. A melting 0.5 s step
has eight such parts, each with its base move on the 1.08 million cells and the runoff of liquid, and one step end;
scaled from the q6p box's melting step (126 markers, 391 238 cells, 162 s) by markers and cells, it should take about 15
to 30 minutes for `dev`, and up to about an hour if most skins turn top-zone as the 10 mm cube's did (Task 10: 3
minutes a step before, 20 after). `production` cells double the skins' cells. So expect, per run, the build (10 to 20
minutes) and two steps of 15 to 60 minutes, and peak memory of 2.5 to 4 GB.

## 2. Rules

- Do not edit anything under `Reentry Model/reentry_model/`, `analysis/`, `tests/` or `docs/`, except creating the
  results directory of Step 7. If you think code is wrong, report it.
- Do not change any setting and do not re-run a measurement with different settings. A run that stops with an error is
  a result to report (its JSON records the error and the steps before it), not something to fix.
- Do not change system settings (sleep, updates, energy). Run nothing else heavy on the machine while these run.
- Git: never use bare `git stash`. Push only the results branch named in Step 7, never `main` or `levelset-stage3`. Do
  not merge anything.
- Quote every path: the repository's path on Asha's machines contains spaces.
- When in doubt, stop and ask Asha. Report in Asha's style: full sentences, every number with what was measured, its
  units and why it matters, terms spelled out, words rather than symbols in prose (Asha's `~/.claude/CLAUDE.md`).

## 3. Get the code

**Before this runbook can be run, the branch `levelset-stage3` must be on `origin`.** As Task 11 left it, the branch
was committed on the laptop and not pushed (the stage's rule: push after the final review). Asha (or the session that
reviews Task 11) pushes it; until then, stop at this step.

- [ ] **Step 1.** In the Mac Studio's clone (`git@github.com:ashajain652/particle_formation.git`):

```bash
cd "<the clone's root>"
git status --short          # must show no modified tracked files; if it does, stop and ask Asha
git fetch origin
git switch levelset-stage3 || git switch -c levelset-stage3 --track origin/levelset-stage3
git pull --ff-only
git merge-base --is-ancestor 4156d4e HEAD && echo "code is at or after Task 11's runoff fix"
git log --oneline -1        # note this commit for the report
cd "Reentry Model"
ls analysis/levelset/stage3_cases.py tests/levelset/test_levelset_stage3_cases.py   # both must exist
```

The `merge-base` line must print its message (`4156d4e` is the runoff fix Task 11 committed just before the script and
this runbook; the `ls` line checks that those are there too). Every command below runs from `Reentry Model/`.

## 4. The Python environment

- [ ] **Step 2.**

```bash
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python      # adjust if the Mac Studio's path differs
export PY
"$PY" -c "import sys, numpy, scipy, skfem, pyamg; print(sys.version.split()[0], numpy.__version__, scipy.__version__, skfem.__version__, pyamg.__version__)"
```

On the laptop this prints `3.12.14 2.5.3 1.18.1 12.0.2 5.3.0` (Python, numpy, scipy, scikit-fem, pyamg). A different
numpy, scipy or pyamg version is worth noting in the report (the timings depend on them).

## 5. Check the machine first (minutes)

- [ ] **Step 3: the level set's tests and the cold plate.**

```bash
PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests/levelset -m "not reference" -q -p no:cacheprovider -n 6
PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests/levelset/test_levelset_stage3_cases.py -m reference -q -p no:cacheprovider \
    -k "cold_plate and dev"
```

Expected: `181 passed, 6 deselected`, then `1 passed` for the cold plate (about a minute; it writes
`output/levelset/stage3_plate_cold_dev.json`). Any failure: stop and report it with its output.

## 6. The two runs

- [ ] **Step 4: launch both, side by side, each on one thread.**

```bash
mkdir -p output/levelset
export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nohup "$PY" analysis/levelset/stage3_cases.py cost --what sphere --preset dev > output/levelset/cost_sphere_dev.log 2>&1 &
nohup "$PY" analysis/levelset/stage3_cases.py cost --what sphere --preset production > output/levelset/cost_sphere_production.log 2>&1 &
```

Each prints a line when its body is built (`sphere dev built in ... s: {...}`), one line per macro step
(`sphere dev t=0.50s wall ...`) and finally `wrote .../stage3_cost_sphere_<preset>.json`.

- [ ] **Step 5: watch.** Check every hour or two:

```bash
tail -n 3 output/levelset/cost_sphere_dev.log output/levelset/cost_sphere_production.log
ps aux | grep "stage3_cases.py cost" | grep -v grep
```

Expected: the build within about 20 minutes, each step within about an hour, peak memory about 2.5 to 4 GB each. A line
containing `FAILED` means the step raised: the run then writes its JSON with the error and stops; that is a result.

- [ ] **Step 6: when both have written their JSON**, check them:

```bash
"$PY" - <<'EOF'
import json
for p in ("dev", "production"):
    d = json.load(open("output/levelset/stage3_cost_sphere_{}.json".format(p)))
    print(p, "built", d["built"]["build_s"], "s,", d["built"]["markers"], "markers,", d["built"]["cells"], "cells; failure:", d["failure"])
    for r in d["steps"]:
        print("  t", r["t"], "wall", r["wall_s"], "s, peak", r["peak_rss_gb"], "GB, parts", r["parts"],
              "books", r["mass_residual"], r["enthalpy_residual"])
        print("   split:", {k: v for k, v in r["cost_split_s"].items() if not isinstance(v, dict)})
EOF
```

## 7. Bring the results back

- [ ] **Step 7.** Push the two JSON files and logs on a results branch (they are small):

```bash
git switch -c levelset-stage3-cost-results
mkdir -p analysis/levelset/macstudio
cp output/levelset/stage3_cost_sphere_dev.json output/levelset/stage3_cost_sphere_production.json \
   output/levelset/cost_sphere_dev.log output/levelset/cost_sphere_production.log analysis/levelset/macstudio/
git add analysis/levelset/macstudio
git commit -m "Level-set stage 3 Task 11: the 50 mm sphere's cost per macro step at dev and production cells, measured on the Mac Studio."
git push -u origin levelset-stage3-cost-results
```

On the laptop the results are taken in by copying the two JSON files into `Reentry Model/output/levelset/` and running
`"$PY" analysis/levelset/stage3_cases.py record`, which gathers every `stage3_*.json` there into the committed record
`analysis/levelset/stage3_recorded.json`.

## 8. Report

Report to Asha: the commit you ran, the environment's versions, the tests' results, each run's build time, markers and
cells, each step's wall time and its split (which part of the step dominates), the peak memory, the books (mass and
enthalpy residuals, which should be below 1e-12 and 1e-8), and any failure with its message. Say plainly what was
measured once and on which machine.
