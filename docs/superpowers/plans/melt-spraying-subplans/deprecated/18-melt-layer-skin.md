# Sub-plan 18 — the melt-layer skin: stage 0 (the grid gate) and stage 1 (skins with melting) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read `00-shared-context.md` first (required reading for every Step 3 sub-plan), then the spec.

**Goal:** Measure whether Step 3's converged droplet population depends on the surface mesh (stage 0), then give every
melting surface patch a one-dimensional skin that resolves the melt layer and feeds the film continuously, coupled to
the 3D conduction through an implicit, energy-exact interface (stage 1).

**Architecture:** A new module `reentry_model/skin.py` holds every skin as arrays (cells of metal mass and enthalpy along
the patch normal) and solves them together on short sub-steps. `MeltingBody` creates skins where melting begins, runs a
trial pass to get each skin's base-flux law, solves the 3D step with that law as an affine boundary condition (both
thermal backends gain it), runs the real pass, books the interface mismatch, and hands skins on at element deaths.
`--surface-model elements` (the default in stage 1) keeps every earlier run bit for bit.

**Tech Stack:** Python 3.12 (`drama_env`), NumPy, SciPy, scikit-fem backend, FEniCSx 0.11 backend (`fenicsx_env`), pytest.

**Spec:** `docs/superpowers/specs/2026-10-07-melt-layer-skin-design.md` (approved and revised 2026-10-07). This
sub-plan covers its stages 0 and 1 (§12.1). Stages 2 (the film and spray on the sub-steps) and 3 (the slurry cascade)
get sub-plan 19, written after stage 1's measurements fix the presets and the interface tolerance this plan measures —
the spec gates each stage on the previous one's measurement, and stage 0 may change whether stages 1–3 proceed at all.

**How this plan's code was checked (2026-10-07).** Every code block and diff block of Tasks 2–8 was first applied to a
scratch copy of `prototype/work-2026-10-07-scheil-default/code/` and tested there (unit tier 252 to 286 passed, the same
five known failures), then rebased onto `prototype/work-2026-10-07-dt-continuum/code/` (facts 97–100), the latest
committed state, where only Task 8's CLI diff needed merging by hand. On that base: the skin module's 17 tests, the 2 new
scikit-fem backend tests, the 2 new FEniCSx tests, the 8 body tests, the 7 CLI tests, and the whole unit tier, which
goes from 258 to 292 passed (17 skin, 2 thermal, 8 body, 7 CLI), with 1 skipped and the same five known failures, and the FEniCSx tier from 12 to 14 passed. On the earlier base the default element model was also checked at flight level: the 100 mm melting flight to 65 s, run from the original and the changed copy with the same seed, agrees bit for bit in all 90 shared history columns over 131 rows, all 58,972 droplet-source rows and every result field but the run time. The same check on the new base (to 55 s, across the step switch) was still running when this plan was committed; Task 8 Step 4 repeats it. The blocks below are generated from the tested files, not retyped, and replaying them in
order on a fresh copy of the dt-continuum copy reproduces the tested copy file for file. Measurements made while checking are quoted
where they set a tolerance or explain a choice; they were each made once, on the cases named.

## Global Constraints

- Work in a throwaway copy of the prototype, never in `prototype/proto3/`, the package `reentry_model/` at the repo
  root, the specs or the monolithic plan. The prototype is git-ignored: there are **no commits during the build**.
- Interpreters: `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`;
  `FX=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/python`, run with
  `FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang`.
- Unit tier, run in the copy: `"$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py`.
  Record its counts before any change and after every task. The only acceptable failures are the five known ones that
  need Task 11's melting SESAM references (the `--dt-continuum` amendment's baseline: 258 passed, 1 skipped, those five;
  FEniCSx tier 12 passed).
- `--surface-model elements` (the `MeltSettings` and CLI default) must reproduce the copy's earlier flights **bit for
  bit** in every history column they share, every source-table row and every result field but the run time. New
  history columns and result fields are added; none is changed.
- A run with skins is **bit-identical to the element model until the first skin is created**.
- Bit-identity is always under the same seed and call sequence: pyamg's AMG set-up draws on numpy's global generator, so
  two solvers stepped from identical input differ by about 2e-11 K unless each step is seeded alike (measured while
  planning). The model's `--seed` (fact 62) does this for a run; tests seed before each compared step.
- Mass books exact (relative 1e-12) and the coupled energy balance (`energy_balance_residual`) below 1e-8 in every test
  that melts; flights below 1e-8.
- Melting default material: `AA7075_scheil` (fact 88). Seed: default 12345 (fact 62). Default macro step after the
  continuum switch: 0.0125 s (`--dt-continuum`, facts 97–100), so a default 100 mm flight to 120 s takes about 2.3 h;
  every skin run states its post-switch step explicitly.
- Exit codes: 0 ok, 1 model failure, 2 bad arguments or a missing optional library. Every new CLI value is validated;
  malformed values exit 2.
- Run names encode the configuration: a skin run appends `_surface-skin` and its non-default settings' suffix (Task 8).
- The CLI default stays `--surface-model elements` in stage 1; it becomes `skin` only when stage 2 passes the spec's
  acceptance (§12.3), so no default run changes before the skin is verified.
- Long runs under `caffeinate -i`, on mains power, the power source logged at start and end (`pmset -g batt`).
- Results are recorded only in `docs/superpowers/plans/melt-spraying-subplans/`: facts in the shared context, tested
  code as diffs in the affected sub-plans (Task 9). Nothing outside that folder is edited.

## Review Focus

1. A skin whose element cannot supply its target thickness (an element near `PHI_DEATH`): the skin is made thin, never
   negative, counted, and the element dies at `PHI_DEATH` if it reaches it. Pinned by
   `test_a_skin_its_element_cannot_supply_is_made_thin_and_counted` (Task 7).
2. Liquid on a patch with no skin (it ran there during the previous step): it rides the 3D nodes until the next step
   creates the patch's skin, and its energy moves onto the skin's top exactly. Pinned by
   `test_liquid_on_a_cold_patch_moves_onto_its_new_skin_with_its_energy` (Task 7).
3. A macro step whose interface mismatch exceeds the tolerance (heating that changes fast, or a skin feeding — the base
   law is good to about 1 % per kelvin while a skin feeds, measured while planning): the step repeats once and the books
   stay exact. Pinned by `test_a_large_interface_mismatch_repeats_the_step` (Task 7).
4. A death whose element exposes no face (orphan), or skins left without mass: the hand-over falls back to the nearest
   patch and no empty skin survives. Pinned by `test_hand_over_falls_back_and_drops_empty_skins` (Task 2) and
   `test_deaths_hand_skins_to_the_faces_they_expose` (Task 7).
5. A whole skin turning liquid in one sub-step with nothing left to draw: no NaN, no negative mass, the liquid goes to
   the patch's liquid, and the next sub-step still solves. Pinned by `test_an_exhausted_skin_empties_without_nan`
   (Task 4).

## Where the work is done, and how it reaches the sub-plans

The build happens in `prototype/work-<date>-skin/code/`, a copy of the latest tested copy (Task 1). Paths in the tasks
are relative to that `code/` directory. Tasks 6–8 give their changes as unified diff blocks (`a/` and `b/` prefixes):
apply each with `patch -p1` from `code/`, or by hand where a later amendment has moved the context lines (the changes
are additive). When stage 1's tests pass and its measurements are made (Task 9), the tested changes are written into
the affected sub-plans as dated amendment sections with diff blocks against the copy as it stood before the change, and
the measurements into the shared context as new facts — the repository's rule for every Step 3 change.

## File structure

| File (in the copy) | Responsibility | Task |
|---|---|---|
| `reentry_model/skin.py` (new) | `SkinSettings`, `SkinBook`, `SkinField`: storage, rezoning, the top node, depths, hand-over, the implicit sub-step, feed and draw, the trial and real passes | 2–5 |
| `tests/test_reentry_model_skin.py` (new) | the skin alone, against exact solutions and its own books | 2–5 |
| `reentry_model/thermal/__init__.py` | `StepResult.Q_interface`; the protocol's `interface`, `state`, `set_state` | 6 |
| `reentry_model/thermal/skfem_backend.py`, `fenicsx_backend.py` | the per-facet affine interface flux, no radiation under skins, the state round trip | 6 |
| `tests/test_reentry_model_thermal.py`, `tests/test_reentry_model_fenicsx.py` | backend tests for the interface; the two-backend skin test | 6, 7 |
| `reentry_model/body.py` | `MeltSettings.surface_model`/`skin`; skins in `MeltingBody`; the skin history columns | 7, 8 |
| `tests/test_reentry_model_melting_skin.py` (new) | the melting body with skins | 7 |
| `reentry_model/coupled.py`, `reentry_model/cli.py` | columns, frame fields, results, flags, run names | 8 |
| `tests/test_reentry_model_cli.py` | flags, names, a short skin run | 8 |

---

### Task 0: Stage 0 — the grid gate on the current model

**Files:**
- Read: the 2.0 mm, 0.0125 s runs of facts 91–95 (`prototype/work-2026-10-07-scheil-default/runs/`, both seeds) and
  the switched-step comparison script beside them (`.../harness/compare.py`)
- Create: `prototype/work-<date>-gridgate/` (a frozen measurement copy and its runs)

**Interfaces:**
- Consumes: the model as committed (facts 88–100). Its default `--dt-continuum` of 0.0125 s reproduces the fact-91
  harness's 0.0125 s run bit for bit (fact 99), so the existing 2.0 mm runs are the comparison and no harness is needed.
- Produces: one fact in the shared context — whether the converged droplet population moves with the surface mesh.

- [ ] **Step 1: Freeze a measurement copy**

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
G="$(pwd)/work-$(date +%F)-gridgate"; mkdir -p "$G"
rsync -a --exclude __pycache__ --exclude reentry_model_output work-2026-10-07-dt-continuum/code/ "$G/meas/"
cp work-2026-10-07-scheil-default/harness/compare.py "$G/"
```

Confirm the copy is the committed state before using it: its unit tier gives the baseline counts above.

- [ ] **Step 2: Launch the two 1.4 mm runs (seeds 12345 and 1), side by side**

The 100 mm flight of facts 91–95 with every default (Scheil, the rule on, `--dt-continuum` 0.0125 s) and
`--h-surface 1.4`:

```bash
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
for SEED in 12345 1; do
  OUT="$G/runs/hs14_seed$SEED"; mkdir -p "$OUT"
  ( cd "$G/meas" && echo "power at start: $(pmset -g batt | head -1)" > "$OUT.log" && \
    caffeinate -i "$PY" -m reentry_model run --diameter 100 --altitude 77.500133 --velocity 7.5 \
      --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics --h-surface 1.4 \
      --seed $SEED --t-max 120 --outdir "$OUT" >> "$OUT.log" 2>&1; echo "exit $?" >> "$OUT.log"; \
    echo "power at end: $(pmset -g batt | head -1)" >> "$OUT.log" ) &
done
```

Expected: both logs end `exit 0`, and each run's JSON records the switch at about 49.5 s. The 2.0 mm default run took
2.3 h (fact 99); the 1.4 mm mesh has roughly 2.9 times the surface elements, so expect something like 7 to 10 h each
with both running at once (an estimate, unmeasured).

- [ ] **Step 3: Compare with the 2.0 mm series**

Run `compare.py` on the two 1.4 mm runs and on the 2.0 mm 0.0125 s runs (both seeds) of facts 91–95, over 49.5–120 s.
Tabulate, for each mesh: droplet count, the thick branch's share of the sprayed mass, median radius by number and by
mass, sprayed mass, mass at 120 s, the runoff measures and the run time. Read every change against the scatter between
the two seeds **on the same mesh**, and against fact 93's scatter at 2.0 mm (0.74 % in count, 0.31 points in thick
share, 0.16 % and 0.07 % in the medians, 0.14 % in sprayed mass, 0.53 % in mass at 120 s).

- [ ] **Step 4: Record the verdict and stop for Asha**

Append a section to `00-shared-context.md`, "Amendment of <date> — the melt-layer skin, stage 0: the grid gate (fact
<n>)", opening the way the 2026-09-27 section does: "Design: `docs/superpowers/specs/2026-10-07-melt-layer-skin-design.md`,
approved 2026-10-07. Plan: sub-plan 18 (stages 0 and 1); stages 2 and 3 follow in sub-plan 19." Then one fact (the next
free number) stating which quantities moved beyond the scatter, by how much, and the verdict: **mesh-dependent** (the
skin is needed for correctness) or **mesh-independent** (the skin is a cost measure). Another agent may be appending to
the shared context; re-read the file immediately before the edit and append only. Report to Asha in her reporting style
and wait for her decision before Task 9's measurements. Tasks 1–8 build code only and may proceed while the runs go.

---

### Task 1: The working copy for stage 1

**Files:**
- Create: `prototype/work-<date>-skin/code/` (copy), `prototype/work-<date>-skin/unit_before.log`

- [ ] **Step 1: Copy the latest tested copy and verify it**

Use the tested copy of the latest committed Step 3 amendment, `prototype/work-2026-10-07-dt-continuum/code/` (facts
97–100); this plan's diff blocks are made against it. If a later amendment has been committed since, use its copy and
apply the blocks by hand where the context moved. Verify it first: re-apply every committed amendment's
diff blocks in date order to a fresh copy of `prototype/proto3/` and compare byte for byte (`diff -rq`, excluding
`__pycache__` and `reentry_model_output`). The amendment authors' block-extraction scripts live beside their copies;
reuse one, extended to find every section with a given date, as the 2026-10-05 (runoff flux) amendment did.

```bash
cd "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype"
W="$(pwd)/work-$(date +%F)-skin"; mkdir -p "$W"
rsync -a --exclude __pycache__ --exclude reentry_model_output work-2026-10-07-dt-continuum/code/ "$W/code/"
```

- [ ] **Step 2: Record the unit tier and the FEniCSx tier before any change**

```bash
cd "$W/code" && "$PY" -m pytest -m "not drama and not reference" -q --ignore=tests/test_integration_drama.py --ignore=tests/test_sphere_reentry.py --ignore=tests/test_sphere_sweep.py | tail -3 | tee ../unit_before.log
FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q | tail -2 | tee -a ../unit_before.log
```

Expected: 258 passed, 1 skipped, 2 failed, 3 errors — the five failures are the known ones that need Task 11's melting SESAM references; FEniCSx tier 12 passed.

---

### Task 2: `skin.py` — settings, storage, rezoning, the top node, depths and hand-over

**Files:**
- Create: `reentry_model/skin.py`
- Test: `tests/test_reentry_model_skin.py`

**Interfaces:**
- Consumes: `material.Material` (`rho`, `enthalpy`, `enthalpy_liquid`, `enthalpy_mixed`, `cp_mixed`, `k`,
  `temperature_from_enthalpy`, `temperature_from_enthalpy_mixed`, `T_feed`, `T_rigid`); `thermal.SIGMA_SB`.
- Produces: `SKIN_PRESETS` (`{"production": 10e-6, "dev": 20e-6}` m), `INTERFACE_TOLERANCE` (1e-3), `NEWTON_TOL`
  (1e-6 K), `NEWTON_MAX` (40), `LINE_SEARCH_MAX` (12), `DZ_FLOOR`; `SkinSettings(thickness=0.4e-3, cell=10e-6,
  substep=0.010)` with `.n_cells`; `SkinBook` (per-skin arrays `E_top, E_rad, E_base, to_film, to_deep, drawn,
  drawn_H, L`); `thomas(sub, diag, sup, rhs)`; `cumulative(m, Q, x)`; `SkinField(material, settings=None)` with `.n`,
  `.nc`, `.face_id`, `.area`, `.m`, `.H`, `.T` and `mass()`, `energy()`, `thickness()`, `state()`, `set_state(state)`,
  `add(face_id, area, T_profile, thickness) -> (mass, H)`, `keep(rows)`, `rows_of(face_ids)` (-1 for none),
  `top_energy(rows, L)`, `equilibrate_top(rows, L, E_top) -> T0`, `book_top(rows, dE, L)`, `refreeze_top(rows, taken,
  L_after)`, `_rezone(rows, start, D, hD, S=None, S_drawn=1.0) -> (m, H, S)`, `_set_cells(rows, m, H, L, E_liq_top)`,
  `remap(rows, L)`, `add_bottom(rows, mass, H, L)`, `depth_above(T_level) -> (depth, through)`,
  `hand_over(src, tgt_faces, weights, face_area) -> (receiving rows, mass handed)`. Tasks 3–5 add `conduct`,
  `feed_and_draw` and `advance` to the same class.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_skin.py` with this content (Tasks 3–5 append to it):

```python
"""reentry_model.skin: the one-dimensional skins of the 2026-10-07 design (Step 3 sub-plan 18) -- storage and rezoning,
the top node, depths, hand-over, the implicit sub-step against exact solutions, feed and draw, and the macro step's
trial and real passes."""
import math

import numpy as np
import pytest
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from scipy.special import erf, erfc

from reentry_model import material, skin


def constant_material(latent=0.0, T_melt=float("inf")):
    """rho 2800, c_p 1000, k 150, no radiation; a +-2 K melting range at T_melt when latent > 0."""
    return material.Material("const", 2800.0, 0.0, np.array([200.0, 3000.0]), np.array([1000.0, 1000.0]),
                             np.array([200.0, 3000.0]), np.array([150.0, 150.0]), latent,
                             T_melt - 2.0 if latent else float("inf"), T_melt + 2.0 if latent else float("inf"),
                             material.LiquidProperties(2400.0, 1.3e-3, 0.86))


def field(mat=None, n=3, nc=8, thickness=0.4e-3, T=None):
    mat = mat or material.Material.from_drama_json("AA7075_scheil")
    f = skin.SkinField(mat, skin.SkinSettings(thickness=thickness, cell=thickness / nc))
    T = np.linspace(900.0, 820.0, nc)[None, :].repeat(n, axis=0) if T is None else T
    f.add(np.arange(10, 10 + n), np.full(n, 2.0e-6), T, np.full(n, thickness))
    return f


def node_total(f, row, L):
    """The skin's metal enthalpy plus the liquid riding its top node [J]."""
    return float(f.H[row].sum() + L * f.material.enthalpy_liquid(f.T[row, 0]))


# -- storage, rezoning, the top node, depths, hand-over (Task 2) ----------------------------------------------------
def test_settings_are_validated():
    assert skin.SkinSettings().n_cells == 40 and skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"]).n_cells == 20
    for bad in (dict(thickness=0.0), dict(cell=-1e-6), dict(substep=0.0), dict(thickness=1e-4, cell=6e-5)):
        with pytest.raises(ValueError):
            skin.SkinSettings(**bad)


def test_added_skins_hold_their_slab_and_its_enthalpy():
    f = field()
    mat = f.material
    assert f.n == 3 and f.m.shape == (3, 8)
    assert f.mass() == pytest.approx(3 * mat.rho * 2.0e-6 * 0.4e-3, rel=1e-14)
    assert f.energy() == pytest.approx(float((f.m * mat.enthalpy(f.T)).sum()), rel=1e-14)
    assert np.allclose(f.thickness(), 0.4e-3, rtol=1e-14)
    assert list(f.rows_of(np.array([11, 99, 10]))) == [1, -1, 0]


def test_remap_conserves_mass_and_enthalpy_and_keeps_the_order():
    f = field()
    f.m[0] *= np.linspace(0.5, 1.5, 8)                  # unequal cells
    f.H[0] = f.m[0] * f.material.enthalpy(f.T[0])
    M = f.m[0].sum()
    L = np.array([3e-9])
    before = node_total(f, 0, L[0])
    f.remap(np.array([0]), L)
    assert f.m[0].sum() == pytest.approx(M, rel=1e-14) and np.allclose(f.m[0], M / 8, rtol=1e-14)
    assert node_total(f, 0, L[0]) == pytest.approx(before, rel=1e-12)
    assert np.all(np.diff(f.T[0, 1:]) <= 1e-9)          # still hottest at the top


def test_the_top_node_keeps_its_energy_when_liquid_is_booked_or_refreezes():
    f = field()
    rows, L = np.array([1]), np.array([5e-9])
    E0 = f.top_energy(rows, L)
    f.book_top(rows, np.array([1e-4]), L)
    assert f.top_energy(rows, L) == pytest.approx(E0 + 1e-4, rel=1e-12)
    before, T1, M = node_total(f, 1, 5e-9), f.T[1, 0], f.m[1].sum()
    f.refreeze_top(rows, np.array([4e-9]), np.array([1e-9]))     # most of the liquid freezes onto cell 0
    assert f.m[1].sum() == pytest.approx(M + 4e-9, rel=1e-12)
    assert node_total(f, 1, 1e-9) == pytest.approx(before, rel=1e-12)
    assert f.T[1, 0] >= T1 - 1e-9                       # freezing releases latent heat


def test_depth_above_interpolates_between_cell_centres():
    mat = constant_material()
    nc, s = 10, 1.0e-3
    T = np.linspace(1000.0, 900.0, nc)[None, :]         # linear: T = 1000 - 100 (z - dz/2)/(s - dz)
    f = field(mat, n=1, nc=nc, thickness=s, T=T)
    dz = s / nc
    depth, through = f.depth_above(950.0)
    assert depth[0] == pytest.approx(dz / 2 + 0.5 * (s - dz), rel=1e-12) and not through[0]
    depth, through = f.depth_above(800.0)
    assert depth[0] == pytest.approx(s) and through[0]
    depth, _ = f.depth_above(1200.0)
    assert depth[0] == 0.0


def test_hand_over_passes_skins_layer_by_layer():
    f = field()
    M, E = f.mass(), f.energy()
    areas = {100: 1.5e-6, 101: 2.5e-6}
    tgt = np.array([[100, 101], [100, -1]])
    w = np.array([[0.4, 0.6], [1.0, 0.0]])
    _, handed = f.hand_over(np.array([0, 1]), tgt, w, lambda faces: np.array([areas[int(x)] for x in faces]))
    assert handed == pytest.approx(float(f.m[[0, 1]].sum()), rel=1e-14)
    f.keep(np.setdiff1d(np.arange(f.n), [0, 1]))
    assert f.mass() == pytest.approx(M, rel=1e-13) and f.energy() == pytest.approx(E, rel=1e-12)
    r = f.rows_of(np.array([100, 101]))
    assert (r >= 0).all() and np.all(np.diff(f.T[r[0]]) <= 1e-9)          # still hottest at the top
    assert np.allclose(f.m[r[0]], f.m[r[0]].mean(), rtol=1e-12)          # equal cells


def test_hand_over_falls_back_and_drops_empty_skins():
    f = field()
    f.m[2] = 0.0
    f.H[2] = 0.0
    _, handed = f.hand_over(np.array([2]), np.array([[105, -1]]), np.array([[1.0, 0.0]]), lambda faces: np.full(len(faces), 1e-6))
    f.keep(np.flatnonzero(f.m.sum(axis=1) > 0.0))
    assert f.rows_of(np.array([105]))[0] == -1 and handed == 0.0 and f.n == 2
```

- [ ] **Step 2: Run the tests to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: FAIL with `ImportError: cannot import name 'skin'`.

- [ ] **Step 3: Write the module's first half**

Create `reentry_model/skin.py`:

```python
"""One-dimensional skins under the melting surface patches (Step 3 sub-plan 18; design:
docs/superpowers/specs/2026-10-07-melt-layer-skin-design.md, sections 5-7).

A skin belongs to one surface patch and has its area: a slab of `n_cells` equal cells along the patch's inward normal
that owns the outermost metal of that patch. Its state is each cell's metal mass `m` [kg] and metal enthalpy `H`
[J above material.T_REF], with the temperature `T` that follows from them. The liquid riding on the patch -- the film
and the deep account, which MeltingBody keeps as `m_f` and `m_d` -- shares the top cell's temperature: cell 0 and that
liquid are one node, of energy m_0 h(T_0) + L h_l(T_0). Every operation conserves mass and enthalpy exactly. Every
array is per skin, and all skins are solved together."""
import math
from dataclasses import dataclass

import numpy as np

from .thermal import SIGMA_SB

SKIN_PRESETS = {"production": 10.0e-6, "dev": 20.0e-6}  # m, the presets' cell sizes (spec section 11.1)
INTERFACE_TOLERANCE = 1.0e-3     # |mismatch| / (|E_top| + |E_base|) above which a macro step repeats (spec 6.2)
NEWTON_TOL = 1.0e-6              # K: a sub-step has converged when no cell moves by more than this
NEWTON_MAX = 40
LINE_SEARCH_MAX = 12           # halvings of a Newton step before it is taken anyway
DZ_FLOOR = 1.0e-12               # m: an emptied cell conducts as if this thin (it holds no heat)


@dataclass
class SkinSettings:
    thickness: float = 0.4e-3              # m, target thickness s
    cell: float = SKIN_PRESETS["production"]
    substep: float = 0.010                 # s, the skins' sub-step

    def __post_init__(self):
        if not (self.thickness > 0.0 and self.cell > 0.0 and self.substep > 0.0):
            raise ValueError("skin thickness, cell and substep must be positive")
        if 2.0 * self.cell > self.thickness * (1.0 + 1e-12):
            raise ValueError("a skin needs at least two cells: cell <= thickness / 2")

    @property
    def n_cells(self):
        return int(round(self.thickness / self.cell))


@dataclass
class SkinBook:
    """What a real pass did, per skin: energies [J] and masses [kg]."""
    E_top: np.ndarray        # convective heat taken in at the top
    E_rad: np.ndarray        # radiated from the top
    E_base: np.ndarray       # passed down into the 3D model through the base
    to_film: np.ndarray      # liquid handed to the film
    to_deep: np.ndarray      # liquid handed to the deep account
    drawn: np.ndarray        # metal drawn from the element below
    drawn_H: np.ndarray      # the enthalpy that metal brought
    L: np.ndarray            # the liquid on each top node after the pass


def thomas(sub, diag, sup, rhs):
    """Row-wise tridiagonal solve: sub and sup (n, k-1), diag and rhs (n, k). Returns x (n, k)."""
    n, k = diag.shape
    c = np.zeros((n, max(k - 1, 1)))
    d = np.zeros((n, k))
    if k > 1:
        c[:, 0] = sup[:, 0] / diag[:, 0]
    d[:, 0] = rhs[:, 0] / diag[:, 0]
    for i in range(1, k):
        denom = diag[:, i] - sub[:, i - 1] * c[:, i - 1]
        if i < k - 1:
            c[:, i] = sup[:, i] / denom
        d[:, i] = (rhs[:, i] - sub[:, i - 1] * d[:, i - 1]) / denom
    x = np.zeros((n, k))
    x[:, -1] = d[:, -1]
    for i in range(k - 2, -1, -1):
        x[:, i] = d[:, i] - c[:, i] * x[:, i + 1]
    return x


def cumulative(m, Q, x):
    """Per row: the integral from the top down to mass coordinate x (n, j) of a quantity Q (n, nc) carried uniformly
    within each cell of mass m (n, nc). Exact for piecewise-uniform cells; one sorted search for all rows."""
    n, nc = m.shape
    cm = np.concatenate([np.zeros((n, 1)), np.cumsum(m, axis=1)], axis=1)
    cQ = np.concatenate([np.zeros((n, 1)), np.cumsum(Q, axis=1)], axis=1)
    M = np.where(cm[:, -1] > 0.0, cm[:, -1], 1.0)[:, None]
    off = 2.0 * np.arange(n)[:, None]
    flat = (cm / M + off).ravel()
    i = np.searchsorted(flat, (x / M + off).ravel(), side="right").reshape(x.shape) - 1
    i = np.clip(i - (nc + 1) * np.arange(n)[:, None], 0, nc - 1)
    r = np.arange(n)[:, None]
    mi = m[r, i]
    frac = np.clip(np.divide(x - cm[r, i], mi, out=np.zeros_like(x, dtype=float), where=mi > 0.0), 0.0, 1.0)
    return cQ[r, i] + frac * Q[r, i]


class SkinField:
    """Every skin of the body, as arrays (module docstring)."""

    def __init__(self, material, settings=None):
        self.material = material
        self.settings = settings or SkinSettings()
        self.nc = self.settings.n_cells
        self.face_id = np.zeros(0, dtype=np.int64)
        self.area = np.zeros(0)
        self.m = np.zeros((0, self.nc))
        self.H = np.zeros((0, self.nc))
        self.T = np.zeros((0, self.nc))

    # -- books ------------------------------------------------------------------------------------------------------
    @property
    def n(self):
        return int(self.face_id.size)

    def mass(self):
        return float(self.m.sum())

    def energy(self):
        """Metal enthalpy held by the skins [J]. The liquid on their top nodes is MeltingBody.film_energy's."""
        return float(self.H.sum())

    def thickness(self):
        return self.m.sum(axis=1) / (self.material.rho * self.area) if self.n else np.zeros(0)

    def state(self):
        return tuple(a.copy() for a in (self.face_id, self.area, self.m, self.H, self.T))

    def set_state(self, state):
        self.face_id, self.area, self.m, self.H, self.T = (a.copy() for a in state)

    # -- membership -------------------------------------------------------------------------------------------------
    def add(self, face_id, area, T_profile, thickness):
        """Append skins with cell temperatures T_profile (k, nc) and thickness [m] (k,). Returns (mass, H) per new skin."""
        area = np.asarray(area, dtype=float)
        dz = np.asarray(thickness, dtype=float) / self.nc
        m = self.material.rho * (area * dz)[:, None] * np.ones((1, self.nc))
        T = np.asarray(T_profile, dtype=float).reshape(len(area), self.nc)
        H = m * self.material.enthalpy(T)
        self.face_id = np.concatenate([self.face_id, np.asarray(face_id, dtype=np.int64)])
        self.area = np.concatenate([self.area, area])
        self.m, self.H, self.T = np.vstack([self.m, m]), np.vstack([self.H, H]), np.vstack([self.T, T])
        return m.sum(axis=1), H.sum(axis=1)

    def keep(self, rows):
        rows = np.asarray(rows, dtype=np.int64)
        self.face_id, self.area = self.face_id[rows], self.area[rows]
        self.m, self.H, self.T = self.m[rows], self.H[rows], self.T[rows]

    def rows_of(self, face_ids):
        """The row of each face id's skin, -1 where it has none."""
        face_ids = np.asarray(face_ids, dtype=np.int64)
        if not self.n:
            return np.full(face_ids.shape, -1, dtype=np.int64)
        order = np.argsort(self.face_id)
        s = self.face_id[order]
        i = np.clip(np.searchsorted(s, face_ids), 0, len(s) - 1)
        return np.where(s[i] == face_ids, order[i], -1)

    # -- the top node -----------------------------------------------------------------------------------------------
    def top_energy(self, rows, L):
        """Energy of the top node of `rows` with liquid L [kg] on it [J]."""
        return self.H[rows, 0] + L * self.material.enthalpy_liquid(self.T[rows, 0])

    def equilibrate_top(self, rows, L, E_top):
        """Put cell 0 of `rows` and the liquid L [kg] on it at one temperature holding E_top [J] together. Returns T_0."""
        mat = self.material
        m0 = self.m[rows, 0]
        tot = m0 + L
        w = np.divide(L, tot, out=np.zeros_like(tot), where=tot > 0.0)
        h = np.divide(E_top, tot, out=np.zeros_like(tot), where=tot > 0.0)
        T0 = np.where(tot > 0.0, mat.temperature_from_enthalpy_mixed(h, w), self.T[rows, 0])
        self.T[rows, 0] = T0
        self.H[rows, 0] = m0 * mat.enthalpy(T0)
        return T0

    def book_top(self, rows, dE, L):
        """Add dE [J] to the top node of `rows` carrying liquid L [kg] and re-solve its temperature."""
        self.equilibrate_top(rows, L, self.top_energy(rows, L) + dE)

    def refreeze_top(self, rows, taken, L_after):
        """Liquid `taken` [kg] freezes onto cell 0 of `rows`; the node's energy is unchanged (spec 7.1 (6)). L_after is
        the liquid left on the node; the skin is remapped onto equal cells afterwards."""
        E = self.top_energy(rows, L_after + taken)
        self.m[rows, 0] += taken
        self.equilibrate_top(rows, L_after, E)
        self.remap(rows, L_after)

    # -- rezoning ---------------------------------------------------------------------------------------------------
    def _rezone(self, rows, start, D, hD, S=None, S_drawn=1.0):
        """Rows' columns below mass coordinate `start` [kg from the top], with D [kg] of metal added at the bottom at
        specific enthalpy hD [J/kg], redistributed onto nc equal cells: within an old cell the specific enthalpy is
        uniform, so the cumulative enthalpy is linear in the cumulative mass and the new cells' enthalpies are exact
        differences of it. S (an intensive per-cell field, optional) is remapped mass-weighted, the drawn metal
        carrying S_drawn. Returns (m, H, S) of the new cells."""
        m, H = self.m[rows], self.H[rows]
        nc = self.nc
        start, D, hD = (np.broadcast_to(np.asarray(v, dtype=float), (len(rows),)) for v in (start, D, hD))
        rest = m.sum(axis=1) - start
        M_new = rest + D
        y = M_new[:, None] * np.linspace(0.0, 1.0, nc + 1)[None, :]
        x = start[:, None] + np.minimum(y, rest[:, None])
        below = np.maximum(y - rest[:, None], 0.0)
        E_top = cumulative(m, H, start[:, None])
        cH = cumulative(m, H, x) - E_top + below * hD[:, None]
        cH[:, 0] = 0.0
        cH[:, -1] = H.sum(axis=1) - E_top[:, 0] + D * hD                      # exact totals
        m_new = np.repeat((M_new / nc)[:, None], nc, axis=1)
        S_new = None
        if S is not None:
            mS = m * S
            cS = cumulative(m, mS, x) - cumulative(m, mS, start[:, None]) + below * S_drawn
            S_new = np.divide(np.diff(cS, axis=1), m_new, out=np.ones_like(m_new) * S_drawn, where=m_new > 0.0)
        return m_new, np.diff(cH, axis=1), S_new

    def _set_cells(self, rows, m, H, L, E_liq_top):
        """Install cell masses and enthalpies on `rows`, the temperatures from them, the top node with its liquid."""
        mat = self.material
        self.m[rows], self.H[rows] = m, H
        h = np.divide(H, m, out=mat.enthalpy(self.T[rows]), where=m > 0.0)
        self.T[rows] = mat.temperature_from_enthalpy(h)
        self.equilibrate_top(rows, L, H[:, 0] + E_liq_top)

    def remap(self, rows, L):
        """Conservative remap of `rows` (any cell masses) onto nc equal cells; L [kg] is the liquid on their top nodes,
        whose energy stays with the top node."""
        rows = np.asarray(rows, dtype=np.int64)
        if not rows.size:
            return
        E_liq = L * self.material.enthalpy_liquid(self.T[rows, 0])
        m, H, _ = self._rezone(rows, 0.0, 0.0, 0.0)
        self._set_cells(rows, m, H, L, E_liq)

    def add_bottom(self, rows, mass, H, L):
        """Metal `mass` [kg] carrying `H` [J] joins the bottom of `rows` (distinct rows; an element's remainder at its
        death, spec section 9); L [kg] is the liquid on their top nodes."""
        rows = np.asarray(rows, dtype=np.int64)
        self.m[rows, -1] += mass
        self.H[rows, -1] += H
        self.remap(rows, L)

    # -- depths -----------------------------------------------------------------------------------------------------
    def depth_above(self, T_level):
        """Per skin: the depth [m] from the top to where the contiguous cells above T_level end, interpolated between
        cell centres, and whether every cell is above it (the layer then continues into the 3D model)."""
        if not self.n:
            return np.zeros(0), np.zeros(0, dtype=bool)
        k = np.cumprod(self.T > T_level, axis=1).sum(axis=1)
        dz = self.thickness() / self.nc
        r = np.arange(self.n)
        Ta = self.T[r, np.clip(k - 1, 0, self.nc - 1)]
        Tb = self.T[r, np.clip(k, 0, self.nc - 1)]
        frac = np.clip(np.divide(Ta - T_level, Ta - Tb, out=np.zeros_like(Ta), where=Ta > Tb), 0.0, 1.0)
        depth = np.where(k == 0, 0.0, np.where(k == self.nc, self.nc * dz, (k - 0.5 + frac) * dz))
        return depth, k == self.nc

    # -- hand-over --------------------------------------------------------------------------------------------------
    def hand_over(self, src, tgt_faces, weights, face_area):
        """Skins of vanished patches pass to the patches their deaths exposed, layer by layer (spec section 9).
        src (k,) rows; tgt_faces (k, j) face ids of the receiving patches (-1: none); weights (k, j), summing to 1 over
        the valid targets; face_area(face_ids) -> areas. Every source has equal cells, so cell-wise weighted sums do too.
        Returns (receiving rows, mass handed over). The caller removes `src` and any skin left without mass; the
        receiving skins' top nodes are not re-equilibrated with the liquid on them (the caller books that liquid's
        change of enthalpy, as the liquid's own hand-over does)."""
        src = np.asarray(src, dtype=np.int64)
        ok = (tgt_faces >= 0) & (weights > 0.0)
        if not src.size or not ok.any():
            return np.zeros(0, dtype=np.int64), 0.0
        faces = np.unique(tgt_faces[ok])
        new = faces[self.rows_of(faces) < 0]
        if new.size:
            k = new.size
            self.face_id = np.concatenate([self.face_id, new])
            self.area = np.concatenate([self.area, np.asarray(face_area(new), dtype=float)])
            self.m, self.H = np.vstack([self.m, np.zeros((k, self.nc))]), np.vstack([self.H, np.zeros((k, self.nc))])
            self.T = np.vstack([self.T, np.zeros((k, self.nc))])
        sj, tj = np.nonzero(ok)
        trows = self.rows_of(tgt_faces[sj, tj])
        w = weights[sj, tj][:, None]
        np.add.at(self.m, trows, w * self.m[src[sj]])
        np.add.at(self.H, trows, w * self.H[src[sj]])
        recv = np.unique(trows)
        full = recv[self.m[recv].sum(axis=1) > 0.0]
        if full.size:
            self.T[full] = self.material.temperature_from_enthalpy(self.H[full] / np.where(self.m[full] > 0.0, self.m[full], 1.0))
        return recv, float(self.m[src].sum())
```

- [ ] **Step 4: Run the tests to make sure they pass**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: 7 passed.

---

### Task 3: `skin.py` — the implicit sub-step

**Files:**
- Modify: `reentry_model/skin.py` (append `SkinField.conduct` to the class)
- Test: `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Consumes: Task 2's `SkinField`, `thomas`, the module constants.
- Produces: `SkinField.conduct(dt, q_top, T_amb, emissivity, T_b, L, T_surface=None, sens=None) -> (E_top, E_rad,
  E_base, sens)`, energies in J per skin; `sens = (S, dEb)`, `S` (n, nc) the derivative of the cell temperatures with
  respect to the base temperature and `dEb` (n,) that of the base energy, carried through the sub-steps of a trial pass.

Why the iteration is plain Newton in the temperature with a line search, not the 3D solver's enthalpy-consistent
update: in 10 µm cells conduction outweighs the heat capacity some five hundred times over a 10 ms sub-step (a cell's
conductance to its neighbour 38 W/K against a capacity rate of 0.07 W/K, measured while planning), so the cells must move
together. Inverting each cell's own h(T) moved neighbouring cells apart by a tenth of a 7 K step wherever the Scheil
curve's effective heat capacity varies, and the iteration stalled at a 0.4 K change on 888 of 2000 skins. The line search
(halving the step per skin until the residual's norm falls) does across a sharp latent ramp what that update did: the
two-phase Neumann test below, with a ±2 K ramp, needs it.

- [ ] **Step 1: Append the failing tests**

```python
# -- the implicit sub-step (Task 3) ---------------------------------------------------------------------------------
def test_thomas_matches_a_banded_solver():
    rng = np.random.default_rng(0)
    n, k = 5, 9
    sub, sup = -rng.random((n, k - 1)), -rng.random((n, k - 1))
    diag = 3.0 + rng.random((n, k))
    rhs = rng.random((n, k))
    x = skin.thomas(sub, diag, sup, rhs)
    for i in range(n):
        ab = np.zeros((3, k))
        ab[0, 1:], ab[1], ab[2, :-1] = sup[i], diag[i], sub[i]
        assert np.allclose(x[i], solve_banded((1, 1), ab, rhs[i]), rtol=1e-13, atol=0.0)


def test_constant_flux_on_a_semi_infinite_solid():
    """Carslaw and Jaeger 2.9: T - T_i = (2 q / k) sqrt(alpha t) ierfc(z / (2 sqrt(alpha t))), at the top cell's centre."""
    mat = constant_material()
    s, nc = 4.0e-3, 400
    f = field(mat, n=1, nc=nc, thickness=s, T=np.full((1, nc), 300.0))
    q, t_end, dt = 2.0e6, 0.05, 2.5e-4
    for _ in range(int(round(t_end / dt))):
        f.conduct(dt, np.array([q]), 0.0, np.array([0.0]), np.array([300.0]), np.array([0.0]))
    alpha = 150.0 / (2800.0 * 1000.0)
    x = 0.5 * s / nc / (2.0 * math.sqrt(alpha * t_end))
    ierfc = math.exp(-x * x) / math.sqrt(math.pi) - x * math.erfc(x)
    exact = 300.0 + 2.0 * q / 150.0 * math.sqrt(alpha * t_end) * ierfc
    assert f.T[0, 0] == pytest.approx(exact, abs=0.01 * (exact - 300.0))


def test_neumanns_melting_front():
    """Two-phase Neumann problem, equal properties: the front at 2 lambda sqrt(alpha t), lambda from
    St/(exp(l^2) erf l) - St/(exp(l^2) erfc l) = l sqrt(pi) with St = c (T_s - T_m)/L = c (T_m - T_i)/L; the front is
    where the material is half liquid."""
    mat = constant_material(latent=4.0e5, T_melt=900.0)
    s, nc = 8.0e-3, 400
    f = field(mat, n=1, nc=nc, thickness=s, T=np.full((1, nc), 800.0))
    t_end, dt = 0.05, 2.5e-4
    for _ in range(int(round(t_end / dt))):
        f.conduct(dt, np.array([0.0]), 0.0, np.array([0.0]), np.array([800.0]), np.array([0.0]), T_surface=1000.0)
    St = 1000.0 * 100.0 / 4.0e5
    lam = brentq(lambda l: St / (math.exp(l * l) * erf(l)) - St / (math.exp(l * l) * erfc(l)) - l * math.sqrt(math.pi), 1e-3, 3.0)
    alpha = 150.0 / (2800.0 * 1000.0)
    front, _ = f.depth_above(900.0)
    assert front[0] == pytest.approx(2.0 * lam * math.sqrt(alpha * t_end), rel=0.03)


def test_a_substep_closes_the_skin_books_with_the_conducted_base_flux():
    f = field()
    L = np.array([0.0, 2e-9, 0.0])
    E0 = sum(node_total(f, i, L[i]) for i in range(3))
    E_top, E_rad, E_base, _ = f.conduct(0.01, np.full(3, 1.5e6), 0.0, np.full(3, 0.3), np.full(3, 800.0), L)
    E1 = sum(node_total(f, i, L[i]) for i in range(3))
    assert E1 - E0 == pytest.approx(float((E_top - E_rad - E_base).sum()), rel=1e-12)
    dz = f.thickness() / f.nc
    conducted = f.area * f.material.k(f.T[:, -1]) / (0.5 * dz) * (f.T[:, -1] - 800.0) * 0.01
    assert np.allclose(E_base, conducted, rtol=1e-6)
```

- [ ] **Step 2: Run the tests to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: the Thomas test passes (Task 2 defined `thomas`); the other three FAIL with
`AttributeError: 'SkinField' object has no attribute 'conduct'`.

- [ ] **Step 3: Append `conduct` to `SkinField`**

```python
    # -- the sub-step -----------------------------------------------------------------------------------------------
    def conduct(self, dt, q_top, T_amb, emissivity, T_b, L, T_surface=None, sens=None):
        """One implicit (backward-Euler) sub-step for every skin (spec section 5.4): the patch's convective flux q_top
        [W/m2] in at the top, radiation from the top cell, conduction through the cells (harmonic conductances), and the
        base held at T_b [K]. The liquid L [kg] on each top node shares cell 0's temperature and heat capacity. Newton
        on the enthalpy with the enthalpy-consistent update the 3D solver uses, damped if the change grows. T_surface
        (tests only) holds the top face at that temperature instead of the flux. `sens = (S, dEb)` carries dT/dT_b and
        the derivative of the base energy through the sub-step (the trial pass). Returns (E_top, E_rad, E_base, sens)
        [J per skin]; E_base is what the books leave, so they close exactly whatever the iteration error."""
        mat, n, nc = self.material, self.n, self.nc
        if n == 0:
            z = np.zeros(0)
            return z, z, z, sens
        A = self.area
        dz = np.maximum(self.m / (mat.rho * A[:, None]), DZ_FLOOR)
        Mt = self.m.copy()
        Mt[:, 0] += L
        w = np.zeros_like(Mt)
        w[:, 0] = np.divide(L, Mt[:, 0], out=np.zeros(n), where=Mt[:, 0] > 0.0)
        T_old = self.T.copy()
        E_old = Mt * mat.enthalpy_mixed(T_old, w)
        eps_s = np.asarray(emissivity, dtype=float) * SIGMA_SB
        q = np.asarray(q_top, dtype=float)
        T_b = np.asarray(T_b, dtype=float)

        def system(T):
            """Residual R [W] of every cell at T, the conductances, and the diagonal of the Newton matrix."""
            k = mat.k(T)
            G = A[:, None] / (0.5 * dz[:, :-1] / k[:, :-1] + 0.5 * dz[:, 1:] / k[:, 1:])
            Gb = A * k[:, -1] / (0.5 * dz[:, -1])
            flux = G * (T[:, :-1] - T[:, 1:])
            net = np.zeros((n, nc))
            diag = Mt * mat.cp_mixed(T, w) / dt
            if T_surface is None:
                net[:, 0] += q * A - eps_s * A * (T[:, 0] ** 4 - T_amb ** 4)
                diag[:, 0] += 4.0 * eps_s * A * T[:, 0] ** 3
            else:
                Gs = A * k[:, 0] / (0.5 * dz[:, 0])
                net[:, 0] += Gs * (T_surface - T[:, 0])
                diag[:, 0] += Gs
            net[:, :-1] -= flux
            net[:, 1:] += flux
            net[:, -1] -= Gb * (T[:, -1] - T_b)
            diag[:, :-1] += G
            diag[:, 1:] += G
            diag[:, -1] += Gb
            return (Mt * mat.enthalpy_mixed(T, w) - E_old) / dt - net, G, Gb, diag

        # Newton in the temperature with a backtracking line search per skin on the residual's norm. Not the 3D solver's
        # enthalpy-consistent update: in 10-micrometre cells conduction outweighs the heat capacity some five hundred
        # times over a 10 ms sub-step, so a cell must move with its neighbours, and inverting each cell's own h(T) moved
        # them apart by a tenth of the step wherever c_p,eff varies through the Scheil range -- the iteration stalled at
        # a 0.4 K change (measured while planning). Across a sharp latent ramp the line search does what that update did.
        T = T_old.copy()
        R, G, Gb, diag = system(T)
        for _ in range(NEWTON_MAX):
            dT = thomas(-G, diag, -G, -R)
            r0 = np.sqrt((R * R).sum(axis=1))
            lam = np.ones(n)
            for _ in range(LINE_SEARCH_MAX):
                T_try = T + lam[:, None] * dT
                R_try, G_try, Gb_try, diag_try = system(T_try)
                moved = np.abs(T_try - T).max(axis=1)
                worse = (np.sqrt((R_try * R_try).sum(axis=1)) > (1.0 - 1e-4 * lam) * r0) & (moved > NEWTON_TOL)
                if not worse.any():
                    break
                lam = np.where(worse, 0.5 * lam, lam)
            T, R, G, Gb, diag = T_try, R_try, G_try, Gb_try, diag_try
            if moved.max() <= NEWTON_TOL:
                break
        else:
            raise RuntimeError("a skin sub-step did not converge in {} iterations (last change {:.2e} K)".format(NEWTON_MAX, moved.max()))
        E_new = Mt * mat.enthalpy_mixed(T, w)
        if T_surface is None:
            E_top = q * A * dt
            E_rad = eps_s * A * (T[:, 0] ** 4 - T_amb ** 4) * dt
        else:
            E_top = A * mat.k(T[:, 0]) / (0.5 * dz[:, 0]) * (T_surface - T[:, 0]) * dt
            E_rad = np.zeros(n)
        E_base = E_top - E_rad - (E_new.sum(axis=1) - E_old.sum(axis=1))
        if sens is not None:
            S, dEb = sens
            rhs = Mt * mat.cp_mixed(T_old, w) / dt * S
            rhs[:, -1] += Gb
            S = thomas(-G, diag, -G, rhs)
            sens = (S, dEb + dt * Gb * (S[:, -1] - 1.0))
        self.T = T
        self.H = self.m * mat.enthalpy(T)
        return E_top, E_rad, E_base, sens
```

- [ ] **Step 4: Run the tests to make sure they pass**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: 11 passed.

---

### Task 4: `skin.py` — feed and draw

**Files:**
- Modify: `reentry_model/skin.py` (append `SkinField.feed_and_draw`)
- Test: `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Consumes: Tasks 2–3 (`_rezone`, `_set_cells`).
- Produces: `SkinField.feed_and_draw(L, room, avail, T_b, S=None) -> (to_film, to_deep, drawn, drawn_H, L_new,
  S_new)`, per skin. `room` [kg] is ρ A δ_m of metal where Girin's closure applies and `inf` elsewhere; `avail` [kg]
  is what the element below can still give; `S` is the trial pass's sensitivity, rezoned with the cells (drawn metal,
  held at the base temperature, carries 1).

- [ ] **Step 1: Append the failing tests**

```python
# -- feed and draw (Task 4) -----------------------------------------------------------------------------------------
def test_feed_takes_only_the_contiguous_liquid_at_the_top_and_the_draw_restores_the_thickness():
    mat = material.Material.from_drama_json("AA7075_scheil")
    Tf = mat.T_feed
    T = np.array([[Tf + 5, Tf + 3, Tf + 1, Tf - 30, Tf + 4, 850.0, 840.0, 830.0]])
    f = field(mat, n=1, nc=8, T=T)
    M0 = f.mass()
    mu = M0 / 8
    L = np.array([1e-9])
    before = node_total(f, 0, L[0])
    to_film, to_deep, drawn, drawn_H, L_new, _ = f.feed_and_draw(L, np.array([1.5 * mu]), np.array([1.0]), np.array([820.0]))
    assert to_film[0] + to_deep[0] == pytest.approx(3 * mu, rel=1e-13) and to_film[0] == pytest.approx(1.5 * mu, rel=1e-13)
    assert drawn[0] == pytest.approx(3 * mu, rel=1e-12) and f.mass() == pytest.approx(M0, rel=1e-13)
    assert L_new[0] == pytest.approx(L[0] + 3 * mu, rel=1e-13)
    assert node_total(f, 0, L_new[0]) == pytest.approx(before + drawn_H[0], rel=1e-12)
    assert drawn_H[0] == pytest.approx(drawn[0] * float(mat.enthalpy(820.0)), rel=1e-14)


def test_an_exhausted_skin_empties_without_nan():
    mat = material.Material.from_drama_json("AA7075_scheil")
    f = field(mat, n=1, nc=8, T=np.full((1, 8), mat.T_feed + 10.0))
    M0 = f.mass()
    to_film, to_deep, drawn, _, L_new, _ = f.feed_and_draw(np.array([0.0]), np.array([np.inf]), np.array([0.0]), np.array([900.0]))
    assert to_film[0] == pytest.approx(M0, rel=1e-13) and drawn[0] == 0.0 and f.mass() == 0.0
    assert np.isfinite(f.T).all() and L_new[0] == pytest.approx(M0, rel=1e-13)
    E_top, E_rad, E_base, _ = f.conduct(0.01, np.array([1e6]), 0.0, np.array([0.3]), np.array([900.0]), L_new)
    assert np.isfinite(f.T).all() and np.isfinite(E_base).all()


def test_steady_ablation_speed():
    """Melt removed as it forms under a constant flux: the steady recession speed is q / (rho (c (T_feed - T_i) + L));
    the base is held at the steady profile's value T_i + (T_feed - T_i) exp(-v s / alpha) and draws at it."""
    mat = constant_material(latent=4.0e5, T_melt=900.0)
    s, nc, q, T_i = 4.0e-3, 200, 2.0e7, 700.0
    v = q / (2800.0 * (1000.0 * (mat.T_feed - T_i) + 4.0e5))
    alpha = 150.0 / (2800.0 * 1000.0)
    T_b = T_i + (mat.T_feed - T_i) * math.exp(-v * s / alpha)
    z = (np.arange(nc) + 0.5) * s / nc
    f = field(mat, n=1, nc=nc, thickness=s, T=(T_i + (mat.T_feed - T_i) * np.exp(-v * z / alpha))[None, :])
    fed, dt = 0.0, 1.0e-3
    for step in range(1000):
        f.conduct(dt, np.array([q]), 0.0, np.array([0.0]), np.array([T_b]), np.array([0.0]))
        to_film, to_deep, _, _, _, _ = f.feed_and_draw(np.array([0.0]), np.array([np.inf]), np.array([1.0]), np.array([T_b]))
        if step >= 700:                                  # the liquid leaves at once: nothing rides the top node
            fed += float(to_film[0] + to_deep[0])
    speed = fed / (0.3 * 2800.0 * f.area[0])
    assert speed == pytest.approx(v, rel=0.02)
```

- [ ] **Step 2: Run the tests to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: the three new tests FAIL with `AttributeError: 'SkinField' object has no attribute 'feed_and_draw'`.

- [ ] **Step 3: Append `feed_and_draw` to `SkinField`**

```python
    # -- feed and draw ----------------------------------------------------------------------------------------------
    def feed_and_draw(self, L, room, avail, T_b, S=None):
        """The contiguous fully liquid cells at the top (T >= T_feed) leave the skin for the liquid on its patch -- the
        part within `room` [kg] (rho A delta_m, measured from the top, spec 5.4) to the film, the rest to the deep
        account -- and the skin draws metal from the element below, at most `avail` [kg], at its base temperature T_b,
        back to its target thickness. The remaining column and the drawn metal are rezoned onto equal cells, exactly in
        mass and enthalpy. S (the trial pass's dT/dT_b, optional) is rezoned with them, drawn metal carrying 1.
        Returns (to_film, to_deep, drawn, drawn_H, L_new, S_new) per skin."""
        mat, n = self.material, self.n
        if n == 0:
            z = np.zeros(0)
            return z, z, z, z, z, S
        contig = np.cumprod(self.T >= mat.T_feed, axis=1).astype(bool)
        R = (self.m * contig).sum(axis=1)
        E_rem = (self.H * contig).sum(axis=1)
        to_film = np.minimum(R, room)
        to_deep = R - to_film
        rest = self.m.sum(axis=1) - R
        D = np.minimum(np.clip(mat.rho * self.area * self.settings.thickness - rest, 0.0, None), np.maximum(avail, 0.0))
        hD = mat.enthalpy(np.asarray(T_b, dtype=float) * np.ones(n))
        E_liq = L * mat.enthalpy_liquid(self.T[:, 0])
        rows = np.arange(n)
        m_new, H_new, S_new = self._rezone(rows, R, D, hD, S, 1.0)
        L_new = L + R
        self._set_cells(rows, m_new, H_new, L_new, E_liq + E_rem)
        return to_film, to_deep, D, D * hD, L_new, S_new
```

- [ ] **Step 4: Run the tests to make sure they pass**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: 14 passed (the steady-ablation test takes about 5 s).

---

### Task 5: `skin.py` — the macro step's trial and real passes, and the base law

**Files:**
- Modify: `reentry_model/skin.py` (append `SkinField.advance`)
- Test: `tests/test_reentry_model_skin.py` (append)

**Interfaces:**
- Consumes: Tasks 2–4.
- Produces: `SkinField.advance(dt, q_top, T_amb, emissivity, T_b0, L, room, avail, T_b1=None)`. With `T_b1=None`, the
  trial pass, on a copy, base held at `T_b0`, returning `(a, b)` per skin: the law q_b(T) = a + b T [W/m²] of heat into
  the 3D model through the skin's facet, its slope clipped at zero. Otherwise the real pass, base held at `T_b1` (the 3D
  model's end-of-step facet temperature, spec §6.2 as revised), returning a `SkinBook`.

What the law is good for, measured while planning on Scheil skins: exact at the trial pass's own base temperature
(relative error 1e-16); a kelvin away, good to between 9e-6 and 5e-5 where nothing feeds; while a skin feeds, good to
about 1 % per kelvin (0.93 % at 1 K and 3.4 % at 5 K on one skin; in a 2000-skin stress case up to 1.9 % per kelvin from
a mushy start and 12 % from a start at the liquidus), because the feed comes in whole cells. That error is what the
interface mismatch, its booking and the step's repeat exist for; Task 9 measures how often a flight repeats. A slope
above zero occurred on a skin whose base was held inside the top of the mushy range (+1.65e5 W/(m² K)); it is clipped,
since a law rising with the base temperature would cost the 3D matrix its diagonal dominance.

- [ ] **Step 1: Append the failing tests**

```python
# -- the macro step's passes (Task 5) -------------------------------------------------------------------------------
def heated_field(q=(2.5e7, 2.5e7), T_top=905.0, T_bottom=860.0, T_b0=(855.0, 860.0)):
    """Two Scheil skins under a hard heating that melts and feeds them (by default); the second has little room within
    the conjugate depth, so part of what it feeds goes to the deep account."""
    mat = material.Material.from_drama_json("AA7075_scheil")
    f = field(mat, n=2, T=np.linspace(T_top, T_bottom, 8)[None, :].repeat(2, axis=0))
    return (f, np.array(q), np.array([0.4, 0.4]), np.array([2e-9, 0.0]), np.array([np.inf, 2e-7]), np.array([1.0, 1.0]),
            np.array(T_b0))


def test_the_trial_pass_changes_nothing():
    f, q, eps, L, room, avail, T_b0 = heated_field()
    before = f.state()
    f.advance(0.5, q, 0.0, eps, T_b0, L, room, avail)
    assert all(np.array_equal(x, y) for x, y in zip(before, f.state()))


def test_the_base_law_predicts_the_real_pass():
    """Exact at the trial pass's own base temperature; a kelvin away, to 1e-4 where nothing feeds (the conduction is
    nearly linear). While a skin feeds, the feed's discrete steps make the law good to about 1 % per kelvin (measured
    while planning) -- which is what the interface mismatch and the step's repeat are for."""
    for case, d, tol in ((dict(), 0.0, 1e-9), (dict(q=(2e7, 1e7), T_top=900.0, T_bottom=820.0, T_b0=(815.0, 820.0)), 1.0, 1e-4)):
        f, q, eps, L, room, avail, T_b0 = heated_field(**case)
        a, b = f.advance(0.5, q, 0.0, eps, T_b0, L, room, avail)
        assert (b <= 0.0).all()                              # a warmer base takes no more heat from the skin
        book = f.advance(0.5, q, 0.0, eps, T_b0, L, room, avail, T_b1=T_b0 + d)
        assert np.allclose(book.E_base, (a + b * (T_b0 + d)) * f.area * 0.5, rtol=tol)


def test_the_real_pass_closes_its_books():
    f, q, eps, L, room, avail, T_b0 = heated_field()
    before = sum(node_total(f, i, L[i]) for i in range(2))
    M0 = f.mass()
    book = f.advance(0.5, q, 0.0, eps, T_b0, L, room, avail, T_b1=T_b0 + np.array([3.0, 1.0]))
    after = sum(node_total(f, i, book.L[i]) for i in range(2))
    assert after - before == pytest.approx(float((book.E_top - book.E_rad - book.E_base + book.drawn_H).sum()), rel=1e-12)
    assert f.mass() + float((book.to_film + book.to_deep).sum()) == pytest.approx(M0 + float(book.drawn.sum()), rel=1e-13)
    assert np.allclose(book.L, L + book.to_film + book.to_deep, rtol=1e-13)
    assert (book.to_film > 0.0).all() and book.to_deep[1] > 0.0 and (book.drawn > 0.0).all()
```

- [ ] **Step 2: Run the tests to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`
Expected: the three new tests FAIL with `AttributeError: 'SkinField' object has no attribute 'advance'`.

- [ ] **Step 3: Append `advance` to `SkinField`**

```python
    # -- the macro step ---------------------------------------------------------------------------------------------
    def advance(self, dt, q_top, T_amb, emissivity, T_b0, L, room, avail, T_b1=None):
        """The macro step's sub-steps (spec section 6). T_b1=None: the trial pass, on a copy, base held at T_b0, returning
        the base law (a, b) -- q_b(T) = a + b T [W/m2] of heat into the 3D model -- from the pass's base energy and its
        derivative with respect to the base temperature. Otherwise the real pass with the base held at T_b1, the 3D
        model's end-of-step facet temperature (what its backward-Euler step assumed for the boundary): feeds, draws and
        returns a SkinBook."""
        n = self.n
        N = max(1, int(math.ceil(dt / self.settings.substep - 1e-9)))
        h = dt / N
        trial = T_b1 is None
        saved = self.state() if trial else None
        T_b = np.asarray(T_b0 if trial else T_b1, dtype=float) * np.ones(n)
        L = np.array(L, dtype=float)
        avail = np.array(avail, dtype=float)
        book = SkinBook(*(np.zeros(n) for _ in range(7)), L)
        sens = (np.zeros((n, self.nc)), np.zeros(n)) if trial else None
        for _ in range(N):
            E_top, E_rad, E_base, sens = self.conduct(h, q_top, T_amb, emissivity, T_b, L, sens=sens)
            to_film, to_deep, drawn, drawn_H, L, S = self.feed_and_draw(L, room, avail, T_b, None if sens is None else sens[0])
            if sens is not None:
                sens = (S, sens[1])
            avail = avail - drawn
            book.E_top += E_top
            book.E_rad += E_rad
            book.E_base += E_base
            book.to_film += to_film
            book.to_deep += to_deep
            book.drawn += drawn
            book.drawn_H += drawn_H
        book.L = L
        if trial:
            self.set_state(saved)
            # the slope is clipped at zero: near the liquidus a warmer base can feed more and pass more heat down, and a
            # law rising with the base temperature would cost the 3D matrix its diagonal dominance; what the clip
            # misses is the interface mismatch's to book
            b = np.minimum(sens[1] / (self.area * dt), 0.0)
            return book.E_base / (self.area * dt) - b * T_b, b
        return book
```

- [ ] **Step 4: Run the tests, then the unit tier**

Run: `"$PY" -m pytest tests/test_reentry_model_skin.py -q`, then the unit tier.
Expected: 17 passed (about 8 s); the unit tier's counts are the baseline's plus these 17.

---

### Task 6: The thermal backends — the interface flux, no radiation under skins, the state round trip

**Files:**
- Modify: `reentry_model/thermal/__init__.py`, `reentry_model/thermal/skfem_backend.py`,
  `reentry_model/thermal/fenicsx_backend.py`
- Test: `tests/test_reentry_model_thermal.py`, `tests/test_reentry_model_fenicsx.py` (append)

**Interfaces:**
- Produces: `StepResult.Q_interface: float = 0.0` [W]; `step(dt, q_conv, T_amb, dirichlet=None, nodal_load=None,
  interface=None)` with `interface = (mask, a, b)` per facet — on masked facets the heat flux a + b T_f [W/m²] enters the
  body implicitly (T_f the facet-mean temperature, its Jacobian on the facet's node pairs like the radiation's) and the
  facet does not radiate; `None` or an all-False mask is the old step bit for bit (the scalar emissivity is then used
  exactly as before); `radiated_power` leaves out the masked facets of the last step (`_radiating`, reset whenever the
  surface is rebuilt); `state()` and `set_state(state)` restore the temperatures, the predictor's previous field and
  the solve counter, so that a step can be repeated.

- [ ] **Step 1: Append the failing tests**

To `tests/test_reentry_model_thermal.py`:

```python
def _interface_solver(mesh_, name="skfem"):
    s = thermal.thermal_solver(name)
    s.setup(mesh_, material.Material.from_drama_json("AA7075_scheil"), 0.3)
    s.set_temperature(700.0)
    return s


def test_interface_flux_replaces_heating_and_radiation_on_masked_facets(coarse_sphere_mesh):
    """Step 3 skins (sub-plan 18): q = a + b T_f enters implicitly on the masked facets, which do not radiate; the
    discrete balance holds; an all-False mask is the plain step bit for bit. (The AMG set-up draws on numpy's global
    generator, so each compared step is seeded alike, as the model's --seed does for a run.)"""
    s0, s1 = _interface_solver(coarse_sphere_mesh), _interface_solver(coarse_sphere_mesh)
    nf = len(s0.areas)
    q = np.full(nf, 2e5)
    np.random.seed(1)
    r0 = s0.step(0.5, q, 0.0)
    np.random.seed(1)
    r1 = s1.step(0.5, q, 0.0, interface=(np.zeros(nf, dtype=bool), np.zeros(nf), np.zeros(nf)))
    assert np.array_equal(r0.T, r1.T) and r1.Q_interface == 0.0 and r0.Q_rad == r1.Q_rad
    s2 = _interface_solver(coarse_sphere_mesh)
    mask = np.arange(nf) % 2 == 0
    a, b = np.full(nf, 3e5), np.full(nf, -2e3)
    E0 = s2.energy()
    r2 = s2.step(0.5, np.where(mask, 0.0, q), 0.0, interface=(mask, a, b))
    Tf = s2.facet_temperature()
    assert r2.Q_interface == pytest.approx(float(((a + b * Tf) * s2.areas)[mask].sum()), rel=1e-12)
    assert s2.energy() - E0 == pytest.approx((r2.Q_conv + r2.Q_interface - r2.Q_rad) * 0.5, rel=1e-6)
    assert r2.Q_rad == pytest.approx(float((0.3 * SIGMA_SB * Tf ** 4 * s2.areas)[~mask].sum()), rel=1e-12)


def test_state_round_trip_repeats_a_step(coarse_sphere_mesh):
    s = _interface_solver(coarse_sphere_mesh)
    q = np.full(len(s.areas), 2e5)
    s.step(0.5, q, 0.0)
    snap = s.state()
    r1 = s.step(0.5, q, 0.0)
    s.set_state(snap)
    r2 = s.step(0.5, q, 0.0)
    assert np.allclose(r1.T, r2.T, rtol=0.0, atol=1e-6)
```

To `tests/test_reentry_model_fenicsx.py`:

```python
def test_backends_agree_with_an_interface_flux(coarse_sphere_mesh):
    """Step 3 skins (sub-plan 18): the same affine interface flux and radiation mask in both backends."""
    from reentry_model import material, thermal
    out = []
    for name in ("skfem", "fenicsx"):
        s = thermal.thermal_solver(name)
        s.setup(coarse_sphere_mesh, material.Material.from_drama_json("AA7075_scheil"), 0.3)
        s.set_temperature(700.0)
        nf = len(s.areas)
        mask = np.arange(nf) % 2 == 0
        r = s.step(0.5, np.where(mask, 0.0, 2e5), 0.0, interface=(mask, np.full(nf, 3e5), np.full(nf, -2e3)))
        out.append((r, s.temperature()))
    (ra, Ta), (rb, Tb) = out
    assert np.allclose(Ta, Tb, rtol=0.0, atol=1e-6) and ra.Q_interface == pytest.approx(rb.Q_interface, rel=1e-9)
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q -k "masked or round_trip"`
Expected: FAIL with `TypeError: ... step() got an unexpected keyword argument 'interface'` and
`AttributeError: ... has no attribute 'state'`.

- [ ] **Step 3: Apply the three diff blocks**

```diff
--- a/reentry_model/thermal/__init__.py
+++ b/reentry_model/thermal/__init__.py
@@ -24,12 +24,13 @@
     iterations: int               # Newton/Picard iterations taken
     Q_extra: float = 0.0          # W, nodal load applied over the step (Step 3 deferred melt energy)
     Q_dropped: float = 0.0        # W, load that fell on pinned (material-free) nodes and was not applied
+    Q_interface: float = 0.0      # W, heat the skins passed in through their facets (Step 3 sub-plan 18)
 
 
 class ThermalSolver(Protocol):
     def setup(self, mesh, material, emissivity) -> None: ...
     def set_temperature(self, T) -> None: ...                       # uniform float or nodal array
-    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None) -> StepResult: ...
+    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None) -> StepResult: ...
     def temperature(self) -> np.ndarray: ...
     def energy(self) -> float: ...                                  # stored enthalpy above material.T_REF [J]
     def element_energies(self, T=None) -> np.ndarray: ...           # phi_e rho V_e h(T_e) per element [J]
@@ -38,6 +39,8 @@
     def set_film_mass(self, mass) -> None: ...                      # nodal mass of the melt film [kg], Step 3
     def nodal_capacity(self): ...                                   # J/K per node, material + film (Step 3)
     def film_weight(self): ...                                      # the film's share of each node's mass, 0-1 (Step 3)
+    def state(self): ...                                            # what set_state needs to repeat a step (sub-plan 18)
+    def set_state(self, state) -> None: ...
 
 
 def thermal_solver(name, **options):
```

```diff
--- a/reentry_model/thermal/skfem_backend.py
+++ b/reentry_model/thermal/skfem_backend.py
@@ -104,6 +104,7 @@
         n = len(self.points)
         self.facet_pattern = _Pattern(np.repeat(self.faces, 3, axis=1).ravel(), np.tile(self.faces, (1, 3)).ravel(), n)
         self.Bf = np.ones((3, 3))[None] * (self.areas / 9.0)[:, None, None]
+        self._radiating = None                            # the facets changed: last step's skin mask no longer applies
         self._update_pinned()
 
     def _update_pinned(self):
@@ -146,6 +147,15 @@
     def temperature(self):
         return self.T.copy()
 
+    def state(self):
+        """What set_state needs to repeat a step: the temperatures, the predictor's previous field, the solve counter."""
+        return (self.T.copy(), None if self._T_prev is None else self._T_prev.copy(), self._solves)
+
+    def set_state(self, state):
+        self.T = state[0].copy()
+        self._T_prev = None if state[1] is None else state[1].copy()
+        self._solves = state[2]
+
     def facet_temperature(self, T=None):
         return (self.T if T is None else T)[self.faces].mean(axis=1)
 
@@ -195,7 +205,11 @@
 
     def radiated_power(self, T_amb, T=None):
         Tf = self.facet_temperature(T)
-        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())
+        q = self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas
+        radiating = getattr(self, "_radiating", None)
+        if radiating is not None and radiating.shape == q.shape:
+            q = np.where(radiating, q, 0.0)                    # a facet under a skin radiates from the skin's top
+        return float(q.sum())
 
     def film_weight(self):
         """Share of each node's mass that is melt film, 0 to 1. The film is liquid: its enthalpy has no latent
@@ -227,9 +241,11 @@
         h = self.material.enthalpy(self.T if T is None else T)
         return self.phi * self.material.rho * self.vol * h[self.tets].mean(axis=1)
 
-    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None):
+    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None):
         """One backward-Euler step with the per-facet convective flux q_conv [W/m2] and an optional nodal load vector
-        [W] (Step 3: the deferred melt energy; loads on pinned nodes are dropped and reported in StepResult.Q_dropped)."""
+        [W] (Step 3: the deferred melt energy; loads on pinned nodes are dropped and reported in StepResult.Q_dropped).
+        `interface = (mask, a, b)` per facet: the skins' base law, heat a + b T_f [W/m2] into the body through the
+        masked facets, which neither take q_conv's place nor radiate (Step 3 sub-plan 18)."""
         T_old = self.T
         # predictor: extrapolate the previous step (saves ~1 Newton iteration per step); plain T_old on the first step
         T_k = T_old + (T_old - self._T_prev) if self._T_prev is not None and dirichlet is None else T_old.copy()
@@ -241,15 +257,28 @@
         last_relative_change, previous_change, damping, decreases = None, None, 1.0, 0
         mat = self.material
         film_w = self.film_weight() if self.film_mass.any() else np.zeros(len(self.points))
+        # the skins' base law (Step 3 sub-plan 18): on the facets under a skin the heat flux a + b T_f enters
+        # implicitly, linearised like the radiation, and the facet does not radiate (its skin's top does)
+        mask = None
+        if interface is not None and np.any(interface[0]):
+            mask = np.asarray(interface[0], dtype=bool)
+            a_int = np.where(mask, np.asarray(interface[1], dtype=float), 0.0)
+            b_int = np.where(mask, np.asarray(interface[2], dtype=float), 0.0)
+        eps_f = self.emissivity if mask is None else np.where(mask, 0.0, self.emissivity)
+        self._radiating = None if mask is None else ~mask
         for iteration in range(1, self.max_iterations + 1):
             K, M, E = self.operators(T_k, T_old)
             Tf = self.facet_temperature(T_k)
-            B = self.facet_pattern.assemble(self.Bf * (4.0 * self.emissivity * SIGMA_SB * Tf ** 3)[:, None, None])
-            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
+            B = self.facet_pattern.assemble(self.Bf * (4.0 * eps_f * SIGMA_SB * Tf ** 3)[:, None, None])
+            F_rad = self.facet_load(eps_f * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
+            F_int = 0.0
+            if mask is not None:
+                B = B + self.facet_pattern.assemble(self.Bf * (-b_int)[:, None, None])
+                F_int = self.facet_load(a_int + b_int * Tf)
             # Newton on R(T) = E(T)/dt + K T - F_conv + F_rad(T) - F_extra with the tangent M/dt + K + B:
             # A T_new = A T_k - R(T_k); with a smooth c_p (c_tan = c_sec) this is the Step 2 iteration unchanged
             A = (M / dt + K + B).tocsr()
-            b = M @ T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra
+            b = M @ T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra + F_int
             if self.pinned.any():                              # material-free nodes keep their temperature (identity rows)
                 A = A + sp.diags(self.pinned.astype(float), format="csr")
                 b[self.pinned] = T_old[self.pinned]
@@ -281,7 +310,9 @@
                 self.max_iterations, last_relative_change, self.newton_tol))
         self._T_prev, self.T = T_old, T_new
         self.last_damping = damping
-        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()), Q_dropped)
+        Q_int = float(((a_int + b_int * self.facet_temperature(T_new)) * self.areas).sum()) if mask is not None else 0.0
+        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()),
+                          Q_dropped, Q_int)
 
     def _solve_with_dirichlet(self, A, b, x0, dirichlet):
         if dirichlet is None:
```

```diff
--- a/reentry_model/thermal/fenicsx_backend.py
+++ b/reentry_model/thermal/fenicsx_backend.py
@@ -100,6 +100,7 @@
         if not np.array_equal(pinned, self.pinned):
             self._structure_changed = True
         self.pinned = pinned
+        self._radiating = None                            # the facets changed: last step's skin mask no longer applies
 
     def set_film_mass(self, mass):
         """Nodal mass [kg] of the melt film riding on the boundary (see the skfem backend)."""
@@ -122,6 +123,15 @@
     def temperature(self):
         return self.T.copy()
 
+    def state(self):
+        """What set_state needs to repeat a step: the temperatures, the predictor's previous field, the solve counter."""
+        return (self.T.copy(), None if self._T_prev is None else self._T_prev.copy(), self._solves)
+
+    def set_state(self, state):
+        self.T = state[0].copy()
+        self._T_prev = None if state[1] is None else state[1].copy()
+        self._solves = state[2]
+
     # -- nodal pieces (mesh numbering) -----------------------------------------------------------------------------
     def facet_load(self, q):
         return np.bincount(self.faces.ravel(), weights=np.repeat(q * self.areas / 3.0, 3), minlength=self.mesh.n_nodes)
@@ -159,7 +169,11 @@
 
     def radiated_power(self, T_amb, T=None):
         Tf = self.facet_temperature(self.T if T is None else T)
-        return float((self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas).sum())
+        q = self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4) * self.areas
+        radiating = getattr(self, "_radiating", None)
+        if radiating is not None and radiating.shape == q.shape:
+            q = np.where(radiating, q, 0.0)                    # a facet under a skin radiates from the skin's top
+        return float(q.sum())
 
     def _assemble_stiffness(self, T):
         from dolfinx.fem.petsc import assemble_matrix
@@ -176,7 +190,7 @@
         m.sum_duplicates()
         return PETSc.Mat().createAIJ(size=m.shape, csr=(m.indptr.astype(np.int32), m.indices.astype(np.int32), m.data), comm=self.msh.comm)
 
-    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None):
+    def step(self, dt, q_conv, T_amb, dirichlet=None, nodal_load=None, interface=None):
         from petsc4py import PETSc
         mat, n = self.material, self.mesh.n_nodes
         T_old = self.T
@@ -191,6 +205,15 @@
             fixed_values[np.asarray(dirichlet[0])] = dirichlet[1]
         vol = self.mesh.element_volumes()
         film_w = self.film_weight() if self.film_mass.any() else np.zeros(n)
+        # the skins' base law (Step 3 sub-plan 18): on the facets under a skin the heat flux a + b T_f enters
+        # implicitly, linearised like the radiation, and the facet does not radiate (its skin's top does)
+        mask = None
+        if interface is not None and np.any(interface[0]):
+            mask = np.asarray(interface[0], dtype=bool)
+            a_int = np.where(mask, np.asarray(interface[1], dtype=float), 0.0)
+            b_int = np.where(mask, np.asarray(interface[2], dtype=float), 0.0)
+        eps_f = self.emissivity if mask is None else np.where(mask, 0.0, self.emissivity)
+        self._radiating = None if mask is None else ~mask
         converged, last_relative_change, previous_change, damping, decreases = False, None, None, 1.0, 0
         T_new, iteration = T_k, 0
         for iteration in range(1, self.max_iterations + 1):
@@ -204,9 +227,13 @@
             M_tan = self.lumped(mat.rho * c_tan, mat.rho * c_liq)
             E = self.lumped(mat.rho * c_sec, mat.rho * c_sec_l) * dT
             Tf = self.facet_temperature(T_k)
-            B = sp.csr_matrix((np.repeat(4.0 * self.emissivity * SIGMA_SB * Tf ** 3 * self.areas / 9.0, 9), (self._facet_rows, self._facet_cols)), shape=(n, n))
-            F_rad = self.facet_load(self.emissivity * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
-            b = M_tan * T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra
+            B = sp.csr_matrix((np.repeat(4.0 * eps_f * SIGMA_SB * Tf ** 3 * self.areas / 9.0, 9), (self._facet_rows, self._facet_cols)), shape=(n, n))
+            F_rad = self.facet_load(eps_f * SIGMA_SB * (Tf ** 4 - T_amb ** 4))
+            F_int = 0.0
+            if mask is not None:
+                B = B + sp.csr_matrix((np.repeat(-b_int * self.areas / 9.0, 9), (self._facet_rows, self._facet_cols)), shape=(n, n))
+                F_int = self.facet_load(a_int + b_int * Tf)
+            b = M_tan * T_k / dt + B @ T_k - E / dt + F_conv - F_rad + F_extra + F_int
             # operator in dof numbering: K (dolfinx) + diag(M_tan/dt) + B
             A = self.K.copy()
             extra = sp.diags(M_tan / dt) + B
@@ -249,4 +276,6 @@
                 self.max_iterations, last_relative_change, self.newton_tol))
         self._T_prev, self.T = T_old, T_new
         self.last_damping = damping
-        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()), Q_dropped)
+        Q_int = float(((a_int + b_int * self.facet_temperature(T_new)) * self.areas).sum()) if mask is not None else 0.0
+        return StepResult(T_new.copy(), float(F_conv.sum()), self.radiated_power(T_amb), iteration, float(F_extra.sum()),
+                          Q_dropped, Q_int)
```

- [ ] **Step 4: Run the tests, the FEniCSx tier and the unit tier**

Run: `"$PY" -m pytest tests/test_reentry_model_thermal.py -q`; then
`FI_PROVIDER=tcp CC=/Users/ashajain/miniforge3/envs/fenicsx_env/bin/clang "$FX" -m pytest tests/test_reentry_model_fenicsx.py -q`;
then the unit tier.
Expected: the thermal file 18 passed (16 before); the FEniCSx tier 13 passed (about 95 s); the unit tier's counts the
previous task's plus 2.

---

### Task 7: The skins in the melting body

**Files:**
- Modify: `reentry_model/body.py`
- Create: `tests/test_reentry_model_melting_skin.py`
- Modify: `tests/test_reentry_model_fenicsx.py` (append one test)

**Interfaces:**
- Consumes: Tasks 2–6 (`skin.SkinField`, `skin.SkinSettings`, `skin.INTERFACE_TOLERANCE`, the backends' `interface`,
  `state`, `set_state`, `StepResult.Q_interface`).
- Produces: `body.SURFACE_MODEL_NAMES`; `MeltSettings.surface_model` (`"elements"` default, or `"skin"`) and
  `MeltSettings.skin` (`skin.SkinSettings` or `None` for its defaults); on `MeltingBody`: `skins` (`None` with
  `elements` or `--removal instant`), `skin_of_patch`, `patch_of_skin`, `skins_created`, `skin_drawn_mass`,
  `skin_handover_mass`, `skin_short_events`, `interface_repeats`, `last_mismatch`, `skin_wall`, `last_skin_liquid`,
  `last_skin_nonrigid`; the methods `_advance_with_skins`, `_advance_elements` (the old `advance`, unchanged),
  `_create_skins(T) -> int`, `_skin_available(patches)`, `_compose_depth(level, depth3d)`, `_defer_to_facets(energy)`,
  `_freeze_into_skins(solid)`, `_death_to_skins(dead, rest, h_e)`, `_hand_over_targets(old_surface, idx)`,
  `_hand_over_skins(old_surface)`, `_refresh_skin_index()`; and overrides of `surface_temperature` and
  `radiated_power` that read the skins' tops.

What the diff does, so a reviewer can check it against the spec: (a) `advance` dispatches; with skins it creates them,
evaluates the surface flow once before the passes (the feed's reach needs δ_m; the melt step reuses that evaluation),
runs the trial pass, the 3D step with the law on the skins' facets and no heating or radiation there, the real pass,
books the mismatch on the facets' nodes, and repeats once, re-linearised about the new facet temperature, if any skin's
mismatch exceeds `INTERFACE_TOLERANCE` of its energy throughput. (b) The melt step applies the skins' draws to their
elements (the element keeps what the drawn metal left, as a deferred load), stops the feed of elements that own a skin
patch, freezes liquid on skin patches into the skin's top cell, turns the molten cascade off, and sends a dying
element's remainder into the skins above it. (c) Liquid on a skin patch shares the skin's top cell temperature: its
enthalpy is h_l(T₀), energy booked on such a patch goes to the top node at once, and it is no longer on the 3D nodes.
(d) Because a booking moves a skin's top within the melt step, the step's cached liquid enthalpy `h_p` is refreshed
before the deep runoff, before the film and spray, after the runoff's booking inside `_film_and_spray`, and before each
death pass's `before` sum — each refresh only with skins, so the element path's arithmetic is untouched. (e) Skins are
handed on at deaths by the liquid's own rule (the dying element's exposed faces by projected area, the nearest patch for
an orphan), layer by layer; empty skins are dropped. (f) At consumption the skins' metal leaves with the liquid.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reentry_model_melting_skin.py`:

```python
"""MeltingBody with skins (Step 3 sub-plan 18): identical to the element model until the first skin; skins at the
solidus drawing from their elements; exact books; deaths handing skins on; the interface repeat; the edge cases of the
plan's Review Focus."""
import numpy as np
import pytest

from reentry_model import body, heating, material, mesh, skin, thermal
from test_reentry_model_coupled import MASS_100MM, simulator

pytest.importorskip("cantera")

DEV = skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"], substep=0.05)


def make_body(the_mesh, surface_model="skin", name="AA7075_scheil", solver="skfem", **settings):
    mat = material.Material.from_drama_json(name)
    m = mesh.VolumeMesh(the_mesh.points, the_mesh.tets, dict(the_mesh.params), the_mesh.element_layer.copy())
    return body.MeltingBody(m, mat, thermal.thermal_solver(solver), MASS_100MM,
                            settings=body.MeltSettings(surface_model=surface_model, skin=DEV, **settings))


def uniform(b, q):
    return heating.HeatingResult(np.full(b.surface.n_patches, q), q, 0.0, 0.0, 0.0)


def start_at(b, T):
    """A body held at T throughout, its books reckoned from here."""
    b.solver.set_temperature(T)
    b.energy0 = b.energy()


def run(b, q, steps, t0=0.0):
    t = t0
    for _ in range(steps):
        t += 0.5
        b.advance(t, 0.5, uniform(b, q))
    return t


def assert_books(b):
    assert b.mass(0.0) == pytest.approx(b.mass0 - b.removed_mass, rel=1e-12)
    assert abs(b.energy_balance_residual()) < 1e-8


def test_settings_validate_the_surface_model():
    with pytest.raises(ValueError):
        body.MeltSettings(surface_model="magic")
    assert body.MeltSettings().surface_model == "elements"


def test_identical_to_the_element_model_until_the_first_skin(coarse_sphere_mesh):
    """Seeded alike before every step (the AMG set-up draws on numpy's generator, as the model's --seed handles)."""
    be, bs = make_body(coarse_sphere_mesh, "elements"), make_body(coarse_sphere_mesh, "skin")
    assert be.skins is None and bs.skins is not None and bs.skins.n == 0
    t, compared = 0.0, 0
    while t < 120.0:
        t += 0.5
        np.random.seed(int(2 * t))
        be.advance(t, 0.5, uniform(be, 4e5))
        np.random.seed(int(2 * t))
        bs.advance(t, 0.5, uniform(bs, 4e5))
        if bs.skins.n:
            break
        assert np.array_equal(be.solver.temperature(), bs.solver.temperature())
        compared += 1
    assert compared >= 3 and bs.skins_created > 0


def test_skins_draw_from_their_elements_and_the_books_close(coarse_sphere_mesh):
    b = make_body(coarse_sphere_mesh)
    start_at(b, 740.0)
    run(b, 1.5e6, 12)
    assert b.skins.n > 0 and b.skin_drawn_mass > 0.0 and (b.m_f + b.m_d).sum() > 0.0
    assert_books(b)
    assert (b.skins.thickness() >= 0.0).all()
    assert (b.patch_of_skin >= 0).all() and int((b.skin_of_patch >= 0).sum()) == b.skins.n
    assert b.solver.film_mass[np.unique(b.surface.faces[b.skin_of_patch >= 0])].sum() == pytest.approx(
        b._film_nodal()[np.unique(b.surface.faces[b.skin_of_patch >= 0])].sum())


def test_deaths_hand_skins_to_the_faces_they_expose(coarse_sphere_mesh):
    b = make_body(coarse_sphere_mesh)
    start_at(b, 860.0)
    n0, t = b.mesh.n_active, 0.0
    while b.mesh.n_active == n0 and t < 30.0:
        t = run(b, 3e6, 1, t)
    assert b.mesh.n_active < n0 and b.skin_handover_mass > 0.0
    assert_books(b)
    assert (b.skins.m.sum(axis=1) > 0.0).all() and (b.patch_of_skin >= 0).all()


def test_a_large_interface_mismatch_repeats_the_step(coarse_sphere_mesh, monkeypatch):
    monkeypatch.setattr(skin, "INTERFACE_TOLERANCE", 0.0)
    b = make_body(coarse_sphere_mesh)
    start_at(b, 760.0)
    run(b, 1e6, 4)
    assert b.skins.n > 0 and b.interface_repeats >= 3
    assert_books(b)


def test_liquid_on_a_cold_patch_moves_onto_its_new_skin_with_its_energy(coarse_sphere_mesh):
    b = make_body(coarse_sphere_mesh)
    b.m_f[5] = 1e-6
    b.solver.set_film_mass(b._film_nodal())
    E, M = b.energy() + b.pending_load.sum(), b.mass(0.0)
    assert b._create_skins(b.solver.temperature()) == 1 and b.skin_of_patch[5] >= 0
    assert b.energy() + b.pending_load.sum() == pytest.approx(E, rel=1e-12) and b.mass(0.0) == pytest.approx(M, rel=1e-14)
    assert b.solver.film_mass.sum() == 0.0                      # the liquid now rides the skin's top, not the 3D nodes


def test_a_skin_its_element_cannot_supply_is_made_thin_and_counted(coarse_sphere_mesh):
    b = make_body(coarse_sphere_mesh)
    e = b.surface.owner[0]
    b.phi[e] = body.PHI_DEATH + 0.01
    b.solver.set_fractions(b.phi)
    b.m_f[0] = 1e-9
    b.solver.set_film_mass(b._film_nodal())
    b.mass0, b.energy0 = b.mass(0.0), b.energy()
    b._create_skins(b.solver.temperature())
    r = b.skin_of_patch[0]
    assert r >= 0 and 0.0 < b.skins.thickness()[r] < 0.5 * DEV.thickness
    assert b.phi[e] >= body.PHI_DEATH - 1e-12
    b.advance(0.5, 0.5, uniform(b, 1e5))
    assert b.skin_short_events >= 1                             # it drew nothing: its element had nothing to give
    assert b.phi[e] >= body.PHI_DEATH - 1e-12 or not b.mesh.active[e]   # ... and died at PHI_DEATH if it reached it
    assert_books(b)


def test_coupled_flight_steps_with_skins_keep_the_books(coarse_sphere_mesh):
    """Physics-mode loads at 71 km, as test_film_spraying_death_and_balances, with skins: film forms and is sprayed."""
    b = make_body(coarse_sphere_mesh)
    sim = simulator(b, t_max=60.0)
    sim.advance(43.0)
    start_at(b, 880.0)
    model = heating.PhysicsHeating()
    for _ in range(8):
        sim.advance(0.5)
        a = sim.aero_state(sim.t, sim.y[:3], sim.y[3:])
        loads = model.evaluate(a, b.theta, b.surface_temperature(), 0.05, T_mean=b.mean_temperature())
        b.advance(sim.t, 0.5, loads, state=a)
    assert b.skins.n > 0 and b.sprayed_mass > 0.0
    assert_books(b)
```

Append to `tests/test_reentry_model_fenicsx.py`:

```python
def test_skins_match_the_skfem_backend(coarse_sphere_mesh):
    """The melt-layer skins (sub-plan 18) in both backends: the body at 800 K under a uniform 1.5 MW/m2 for six steps
    of 0.5 s, every patch skinned at the first step. Skins, mass, film and temperatures agree and both energy balances
    are exact."""
    pytest.importorskip("cantera")
    from reentry_model import body, heating, mesh, skin
    from test_reentry_model_coupled import MASS_100MM
    out = {}
    for name in ("skfem", "fenicsx"):
        np.random.seed(12345)
        m = mesh.VolumeMesh(coarse_sphere_mesh.points, coarse_sphere_mesh.tets, dict(coarse_sphere_mesh.params))
        b = body.MeltingBody(m, material.Material.from_drama_json("AA7075_scheil"), thermal.thermal_solver(name), MASS_100MM,
                             settings=body.MeltSettings(surface_model="skin",
                                                        skin=skin.SkinSettings(cell=skin.SKIN_PRESETS["dev"], substep=0.05)))
        b.solver.set_temperature(800.0)
        b.energy0 = b.energy()
        loads = heating.HeatingResult(np.full(b.surface.n_patches, 1.5e6), 1.5e6, 0.0, 0.0, 0.0)
        for k in range(6):
            b.advance(0.5 * (k + 1), 0.5, loads)
        out[name] = b
    s, f = out["skfem"], out["fenicsx"]
    assert s.skins.n == f.skins.n == s.surface.n_patches and np.array_equal(s.mesh.active, f.mesh.active)
    assert f.mass(0.0) == pytest.approx(s.mass(0.0), rel=1e-8) and f.m_f.sum() == pytest.approx(s.m_f.sum(), rel=1e-4)
    assert np.abs(f.solver.temperature() - s.solver.temperature()).max() < 1e-3
    assert np.abs(f.skins.T - s.skins.T).max() < 1e-3
    assert abs(s.energy_balance_residual()) < 1e-8 and abs(f.energy_balance_residual()) < 1e-8
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_melting_skin.py -q`
Expected: FAIL with `TypeError: MeltSettings.__init__() got an unexpected keyword argument 'surface_model'`.

- [ ] **Step 3: Apply the diff block**

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -7,7 +7,9 @@
 
 import numpy as np
 
+from .thermal import SIGMA_SB
 
+
 class Body(Protocol):
     def mass(self, t) -> float: ...
     def advance(self, t, dt, loads) -> None: ...          # loads: heating.HeatingResult applied over [t - dt, t]
@@ -148,6 +150,7 @@
                                  # next step, as every exposed element did before the cascade (amendment of 2026-10-03)
 NONRIGID_MAX_CROSSINGS = 64      # elements the non-rigid depth's line may cross before it stops (amendment of 2026-10-06)
 REMOVAL_NAMES = ("girin", "instant")
+SURFACE_MODEL_NAMES = ("elements", "skin")    # skin: a one-dimensional skin under every melting patch (sub-plan 18)
 
 
 @dataclass
@@ -166,8 +169,13 @@
     rigid_substrate: bool = True      # the thin branch needs a rigid substrate: the regime test reads the film plus the slurry
                                       # (above 50 % liquid) beneath it, and off Girin's closure a film on slurry does not
                                       # spray (amendment of 2026-10-06)
+    surface_model: str = "elements"   # skin: a one-dimensional skin of fine cells under every melting patch resolves the
+                                      # melt layer and feeds the film continuously (sub-plan 18); elements: as before
+    skin: object = None               # skin.SkinSettings for surface_model "skin"; None: its defaults
 
     def __post_init__(self):
+        if self.surface_model not in SURFACE_MODEL_NAMES:
+            raise ValueError("surface_model must be one of {}, got {!r}".format(SURFACE_MODEL_NAMES, self.surface_model))
         if self.removal not in REMOVAL_NAMES:
             raise ValueError("removal must be one of {}, got {!r}".format(REMOVAL_NAMES, self.removal))
         if self.size_feedback not in SIZE_FEEDBACK_NAMES:
@@ -269,6 +277,16 @@
         self.consumed = False
         self._refresh_geometry()
         self.newtonian_drag_sphere = self.newtonian_drag()[0]      # the meshed sphere's own value: the normalisation
+        self.skins = None                                          # sub-plan 18: the skins, with --surface-model skin
+        if self.settings.surface_model == "skin" and self.settings.removal != "instant":
+            from . import skin as skin_mod
+            self.skins = skin_mod.SkinField(material, self.settings.skin or skin_mod.SkinSettings())
+        self._draw_m, self._draw_H = np.zeros(mesh.n_elements), np.zeros(mesh.n_elements)   # metal the skins drew this step
+        self.skins_created = self.skin_short_events = self.interface_repeats = 0
+        self.skin_drawn_mass = self.skin_handover_mass = self.last_mismatch = self.skin_wall = 0.0
+        self.last_skin_liquid = self.last_skin_nonrigid = None
+        self._step_flow = None
+        self._refresh_skin_index()
 
     # -- geometry -----------------------------------------------------------------------------------------------
     def _refresh_geometry(self):
@@ -285,6 +303,7 @@
         self.patch_of_face = np.full(len(self.mesh._face_nodes), -1, dtype=np.int64)     # face id -> patch index (-1: not a patch)
         self.patch_of_face[self.surface.face_ids] = np.arange(self.surface.n_patches)
         self._fit_nose()
+        self._refresh_skin_index()
 
     def _fit_nose(self):
         """Mass centre, transverse radius and the windward-cap radius of the current surface (module docstring)."""
@@ -350,7 +369,8 @@
         return 2.0 * self.equivalent_radius()
 
     def mass(self, t):
-        return float((self.phi * self.element_mass).sum() + self.m_f.sum() + self.m_d.sum())
+        m = float((self.phi * self.element_mass).sum() + self.m_f.sum() + self.m_d.sum())
+        return m + self.skins.mass() if getattr(self, "skins", None) is not None else m
 
     def reference_area(self):
         """Drag reference area: the hull's projected area, paired with the C_D of drag_shape_factor (the raw surface's
@@ -363,7 +383,8 @@
         return (3.0 * self.mass(0.0) / (4.0 * math.pi * self.material.rho)) ** (1.0 / 3.0)
 
     def energy(self):
-        return self.solver.energy() + self.film_energy()
+        E = self.solver.energy() + self.film_energy()
+        return E + self.skins.energy() if getattr(self, "skins", None) is not None else E
 
     def mean_temperature(self):
         m = self.mass(0.0)
@@ -380,6 +401,87 @@
 
     # -- the step --------------------------------------------------------------------------------------------------
     def advance(self, t, dt, loads, state=None):
+        if self.skins is None:
+            return self._advance_elements(t, dt, loads, state)
+        return self._advance_with_skins(t, dt, loads, state)
+
+    def _advance_with_skins(self, t, dt, loads, state):
+        """The macro step with skins (sub-plan 18; spec sections 4 and 6): create skins; the trial pass gives each skin's
+        base law; the 3D step takes that law on the skins' facets (no heating, no radiation there) with the deferred
+        loads; the real pass, base held at the new facet temperature, feeds, draws and books; the mismatch goes to the 3D
+        nodes, and a step whose mismatch exceeds the tolerance repeats once, re-linearised about that temperature. Then
+        the rest of the melt step. Until the first skin exists this is the element model's step exactly."""
+        import dataclasses
+        import time
+        from . import skin as skin_mod
+        from . import spray as spray_mod
+        started = time.perf_counter()
+        mat, liq, sk = self.material, self.liquid, self.skins
+        self._create_skins(self.solver.temperature())
+        if not sk.n:
+            self._step_flow = None
+            return self._advance_elements(t, dt, loads, state)
+        on = self.skin_of_patch >= 0
+        p = self.patch_of_skin
+        flow = delta_m = None
+        if state is not None:                     # evaluated before the passes: the feed's reach needs delta_m
+            flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
+            delta_m, _ = spray_mod.melt_layer(flow, liq)
+        self._step_flow = (flow, delta_m)
+        dm_p = delta_m[p] if delta_m is not None else np.full(sk.n, np.nan)
+        room = np.where(np.isfinite(dm_p), mat.rho * sk.area * np.nan_to_num(dm_p), np.inf)
+        avail = self._skin_available(p)
+        L = (self.m_f + self.m_d)[p]
+        q_top = loads.q_conv[p]
+        eps = np.full(sk.n, self.emissivity)
+        T_b0 = self.surface.facet_mean(self.solver.temperature())[p]
+        q3 = np.where(on, 0.0, loads.q_conv)
+        capacity = LOAD_DT_MAX * self.solver.nodal_capacity()
+        applied = np.clip(self.pending_load, -capacity, capacity)
+        solver_state, skin_state = self.solver.state(), sk.state()
+        T_lin = T_b0
+        for attempt in (0, 1):
+            a, b = sk.advance(dt, q_top, self.T_ambient, eps, T_lin, L, room, avail)
+            A_f, B_f = np.zeros(self.surface.n_patches), np.zeros(self.surface.n_patches)
+            A_f[p], B_f[p] = a, b
+            res = self.solver.step(dt, q3, self.T_ambient, nodal_load=applied / dt, interface=(on, A_f, B_f))
+            T_b1 = self.surface.facet_mean(self.solver.temperature())[p]
+            book = sk.advance(dt, q_top, self.T_ambient, eps, T_b0, L, room, avail, T_b1=T_b1)
+            mismatch = book.E_base - (a + b * T_b1) * sk.area * dt
+            scale = np.abs(book.E_top) + np.abs(book.E_base) + 1e-30
+            if attempt == 0 and (np.abs(mismatch) > skin_mod.INTERFACE_TOLERANCE * scale).any():
+                self.solver.set_state(solver_state)
+                sk.set_state(skin_state)
+                T_lin = T_b1
+                self.interface_repeats += 1
+                continue
+            break
+        self.pending_load = self.pending_load - applied
+        full = np.zeros(self.surface.n_patches)
+        full[p] = mismatch
+        self._defer_to_facets(full)
+        self.last_mismatch = float(np.abs(mismatch).sum())
+        E_top, E_rad = float(book.E_top.sum()), float(book.E_rad.sum())
+        self.integrated_heat += res.Q_conv * dt + E_top
+        self.absorbed_heat += (res.Q_conv - res.Q_rad) * dt + E_top - E_rad
+        self.radiated_heat += res.Q_rad * dt + E_rad
+        self.iterations.append(res.iterations)
+        self.last = dataclasses.replace(res, Q_conv=res.Q_conv + E_top / dt, Q_rad=res.Q_rad + E_rad / dt)
+        self._applied_total = getattr(self, "_applied_total", 0.0) + res.Q_extra * dt
+        self._dropped_total = getattr(self, "_dropped_total", 0.0) + res.Q_dropped * dt
+        self.m_f[p] = self.m_f[p] + book.to_film
+        self.m_d[p] = self.m_d[p] + book.to_deep
+        e = self.surface.owner[p]
+        np.add.at(self._draw_m, e, book.drawn)
+        np.add.at(self._draw_H, e, book.drawn_H)
+        self.skin_drawn_mass += float(book.drawn.sum())
+        self.skin_short_events += int((sk.thickness() < 0.5 * sk.settings.thickness).sum())
+        if book.to_film.sum() + book.to_deep.sum() > 0.0 and self.melt_onset is None:
+            self.melt_onset = t
+        self.melt_step(t, dt, state)
+        self.skin_wall += time.perf_counter() - started
+
+    def _advance_elements(self, t, dt, loads, state=None):
         # A deferred melt load is energy the transferred mass delivered to the nodes it arrived at; a node that has
         # since melted away has nothing to heat with it, so no node is asked to move more than LOAD_DT_MAX in one
         # step and the remainder waits (it is applied as soon as the node has the capacity, dropped into Q_dropped
@@ -422,7 +524,12 @@
         # receded to it. Feeding from any depth let the interior drain through an intact skin -- by 82 s of the 100 mm
         # flight the core was being consumed faster than the outermost shell (measured 2026-09-22).
         flow = delta_m = None
-        if state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
+        if self.skins is not None and self._step_flow is not None:      # with skins: the step's own evaluation
+            flow, delta_m = self._step_flow
+            self._step_flow = None
+            if flow is not None:
+                self.last_face_ids = self.surface.face_ids.copy()
+        elif state is not None and s.removal != "instant" and self.surface.n_patches and self.m_f.size:
             flow = self.flow.evaluate(state, self.theta, self.nose_radius(), liq.rho, self.surface_temperature())
             delta_m, _ = spray_mod.melt_layer(flow, liq)
             self.last_face_ids = self.surface.face_ids.copy()      # the surface delta_m belongs to (the deaths come later)
@@ -431,6 +538,10 @@
         h_e = h_node[tets].mean(axis=1)                          # as solver.energy() weighs the solid
         h_liq = mat.enthalpy_liquid(T)
         h_p = self.film_enthalpy(h_liq)                          # as film_energy() weighs the film: liquid, with L_f
+        if self.skins is not None and self._draw_m.any():        # the metal the skins drew leaves its elements
+            self.phi = self.phi - self._draw_m / self.element_mass
+            self._defer_to_elements(self._draw_m * h_e - self._draw_H)     # the element keeps what the drawn metal left
+            self._draw_m[:], self._draw_H[:] = 0.0, 0.0
         # (i) feed and (iv) re-solidification, netted. The feed hands the molten fraction of what each element still
         # holds to the film (f_feed of phi_e, so an element that is a quarter molten hands over a quarter of its
         # remainder); the mirror rule hands back the fraction of the film that has fallen below the ramp. What stops
@@ -446,11 +557,19 @@
         reach = self._shear_reaches(delta_m) if s.feed_depth == "conjugate" and s.removal != "instant" else None
         if reach is not None:
             f = f * reach
+        if self.skins is not None and self.skins.n:              # a patch's skin feeds its film: its element does not
+            skin_owner = np.zeros(self.mesh.n_elements, dtype=bool)
+            skin_owner[self.surface.owner[self.skin_of_patch >= 0]] = True
+            f = np.where(skin_owner, 0.0, f)
         h_hot = (fn * h_node[tets]).mean(axis=1) * self.mesh.active     # what the molten part carries, per kg of element
         cap = np.where(self.owner_area > 0.0, self.phi, np.maximum(self.phi - PHI_MIN, 0.0))
         gross = np.minimum(f * self.phi, cap) * self.element_mass
         # the liquid on each patch below the ramp re-solidifies: the film and the deep liquid beneath it alike
         solid = (1.0 - mat.feed_fraction(self.film_temperature())) * (self.m_f + self.m_d) if s.removal != "instant" else np.zeros(0)
+        skin_solid = None
+        if self.skins is not None and self.skins.n and solid.size:      # on a skin it freezes into the skin's top cell
+            on = self.skin_of_patch >= 0
+            skin_solid, solid = np.where(on, solid, 0.0), np.where(on, 0.0, solid)
         want = np.bincount(self.surface.owner, solid, self.mesh.n_elements) if solid.size else np.zeros(self.mesh.n_elements)
         net = gross - want
         fed, want = np.maximum(net, 0.0), np.maximum(-net, 0.0)
@@ -476,16 +595,22 @@
             delta, carried = self._add_to_film(fed, fed_h)
             self._defer_to_patches(carried - delta * h_p)        # ... and the melt arrives at its patch's temperature
             if s.deep_runoff and flow is not None:
+                if self.skins is not None:                       # a skin's top moves when energy is booked on it
+                    h_p = self.film_enthalpy(h_liq)
                 deep_liquid = self._deep_runoff(dt, flow, delta_m, reach, T, h_node, h_e, h_p)
+            if self.skins is not None:
+                h_p = self.film_enthalpy(h_liq)
             released = self._film_and_spray(t, dt, state, h_p, flow, delta_m)
         frozen = 0.0 if s.removal == "instant" else self._freeze_back(want, solid, h_e, h_p)
+        if skin_solid is not None:
+            frozen += self._freeze_into_skins(skin_solid)
         # (v) death of consumed patch owners; a death exposes its neighbours, which die in turn if they are consumed --
         # and, with the molten cascade (amendment of 2026-10-03), are fed within the step if the death exposed them
         # fully molten (`_feed_exposed`), so that they die in the next pass and expose the next
         n_dead, passes, cascade, capped = 0, 0, 0.0, False
         # fully molten as `molten_depth` has it -- the element's mean nodal temperature at or above T_feed -- so the cascade
         # removes exactly the contiguous molten layer the branch test reads, wherever a death exposes it
-        full = T[tets].mean(axis=1) >= mat.T_feed if s.molten_cascade and s.removal != "instant" else None
+        full = T[tets].mean(axis=1) >= mat.T_feed if s.molten_cascade and s.removal != "instant" and self.skins is None else None
         while not self.consumed:
             dead = np.flatnonzero((self.owner_area > 0.0) & self.mesh.active & (self.phi <= PHI_DEATH))
             if dead.size == 0:
@@ -496,11 +621,15 @@
                 self.removed_enthalpy += float(rest.sum()) * mat.h_liquid
                 self._defer_to_elements(rest * (h_e[dead] - mat.h_liquid), dead)
             else:
+                # with skins, an element under a skin passes its remainder to that skin's base (spec section 9)
+                went = self._death_to_skins(dead, rest, h_e) if self.skins is not None and self.skins.n else np.zeros(dead.size, dtype=bool)
                 extra, extra_h = np.zeros(self.mesh.n_elements), np.zeros(self.mesh.n_elements)
-                extra[dead], extra_h[dead] = rest, rest * h_e[dead]
+                extra[dead[~went]], extra_h[dead[~went]] = rest[~went], rest[~went] * h_e[dead[~went]]
                 delta, carried = self._add_to_film(extra, extra_h)
                 self._defer_to_patches(carried - delta * h_p)
             self.phi[dead] = 0.0
+            if self.skins is not None:
+                h_p = self.film_enthalpy(h_liq)
             before = float(((self.m_f + self.m_d) * h_p).sum())
             was_owner = self.owner_area > 0.0
             self._kill(dead)                                     # the liquid of a dead patch moves to another patch...
@@ -607,8 +736,17 @@
             self._defer_to_patches(energy * weights / total)
 
     def _defer_to_patches(self, energy):
-        """Book a per-patch energy [J] on the patches' nodes, to be applied over the next step."""
+        """Book a per-patch energy [J]: on the nodes of a patch without a skin, to be applied over the next step; on a
+        skin's top node at once (sub-plan 18), shared with the liquid that node carries."""
         e = np.asarray(energy, dtype=float)
+        if not e.any():
+            return
+        if getattr(self, "skins", None) is not None and self.skins.n:
+            on = (self.skin_of_patch >= 0) & (e != 0.0)
+            if on.any():
+                p = np.flatnonzero(on)
+                self.skins.book_top(self.skin_of_patch[p], e[p], (self.m_f + self.m_d)[p])
+                e = np.where(self.skin_of_patch >= 0, 0.0, e)
         if e.any():
             np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))
 
@@ -618,7 +756,10 @@
         beneath the film rides the same nodes in the same way (amendment of 2026-10-02)."""
         if not self.m_f.size:
             return np.zeros(self.mesh.n_nodes)
-        return np.bincount(self.surface.faces.ravel(), np.repeat((self.m_f + self.m_d) / 3.0, 3), self.mesh.n_nodes)
+        lm = self.m_f + self.m_d
+        if getattr(self, "skins", None) is not None:              # liquid on a skin rides the skin's top, not the nodes
+            lm = np.where(self.skin_of_patch >= 0, 0.0, lm)
+        return np.bincount(self.surface.faces.ravel(), np.repeat(lm / 3.0, 3), self.mesh.n_nodes)
 
     def _deep_runoff(self, dt, flow, delta_m, reach, T, h_node, h_e, h_p):
         """(ii-a) The deep runoff (amendment of 2026-10-02): the contiguous liquid below the conjugate depth runs off
@@ -747,6 +888,11 @@
         # instead of the fully molten material alone. The Rayleigh-Taylor criteria and the molten region the wave-fits
         # test measures against keep the liquid.
         nonrigid = self.nonrigid_depth() if s.rigid_substrate else None
+        if self.skins is not None and self.skins.n:              # the skins resolve both depths (spec section 5.5)
+            molten = self._compose_depth(mat.T_feed, molten)
+            if nonrigid is not None:
+                nonrigid = self._compose_depth(mat.T_rigid, nonrigid)
+            self.last_skin_liquid, self.last_skin_nonrigid = molten, nonrigid
         under = molten if nonrigid is None else nonrigid
         # (ii) lubrication and runoff
         n_sub, moved = 0, 0.0
@@ -760,6 +906,8 @@
             # over the patches that gained film, which is exact in total and second order in the attribution)
             d = self.m_f - before
             self._spread_to_patches(-float((d * h_p).sum()), np.maximum(d, 0.0))
+            if self.skins is not None and self.skins.n:          # the booking moved skin tops: the droplets leave at theirs
+                h_p = self.film_enthalpy()
         b = self.m_f / (liq.rho * areas)
         layer = b + molten + deep
         regime_layer = layer if nonrigid is None else b + nonrigid + deep
@@ -858,13 +1006,20 @@
             self.consumed = True                                 # the liquid still on it leaves with the body
             self.removed_mass += float(self.m_f.sum() + self.m_d.sum())
             self.removed_enthalpy += float(((self.m_f + self.m_d) * self.film_enthalpy()).sum())
+            if self.skins is not None and self.skins.n:          # ... and so does the skins' metal
+                self.removed_mass += self.skins.mass()
+                self.removed_enthalpy += self.skins.energy()
+                self.skins.keep(np.zeros(0, dtype=np.int64))
             self.m_f, self.m_d = np.zeros(0), np.zeros(0)
+            self.skin_of_patch, self.patch_of_skin = np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
             return
         self._refresh_geometry()
         keep = self.patch_of_face[old_surface.face_ids]
         kept = keep >= 0
         self.m_f = self._hand_over(old_surface, old_m_f, keep, kept)
         self.m_d = self._hand_over(old_surface, old_m_d, keep, kept)   # the deep liquid beneath it goes the same way
+        if self.skins is not None and self.skins.n:
+            self._hand_over_skins(old_surface)
 
     def _hand_over(self, old_surface, old, keep, kept):
         """A per-patch liquid mass of the surface before the deaths, carried onto the current surface: kept patches keep
@@ -905,6 +1060,156 @@
                 np.add.at(m_f, self._patch_tree.query(old_surface.centroids[orphan])[1], old[orphan])
         return m_f
 
+    # -- skins (sub-plan 18) -----------------------------------------------------------------------------------------
+    def _refresh_skin_index(self):
+        """skin_of_patch: the skin row under each patch of the current surface (-1: none); patch_of_skin the reverse
+        (-1 for a skin whose patch has vanished, until the death loop hands it over)."""
+        self.skin_of_patch = np.full(self.surface.n_patches, -1, dtype=np.int64)
+        sk = getattr(self, "skins", None)
+        if sk is None or not sk.n:
+            self.patch_of_skin = np.zeros(0, dtype=np.int64)
+            return
+        self.patch_of_skin = self.patch_of_face[sk.face_id]
+        ok = self.patch_of_skin >= 0
+        self.skin_of_patch[self.patch_of_skin[ok]] = np.flatnonzero(ok)
+
+    def _skin_available(self, patches):
+        """Metal [kg] a skin on each of `patches` may take from its element: its share, by area among the element's
+        patches, of what the element holds above PHI_DEATH."""
+        e = self.surface.owner[patches]
+        share = self.surface.areas[patches] / self.owner_area[e]
+        return np.maximum(self.phi[e] - PHI_DEATH, 0.0) * self.element_mass[e] * share
+
+    def _create_skins(self, T):
+        """Patches whose surface has reached the solidus, or that hold liquid, get a skin (spec 5.2): a slab of the
+        target thickness -- less where the element cannot give that above PHI_DEATH -- taken from the patch's element,
+        its cells filled from the P1 field along the inward normal. The element's mass leaves at its own mean nodal
+        enthalpy and the difference to what the cells hold is booked on its nodes; liquid already on the patch moves
+        onto the skin's top with its energy. Returns how many skins were made."""
+        sk, s, mat = self.skins, self.surface, self.material
+        if sk is None or not s.n_patches:
+            return 0
+        Tf = s.facet_mean(T)
+        idx = np.flatnonzero((self.skin_of_patch < 0) & ((Tf >= mat.T_solidus) | (self.m_f + self.m_d > 0.0)))
+        if not idx.size:
+            return 0
+        thick = np.minimum(sk.settings.thickness, self._skin_available(idx) / (mat.rho * s.areas[idx]))
+        idx, thick = idx[thick > 0.0], thick[thick > 0.0]
+        if not idx.size:
+            return 0
+        e = s.owner[idx]
+        z = (np.arange(sk.nc) + 0.5)[None, :] * (thick / sk.nc)[:, None]
+        x = s.centroids[idx][:, None, :] - s.normals[idx][:, None, :] * z[:, :, None]
+        P = self.mesh.points[self.mesh.tets[e]]
+        inv = np.linalg.inv(np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]], axis=2))
+        lam = np.einsum("kij,kcj->kci", inv, x - P[:, None, 0])
+        lam = np.concatenate([1.0 - lam.sum(axis=2, keepdims=True), lam], axis=2)
+        lam = np.clip(lam, 0.0, None)
+        lam = lam / lam.sum(axis=2, keepdims=True)            # a point past the element takes its nearest face's value
+        Tp = np.einsum("kci,ki->kc", lam, T[self.mesh.tets[e]])
+        L = (self.m_f + self.m_d)[idx]
+        h_nodes = self.film_enthalpy()[idx]                   # the liquid's enthalpy on the 3D nodes, before the skins
+        mass, H = sk.add(s.face_ids[idx], s.areas[idx], Tp, thick)
+        h_e = mat.enthalpy(T)[self.mesh.tets[e]].mean(axis=1)
+        np.subtract.at(self.phi, e, mass / self.element_mass[e])
+        self._defer_to_elements(np.bincount(e, mass * h_e - H, self.mesh.n_elements))
+        self._refresh_skin_index()
+        if (L > 0.0).any():
+            rows = self.skin_of_patch[idx]
+            sk.book_top(rows, L * (h_nodes - mat.enthalpy_liquid(sk.T[rows, 0])), L)
+        self.solver.set_fractions(self.phi)
+        self.solver.set_film_mass(self._film_nodal())
+        self.skins_created += int(idx.size)
+        return int(idx.size)
+
+    def _compose_depth(self, level, depth3d):
+        """Per patch: the depth of the layer above `level` read from the skin where there is one, continued into the 3D
+        model's own depth (measured from the skin's base) when the whole skin is above it (spec 5.5)."""
+        out = np.array(depth3d, dtype=float, copy=True)
+        if self.skins is None or not self.skins.n:
+            return out
+        d, through = self.skins.depth_above(level)
+        on = self.skin_of_patch >= 0
+        r = self.skin_of_patch[on]
+        out[on] = d[r] + np.where(through[r], out[on], 0.0)
+        return out
+
+    def _defer_to_facets(self, energy):
+        """Book a per-patch energy [J] on the 3D nodes of the patches' facets, skins or not (the interface mismatch)."""
+        e = np.asarray(energy, dtype=float)
+        if e.any():
+            np.add.at(self.pending_load, self.surface.faces.ravel(), np.repeat(e / 3.0, 3))
+
+    def _freeze_into_skins(self, solid):
+        """Liquid on a skin patch whose top has fallen below the feed ramp refreezes into the skin's top cell, the node's
+        energy unchanged (spec 7.1 (6)); the deep liquid, against the solid, goes first. Returns the mass [kg]."""
+        p = np.flatnonzero(solid > 0.0)
+        if not p.size:
+            return 0.0
+        from_deep = np.minimum(self.m_d[p], solid[p])
+        from_film = np.minimum(self.m_f[p], solid[p] - from_deep)
+        taken = from_deep + from_film
+        self.m_d[p] = self.m_d[p] - from_deep
+        self.m_f[p] = self.m_f[p] - from_film
+        self.skins.refreeze_top(self.skin_of_patch[p], taken, (self.m_f + self.m_d)[p])
+        return float(taken.sum())
+
+    def _death_to_skins(self, dead, rest, h_e):
+        """Dying elements under skins pass their remainder `rest` [kg] into those skins' bottoms, shared by area among
+        the element's skin patches, at the element's own mean nodal enthalpy (spec 9). Returns which of `dead` did."""
+        on = self.skin_of_patch >= 0
+        p = np.flatnonzero(on & np.isin(self.surface.owner, dead))
+        went = np.isin(dead, self.surface.owner[p])
+        if p.size:
+            e = self.surface.owner[p]
+            skin_area = np.bincount(e, self.surface.areas[p], self.mesh.n_elements)
+            full = np.zeros(self.mesh.n_elements)
+            full[dead] = rest
+            mass = full[e] * self.surface.areas[p] / skin_area[e]
+            self.skins.add_bottom(self.skin_of_patch[p], mass, mass * h_e[e], (self.m_f + self.m_d)[p])
+        return went
+
+    def _hand_over_targets(self, old_surface, idx):
+        """For vanished patches `idx` of `old_surface`: (targets (k, 4), patch indices of the current surface or -1;
+        weights (k, 4)) -- the faces their own element exposed, weighted by area projected on the vanished patch's
+        normal (raw area where none faces that way), or the nearest patch where the element exposed nothing. The rule
+        of `_hand_over`, kept apart so that the liquid's hand-over keeps its exact arithmetic."""
+        targets = self.patch_of_face[self.mesh.element_faces(old_surface.owner[idx])]
+        ok = targets >= 0
+        safe = np.where(ok, targets, 0)
+        cosine = np.clip(np.einsum("nij,nj->ni", self.surface.normals[safe], old_surface.normals[idx]), 0.0, None)
+        area = np.where(ok, self.surface.areas[safe] * cosine, 0.0)
+        total = area.sum(axis=1)
+        flat = ok & (total[:, None] <= 0.0)
+        if flat.any():
+            area = np.where(flat, self.surface.areas[safe], area)
+            total = area.sum(axis=1)
+        has = total > 0.0
+        weights = np.where(has[:, None], area / np.where(has, total, 1.0)[:, None], 0.0)
+        targets = np.where(ok & has[:, None], targets, -1)
+        if (~has).any():
+            orphan = np.flatnonzero(~has)
+            targets[orphan, 0] = self._patch_tree.query(old_surface.centroids[idx[orphan]])[1]
+            weights[orphan, 0] = 1.0
+        return targets, weights
+
+    def _hand_over_skins(self, old_surface):
+        """Skins whose patch vanished pass to the faces the deaths exposed, layer by layer (spec 9); empty skins go."""
+        sk = self.skins
+        lost = np.flatnonzero(self.patch_of_face[sk.face_id] < 0)
+        if lost.size:
+            old_patch = np.full(len(self.patch_of_face), -1, dtype=np.int64)
+            old_patch[old_surface.face_ids] = np.arange(old_surface.n_patches)
+            targets, weights = self._hand_over_targets(old_surface, old_patch[sk.face_id[lost]])
+            tgt_faces = np.where(targets >= 0, self.surface.face_ids[np.maximum(targets, 0)], -1)
+            _, handed = sk.hand_over(lost, tgt_faces, weights, lambda f: self.surface.areas[self.patch_of_face[f]])
+            self.skin_handover_mass += handed
+        keep = np.ones(sk.n, dtype=bool)
+        keep[lost] = False
+        keep &= sk.m.sum(axis=1) > 0.0
+        sk.keep(np.flatnonzero(keep))
+        self._refresh_skin_index()
+
     # -- reporting -----------------------------------------------------------------------------------------------
     def on_current_surface(self, values, fill):
         """A per-patch array from the last flow and spray evaluation (`last_flow`, `last_spray`), carried onto the
@@ -920,7 +1225,22 @@
         index[self.last_face_ids] = np.arange(self.last_face_ids.size)
         j = index[self.surface.face_ids]
         return np.where(j >= 0, np.asarray(values, dtype=float)[np.maximum(j, 0)], fill)
+
+    def surface_temperature(self):
+        """Wall temperature per patch: the skin's top cell where a patch has a skin (sub-plan 18), else the facet's."""
+        Tf = self.surface.facet_mean(self.solver.temperature())
+        if getattr(self, "skins", None) is not None and self.skins.n:
+            on = self.skin_of_patch >= 0
+            Tf[on] = self.skins.T[self.skin_of_patch[on], 0]
+        return Tf
 
+    def radiated_power(self):
+        P = self.solver.radiated_power(self.T_ambient)
+        if getattr(self, "skins", None) is not None and self.skins.n:   # skins radiate from their tops, at patch area
+            T0 = self.skins.T[:, 0]
+            P += float((self.emissivity * SIGMA_SB * (T0 ** 4 - self.T_ambient ** 4) * self.skins.area).sum())
+        return P
+
     def film_temperature(self):
         """The film's temperature per patch: the surface's own (the film is thermally thin, class docstring)."""
         return self.surface_temperature()
@@ -933,7 +1253,13 @@
         coupled balance), and the mixture enthalpy h(T) (which makes melting free on the ramp)."""
         if h_node is None:
             h_node = self.material.enthalpy_liquid(self.solver.temperature())
-        return h_node[self.surface.faces].mean(axis=1) if self.m_f.size else np.zeros(0)
+        if not self.m_f.size:
+            return np.zeros(0)
+        h = h_node[self.surface.faces].mean(axis=1)
+        if getattr(self, "skins", None) is not None and self.skins.n:   # liquid on a skin shares its top cell's temperature
+            on = self.skin_of_patch >= 0
+            h[on] = self.material.enthalpy_liquid(self.skins.T[self.skin_of_patch[on], 0])
+        return h
 
     def film_energy(self):
         """Enthalpy held by the film [J]. Identical to the nodal sum the solver carries (`set_film_mass` gives each
@@ -1165,7 +1491,8 @@
 
     def demised(self):
         """True once the body's material (film excluded) is below the demise fraction of the initial mass."""
-        return self.consumed or float((self.phi * self.element_mass).sum()) < self.settings.demise_fraction * self.mass0
+        metal = float((self.phi * self.element_mass).sum()) + (self.skins.mass() if self.skins is not None else 0.0)
+        return self.consumed or metal < self.settings.demise_fraction * self.mass0
 
     def energy_balance_residual(self):
         """(E_body + E_film - E0 - absorbed heat + removed enthalpy + pending load + dropped load) / absorbed heat:
```

- [ ] **Step 4: Run the tests, the FEniCSx tier and the unit tier**

Run: `"$PY" -m pytest tests/test_reentry_model_melting_skin.py -q --durations=10`; then the FEniCSx tier; then the unit
tier.
Expected: 8 passed in about 60 s (the coupled flight steps about 22 s, the drawing test about 19 s); FEniCSx 14 passed;
the unit tier's counts the previous task's plus 8, every earlier test unchanged — every one of them runs
`surface_model="elements"`.

---

### Task 8: CLI flags, run names, history columns, results and frame fields

**Files:**
- Modify: `reentry_model/cli.py`, `reentry_model/coupled.py`, `reentry_model/body.py`
- Test: `tests/test_reentry_model_cli.py` (append)

**Interfaces:**
- Consumes: Task 7.
- Produces: CLI `--surface-model {elements,skin}` (default `elements` in stage 1), `--skin-thickness` [mm] (0.4),
  `--skin-cell` [µm] (default: the preset's), `--skin-substep` [ms] (10), `--skin-preset {production,dev}`;
  `cli.skin_settings(args)`; `cli.skin_name_suffix(thickness_mm, cell_um, substep_ms, preset) -> str`;
  `model_run_name(..., surface_model="elements", skin_suffix="")`; `coupled.SKIN_COLUMNS` (in `MELT_COLUMNS`) filled by
  `MeltingBody._skin_stats()` through `melt_stats()`; result fields `skins_created`, `skin_drawn_mass_kg`,
  `skin_handover_mass_kg`, `interface_repeats`, `skin_wall_s`; frame fields `skin_T_top`, `skin_liquid_depth`,
  `skin_nonrigid_depth`, `skin_thickness` (nan on a patch without a skin); the settings record `surface_model` and
  `skin`.

- [ ] **Step 1: Append the failing tests**

To `tests/test_reentry_model_cli.py`:

```python
# ---------------------------------------------------------------------------------------------------------------
# Step 3 sub-plan 18: the melt-layer skins' flags, run names, columns and frame fields

def test_skin_settings_reach_the_run_name():
    name = cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin")
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin",
                              surface_model="skin") == name + "_surface-skin"
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", "girin", surface_model="skin",
                              skin_suffix="_skin-dev", seed=7) == name + "_surface-skin_skin-dev_seed-7"
    assert cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table", "none", "physics", None,
                              surface_model="skin") == cli.model_run_name(0.1, 7500.0, 77500.133, "us76", "sesam-table",
                                                                          "none", "physics", None)   # no melt, no skins
    assert cli.skin_name_suffix(0.4, None, 10.0, "production") == "" and cli.skin_name_suffix(0.4, 10.0, 10.0, "production") == ""
    assert cli.skin_name_suffix(0.4, None, 10.0, "dev") == "_skin-dev"
    assert cli.skin_name_suffix(0.8, 5.0, 20.0, "production") == "_skin0.8mm_cell5um_sub20ms"


def test_a_skin_run_writes_its_columns_frames_and_settings(tmp_path):
    """3 s from 71 km with a warm body and skins (dev preset), the 0.5 s step throughout (`--dt-continuum off`, as the
    short melting run passes its own step): exit 0, the run name, the skin columns, the frame fields."""
    pytest.importorskip("cantera")
    import glob
    import pyvista as pv
    argv = [a for a in MELT if a not in ("--heating", "sesam")] + [
        "--heating", "physics", "--altitude", "71", "--velocity", "7.24", "--temperature", "800", "--t-max", "3",
        "--dt-continuum", "off", "--surface-model", "skin", "--skin-preset", "dev", "--frames-every", "2", "--outdir", str(tmp_path)]
    assert cli.main(argv) == 0
    (path,) = glob.glob(str(tmp_path / "*.json"))
    doc = json.load(open(path))
    assert doc["run_name"].endswith("_melt-girin_dtcontinuum-off_surface-skin_skin-dev")
    s, r = doc["settings"], doc["results"]
    assert s["surface_model"] == "skin" and s["skin"] == {"thickness_mm": 0.4, "cell_um": 20.0, "substep_ms": 10.0, "preset": "dev"}
    assert r["skins_created"] > 0 and r["skin_drawn_mass_kg"] >= 0.0 and abs(r["melt_energy_balance_residual"]) < 1e-6
    rows = list(csv.DictReader(open(path[:-5] + ".csv")))
    assert all(k in rows[0] for k in coupled.SKIN_COLUMNS) and float(rows[-1]["skin_count"]) > 0.0
    surf = sorted(glob.glob(str(tmp_path / "**" / "surface_*.vtp"), recursive=True))
    assert surf and all(k in pv.read(surf[-1]).cell_data for k in ("skin_T_top", "skin_liquid_depth", "skin_nonrigid_depth",
                                                                    "skin_thickness"))


@pytest.mark.parametrize("argv", [
    MELT + ["--surface-model", "magic"],
    MELT + ["--surface-model", "skin", "--skin-thickness", "0"],
    MELT + ["--surface-model", "skin", "--skin-substep", "-1"],
    MELT + ["--surface-model", "skin", "--skin-cell", "300"],          # more than half of the 0.4 mm skin
    MELT + ["--surface-model", "skin", "--skin-preset", "fine"],
])
def test_bad_skin_arguments_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q -k skin`
Expected: the run-name test FAILS (`model_run_name() got an unexpected keyword argument 'surface_model'`) and so does
the short run (argparse rejects `--surface-model`); the five bad-argument cases already pass, because argparse rejects the
unknown flags with exit 2 — once the flags exist they pin the validation instead.

- [ ] **Step 3: Apply the three diff blocks**

```diff
--- a/reentry_model/cli.py
+++ b/reentry_model/cli.py
@@ -19,7 +19,7 @@
 
 import numpy as np
 
-from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, spray, surface_flow, thermal, viz
+from . import __version__, aero, atmosphere, body, compare, coupled, dispersion, fap, heating, material, mesh, sesam_io, skin, spray, surface_flow, thermal, viz
 from . import trajectory as tj
 from .earth import GRAVITY_MODELS
 from .thermal import MissingBackend
@@ -90,13 +90,15 @@
 
 
 def model_run_name(diameter_m, velocity_ms, altitude_m, atmosphere_name, bridging_name, wind_name, heating_name=None, melt=None,
-                   deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED, rigid_substrate=True, dt_continuum="default"):
+                   deep_runoff=True, molten_cascade=True, seed=DEFAULT_SEED, rigid_substrate=True, dt_continuum="default",
+                   surface_model="elements", skin_suffix=""):
     """The run name encodes the configuration, so runs of different settings can share an output directory. A melting
     run with the deep runoff off (amendment of 2026-10-02) ends in `_deeprunoff-off`, one with the molten cascade off
     (amendment of 2026-10-03) in `_moltencascade-off`, one with the rigid substrate off (amendment of 2026-10-06) in
     `_rigidsubstrate-off`, one whose continuum step is not its removal mode's default (amendment of 2026-10-07) in
-    `_dtcontinuum-off` or `_dtcontinuum-<s>`, and any run with a seed other than DEFAULT_SEED (amendment of 2026-10-05)
-    in `_seed-<n>`; the defaults keep the old name."""
+    `_dtcontinuum-off` or `_dtcontinuum-<s>`, one with skins (sub-plan 18) in `_surface-skin` followed by `skin_suffix`
+    (the skins' non-default settings, `skin_name_suffix`), and any run with a seed other than DEFAULT_SEED (amendment of
+    2026-10-05) in `_seed-<n>`; the defaults keep the old name."""
     dtc_default = DEFAULT_DT_CONTINUUM if melt == "girin" else None
     dtc = dtc_default if dt_continuum == "default" else dt_continuum
     name = "model_d{:06.2f}mm_v{:08.5f}kms_h{:07.3f}km_{}_{}_{}".format(
@@ -105,9 +107,24 @@
             + ("_deeprunoff-off" if melt and not deep_runoff else "") + ("_moltencascade-off" if melt and not molten_cascade else "")
             + ("_rigidsubstrate-off" if melt and not rigid_substrate else "")
             + ("" if not melt or dtc == dtc_default else "_dtcontinuum-off" if dtc is None else "_dtcontinuum-{:g}".format(dtc))
+            + ("_surface-skin" + skin_suffix if melt and surface_model == "skin" else "")
             + ("_seed-{}".format(seed) if seed != DEFAULT_SEED else ""))
 
 
+def skin_name_suffix(thickness_mm, cell_um, substep_ms, preset):
+    """The run-name suffix of a skin run's non-default settings (sub-plan 18): `_skin-dev` for the dev preset,
+    `_skin<mm>mm` for another thickness, `_cell<um>um` for an explicit cell other than the preset's, `_sub<ms>ms` for
+    another sub-step; empty for the defaults."""
+    out = "_skin-dev" if preset == "dev" else ""
+    if thickness_mm != 0.4:
+        out += "_skin{:g}mm".format(thickness_mm)
+    if cell_um is not None and abs(cell_um * 1e-6 - skin.SKIN_PRESETS[preset]) > 1e-12:
+        out += "_cell{:g}um".format(cell_um)
+    if substep_ms != 10.0:
+        out += "_sub{:g}ms".format(substep_ms)
+    return out
+
+
 def make_atmosphere(spec, epoch, wind_name):
     """(atmosphere object, short name, provenance dict) for an --atmosphere value."""
     wind = atmosphere.NoWind() if wind_name == "none" else atmosphere.StaticProfileWind()
@@ -233,6 +250,15 @@
                          "with the film plus everything more than half liquid beneath it (slurry included), and where there "
                          "is no conjugate depth a film on slurry does not spray (default on). off: the regime test reads the "
                          "film plus the fully molten material, as in every run before the amendment of 2026-10-06")
+    me.add_argument("--surface-model", choices=body.SURFACE_MODEL_NAMES, default="elements",
+                    help="skin: a one-dimensional skin of fine cells under every melting patch resolves the melt layer and "
+                         "feeds the film continuously, coupled to the 3D conduction through an implicit interface (Step 3 "
+                         "sub-plan 18). elements: the surface elements melt and feed the film, as in every run before it (default)")
+    me.add_argument("--skin-thickness", type=float, default=0.4, help="skin target thickness [mm] (default %(default)s)")
+    me.add_argument("--skin-cell", type=float, default=None, help="skin cell size [um] (default: the preset's)")
+    me.add_argument("--skin-substep", type=float, default=10.0, help="skin sub-step [ms] (default %(default)s)")
+    me.add_argument("--skin-preset", choices=tuple(skin.SKIN_PRESETS), default="production",
+                    help="skin cell-size preset: production 10 um, dev 20 um (default %(default)s)")
     me.add_argument("--rt-spray", choices=("on", "off"), default="on",
                     help="apply Girin & Kopyt's front-surface Rayleigh-Taylor mode as a release mechanism, on the patches where it "
                          "grows faster than the shear mode (default on). off: evaluate and report it but release nothing through it, "
@@ -296,7 +322,8 @@
         melt_settings = body.MeltSettings(removal=args.removal, runoff=args.runoff == "on", demise_fraction=args.demise_fraction,
                                           particles=args.particles, size_feedback=size_feedback,
                                           deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on",
-                                          rigid_substrate=args.rigid_substrate == "on")
+                                          rigid_substrate=args.rigid_substrate == "on", surface_model=args.surface_model,
+                                          skin=skin_settings(args))
         the_body = body.MeltingBody(the_mesh, mat, solver, mass, flow, spray_model, melt_settings, T0=args.temperature,
                                     emissivity=args.emissivity, T_ambient=args.t_ambient)
     else:
@@ -318,12 +345,21 @@
     if melting:
         info.update({"removal": args.removal, "runoff": args.runoff, "deep_runoff": args.deep_runoff,
                      "molten_cascade": args.molten_cascade, "rigid_substrate": args.rigid_substrate, "rt_spray": args.rt_spray,
+                     "surface_model": args.surface_model,
+                     "skin": {"thickness_mm": args.skin_thickness, "cell_um": skin_settings(args).cell * 1e6,
+                              "substep_ms": args.skin_substep, "preset": args.skin_preset},
                      "rarefied_shear": args.rarefied_shear, "we_critical": args.we_critical,
                      "k_r": args.kr, "k_t": args.kt, "demise_fraction": args.demise_fraction, "particles": args.particles,
                      "size_feedback": size_feedback, "gamma_pm": args.gamma_pm, "kn_body_shock": args.kn_body_shock,
                      "liquid": {"rho": mat.liquid.rho, "mu": mat.liquid.mu, "sigma": mat.liquid.sigma},
                      "T_solidus_K": mat.T_solidus, "T_liquidus_K": mat.T_liquidus, "latent_heat_Jkg": mat.latent_heat})
     return the_body, heating_model, info
+
+
+def skin_settings(args):
+    """The skins' settings from the command line (sub-plan 18)."""
+    cell = args.skin_cell * 1e-6 if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset]
+    return skin.SkinSettings(thickness=args.skin_thickness * 1e-3, cell=cell, substep=args.skin_substep * 1e-3)
 
 
 def cmd_run(args, parser):
@@ -363,6 +399,13 @@
     dt_continuum = continuum_step(args)
     if dt_continuum is not None and dt_continuum > args.dt:
         parser.error("--dt-continuum must be no longer than --dt")
+    for label, value in (("--skin-thickness", args.skin_thickness), ("--skin-substep", args.skin_substep),
+                         ("--skin-cell", 1.0 if args.skin_cell is None else args.skin_cell)):
+        if value <= 0.0:
+            parser.error("{} must be > 0".format(label))
+    cell_um = args.skin_cell if args.skin_cell is not None else skin.SKIN_PRESETS[args.skin_preset] * 1e6
+    if 2.0 * cell_um * 1e-3 > args.skin_thickness * (1.0 + 1e-12):
+        parser.error("--skin-cell must be at most half of --skin-thickness")
 
     initial = tj.InitialState(velocity=args.velocity * 1e3, altitude=args.altitude * 1e3,
                               flight_path=math.radians(args.flight_path_angle), heading=math.radians(args.heading),
@@ -375,7 +418,10 @@
     name = args.name or model_run_name(settings.diameter, initial.velocity, initial.altitude, atm_name, args.bridging, args.wind,
                                        args.heating if args.thermal == "fem" else None, args.removal if args.melt == "on" else None,
                                        deep_runoff=args.deep_runoff == "on", molten_cascade=args.molten_cascade == "on", seed=args.seed,
-                                       rigid_substrate=args.rigid_substrate == "on", dt_continuum=dt_continuum)
+                                       rigid_substrate=args.rigid_substrate == "on", dt_continuum=dt_continuum,
+                                       surface_model=args.surface_model,
+                                       skin_suffix=skin_name_suffix(args.skin_thickness, args.skin_cell, args.skin_substep,
+                                                                    args.skin_preset))
     run_dir = os.path.join(args.outdir, name)
     os.makedirs(args.outdir, exist_ok=True)
     thermal_info, the_body = {}, body.ConstantBody(mass, args.temperature)
```

```diff
--- a/reentry_model/coupled.py
+++ b/reentry_model/coupled.py
@@ -30,6 +30,10 @@
                 "delta_m_mean_um", "thick_branch_fraction", "n_dead_elements", "deep_liquid_kg", "deep_mass_kg",
                 "deep_runoff_mass_kg", "deep_surfaced_mass_kg", "deep_blob_fraction", "cascade_passes", "cascade_mass_kg",
                 "nonrigid_depth_mean_mm", "slurry_thick_fraction", "rigid_thin_fraction", "slurry_held_mass_kg"]
+SKIN_COLUMNS = ["skin_count", "skin_mass_kg", "skin_thickness_min_mm", "skin_thickness_median_mm", "skin_drawn_mass_kg",
+                "skin_handover_mass_kg", "skin_short_events", "interface_mismatch_J", "interface_repeats",
+                "liquid_depth_mean_um", "liquid_depth_max_um", "nonrigid_depth_skin_mean_um", "skin_wall_s"]   # sub-plan 18
+MELT_COLUMNS = MELT_COLUMNS + SKIN_COLUMNS
 PVD_TEMPLATE = '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n{}</Collection>\n</VTKFile>\n'
 
 
@@ -143,6 +147,9 @@
                 "r_median_um": float(np.nanmedian(c["r_median_um"])) if np.isfinite(c["r_median_um"]).any() else None,
                 "n_dead_elements": int(body.mesh.n_elements - body.mesh.n_active), "n_source_rows": len(body.source_rows),
                 "removed_enthalpy_J": body.removed_enthalpy, "frozen_mass_kg": body.frozen_mass,
+                "skins_created": body.skins_created, "skin_drawn_mass_kg": body.skin_drawn_mass,
+                "skin_handover_mass_kg": body.skin_handover_mass, "interface_repeats": body.interface_repeats,
+                "skin_wall_s": body.skin_wall,
                 "melt_energy_balance_residual": body.energy_balance_residual()}
 
     def thermal_results(self, history):
@@ -234,6 +241,19 @@
         # the depth the regime test read under each patch: down to the first rigid point, slurry included (amendment of
         # 2026-10-06); the spray step's own, carried like delta_m, nan where it was not evaluated
         poly.cell_data["nonrigid_depth"] = carry(body.last_nonrigid, np.nan)
+        # the skins (sub-plan 18): top temperature, resolved liquid and non-rigid depths, thickness; nan without a skin
+        sk = getattr(body, "skins", None)
+        if sk is not None:
+            on = body.skin_of_patch >= 0
+            r = body.skin_of_patch[on]
+            d_liq, _ = sk.depth_above(body.material.T_feed)
+            d_nr, _ = sk.depth_above(body.material.T_rigid)
+            for key, v in (("skin_T_top", sk.T[:, 0] if sk.n else None), ("skin_liquid_depth", d_liq),
+                           ("skin_nonrigid_depth", d_nr), ("skin_thickness", sk.thickness())):
+                out = np.full(surface.n_patches, np.nan)
+                if sk.n:
+                    out[on] = v[r]
+                poly.cell_data[key] = out
     poly.save(os.path.join(output_dir, "surface_{}.vtp".format(k)))
 
 
```

```diff
--- a/reentry_model/body.py
+++ b/reentry_model/body.py
@@ -1552,4 +1552,23 @@
                 "film_thickness_max_mm": self.film_thickness_max() * 1e3, "film_thickness_mean_mm": self.film_thickness_mean() * 1e3,
                 "nose_radius_mm": self.nose_radius() * 1e3, "transverse_radius_mm": self.transverse_radius * 1e3,
                 "fitted_nose_radius_mm": self.fitted_nose_radius * 1e3,
-                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active)}
+                "n_dead_elements": float(self.mesh.n_elements - self.mesh.n_active), **self._skin_stats()}
+
+    def _skin_stats(self):
+        """The skins' history columns (sub-plan 18); zero or nan without skins."""
+        sk = self.skins
+        th = sk.thickness() if sk is not None and sk.n else np.zeros(0)
+        on = self.skin_of_patch >= 0
+        liq = self.last_skin_liquid if self.last_skin_liquid is not None and len(self.last_skin_liquid) == len(on) else None
+        nr = self.last_skin_nonrigid if self.last_skin_nonrigid is not None and len(self.last_skin_nonrigid) == len(on) else None
+        some = bool(on.any())
+        return {"skin_count": float(sk.n) if sk is not None else 0.0, "skin_mass_kg": sk.mass() if sk is not None else 0.0,
+                "skin_thickness_min_mm": float(th.min() * 1e3) if th.size else float("nan"),
+                "skin_thickness_median_mm": float(np.median(th) * 1e3) if th.size else float("nan"),
+                "skin_drawn_mass_kg": self.skin_drawn_mass, "skin_handover_mass_kg": self.skin_handover_mass,
+                "skin_short_events": float(self.skin_short_events), "interface_mismatch_J": self.last_mismatch,
+                "interface_repeats": float(self.interface_repeats),
+                "liquid_depth_mean_um": float(liq[on].mean() * 1e6) if liq is not None and some else float("nan"),
+                "liquid_depth_max_um": float(liq[on].max() * 1e6) if liq is not None and some else float("nan"),
+                "nonrigid_depth_skin_mean_um": float(nr[on].mean() * 1e6) if nr is not None and some else float("nan"),
+                "skin_wall_s": self.skin_wall}
```

- [ ] **Step 4: Run the tests and the unit tier; confirm the element model's bits**

Run: `"$PY" -m pytest tests/test_reentry_model_cli.py -q -k skin`, then the unit tier.
Expected: 7 passed in about 70 s (the short skin run takes most of it); the unit tier 286 passed, 1 skipped, 2 failed, 3 errors — 34 more passing than the baseline (17 skin, 2 thermal, 8 body, 7 CLI) and the same five known failures; about 5 minutes.

Then confirm bit-identity with the copy before the change (`work-2026-10-07-dt-continuum/code`): the 50 mm default
melting flight to the ground (it never switches step), and the 100 mm default flight to 55 s (it switches to 0.0125 s at
49.5 s), each run from both copies with the same arguments and the default seed, agree in every history column they
share, every source-table row (`particles.npz`) and every result field but the run time. The comparison script used
while planning is `prototype/work-2026-10-07-skinplan/bitid/compare.py`; record its output.

---

### Task 9: Stage 1 measurements, and the amendment written into the sub-plans

**Files:**
- Create: `prototype/work-<date>-skin/meas/` (frozen copy), `.../runs/`, `.../analysis/`
- Modify (documents only): `00-shared-context.md`, `03-thermal-core.md`, `09-melting-body.md`, `10-coupled-loop.md`,
  `13-cli-wiring.md`, `14-verification-runs.md`, `15-documentation.md`, and this sub-plan's status line.

- [ ] **Step 1: Freeze the measurement copy and record its manifest**

`rsync` the tested `code/` to `meas/`; write `meas_manifest.md5` (`find meas -name '*.py' | sort | xargs md5`); never
edit `meas/` while a run uses it.

- [ ] **Step 2: The box against the fine box (spec §12.2)**

Write `analysis/box_ablation.py`: an 8 × 2 × 2 mm box (`mesh.box_mesh(8e-3, 2e-3, 2e-3, h)`) heated by a uniform
2 MW/m² on its x = 0 face only (the face its analytic tests heat), no heating on the other faces and radiation off
everywhere (`MeltingBody(..., emissivity=0.0)`), so that the test is one-dimensional; Scheil material, starting at 700 K,
for 2 s. Liquid is removed as soon as it reaches the patch's liquid (the script zeroes `m_f` and `m_d` after each step
and adds their mass and enthalpy to `removed_mass` and `removed_enthalpy`, so the books still close). Compare (i) the
coarse box, `h = 2 mm`, with skins at the dev and production presets and macro steps of 0.05, 0.1, 0.25 and 0.5 s, with
(ii) the fine box, `h = 0.1 mm`, element model, at 0.005 s. Report the recession (mass removed per unit area) and the
top temperature against time, the mismatch and repeats per step, and the energy balances. The skin must recover the
fine box's recession within 5 % at every macro step, with no oscillation in the top temperature.

- [ ] **Step 3: The flights**

All on the 100 mm physics flight to 120 s (Scheil default, seed 12345, `--surface-model skin`). The macro step after
the continuum switch is set explicitly on every run: `--dt-continuum off` keeps 0.5 s throughout, `--dt-continuum 0.25`
and `0.125` shorten it after the switch.
(a) macro step after the switch 0.5 s (`off`), 0.25 s and 0.125 s, production preset, sub-step 10 ms;
(b) sub-step 20, 10 and 5 ms, and cell 20, 10 and 5 µm, with `--dt-continuum off`;
(c) skin thickness 0.4 against 0.8 mm, with `--dt-continuum off`;
(d) `--h-surface` 2.0 against 1.4 mm, with `--dt-continuum off`;
(e) seed 1 beside the default, with `--dt-continuum off`, for the scatter;
(f) skins at the default 0.0125 s step after the switch, against the element model's default runs of stage 0 (the
converged reference of facts 92–94), if (a)'s measured cost makes it affordable — it is the comparison spec §12.3 calls
"against the converged reference".
Plus the 50 mm whole flight with skins (it never switches). Read every comparison against the scatter at its own step. Report the
droplet population and branch split (expected still step-dependent in stage 1: the film and spray still run once per
macro step, which is stage 2's to change), the masses, the interface mismatch and repeat rate, the skin thickness
statistics and short events, the books, and the cost per macro step and peak memory at each preset. Compare with stage
0's verdict and facts 92–94. The cost measured while planning, for scale: 2000 skins take 1.9–2.7 s for a trial pass and
1.6–2.5 s for a real pass at 10 µm, 1.0–1.7 s and 0.9–1.5 s at 20 µm (one macro step of 0.5 s, single runs).

- [ ] **Step 4: Bit-identity of the element model**

Re-run the copy's earlier 50 mm and 100 mm (to 120 s) default flights with `--surface-model elements` from `meas/` and
show they reproduce the pre-change runs in every shared history column, every source row and every result field but the
run time (Task 8's check, now from the frozen copy).

- [ ] **Step 5: Write the amendment**

In each affected sub-plan, a section "Amendment of <date> — the melt-layer skin, stage 1" with the tested changes as
diff blocks against the copy before the change: 03 the backends' interface (Task 6); 09 `skin.py` in full as a new file
and the body's changes (Tasks 2–5, 7, 8's body part); 10 the columns, results and frame fields; 13 the flags and names;
14 the runs and their results; 15 README, assumptions and the spec's measured values. Add the measurements as new facts
in `00-shared-context.md` (next free numbers) and verify that every committed amendment's blocks plus these reproduce
`meas/` byte for byte from `prototype/proto3/`. Update this sub-plan's status line to "stage 1 built and measured;
stages 2–3 in sub-plan 19". Report to Asha in her reporting style, with the decisions stage 1 raises: the presets, the
interface tolerance (and whether the repeat rate makes it too tight), and whether to proceed to stage 2.

---

## Self-review

**Spec coverage.** §5.1–5.3 (the skin, creation, ownership): Tasks 2 and 7. §5.4 (the sub-step, feed, draw, rezoning,
bounds): Tasks 3, 4 and 7. §5.5 (resolved depths): Tasks 2 (`depth_above`) and 7 (`_compose_depth`). §6 (two clocks,
the interface, the repeat, the books): Tasks 5–7. §7 (the film and spray on the sub-steps): **stage 2, sub-plan 19** —
in stage 1 the film and spray run once per macro step, fed by the skins (spec §12.1, stage 1). §8 (the slurry cascade):
**stage 3, sub-plan 19**. §9 (deaths, hand-over, end of flight): Tasks 2 and 7; remeshing waits for sub-plan 17. §10
(retired and kept): Task 7 — feed off for skin owners, cascade off, remainder to the skin, film on the skin's top,
freeze-back into the skin, resolved depths in the branch test; interior feed, deep runoff, the rigid-substrate rule,
`PHI_DEATH` and the once-per-step surface flow kept. §11 (settings, outputs, code): Task 8, with the CLI default left at
`elements` until stage 2's acceptance. §12.1 stage 0: Task 0; stage 1: Tasks 1–9. §12.2 tests: Tasks 2–7 and Task 9
Step 2. §12.3 measurements: Task 9 Step 3 (acceptance proper is stage 2's).

**Corrections made while planning and checking the code, all recorded in the spec (2026-10-07):** the real pass holds
the base at the 3D model's end-of-step facet temperature, not a ramp (§6.2); the conjugate-depth reach is measured from
the skin's top in metal (§5.4); the skin iterates in the temperature with a line search rather than the 3D solver's
enthalpy-consistent update (§5.4); the base law's slope is clipped at zero, and the interface tolerance is measured
against each patch's energy throughput rather than its base energy alone (§6.2). The surface flow is evaluated before
the passes, at the start-of-step wall temperatures, as the spec's step diagram already has it; the cached liquid
enthalpy's refreshes (Task 7 (d)) are an implementation detail with no spec counterpart.

**Placeholder scan.** None: every code step carries its code or its diff, generated from the tested files.

**Type consistency.** `feed_and_draw` returns six arrays in the order `advance` unpacks; `advance` returns `(a, b)` in
the trial pass and a `SkinBook` otherwise, as `_advance_with_skins` uses them; `hand_over` returns `(rows, mass)`, of
which `_hand_over_skins` books the mass; `_skin_available` serves creation and the step alike; `SKIN_PRESETS` keys are
the CLI's `--skin-preset` choices.

**Review Focus.** Each of the five lines has its test in the owning task (Tasks 2, 4 and 7).
