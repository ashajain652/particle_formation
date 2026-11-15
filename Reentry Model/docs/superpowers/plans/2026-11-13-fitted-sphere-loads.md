# The skins' windward loads from one fitted sphere: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The skins' windward loads come from one sphere fitted to the exposed top corners: those ahead of the
mass-centroid plane whose upstream ray does not re-enter the tops. Each exposed corner is read at the angle of the ray
from the sphere's centre. The old rule stays available as `--windward-loads normals`.

**Architecture:** A new skyscrapers module, `shelter.py`, casts the upstream rays against the top surface on a grid.
`SkinBody` gains a `windward_loads` setting and a lazily made `SphereFit`. These feed its existing seams:
`load_angles`, `nose_fit`/`nose_radius` and `windward_corners`, plus a new `load_meridian` that points the push on
the melt. The command line adds the option to `run` and `resume`. No consumer's interface changes.

**Tech Stack:** Python 3.12 (drama_env), numpy, pytest with xdist.

**Spec:** `docs/superpowers/specs/2026-11-13-fitted-sphere-loads-design.md`. It was approved by Asha on 2026-11-13.
Read it, with its figures, before Task 1.

## Global Constraints

- **Branch.** Work on `fitted-sphere-loads`, cut from main at 4e857e1, in the worktree
  `.claude/worktrees/geometry-stop`. Commit locally only, and never push.
  - Never touch the `interim-gas-loads` branch or its worktree.
  - Merge into main only after `interim-gas-loads` has landed there (Task 7).
- **Commands.** Run from `Reentry Model/` with `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`. The unit
  tier is `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4`.
- **The option.** `--windward-loads sphere|normals` is a `run` and `resume` option.
  - Its default is `sphere`, on the command line and in the `SkinBody` constructor.
  - The class attribute is `windward_loads = "normals"`, so states pickled before the change fly as they were.
  - `normals` marks the run name `_windward-normals`. The default is unmarked.
- **The candidates.** A corner is a candidate when `(points - centroid) @ v_hat > 0` (`SkinBody._leeward`).
- **The shelter test.**
  - The ray starts at `points + normals * 0.05 * sqrt(areas)`, using the top surface's corner normals and areas.
  - It runs along `+v_hat`, the direction of motion.
  - A corner is sheltered when the ray crosses any triangle of the top surface at a positive distance. Edges and
    corners of triangles count, and the corner's own triangles count too.
  - The grid version must give exactly the plain search's answer.
- **The sphere.**
  - It is `body.fit_sphere` on the exposed corners (candidates that are not sheltered).
  - The sphere is not valid when there are fewer than 4 exposed corners, the radius is not finite and positive,
    the exposed corners are coplanar to rounding (smallest singular value below 1e-9 of the largest), or the fit is
    concave (`(centre - exposed.mean(0)) @ v_hat > 0`).
- **The radius used.**
  - With a valid sphere, it is `min(fitted, bound)` when the bound is finite and positive, and `fitted` otherwise.
  - With no valid sphere, it is the bound, or the mesh radius `mesh.params["radius_m"]` where there is no bound.
  - The bound is today's: `NOSE_CAP_FACTOR` times the windward transverse radius, taken over the corners whose base
    normal faces the flow and that are ahead of the plane.
- **The angles.** An exposed corner reads its ray angle from the fitted centre, or its own `theta` when there is no
    valid sphere. Every other corner reads `max(theta, pi/2 + 1e-6)`.
- **The push direction under `sphere`.** For exposed corners with a valid sphere, it is the fitted sphere's meridian
  away from its stagnation point, `-(v - (r.v) r)` with `r` the unit ray from the centre, projected onto the base's
  tangent plane and normalised. Every other corner keeps `skin.along_meridian`'s direction.
- **The history.** Under `normals` the history and the run JSON are today's. The JSON gains one key, `windward_loads`.
  Under `sphere`, the nose columns carry the new fit, and
  `SPHERE_COLUMNS = ("sheltered_corners", "sphere_centre_mm", "sphere_fit_rms_mm")` follow the two leeward columns.
- **Byte-identity.** `--windward-loads normals` is byte-identical to main at 4e857e1, apart from the JSON key.
- **Speed.** The ray test takes under 0.25 s on the 100 mm saved state (9,417 corners).

## Review Focus

1. **A state pickled before the change flies `normals` exactly.** The class attribute covers it. Tested in Task 4.
2. **A flat or dished exposed face (a pancake).** It falls back to normal angles at the bound, marked in
   `nose_cap_concave`. Tested in Task 2.
3. **Rays that graze a triangle's edge or corner, or that fall on a grid-cell boundary.** A triangle must be listed in
   every cell its projection overlaps, and edges count as hits. The grid must equal the plain search on a ringed
   sphere and on a jittered one. Tested in Task 1.
4. **A flight direction other than +x.** The test bodies fly along -x. The ray test and the meridian must use `v_hat`,
   whatever it is. Tested in Tasks 1 and 3 with -x and with an oblique direction.
5. **A resume that switches the rule.** It must add or drop `SPHERE_COLUMNS` in the saved row, so that the history
   has one header. Tested in Task 4.

---

### Task 1: The upstream-ray shelter test

**Files:**
- Create: `reentry_model/skyscrapers/shelter.py`
- Test: `tests/test_reentry_model_shelter.py`

**Interfaces:**
- Produces:
  - `OFFSET_FRACTION = 0.05`.
  - `upstream_hits(points, faces, normals, areas, v_hat, candidates=None) -> np.ndarray`: bool `(n,)`. It is True where
    the corner's ray (Global Constraints) crosses a triangle of `(points, faces)` at a positive distance along
    `v_hat`. Only the `candidates` mask is tested, and every other corner is False. It is grid-accelerated.
  - `upstream_hits_plain(...)`: the same signature, the plain search over every triangle, kept as the reference.

The approach: project the triangles onto the plane across `v_hat`, using an orthonormal basis `(e1, e2)`. Lay a
uniform grid with cells about one median edge length wide, and list each triangle in every cell its projected
bounding box overlaps. Each ray tests only its own cell's triangles: an inclusive barycentric test, then the hit's
position along `v_hat` must exceed the start's. The investigation's version is in
`Reentry Model/output/shape_feedback/scripts/shelter.py` in the main checkout. That version also skipped the corner's
own triangles, which this design does not.

- [ ] **Step 1: Write the failing tests** (on `coarse_sphere_mesh`, building the top surface as
  `columns.CornerSurface.from_points(points, faces)` from the mesh's surface corners moved out 0.8 mm):

```python
def test_no_corner_of_a_sphere_is_sheltered(coarse_sphere_mesh): ...
    # for v_hat in (-1,0,0), (1,0,0) and an oblique unit vector: candidates = the corners with (p - c) . v > 0
    assert not shelter.upstream_hits(P, F, N, A, v, candidates=ahead).any()

def test_a_raised_ring_shelters_the_band_behind_it_and_exposes_its_crest(coarse_sphere_mesh): ...
    # the tops of corners whose base normal is 45-52 degrees from v_hat raised 4 mm along their base normals
    # (raise further if no ring corner faces away; assert that one does)
    hit = shelter.upstream_hits(P, F, N, A, v, candidates=ahead)
    assert not hit[base_angle < 30].any()                      # the nose is untouched
    assert hit[(base_angle > 52) & (base_angle < 65)].any()     # the band just behind the ring is sheltered
    assert (base_angle[hit] >= 45).all() and (base_angle[hit] <= 75).all()
    assert (ring & (theta > np.pi / 2) & ~hit).any()            # a crest corner faces away and is exposed

def test_the_grid_gives_the_plain_search_s_answer(coarse_sphere_mesh): ...
    # the ringed sphere and a sphere whose tops are jittered by +-0.3 mm (seeded), for v_hat -x and oblique
    assert np.array_equal(shelter.upstream_hits(*args, candidates=ahead), shelter.upstream_hits_plain(*args, candidates=ahead))

def test_only_candidates_are_tested(coarse_sphere_mesh): ...
    assert not shelter.upstream_hits(*ringed_args, candidates=np.zeros(n, bool)).any()
```

- [ ] **Step 2: Run** `"$PY" -m pytest tests/test_reentry_model_shelter.py -q`. Expect FAIL, because there is no
  module `shelter`.
- [ ] **Step 3: Implement** `shelter.py` as above, with a module docstring that cites the design and the level-set
  test `levelset/shadow.py`.
- [ ] **Step 4: Run** the same command. Expect PASS. Also time `upstream_hits` against `upstream_hits_plain` on the
  coarse ringed sphere and put both times in the report.
- [ ] **Step 5: Commit** with the message `reentry_model: the skins' upstream-ray shelter test (a corner is sheltered
  where its ray up the flow re-enters the tops), grid-accelerated with a plain-search reference`.

### Task 2: The sphere rule in `SkinBody`

**Files:**
- Modify: `reentry_model/skyscrapers/skinbody.py`. The parts are: the constants near `NOSE_COLUMNS` (:206) and
  `LEEWARD_COLUMNS` (:219), `NoseFit` (:231), the class attributes (:748-761), `__init__` (:763), `_refresh_geometry`
  (:876-889), `windward_corners` (:916), `load_angles` (:922), `_fit_nose` (:933, with its bound extracted), `nose_fit`
  (:969), and `nose_stats` (:1409).
- Modify: `reentry_model/skyscrapers/coupled.py`, in `write_skin_frame` (:513-534).
- Modify: `tests/test_reentry_model_skinbody.py`. This covers the `with_tops` helper (:2066, which also resets
  `_sphere_raw`) and the existing tests of the old rule.
- Test: `tests/test_reentry_model_skinbody.py` (new tests) and `tests/test_reentry_model_coupled.py` (frames).

**Interfaces:**
- Consumes: `shelter.upstream_hits` (Task 1).
- Produces:
  - `WINDWARD_LOADS = ("sphere", "normals")`.
  - `SPHERE_COLUMNS` (Global Constraints).
  - `@dataclass SphereFit`. Its fields are: `centre` ((3,) m, nan where there is no valid sphere), `fitted` (m, nan),
    `bound` (m), `used` (m), `exposed` and `sheltered` (bool (n,)), `points` (int, the exposed count), `concave` (bool,
    True where there is no valid sphere) and `rms` (m, nan).
  - `SkinBody(..., windward_loads="sphere")`, validated against `WINDWARD_LOADS`. The class attribute is `"normals"`,
    and the class also has `_sphere_raw = None`.
  - `SkinBody.sphere_fit` (property, lazy, reset in `_refresh_geometry` beside `_nose_raw`).
  - Under `sphere`:
    - `windward_corners()` is `exposed & (load_angles() <= pi/2)`;
    - `load_angles()` follows Global Constraints;
    - `nose_fit` is `NoseFit(fitted, bound, points, concave, used)` from the `SphereFit`, so `nose_radius()` and the
      nose columns carry the new fit;
    - `nose_stats()` appends `SPHERE_COLUMNS`. These are `sheltered_corners` as an int, `sphere_centre_mm` as
      `centre @ v_hat * 1e3` (from the mesh's origin, the original centre; nan without a sphere) and
      `sphere_fit_rms_mm`.
  - Under `normals`, every one of these is today's code path, untouched. `leeward_test` and `size_feedback` are read
    only under `normals`.
  - `write_skin_frame` adds the point data `load_angle_deg` (always) and `sheltered` (under `sphere`).

- [ ] **Step 1: Write the failing tests:**

```python
def test_an_untouched_sphere_fits_its_own_sphere_and_reads_its_normals_angles(coarse_sphere_mesh): ...
    b = make(coarse_sphere_mesh)                         # the library default is now sphere
    sf = b.sphere_fit
    assert b.windward_loads == "sphere" and not sf.concave and not sf.sheltered.any()
    assert abs(sf.fitted / (0.05 + 0.0008) - 1) < 0.01 and np.linalg.norm(sf.centre) < 0.01 * 0.05
    assert np.abs(b.load_angles()[sf.exposed] - b.theta[sf.exposed]).max() < TOL   # TOL measured first; record it
def test_moved_tops_move_the_sphere_s_centre_and_scaled_tops_scale_its_radius(coarse_sphere_mesh): ...
    # tops translated 5 mm along v_hat: the centre moves 5 mm; tops scaled 1.2 about the origin: the radius scales 1.2
def test_a_sheltered_corner_is_leeward_and_a_crest_facing_away_reads_its_ray_angle(coarse_sphere_mesh): ...
    # the raised ring of Task 1 through with_tops: sheltered corners read max(theta, pi/2 + 1e-6), and the away-facing
    # crest corner reads its ray angle (< pi/2) and is in windward_corners()
def test_the_radius_used_is_the_fit_capped_at_the_flat_face_bound(coarse_sphere_mesh): ...
    # flatten_front(factor=2.0) across the whole front (rho_max=1.0): fitted > bound and nose_radius() == bound,
    # angles still from the fitted centre
def test_a_flat_exposed_face_falls_back_to_normal_angles_at_the_bound(coarse_sphere_mesh): ...
    # every exposed top projected onto one plane across v_hat: sf.concave, load_angles() == theta on the exposed
    # corners, nose_radius() == bound, skin_stats()["nose_cap_concave"] == 1.0
def test_the_history_has_the_sphere_columns_only_under_sphere(coarse_sphere_mesh): ...
    assert set(SPHERE_COLUMNS) <= set(make(m).skin_stats())
    assert not set(SPHERE_COLUMNS) & set(make(m, windward_loads="normals").skin_stats())
def test_a_state_without_the_attribute_flies_normals(coarse_sphere_mesh): ...
    b = make(m); del b.__dict__["windward_loads"]; assert b.windward_loads == "normals"
def test_the_frame_carries_the_load_angle_and_the_shelter(...):  # coupled: both point arrays present under sphere
```

- [ ] **Step 2: Run** `"$PY" -m pytest tests/test_reentry_model_skinbody.py -q -n 4 -k "sphere or sheltered or crest or flat_exposed or radius_used or frame"`.
  Expect FAIL.
- [ ] **Step 3: Implement** as in the Interfaces block. `_fit_nose` keeps its body for `normals`. Its bound computation
  moves into `_windward_bound()`, which both rules use.
- [ ] **Step 4: Re-pin the existing tests.** The library default is now `sphere`, so run
  `tests/test_reentry_model_skinbody.py tests/test_reentry_model_coupled.py tests/test_reentry_model_melting.py tests/test_reentry_model_restart.py -n 4`.
  Treat each failure as one of two kinds:
  - **A test of the old rule's mechanics** (the 30-degree cap, `leeward_test`, `size_feedback`, for example at
    :2093-2140, :2472-2583, :2616, :2640). Give it `windward_loads="normals"`.
  - **A pinned value that moved with the rule.** Give it `windward_loads="normals"` and say so in its docstring, as the
    14 bridge-law tests say they belong to the bridge law. Do not change the value.

  List every re-pinned test in the report.
- [ ] **Step 5: Run** the four modules again. Expect PASS.
- [ ] **Step 6: Commit** with the message `reentry_model: under --windward-loads sphere (the new default) the skins
  read their windward loads from one sphere fitted to the exposed tops, each corner at the angle of the ray from its
  centre; sheltered corners are leeward`.

### Task 3: The push on the melt follows the fitted sphere's meridian

**Files:**
- Modify: `reentry_model/skyscrapers/skin.py`. The parts are `push` (:1253, its pressure term at :1284), `Flow`
  (:1907-1940) and the call in `FlowBook` (:2042).
- Modify: `reentry_model/skyscrapers/skinbody.py`, the `Flow` built in `_advance_part` (:1102-1104).
- Test: `tests/test_reentry_model_skin.py` and `tests/test_reentry_model_skinbody.py`.

**Interfaces:**
- Consumes: `SkinBody.sphere_fit` (Task 2).
- Produces:
  - `push(..., wet=None, meridian=None)`. Given `meridian` ((n, 3) unit vectors), the pressure term is
    `-meridian * dp_ds[:, None]`. Otherwise it is today's `-along_meridian(normals, v_hat, dp_ds)`, bit for bit.
  - `Flow.meridian: np.ndarray = None`, passed by `FlowBook` to `push`.
  - `SkinBody.load_meridian()`, which returns `(n, 3)` under `sphere` (Global Constraints) and `None` under `normals`.
  - In `_advance_part`, the shear is `mer * f.tau[:, None]` when `mer` is not None. Otherwise it stays today's
    `along_meridian(...)`.

- [ ] **Step 1: Write the failing tests:**
  - `push` with `meridian=None` equals today's to the last bit.
  - With a given `meridian`, the pressure term equals `-meridian * dp_ds[:, None]`.
  - On the untouched coarse sphere, `load_meridian()` is within 2 degrees of `along_meridian` on every exposed corner
    with `theta > 5` degrees.
  - With the tops translated 10 mm across the flight direction, it points away from the moved stagnation point, both
    for v_hat -x and for an oblique direction.
  - Under `normals`, `load_meridian()` is None and the `Flow` built is today's.
- [ ] **Step 2: Run.** Expect FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `tests/test_reentry_model_skin.py tests/test_reentry_model_skinbody.py -n 4`. Expect PASS.
- [ ] **Step 5: Commit** with the message `reentry_model: under --windward-loads sphere the pressure gradient and the gas
  shear push the melt along the fitted sphere's meridian, away from its stagnation point`.

### Task 4: The command line, run and resume

**Files:**
- Modify: `reentry_model/cli.py`. The parts are: the helpers near `skin_leeward_test_of` (:171), `skin_name_suffix`
  (:286, the marks at :350-355), the `resume` options (:403-415), the `run` options (:603-627, after `--geometry-stop`),
  the `SkinBody` built (:680-687), the run JSON (:720-730), the run-name call (:786-795), the checks (:894-897),
  `cmd_resume`'s checks (:1204-1215), how it applies them (:1283-1289), the saved row (:1306-1314), the saved
  arguments (:1331-1336), the printed line (:1340-1356) and `thermal_info` (:1378-1383). Also the usage docstring
  (:8-19).
- Test: `tests/test_reentry_model_cli.py` (`BASE_RUN_OPTION_STRINGS` :1302 gains `'--windward-loads'` after
  `'--geometry-stop'`) and `tests/test_reentry_model_restart.py`.

**Interfaces:**
- Consumes: `skinbody.WINDWARD_LOADS` and `SPHERE_COLUMNS`, and `SkinBody.windward_loads` (Task 2).
- Produces:
  - `skin_windward_loads_of(args) -> str`, which is `args.windward_loads or "sphere"`.
  - `skin_name_suffix(..., windward_loads="sphere")`, which appends `_windward-normals` for `normals`.
  - `run`:
    - `--windward-loads` needs `--surface-model skin` (exit 2).
    - Under a resolved `sphere`, an explicit `--leeward-test`, or `--size-feedback` with the skins, exits 2 with
      "applies only with --windward-loads normals".
    - The JSON's `skin.windward_loads` is always written. Under `sphere`, `skin.leeward_test` and
      `skin.size_feedback` are not written.
  - `resume`:
    - `--windward-loads` needs a saved flight of skins (exit 2).
    - The resolved rule is the argument, or the body's own (`normals` for an old state). The same refusals apply
      under `sphere`.
    - A change sets `run.body.windward_loads` and clears `_sphere_raw` and `_nose_raw`.
    - The saved row pops `SPHERE_COLUMNS` and then takes `nose_stats()`, with `redo` true for this option.
    - It is recorded in the saved arguments, in the printed line ("windward loads {}") and in
      `thermal_info["windward_loads"]`.

- [ ] **Step 1: Write the failing tests:**
  - The option list includes `--windward-loads`.
  - `skin_name_suffix` marks only `normals`.
  - Each refusal exits 2.
  - A skins run's JSON has `windward_loads: "sphere"`.
  - The `flown` state resumes under `sphere` with `SPHERE_COLUMNS` in its CSV header.
  - `resume --windward-loads normals` drops them and records `normals`.
  - A state whose body lost the attribute (`del body.__dict__["windward_loads"]` before saving) resumes `normals`,
    and its history header equals a `normals` run's.
  - Existing tests that pass `--leeward-test` or `--size-feedback` with the skins get `--windward-loads normals`.
- [ ] **Step 2: Run** `tests/test_reentry_model_cli.py tests/test_reentry_model_restart.py -n 4`. Expect FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the same command. Expect PASS.
- [ ] **Step 5: Commit** with the message `reentry_model: run and resume take --windward-loads sphere|normals (sphere by
  default; normals is the earlier rule, named _windward-normals, and states saved before the option resume on it)`.

### Task 5: Verification on the real bodies and against main

**Files:**
- Create: `analysis/fitted_sphere_check.py`. It defaults to the two saved states in
  `.claude/worktrees/subplan20/prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_compare/`:
  `d100__handdown-off.state.pkl` and `d050__handdown-off.state.pkl`. For each state it prints the sheltered and
  exposed counts, the fit, the radius used, the centre and the rms, and times `upstream_hits` against
  `upstream_hits_plain`. It sets `windward_loads = "sphere"` on the loaded body.

- [ ] **Step 1: Run the check.** Run `"$PY" analysis/fitted_sphere_check.py`. Expect:
  - 100 mm body: 568 sheltered of 4,737 candidates, and a fit of 50.89 mm within 0.05 mm.
  - 50 mm body: 233 of 1,255, and 49.03 mm within 0.05 mm.
  - `upstream_hits` under 0.25 s on the 100 mm body.

  If the counts differ, stop and report: the investigation skipped own triangles, but measured no difference.
- [ ] **Step 2: Check byte-identity with main.** Make a detached checkout of 4e857e1
  (`git worktree add --detach "$SCRATCH/main-4e857e1" 4e857e1`).
  - Run `analysis/skyscrapers_port/same_results.py run` for the `skins` case. Use main's code, then this branch with
    `--append "--windward-loads normals"`. The element and non-melting cases run without the option, which they refuse.
  - Run `compare --report`. The only difference must be the JSON key `windward_loads`.
  - Remove the detached worktree afterwards.
- [ ] **Step 3: Fly the default rule.** Run the short melting smoke flight from `CLAUDE.md` with no
  `--windward-loads`. It must exit 0. Report its wall time per macro step against the same flight under `normals`.
- [ ] **Step 4: Commit** `analysis/fitted_sphere_check.py` with a message that states the measured counts and times.

### Task 6: Documents and the level-set design

**Files:**
- Modify: `docs/model_assumptions.md`. Add a new section, "The skins' windward loads: the fitted sphere and the
  shelter test", covering:
  - the rule, and each decision with its date;
  - the measured counts from Task 5;
  - the straight-line shadow, which leaves out reattachment heating behind a step;
  - the back slope beside a heated crest, which is under-resolved;
  - that the element model is unchanged.

  Number it after the existing sections. The interim gas loads add a section too, so the numbers are settled at
  Task 7.
- Modify: the root `README.md` (the option list) and `Reentry Model/CLAUDE.md` (the skins' description, beside
  `--geometry-stop`).
- Modify: `docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`, as in the design's
  level-set section:
  - in section 20.8, the windward rule, the fitted sphere and the ray angle;
  - in section 20.14, `--leeward` gains `sphere`, the default;
  - in sections 20.10 and 20.11 (item 5), the parity flights use `--leeward sphere` against `--windward-loads sphere`;
  - a dated line in its decisions log.
- Modify: the design's status line, to say it is built.

- [ ] **Step 1: Edit** the documents. Every number quoted must come from Tasks 1-5's reports.
- [ ] **Step 2: Commit** with the message `docs: the skins' windward loads from the fitted sphere, with the shelter
  test, in the model assumptions, README and CLAUDE.md; the level-set design takes the same rule`.

### Task 7: Merge after the interim gas loads have landed

- [ ] **Step 1: Wait.** Watch `main` and `origin/main` until `interim-gas-loads` has been merged and the session
  building it is idle. Then tell Asha before going on.
- [ ] **Step 2: Merge** `main` into `fitted-sphere-loads` in this worktree.
  - In `SkinBody.load_angles`, the corners behind the plane and the sheltered corners read `np.pi`, as their Task 7
    gives the separated wake. Exposed corners keep their ray angle, even past 90 degrees.
  - Under `normals`, follow their code.
  - Settle the `model_assumptions` section numbers.
- [ ] **Step 3: Run** the unit tier on the merged tree. Expect PASS. Report the count against both branches' counts.
- [ ] **Step 4: Fast-forward main.** Check that main has not moved, fast-forward local main, and report to Asha.
  Pushing is her decision.
