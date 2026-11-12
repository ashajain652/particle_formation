# Geometry stop and end state Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A skin flight ends itself when its geometry stops describing the metal, and every `--thermal fem` flight of
`run` and `resume` saves a state when it ends, the last completed step's when it fails or is stopped by hand.

**Architecture:** Two measures on `SkinBody` (`geometry_stats`, `geometry_breach`) feed a second body-level end
condition in `CoupledRun.run`, beside demise. A new `CoupledSettings.end_path` makes the run keep an in-memory pickle of
each completed step and write it, or the final state, to `<name>.end.state.pkl` at every end; the command line always
sets it and turns `SIGTERM` into a catchable stop. With `end_path` unset (the library default) the run behaves exactly as
today.

**Tech Stack:** Python 3.12 (`drama_env`), numpy, pytest with pytest-xdist; `restart`'s pickler (protocol 5).

**Spec:** `docs/superpowers/specs/2026-11-12-geometry-stop-and-end-state-design.md` (approved by Asha 2026-11-12).

## Global Constraints

- Run everything from `Reentry Model/` with `PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python`.
- `GEOMETRY_VOLUME_RATIO = 1.25` (Asha's); the centre test fires when `core_front_mm < GEOMETRY_FRONT_LIMIT_MM = 0.0`
  (the sphere's original centre is the origin: `mesh.py` meshes `addSphere(0, 0, 0, radius)`).
- Exit codes: 0 for every normal end, the geometry stop included; 1 for a failure (and escape, as today); 2 for bad
  arguments; 143 for a stop by `SIGTERM` and 130 for `SIGINT` (128 plus the signal number).
- End reasons: `geometry`, `failed`, `stopped` join `ground`, `escape`, `t_max`, `demise`.
- The end state is `<outdir>/<name>.end.state.pkl` (`END_STATE_SUFFIX = ".end.state.pkl"` in `cli.py`); it never
  replaces `--save-state-at`'s `<name>.state.pkl`. `restart.FORMAT` stays 1; a state gains the key `"end"` (the end
  reason it was saved at, or None for a mid-flight save; a missing key reads as None).
- `--thermal none` flights (`sim.run`, no `CoupledRun`) have no state format and are untouched.
- A save that fails is recorded and printed and never changes a flight's end or exit code.
- `--geometry-stop on|off`: skin flights only (exit 2 otherwise), default on for `run` and `resume`, an old state included;
  off is named `_geomstop-off`.
- Commit messages: a scope prefix and a sentence on the behaviour that changed, ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Before the branch is offered for merging the unit tier must pass:
  `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4`
  (main f1ab736 gave 856 passed, 1 skipped, about 74 minutes under load).

## Review Focus

1. A failure after the step's droplet records went to the spool (in `row()` or `_frame`): the end state's records must
   stop at its own step, so a resume does not duplicate them. Test in Task 4.
2. Control-C pressed again while the end state is being written: no partial `.end.state.pkl` is left, the error is
   recorded, the exit code is still 130. Test in Task 4.
3. A failure in the first step of a resume of a state saved before this change: the end state is the resumed step
   itself. Test in Task 4.
4. The geometry test breached on the same step that another end (`t_max`, demise, ground) is reached: that other end
   wins and `geometry_stop` is None. Test in Task 2.
5. A top surface whose triangles are wound inward, or folded: the enclosed volume is the absolute value of the signed
   sum, so the ratio does not change sign with the winding. Test in Task 1.

---

### Task 1: The two geometry measures on the skins

**Files:**
- Modify: `reentry_model/skyscrapers/skinbody.py` (constants beside `LEEWARD_COLUMNS` at about line 210; methods near
  `nose_stats` at about line 1311; `skin_stats` at about line 1346)
- Test: `tests/test_reentry_model_skinbody.py`

**Interfaces:**
- Produces: `GEOMETRY_COLUMNS = ("core_front_mm", "tops_volume_ratio")`; `GEOMETRY_VOLUME_RATIO = 1.25`;
  `GEOMETRY_FRONT_LIMIT_MM = 0.0`; `SkinBody.geometry_stats() -> dict` (the two columns, floats);
  `SkinBody.geometry_breach() -> tuple[str, float] | None`, `("volume", ratio)` when the ratio exceeds
  `GEOMETRY_VOLUME_RATIO`, else `("centre", core_front_mm)` when the front is below `GEOMETRY_FRONT_LIMIT_MM`, else None
  (both read the module constants at call time, so tests can monkeypatch them); `skin_stats()` includes
  `geometry_stats()` right after `nose_stats()`.

- [ ] **Step 1: Write the failing tests** (use the file's `make(coarse_sphere_mesh)` factory)

```python
def test_an_undeformed_body_has_its_front_at_the_radius_and_tops_enclosing_its_metal(coarse_sphere_mesh):
    b = make(coarse_sphere_mesh)
    g = b.geometry_stats()
    assert g["core_front_mm"] == pytest.approx(coarse_sphere_mesh.params["radius_m"] * 1e3, rel=5e-3)
    assert abs(g["tops_volume_ratio"] - 1.0) < 0.01
    assert b.geometry_breach() is None
    assert set(skinbody.GEOMETRY_COLUMNS) <= set(b.skin_stats())

def test_the_breach_names_the_volume_test_first_then_the_centre(coarse_sphere_mesh, monkeypatch):
    b = make(coarse_sphere_mesh)
    monkeypatch.setattr(skinbody, "GEOMETRY_FRONT_LIMIT_MM", 1e3)
    assert b.geometry_breach()[0] == "centre"
    monkeypatch.setattr(skinbody, "GEOMETRY_VOLUME_RATIO", 0.5)
    test, value = b.geometry_breach()
    assert test == "volume" and value == pytest.approx(b.geometry_stats()["tops_volume_ratio"])

def test_the_enclosed_volume_does_not_depend_on_the_winding_of_the_top_triangles(coarse_sphere_mesh):
    b = make(coarse_sphere_mesh)
    before = b.geometry_stats()["tops_volume_ratio"]
    b.surface.faces = b.surface.faces[:, ::-1].copy()
    assert b.geometry_stats()["tops_volume_ratio"] == pytest.approx(before, rel=1e-12)
```

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_skinbody.py -q -k "front_at_the_radius or breach_names or winding"`
Expected: FAIL with `AttributeError: 'SkinBody' object has no attribute 'geometry_stats'`

- [ ] **Step 3: Implement the constants and the two methods in `skinbody.py`**

`core_front_mm`: the largest `self.columns.points @ self.v_hat`, times 1e3 (the base corners are the 3D core's surface).
`tops_volume_ratio`: the absolute value of the signed volume the top surface (`self.surface.points`,
`self.surface.faces`) encloses, by tetrahedra from the points' mean as `spray.enclosed_centroid` builds them, over the
metal's volume, the 3D model's (`self.mesh.volume()`, cached as `mass()` caches it) plus `(self.skins.dz() *
self.skins.area[:, None]).sum()`. Docstrings say what each measures, why (open item 4, Asha 2026-11-12) and the slip-band-on
flight's figure values (1.03 at 179 s, 1.17 at 183 s, 1.42 at 185.5 s).

- [ ] **Step 4: Run the tests to see them pass**

Run: same command. Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add reentry_model/skyscrapers/skinbody.py tests/test_reentry_model_skinbody.py
git commit -m "reentry_model: the skins report every step how far the core's front lies ahead of the sphere's original centre and how many times the metal's volume the column tops enclose (open item 4)"
```

### Task 2: The geometry stop in the run and on the command line

**Files:**
- Modify: `reentry_model/skyscrapers/coupled.py` (the loop at about line 185, `skin_results`), `reentry_model/skyscrapers/skinbody.py`
  (constructor keyword and class default), `reentry_model/cli.py` (`skin_name_suffix`, the run and resume parsers beside
  `--leeward-test`, the skin-only checks at about line 806, `build_thermal`, `cmd_resume`'s changed options and saved-row
  columns at about line 1120, `finish_run`'s print)
- Modify: `README.md` (root) and `Reentry Model/CLAUDE.md`, where they list end reasons and skin options
- Test: `tests/test_reentry_model_restart.py` (it has the small spraying flight `FLIGHT`), `tests/test_reentry_model_cli.py`

**Interfaces:**
- Consumes: Task 1's `geometry_breach()`, `GEOMETRY_COLUMNS`.
- Produces: `SkinBody(..., geometry_stop=True)` and class attribute `geometry_stop = True` (old states fly with it on);
  in `CoupledRun.run`, after the demise check and only while `sim.end_reason is None`:
  `sim.end_reason = "geometry"` and `self.geometry_stopped = (test, t, value)`; `history.results["geometry_stop"]` is
  `{"test": test, "time_s": t, "value": value}` or None for skin runs; `skin_name_suffix(..., geometry_stop=True)`
  adds `_geomstop-off` when False; `thermal_info["geometry_stop"] = "on" | "off"`; the printed line
  `geometry stop at {t:.6g} s: the column tops enclose {v:.3g} times the metal's volume (limit 1.25)` or
  `geometry stop at {t:.6g} s: the core's front is {-v:.3g} mm behind the sphere's original centre`, printed even
  under `--quiet`.

- [ ] **Step 1: Write the failing tests**

```python
def test_the_volume_test_ends_a_skin_flight_normally(tmp_path, monkeypatch):
    monkeypatch.setattr(skinbody, "GEOMETRY_VOLUME_RATIO", 0.5)
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "g"]) == 0
    res = json.load(open(str(tmp_path / "g.json")))["results"]
    assert res["end_reason"] == "geometry" and res["geometry_stop"]["test"] == "volume"
    assert res["geometry_stop"]["time_s"] == pytest.approx(0.025) and len(rows_of(str(tmp_path / "g.csv"))) == 2

def test_the_centre_test_ends_a_skin_flight_and_says_so(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(skinbody, "GEOMETRY_FRONT_LIMIT_MM", 1e3)
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "c"]) == 0
    assert json.load(open(str(tmp_path / "c.json")))["results"]["geometry_stop"]["test"] == "centre"
    assert "geometry stop at 0.025 s: the core's front is" in capsys.readouterr().out

def test_another_end_on_the_same_step_wins_over_the_geometry_stop(tmp_path, monkeypatch):
    calls = {"n": 0}
    real = skinbody.SkinBody.geometry_breach
    def late(self):
        calls["n"] += 1
        return ("volume", 9.0) if calls["n"] >= 2 else real(self)
    monkeypatch.setattr(skinbody.SkinBody, "geometry_breach", late)
    assert cli.main(FLIGHT + ["--t-max", "0.0375", "--outdir", str(tmp_path), "--name", "t"]) == 0
    res = json.load(open(str(tmp_path / "t.json")))["results"]
    assert res["end_reason"] == "t_max" and res["geometry_stop"] is None

def test_the_stop_can_be_turned_off_and_is_then_named(tmp_path, monkeypatch):
    monkeypatch.setattr(skinbody, "GEOMETRY_VOLUME_RATIO", 0.5)
    assert cli.main(FLIGHT + ["--t-max", "0.05", "--outdir", str(tmp_path), "--name", "o", "--geometry-stop", "off"]) == 0
    assert json.load(open(str(tmp_path / "o.json")))["results"]["end_reason"] == "t_max"
    assert cli.skin_name_suffix(0.8, None, 1.10, 10.0, "production", 15.0, geometry_stop=False) == "_geomstop-off"

def test_a_resume_can_turn_the_stop_off_and_old_states_fly_with_it_on(flown, tmp_path, monkeypatch):
    out, flight, records, state_path = flown
    monkeypatch.setattr(skinbody, "GEOMETRY_VOLUME_RATIO", 0.5)
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "on"]) == 0
    assert json.load(open(str(tmp_path / "on.json")))["results"]["end_reason"] == "geometry"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "off", "--geometry-stop", "off"]) == 0
    doc = json.load(open(str(tmp_path / "off.json")))
    assert doc["results"]["end_reason"] == "t_max" and doc["resumed_from"]["changed"] == {"--geometry-stop": "off"}
```

In `tests/test_reentry_model_cli.py`, a parametrized case: `--geometry-stop off` with `--surface-model elements`
exits 2 (follow the file's existing exit-2 cases for skin-only options).

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_restart.py tests/test_reentry_model_cli.py -q -k "geometry or stop_can or same_step"`
Expected: FAIL (`unrecognized arguments: --geometry-stop`; end reason `t_max` where `geometry` is expected)

- [ ] **Step 3: Implement the stop in `coupled.py`, the attribute in `skinbody.py` and the option in `cli.py`**

On resume, `--geometry-stop` given sets `run.body.geometry_stop` and enters `changed`; the saved row gains the
geometry columns where it lacks them, as the nose columns are added (cli.py about line 1120). Update the root
`README.md` and `CLAUDE.md` where they list the end reasons and the skin options.

- [ ] **Step 4: Run the tests to see them pass**

Run: same command. Expected: all selected pass.

- [ ] **Step 5: Commit**

```bash
git add reentry_model tests README.md CLAUDE.md
git commit -m "reentry_model: a skin flight ends normally, with a message, once its column tops enclose more than 1.25 times the metal's volume or its core's front passes the sphere's original centre, unless --geometry-stop off (Asha 2026-11-12)"
```

(`README.md` is at the repository root: `git add ../README.md`.)

### Task 3: A saved state at every normal end, and resuming one

**Files:**
- Modify: `reentry_model/skyscrapers/coupled.py` (`CoupledSettings`, `_save`, the end of `run`),
  `reentry_model/skyscrapers/restart.py`, `reentry_model/cli.py` (`cmd_run`, `cmd_resume`, `finish_run`)
- Test: `tests/test_reentry_model_restart.py`

**Interfaces:**
- Produces: `CoupledSettings.end_path: str = None`; `CoupledRun._save(step, row, path=None, end=None)` (the state dict
  gains `"end": end`; `path` defaults to `settings.save_path`); after the loop, with `end_path` set, the final state is
  saved there with `end=sim.end_reason`, recorded as `self.end_saved_at` / `self.end_save_error` and in
  `history.results` as `end_state_saved_at_s` / `end_state_error`; `restart.dumps(state) -> bytes` and
  `restart.write_bytes(path, data) -> None` (the same pickler and the same part-file-then-rename as `save`, which
  becomes `write_bytes(path, dumps(state))`); `cli.END_STATE_SUFFIX = ".end.state.pkl"`; `cmd_run` and `cmd_resume`
  always set `end_path` and `save_meta`; `doc["files"]["end_state"]` when it was saved.
- Resume of an end state (`state.get("end")`): `t_max` with a later `--t-max`, and `geometry` with
  `--geometry-stop off`, clear `run.sim.end_reason` and fly on; `t_max` without a later `--t-max`, `geometry` with the
  stop still on, and `ground`, `escape`, `demise` are refused with exit 2 and a message naming the end.

- [ ] **Step 1: Write the failing tests**

```python
def test_every_normal_end_saves_its_state_beside_the_mid_flight_one(tmp_path):
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "e", "--save-state-at", "0.05"]) == 0
    doc = json.load(open(str(tmp_path / "e.json")))
    assert doc["files"]["end_state"] == os.path.abspath(str(tmp_path / "e.end.state.pkl"))
    assert restart.load(str(tmp_path / "e.end.state.pkl"))["end"] == "t_max"
    assert restart.load(str(tmp_path / "e.state.pkl"))["end"] is None
    assert doc["results"]["end_state_saved_at_s"] == pytest.approx(0.0625)

def test_an_end_state_at_t_max_resumes_to_a_later_t_max_as_the_unbroken_flight(tmp_path):
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "a"]) == 0
    assert cli.main(FLIGHT + ["--t-max", "0.0875", "--outdir", str(tmp_path), "--name", "u"]) == 0
    assert cli.main(["resume", str(tmp_path / "a.end.state.pkl"), "--t-max", "0.0875",
                     "--outdir", str(tmp_path), "--name", "r"]) == 0
    assert rows_of(str(tmp_path / "r.csv")) == rows_of(str(tmp_path / "u.csv"))[4:]

@pytest.mark.parametrize("extra", [[], ["--t-max", "0.0625"]])
def test_an_end_state_is_not_resumed_past_an_end_that_still_holds(tmp_path, extra):
    assert cli.main(FLIGHT + ["--t-max", "0.05", "--outdir", str(tmp_path), "--name", "a"]) == 0
    path = str(tmp_path / "a.end.state.pkl")
    if extra:
        state = restart.load(path); state["end"] = "demise"; restart.save(path, state)
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", path, "--outdir", str(tmp_path), "--name", "r"] + extra)
    assert exc.value.code == 2
```

Also: the element model's melting flight (`MELTING_ELEMENTS`, `--t-max 0.05`) writes `<name>.end.state.pkl`; and the
existing test whose monkeypatched `restart.save` raises `OSError` still exits 0, now with `end_state` absent from
`doc["files"]` and `end_state_error` set as well.

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_restart.py -q -k "end_state or normal_end or still_holds"`
Expected: FAIL (`KeyError: 'end_state'`)

- [ ] **Step 3: Implement `end_path`, the end save, `dumps`/`write_bytes` and the resume rules**

- [ ] **Step 4: Run the whole restart module**

Run: `"$PY" -m pytest tests/test_reentry_model_restart.py -q -n 4`
Expected: all pass (the existing resume tests now also leave end states in their folders; none of them should change)

- [ ] **Step 5: Commit**

```bash
git add reentry_model tests
git commit -m "reentry_model: every flight of run and resume saves its whole state at its end as <name>.end.state.pkl, and an end state resumes past an end that no longer holds (a later --t-max, --geometry-stop off) (Asha 2026-11-12)"
```

### Task 4: The last completed step on a failure or a stop by hand

**Files:**
- Modify: `reentry_model/skyscrapers/coupled.py` (`run`: the copy and the catch; the post-loop part becomes
  `_history`), `reentry_model/cli.py` (signal handlers around `run.run`, exit codes, `finish_run`),
  `analysis/flights/README.md` (what a flight leaves when it ends), root `README.md` and `CLAUDE.md` (exit codes 130 and
  143, the end reasons `failed` and `stopped`, failures now writing their results files)
- Test: `tests/test_reentry_model_restart.py`

**Interfaces:**
- Consumes: Task 3's `end_path`, `restart.dumps`, `restart.write_bytes`, `_save`'s state dict.
- Produces: `class FlightStopped(Exception)` in `coupled.py` with attribute `signum: int`;
  `CoupledRun._history(rows, states, started) -> tj.History` (today's post-loop code, unchanged in what it computes);
  with `end_path` set, `run()` keeps `self._last = (bytes, n_rows)` from before the loop and after every completed step
  (after the demise, geometry and `save_at` checks; the copy's `"end"` is None), replacing it only once the new pickle is
  complete; on `Exception`, `KeyboardInterrupt` or `FlightStopped` it loads `self._last` back (the consistent run of the
  last completed step), sets the loaded state's `"end"` to `failed` or `stopped` and saves it to `end_path` with
  `restart.save` (any exception of that save, `KeyboardInterrupt` and `FlightStopped` included, goes to
  `end_save_error` and changes nothing else), sets the restored run's `sim.end_reason` likewise, and returns that run's
  `_history(rows[:n], states[:n], started)` with `history.run` (the restored run, whose body `finish_run` must use) and
  `history.results["failure"] = {"kind": "failed" | "stopped", "error": "<Type>: <message>", "time_s": <last completed
  step's time>}`. Without `end_path` it raises as today. The CLI installs a `SIGTERM` handler raising
  `FlightStopped(signal.SIGTERM)` around `run.run` and restores the previous handler afterwards; it prints
  `ERROR: <message>` to standard error as `main` does today and returns 1, `128 + signum` (143), or 130 for
  `KeyboardInterrupt`.

- [ ] **Step 1: Write the failing tests**

```python
def fail_on_call(monkeypatch, n, exc):
    real, calls = skinbody.SkinBody.advance, {"n": 0}
    def advance(self, *a, **k):
        calls["n"] += 1
        if calls["n"] == n:
            raise exc
        return real(self, *a, **k)
    monkeypatch.setattr(skinbody.SkinBody, "advance", advance)

def test_a_failure_saves_the_last_completed_step_and_writes_its_results(tmp_path, monkeypatch):
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "u"]) == 0
    fail_on_call(monkeypatch, 3, RuntimeError("a corner move degrades the mesh"))
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "f"]) == 1
    doc = json.load(open(str(tmp_path / "f.json")))
    assert doc["results"]["end_reason"] == "failed" and "degrades the mesh" in doc["results"]["failure"]["error"]
    assert rows_of(str(tmp_path / "f.csv")) == rows_of(str(tmp_path / "u.csv"))[:3]
    state = restart.load(str(tmp_path / "f.end.state.pkl"))
    assert state["end"] == "failed" and state["step"] == 2 and state["time_s"] == pytest.approx(0.0375)
    monkeypatch.undo()
    assert cli.main(["resume", str(tmp_path / "f.end.state.pkl"), "--outdir", str(tmp_path), "--name", "r"]) == 0
    assert rows_of(str(tmp_path / "r.csv")) == rows_of(str(tmp_path / "u.csv"))[2:]

def test_a_sigterm_stops_the_flight_with_the_last_completed_step_and_exit_143(tmp_path, monkeypatch):
    before = signal.getsignal(signal.SIGTERM)
    real, calls = skinbody.SkinBody.advance, {"n": 0}
    def advance(self, *a, **k):
        calls["n"] += 1
        if calls["n"] == 2:
            os.kill(os.getpid(), signal.SIGTERM)
        return real(self, *a, **k)
    monkeypatch.setattr(skinbody.SkinBody, "advance", advance)
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "s"]) == 143
    assert json.load(open(str(tmp_path / "s.json")))["results"]["end_reason"] == "stopped"
    assert restart.load(str(tmp_path / "s.end.state.pkl"))["step"] == 1
    assert signal.getsignal(signal.SIGTERM) == before

def test_control_c_during_the_end_write_leaves_no_partial_state(tmp_path, monkeypatch):
    fail_on_call(monkeypatch, 2, KeyboardInterrupt())
    def interrupted(path, data):
        with open(path + restart.PART, "wb") as fh:
            fh.write(data[:10])
        raise KeyboardInterrupt
    monkeypatch.setattr(restart, "write_bytes", interrupted)
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "k"]) == 130
    assert not os.path.exists(str(tmp_path / "k.end.state.pkl"))
    assert json.load(open(str(tmp_path / "k.json")))["results"]["end_state_error"]

def test_a_failure_in_the_first_step_of_a_resume_saves_the_resumed_step(flown, tmp_path, monkeypatch):
    out, flight, records, state_path = flown
    fail_on_call(monkeypatch, 1, RuntimeError("boom"))
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "r"]) == 1
    assert restart.load(str(tmp_path / "r.end.state.pkl"))["step"] == 14

def test_a_failure_after_the_droplet_records_were_spooled_resumes_without_repeating_them(tmp_path, monkeypatch):
    assert cli.main(FLIGHT + ["--t-max", "0.0875", "--outdir", str(tmp_path), "--name", "u"]) == 0
    real, calls = coupled.CoupledRun._frame, {"n": 0}
    def frame(self, *a, **k):
        calls["n"] += 1
        if calls["n"] == 6:                        # the frame of step 5 (0.075 s), after its records were spooled
            raise RuntimeError("frame writer failed")
        return real(self, *a, **k)
    monkeypatch.setattr(coupled.CoupledRun, "_frame", frame)
    assert cli.main(FLIGHT + ["--t-max", "0.0875", "--outdir", str(tmp_path), "--name", "f"]) == 1
    monkeypatch.undo()
    assert cli.main(["resume", str(tmp_path / "f.end.state.pkl"), "--outdir", str(tmp_path), "--name", "r"]) == 0
    whole = np.load(str(tmp_path / "u" / "particles.npz"))
    got = np.load(str(tmp_path / "r" / "particles.npz"))
    after = whole["time_s"] > restart.load(str(tmp_path / "f.end.state.pkl"))["time_s"] + 1e-9
    assert all(np.array_equal(got[k], whole[k][after], equal_nan=True) for k in whole.files)
```

(The `_frame` call count includes the call before the loop; check the count against the loop when implementing and
adjust the number so the failure falls after step 5's records are spooled. If the restored run's spool count and the
spool file disagree, the fix is to cut the spool back to the restored count when restoring.)

- [ ] **Step 2: Run them to see them fail**

Run: `"$PY" -m pytest tests/test_reentry_model_restart.py -q -k "failure or sigterm or control_c"`
Expected: FAIL (exit codes 1 from `main`'s `RuntimeError` branch with no files; 143 not returned)

- [ ] **Step 3: Implement the copy, the catch, `_history`, `FlightStopped`, the handlers and the exit codes**

- [ ] **Step 4: Run the restart and command-line modules**

Run: `"$PY" -m pytest tests/test_reentry_model_restart.py tests/test_reentry_model_cli.py -q -n 4`
Expected: all pass; any existing test that expected exit 2 from a `ValueError` raised inside a flight now sees 1 (a
model failure): report each such test to the controller rather than changing it.

- [ ] **Step 5: Run the unit tier**

Run: `PYTHONDONTWRITEBYTECODE=1 "$PY" -m pytest tests -m "not drama and not reference" -q -p no:cacheprovider -n 4`
Expected: 856 + the new tests passed, 1 skipped.

- [ ] **Step 6: Commit**

```bash
git add reentry_model tests analysis/flights/README.md CLAUDE.md ../README.md
git commit -m "reentry_model: a flight that fails or is stopped by SIGTERM or Control-C saves its last completed step as <name>.end.state.pkl, writes its results with end reason failed or stopped, and exits 1, 143 or 130 (Asha 2026-11-12)"
```
