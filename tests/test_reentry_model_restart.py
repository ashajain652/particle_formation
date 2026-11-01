"""restart.py and the coupled run's save and resume (sub-plan 22 Task 5a; Asha 2026-11-01): a flight saves its whole
state once, at the end of the first macro step at or after a given flight time, without changing anything it does, and
a saved flight resumes from it, in an object graph of its own, to the last bit."""
import csv
import json
import os

import numpy as np
import pytest

from reentry_model import cli, restart, spray
from reentry_model import trajectory as tj

# the spraying skin flight of the command line's test (sub-plan 22 Task 3): 65 km, a body at 893 K under physics heating,
# the dev skins on the 4 mm mesh; a 25 ms first step, then 12.5 ms ones, the tops spraying from the step ending 0.0625 s
FLIGHT = ["run", "--diameter", "100", "--velocity", "7.24", "--altitude", "65", "--flight-path-angle", "-0.959331",
          "--atmosphere", "us76", "--thermal", "fem", "--heating", "physics", "--h-surface", "4", "--h-core", "20",
          "--melt", "on", "--prism-layers", "0", "--temperature", "893", "--dt", "0.025", "--skin-preset", "dev",
          "--fine-size", "2", "--skin-flow-cycles", "2", "--quiet"]


def rows_of(path):
    """A history CSV's rows as text, the wall-time columns left out (they time the run, not the flight)."""
    with open(path) as fh:
        return [{k: v for k, v in r.items() if not k.startswith("wall_")} for r in csv.DictReader(fh)]


def test_a_coarse_grid_solver_is_saved_with_its_factor(tmp_path):
    """pyamg's coarse-grid solver is a class made inside a function, which pickle cannot name; it is saved as its
    solver's name and the factor it has cached, and comes back as the same solver holding the same factor, so a
    hierarchy built before a save solves after it exactly as before."""
    pyamg = pytest.importorskip("pyamg")
    from pyamg.gallery import poisson
    A = poisson((12, 12), format="csr")
    np.random.seed(3)
    ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
    b = np.random.rand(A.shape[0])
    x = ml.solve(b, tol=1e-12)                     # the coarse solver caches its pseudo-inverse on the way
    assert hasattr(ml.coarse_solver, "P")
    path = str(tmp_path / "s.pkl")
    restart.save(path, {"ml": ml})
    assert not os.path.exists(path + restart.PART)
    got = restart.load(path)["ml"]
    assert got is not ml and repr(got.coarse_solver) == repr(ml.coarse_solver)
    np.testing.assert_array_equal(got.coarse_solver.P, ml.coarse_solver.P)
    np.testing.assert_array_equal(got.solve(b, tol=1e-12), x)


def test_a_coarse_grid_solver_with_arguments_is_refused(tmp_path):
    """A coarse solver given arguments keeps them where a save cannot reach them, so saving one is refused rather than
    bringing back a different solver."""
    pyamg = pytest.importorskip("pyamg")
    from pyamg.multilevel import coarse_grid_solver
    with pytest.raises(TypeError, match="coarse"):
        restart.save(str(tmp_path / "s.pkl"), coarse_grid_solver(("pinv", {"rtol": 1e-6})))
    assert not os.path.exists(str(tmp_path / "s.pkl"))


def test_a_flight_saved_on_its_way_is_unchanged_by_the_save_and_resumes_to_the_last_bit(tmp_path):
    """The flight to 0.2125 s twice, once saving its state at 0.1875 s: the same history and droplet records, so a save
    changes nothing in the flight that makes it. The saved state, loaded as a new object graph with the random state it
    carries and a new folder for its droplet records, resumes from the end of its 14th step: the history it makes
    starts with the saved step's row and repeats the flight's from there in every column, and its droplet records are
    the flight's after the save, to the last bit. The heat solver solves twice a step and builds a new multigrid
    hierarchy every 30 solves, drawing from numpy's random state, so the saved one is built before the save and used
    after it, and the 16th step builds a new one: both the saved hierarchy and the saved random state are tested (with
    the generator re-seeded instead, the history differs)."""
    pytest.importorskip("cantera")
    plain, saved = str(tmp_path / "plain"), str(tmp_path / "saved")
    assert cli.main(FLIGHT + ["--t-max", "0.2125", "--outdir", plain, "--name", "f"]) == 0
    assert cli.main(FLIGHT + ["--t-max", "0.2125", "--outdir", saved, "--name", "f", "--save-state-at", "0.1875"]) == 0
    flight = rows_of(os.path.join(plain, "f.csv"))
    assert len(flight) == 17 and rows_of(os.path.join(saved, "f.csv")) == flight
    a, b = np.load(os.path.join(plain, "f", "particles.npz")), np.load(os.path.join(saved, "f", "particles.npz"))
    assert a.files == b.files and len(a["dm_kg"]) > 0 and all(np.array_equal(a[k], b[k], equal_nan=True) for k in a.files)
    doc = json.load(open(os.path.join(saved, "f.json")))
    assert doc["settings"]["save_state_at_s"] == 0.1875 and doc["results"]["state_saved_at_s"] == pytest.approx(0.1875)
    assert doc["results"]["state_save_error"] is None and doc["files"]["state"] == os.path.join(saved, "f.state.pkl")
    assert "save_state_at_s" not in json.load(open(os.path.join(plain, "f.json")))["settings"]

    state = restart.load(doc["files"]["state"])
    assert state["format"] == restart.FORMAT and state["step"] == 14 and state["time_s"] == pytest.approx(0.1875)
    assert state["meta"]["name"] == "f" and state["meta"]["args"]["save_state_at"] == 0.1875
    run = state["run"]
    assert run.body.solver._solves == 28 and run.body.solver.amg_rebuild_every == 30
    np.random.set_state(state["rng"])
    run.body.source_rows = spray.SourceSpool(str(tmp_path / "resumed_rows"))
    run.settings.save_at = None
    history = run.run(resume=state)
    assert run.body.solver._solves == 32                          # the 31st solve built a new hierarchy
    tj.write_history_csv(history, str(tmp_path / "resumed.csv"))
    assert rows_of(str(tmp_path / "resumed.csv")) == flight[14:]
    after = a["time_s"] > state["time_s"]          # the saved step's records carry its end: 0.1875 and a round-off
    got = run.body.source_rows.columns()
    assert after.any() and all(np.array_equal(got[k], a[k][after], equal_nan=True) for k in a.files)


def test_a_save_that_fails_leaves_the_flight_flying(tmp_path, monkeypatch):
    """A save that fails (a full disk, say) is recorded and the flight goes on: the state is worth less than the flight."""
    pytest.importorskip("cantera")

    def broken(path, state):
        raise OSError("no space left on device")

    monkeypatch.setattr(restart, "save", broken)
    assert cli.main(FLIGHT + ["--t-max", "0.0625", "--outdir", str(tmp_path), "--name", "f", "--save-state-at", "0.05"]) == 0
    doc = json.load(open(str(tmp_path / "f.json")))
    assert doc["results"]["end_reason"] == "t_max" and doc["results"]["state_saved_at_s"] is None
    assert "no space left on device" in doc["results"]["state_save_error"] and "state" not in doc["files"]
    assert len(rows_of(str(tmp_path / "f.csv"))) == 5


@pytest.mark.parametrize("argv", [
    FLIGHT + ["--t-max", "0.0875", "--save-state-at", "0"],
    FLIGHT + ["--t-max", "0.0875", "--save-state-at", "-0.05"],
    FLIGHT + ["--t-max", "0.0875", "--save-state-at", "0.0875"],     # nothing would be left to resume
    FLIGHT + ["--t-max", "0.0875", "--save-state-at", "nan"],
    ["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.5", "--t-max", "10", "--save-state-at", "5"],
])
def test_bad_save_times_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli.main(argv + ["--outdir", str(tmp_path)])
    assert exc.value.code == 2
