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


@pytest.fixture(scope="module")
def flown(tmp_path_factory):
    """The test flight to 0.2125 s once, saving its state at 0.1875 s (the end of its 14th step), for the resume command's
    tests: (the flight's directory, its history rows as text, its droplet records, the saved state's path)."""
    pytest.importorskip("cantera")
    out = str(tmp_path_factory.mktemp("flown"))
    assert cli.main(FLIGHT + ["--t-max", "0.2125", "--outdir", out, "--name", "f", "--save-state-at", "0.1875"]) == 0
    rows = rows_of(os.path.join(out, "f.csv"))
    assert len(rows) == 17
    return out, rows, dict(np.load(os.path.join(out, "f", "particles.npz"))), os.path.join(out, "f.state.pkl")


def test_the_resume_command_repeats_the_flight_to_the_last_bit_and_writes_its_files(flown, tmp_path):
    """`resume` on the state saved at 0.1875 s, with nothing changed: it writes the files `run` writes under its own name
    -- the history, the JSON, the droplet files -- the history starting with the saved step's row and repeating the
    unbroken flight's from there in every column, the droplet records those of the flight after the save, to the last bit,
    and the JSON naming the state it came from with nothing changed (sub-plan 22 Task 5b, Asha 2026-11-06)."""
    out, flight, records, state_path = flown
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "r"]) == 0
    assert rows_of(str(tmp_path / "r.csv")) == flight[14:]
    after = records["time_s"] > 0.1875 + 1e-9        # the saved step's records carry its end: 0.1875 and a round-off
    got = np.load(str(tmp_path / "r" / "particles.npz"))
    assert after.any() and got.files == list(records) and all(np.array_equal(got[k], records[k][after], equal_nan=True) for k in records)
    assert os.path.exists(str(tmp_path / "r" / "particles_summary.csv")) and not os.path.exists(str(tmp_path / "r" / cli.SOURCE_SPOOL))
    doc = json.load(open(str(tmp_path / "r.json")))
    assert doc["run_name"] == "r" and doc["settings"]["t_max_s"] == 0.2125
    assert doc["resumed_from"] == {"state": os.path.abspath(state_path), "step": 14, "time_s": pytest.approx(0.1875),
                                   "run_name": "f", "changed": {}}
    assert "resumed_from" not in json.load(open(os.path.join(out, "f.json")))
    assert doc["results"]["end_reason"] == "t_max" and "state_saved_at_s" not in doc["results"]


def test_the_resume_command_changes_what_it_is_told(flown, tmp_path):
    """Each option changes the resumed flight and is recorded in the JSON's "changed": --t-max flies further (two more
    12.5 ms steps than the saved flight's end), a halved --dt-continuum takes twice the steps over the same time (4 from
    0.1875 s to 0.2125 s instead of 2), and --seed seeds numpy's generator instead of resuming the saved random state, so
    the multigrid hierarchy the heat solver builds in the 16th step is another and the history differs from the flight's."""
    out, flight, records, state_path = flown
    run = lambda name, *options: (cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", name, *options]),
                                  rows_of(str(tmp_path / (name + ".csv"))), json.load(open(str(tmp_path / (name + ".json")))))
    code, longer, doc = run("longer", "--t-max", "0.2375")
    assert code == 0 and len(longer) == 5 and longer[:2] == flight[14:16] and doc["resumed_from"]["changed"] == {"--t-max": 0.2375}
    assert doc["results"]["final_time_s"] == pytest.approx(0.2375)       # (its 3rd row is mid-flight, the flight's the last: only its momentum_residual differs)
    code, finer, doc = run("finer", "--dt-continuum", "0.00625")
    assert code == 0 and len(finer) == 5 and doc["resumed_from"]["changed"] == {"--dt-continuum": 0.00625}
    assert doc["results"]["final_time_s"] == pytest.approx(0.2125) and doc["settings"]["thermal"] == "fem"
    code, reseeded, doc = run("reseeded", "--seed", "7")
    assert code == 0 and len(reseeded) == 3 and reseeded != flight[14:] and doc["resumed_from"]["changed"] == {"--seed": 7}
    assert doc["settings"]["seed"] == 7
    code, same, doc = run("same")
    assert same == flight[14:] and doc["settings"]["seed"] == cli.DEFAULT_SEED


def test_resuming_with_a_limit_adds_the_spreading_columns_and_off_drops_them(flown, tmp_path):
    """The rows of a run all take the first row's columns, the saved one's in a resume, so a limit given to a resume adds the
    spreading columns to the saved row, or the later rows' would be lost, and `off` drops them from a state saved with a
    limit, or the later rows would raise KeyError (sub-plan 22 Task 5c). The test flight's state saved without a limit, as the
    two real states were, resumed with one: the columns are in every row, zero in the saved one, the JSON records the limit in
    its settings, its results and what was changed, and every other column is the flight's, the limit binding nowhere. Then
    the flight saved with a limit: the resume keeps it unless told, and `off` gives back the flight's column set."""
    out, flight, records, state_path = flown
    names = ["spread_mass_kg", "spread_unplaced_kg", "spread_edges_over", "spread_jump_max_K"]
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "on", "--spread-limit", "60"]) == 0
    rows = rows_of(str(tmp_path / "on.csv"))
    assert len(rows) == 3 and list(rows[0]) == list(flight[14]) + names and all(list(r) == list(rows[0]) for r in rows)
    assert all(float(rows[0][k]) == 0.0 for k in names)
    assert [{k: v for k, v in r.items() if k not in names} for r in rows] == flight[14:]
    doc = json.load(open(str(tmp_path / "on.json")))
    assert doc["resumed_from"]["changed"] == {"--spread-limit": 60.0} and doc["settings"]["spread_limit_deg"] == 60.0
    assert doc["results"]["spread_limit_deg"] == 60.0 and doc["results"]["spread_mass_kg"] == 0.0

    limited = str(tmp_path / "limited")
    assert cli.main(FLIGHT + ["--t-max", "0.2125", "--outdir", limited, "--name", "f", "--save-state-at", "0.1875",
                              "--spread-limit", "60"]) == 0
    saved = os.path.join(limited, "f.state.pkl")
    doc = json.load(open(os.path.join(limited, "f.json")))
    assert doc["settings"]["skin"]["spread_limit_deg"] == 60.0 and doc["results"]["spread_limit_deg"] == 60.0
    assert list(rows_of(os.path.join(limited, "f.csv"))[0]) == list(flight[0]) + names

    assert cli.main(["resume", saved, "--outdir", str(tmp_path), "--name", "kept"]) == 0
    assert list(rows_of(str(tmp_path / "kept.csv"))[0]) == list(flight[0]) + names
    assert json.load(open(str(tmp_path / "kept.json")))["results"]["spread_limit_deg"] == 60.0

    assert cli.main(["resume", saved, "--outdir", str(tmp_path), "--name", "off", "--spread-limit", "off"]) == 0
    assert rows_of(str(tmp_path / "off.csv")) == flight[14:]
    doc = json.load(open(str(tmp_path / "off.json")))
    assert doc["resumed_from"]["changed"] == {"--spread-limit": "off"}
    assert "spread_limit_deg" not in doc["results"] and "spread_mass_kg" not in doc["results"]
    assert "spread_limit_deg" not in doc["settings"]


@pytest.mark.parametrize("case", ["missing", "format", "t_max", "dt_continuum", "t_max_nan", "dt_continuum_nan",
                                  "spread_0", "spread_180", "spread_nan", "save_not_after", "save_at_end", "save_nan"])
def test_the_resume_command_refuses_what_it_cannot_resume(case, flown, tmp_path):
    """A missing state file, a state of another format, an end of flight not after the saved time, a continuum step outside
    (0, the saved flight's --dt of 0.025 s], a spread limit outside (0, 180) degrees and a value that is not a number are
    refused with exit 2 (bad input), and nothing is written."""
    out, flight, records, state_path = flown
    bad_format = str(tmp_path / "other.pkl")
    restart.save(bad_format, {"format": restart.FORMAT + 1})
    state, options = {"missing": (str(tmp_path / "none.pkl"), []), "format": (bad_format, []),
                      "t_max": (state_path, ["--t-max", "0.1875"]), "dt_continuum": (state_path, ["--dt-continuum", "0.05"]),
                      "t_max_nan": (state_path, ["--t-max", "nan"]), "dt_continuum_nan": (state_path, ["--dt-continuum", "nan"]),
                      "spread_0": (state_path, ["--spread-limit", "0"]), "spread_180": (state_path, ["--spread-limit", "180"]),
                      "spread_nan": (state_path, ["--spread-limit", "nan"]),
                      "save_not_after": (state_path, ["--save-state-at", "0.1875"]),
                      "save_at_end": (state_path, ["--save-state-at", "0.2125"]),
                      "save_nan": (state_path, ["--save-state-at", "nan"])}[case]
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", *options])
    assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))


def test_a_spread_limit_is_refused_on_a_flight_that_does_not_melt(tmp_path):
    """A flight that does not melt has nothing re-frozen to spread, so a limit given to its resume exits 2, as `run` refuses
    it; `off` asks nothing of the flight, so that goes through."""
    pytest.importorskip("cantera")
    k = FLIGHT.index("--melt")
    argv = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--quiet"]              # without --melt on and --skin-flow-cycles 2
    assert cli.main(argv + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "dry", "--save-state-at", "0.05"]) == 0
    state = str(tmp_path / "dry.state.pkl")
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--spread-limit", "60"])
    assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
    assert cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--spread-limit", "off"]) == 0


def test_a_resumed_flight_saves_its_state_on_its_way_and_that_state_resumes(flown, tmp_path):
    """A resumed flight can save its state as the first did (final review, 2026-11-07): `--save-state-at 0.2` on the resume
    of the state saved at 0.1875 s writes <outdir>/<name>.state.pkl at the end of the 15th step, though the original
    flight's own save had set the run's record of one; the JSON records the time and the file, the new state carries the
    resumed run's name and arguments, and resumed in its turn it repeats the unbroken flight from its step to the last bit.
    The same flight resumed without the option saves nothing."""
    out, flight, records, state_path = flown
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "r", "--save-state-at", "0.2"]) == 0
    doc = json.load(open(str(tmp_path / "r.json")))
    new = str(tmp_path / "r.state.pkl")
    assert doc["settings"]["save_state_at_s"] == 0.2 and doc["results"]["state_saved_at_s"] == pytest.approx(0.2)
    assert doc["results"]["state_save_error"] is None and doc["files"]["state"] == new
    assert doc["resumed_from"]["changed"] == {"--save-state-at": 0.2}
    state = restart.load(new)
    assert state["step"] == 15 and state["meta"]["name"] == "r" and state["meta"]["args"]["save_state_at"] == 0.2
    assert state["meta"]["args"]["outdir"] == str(tmp_path)
    assert cli.main(["resume", new, "--outdir", str(tmp_path), "--name", "again"]) == 0
    assert rows_of(str(tmp_path / "again.csv")) == flight[15:]
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "plain"]) == 0
    assert not os.path.exists(str(tmp_path / "plain.state.pkl"))
    assert "state" not in json.load(open(str(tmp_path / "plain.json")))["files"]


def test_a_resume_says_what_it_is_doing_before_it_flies(flown, tmp_path, capsys):
    """One line at the start of a resume, before the flight, names the state, its saved time, the end of flight and every
    option in force, the spread limit or off included, so a flight's log says what it is without waiting for its JSON."""
    out, flight, records, state_path = flown
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "r", "--t-max", "0.225", "--seed", "7",
                     "--spread-limit", "60", "--save-state-at", "0.2"]) == 0
    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")]
    assert len(line) == 1
    for text in (os.path.abspath(state_path), "saved at 0.1875 s", "step 14", "to 0.225 s", "seed 7", "spread limit 60 deg",
                 "save state at 0.2 s", "dt_continuum 0.0125", "slurry spray off", os.path.join(str(tmp_path), "r")):
        assert text in line[0]
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "q"]) == 0
    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")][0]
    assert "saved random state" in line and "spread limit off" in line and "save state at no" in line


def test_resuming_with_the_melting_rule_sets_it_on_the_body_and_records_it(flown, tmp_path, capsys):
    """`--part-count-from melting` on a resume (sub-plan 22 Task 5d) sets the rule on the saved flight's skin body, which the
    one line before the flight names, records it in the JSON's settings and in what was changed, and carries it into the
    state the resumed flight saves, so a resume of that state keeps it unless told; `moves` asked of it gives back the
    saved flight's record (no key), and a flight resumed without the option records none."""
    out, flight, records, state_path = flown
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "m", "--part-count-from", "melting",
                     "--save-state-at", "0.2"]) == 0
    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")][0]
    assert "part count from melting" in line
    doc = json.load(open(str(tmp_path / "m.json")))
    assert doc["resumed_from"]["changed"] == {"--part-count-from": "melting", "--save-state-at": 0.2}
    assert doc["settings"]["part_count_from"] == "melting"
    assert restart.load(str(tmp_path / "m.state.pkl"))["meta"]["args"]["part_count_from"] == "melting"
    assert cli.main(["resume", str(tmp_path / "m.state.pkl"), "--outdir", str(tmp_path), "--name", "kept"]) == 0
    assert "part count from melting" in capsys.readouterr().out
    assert json.load(open(str(tmp_path / "kept.json")))["settings"]["part_count_from"] == "melting"
    assert cli.main(["resume", str(tmp_path / "m.state.pkl"), "--outdir", str(tmp_path), "--name", "back",
                     "--part-count-from", "moves"]) == 0
    doc = json.load(open(str(tmp_path / "back.json")))
    assert "part_count_from" not in doc["settings"] and doc["resumed_from"]["changed"] == {"--part-count-from": "moves"}
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "plain"]) == 0
    assert "part_count_from" not in json.load(open(str(tmp_path / "plain.json")))["settings"]
    assert cli.main(FLIGHT + ["--t-max", "0.05", "--outdir", str(tmp_path), "--name", "run", "--part-count-from", "melting"]) == 0
    assert json.load(open(str(tmp_path / "run.json")))["settings"]["skin"]["part_count_from"] == "melting"   # run's record


def test_the_melting_rule_is_refused_on_a_resume_of_a_flight_that_does_not_melt(tmp_path):
    """A flight that does not melt has no draws to set a part count by, so the melting rule asked of its resume exits 2,
    as `run` refuses it; `moves` asks nothing of the flight, so that goes through (sub-plan 22 Task 5d)."""
    pytest.importorskip("cantera")
    k = FLIGHT.index("--melt")
    argv = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--quiet"]              # without --melt on and --skin-flow-cycles 2
    assert cli.main(argv + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "dry", "--save-state-at", "0.05"]) == 0
    state = str(tmp_path / "dry.state.pkl")
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--part-count-from", "melting"])
    assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
    assert cli.main(["resume", state, "--outdir", str(tmp_path / "o"), "--name", "r", "--part-count-from", "moves"]) == 0
