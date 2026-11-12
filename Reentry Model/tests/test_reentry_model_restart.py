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
# the dev skins on the 4 mm mesh; a 25 ms first step, then 12.5 ms ones, the tops spraying from the step ending 0.0625 s.
# With the options not given it is a flight of melting skins as they now default: the hand-down off, the part count by the melting
FLIGHT = ["run", "--diameter", "100", "--velocity", "7.24", "--altitude", "65", "--flight-path-angle", "-0.959331",
          "--atmosphere", "us76", "--thermal", "fem", "--heating", "physics", "--h-surface", "4", "--h-core", "20",
          "--melt", "on", "--prism-layers", "0", "--temperature", "893", "--dt", "0.025", "--skin-preset", "dev",
          "--fine-size", "2", "--skin-flow-cycles", "2", "--quiet"]
DEFAULT_FLIGHT = FLIGHT
# The old flight in its hand-down and its part count only (handing down and counting parts by the moves, as every flight before
# 2026-11-10 did), for the resume tests of the spread limit (spreading moves only handed-down metal), of the melting rule's
# change and of the state's `flown` saved options, which need those settings and not the new defaults (Task 5h); it flies
# `--spray-gates 2026-11-10`, the default, as every other flight here does
FLIGHT_ON = FLIGHT[:-3] + ["--hand-down", "on", "--part-count-from", "moves"] + FLIGHT[-3:]
# The element model's melting flight (`--surface-model elements`): the same 65 km, 893 K body and 25 ms first step, with no
# prism-layer, skin or fine-band option (the port's `element` case of analysis/skyscrapers_port/same_results.py)
MELTING_ELEMENTS = FLIGHT[:FLIGHT.index("--prism-layers")] + ["--temperature", "893", "--dt", "0.025",
                                                               "--surface-model", "elements", "--quiet"]


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
    assert cli.main(FLIGHT_ON + ["--t-max", "0.2125", "--outdir", out, "--name", "f", "--save-state-at", "0.1875"]) == 0
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
    assert cli.main(FLIGHT_ON + ["--t-max", "0.2125", "--outdir", limited, "--name", "f", "--save-state-at", "0.1875",
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


def test_resuming_with_the_hand_down_off_from_a_state_with_a_spread_limit_needs_the_limit_off(tmp_path, capsys):
    """Spreading moves only metal being handed down (Asha, 2026-11-08; sub-plan 22 Task 5e), so `--hand-down off` on a resume of
    a state saved with a spread limit (the spreading states carry 60) exits 2 unless `--spread-limit off` is given too, and
    flies with it: recorded in the JSON's settings and in what was changed, named in the line before the flight, and kept by
    a resume of the state it saves. The reverse, a limit given to a state saved with the hand-down off, exits 2 as well."""
    pytest.importorskip("cantera")
    limited = str(tmp_path / "limited")
    assert cli.main(FLIGHT_ON + ["--t-max", "0.075", "--outdir", limited, "--name", "f", "--save-state-at", "0.05",
                              "--spread-limit", "60"]) == 0
    state, out = os.path.join(limited, "f.state.pkl"), str(tmp_path / "o")
    capsys.readouterr()
    for options in (["--hand-down", "off"], ["--hand-down", "off", "--spread-limit", "60"]):
        with pytest.raises(SystemExit) as exc:
            cli.main(["resume", state, "--outdir", out, "--name", "r"] + options)
        assert exc.value.code == 2 and not os.path.exists(out)
        err = capsys.readouterr().err
        assert "spreading moves only metal being handed down" in err and "--spread-limit off" in err
    assert cli.main(["resume", state, "--outdir", out, "--name", "r", "--hand-down", "off", "--spread-limit", "off",
                     "--save-state-at", "0.0625"]) == 0
    assert "hand-down off" in capsys.readouterr().out
    doc = json.load(open(os.path.join(out, "r.json")))
    assert doc["resumed_from"]["changed"] == {"--hand-down": "off", "--spread-limit": "off", "--save-state-at": 0.0625}
    assert doc["settings"]["hand_down"] == "off" and "spread_limit_deg" not in doc["settings"]
    kept = os.path.join(out, "r.state.pkl")
    assert restart.load(kept)["meta"]["args"]["hand_down"] == "off"
    with pytest.raises(SystemExit) as exc:                  # a limit on a flight that hands nothing down
        cli.main(["resume", kept, "--outdir", out, "--name", "again", "--spread-limit", "60"])
    assert exc.value.code == 2 and "--hand-down on" in capsys.readouterr().err
    assert cli.main(["resume", kept, "--outdir", out, "--name", "same"]) == 0
    assert json.load(open(os.path.join(out, "same.json")))["settings"]["hand_down"] == "off"
    assert cli.main(["resume", kept, "--outdir", out, "--name", "back", "--hand-down", "on"]) == 0
    doc = json.load(open(os.path.join(out, "back.json")))
    assert "hand_down" not in doc["settings"] and doc["resumed_from"]["changed"] == {"--hand-down": "on"}


def test_resuming_with_height_pairing_sets_it_on_the_body_and_records_it(flown, tmp_path, capsys):
    """`--sideways-heat height` on a resume (sub-plan 22 Task 5f) sets the pairing on the saved flight's skin body, which pairs
    the heat otherwise from the next step on and which the one line before the flight names; it is recorded in the JSON's
    settings and in what was changed and carried into the state the resumed flight saves, so a resume of that state keeps
    it unless told, and `depth` asked of it gives back the saved flight's record (no key). `run` records it in the skins'
    settings only when given, and the resume of a flight without skins refuses it (exit 2)."""
    out, flight, records, state_path = flown
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "h", "--sideways-heat", "height",
                     "--save-state-at", "0.2"]) == 0
    line = [l for l in capsys.readouterr().out.splitlines() if l.startswith("resuming ")][0]
    assert "sideways heat height" in line
    doc = json.load(open(str(tmp_path / "h.json")))
    assert doc["resumed_from"]["changed"] == {"--sideways-heat": "height", "--save-state-at": 0.2}
    assert doc["settings"]["sideways_heat"] == "height"
    rows = rows_of(str(tmp_path / "h.csv"))
    assert rows[0] == flight[14] and rows[1:] != flight[15:]
    saved = restart.load(str(tmp_path / "h.state.pkl"))
    assert saved["meta"]["args"]["sideways_heat"] == "height" and saved["run"].body.sideways_pairing == "height"
    assert cli.main(["resume", str(tmp_path / "h.state.pkl"), "--outdir", str(tmp_path), "--name", "kept"]) == 0
    assert "sideways heat height" in capsys.readouterr().out
    assert json.load(open(str(tmp_path / "kept.json")))["settings"]["sideways_heat"] == "height"
    assert cli.main(["resume", str(tmp_path / "h.state.pkl"), "--outdir", str(tmp_path), "--name", "back",
                     "--sideways-heat", "depth"]) == 0
    doc = json.load(open(str(tmp_path / "back.json")))
    assert "sideways_heat" not in doc["settings"] and doc["resumed_from"]["changed"] == {"--sideways-heat": "depth"}
    assert "sideways_heat" not in json.load(open(os.path.join(out, "f.json")))["settings"]["skin"]
    assert cli.main(FLIGHT + ["--t-max", "0.05", "--outdir", str(tmp_path), "--name", "run", "--sideways-heat", "height"]) == 0
    assert json.load(open(str(tmp_path / "run.json")))["settings"]["skin"]["sideways_heat"] == "height"
    k = FLIGHT.index("--melt")
    elements = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--surface-model", "elements", "--quiet"]
    assert cli.main(elements + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "el", "--save-state-at", "0.05"]) == 0
    capsys.readouterr()
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", str(tmp_path / "el.state.pkl"), "--outdir", str(tmp_path / "o"), "--name", "r",
                  "--sideways-heat", "height"])
    assert exc.value.code == 2 and not os.path.exists(str(tmp_path / "o"))
    assert "--sideways-heat height needs a saved flight of skins" in capsys.readouterr().err


def test_resuming_with_size_feedback_current_fits_the_cap_and_records_it(flown, tmp_path, capsys):
    """`--size-feedback` on a resume (sub-plan 22 Task 5g, Asha 2026-11-10). The test flight is `current` by default (its
    JSON says so in the skins' settings); a resume keeps that, names it in the line before the flight and records it in
    the JSON's settings, and `initial` asked of it flies the mesh's 50 mm radius, recorded in what was changed with no
    settings key, and is kept by the state it saves. A state saved before Task 5g (no setting, no fit, no nose columns on
    its row) resumes as `initial` with the nose columns filled, and `current` asked of it fits the cap at the first load."""
    out, flight, records, state_path = flown
    assert json.load(open(os.path.join(out, "f.json")))["settings"]["skin"]["size_feedback"] == "current"
    assert all(k in flight[0] for k in ("nose_radius_mm", "fitted_nose_radius_mm", "nose_bound_mm", "nose_cap_points",
                                        "nose_cap_concave")) and float(flight[-1]["nose_radius_mm"]) != 50.0
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "same"]) == 0
    assert "size feedback current" in capsys.readouterr().out
    assert json.load(open(str(tmp_path / "same.json")))["settings"]["size_feedback"] == "current"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "i", "--size-feedback", "initial",
                     "--save-state-at", "0.2"]) == 0
    assert "size feedback initial" in capsys.readouterr().out
    doc = json.load(open(str(tmp_path / "i.json")))
    assert doc["resumed_from"]["changed"] == {"--size-feedback": "initial", "--save-state-at": 0.2}
    assert "size_feedback" not in doc["settings"]
    rows = rows_of(str(tmp_path / "i.csv"))
    assert all(float(r["nose_radius_mm"]) == 50.0 for r in rows)       # the saved row too (Task 5k, minor 2), not the old setting's
    assert {k: v for k, v in rows[0].items() if k != "nose_radius_mm"} == {k: v for k, v in flight[14].items() if k != "nose_radius_mm"}
    assert flight[14]["nose_radius_mm"] != rows[0]["nose_radius_mm"]
    assert restart.load(str(tmp_path / "i.state.pkl"))["run"].body.size_feedback == "initial"
    # a state saved before Task 5g
    state = restart.load(state_path)
    body = state["run"].body
    del vars(body)["size_feedback"], vars(body)["_nose_raw"]
    for key in ("nose_radius_mm", "fitted_nose_radius_mm", "nose_bound_mm", "nose_cap_points", "nose_cap_concave"):
        del state["row"][key]
    old = str(tmp_path / "old.state.pkl")
    restart.save(old, state)
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "o"]) == 0
    assert "size feedback initial" in capsys.readouterr().out
    rows = rows_of(str(tmp_path / "o.csv"))
    assert "size_feedback" not in json.load(open(str(tmp_path / "o.json")))["settings"]
    assert all(float(r["nose_radius_mm"]) == 50.0 and float(r["nose_cap_points"]) > 0 for r in rows)
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "c", "--size-feedback", "current"]) == 0
    assert "size feedback current" in capsys.readouterr().out
    doc = json.load(open(str(tmp_path / "c.json")))
    assert doc["resumed_from"]["changed"] == {"--size-feedback": "current"} and doc["settings"]["size_feedback"] == "current"
    rows = rows_of(str(tmp_path / "c.csv"))
    assert all(float(r["nose_radius_mm"]) == float(r["fitted_nose_radius_mm"]) != 50.0 for r in rows[1:])
    k = FLIGHT.index("--melt")
    elements = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--surface-model", "elements", "--quiet"]
    assert cli.main(elements + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "el", "--save-state-at", "0.05"]) == 0
    capsys.readouterr()
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", str(tmp_path / "el.state.pkl"), "--outdir", str(tmp_path / "o2"), "--name", "r",
                  "--size-feedback", "current"])
    assert exc.value.code == 2 and "--size-feedback needs a saved flight of skins" in capsys.readouterr().err


def test_resuming_with_a_leeward_test_keeps_the_saved_one_changes_it_and_reads_an_old_state_as_facing(flown, tmp_path, capsys):
    """`--leeward-test` on a resume (sub-plan 22 Task 5j, Asha 2026-11-10). The test flight is `centroid` by default; a resume
    keeps it, names it in the line before the flight and records it in the settings, and `facing` asked of it is recorded in
    what was changed and kept by the state it saves. A state saved before Task 5j (no setting on its body, no leeward
    columns on its row) resumes as `facing` with the two columns filled, and `centroid` asked of it is applied. The element
    model's saved flight refuses it (exit 2)."""
    out, flight, records, state_path = flown
    assert json.load(open(os.path.join(out, "f.json")))["settings"]["skin"]["leeward_test"] == "centroid"
    assert all(k in flight[0] for k in ("centroid_offset_mm", "leeward_corners"))
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "same"]) == 0
    assert "leeward test centroid" in capsys.readouterr().out
    assert json.load(open(str(tmp_path / "same.json")))["settings"]["leeward_test"] == "centroid"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "f", "--leeward-test", "facing",
                     "--save-state-at", "0.2"]) == 0
    assert "leeward test facing" in capsys.readouterr().out
    doc = json.load(open(str(tmp_path / "f.json")))
    assert doc["resumed_from"]["changed"] == {"--leeward-test": "facing", "--save-state-at": 0.2}
    assert "leeward_test" not in doc["settings"]
    assert restart.load(str(tmp_path / "f.state.pkl"))["run"].body.leeward_test == "facing"
    fresh = restart.load(state_path)["run"].body                  # Task 5k, minor 2: the saved row's nose and leeward columns are
    fresh.leeward_test = "facing"                                 # the saved step's surface read under the new test, not the old one's
    fresh._leeward_raw = fresh._nose_raw = None
    first = rows_of(str(tmp_path / "f.csv"))[0]
    assert all(float(first[k]) == pytest.approx(v, rel=1e-7) for k, v in fresh.nose_stats().items())      # the CSV keeps 9 digits
    state = restart.load(state_path)
    body = state["run"].body
    del vars(body)["leeward_test"], vars(body)["_leeward_raw"]
    for key in ("centroid_offset_mm", "leeward_corners"):
        del state["row"][key]
    n40, m40 = body.hist_n[:40].copy(), body.hist_m[:40].copy()     # and a body saved with 40 size bins, as before Task 5j: padded on resume
    assert body.hist_n[40:].sum() == 0.0 and n40.sum() > 0.0
    body.hist_n, body.hist_m = n40.copy(), m40.copy()
    old = str(tmp_path / "old.state.pkl")
    restart.save(old, state)
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "o", "--save-state-at", "0.2"]) == 0
    assert "leeward test facing" in capsys.readouterr().out
    after = restart.load(str(tmp_path / "o.state.pkl"))["run"].body       # the flight went on, sprayed, and kept the first 40 bins
    assert len(after.hist_n) == len(after.hist_m) == 50
    assert (after.hist_n[:40] >= n40).all() and (after.hist_m[:40] >= m40).all() and after.hist_n.sum() >= n40.sum()
    rows = rows_of(str(tmp_path / "o.csv"))
    assert "leeward_test" not in json.load(open(str(tmp_path / "o.json")))["settings"]
    assert all(k in r for r in rows for k in ("centroid_offset_mm", "leeward_corners"))
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "c", "--leeward-test", "centroid"]) == 0
    assert "leeward test centroid" in capsys.readouterr().out
    doc = json.load(open(str(tmp_path / "c.json")))
    assert doc["resumed_from"]["changed"] == {"--leeward-test": "centroid"} and doc["settings"]["leeward_test"] == "centroid"
    k = FLIGHT.index("--melt")
    elements = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--surface-model", "elements", "--quiet"]
    assert cli.main(elements + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "el", "--save-state-at", "0.05"]) == 0
    capsys.readouterr()
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", str(tmp_path / "el.state.pkl"), "--outdir", str(tmp_path / "o2"), "--name", "r",
                  "--leeward-test", "centroid"])
    assert exc.value.code == 2 and "--leeward-test needs a saved flight of skins" in capsys.readouterr().err


def test_a_melting_element_flight_saves_its_state_and_resumes_to_the_last_bit(tmp_path):
    """Defect 1 of the 2026-11-11 handoff: the element model's melting body (`body.MeltingBody`) held the film module as an
    attribute, which pickle cannot save, so its flight recorded `state_save_error` "TypeError: cannot pickle 'module' object"
    and saved nothing. The flight to 0.125 s twice, once saving at 0.05 s: the save succeeds and changes nothing in the
    flight, and the resumed flight repeats the unbroken one's rows from the saved step on, to the last bit."""
    pytest.importorskip("cantera")
    plain, saved = str(tmp_path / "plain"), str(tmp_path / "saved")
    assert cli.main(MELTING_ELEMENTS + ["--t-max", "0.125", "--outdir", plain, "--name", "f"]) == 0
    assert cli.main(MELTING_ELEMENTS + ["--t-max", "0.125", "--outdir", saved, "--name", "f", "--save-state-at", "0.05"]) == 0
    results = json.load(open(os.path.join(saved, "f.json")))["results"]
    assert results["state_save_error"] is None and results["state_saved_at_s"] == pytest.approx(0.05)
    flight = rows_of(os.path.join(plain, "f.csv"))
    assert rows_of(os.path.join(saved, "f.csv")) == flight
    assert cli.main(["resume", os.path.join(saved, "f.state.pkl"), "--outdir", str(tmp_path), "--name", "r"]) == 0
    resumed = rows_of(str(tmp_path / "r.csv"))
    assert len(flight) == 11 and len(resumed) == 8 and resumed == flight[-8:]


def test_a_non_melting_element_flight_resumes_to_the_last_bit(tmp_path):
    """Defect 2 of the 2026-11-11 handoff: a non-melting element flight's body is a plain `ThermalBody`, which has no
    `source_rows`, and `cmd_resume` read them unconditionally (AttributeError). Saved at 0.05 s of a 0.1 s flight, it
    resumes and repeats the unbroken flight's rows from the saved step on, to the last bit."""
    pytest.importorskip("cantera")
    k = FLIGHT.index("--melt")
    elements = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--surface-model", "elements", "--quiet"]
    assert cli.main(elements + ["--t-max", "0.1", "--outdir", str(tmp_path), "--name", "f", "--save-state-at", "0.05"]) == 0
    assert cli.main(["resume", str(tmp_path / "f.state.pkl"), "--outdir", str(tmp_path), "--name", "r"]) == 0
    flight, resumed = rows_of(str(tmp_path / "f.csv")), rows_of(str(tmp_path / "r.csv"))
    assert len(flight) == 5 and len(resumed) == 3 and resumed == flight[-3:]


def test_an_element_model_state_saved_with_40_size_bins_is_padded_to_50_on_resume(tmp_path):
    """The parked gap of Task 5j's re-review (Task 5k): the skins' state with 40-bin histogram accumulators resumes padded to
    the 50 bins (`cmd_resume`, the first 40 edges unchanged), and so must an element model's, whose melting body
    (`body.MeltingBody`) holds `hist_n` and `hist_m` as the skins do. A melting element flight's saved state is given 40-bin
    accumulators holding known counts, as saved before Task 5j, saved again and resumed: the state the resumed flight saves
    has them 50 long, the first 40 bins exactly as they were (the flight sprays nothing in this time) and the added ten
    empty."""
    pytest.importorskip("cantera")
    assert cli.main(MELTING_ELEMENTS + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "el",
                                        "--save-state-at", "0.05"]) == 0
    state = restart.load(str(tmp_path / "el.state.pkl"))
    body = state["run"].body
    assert not hasattr(body, "skins") and len(body.hist_n) == spray.N_BINS
    body.hist_n, body.hist_m = np.arange(1.0, 41.0), 1e-9 * np.arange(1.0, 41.0)         # as saved with 40 bins before Task 5j
    old = str(tmp_path / "old.state.pkl")
    restart.save(old, state)
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "o", "--t-max", "0.125",
                     "--save-state-at", "0.075"]) == 0
    after = restart.load(str(tmp_path / "o.state.pkl"))["run"].body
    assert after.n_released == 0.0
    assert len(after.hist_n) == len(after.hist_m) == spray.N_BINS == 50
    np.testing.assert_array_equal(after.hist_n[:40], np.arange(1.0, 41.0))
    np.testing.assert_array_equal(after.hist_m[:40], 1e-9 * np.arange(1.0, 41.0))
    assert after.hist_n[40:].sum() == 0.0 and after.hist_m[40:].sum() == 0.0


def test_a_state_saved_before_the_bottom_cell_count_resumes_with_it_continuing_from_the_saved_fault_count(flown, tmp_path):
    """Task 5k, fix round 1 (controller's ruling, 2026-11-11): a state saved before `jump_faults_bottom_cell` existed (the attribute
    deleted from its body, the column from its saved row, as the old-state tests of this file do) resumes, takes base moves, and
    reports the old-rule count as the saved `jump_faults` plus the old-rule faults of the new steps; `jump_faults` continues under the
    new rule from the resume. Before the fix the first base move raised AttributeError."""
    out, flight, records, state_path = flown
    state = restart.load(state_path)
    body = state["run"].body
    del vars(body)["jump_faults_bottom_cell"]
    del state["row"]["jump_faults_bottom_cell"]
    body.jump_faults = state["row"]["jump_faults"] = 7                 # a saved count, made by the old rule
    old = str(tmp_path / "old.state.pkl")
    restart.save(old, state)
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "o"]) == 0
    rows = rows_of(str(tmp_path / "o.csv"))
    assert len(rows) == 3 and rows[0]["jump_faults_bottom_cell"] == "7" and rows[0]["jump_faults"] == "7"
    last = rows[-1]
    assert float(last["jump_faults_bottom_cell"]) >= 7 and float(last["jump_faults"]) >= 7
    assert float(last["jump_faults_bottom_cell"]) - 7 >= float(last["jump_faults"]) - 7      # the new rule never counts more
    doc = json.load(open(str(tmp_path / "o.json")))["results"]
    assert doc["jump_faults_bottom_cell"] == int(float(last["jump_faults_bottom_cell"]))


def test_a_melting_skin_flight_with_no_options_hands_down_nothing_and_counts_parts_by_the_melting(tmp_path):
    """Asha, 2026-11-10 (sub-plan 22 Task 5h): with `--hand-down` and `--part-count-from` not given, a melting skin flight
    runs with the hand-down off and the part count by the melting, records both in the skins' settings and in the
    arguments its saved state carries, and its name carries no suffix for either. `--spray-gates` is `2026-11-10` and
    recorded; a resume keeps all three, and `--hand-down on` asked of it is the one that changes."""
    pytest.importorskip("cantera")
    assert "--hand-down" not in DEFAULT_FLIGHT and "--part-count-from" not in DEFAULT_FLIGHT
    assert cli.main(DEFAULT_FLIGHT + ["--t-max", "0.075", "--outdir", str(tmp_path), "--save-state-at", "0.05"]) == 0
    (doc_path,) = [p for p in os.listdir(str(tmp_path)) if p.endswith(".json")]
    doc = json.load(open(os.path.join(str(tmp_path), doc_path)))
    assert "handdown" not in doc["run_name"] and "parts-" not in doc["run_name"]
    assert doc["settings"]["skin"]["hand_down"] == "off" and doc["settings"]["skin"]["part_count_from"] == "melting"
    assert doc["settings"]["spray_gates"] == "2026-11-10"
    state_path = os.path.join(str(tmp_path), doc["run_name"] + ".state.pkl")
    state = restart.load(state_path)
    assert state["meta"]["args"]["hand_down"] == "off" and state["meta"]["args"]["part_count_from"] == "melting"
    body = state["run"].body
    assert body.hand_down is False and body.part_count_from == "melting" and body.spray.spray_gates == "2026-11-10"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "same"]) == 0
    doc = json.load(open(str(tmp_path / "same.json")))
    assert doc["settings"]["spray_gates"] == "2026-11-10" and doc["settings"]["part_count_from"] == "melting"
    assert "hand_down" in doc["settings"] and doc["settings"]["hand_down"] == "off"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "on", "--hand-down", "on"]) == 0
    doc = json.load(open(str(tmp_path / "on.json")))
    assert doc["resumed_from"]["changed"] == {"--hand-down": "on"} and "hand_down" not in doc["settings"]
    assert doc["settings"]["spray_gates"] == "2026-11-10" and doc["settings"]["part_count_from"] == "melting"


def test_resuming_with_spray_gates_before_sets_them_on_the_spray_model_and_a_state_without_them_resumes_as_before(
        flown, tmp_path, capsys):
    """`--spray-gates` on a resume (sub-plan 22 Task 5h). The test flight is `before` through its saved state only if the
    state carries no setting: a state saved before the option (the setting removed from its spray model, which then reads the
    class default) resumes as `before`, with `spray_gates: before` in the JSON's settings (written always, as `run` writes it:
    Task 5k, final review minor 1), named in the line before the
    flight; `2026-11-10` asked of it is recorded in the JSON's settings and in what was changed, and kept by the state the
    resumed flight saves. A flight that does not spray refuses the option."""
    out, flight, records, state_path = flown
    assert json.load(open(os.path.join(out, "f.json")))["settings"]["spray_gates"] == "2026-11-10"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "same"]) == 0
    assert "spray gates 2026-11-10" in capsys.readouterr().out
    assert json.load(open(str(tmp_path / "same.json")))["settings"]["spray_gates"] == "2026-11-10"
    state = restart.load(state_path)
    del vars(state["run"].body.spray)["spray_gates"]                       # a model pickled before the setting existed
    old = str(tmp_path / "old.state.pkl")
    restart.save(old, state)
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "o"]) == 0
    assert "spray gates before" in capsys.readouterr().out
    assert json.load(open(str(tmp_path / "o.json")))["settings"]["spray_gates"] == "before"      # recorded always (Task 5k, minor 1)
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "n", "--spray-gates", "2026-11-10",
                     "--save-state-at", "0.2"]) == 0
    doc = json.load(open(str(tmp_path / "n.json")))
    assert doc["resumed_from"]["changed"] == {"--spray-gates": "2026-11-10", "--save-state-at": 0.2}
    assert doc["settings"]["spray_gates"] == "2026-11-10"
    assert restart.load(str(tmp_path / "n.state.pkl"))["run"].body.spray.spray_gates == "2026-11-10"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "b", "--spray-gates", "before"]) == 0
    doc = json.load(open(str(tmp_path / "b.json")))
    assert doc["settings"]["spray_gates"] == "before" and doc["resumed_from"]["changed"] == {"--spray-gates": "before"}
    k = FLIGHT.index("--melt")
    dry = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--quiet"]
    assert cli.main(dry + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "dry", "--save-state-at", "0.05"]) == 0
    capsys.readouterr()
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", str(tmp_path / "dry.state.pkl"), "--outdir", str(tmp_path / "o2"), "--name", "r",
                  "--spray-gates", "before"])
    assert exc.value.code == 2 and "--spray-gates needs a saved flight with spraying" in capsys.readouterr().err


def test_resuming_with_the_slip_band_sets_it_on_the_spray_model_and_a_state_without_it_resumes_with_it_off(flown, tmp_path, capsys):
    """`--slip-band` on a resume (sub-plan 22 Task 5i). The test flight has it off, recorded as `off` in the JSON's settings of the
    flight and of a plain resume of it, `slip band off` in the line before the flight; `on` asked of it is recorded in the settings and in what
    was changed and kept by the state the resumed flight saves; a state saved before the option (the attribute removed from the
    spray model, which reads the class default) resumes with it off, recorded as such; a flight that does not spray refuses the option."""
    out, flight, records, state_path = flown
    assert json.load(open(os.path.join(out, "f.json")))["settings"]["slip_band"] == "off"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "same"]) == 0
    assert "slip band off" in capsys.readouterr().out
    assert json.load(open(str(tmp_path / "same.json")))["settings"]["slip_band"] == "off"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "on", "--slip-band", "on",
                     "--save-state-at", "0.2"]) == 0
    assert "slip band on" in capsys.readouterr().out
    doc = json.load(open(str(tmp_path / "on.json")))
    assert doc["resumed_from"]["changed"] == {"--slip-band": "on", "--save-state-at": 0.2} and doc["settings"]["slip_band"] == "on"
    saved = restart.load(str(tmp_path / "on.state.pkl"))
    assert saved["run"].body.spray.slip_band is True and saved["meta"]["args"]["slip_band"] == "on"
    assert cli.main(["resume", str(tmp_path / "on.state.pkl"), "--outdir", str(tmp_path), "--name", "kept"]) == 0
    assert json.load(open(str(tmp_path / "kept.json")))["settings"]["slip_band"] == "on"      # a resume keeps the saved setting
    state = restart.load(state_path)
    state["run"].body.spray.slip_band = True
    vars(state["run"].body.spray).pop("slip_band")                                           # pickled before the option existed
    old = str(tmp_path / "old.state.pkl")
    restart.save(old, state)
    assert cli.main(["resume", old, "--outdir", str(tmp_path), "--name", "o"]) == 0
    assert "slip band off" in capsys.readouterr().out and json.load(open(str(tmp_path / "o.json")))["settings"]["slip_band"] == "off"
    assert cli.main(["resume", state_path, "--outdir", str(tmp_path), "--name", "off", "--slip-band", "off"]) == 0
    assert json.load(open(str(tmp_path / "off.json")))["resumed_from"]["changed"] == {"--slip-band": "off"}
    k = FLIGHT.index("--melt")
    dry = FLIGHT[:k] + FLIGHT[k + 2:-3] + ["--quiet"]
    assert cli.main(dry + ["--t-max", "0.075", "--outdir", str(tmp_path), "--name", "dry", "--save-state-at", "0.05"]) == 0
    capsys.readouterr()
    with pytest.raises(SystemExit) as exc:
        cli.main(["resume", str(tmp_path / "dry.state.pkl"), "--outdir", str(tmp_path / "o2"), "--name", "r", "--slip-band", "on"])
    assert exc.value.code == 2 and "--slip-band needs a saved flight with spraying" in capsys.readouterr().err
