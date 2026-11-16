"""analysis/flights/fly.py, the flight launcher (ported from the prototype's planning copy on 2026-11-12): it runs the
package's command line in its own process with a hook on `CoupledRun.row` that only reads each history row, keeping every
row in <outdir>/<name>.rows.csv as it is made and a progress line every PROGRESS_EVERY rows in <outdir>/<name>.progress,
so a flight that stops early keeps its rows to its last completed macro step."""
import csv
import importlib.util
import json
import os

import pytest

from reentry_model import PROJECT_ROOT
from reentry_model.skyscrapers import body, coupled

FLY = os.path.join(PROJECT_ROOT, "analysis", "flights", "fly.py")
# the element model's non-melting flight on the command-line tests' coarse mesh, 2 s in four 0.5 s macro steps
ARGS = ["--diameter", "100", "--velocity", "7.5", "--altitude", "77.500133", "--flight-path-angle", "-0.959331",
        "--atmosphere", "us76", "--thermal", "fem", "--heating", "sesam", "--h-surface", "4", "--h-core", "20",
        "--surface-model", "elements", "--t-max", "2", "--quiet"]


@pytest.fixture
def fly():
    spec = importlib.util.spec_from_file_location("fly", FLY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def test_a_flight_keeps_every_history_row_as_it_is_made_and_a_progress_line(fly, tmp_path):
    """The rows the hook kept are the run's own history, row for row and column for column (the history CSV prints nine
    significant digits, the rows file every digit); the progress file has its line for row 0; the hook is undone after the
    flight. The resume form keeps the rows made after the saved step."""
    original = coupled.CoupledRun.row
    assert fly.main([str(tmp_path), "f"] + ARGS + ["--save-state-at", "1"]) == 0
    assert coupled.CoupledRun.row is original
    history, rows = read(tmp_path / "f.csv"), read(tmp_path / "f.rows.csv")
    assert len(history) == len(rows) == 5
    assert all(set(r) == set(h) and all("{:.9g}".format(float(r[k])) == h[k] for k in h) for r, h in zip(rows, history))
    progress = (tmp_path / "f.progress").read_text().splitlines()
    assert len(progress) == 1 and progress[0].startswith("step      0  t    0.000 s")
    assert fly.main([str(tmp_path), "r", "resume", str(tmp_path / "f.state.pkl")]) == 0
    assert coupled.CoupledRun.row is original
    history, rows = read(tmp_path / "r.csv"), read(tmp_path / "r.rows.csv")
    assert len(history) == 3 and len(rows) == 2           # the saved step's row is the state's, not made again
    assert all("{:.9g}".format(float(r[k])) == h[k] for r, h in zip(rows, history[1:]) for k in h)


def test_a_flight_that_raises_part_way_returns_1_and_keeps_the_rows_made_before(fly, tmp_path, monkeypatch, capsys):
    """An exception raised in a macro step (here the third macro step's heat solve) is printed with its traceback and gives
    exit code 1; the rows file holds the three rows made before the failure, and the hook is undone. Since Task 4 of the
    geometry stop (2026-11-12) the command line also writes the history to the last completed step, the same three rows,
    with end reason `failed`, and saves that step as the end state; before it, it wrote no history."""
    advance, calls = body.ThermalBody.advance, []

    def failing(self, *args, **kwargs):
        calls.append(1)
        if len(calls) == 3:
            raise ArithmeticError("a failed solve")
        return advance(self, *args, **kwargs)
    monkeypatch.setattr(body.ThermalBody, "advance", failing)
    original = coupled.CoupledRun.row
    assert fly.main([str(tmp_path), "f"] + ARGS) == 1
    assert coupled.CoupledRun.row is original
    assert "ArithmeticError: a failed solve" in capsys.readouterr().err
    assert [float(r["time_s"]) for r in read(tmp_path / "f.rows.csv")] == [0.0, 0.5, 1.0]
    assert [float(r["time_s"]) for r in read(tmp_path / "f.csv")] == [0.0, 0.5, 1.0]
    assert json.load(open(str(tmp_path / "f.json")))["results"]["end_reason"] == "failed"
    assert (tmp_path / "f.end.state.pkl").exists()
