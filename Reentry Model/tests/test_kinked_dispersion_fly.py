"""analysis/kinked_dispersion/fly.py, the flight launcher that keeps the history as it goes, checked without a flight:
`CoupledRun.row` is replaced by a stand-in for the test."""
import csv

import pytest

from analysis.kinked_dispersion import fly
from reentry_model.skyscrapers import coupled


@pytest.fixture
def stand_in(monkeypatch):
    """CoupledRun.row replaced by a function returning a row of awkward floats; returns the list of rows it made."""
    made = []

    def row(self, t, y, a, loads):
        i = len(made)
        r = {"time_s": 0.1 * (i + 1), "mass_kg": 1.0 / 3.0 + i, "released_mass_kg": 1e-300 * (i + 1),
             "floor_like_share": 0.1 + 0.2, "extra": float(i)}
        made.append(r)
        return r

    monkeypatch.setattr(coupled.CoupledRun, "row", row)
    return made


def test_keep_writes_header_and_each_row_as_made(stand_in, tmp_path):
    undo = fly.keep(str(tmp_path), "t")
    try:
        run = coupled.CoupledRun.__new__(coupled.CoupledRun)
        path = tmp_path / "t.rows.csv"
        for i in range(3):
            run.row(0, None, None, None)
            with open(path, newline="") as fh:          # readable before the run ends
                rows = list(csv.DictReader(fh))
            assert len(rows) == i + 1
            assert list(rows[0]) == list(stand_in[0])
    finally:
        undo()


def test_progress_line_on_rows_0_20_40(stand_in, tmp_path):
    undo = fly.keep(str(tmp_path), "t")
    try:
        run = coupled.CoupledRun.__new__(coupled.CoupledRun)
        for _ in range(41):
            run.row(0, None, None, None)
    finally:
        undo()
    lines = (tmp_path / "t.progress").read_text().splitlines()
    assert [int(l.split()[1]) for l in lines] == [0, 20, 40]
    assert "floor-like share 0.300" in lines[0] and "released" in lines[0]
    assert "T_surface_max" not in lines[0]               # keys the row lacks are left out


def test_undo_restores_the_method_and_closes_the_file(stand_in, tmp_path):
    before = coupled.CoupledRun.row
    undo = fly.keep(str(tmp_path), "t")
    assert coupled.CoupledRun.row is not before
    coupled.CoupledRun.__new__(coupled.CoupledRun).row(0, None, None, None)
    handle = None
    for cell in coupled.CoupledRun.row.__closure__:      # the hook's open file
        if hasattr(cell.cell_contents, "get") and "fh" in cell.cell_contents:
            handle = cell.cell_contents["fh"]
    assert handle is not None and not handle.closed
    undo()
    assert coupled.CoupledRun.row is before
    assert handle.closed


def test_values_round_trip_bit_for_bit(stand_in, tmp_path):
    undo = fly.keep(str(tmp_path), "t")
    try:
        run = coupled.CoupledRun.__new__(coupled.CoupledRun)
        for _ in range(3):
            run.row(0, None, None, None)
    finally:
        undo()
    with open(tmp_path / "t.rows.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    for written, made in zip(rows, stand_in):
        assert {k: float(v) for k, v in written.items()} == made


def test_exit_codes(tmp_path, monkeypatch):
    """A bad argument comes back as argparse's own 2 (the real parser here), a normal return as its code, an ordinary
    exception as 1 after its traceback; the hook is undone each time."""
    before = coupled.CoupledRun.row
    bad = ["--melt", "off", "--dispersion", "table", "--diameter", "1", "--velocity", "1", "--altitude", "1"]
    assert fly.main([str(tmp_path), "bad"] + bad) == 2          # --dispersion is refused on a run that does not spray
    assert coupled.CoupledRun.row is before
    monkeypatch.setattr(fly.cli, "main", lambda argv: 7)
    assert fly.main([str(tmp_path), "ok"]) == 7

    def boom(argv):
        raise RuntimeError("model failure")
    monkeypatch.setattr(fly.cli, "main", boom)
    assert fly.main([str(tmp_path), "boom"]) == 1
    assert coupled.CoupledRun.row is before
