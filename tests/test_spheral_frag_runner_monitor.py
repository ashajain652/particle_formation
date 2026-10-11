"""spheral_frag.runner.run replay: the crossings-error monitor (Spheral M2 plan, decision 7; Task 6b).

On the monitored intervals (the first, and every --monitor-every M-th after it) the replay reruns the accepted
interval from its checkpoint at twice the crossings and estimates the error at N as
2 max vM(sigma_N - sigma_2N) / max vM(sigma_N) (`run.crossings_error`) and its pressure part 2 max |dP| / max
vM(sigma_N) (`run.pressure_error`); when the larger of the two (`run.limit_error`) exceeds --monitor-limit it takes the
2 N state and continues at 2 N. On the synthetic receding sphere of tests/test_spheral_frag_runner_replay.py (four intervals):

  - `--monitor-every 2` monitors intervals 0, 1 and 3, and leaves every accepted state, check, removal and
    intervals.csv value bitwise identical to `--monitor-every 0` (the estimates are below the limit);
  - `--monitor-limit 1e-12` forces the switch: the switched interval's accepted state equals a plain run at 2 N
    bitwise, the following interval runs at 2 N and is monitored against 4 N, the run name carries the suffixes,
    and summary.json lists the switches; with --max-crossings too low to double again the run fails (exit 1,
    failed.json);
  - a monitored run stopped and resumed equals the uninterrupted one bitwise, also across a switch (the run's
    crossings are carried in the checkpoints).

Measured 2026-10-10 on the synthetic sphere (10 mm spacing, 10 crossings against 20, intervals 0, 1, 3): 3D 6.0 %,
3.8 %, 5.6 % (RMS-based 1.2-1.9 %, p99 3.0-4.6 %), RZ 3.2 %, 3.8 %, 3.7 % (RMS 1.6-1.9 %, p99 3.1-3.7 %): at 10 mm,
five particles per radius, the coarse body exceeds the default 5 % in 3D, so the bitwise test runs at
--monitor-limit 0.5. Each rerun costs twice its interval (3D 3.5-4.3 s against 1.7-2.2 s; RZ 0.66-0.74 s against
0.33-0.38 s)."""
import json
import os

import numpy as np
import pytest

from helpers import REPO_ROOT
from test_spheral_frag_runner_replay import NPROC, inputs, intervals, launch, rel  # noqa: F401 (fixture)

from spheral_frag import naming, record
from spheral_frag.runner import checkpoint

MON_COLUMNS = ("monitored", "crossings_err_est", "crossings_err_rms", "crossings_err_p99", "crossings_err_pressure",
               "monitor_crossings", "monitor_wall_s", "monitor_ke_internal_ratio", "crossings_switch",
               "crossings_err_limit", "crossings_err_limit_basis")
WALL = ("wall_s", "wall_per_step_s", "monitor_wall_s")
LIMIT = 0.5     # above every estimate on the coarse synthetic sphere (measured: test below), so no switch


def run_dir_of(out):
    """The one run under `out`, at <out>/<prepared>/<run suffix> (naming.run_dir)."""
    names = os.listdir(out)
    assert len(names) == 1, names
    suffixes = os.listdir(os.path.join(out, names[0]))
    assert len(suffixes) == 1, suffixes
    d = os.path.join(out, names[0], suffixes[0])
    assert naming.run_dir(out, naming.run_name_from_dir(d)) == d
    return d


def same_checkpoints(a, b, ks, skip_extra=()):
    for k in ks:
        x = checkpoint.read(os.path.join(a, "checkpoints", checkpoint.filename(k)))
        y = checkpoint.read(os.path.join(b, "checkpoints", checkpoint.filename(k)))
        assert set(x) == set(y)
        for key in x:
            if key == "run":
                continue
            if key == "extra":
                for e in x[key]:
                    if e not in skip_extra:
                        assert np.array_equal(x[key][e], y[key][e]), (k, e)
            elif isinstance(x[key], np.ndarray):
                assert np.array_equal(x[key], y[key]), (k, key)
            else:
                assert x[key] == y[key], (k, key)


def same_records(a, b):
    ca, cb = record.list_checks(a), record.list_checks(b)
    assert [n for n, _ in ca] == [n for n, _ in cb]
    for (_, pa), (_, pb) in zip(ca, cb):
        with np.load(pa) as za, np.load(pb) as zb:
            for key in za.files:
                assert np.array_equal(za[key], zb[key]), (pa, key)
    ra, rb = (record.read_removed(os.path.join(r, "removed.npz")) for r in (a, b))
    for f in ("id", "t", "check", "mass", "h", "reason"):
        assert np.array_equal(getattr(ra, f), getattr(rb, f)), f


@pytest.fixture(scope="module")
def monitored(inputs):     # noqa: F811
    """form -> (unmonitored run, monitored run (--monitor-every 2), monitored stopped at 2 and resumed)."""
    d, prep, bodies = inputs
    cache = {}

    def get(form):
        if form not in cache:
            common = ["replay", "--body", rel(bodies[form]), "--prepared", rel(prep), "--frames", "0-4"]
            o0, o2, o2r = (os.path.join(d, f"mon_{form}_{t}") for t in ("off", "every2", "every2_resumed"))
            mon = ["--monitor-every", "2", "--monitor-limit", repr(LIMIT)]
            launch(*common, "--monitor-every", "0", "--out", rel(o0), nproc=NPROC[form])
            full = launch(*common, *mon, "--out", rel(o2), nproc=NPROC[form])
            launch(*common, *mon, "--out", rel(o2r), "--stop-at", "2", nproc=NPROC[form])
            launch(*common, *mon, "--out", rel(o2r), "--resume", nproc=NPROC[form])
            cache[form] = (run_dir_of(o0), run_dir_of(o2), run_dir_of(o2r), full.stdout)
        return cache[form]
    return get


@pytest.mark.spheral
@pytest.mark.parametrize("form", ("3d", "rz"))
def test_monitor_estimates_on_the_right_intervals_and_leaves_the_run_bitwise(monitored, form):
    """Estimates on intervals 0, 1 and 3 only, each with its rerun's record; every accepted state, check, removal
    and non-monitor column as in the unmonitored run, bitwise (measured estimates: module docstring)."""
    off, on, _, out = monitored(form)
    assert os.path.basename(off).endswith("_monitorevery-0")
    assert os.path.basename(on).endswith("_monitorevery-2_monitorlimit-0.5")
    rows_off, rows_on = intervals(off), intervals(on)
    print(form, [(r["k0"], r["crossings_err_est"], r["crossings_err_rms"], r["crossings_err_p99"], r["wall_s"],
                  r["monitor_wall_s"]) for r in rows_on])
    assert [int(r["monitored"]) for r in rows_on] == [1, 1, 0, 1]
    assert [int(r["monitored"]) for r in rows_off] == [0, 0, 0, 0]
    for r in rows_on:
        if r["monitored"] == "1":
            est, rms, p99 = (float(r[c]) for c in ("crossings_err_est", "crossings_err_rms", "crossings_err_p99"))
            assert 0.0 < rms <= p99 <= est < LIMIT, r
            pr, lim = float(r["crossings_err_pressure"]), float(r["crossings_err_limit"])
            assert lim == max(est, pr) < LIMIT
            assert r["crossings_err_limit_basis"] == ("von_mises" if est >= pr else "pressure")
            assert int(r["monitor_crossings"]) == 20 and float(r["monitor_wall_s"]) > 0.0
            assert int(r["crossings_switch"]) == 0
            with np.load(os.path.join(on, "monitor", "interval_{:05d}.npz".format(int(r["k0"])))) as z:
                assert float(z["est"]) == float(r["crossings_err_est"])
                assert float(z["est"]) == pytest.approx(2.0*z["vm_diff"].max()/z["vm_n"].max(), rel=1e-15)
        else:
            assert r["crossings_err_est"] == "nan" and int(r["monitor_crossings"]) == 0
            assert r["crossings_err_limit"] == "nan" and r["crossings_err_limit_basis"] == ""
    assert all(int(r["run_crossings"]) == 10 and int(r["crossings"]) == 10 for r in rows_on)
    # the accepted states are the unmonitored run's, bitwise
    same_checkpoints(off, on, range(5))
    same_records(off, on)
    for x, y in zip(rows_off, rows_on):
        for c in x:
            if c not in WALL and c not in MON_COLUMNS:
                assert x[c] == y[c], c
    with open(os.path.join(on, "summary.json")) as fh:
        mon = json.load(fh)["crossings_monitor"]
    assert mon["intervals_monitored"] == [[0, 1], [1, 2], [3, 4]] and mon["switches"] == []
    assert mon["max_crossings_err_est"] == max(float(r["crossings_err_est"]) for r in rows_on if r["monitored"] == "1")
    assert mon["max_crossings_err_limit"] == max(float(r["crossings_err_limit"]) for r in rows_on
                                                 if r["monitored"] == "1")
    assert mon["max_crossings_err_limit_basis"] in ("von_mises", "pressure")
    assert mon["final_crossings"] == 10


@pytest.mark.spheral
@pytest.mark.parametrize("form", ("3d", "rz"))
def test_monitored_run_stopped_and_resumed_is_bitwise(monitored, form):
    _, on, resumed, _ = monitored(form)
    same_checkpoints(on, resumed, range(5))
    same_records(on, resumed)
    for x, y in zip(intervals(on), intervals(resumed)):
        for c in x:
            if c not in WALL:
                assert x[c] == y[c], c
    for k in (0, 1, 3):
        with np.load(os.path.join(on, "monitor", f"interval_{k:05d}.npz")) as za, \
                np.load(os.path.join(resumed, "monitor", f"interval_{k:05d}.npz")) as zb:
            for key in za.files:
                assert np.array_equal(za[key], zb[key]), (k, key)


@pytest.mark.spheral
def test_a_tiny_limit_switches_to_twice_the_crossings(inputs):    # noqa: F811
    """--monitor-limit 1e-12 on frames 0-2, --monitor-every 2 (intervals 0 and 1 monitored): interval 0 is accepted
    at 20 crossings (the monitor's rerun), interval 1 runs at 20 and is switched to 40. The switched interval 0's
    checkpoint equals a plain 20-crossing run's bitwise; a run stopped after interval 0 and resumed continues at 20
    and equals the uninterrupted one; with --max-crossings 40 the second switch is impossible and the run fails."""
    d, prep, bodies = inputs
    common = ["replay", "--body", rel(bodies["rz"]), "--prepared", rel(prep)]
    tiny = ["--monitor-every", "2", "--monitor-limit", "1e-12"]
    o_t, o_tr, o_20, o_f = (os.path.join(d, f"tiny_{t}") for t in ("full", "resumed", "plain20", "fail"))
    r = launch(*common, "--frames", "0-2", *tiny, "--out", rel(o_t))
    assert "takes the 20-crossing state" in r.stdout and "continues at 40 crossings" in r.stdout
    launch(*common, "--frames", "0-2", *tiny, "--out", rel(o_tr), "--stop-at", "1")
    launch(*common, "--frames", "0-2", *tiny, "--out", rel(o_tr), "--resume")
    launch(*common, "--frames", "0-1", "--crossings", "20", "--monitor-every", "0", "--out", rel(o_20))
    run_t, run_tr, run_20 = run_dir_of(o_t), run_dir_of(o_tr), run_dir_of(o_20)
    name = naming.run_name_from_dir(run_t)
    assert name.endswith("_monitorevery-2_monitorlimit-1e-12")
    assert naming.parse_run_name(name)["brackets"]["monitor_limit"] == 1e-12
    rows = intervals(run_t)
    assert [int(x["crossings"]) for x in rows] == [20, 40]
    assert [int(x["run_crossings"]) for x in rows] == [20, 40]
    assert [int(x["crossings_switch"]) for x in rows] == [1, 1]
    assert [int(x["monitor_crossings"]) for x in rows] == [20, 40]
    assert all(float(x["crossings_err_limit"]) > 1e-12 for x in rows)
    # the switched interval's accepted state is the 2 N rerun: a plain run at 20 crossings, bitwise
    same_checkpoints(run_20, run_t, [1], skip_extra=("run_crossings",))
    ck = checkpoint.read(os.path.join(run_t, "checkpoints", checkpoint.filename(1)))
    assert ck["extra"]["run_crossings"] == 20 and ck["extra"]["crossings"] == 20
    with open(os.path.join(run_t, "summary.json")) as fh:
        mon = json.load(fh)["crossings_monitor"]
    assert [(s["k0"], s["to_crossings"]) for s in mon["switches"]] == [(0, 20), (1, 40)]
    assert all(s["limit"] > 1e-12 and s["basis"] in ("von_mises", "pressure") for s in mon["switches"])
    assert mon["final_crossings"] == 40
    # resumed across the switch
    same_checkpoints(run_t, run_tr, range(3))
    same_records(run_t, run_tr)
    for x, y in zip(rows, intervals(run_tr)):
        for c in x:
            if c not in WALL:
                assert x[c] == y[c], c
    # no room to double again: exit 1 with failed.json
    bad = launch(*common, "--frames", "0-4", "--monitor-every", "1", "--monitor-limit", "1e-12", "--max-crossings",
                 "40", "--out", rel(o_f), check=False)
    assert bad.returncode == 1 and "crossings error" in bad.stderr
    fd = run_dir_of(o_f)
    with open(os.path.join(fd, "failed.json")) as fh:
        f = json.load(fh)
    assert f["reason"] == "crossings_error" and (f["k0"], f["crossings"], f["monitor_crossings"]) == (1, 20, 40)
    assert f["crossings_err_limit"] == max(f["crossings_err_est"], f["crossings_err_pressure"])
    assert f["crossings_err_limit_basis"] in ("von_mises", "pressure")
    assert not os.path.isfile(os.path.join(fd, "summary.json"))


# ---------------------------------------------------------------------------------------------------------------
# Numpy helpers
# ---------------------------------------------------------------------------------------------------------------

def _state(pid, sigma):
    return {"pid": np.asarray(pid), "sigma": np.asarray(sigma, dtype=float)}


def test_crossings_error_definition():
    from spheral_frag.runner import run
    rng = np.random.default_rng(3)
    n = 50
    a = rng.normal(size=(n, 3, 3))
    a = 0.5*(a + a.transpose(0, 2, 1))*1e7
    b = a.copy()
    b[7] += np.diag([1e6, 0.0, 0.0])                      # one particle differs by a uniaxial 1 MPa: vM = 1 MPa
    est, rms, p99 = run.crossings_error(_state(np.arange(n), a), _state(np.arange(n), b))
    top = run.von_mises(a).max()
    assert est == pytest.approx(2.0*1e6/top, rel=1e-12)
    assert rms == pytest.approx(2.0*1e6/np.sqrt(n)/top, rel=1e-12)
    assert 0.0 < p99 < est
    # matched by pid, not by row
    perm = rng.permutation(n)
    assert run.crossings_error(_state(np.arange(n), a), _state(np.arange(n)[perm], b[perm])) == (est, rms, p99)
    # a pure pressure difference does not enter (von Mises of the difference)
    c = a - 3e6*np.eye(3)[None]
    assert run.crossings_error(_state(np.arange(n), a), _state(np.arange(n), c))[0] == pytest.approx(0.0, abs=1e-12)
    # equal states: zero; zero stress with a difference: inf; different particles: an error
    assert run.crossings_error(_state(np.arange(n), a), _state(np.arange(n), a)) == (0.0, 0.0, 0.0)
    z = np.zeros((2, 3, 3))
    assert run.crossings_error(_state([1, 2], z), _state([1, 2], z)) == (0.0, 0.0, 0.0)
    assert run.crossings_error(_state([1, 2], z), _state([1, 2], z + np.diag([1.0, 0, 0])))[0] == np.inf
    with pytest.raises(ValueError):
        run.crossings_error(_state([1, 2], z), _state([1, 3], z))
    # von Mises: uniaxial s -> |s|; pure shear t -> sqrt(3) t
    assert run.von_mises(np.diag([5.0, 0, 0])[None])[0] == pytest.approx(5.0)
    sh = np.zeros((1, 3, 3))
    sh[0, 0, 1] = sh[0, 1, 0] = 2.0
    assert run.von_mises(sh)[0] == pytest.approx(2.0*np.sqrt(3.0))


def test_monitor_schedule_and_names():
    from spheral_frag.runner import run
    assert [j for j in range(25) if run.monitor_due(j, 10)] == [0, 9, 19]
    assert [j for j in range(6) if run.monitor_due(j, 2)] == [0, 1, 3, 5]
    assert [j for j in range(4) if run.monitor_due(j, 1)] == [0, 1, 2, 3]
    assert not any(run.monitor_due(j, 0) for j in range(20))
    assert naming.BRACKET_DEFAULTS["monitor_every"] == 10 and naming.BRACKET_DEFAULTS["monitor_limit"] == 0.05
    prep = naming.prepare_name("synthetic_receding_sphere_d100.00mm_seed12345", 0, 4)
    default = naming.run_name(prep, "replay", (0, 4), "rz", 2.2, naming.canonical_brackets({}), 1)
    assert default.endswith("_seed1")                                   # existing names unchanged
    b = naming.canonical_brackets({"crossings": 20, "monitor_every": 2, "monitor_limit": 1e-12})
    name = naming.run_name(prep, "replay", (0, 4), "rz", 2.2, b, 1)
    assert name.endswith("_crossings-20_monitorevery-2_monitorlimit-1e-12")
    assert naming.parse_run_name(name)["brackets"] == b
    # the summary's account
    rows = [dict(k0="0", k1="1", crossings="10", monitored="1", crossings_err_est="0.01", crossings_err_rms="0.001",
                 crossings_err_p99="0.005", monitor_crossings="20", monitor_wall_s="2.0", crossings_switch="0",
                 run_crossings="10"),
            dict(k0="1", k1="2", crossings="10", monitored="0", crossings_err_est="nan", crossings_err_rms="nan",
                 crossings_err_p99="nan", monitor_crossings="0", monitor_wall_s="0.0", crossings_switch="0",
                 run_crossings="10"),
            dict(k0="2", k1="3", crossings="20", monitored="1", crossings_err_est="0.07", crossings_err_rms="0.01",
                 crossings_err_p99="0.03", monitor_crossings="20", monitor_wall_s="3.0", crossings_switch="1",
                 run_crossings="20"),
            dict(k0="3", k1="4", crossings="20", monitored="1", crossings_err_est="0.01", crossings_err_rms="0.001",
                 crossings_err_p99="0.005", crossings_err_pressure="0.09", crossings_err_limit="0.09",
                 crossings_err_limit_basis="pressure", monitor_crossings="40", monitor_wall_s="4.0",
                 crossings_switch="1", run_crossings="40")]

    class A:
        monitor_every, monitor_limit, crossings = 2, 0.05, 10
    s = run.monitor_summary(rows, A)
    assert s["intervals_monitored"] == [[0, 1], [2, 3], [3, 4]] and s["max_crossings_err_est"] == 0.07
    assert s["max_crossings_err_limit"] == 0.09 and s["max_crossings_err_limit_basis"] == "pressure"
    assert s["max_at"] == [3, 4]
    assert s["switches"] == [{"k0": 2, "k1": 3, "to_crossings": 20, "est": 0.07, "limit": 0.07, "basis": "von_mises"},
                             {"k0": 3, "k1": 4, "to_crossings": 40, "est": 0.01, "limit": 0.09, "basis": "pressure"}]
    assert s["final_crossings"] == 40 and s["monitor_wall_s"] == 9.0


def test_pressure_error():
    from spheral_frag.runner import run
    a = np.zeros((3, 3, 3))
    a[:, 0, 0] = [1e7, 2e7, 4e7]                          # uniaxial: vM 10, 20, 40 MPa
    b = a - 2e6*np.eye(3)[None]                           # P differs by 2 MPa everywhere, vM unchanged
    pr = run.pressure_error(_state([0, 1, 2], a), _state([0, 1, 2], b))
    assert pr == pytest.approx(2*2e6/4e7, rel=1e-12)
    est = run.crossings_error(_state([0, 1, 2], a), _state([0, 1, 2], b))[0]
    assert est == pytest.approx(0.0, abs=1e-12)
    # the limit keys on the larger of the two: a pure pressure difference counts (controller, after Task 6b)
    assert run.limit_error(est, pr) == (pr, "pressure")
    assert run.limit_error(0.03, 0.01) == (0.03, "von_mises")
    assert run.limit_error(0.02, 0.02) == (0.02, "von_mises")
    assert run.limit_error(0.02, float("nan")) == (0.02, "von_mises")


def test_measured_table():
    """data/spheral/m2_monitor.csv (analysis/spheral_m2_monitor.py, measured 2026-10-10): the monitor's own estimate
    on check 4's RZ 2.2 mm interval 0 -> 1 equals the estimator applied to check 4's 10- and 20-crossing runs bitwise
    (the same two states); the real flight's estimates (RZ 2.2 mm, frames 0-20, every 2nd interval) stay below the
    5 % limit [largest 1.45 %, k 7 -> 8], so no switch, also under the larger-of-the-two rule; check 4's RZ 2.2 mm
    ramp interval exceeds it on its pressure part (6.1 %) and would switch."""
    import csv
    from spheral_frag.runner import run
    with open(os.path.join(REPO_ROOT, "data", "spheral", "m2_monitor.csv"), newline="") as fh:
        rows = list(csv.DictReader(fh))
    v = [r for r in rows if r["case"] == "check4" and r["form"] == "rz" and r["dx_mm"] == "2.2" and r["k1"] == "1"]
    assert len(v) == 1 and float(v[0]["est"]) == float(v[0]["monitor_run_est"])
    fl = [r for r in rows if r["case"] == "flight" and r["monitor_crossings"] != "0"]
    assert len(fl) == 11 and all(int(r["crossings"]) == 10 for r in fl)
    assert max(float(r["est"]) for r in fl) < 0.02          # [1.446e-2]
    # the limit keys on the larger of the von Mises and pressure estimates (controller, after Task 6b)
    for r in rows:
        if r["est"] != "nan":
            assert (float(r["limit"]), r["limit_basis"]) == run.limit_error(float(r["est"]), float(r["pressure"]))
    assert max(float(r["limit"]) for r in fl) < 0.02        # [1.446e-2, von Mises]
    assert float(v[0]["limit"]) > 0.05 and v[0]["limit_basis"] == "pressure"    # [6.12e-2]: the ramp would switch
