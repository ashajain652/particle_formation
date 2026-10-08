"""Spheral M1 on a real finite-element flight (plan 2026-10-07-spheral-m1, Task 10): the thresholds as measured.

The flight is the 2026-10-08 Scheil rerun of the reconstructed Step 3 prototype (M1 decision 13), the MVP input:

    reentry_model_output/FLIGHT  (1,255 frames, k 0-1254, 0.5 s step, AA7075_scheil, seed 12345, molten cascade on,
                                  deep runoff off; package prototype/work-2026-10-08-spheral-mvp/code, sha256 9629f39c...)

prepared whole on 2026-10-08 (1,192 s of frame processing, 791 MB, median frame 0.51 MB) and measured by
`analysis/spheral_m1_flight.py` into data/spheral/m1_frames.csv. Every threshold below is that flight's worst case
rounded up to one significant figure, or the a priori one of Tasks 3-8 where that is tighter; the comment says which.

`test_committed_measurement_meets_the_thresholds` reads only the committed CSV and runs in the normal loop. The
`fe_flight` tests need the flight itself (git-ignored): set SPHERAL_FRAG_FE_RUN to its run directory and
SPHERAL_FRAG_FE_PACKAGE to the package, e.g.

    SPHERAL_FRAG_FE_RUN=reentry_model_output/FLIGHT SPHERAL_FRAG_FE_PACKAGE=prototype/work-2026-10-08-spheral-mvp/code \\
        "$PY" -m pytest tests/test_spheral_frag_flight.py -q

They prepare six frames afresh (frame 0, the first film, the worst stagnation mismatch, the thinnest patch, peak
dynamic pressure, the last frame) and require the committed CSV's rows back, then run a synthetic Spheral run on the
hot-phase interval and close its accounts.

The drag does not meet spec §7.3's "a few percent" (5 %): the frame's own loads on the derived normals differ from the
history's drag by +10 % before the first film, -36 % to +9 % in the shock-layer phase and a constant -11 % on the
late subsonic branch. That is M1's answer (Review focus 6), reported, not fitted away; the threshold only pins it."""
import csv
import importlib.util
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from helpers import REPO_ROOT

FLIGHT = ("model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-AA7075_scheil_"
          "deeprunoff-off")
PACKAGE_SHA256 = "9629f39c77d1f76c317b8640ae046ae88a0fd576f29fee671df16009b9c3687f"
CSV_PATH = os.path.join(REPO_ROOT, "data", "spheral", "m1_frames.csv")
N_FRAMES = 1255
SAMPLE_FRAMES = (0, 48, 64, 166, 182, 1254)

# measured on FLIGHT, 2026-10-08 (worst case in brackets)
THRESHOLDS = {
    "n_open_directed_edges": 0,          # a priori (plan fact 7) [0]
    "n_nonmanifold_edges": 300,          # measured [288, on 1,204 frames]: the staircase surface pinches
    "n_nonmanifold_vertices": 20,        # measured [15, on 131 frames]
    "volume_rel_diff": 1e-12,            # a priori, round-off [2.2e-16]
    "mass_rel_diff": 5e-9,               # a priori, the history's 9 digits [4.2e-9]
    "film_rel_diff": 5e-9,               # measured, the history's 9 digits [4.1e-9]
    "max_abs_dfl": 0.0,                  # a priori, bitwise [0]
    "n_delta_m_mismatch": 0,             # a priori [0]
    "n_nan_forbidden": 0,                # a priori [0]
    "rel_table_patch": 8e-4,             # measured, binning by 1 deg [7.85e-4]
    "rel_smooth_hist": 0.4,              # measured [0.390 at 192.5 s]; spec §7.3's 5 % is NOT met (see above)
    "rel_p_w_max_stag": 0.2,             # measured [0.129 at k 64, where the stagnation patch died; median 6e-10]
    "release_step_low": -0.3,            # measured: sum release_rate A / sprayed increment - 1 >= [-0.295]
    "release_step_high": 0.04,           # measured: <= [+0.038]
    "release_total": 0.2,                # measured: |flight total| [0.123: 1.025 of 1.169 kg]
    "frame0_thickness_vs_2R": 6e-4,      # measured, facetted sphere [5.8e-4]
    "thickness_nan": 0,                  # a priori [0]
    "n_wrong_way": 0,                    # a priori [0]
    "accounts": 1e-14,                   # a priori (plan Task 8) [4e-17 on frames 100-110 and 800-810]
}
FRAMES_WITHOUT_LOADS = [0]               # flow_eval 0 on frame 0 only (decision 4)


def _rows():
    with open(CSV_PATH, newline="") as fh:
        return list(csv.DictReader(fh))


def _col(rows, name):
    return np.array([np.nan if r[name] == "" else float(r[name]) for r in rows])


def _analysis():
    spec = importlib.util.spec_from_file_location("spheral_m1_flight",
                                                  os.path.join(REPO_ROOT, "analysis", "spheral_m1_flight.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _check_rows(rows):
    T = THRESHOLDS
    for name in ("n_open_directed_edges", "n_nonmanifold_edges", "n_nonmanifold_vertices", "volume_rel_diff",
                 "mass_rel_diff", "film_rel_diff", "max_abs_dfl", "n_delta_m_mismatch", "n_nan_forbidden",
                 "thickness_nan", "n_wrong_way"):
        v = np.abs(_col(rows, name))
        assert np.nanmax(v) <= T[name], (name, float(np.nanmax(v)))
    assert all(r["fl_bitwise"] == "1" for r in rows)
    valid = _col(rows, "loads_valid") == 1
    assert [int(r["k"]) for r, ok in zip(rows, valid) if not ok] == [k for k in FRAMES_WITHOUT_LOADS
                                                                    if any(int(r["k"]) == k for r in rows)]
    for name in ("rel_table_patch", "rel_smooth_hist"):
        assert np.nanmax(np.abs(_col(rows, name)[valid])) <= T[name], name
    assert np.nanmax(np.abs(_col(rows, "rel_p_w_max_stag"))) <= T["rel_p_w_max_stag"]
    rs, di = _col(rows, "release_sum_kg"), _col(rows, "sprayed_increment_kg")
    spr = di > 0
    if spr.any():
        rel = rs[spr] / di[spr] - 1.0
        assert rel.min() >= T["release_step_low"] and rel.max() <= T["release_step_high"]
    assert not np.any((rs > 0) & (di <= 0))               # nothing released on a step the history sprayed nothing


def test_committed_measurement_meets_the_thresholds():
    rows = _rows()
    assert len(rows) == N_FRAMES and [int(r["k"]) for r in rows] == list(range(N_FRAMES))
    _check_rows(rows)
    rs, di = _col(rows, "release_sum_kg"), _col(rows, "sprayed_increment_kg")
    assert abs(rs.sum() / di.sum() - 1.0) <= THRESHOLDS["release_total"]
    assert np.all(np.diff(_col(rows, "time_s")) > 0)


# ---------------------------------------------------------------------------------------------------- the flight
@pytest.fixture(scope="module")
def prepared_sample(tmp_path_factory):
    run_dir = os.path.abspath(os.environ["SPHERAL_FRAG_FE_RUN"])
    package = os.path.abspath(os.environ["SPHERAL_FRAG_FE_PACKAGE"])
    assert os.path.basename(run_dir.rstrip(os.sep)) == FLIGHT, "the thresholds are this flight's"
    out = tmp_path_factory.mktemp("prepare")
    dirs = {}
    for k in SAMPLE_FRAMES:
        cmd = [sys.executable, "-m", "spheral_frag", "prepare", "--fe-run", run_dir, "--fe-package", package,
               "--frames", "{0}:{0}".format(k), "--outdir", str(out), "--quiet"]
        p = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
        assert p.returncode == 0, p.stdout + p.stderr
        (d,) = [os.path.join(out, n) for n in os.listdir(out) if n.endswith("_k{0:05d}-{0:05d}".format(k))]
        dirs[k] = d
    return run_dir, dirs


@pytest.mark.fe_flight
def test_flight_contract_and_gating_checks(prepared_sample):
    from spheral_frag import contract
    _, dirs = prepared_sample
    for k, d in dirs.items():
        with open(os.path.join(d, "prepare.json")) as fh:
            prep = json.load(fh)
        assert prep["fe_run"]["package"]["sha256"] == PACKAGE_SHA256
        assert prep["flight"]["material"] == "AA7075_scheil" and prep["flight"]["seed"] == 12345
        assert prep["contract"]["unconfirmed"] == [] and prep["contract"]["absent_optional"] == []
        assert not prep["contract"].get("missing")
        for name, c in prep["checks"].items():
            assert c["passed"] in (True, None), (k, name, c)
        e = prep["frames"][0]
        assert set(e["absent_optional"]) == set()
        assert e["n_derived"]["n_against_facet"] == 0 and e["n_derived"]["max_norm_error"] < 1e-15
    assert {f.key for f in contract.FE_FIELDS} >= {"n_derived", "flow_eval", "p_w_stag_step_Pa"}


@pytest.mark.fe_flight
def test_flight_frames_reproduce_the_committed_rows(prepared_sample):
    from spheral_frag import frames
    mod = _analysis()
    run_dir, dirs = prepared_sample
    history = frames.read_fe_run(run_dir).history
    committed = {int(r["k"]): r for r in _rows()}
    got = []
    for k, d in dirs.items():
        with open(os.path.join(d, "prepare.json")) as fh:
            (e,) = json.load(fh)["frames"]
        row = mod.frame_row(e, history, d)
        got.append(row)
        want = committed[k]
        for name in mod.COLUMNS:
            if name.startswith("t_") or name == "size_bytes":
                continue                                   # timings and the deflated size are not reproducible
            v = row.get(name)
            s = "" if v is None else "{:.9g}".format(v) if isinstance(v, float) else str(v)
            assert s == want[name], (k, name, s, want[name])
    _check_rows([{c: ("" if r.get(c) is None else
                      "{:.9g}".format(r[c]) if isinstance(r.get(c), float) else str(r[c])) for c in mod.COLUMNS}
                 for r in got])
    prep0 = dirs[0]
    thick = mod.frame0_thickness(prep0, [{"k": 0, "file": json.load(open(os.path.join(prep0, "prepare.json")))
                                          ["frames"][0]["file"]}], 0.1)
    assert thick["max_abs_rel"] <= THRESHOLDS["frame0_thickness_vs_2R"] and thick["n_nan"] == 0


@pytest.mark.fe_flight
def test_flight_synthetic_run_closes_its_accounts(tmp_path):
    """The hot-phase interval (frames 100-110, 50-55 s) prepared and run through fake_run and analyse, 3D, dx 2.2 mm."""
    import spheral_frag_synthetic as syn
    run_dir = os.path.abspath(os.environ["SPHERAL_FRAG_FE_RUN"])
    package = os.path.abspath(os.environ["SPHERAL_FRAG_FE_PACKAGE"])
    p = subprocess.run([sys.executable, "-m", "spheral_frag", "prepare", "--fe-run", run_dir, "--fe-package",
                        package, "--frames", "100:110", "--outdir", str(tmp_path / "prepare"), "--quiet"],
                       cwd=REPO_ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    (prepared,) = [str(x) for x in (tmp_path / "prepare").iterdir()]
    run, answers = syn.fake_run(prepared, 100, 110, 2.2e-3, seed=1, runs_dir=str(tmp_path / "runs"))
    p = subprocess.run([sys.executable, "-m", "spheral_frag", "analyse", "--run", run, "--prepared", prepared,
                        "--outdir", str(tmp_path / "analyse")], cwd=REPO_ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    with open(os.path.join(tmp_path, "analyse", os.path.basename(run), "analyse.json")) as fh:
        a = json.load(fh)
    assert a["passed"] and a["max_abs_rel_residual"] <= THRESHOLDS["accounts"]
    assert a["counts"]["n_fragments"] == 1 and a["counts"]["n_debris_rows"] == syn.FAKE_N_DEBRIS
    assert math.isclose(a["counts"]["fragment_mass_kg"], answers["ring_mass_kg"], rel_tol=1e-14)
