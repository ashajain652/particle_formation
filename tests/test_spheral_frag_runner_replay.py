"""spheral_frag.runner.run replay: the mass-scaled replay driver (Spheral M2 plan, Task 5; spec §9.1, §11.3).

A synthetic receding sphere, five frames 0.5 s apart (four replay intervals), built here in drama_env with the M1
synthetic helpers (`tests/spheral_frag_synthetic.py`: the coarse 100 mm tetrahedral sphere, the synthetic material,
modified-Newtonian loads): the temperature is T = 300 + (T_wall - 300) (r / R)^2 with T_wall = RECEDE_T_WALL per
frame (below the synthetic solidus, 750 K: no melting), and the tetrahedra whose centroid lies upstream of
RECEDE_X[k] are dead from frame k on, so the nose recedes by 10 mm per frame from frame 2. It is prepared
(`spheral_frag prepare`) and turned into 10 mm bodies (`spheral_frag body`, 3D and RZ, frames 0-4), and replayed at
the defaults under Spheral by `spheral_frag/runner/run.py`:

  - the run goes through the four intervals and removes exactly the particles whose `remove_k` says, at that frame,
    with their body-stage masses and the enthalpy of their held temperature, into the film account;
  - `spheral_frag analyse` reads the run directory and closes its accounts to 1e-12 m0 (spec §11.3; the
    threshold is analyse's own ACCOUNTS_RTOL);
  - a run stopped after interval 2 (`--stop-at 2`) and resumed (`--resume`) equals the uninterrupted run bitwise in
    every check, every checkpoint array, the removals and every column of intervals.csv but the wall times;
  - the run name parses back to the run's settings.

Measured 2026-10-10. Synthetic, 10 mm: 3D 487 particles (105 removed: 7, 42, 56 at frames 2-4) on 2 processes,
RZ 39 rings (9 removed: 3, 1, 5) on 1; s 1.16e7-1.46e7 at 10 crossings, 289-333 steps per interval, KE/E at most
4.1e-3 (3D) and 3.5e-3 (RZ), no repeats; 6.7 ms and 1.1 ms per step; the accounts close to 3.8e-17 m0. The net
axial momentum at an interval's end is round-off in 3D (3e-14 of the interval's drag impulse) and 6 to 154 times
the impulse in RZ (SolidSPHRZ, Task 4). Real flight (MVP frames, k 0-4, defaults): RZ 2.2 mm (813 rings) s 9.63e6,
1,649-1,651 steps and 6.5-6.6 s per interval (4.0 ms/step, 1 process), KE/E 6.0e-3-8.2e-3 of which 87-97 % is
the centre-of-mass drift (internal KE/E 1.9e-4-9.9e-4), drag -2.5 % against the prepared tables' integral, max
displacement within an interval 0.018 dx; 3D 2.2 mm (49,173 particles, k 0-1, 18 processes) 1,609 steps in 178 s
(111 ms/step), KE/E 9.2e-4, drag -0.11 %, momentum 3e-15 of the impulse."""
import csv
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from conftest import spheral_launcher
from helpers import REPO_ROOT

from spheral_frag import naming, record
from spheral_frag.runner import checkpoint

RUN = "spheral_frag/runner/run.py"
RECEDE_TIMES = (0.0, 0.5, 1.0, 1.5, 2.0)
RECEDE_T_WALL = (300.0, 500.0, 600.0, 650.0, 700.0)     # K at r = R
RECEDE_X = (np.inf, np.inf, 0.040, 0.030, 0.020)        # m: tetrahedra with centroid x above it are dead
RECEDE_RUN_NAME = "synthetic_receding_sphere_d100.00mm_seed12345"
DX_MM = 10.0
FORMS = ("3d", "rz")
NPROC = {"3d": 2, "rz": 1}
WALL_COLUMNS = ("wall_s", "wall_per_step_s", "monitor_wall_s")


def receding_run(directory):
    """Write the synthetic receding sphere's finite-element run (module docstring); returns its directory."""
    points, tets = syn.sphere_mesh()
    cent = points[tets].mean(axis=1)
    frames, drags, states = [], [], []
    for k, (t, T_wall, x_dead) in enumerate(zip(RECEDE_TIMES, RECEDE_T_WALL, RECEDE_X)):
        state = syn.flight_state(t)
        T = syn.sphere_temperature(points, T_wall)
        tk = tets[cent[:, 0] <= x_dead]
        frame = syn.make_frame(k, t, points, tk, T, np.ones(len(tk)), p_stag=state["p_w_stag_step_Pa"],
                               evaluated=k > 0)
        frames.append(frame)
        states.append(state)
        if k == 0:
            ev = syn.make_frame(k, t, points, tk, T, np.ones(len(tk)), p_stag=state["p_w_stag_step_Pa"],
                                evaluated=True)
            drags.append(syn.frame_drag(ev))
        else:
            drags.append(syn.frame_drag(frame))
    history = syn.history_rows(frames, states, drags, [syn.sprayed_this_step(f) for f in frames])
    synthetic = {"device": "receding_sphere", "R_m": syn.SPHERE_R, "T_wall_K": list(RECEDE_T_WALL),
                 "dead_upstream_of_x_m": [None if not math.isfinite(x) else x for x in RECEDE_X]}
    doc = syn.run_json(RECEDE_RUN_NAME, 2e3*syn.SPHERE_R, syn.RHO*4.0/3.0*math.pi*syn.SPHERE_R**3,
                       syn.frame_mass(frames[0]), len(frames), synthetic)
    return syn.write_fe_run(directory, frames, history, doc)


def rel(p):
    return os.path.relpath(p, REPO_ROOT)


def launch(*args, nproc=1, timeout=1800, check=True):
    cmd = [spheral_launcher()] + (["-n", str(nproc)] if nproc > 1 else []) + [RUN, *args]
    r = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=timeout)
    if check:
        assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-4000:]
    return r


@pytest.fixture(scope="module")
def inputs():
    """(work dir, prepared dir, {form: body dir}) in spheral_output/tests (the container sees only the repo)."""
    pytest.importorskip("pyvista")
    from spheral_frag.__main__ import main
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")
    os.makedirs(base, exist_ok=True)
    d = tempfile.mkdtemp(prefix="runner_replay_", dir=base)
    fe_dir = receding_run(os.path.join(d, "fe"))
    tdir = os.path.join(d, "material")
    os.makedirs(tdir)
    syn.synthetic_material_table(os.path.join(tdir, "material_table.npz"))
    assert main(["prepare", "--fe-run", fe_dir, "--material-table", tdir, "--outdir", os.path.join(d, "prepare"),
                 "--quiet"]) == 0
    prep = os.path.join(d, "prepare", naming.prepare_name(RECEDE_RUN_NAME, 0, 4))
    bodies = {}
    for form in FORMS:
        assert main(["body", "--prepared", prep, "--frames", "0-4", "--form", form, "--dx-mm", str(DX_MM),
                     "--out", os.path.join(d, "body"), "--quiet"]) == 0
        bodies[form] = naming.body_dir(os.path.join(d, "body"), naming.body_name(os.path.basename(prep), 0, 4, form,
                                                                                 DX_MM))
    yield d, prep, bodies
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture(scope="module")
def runs(inputs):
    """form -> (uninterrupted run dir, stopped-and-resumed run dir), launched once each."""
    d, prep, bodies = inputs
    cache = {}

    def get(form):
        if form not in cache:
            common = ["replay", "--body", rel(bodies[form]), "--prepared", rel(prep), "--frames", "0-4"]
            full = launch(*common, "--out", rel(os.path.join(d, "full")), nproc=NPROC[form])
            stop = launch(*common, "--out", rel(os.path.join(d, "resumed")), "--stop-at", "2", nproc=NPROC[form])
            assert "stopped at frame 2" in stop.stdout
            resumed = launch(*common, "--out", rel(os.path.join(d, "resumed")), "--resume", nproc=NPROC[form])
            assert "resuming from" in resumed.stdout
            run_name = [ln for ln in full.stdout.splitlines() if ln.startswith("replay ")][0].split()[1]
            cache[form] = (naming.run_dir(os.path.join(d, "full"), run_name),
                           naming.run_dir(os.path.join(d, "resumed"), run_name), full.stdout, stop.stdout)
        return cache[form]
    return get


def only_run(out):
    """The one run directory under the runs directory `out` (<out>/<prepared>/<run suffix>, naming.run_dir)."""
    (prepared,) = os.listdir(out)
    (suffix,) = os.listdir(os.path.join(out, prepared))
    return os.path.join(out, prepared, suffix)


def intervals(run_dir):
    with open(os.path.join(run_dir, "intervals.csv"), newline="") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------------------------------------------
# Under Spheral
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_replay_runs_and_removes_exactly_the_body_stage_removals(inputs, runs, form):
    d, prep, bodies = inputs
    full, _, out, _ = runs(form)
    assert os.path.isfile(os.path.join(full, "summary.json"))
    with np.load(os.path.join(bodies[form], "body.npz")) as z:
        b = {k: z[k] for k in z.files}
    rows = intervals(full)
    assert [(int(r["k0"]), int(r["k1"])) for r in rows] == [(0, 1), (1, 2), (2, 3), (3, 4)]
    assert all(float(r["ke_internal_ratio"]) < 0.05 for r in rows)          # the test's key (KE_FRACTION)
    assert all(int(r["repeats"]) == 0 for r in rows)
    removed = record.read_removed(os.path.join(full, "removed.npz"))
    gone = b["remove_k"] >= 0
    assert gone.any() and (b["remove_k"] == -1).any()
    assert np.array_equal(np.sort(removed.id), np.sort(b["pid"][gone]))
    order = np.argsort(b["pid"])
    rows_of = order[np.searchsorted(b["pid"], removed.id, sorter=order)]
    assert np.array_equal(removed.t, b["frames_time_s"][b["remove_k"][rows_of]])     # removed at its frame
    assert np.array_equal(removed.mass, b["m"][rows_of])                             # physical, body-stage mass
    assert np.all(removed.reason == record.REASON_CENTRE_CROSSING)
    for k in range(1, 5):                                                            # per interval
        assert int(rows[k - 1]["n_removed"]) == int(np.count_nonzero(b["remove_k"] == k))
    ck = checkpoint.read(os.path.join(full, "checkpoints", checkpoint.filename(4)))
    assert np.array_equal(ck["pid"], np.sort(b["pid"][~gone]))
    assert ck["scale"] > 1e6                                                         # mass-scaled


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_analyse_closes_the_accounts(inputs, runs, form, tmp_path):
    from spheral_frag.__main__ import main
    import json
    d, prep, _ = inputs
    full, _, _, _ = runs(form)
    assert main(["analyse", "--run", full, "--prepared", prep, "--outdir", str(tmp_path), "--quiet"]) == 0
    name = naming.run_name_from_dir(full)
    with open(os.path.join(naming.analyse_dir(str(tmp_path), name), "analyse.json")) as fh:
        doc = json.load(fh)
    assert doc["passed"] and doc["max_abs_rel_residual"] <= 1e-12        # analyse.ACCOUNTS_RTOL
    assert doc["counts"]["n_checks"] == 5 and doc["counts"]["n_fragments"] == 0
    last = doc["accounts"][-1]
    assert last["mass_film_kg"] > 0.0 and last["n_removals"] == doc["counts"]["n_removals"]


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_stopped_and_resumed_equals_uninterrupted_bitwise(runs, form):
    full, resumed, _, _ = runs(form)
    for k in range(5):
        a = checkpoint.read(os.path.join(full, "checkpoints", checkpoint.filename(k)))
        b = checkpoint.read(os.path.join(resumed, "checkpoints", checkpoint.filename(k)))
        assert set(a) == set(b)
        for key in a:
            if key == "extra":
                assert a[key].keys() == b[key].keys()
                for e in a[key]:
                    assert np.array_equal(a[key][e], b[key][e]), (k, e)
            elif isinstance(a[key], np.ndarray):
                assert np.array_equal(a[key], b[key]), (k, key)
            else:
                assert a[key] == b[key], (k, key)
    ca, cb = record.list_checks(full), record.list_checks(resumed)
    assert [n for n, _ in ca] == [n for n, _ in cb] == list(range(5))
    for (_, pa), (_, pb) in zip(ca, cb):
        with np.load(pa) as za, np.load(pb) as zb:
            assert za.files == zb.files
            for key in za.files:
                assert np.array_equal(za[key], zb[key]), (pa, key)
    ra, rb = (record.read_removed(os.path.join(r, "removed.npz")) for r in (full, resumed))
    for f in ("id", "t", "check", "mass", "h", "reason"):
        assert np.array_equal(getattr(ra, f), getattr(rb, f)), f
    ia, ib = intervals(full), intervals(resumed)
    assert len(ia) == len(ib) == 4
    for x, y in zip(ia, ib):
        for c in x:
            if c not in WALL_COLUMNS:
                assert x[c] == y[c], c
    # the run did something worth comparing: stressed, moving, removed
    last = checkpoint.read(os.path.join(full, "checkpoints", checkpoint.filename(4)))
    assert np.abs(last["S"]).max() > 1e6 and np.abs(last["v"]).max() > 0.0


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_run_name_parses_back(inputs, runs, form):
    full, _, _, _ = runs(form)
    name = naming.run_name_from_dir(full)
    assert full == naming.run_dir(os.path.dirname(os.path.dirname(full)), name)    # <out>/<prepared>/<suffix>
    p = naming.parse_run_name(name)
    _, prep, _ = inputs
    assert p["prepared"] == os.path.basename(prep) and p["mode"] == "replay" and p["frames"] == (0, 4)
    assert p["form"] == form and p["dx_mm"] == DX_MM and p["seed"] == 1
    assert p["brackets"] == naming.BRACKET_DEFAULTS                       # every default: no suffix
    meta = record.read_meta(os.path.join(full, "meta.json"))
    assert meta["run_name"] == name


@pytest.mark.spheral
@pytest.mark.parametrize("form", FORMS)
def test_body_is_recentred_at_every_boundary(inputs, runs, form):
    """RZ: after each boundary the mass-weighted mean axial displacement from x_ref is zero to round-off (the shift
    is recorded). 3D: momentum-conserving, its mean displacement is round-off (1.3e-16 dx measured) until particles
    are removed, and is not touched then; a removal leaves the remaining particles with a mean displacement of their
    own (the removed nose had expanded forward: 1.6e-4 to 1.6e-3 dx measured 2026-10-10), which exceeds
    RECENTRE_3D_DX and is removed like RZ's."""
    from spheral_frag.runner import run
    full, _, _, _ = runs(form)
    dx = DX_MM*1e-3
    rows = intervals(full)
    for k in range(1, 5):
        ck = checkpoint.read(os.path.join(full, "checkpoints", checkpoint.filename(k)))
        m = ck["m"]
        d = (m[:, None]*(ck["x"] - ck["x_ref"])).sum(axis=0)/m.sum()
        applied = int(rows[k - 1]["recentre_applied"]) == 1
        if form == "rz" or applied:
            # positions ~0.05 m carry round-off ~1e-17 m; threshold 1e-12 dx = 1e-14 m
            assert abs(d[0]) <= 1e-12*dx and (form == "rz" or np.linalg.norm(d) <= 1e-12*dx), (k, d/dx)
        else:
            assert np.linalg.norm(d) <= run.RECENTRE_3D_DX*dx
    if form == "rz":
        assert all(int(r["recentre_applied"]) == 1 for r in rows)
        assert any(abs(float(r["recentre_shift_dx"])) > 1e-6 for r in rows)    # SolidSPHRZ's drift was removed
    else:
        assert int(rows[0]["n_removed"]) == 0 and int(rows[0]["recentre_applied"]) == 0
        assert float(rows[0]["recentre_shift_dx"]) <= run.RECENTRE_3D_DX


@pytest.mark.spheral
def test_the_internal_ratio_drives_the_repeat(inputs, runs):
    """A threshold between interval 0->1's internal and whole KE ratios passes (the internal one is the key); a
    threshold below both repeats the interval at 2 N_c and then fails (exit 1) at --max-crossings, with failed.json
    naming the internal ratio."""
    import json
    d, prep, bodies = inputs
    full, _, _, _ = runs("rz")
    r0 = intervals(full)[0]
    # the first attempt, at 10 crossings (the default run's monitor may have switched interval 0 to 20 since the
    # limit keys on the larger of the von Mises and pressure estimates); the runs below are unmonitored
    n, ratios = r0["attempts_ke_ratio"].split(";")[0].split(":")
    internal, whole = (float(x) for x in ratios.split("/"))
    assert int(n) == 10 and internal < whole
    common = ["replay", "--body", rel(bodies["rz"]), "--prepared", rel(prep), "--frames", "0-1", "--monitor-every",
              "0"]
    between = math.sqrt(internal*whole)
    ok = launch(*common, "--out", rel(os.path.join(d, "between")), "--ke-fraction", repr(between))
    run_dir = only_run(os.path.join(d, "between"))
    assert int(intervals(run_dir)[0]["repeats"]) == 0, ok.stdout[-2000:]
    bad = launch(*common, "--out", rel(os.path.join(d, "tiny")), "--ke-fraction", "1e-12", "--max-crossings", "20",
                 check=False)
    assert bad.returncode == 1 and "internal KE" in bad.stderr
    fail_dir = only_run(os.path.join(d, "tiny"))
    with open(os.path.join(fail_dir, "failed.json")) as fh:
        f = json.load(fh)
    assert f["keyed_on"] == "ke_internal_ratio" and [x["crossings"] for x in f["attempts"]] == [10, 20]
    assert all(1e-12 <= x["ke_internal_ratio"] < x["ke_ratio"] for x in f["attempts"])
    assert not os.path.isfile(os.path.join(fail_dir, "summary.json"))


# ---------------------------------------------------------------------------------------------------------------
# Command line and numpy helpers (no long runs)
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.spheral
def test_command_line_exit_codes(inputs):
    d, prep, bodies = inputs
    common = ["replay", "--body", rel(bodies["rz"]), "--prepared", rel(prep), "--out", rel(os.path.join(d, "bad"))]
    assert launch(*common, "--frames", "0-9", check=False).returncode == 2           # frame outside the body
    assert launch(*common, "--frames", "0-4", "--scheme", "crksph", check=False).returncode == 2   # Task 9
    assert launch(*common, "--frames", "0-4", "--strength", "rubber", check=False).returncode == 2
    assert launch("replay", "--body", rel(bodies["rz"]), check=False).returncode == 2               # missing args
    r = launch("replay", "--body", rel(os.path.join(d, "nowhere")), "--prepared", rel(prep), "--frames", "0-4",
               check=False)
    assert r.returncode == 2 and "missing" in r.stderr


def test_non_default_brackets_name_the_run():
    prep = naming.prepare_name(RECEDE_RUN_NAME, 0, 4)
    b = naming.canonical_brackets({"crossings": 20, "damping": 0.2, "strength": "constant:100e6",
                                   "rate_min": 1e-5, "sampler": "radial"})
    name = naming.run_name(prep, "replay", (0, 4), "rz", 2.2, b, 1)
    assert name.endswith("_crossings-20_damping-0.2_strength-constant-100000000_ratemin-1e-05_sampler-radial")
    assert naming.parse_run_name(name)["brackets"] == b
    assert naming.run_name(prep, "replay", (0, 4), "rz", 2.2, naming.canonical_brackets({}), 1).endswith("seed1")
    with pytest.raises(ValueError):
        naming.canonical_brackets({"strength": "Flight"})
    with pytest.raises(ValueError):
        naming.parse_run_name(name.replace("constant-100000000", "constant-1e8"))     # not canonical


def test_scale_factor_and_energies():
    from spheral_frag import material
    from spheral_frag.runner import run
    # s = (c_l dt / (N_c D))^2: 10 crossings of 0.1 m in 0.5 s at 6,000 m/s -> (3,000)^2
    assert run.scale_factor(6000.0, 0.5, 10, 0.1) == pytest.approx(9.0e6, rel=1e-15)
    table = material.MaterialTable.load(os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag", "material",
                                                  "material_AA7075_scheil.npz"))
    n = 4
    rho = table.rho_free(np.full(n, 400.0))
    st = dict(m=np.full(n, 2.0), v=np.zeros((n, 2)), u=table.enthalpy(np.full(n, 400.0)), rho=rho*1.001,
              S=np.zeros((n, 2, 2)), S_zz=np.full(n, 1e6))
    ke, el, ke_com = run.energies(st, table, 1.0)
    assert ke == 0.0 and ke_com == 0.0
    K = table.bulk_modulus(np.array(400.0))
    G = table.shear_modulus(np.array(400.0))
    P = K*0.001
    V = 2.0/(rho*1.001)
    assert el == pytest.approx(float(np.sum(V*(1e12/(4*G) + P*P/(2*K)))), rel=1e-12)
    assert run.ke_ratio(1.0, 0.0) == math.inf and run.ke_ratio(0.0, 0.0) == 0.0
    # the slowness test keys on the internal KE: a rigid drift alone passes, ringing fails
    assert run.slowness(10.0, 9.9, 100.0) == (True, pytest.approx(1e-3), 0.1)
    assert run.slowness(10.0, 0.0, 100.0)[0] is False
    # re-centring: RZ shifts the axis by minus the mean axial displacement; 3D leaves round-off alone
    m = np.array([1.0, 3.0])
    xr = np.array([[0.0, 0.01], [0.02, 0.02]])
    x, sh, ap = run.recentre("rz", m, xr + [0.004, 0.001], xr, 0.01)
    assert ap and sh == pytest.approx(-0.4) and abs((m*(x[:, 0] - xr[:, 0])).sum()) < 1e-18
    assert np.array_equal(x[:, 1], xr[:, 1] + 0.001)
    x3 = np.zeros((2, 3))
    assert run.recentre("3d", m, x3 + 1e-14, x3, 0.01)[2] is False
    assert run.recentre("3d", m, x3 + 1e-6, x3, 0.01)[2] is True
