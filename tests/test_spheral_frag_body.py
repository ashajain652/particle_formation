"""The body stage (Spheral M2 plan, Task 2): particles from a frame, their temperatures at every frame, check 3.

The normal-loop tests run on the committed synthetic sphere run (`tests/fixtures/spheral_frag/sphere_run`, prepared
here with the synthetic material table): the lattice count within the faceting bound of V/dx^3 and of 1/2 pi R^2/dx^2
(RZ); the mass exact to 1e-14; a uniform temperature's heat content exact to 1e-14; a linear nodal field reproduced at
every particle to round-off; the removal frames on the receding sphere (frame 2 loses ten nose elements) exactly
those of the points inside the deactivated elements, the elements' own centroids among them; the RZ azimuth spread of
a field linear in y; the files rebuilt byte-identical, a second build skipped, the exit codes; the names parsed back;
and the committed check-3 table (data/spheral/m2_check3.csv) against the tolerances below.

The `fe_flight` test rebuilds two check-3 cases on the real flight's prepared input (git-ignored) and requires the
committed rows back."""
import csv
import math
import os
from dataclasses import replace

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from helpers import REPO_ROOT
from spheral_frag import body, frames, geometry, material, naming
from spheral_frag.__main__ import main

pytest.importorskip("pyvista")

FIXTURE_RUN = os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag", "sphere_run", syn.SPHERE_RUN_NAME)
PREP_NAME = naming.prepare_name(syn.SPHERE_RUN_NAME, 0, 2)
CHECK3_CSV = os.path.join(REPO_ROOT, "data", "spheral", "m2_check3.csv")

# Check 3's tolerances (plan Task 2, Step 5): the worst case over k = 0, 48, 100, 182, 300, 600 of the MVP flight,
# measured 2026-10-09 by analysis/spheral_m2_check3.py into data/spheral/m2_check3.csv, rounded up to one significant
# figure (worst case in brackets). The heat-content difference does not shrink with dx: it is the mass weighting the
# plan prescribes, not sampling -- m_i = rho_free(T_i) V_i makes hot particles lighter than the finite element's
# constant-density mass (-0.56 % at k = 48), and a particle in a partly consumed element (phi < 1) carries the whole
# element's volume (+1.5 % at k = 182, with rho_free -0.46 % against it); measured on 3D 1.5 mm, plan report.
HEAT_RTOL = {("rz", 2.2): 2e-2,    # [1.053e-2, k 182]
             ("rz", 1.1): 1e-2,    # [9.845e-3, k 182]
             ("rz", 0.55): 1e-2,   # [9.097e-3, k 182]
             ("3d", 3.0): 2e-2,    # [1.057e-2, k 182]
             ("3d", 2.2): 2e-2,    # [1.084e-2, k 182]
             ("3d", 1.5): 2e-2}    # [1.086e-2, k 182]
MASS_RTOL = 1e-14                  # a priori, exact by c up to round-off [2.2e-16]
VOLUME_RTOL = 0.2                  # measured, reported not removed [0.1019, k 300: phi < 1 elements and expansion]
AZIMUTH_SPREAD_K = {2.2: 8.0, 1.1: 8.0, 0.55: 9.0}   # measured, RZ [7.619, 7.870, 8.895 K, all k 182]


# ------------------------------------------------------------------------------------------------------ fixtures
@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    tdir = tmp_path_factory.mktemp("material")
    syn.synthetic_material_table(str(tdir / "material_table.npz"))
    out = tmp_path_factory.mktemp("prepare")
    assert main(["prepare", "--fe-run", FIXTURE_RUN, "--material-table", str(tdir), "--outdir", str(out),
                 "--quiet"]) == 0
    d = os.path.join(str(out), PREP_NAME)
    return d, material.MaterialTable.load(d)


@pytest.fixture(scope="module")
def loaded(prepared):
    d, table = prepared
    ld = body.frame_loader(d)
    return {k: ld(k) for k in (0, 1, 2)}, table


def with_T(frame, T):
    return replace(frame, node={**frame.node, "T": np.ascontiguousarray(T, dtype=np.float64)})


def inradius(frame):
    """Smallest distance from the origin to a surface face's plane (the faceted sphere lies outside it)."""
    a, b, c = (frame.points[frame.faces[:, i]] for i in range(3))
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    return float(np.abs(np.einsum("ij,ij->i", n, a)).min())


# ------------------------------------------------------------------------------------------------------ lattice
@pytest.mark.parametrize("form,dx", [("3d", 5e-3), ("3d", 2.5e-3), ("rz", 5e-3), ("rz", 1.25e-3)])
def test_lattice_count_within_the_faceting_bound(loaded, form, dx):
    fr, table = loaded
    b = body.build_body(fr.get, [0], form, dx, table)
    R, r_in = syn.SPHERE_R, inradius(fr[0])
    assert 0.95 * R < r_in < R
    if form == "3d":
        s = math.sqrt(3.0) / 2.0 * dx         # each point's cell lies within s of it: count bounded by the balls
        lo, hi = 4.0 / 3.0 * math.pi * (r_in - s) ** 3 / dx ** 3, 4.0 / 3.0 * math.pi * (R + s) ** 3 / dx ** 3
        assert np.array_equal(b.x_ref, b.x_sample)
        on = np.round(b.x_ref / dx)
        assert np.allclose(b.x_ref, on * dx, rtol=0, atol=1e-15)           # a node at the origin
    else:
        s = math.sqrt(2.0) / 2.0 * dx
        lo, hi = 0.5 * math.pi * (r_in - s) ** 2 / dx ** 2, 0.5 * math.pi * (R + s) ** 2 / dx ** 2
        j = b.x_ref[:, 1] / dx - 0.5
        assert np.allclose(j, np.round(j), atol=1e-12) and b.x_ref[:, 1].min() == 0.5 * dx
        assert np.all(b.x_sample[:, 2] == 0.0) and np.array_equal(b.x_sample[:, 1], b.x_ref[:, 1])
    assert lo <= len(b.pid) <= hi, (len(b.pid), lo, hi)
    assert np.all(geometry.TetLocator.from_frame(fr[0]).contains(b.x_sample))
    assert np.array_equal(b.pid, np.arange(len(b.pid)))


def test_rz_azimuth_rotates_the_half_plane():
    _, xs0 = body.lattice("rz", 1e-3, [-0.01] * 3, [0.01] * 3, 0.0)
    xr, xs = body.lattice("rz", 1e-3, [-0.01] * 3, [0.01] * 3, 90.0)
    assert np.allclose(xs[:, 2], xr[:, 1], rtol=0, atol=1e-18) and np.allclose(xs[:, 1], 0.0, atol=1e-18)
    assert np.array_equal(xs0[:, 0], xs[:, 0])


# ------------------------------------------------------------------------------------------- masses and heat
@pytest.mark.parametrize("form", ["3d", "rz"])
def test_mass_is_exact_and_the_volume_is_reported(loaded, form):
    fr, table = loaded
    b = body.build_body(fr.get, [1], form, 4e-3, table)           # frame 1: 300-900 K, so rho_free varies
    M = body.fe_solid_mass(fr[1], table.rho_solid)
    assert abs(math.fsum(b.m) / M - 1.0) <= MASS_RTOL and abs(b.check["mass_rel_diff"]) <= MASS_RTOL
    assert np.array_equal(b.rho, table.rho_free(b.T0)) and np.allclose(b.m, b.rho * b.V, rtol=1e-15, atol=0)
    if form == "3d":
        assert np.all(b.V == b.V[0]) and math.isclose(b.V[0], b.c * 4e-3 ** 3, rel_tol=1e-15)
    else:
        assert np.allclose(b.V, b.c * 2 * math.pi * b.x_ref[:, 1] * 4e-3 ** 2, rtol=1e-15, atol=0)
    vt, vphi = body.fe_volumes(fr[1])
    assert b.check["volume_rel_diff"] == math.fsum(b.V) / vphi - 1.0
    rho_mean = math.fsum(b.m) / math.fsum(b.V)                   # the particles are expanded: more volume
    assert b.check["volume_rel_diff"] > 0 and rho_mean < table.rho_solid


@pytest.mark.parametrize("form", ["3d", "rz"])
def test_uniform_temperature_heat_content_is_exact(loaded, form):
    fr, table = loaded
    f = with_T(fr[1], np.full(len(fr[1].points), 650.0))
    b = body.build_body({1: f}.get, [1], form, 4e-3, table)
    hp = b.per_frame[0]["heat_particles_J"]
    hfe = body.fe_heat_content(f, table)
    assert hfe > 0 and abs(hp / hfe - 1.0) <= 1e-14 and abs(b.per_frame[0]["heat_rel_diff"]) <= 1e-14


@pytest.mark.parametrize("form", ["3d", "rz"])
def test_linear_nodal_field_is_reproduced_at_every_particle(loaded, form):
    fr, table = loaded
    g = np.array([2000.0, -1500.0, 800.0])          # K/m
    lin = lambda x: 600.0 + x @ g                                     # noqa: E731
    fs = {k: with_T(fr[k], lin(fr[k].points) + 10.0 * k) for k in (0, 1, 2)}
    b = body.build_body(fs.get, [0, 1, 2], form, 4e-3, table, azimuth_deg=30.0 if form == "rz" else 0.0)
    for j, k in enumerate((0, 1, 2)):
        live = (b.remove_k < 0) | (b.remove_k > k)
        want = lin(b.x_sample[live]) + 10.0 * k
        assert np.max(np.abs(b.temperature[live, j] - want)) <= 1e-12 * np.abs(want).max()


def test_rz_azimuth_spread_of_a_field_linear_in_y(loaded):
    fr, table = loaded
    G = 1000.0                                        # K/m along y: T(phi) - T(0) = G r (cos(phi) - 1)
    f = with_T(fr[1], 500.0 + G * fr[1].points[:, 1])
    b = body.build_body({1: f}.get, [1], "rz", 4e-3, table)
    row = b.per_frame[0]
    rot = body.rotate_about_axis(b.x_sample, 180.0)
    inside = geometry.TetLocator.from_frame(f).contains(rot)
    want = 2.0 * G * b.x_ref[inside, 1].max()
    assert math.isclose(row["azimuth_spread_K"], want, rel_tol=1e-9)
    assert row["n_outside_other_azimuths"] >= 0


# ------------------------------------------------------------------------------------------------------ removal
def dead_tets(fr):
    """Frame 0's tetrahedra (compact coordinates) that frame 2 no longer has (by finite-element node indices)."""
    def keys(f):
        t = np.sort(f.node["node_ids"][f.tets], axis=1)
        return {tuple(r) for r in t}
    gone = keys(fr[0]) - keys(fr[2])
    t0 = np.sort(fr[0].node["node_ids"][fr[0].tets], axis=1)
    return np.array([i for i, r in enumerate(t0) if tuple(r) in gone])


def barycentric_inside(points, tets, x, tol):
    """(n_x,) True where x lies in one of `tets` with every barycentric weight >= -tol (brute force)."""
    out = np.zeros(len(x), bool)
    for t in tets:
        p = points[t]
        lam = np.linalg.solve((p[1:] - p[0]).T, (x - p[0]).T).T
        lam = np.concatenate([1.0 - lam.sum(axis=1, keepdims=True), lam], axis=1)
        out |= (lam >= -tol).all(axis=1)
    return out


def test_receding_sphere_removal_is_exactly_the_deactivated_elements(loaded):
    fr, table = loaded
    dead = dead_tets(fr)
    assert len(dead) == syn.SPHERE_N_DEAD
    # (a) particles at every frame-0 element's centroid: removed at frame 2 exactly where the element died
    T, rk, _ = sample_points(fr, table, fr[0].points[fr[0].tets].mean(axis=1))
    want = np.full(len(fr[0].tets), -1)
    want[dead] = 2
    assert np.array_equal(rk, want)
    # (b) a lattice body: removed = inside a dead element (strictly inside -> removed; removed -> inside, with slack)
    b = body.build_body(fr.get, [0, 1, 2], "3d", 2e-3, table)
    removed = b.remove_k >= 0
    assert set(np.unique(b.remove_k[removed])) == {2} and removed.sum() > 10
    strict = barycentric_inside(fr[0].points, fr[0].tets[dead], b.x_sample, -1e-9)
    loose = barycentric_inside(fr[0].points, fr[0].tets[dead], b.x_sample, 1e-9)
    assert np.all(removed[strict]) and np.all(loose[removed])
    # held at the removal frame, NaN never before it (no NaN after it: the range ends at frame 2)
    assert np.array_equal(b.temperature[removed, 2], b.temperature[removed, 1])
    assert np.all(np.isfinite(b.temperature))
    assert b.per_frame[2]["n_removed"] == removed.sum() and b.per_frame[2]["n_present"] == (~removed).sum()


def sample_points(fr, table, x):
    """build_body's frame loop on arbitrary points: (temperature, remove_k, per-frame rows)."""
    b = body.build_body(fr.get, [0, 1, 2], "3d", 1.0, table, points=(x, x))
    assert len(b.pid) == len(x)
    return b.temperature, b.remove_k, b.per_frame


def test_a_removed_particle_stays_removed_and_reads_nan_after(loaded):
    fr, table = loaded
    seq = {0: fr[0], 1: replace(fr[2], k=1), 2: replace(fr[0], k=2), 3: replace(fr[2], k=3)}   # dies, revives
    b = body.build_body(seq.get, [0, 1, 2, 3], "3d", 2e-3, table)
    gone = b.remove_k >= 0
    assert gone.any() and set(np.unique(b.remove_k[gone])) == {1}
    assert np.array_equal(b.temperature[gone, 1], b.temperature[gone, 0])
    assert np.all(np.isnan(b.temperature[gone, 2:])) and np.all(np.isfinite(b.temperature[~gone]))
    assert b.per_frame[2]["n_present"] == (~gone).sum()


# ------------------------------------------------------------------------------------------------ the command
def run_body(prep_dir, out, *extra, frames_="0-2", form="3d", dx="5"):
    return main(["body", "--prepared", prep_dir, "--frames", frames_, "--form", form, "--dx-mm", dx,
                 "--out", str(out), "--quiet", *extra])


def test_command_writes_the_files_rebuilds_byte_identical_and_resumes(prepared, tmp_path, capsys):
    d, table = prepared
    assert run_body(d, tmp_path, form="rz", dx="2.5") == 0
    name = naming.body_name(PREP_NAME, 0, 2, "rz", 2.5)
    bd = tmp_path / name
    files = {f: (bd / f).read_bytes() for f in ("body.npz", "temperature.npy")}
    import json
    doc = json.loads((bd / "body.json").read_text())
    assert doc["name"] == name and doc["passed"] and doc["n_particles"] > 0 and doc["frames"]["k"] == [0, 1, 2]
    assert doc["temperature"]["order"] == "F" and len(doc["per_frame"]) == 3
    arrays, T = body.load_body(str(bd))
    n = len(arrays["pid"])
    assert T.shape == (n, 3) and T.dtype == np.float64 and T.flags.f_contiguous and isinstance(T, np.memmap)
    assert arrays["x_ref"].shape == (n, 2) and arrays["x_sample"].shape == (n, 3)
    assert set(arrays) >= {"pid", "x_ref", "V", "m", "rho", "T0", "remove_k", "frames_k", "frames_time_s"}
    assert np.array_equal(arrays["T0"], T[:, 0]) and list(arrays["frames_k"]) == [0, 1, 2]
    assert abs(math.fsum(arrays["m"]) / doc["check3"]["fe_solid_mass_k0_kg"] - 1.0) <= MASS_RTOL
    # the same body in memory
    b = body.build_body(body.frame_loader(d), [0, 1, 2], "rz", 2.5e-3, table)
    assert np.array_equal(b.temperature, np.asarray(T)) and np.array_equal(b.remove_k, arrays["remove_k"])
    del T
    capsys.readouterr()
    assert run_body(d, tmp_path, form="rz", dx="2.5") == 0
    assert "built already" in capsys.readouterr().out
    assert run_body(d, tmp_path, "--force", form="rz", dx="2.5") == 0
    assert {f: (bd / f).read_bytes() for f in files} == files
    (bd / "temperature.npy").write_bytes(b"x")              # a damaged file is rebuilt, not skipped
    capsys.readouterr()
    assert run_body(d, tmp_path, form="rz", dx="2.5") == 0
    assert "built already" not in capsys.readouterr().out and (bd / "temperature.npy").read_bytes() == \
        files["temperature.npy"]


def test_command_exit_codes(prepared, tmp_path):
    d, _ = prepared
    assert run_body(str(tmp_path / "nowhere"), tmp_path) == 2                       # no prepared input
    assert run_body(d, tmp_path, frames_="0-7") == 2                                 # no such frame
    assert run_body(d, tmp_path, frames_="2-0") == 2                                 # reversed
    assert run_body(d, tmp_path, dx="-1") == 2
    assert run_body(d, tmp_path, "--azimuth-deg", "90") == 2                         # azimuth is RZ only
    assert main(["body", "--prepared", d]) == 2                                      # missing arguments
    assert run_body(d, tmp_path, form="rz", dx="500") == 1                           # no lattice point inside
    import shutil
    broken = tmp_path / "broken" / PREP_NAME
    shutil.copytree(d, broken)
    (broken / "frames" / "frame_00001.npz").write_bytes(b"not a frame")
    assert run_body(str(broken), tmp_path / "out", frames_="0-1") == 1              # a frame that cannot be read
    assert not any(p.name == "body.json" for p in (tmp_path / "out").rglob("*"))
    assert run_body(d, tmp_path, form="rz", dx="5", frames_="1") == 0
    assert run_body(d, tmp_path, "--azimuth-deg", "45", form="rz", dx="5", frames_="1") == 0
    assert (tmp_path / naming.body_name(PREP_NAME, 1, 1, "rz", 5.0, 45.0)).is_dir()


# ------------------------------------------------------------------------------------------------------ names
def test_body_names_parse_back():
    name = naming.body_name(PREP_NAME, 0, 47, "rz", 2.2)
    assert name == "body_{}_k00000-00047_rz_dx2.20mm".format(PREP_NAME)
    p = naming.parse_body_name(name)
    assert p == {"prepared": PREP_NAME, "frames": (0, 47), "form": "rz", "dx_mm": 2.2, "azimuth_deg": 0.0,
                 "synthetic": None, "prepare": naming.parse_prepare_name(PREP_NAME)}
    for args in [(PREP_NAME, 3, 9, "3d", 1.5), (PREP_NAME, 0, 0, "rz", 0.55, 90.0),
                 (PREP_NAME, 0, 4, "rz", 1.1, -22.5, "radial"), ("other", 1, 2, "3d", 3.0, 0.0, "radial")]:
        n = naming.body_name(*args)
        q = naming.parse_body_name(n)
        assert naming.body_name(q["prepared"], *q["frames"], q["form"], q["dx_mm"], q["azimuth_deg"],
                                q["synthetic"]) == n
    for bad in ("body_x_k0-1_rz_dx2.20mm", name + "_azimuth-0deg", name.replace("rz", "2d"), "x" + name,
                name + "_synthetic-Radial"):
        with pytest.raises(ValueError):
            naming.parse_body_name(bad)
    with pytest.raises(ValueError):
        naming.body_name(PREP_NAME, 0, 1, "rz", 2.2, synthetic="a b")


# ----------------------------------------------------------------------------------------- the flight's check 3
def _check3_rows():
    with open(CHECK3_CSV, newline="") as fh:
        return list(csv.DictReader(fh))


def test_committed_check3_meets_the_tolerances():
    rows = _check3_rows()
    cases = {(r["form"], float(r["dx_mm"])) for r in rows}
    assert cases == set(HEAT_RTOL) and sorted({int(r["k"]) for r in rows}) == [0, 48, 100, 182, 300, 600]
    for r in rows:
        key = (r["form"], float(r["dx_mm"]))
        assert abs(float(r["mass_rel_diff"])) <= MASS_RTOL
        assert abs(float(r["volume_rel_diff"])) <= VOLUME_RTOL
        h = float(r["heat_rel_diff"])
        if int(r["k"]) == 0:              # a uniform 300 K: both heat contents zero, no relative number
            assert math.isnan(h) and abs(float(r["heat_abs_diff_J"])) < 1e-6
        else:
            assert abs(h) <= HEAT_RTOL[key], (r["k"], key, h)
        if r["form"] == "rz":
            assert float(r["azimuth_spread_K"]) <= AZIMUTH_SPREAD_K[key[1]]


@pytest.mark.fe_flight
def test_flight_check3_rows_rebuild(tmp_path):
    flight = ("model_d100.00mm_v07.50000kms_h077.500km_us76_sesam-table_none_fem-physics_melt-girin_material-"
              "AA7075_scheil_deeprunoff-off")
    prep = os.environ.get("SPHERAL_FRAG_PREPARED", os.path.join(REPO_ROOT, "spheral_output", "prepare",
                                                                 "prep_{}_k00000-01254".format(flight)))
    if not os.path.isfile(os.path.join(prep, "prepare.json")):
        pytest.skip("the flight's prepared input is not at {} (set SPHERAL_FRAG_PREPARED)".format(prep))
    import json
    rows = {(int(r["k"]), r["form"], float(r["dx_mm"])): r for r in _check3_rows()}
    for k, form, dx in ((182, "rz", 2.2), (182, "3d", 3.0), (48, "rz", 1.1)):
        assert run_body(prep, tmp_path, frames_="{0}-{0}".format(k), form=form, dx=str(dx)) == 0
        name = naming.body_name(os.path.basename(prep.rstrip(os.sep)), k, k, form, dx)
        doc = json.loads((tmp_path / name / "body.json").read_text())
        want = rows[(k, form, dx)]
        f0 = doc["per_frame"][0]
        assert doc["n_particles"] == int(want["n_particles"])
        for key, got in (("c", doc["c"]), ("mass_rel_diff", doc["check3"]["mass_rel_diff"]),
                         ("volume_rel_diff", doc["check3"]["volume_rel_diff"]),
                         ("heat_rel_diff", f0["heat_rel_diff"])):
            assert "{:.9g}".format(got) == want[key], (k, form, dx, key)
        assert abs(f0["heat_rel_diff"]) <= HEAT_RTOL[(form, dx)]
