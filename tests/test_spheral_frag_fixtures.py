"""The synthetic finite-element runs of tests/fixtures/spheral_frag/ (Spheral M1 plan, Task 2).

They carry every item of the frame contract under its finite-element name (`contract.FE_FIELDS`), they are what
`spheral_frag_synthetic` builds (the regenerated sphere run equals the committed one bitwise), and their analytic
content holds, so Tasks 4-9 can test the core against the answers in the fixtures' README."""
import csv
import json
import math
import os
import re

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from helpers import REPO_ROOT
from spheral_frag import contract

pv = pytest.importorskip("pyvista")

FIXTURE_DIR = os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag")
FIXTURES = {"sphere_run": syn.SPHERE_RUN_NAME, "dumbbell_frame": syn.DUMBBELL_RUN_NAME,
            "slab_frame": syn.SLAB_RUN_NAME}
CSV_DIGITS_RTOL = 5e-9   # a value through "{:.9g}" is within half a unit of its 9th digit: <= 5e-9 relative
CSV_RTOL = 1e-8          # a product of two such values


def run_paths(fixture):
    out = os.path.join(FIXTURE_DIR, fixture)
    name = FIXTURES[fixture]
    return out, name, os.path.join(out, name)


def read_history(fixture):
    out, name, _ = run_paths(fixture)
    with open(os.path.join(out, name + ".csv")) as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def read_json(fixture):
    out, name, _ = run_paths(fixture)
    with open(os.path.join(out, name + ".json")) as fh:
        return json.load(fh)


def json_path(doc, dotted):
    for part in dotted.split("."):
        doc = doc[part]
    return doc


def pvd_entries(run_dir, kind):
    with open(os.path.join(run_dir, "vtk", kind + ".pvd")) as fh:
        return re.findall(r'<DataSet timestep="([^"]+)" file="([^"]+)"/>', fh.read())


def committed_frames(fixture):
    _, _, run_dir = run_paths(fixture)
    return [syn.read_frame_pyvista(run_dir, k)[0] for k in range(len(pvd_entries(run_dir, "field")))]


def directed_edges_without_reverse(faces):
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    have = set(map(tuple, e))
    return sum((b, a) not in have for a, b in have) + (len(e) - len(have))


def divergence_volume(points, faces):
    a, b, c = (points[faces[:, i]] for i in range(3))
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


@pytest.mark.parametrize("fixture", sorted(FIXTURES))
def test_fixture_carries_every_contract_item_under_its_fe_name(fixture):
    out, name, run_dir = run_paths(fixture)
    history, doc = read_history(fixture), read_json(fixture)
    field_entries, surface_entries = pvd_entries(run_dir, "field"), pvd_entries(run_dir, "surface")
    assert [t for t, _ in field_entries] == [t for t, _ in surface_entries]
    assert [f for _, f in field_entries] == ["field_{}.vtu".format(k) for k in range(len(field_entries))]
    for k, (t, _) in enumerate(field_entries):
        grid = pv.read(os.path.join(run_dir, "vtk", "field_{}.vtu".format(k)))
        poly = pv.read(os.path.join(run_dir, "vtk", "surface_{}.vtp".format(k)))
        for f in contract.data_fields("node"):
            assert f.name in grid.point_data, (fixture, k, f.key)
        for f in contract.data_fields("tet"):
            assert f.name in grid.cell_data, (fixture, k, f.key)
        for f in contract.data_fields("patch"):                    # required and optional: the export writes all
            assert f.name in poly.cell_data, (fixture, k, f.key)
        assert set(grid.point_data) | set(grid.cell_data) == {f.name for f in contract.data_fields("node")
                                                              + contract.data_fields("tet")}
        assert set(poly.cell_data) == {f.name for f in contract.data_fields("patch")}
        assert set(poly.point_data) == {contract.fe_name("T")}
        # one node array, every mesh node in both files; vtu cells all tetrahedra, vtp cells all triangles
        assert np.asarray(grid.points).dtype == np.float64
        assert np.asarray(grid.points).tobytes() == np.asarray(poly.points).tobytes()
        assert set(np.unique(grid.celltypes)) == {pv.CellType.TETRA}
        assert (np.asarray(poly.faces).reshape(-1, 4)[:, 0] == 3).all()
        # the frame time is the history row's within the contract's tolerance (the pvd prints 6 digits)
        row = int(np.argmin(np.abs(history[contract.fe_name("time_s")] - float(t))))
        assert abs(history["time_s"][row] - float(t)) <= contract.FRAME_TIME_RTOL * max(1.0, float(t))
        assert t == "{:.6g}".format(history["time_s"][row])
    for f in contract.fields("history"):
        assert f.name in history, f.key
    for f in contract.fields("run"):
        json_path(doc, f.name)                                      # KeyError if absent
    assert json_path(doc, contract.fe_name("run_name")) == name
    assert json_path(doc, contract.fe_name("frames_every")) == 1
    assert len(history["time_s"]) == len(field_entries)
    assert json_path(doc, contract.fe_name("v_hat")) == list(contract.V_HAT)


def test_fixtures_stay_below_the_size_budget():
    total = sum(os.path.getsize(os.path.join(d, f)) for d, _, files in os.walk(FIXTURE_DIR) for f in files)
    assert total < 3 * 1024 * 1024          # M1 plan, Global constraints


def test_regenerated_sphere_run_equals_committed_bitwise(tmp_path):
    pytest.importorskip("gmsh")
    frames, history, doc, _ = syn.sphere_run_parts()
    for built, read in zip(frames, committed_frames("sphere_run"), strict=True):
        for attr in ("points", "tets", "faces"):
            a, b = getattr(built, attr), getattr(read, attr).astype(getattr(built, attr).dtype)
            assert a.shape == b.shape and a.tobytes() == b.tobytes(), (built.k, attr)
        for where in ("node", "tet", "patch"):
            got, want = getattr(read, where), getattr(built, where)
            assert set(got) == set(want), (built.k, where)
            for key in want:
                assert got[key].dtype == np.float64 and got[key].tobytes() == want[key].tobytes(), (built.k, key)
    # the text files are byte-identical when rewritten
    run_dir = syn.write_fe_run(str(tmp_path), frames, history, doc)
    out, name, committed = run_paths("sphere_run")
    for rel in (name + ".csv", name + ".json"):
        with open(os.path.join(out, rel), "rb") as a, open(os.path.join(str(tmp_path), rel), "rb") as b:
            assert a.read() == b.read(), rel
    for rel in ("field.pvd", "surface.pvd"):
        with open(os.path.join(committed, "vtk", rel), "rb") as a, open(os.path.join(run_dir, "vtk", rel), "rb") as b:
            assert a.read() == b.read(), rel


def test_sphere_run_analytic_content():
    frames, history, doc = committed_frames("sphere_run"), read_history("sphere_run"), read_json("sphere_run")
    g = [syn.frame_geometry(f) for f in frames]
    # every frame: closed as an oriented surface, outward (divergence volume = tetrahedra volume), all nodes kept
    for f in frames:
        assert directed_edges_without_reverse(f.faces) == 0
        vt = syn.tet_volumes(f.points, f.tets).sum()
        assert divergence_volume(f.points, f.faces) == pytest.approx(vt, rel=1e-12)
        assert len(f.points) == len(frames[0].points)
        assert np.array_equal(f.node["f_l"], syn.synthetic_liquid_fraction(f.node["T"]))
    # frame 0: written before any step, nothing evaluated
    f0 = frames[0]
    assert (f0.patch["p_w"] == 0).all() and (f0.patch["tau"] == 0).all() and (f0.patch["closure"] == 1).all()
    assert np.isnan(f0.patch["delta_m"]).all() and (f0.patch["film_thickness"] == 0).all()
    assert (f0.patch["release_rate"] == 0).all() and (f0.node["T"] == 300.0).all()
    # frames 1-2: the analytic loads in each patch's own inclination; p_stag is the CSV's value exactly
    for k in (1, 2):
        f, (area, normal, centroid, theta) = frames[k], g[k]
        p_stag = history["p_w_stag_Pa"][k]
        p, tau = syn.newtonian_loads(theta, p_stag)
        assert np.array_equal(f.patch["p_w"], p) and np.array_equal(f.patch["tau"], tau)
        assert f.patch["p_w"].max() <= p_stag and (f.patch["p_w"][np.degrees(theta) >= syn.THETA_EVAL_DEG] == 0).all()
        # delta_m finite exactly where the closure is Girin's (0)
        assert np.array_equal(np.isfinite(f.patch["delta_m"]), f.patch["closure"] == 0.0)
        # the history's drag is the frame's own: mass_kg * load_factor_g * g0 = sum A (p cos + tau sin)
        D_hist = history["mass_kg"][k] * history["load_factor_g"][k] * syn.G0
        assert D_hist == pytest.approx(syn.drag(theta, area, f.patch["p_w"], f.patch["tau"]), rel=CSV_RTOL)
    # mass: sum phi rho V + film + deep against the CSV's 9 digits; the melt accounts and the sprayed increments
    sprayed = 0.0
    for k, f in enumerate(frames):
        solid, film, deep = syn.frame_masses(f)
        assert history["mass_kg"][k] == pytest.approx(solid + film + deep, rel=CSV_DIGITS_RTOL)
        assert history["film_mass_kg"][k] == pytest.approx(film, rel=CSV_DIGITS_RTOL, abs=0.0)
        assert history["deep_mass_kg"][k] == 0.0 and deep == 0.0
        sprayed += math.fsum(f.patch["release_rate"] * g[k][0])
        assert history["sprayed_mass_kg"][k] == pytest.approx(sprayed, rel=CSV_DIGITS_RTOL, abs=0.0)
    assert doc["results"]["initial_mass_kg"] == syn.frame_mass(frames[0])
    # frame 2: ten elements dead (the staircase), two partly consumed
    assert len(frames[1].tets) - len(frames[2].tets) == syn.SPHERE_N_DEAD
    assert sorted(frames[2].tet["phi"][frames[2].tet["phi"] < 1.0]) == sorted(syn.SPHERE_PHI)
    dead = doc["synthetic"]["dead_elements_frame2"]
    assert len(set(dead)) == syn.SPHERE_N_DEAD
    full = {tuple(t) for t in frames[1].tets}
    alive = {tuple(t) for t in frames[2].tets}
    assert alive < full and full - alive == {tuple(frames[1].tets[i]) for i in dead}


def test_slab_columns_have_their_exact_slurry_depths():
    (f,), history = committed_frames("slab_frame"), read_history("slab_frame")
    H, h = syn.SLAB["H"], syn.SLAB["h"]
    col = syn.slab_columns(f.points, syn.SLAB["L"], h)
    depth = H - f.points[:, 2]
    want = 0.5 + (np.asarray(syn.SLAB_SLURRY_DEPTH)[col] - depth) / syn.SLAB_GRADIENT
    np.testing.assert_allclose(f.node["f_l"], want, rtol=0, atol=1e-15)       # linear in depth, inside (0, 1)
    assert f.node["f_l"].min() > 0.0 and f.node["f_l"].max() < 1.0
    top, pcol, exact = syn.slab_exact_patches(f, h)
    assert np.bincount(pcol[exact]).tolist() == [24, 24, 24, 24, 32]
    for a in syn.slab_answers():
        sel = exact & (pcol == a["column"])
        assert (f.patch["film_thickness"][sel] == syn.SLAB_FILM).all()
        assert (f.patch["deep_thickness"][sel] == syn.SLAB_DEEP).all()
        if math.isnan(a["delta_m_m"]):
            assert np.isnan(f.patch["delta_m"][sel]).all() and (f.patch["closure"][sel] == 1).all()
        else:
            assert (f.patch["delta_m"][sel] == a["delta_m_m"]).all() and (f.patch["closure"][sel] == 0).all()
        # P1 along the vertical below an exact patch: f_l = 0.5 at the column's slurry depth (to round-off)
        nodes = np.flatnonzero(col == a["column"])
        z, fl = f.points[nodes, 2], f.node["f_l"][nodes]
        slope, intercept = np.polyfit(H - z, fl, 1)
        assert (0.5 - intercept) / slope == pytest.approx(a["slurry_depth_m"], abs=1e-15)
    assert not top[~np.isfinite(f.patch["delta_m"]) & top & (pcol < 4)].any()
    assert (f.patch["film_thickness"][~top] == 0).all() and np.isnan(f.patch["delta_m"][~top]).all()
    assert directed_edges_without_reverse(f.faces) == 0
    solid, film, deep = syn.frame_masses(f)
    assert history["mass_kg"][0] == pytest.approx(solid + film + deep, rel=CSV_DIGITS_RTOL)
    assert solid == pytest.approx(syn.RHO * syn.SLAB["L"] * syn.SLAB["W"] * H, rel=1e-12)


def test_dumbbell_neck_and_spheres():
    (f,) = committed_frames("dumbbell_frame")
    d = syn.DUMBBELL
    neck = syn.dumbbell_neck_patches(f, d["R"], d["r_neck"], d["L_neck"])
    spheres = syn.dumbbell_sphere_patches(f, d["R"], d["r_neck"], d["L_neck"])
    assert neck.sum() > 400 and spheres.sum() > 1000 and not (neck & spheres).any()
    # the neck's vertices lie on the cylinder r = r_neck (gmsh places them on the CAD surface)
    v = f.points[f.faces[neck]].reshape(-1, 3)
    np.testing.assert_allclose(np.linalg.norm(v[:, 1:], axis=1), d["r_neck"], rtol=0, atol=1e-9)
    assert directed_edges_without_reverse(f.faces) == 0
    vt = syn.tet_volumes(f.points, f.tets).sum()
    assert divergence_volume(f.points, f.faces) == pytest.approx(vt, rel=1e-12)
    assert vt == pytest.approx(syn.dumbbell_volume(d["R"], d["r_neck"], d["L_neck"]), rel=0.03)   # faceting
    assert (f.patch["film_thickness"] == 0).all() and np.isnan(f.patch["delta_m"]).all()
