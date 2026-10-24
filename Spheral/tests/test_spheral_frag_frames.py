"""Frame import and the compact prepared frame (Spheral M1 plan, Task 4; spec §12 check "Frame import").

On the committed synthetic runs (tests/fixtures/spheral_frag/, README): every array `frames.read_frame` returns is
the array pyvista reads directly, bitwise, and survives compact -> write -> load bitwise per original node; the
surfaces are closed, manifold and outward with divergence volume = tetrahedra volume; the mass sum phi rho V + film
+ deep matches the history's `mass_kg` to the CSV's 9 digits; and the checks detect what they exist to detect
(a removed triangle, inward-wound faces, non-manifold edges and vertices, a missing field, an unmatched time)."""
import json
import os
import shutil

from dataclasses import replace

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from helpers import REPO_ROOT
from spheral_frag import contract, frames
from spheral_frag.contract import ContractError, Frame

pytest.importorskip("pyvista")

FIXTURE_DIR = os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag")
FIXTURES = {"sphere_run": syn.SPHERE_RUN_NAME, "dumbbell_frame": syn.DUMBBELL_RUN_NAME,
            "slab_frame": syn.SLAB_RUN_NAME}

# Thresholds (each with its source):
# - CSV_DIGITS_RTOL: a value written "{:.9g}" is within half a unit of its 9th significant digit, i.e. within
#   5e-9 relative (the bound at a leading digit 1). The plan's "2e-9 (the CSV's 9 digits)" is not a bound: the
#   sphere's row 1 mass differs from the exact sum by 3.4e-9 (measured, below).
# - VOLUME_RTOL: plan Task 4, divergence volume against the tetrahedra's volume (a priori: round-off).
CSV_DIGITS_RTOL = 5e-9
VOLUME_RTOL = 1e-12


def run_of(fixture):
    return frames.read_fe_run(os.path.join(FIXTURE_DIR, fixture, FIXTURES[fixture]))


@pytest.fixture(scope="module")
def sphere():
    run = run_of("sphere_run")
    return run, [frames.read_frame(run, k) for k, _, _ in run.frames]


FRAME_CASES = [(fixture, k) for fixture in sorted(FIXTURES) for k in ((0, 1, 2) if fixture == "sphere_run" else (0,))]


def copy_run(tmp_path, fixture="sphere_run"):
    dst = tmp_path / fixture
    shutil.copytree(os.path.join(FIXTURE_DIR, fixture), dst)
    return dst, FIXTURES[fixture]


# ------------------------------------------------------------------------------------------------------- reading
def test_read_fe_run_accepts_run_dir_json_and_outdir():
    out = os.path.join(FIXTURE_DIR, "sphere_run")
    name = syn.SPHERE_RUN_NAME
    runs = [frames.read_fe_run(p) for p in (os.path.join(out, name), os.path.join(out, name + ".json"), out)]
    for run in runs:
        assert run.name == name and run.directory == os.path.join(out, name)
        assert run.frames == [(0, 0.0, 0), (1, 0.5, 1), (2, 1.0, 2)]
        assert run.history["time_s"].tolist() == [0.0, 0.5, 1.0]
        assert set(f.name for f in contract.fields("history")) <= set(run.history)
        assert "surface_T_max_K" in run.history                 # extra columns are kept, not required
        assert frames.run_value(run, "seed") == 12345 and frames.run_value(run, "material") == syn.SYNTHETIC_MATERIAL
        np.testing.assert_array_equal(frames.run_v_hat(run), [1.0, 0.0, 0.0])
    assert frames.history_row(runs[0], 1)["p_w_stag_step_Pa"] == 1715.32034


def test_read_fe_run_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        frames.read_fe_run(tmp_path / "nowhere")
    with pytest.raises(FileNotFoundError):
        frames.read_fe_run(tmp_path)                             # an empty directory holds no run
    dst, name = copy_run(tmp_path)
    # a history column of the contract missing -> ContractError naming it
    csv_path = dst / (name + ".csv")
    lines = csv_path.read_text().splitlines()
    cols = lines[0].split(",")
    i = cols.index("load_factor_g")
    csv_path.write_text("\n".join(",".join(c for j, c in enumerate(l.split(",")) if j != i) for l in lines) + "\n")
    with pytest.raises(ContractError, match="load_factor_g"):
        frames.read_fe_run(dst)
    shutil.rmtree(dst)
    dst, name = copy_run(tmp_path)
    doc = json.loads((dst / (name + ".json")).read_text())
    del doc["settings"]["v_hat_body"], doc["settings"]["seed"]
    (dst / (name + ".json")).write_text(json.dumps(doc))
    with pytest.raises(ContractError, match=r"v_hat.*|seed"):
        frames.read_fe_run(dst)


@pytest.mark.parametrize("fixture,k", FRAME_CASES)
def test_read_frame_equals_pyvista_bitwise(fixture, k):
    run = run_of(fixture)
    got = frames.read_frame(run, k)
    want, vtp_points = syn.read_frame_pyvista(run.directory, k)
    assert got.k == k and got.time_s == run.frame_time(k)
    for attr in ("points", "tets", "faces"):
        a, b = getattr(got, attr), getattr(want, attr)
        assert a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes(), attr
    assert vtp_points.tobytes() == got.points.tobytes()        # one node array
    for where in ("node", "tet", "patch"):
        g, w = getattr(got, where), getattr(want, where)
        assert set(g) == set(w)
        for key in w:
            assert g[key].dtype == w[key].dtype and g[key].tobytes() == w[key].tobytes(), (where, key)
    assert frames.missing(got) == [] and sorted(frames.absent_optional(got)) == sorted(contract.FIXTURE_OMITTED)


def test_frame_times_map_to_their_rows(tmp_path):
    dst, name = copy_run(tmp_path)
    vtk = dst / name / "vtk"
    for pvd in ("field.pvd", "surface.pvd"):                    # within 1e-6 s * max(1, t): still row 1
        p = vtk / pvd
        p.write_text(p.read_text().replace('timestep="0.5"', 'timestep="0.5000009"'))
    run = frames.read_fe_run(dst)
    assert run.frames[1] == (1, 0.5000009, 1)
    for pvd in ("field.pvd", "surface.pvd"):                    # 1.1e-6 s off: no row
        p = vtk / pvd
        p.write_text(p.read_text().replace('timestep="0.5000009"', 'timestep="0.5000011"'))
    with pytest.raises(ContractError, match="no history row"):
        frames.read_fe_run(dst)
    p = vtk / "surface.pvd"                                     # the two collections disagree
    p.write_text(p.read_text().replace('timestep="0.5000011"', 'timestep="0.5"'))
    with pytest.raises(ContractError, match="field.pvd time"):
        frames.read_fe_run(dst)
    assert frames.match_row(np.array([0.0, 0.5, 1.0]), 1.0 + 0.9e-6) == 2


def test_missing_required_field_is_reported(sphere, tmp_path):
    run, fr = sphere
    f = fr[1]
    patch = {k: v for k, v in f.patch.items() if k != "p_w"}
    frame = Frame(k=1, time_s=0.5, points=f.points, tets=f.tets, faces=f.faces, node=f.node, tet=f.tet, patch=patch)
    assert frames.missing(frame) == ["p_w"]
    # on disk: a run whose vtp lacks p_w and q_rad (optional) reads, and `missing` names exactly p_w
    lack = [Frame(k=g.k, time_s=g.time_s, points=g.points, tets=g.tets, faces=g.faces, node=g.node, tet=g.tet,
                  patch={k: v for k, v in g.patch.items() if k not in ("p_w", "q_rad")}) for g in fr]
    doc = json.loads(open(run.json_path).read())
    syn.write_fe_run(str(tmp_path), lack, run.history, doc)
    run2 = frames.read_fe_run(tmp_path / run.name)
    g = frames.read_frame(run2, 2)
    assert frames.missing(g) == ["p_w"] and sorted(frames.absent_optional(g)) == sorted(["q_rad", *contract.FIXTURE_OMITTED])
    assert frames.missing(Frame(k=0, time_s=0.0, points=f.points, tets=f.tets[:0], faces=f.faces)) == \
        ["tets"] + [x.key for loc in ("node", "tet", "patch") for x in contract.data_fields(loc) if x.required]


# ------------------------------------------------------------------------------------------------------- checks
@pytest.mark.parametrize("fixture,k", FRAME_CASES)
def test_surfaces_closed_manifold_outward(fixture, k):
    run = run_of(fixture)
    s = frames.surface_checks(frames.read_frame(run, k))
    assert s["n_open_directed_edges"] == 0 and s["n_boundary_edges"] == 0
    assert s["n_nonmanifold_edges"] == 0 and s["n_nonmanifold_vertices"] == 0
    assert s["faces_on_active_nodes"] and s["surface_matches_tets"] and s["n_faces_interior"] == 0
    assert s["rel_diff"] <= VOLUME_RTOL                          # measured 0.0 on every fixture frame
    assert s["n_inward_faces"] == 0 and s["outward"]


def test_staircase_frame_2(sphere):
    run, fr = sphere
    f = fr[2]
    s = frames.surface_checks(f)
    assert s["n_open_directed_edges"] == 0 and s["outward"] and s["surface_matches_tets"]
    assert (s["n_tets"], s["n_faces"], s["n_nodes"], s["n_nodes_active"]) == (2487, 1198, 753, 752)
    assert s["volume_div"] < frames.surface_checks(fr[1])["volume_div"]
    m = frames.mass_checks(f, syn.RHO, syn.RHO_LIQUID, frames.history_row(run, 2))
    assert m["n_partial_tets"] == 2 and m["phi_min"] == 0.25
    assert m["volume_phi"] < m["volume_tets"]                     # phi < 1: the tets' union is not the body
    vol = frames.tet_volumes(f.points, f.tets)
    assert m["solid_mass"] == pytest.approx(syn.RHO * (vol.sum() - vol[f.tet["phi"] < 1] @ (1 - f.tet["phi"][f.tet["phi"] < 1])),
                                            rel=1e-14)


@pytest.mark.parametrize("fixture,k", FRAME_CASES)
def test_mass_against_history(fixture, k):
    run = run_of(fixture)
    f = frames.read_frame(run, k)
    row = frames.history_row(run, run.frame_row(k))
    m = frames.mass_checks(f, syn.RHO, syn.RHO_LIQUID, row)
    # the CSV's 9 digits; measured sphere 1.39e-9, 3.42e-9, 1.51e-9; dumbbell 8.0e-10; slab 2.2e-16
    assert m["mass_rel_diff"] <= CSV_DIGITS_RTOL
    assert m["fe_mass"] == pytest.approx(syn.frame_mass(f), rel=1e-14)   # the fixture builder's sum (round-off)
    for part in ("film", "deep"):
        assert m[part + "_rel_diff"] <= CSV_DIGITS_RTOL


def test_one_triangle_removed_gives_three_open_edges(sphere):
    _, fr = sphere
    for f in (fr[0], fr[2]):
        g = Frame(k=f.k, time_s=f.time_s, points=f.points, tets=f.tets, faces=np.delete(f.faces, 17, axis=0))
        s = frames.surface_checks(g)
        assert s["n_open_directed_edges"] == 3 and s["n_boundary_edges"] == 3
        assert s["n_boundary_faces_missing"] == 1 and not s["surface_matches_tets"]
        assert s["n_nonmanifold_edges"] == 0


def test_inward_faces_detected_and_reoriented(sphere, tmp_path):
    """The real Step 3 export may wind some triangles inward after deaths (Task 2's finding: 4 of 1198 on a
    10-death frame). Four isolated faces wound inward: reported face by face, 6 open directed edges each, and
    orient_outward restores the outward winding exactly, flipping only those."""
    _, fr = sphere
    f = fr[2]
    edges = [set(map(frozenset, ((a, b), (b, c), (c, a)))) for a, b, c in f.faces]
    pick = []
    for i in range(0, len(f.faces), 97):                          # mutually non-adjacent faces
        if all(not (edges[i] & edges[j]) for j in pick):
            pick.append(i)
        if len(pick) == 4:
            break
    bad = f.faces.copy()
    bad[pick] = bad[pick][:, [0, 2, 1]]
    g = Frame(k=2, time_s=1.0, points=f.points, tets=f.tets, faces=bad, node=f.node, tet=f.tet, patch=f.patch)
    s = frames.surface_checks(g)
    assert s["n_inward_faces"] == 4 and s["inward_faces"] == sorted(pick) and not s["outward"]
    assert s["n_open_directed_edges"] == 24 and s["surface_matches_tets"]
    o, flipped = frames.orient_outward(g)
    assert flipped.tolist() == sorted(pick)
    assert o.faces.tobytes() == f.faces.tobytes() and o.patch is g.patch
    assert frames.surface_checks(o)["outward"]
    same, none = frames.orient_outward(f)
    assert same is f and len(none) == 0
    with pytest.raises(ContractError, match="orient"):
        cg, ids = frames.compact(g)
        frames.write_prepared_frame(str(tmp_path / "x.npz"), cg, ids, frames.derived_patch_arrays(cg),
                                    frames.mesh_from_frame(g, 0))
    h = Frame(k=2, time_s=1.0, points=f.points, tets=f.tets[1:], faces=f.faces)
    with pytest.raises(ContractError, match="cannot orient"):
        frames.orient_outward(h)


def _two_tets(shared):
    """Two unit-ish tetrahedra sharing one vertex (shared=1) or one edge (shared=2), and their boundary."""
    p = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    if shared == 1:
        points = np.vstack([p, -1.3 * p[1:]])                               # tet B = (0, 4, 5, 6): touches A at node 0
        tets = np.array([[0, 1, 2, 3], [0, 4, 6, 5]])
    else:
        r = np.array([[0.5, -1.0, -1.0], [-1.0, -1.0, 0.3]])
        points = np.vstack([p, r])                                # tet B = (0, 1, 4, 5): touches A along edge 0-1
        tets = np.array([[0, 1, 2, 3], [0, 1, 4, 5]])
    faces, _ = syn.outward_surface(points, tets)
    return Frame(k=0, time_s=0.0, points=points, tets=tets, faces=faces)


def test_nonmanifold_vertex_and_edge_counted():
    bowtie = frames.surface_checks(_two_tets(1))
    assert bowtie["n_nonmanifold_vertices"] == 1 and bowtie["nonmanifold_vertices"] == [0]
    assert bowtie["n_nonmanifold_edges"] == 0 and bowtie["n_open_directed_edges"] == 0 and bowtie["outward"]
    hinge = frames.surface_checks(_two_tets(2))
    assert hinge["n_nonmanifold_edges"] == 1 and hinge["n_nonmanifold_vertices"] == 0
    assert hinge["n_open_directed_edges"] == 0 and hinge["outward"]
    for s in (bowtie, hinge):
        assert s["rel_diff"] <= VOLUME_RTOL


class LinearTable:
    """Stands in for material.MaterialTable (Task 3): the fixtures' synthetic linear f_l(T)."""
    def liquid_fraction(self, T):
        return syn.synthetic_liquid_fraction(T)


def test_consistency_checks(sphere):
    _, fr = sphere
    for f in fr:
        c = frames.consistency_checks(f, LinearTable())
        assert c["max_abs_dfl"] == 0.0 and c["fl_bitwise"] and c["n_delta_m_mismatch"] == 0
        assert c["nan_forbidden"] == {}
    assert [frames.consistency_checks(f, None)["n_girin"] for f in fr] == [0, 306, 304]
    f = fr[1]
    patch = dict(f.patch, delta_m=f.patch["delta_m"].copy(), p_w=f.patch["p_w"].copy())
    i = int(np.flatnonzero(f.patch["closure"] == 0)[0])
    patch["delta_m"][i] = np.nan
    patch["p_w"][:3] = np.nan
    node = dict(f.node, f_l=np.nextafter(f.node["f_l"], 2.0))
    g = Frame(k=1, time_s=0.5, points=f.points, tets=f.tets, faces=f.faces, node=node, tet=f.tet, patch=patch)
    c = frames.consistency_checks(g, LinearTable())
    assert c["n_delta_m_mismatch"] == 1 and c["nan_forbidden"] == {"p_w": 3}
    assert not c["fl_bitwise"] and 0 < c["max_abs_dfl"] < 1e-15


# ------------------------------------------------------------------------------------------------------- compact
def test_compact_drops_unreferenced_nodes(sphere):
    _, fr = sphere
    c, ids = frames.compact(fr[2])
    assert len(ids) == 752 and len(c.points) == 752 and (np.diff(ids) > 0).all()
    assert c.node["T"].tobytes() == fr[2].node["T"][ids].tobytes()
    assert (ids[c.tets] == fr[2].tets).all() and (ids[c.faces] == fr[2].faces).all()
    s0, s1 = frames.surface_checks(fr[2]), frames.surface_checks(c)
    assert s1["volume_div"] == s0["volume_div"] and s1["outward"] and s1["n_nodes"] == s1["n_nodes_active"]
    c0, ids0 = frames.compact(fr[0])
    assert ids0.tolist() == list(range(753))


def test_derived_patch_arrays(sphere):
    _, fr = sphere
    f = fr[1]
    d = frames.derived_patch_arrays(f)
    area, normal, centroid, theta = syn.patch_geometry(f.points, f.faces)
    np.testing.assert_allclose(d["area"], area, rtol=1e-15)
    np.testing.assert_allclose(d["normal"], normal, rtol=0, atol=1e-15)
    np.testing.assert_array_equal(d["centroid"], centroid)
    np.testing.assert_allclose(d["theta"], theta, rtol=0, atol=1e-12)
    # closed surface: sum A n = 0 to round-off; theta in [0, pi]; the nose patches near 0
    assert np.abs((d["area"][:, None] * d["normal"]).sum(axis=0)).max() <= 1e-15 * d["area"].sum()
    assert d["theta"].min() < np.radians(10) and d["theta"].max() > np.radians(170)


def test_prepared_frame_round_trip_bitwise(sphere, tmp_path):
    """compact -> write (reduced, on mesh version 0 from frame 0) -> load equals the frame read from the files,
    bitwise, per original node via node_ids, less the dropped informational fields; the patch geometry recomputed on
    loading is prepare's bitwise; a rewrite is byte-identical; the files are numeric npz without pickles and follow
    the stored schema; frame 2 (ten nose elements dead, two at phi < 1) stores only its mask and its two phis."""
    _, fr = sphere
    mesh = frames.mesh_from_frame(fr[0], 0)
    frames.write_mesh(str(tmp_path), mesh)
    for f in fr:
        index = mesh.fits(f)
        assert index is not None and np.all(np.diff(index) > 0)        # the export keeps the mesh's order
        c, ids = frames.compact(f)
        d = frames.derived_patch_arrays(c)
        path = tmp_path / frames.prepared_frame_name(f.k)
        frames.write_prepared_frame(str(tmp_path), c, ids, d, mesh)
        assert path.exists()
        g = frames.load_prepared_frame(str(path))
        assert g.k == f.k and g.time_s == f.time_s
        ids_back = g.node["node_ids"]
        assert ids_back.tobytes() == ids.tobytes()
        assert g.points.tobytes() == f.points[ids_back].tobytes()
        assert (ids_back[g.tets] == f.tets).all() and (ids_back[g.faces] == f.faces).all()
        for key, v in f.node.items():
            assert g.node[key].tobytes() == v[ids_back].tobytes(), key
        for key, v in f.tet.items():
            assert g.tet[key].tobytes() == v.tobytes(), key
        for key, v in f.patch.items():
            if key in contract.PREPARED_DROPPED:
                assert key not in g.patch, key
            else:
                assert g.patch[key].tobytes() == v.tobytes(), key
        for key, v in d.items():
            assert g.patch[key].tobytes() == v.tobytes(), key
        # the schemas: the expanded arrays as contract.PREPARED_FRAME_ARRAYS says, the file as PREPARED_FRAME_STORED
        full = frames.load_prepared_arrays(str(path))
        sizes = {"n": len(ids), "m": len(f.tets), "f": len(f.faces)}
        for name, a in full.items():
            dtype, shape, _ = contract.PREPARED_FRAME_ARRAYS[name]
            assert a.dtype == np.dtype(dtype) and a.dtype.kind != "O", name
            assert a.shape == tuple(sizes[s] if isinstance(s, str) else s for s in shape), name
        assert set(contract.PREPARED_FRAME_ARRAYS) - set(full) == {
            "patch_thickness", "patch_slurry_depth", "patch_liquid_depth",
            *("patch_" + key for key in contract.FIXTURE_OMITTED if key not in contract.PREPARED_DROPPED)}
        raw = frames.load_stored_arrays(str(path))
        assert not set(raw) & set(contract.PREPARED_FRAME_RECOMPUTED)
        for name, a in raw.items():
            assert a.dtype == np.dtype(contract.PREPARED_FRAME_STORED[name][0]), name
        assert len(raw["tet_alive"]) == (len(mesh.tets) + 7) // 8 and len(raw["moved_node"]) == 0
        assert len(raw["phi_index"]) == int(np.sum(f.tet["phi"] != 1.0)) == (2 if f.k == 2 else 0)
        # rebuilt: byte-identical file
        first = path.read_bytes()
        frames.write_prepared_frame(str(path), c, ids, frames.derived_patch_arrays(c), mesh)
        assert path.read_bytes() == first


def test_moved_nodes_and_a_new_mesh_version(sphere, tmp_path):
    """A receding surface (three nodes moved) is stored as those nodes alone and loads bitwise; a tetrahedron the
    version does not hold, or another node count, is not the version's (prepare starts a new one); a frame whose
    tetrahedra are out of the version's order is put in it, fields with them, and is refused by the writer if not."""
    _, fr = sphere
    mesh = frames.mesh_from_frame(fr[0], 0)
    frames.write_mesh(str(tmp_path), mesh)
    f = fr[1]
    pts = f.points.copy()
    moved = np.unique(f.faces[:3].ravel())[:3]
    pts[moved] *= 0.999
    g = replace(f, points=pts)
    c, ids = frames.compact(g)
    frames.write_prepared_frame(str(tmp_path), c, ids, frames.derived_patch_arrays(c), mesh)
    path = str(tmp_path / frames.prepared_frame_name(1))
    assert sorted(frames.load_stored_arrays(path)["moved_node"]) == sorted(moved)
    back = frames.load_prepared_frame(path)
    assert back.points.tobytes() == pts[back.node["node_ids"]].tobytes()
    assert back.patch["area"].tobytes() == frames.derived_patch_arrays(c)["area"].tobytes()
    # re-gridding
    other = replace(f, tets=np.vstack([f.tets[1:], f.tets[:1, [1, 0, 2, 3]]]))
    assert mesh.fits(other) is None
    assert mesh.fits(replace(f, points=np.vstack([f.points, f.points[:1]]))) is None
    # order: reversed tetrahedra are put back in the version's order, fields with them
    rev = replace(f, tets=f.tets[::-1], tet={key: v[::-1] for key, v in f.tet.items()})
    put, index = frames.in_mesh_order(rev, mesh.fits(rev))
    assert np.array_equal(put.tets, f.tets) and np.all(np.diff(index) > 0)
    assert all(np.array_equal(put.tet[key], f.tet[key]) for key in f.tet)
    c, ids = frames.compact(rev)
    with pytest.raises(ValueError, match="order"):
        frames.write_prepared_frame(str(tmp_path / "x.npz"), c, ids, frames.derived_patch_arrays(c), mesh)


def test_prepared_frame_schema_refusals(sphere, tmp_path):
    _, fr = sphere
    c, ids = frames.compact(fr[1])
    d = frames.derived_patch_arrays(c)
    with pytest.raises(ValueError, match="not in contract"):
        frames.prepared_arrays(c, ids, dict(d, wobble=d["area"]))
    with pytest.raises(ValueError, match="shape"):
        frames.prepared_arrays(c, ids, dict(d, thickness=d["area"][:-1]))
    with pytest.raises(ValueError, match="lacks required"):
        frames.prepared_arrays(c, ids, {k: v for k, v in d.items() if k != "theta"})
    a = frames.prepared_arrays(c, ids, dict(d, thickness=np.full(len(c.faces), 0.1)))
    assert a["patch_thickness"].dtype == np.float64 and a["tets"].dtype == np.int32
