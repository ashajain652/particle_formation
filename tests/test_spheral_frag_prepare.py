"""`python -m spheral_frag prepare` (Spheral M1 plan, Task 9, Steps 1 and 4) on the committed synthetic sphere run.

prepare exits 0 and writes every file, its JSON validates against the schema below, a second run is skipped, --force
rebuilds byte-identical frame, load and flight files, an interrupted prepare rebuilds only the frames it lacks, a run
without `p_w` exits 1 naming it, bad arguments and a missing directory exit 2, and the layer depths follow decision 9:
along the derived surface's normals (`n_derived`) when the frame carries them, with the facet-normal depths kept as
a diagnostic, and along the facet normals otherwise (`depth_normals: "facet"`).

The fixture's material (`synthetic_linear_fl`) is no finite-element material, so the runs pass the synthetic material
table (`spheral_frag_synthetic.synthetic_material_table`) with --material-table; one test checks that without it the
missing material exits 2."""
import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

import spheral_frag_synthetic as syn
from helpers import REPO_ROOT
from spheral_frag import contract, frames, loads, material, naming, record
from spheral_frag.__main__ import main

pytest.importorskip("pyvista")
jsonschema = pytest.importorskip("jsonschema")

FIXTURE = os.path.join(REPO_ROOT, "tests", "fixtures", "spheral_frag", "sphere_run")
RUN_DIR = os.path.join(FIXTURE, syn.SPHERE_RUN_NAME)
PREP_NAME = naming.prepare_name(syn.SPHERE_RUN_NAME, 0, 2)

_num = {"type": ["number", "null"]}
_int = {"type": "integer"}
FRAME_SCHEMA = {
    "type": "object",
    "required": ["k", "time_s", "history_row", "file", "sha256", "n_nodes", "n_tets", "n_faces", "volume_m3",
                 "fe_mass_kg", "mass_rel_diff", "fe_heat_content_J", "surface", "loads_valid", "drag", "thickness",
                 "zones", "depths", "timing_s"],
    "properties": {
        "k": _int, "time_s": {"type": "number"}, "history_row": _int, "file": {"type": "string"},
        "sha256": {"type": "string", "pattern": "^[0-9a-f]{64}$"}, "n_nodes": _int, "n_tets": _int, "n_faces": _int,
        "volume_m3": {"type": "number"}, "fe_mass_kg": {"type": "number"}, "mass_rel_diff": {"type": "number"},
        "fe_heat_content_J": {"type": "number"},
        "surface": {"type": "object", "required": ["n_open_directed_edges", "n_nonmanifold_edges",
                                                   "n_nonmanifold_vertices", "volume_rel_diff",
                                                   "n_inward_faces_as_written", "n_flipped"]},
        "loads_valid": {"type": "boolean"},
        "drag": {"type": "object", "required": ["smooth_N", "normals", "rel_smooth_hist", "patch_N", "table_N", "history_N", "lee_share"],
                 "properties": {"patch_N": {"type": "number"}, "table_N": _num, "history_N": {"type": "number"},
                                "lee_share": _num}},
        "thickness": {"type": "object", "required": ["min_m", "n_nan"]},
        "zones": {"type": "object", "required": ["bulk_area_m2_at_2mm", "max_layer_m"]},
        "depths": {"type": "object", "required": ["normals", "n_marched", "n_wrong_way", "n_left_body"],
                   "properties": {"normals": {"enum": ["derived", "facet"]}}},
    },
}
CHECK_SCHEMA = {"type": "object", "required": ["value", "threshold", "source", "passed"],
                "properties": {"source": {"type": "string"}, "passed": {"type": ["boolean", "null"]}}}
PREPARE_SCHEMA = {
    "type": "object",
    "required": ["schema", "schema_version", "contract_version", "name", "created_utc", "spheral_frag_version",
                 "git_commit", "fe_run", "flight", "contract", "material_table", "frames", "checks", "timing_s",
                 "size_bytes", "depth_normals", "passed"],
    "properties": {
        "schema": {"const": "spheral_frag.prepare"}, "schema_version": {"const": 1},
        "contract_version": {"const": contract.CONTRACT_VERSION}, "name": {"type": "string"},
        "created_utc": {"type": "string", "pattern": r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$"},
        "spheral_frag_version": {"type": "string"}, "git_commit": {"type": "string"},
        "fe_run": {"type": "object", "required": ["name", "dir", "json_sha256", "csv_sha256", "package", "settings"],
                   "properties": {"package": {"type": "object", "required": ["path", "file", "sha256"]},
                                  "settings": {"type": "object"}}},
        "flight": {"type": "object", "required": ["v_hat", "diameter_m", "initial_mass_kg", "material",
                                                  "macro_step_s", "frames_every", "seed"],
                   "properties": {"v_hat": {"type": "array", "items": {"type": "number"}, "minItems": 3,
                                            "maxItems": 3}}},
        "contract": {"type": "object", "required": ["present", "absent_optional", "unconfirmed"],
                     "properties": {k: {"type": "array", "items": {"type": "string"}}
                                    for k in ("present", "absent_optional", "unconfirmed")}},
        "material_table": {"type": "object", "required": ["file", "fe_material", "provisional",
                                                          "max_abs_dh_J_per_kg", "fl_bitwise"]},
        "frames": {"type": "array", "items": FRAME_SCHEMA},
        "checks": {"type": "object", "additionalProperties": CHECK_SCHEMA},
        "timing_s": {"type": "object", "required": ["total", "per_frame_median"]},
        "size_bytes": _int, "depth_normals": {"enum": ["derived", "facet", "mixed", None]},
        "passed": {"type": "boolean"},
    },
}


@pytest.fixture(scope="module")
def table_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("material")
    syn.synthetic_material_table(str(d / "material_table.npz"))
    return str(d)


def run_prepare(fe_run, outdir, table_dir, *extra):
    return main(["prepare", "--fe-run", str(fe_run), "--material-table", table_dir, "--outdir", str(outdir),
                 "--quiet", *extra])


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, table_dir):
    out = tmp_path_factory.mktemp("prepare")
    assert run_prepare(RUN_DIR, out, table_dir) == 0
    d = os.path.join(str(out), PREP_NAME)
    with open(os.path.join(d, "prepare.json")) as fh:
        return d, json.load(fh)


def file_bytes(d, names):
    return {n: open(os.path.join(d, n), "rb").read() for n in names}


PRODUCTS = ["frames/frame_00000.npz", "frames/frame_00001.npz", "frames/frame_00002.npz", "loads.npz", "flight.npz"]


# --------------------------------------------------------------------------------------------------- the outputs
def test_prepare_writes_every_file_and_a_valid_json(prepared):
    d, doc = prepared
    for name in PRODUCTS + ["material_table.npz", "material_table.json", "loads.json", "prepare.json"]:
        assert os.path.isfile(os.path.join(d, name)), name
    jsonschema.validate(doc, PREPARE_SCHEMA)
    assert doc["name"] == PREP_NAME and doc["passed"] is True
    assert [f["k"] for f in doc["frames"]] == [0, 1, 2] and [f["time_s"] for f in doc["frames"]] == [0.0, 0.5, 1.0]
    assert all(c["passed"] is not False for c in doc["checks"].values())
    assert doc["flight"] == {**doc["flight"], "v_hat": [1.0, 0.0, 0.0], "diameter_m": 0.1, "frames_every": 1,
                             "macro_step_s": 0.5, "seed": syn.SEED, "material": syn.SYNTHETIC_MATERIAL}
    assert doc["contract"]["missing"] == [] and doc["contract"]["absent_optional"] == ["n_derived"]
    assert doc["depth_normals"] == "facet" and all(f["depths"]["normals"] == "facet" for f in doc["frames"])
    for f in doc["frames"]:
        with open(os.path.join(d, f["file"]), "rb") as fh:
            import hashlib
            assert hashlib.sha256(fh.read()).hexdigest() == f["sha256"]
    # the numbers the fixture's README gives (tests/fixtures/spheral_frag/README.md)
    f0, f1, f2 = doc["frames"]
    assert (f0["n_tets"], f0["n_faces"], f2["n_tets"], f2["n_faces"], f2["n_nodes"]) == (2497, 1194, 2487, 1198, 752)
    assert f0["loads_valid"] is False and f1["loads_valid"] is True
    assert abs(f1["drag"]["patch_N"] - 6.866822220) < 5e-9 and abs(f2["drag"]["patch_N"] - 7.311101107) < 5e-9
    assert all(f["drag"]["normals"] == "facet" for f in doc["frames"])          # no n_derived: the table at theta
    assert math.isclose(f1["drag"]["smooth_N"], f1["drag"]["table_lee_N"], rel_tol=1e-12)
    assert max(f["mass_rel_diff"] for f in doc["frames"]) <= 5e-9       # the CSV's 9 digits (Task 4)
    assert f2["surface"]["n_open_directed_edges"] == 0 and f2["surface"]["volume_rel_diff"] <= 1e-12
    assert doc["checks"]["material_table_h"]["passed"] is None          # --material-table: no FE material to check


def test_prepared_frames_loads_and_flight_read_back(prepared):
    d, doc = prepared
    run = frames.read_fe_run(RUN_DIR)
    fr = frames.load_prepared_frame(os.path.join(d, "frames", "frame_00002.npz"))
    assert {"thickness", "slurry_depth", "liquid_depth", "area", "normal", "centroid", "theta"} <= set(fr.patch)
    assert "n_derived" not in fr.patch and "slurry_depth_facet" not in fr.patch    # no derived normals: no diagnostic
    z = dict(np.load(os.path.join(d, "loads.npz")))
    assert z["p"].shape == z["tau"].shape == z["area"].shape == (3, 180) and list(z["k"]) == [0, 1, 2]
    tables = loads.unstack(z)
    assert not tables[0].has_loads and tables[1].has_loads and tables[2].theta_last_deg == 149.5
    # the table of frame 1 is the one loads.build_table makes of the prepared frame's own patches
    f1 = frames.load_prepared_frame(os.path.join(d, "frames", "frame_00001.npz"))
    want = loads.build_table(f1.patch["theta"], f1.patch["area"], f1.patch["p_w"], f1.patch["tau"],
                             run.history["p_w_stag_Pa"][1])
    assert np.array_equal(tables[1].p, want.p, equal_nan=True) and np.array_equal(tables[1].tau, want.tau,
                                                                                  equal_nan=True)
    hdr = json.load(open(os.path.join(d, "loads.json")))
    assert hdr["lee_label"] == loads.LEE_LABEL and hdr["lee_params"] is None
    fl = dict(np.load(os.path.join(d, "flight.npz")))
    for f in contract.fields("history"):
        assert fl[f.key].tobytes() == run.history[f.name].tobytes(), f.key
    assert np.array_equal(fl["deceleration_ms2"], run.history["load_factor_g"] * loads.G0)
    assert np.array_equal(fl["drag_N"], run.history["mass_kg"] * run.history["load_factor_g"] * loads.G0)
    assert np.all((fl["T_air_K"] > 180.0) & (fl["T_air_K"] < 220.0))   # US76 near 77.5 km: about 200 K
    flight = record.FlightTable.load(os.path.join(d, "flight.npz"))    # the record's reader takes it
    assert flight.at(0.25)["T_air_K"] == pytest.approx(0.5 * (fl["T_air_K"][0] + fl["T_air_K"][1]), rel=1e-15)
    t = material.MaterialTable.load(d)
    assert t.name == syn.SYNTHETIC_MATERIAL and doc["material_table"]["fe_material"] == syn.SYNTHETIC_MATERIAL


# --------------------------------------------------------------------------------------------------- resume
def test_second_run_is_skipped_and_force_rebuilds_byte_identical(prepared, table_dir, capsys):
    d, doc = prepared
    out = os.path.dirname(d)
    before = file_bytes(d, PRODUCTS)
    mtime = os.path.getmtime(os.path.join(d, "prepare.json"))
    capsys.readouterr()
    assert run_prepare(RUN_DIR, out, table_dir) == 0
    assert "prepared already" in capsys.readouterr().out
    assert os.path.getmtime(os.path.join(d, "prepare.json")) == mtime
    assert run_prepare(RUN_DIR, out, table_dir, "--force") == 0
    assert file_bytes(d, PRODUCTS) == before
    doc2 = json.load(open(os.path.join(d, "prepare.json")))
    assert doc2["timing_s"]["n_frames_reused"] == 0 and doc2["created_utc"] >= doc["created_utc"]


def test_interrupted_prepare_rebuilds_only_the_missing_frames(tmp_path, prepared, table_dir):
    d0, _ = prepared
    out = tmp_path / "prep"
    shutil.copytree(d0, out / PREP_NAME)
    d = str(out / PREP_NAME)
    want = file_bytes(d, PRODUCTS)
    os.remove(os.path.join(d, "prepare.json"))                     # interrupted: no summary,
    os.remove(os.path.join(d, "frames", "frame_00001.npz"))        # one frame never written
    with open(os.path.join(d, "frames", "frame_00002.npz"), "ab") as fh:
        fh.write(b"x")                                             # and one whose sha256 no longer matches
    assert run_prepare(RUN_DIR, out, table_dir) == 0
    doc = json.load(open(os.path.join(d, "prepare.json")))
    assert doc["timing_s"]["n_frames_reused"] == 1 and doc["timing_s"]["n_frames_built"] == 2
    assert file_bytes(d, PRODUCTS) == want


# --------------------------------------------------------------------------------------------------- failures
def test_run_without_p_w_exits_1_naming_it(tmp_path, table_dir, capsys):
    run = frames.read_fe_run(RUN_DIR)
    lack = []
    for k, _, _ in run.frames:
        f = frames.read_frame(run, k)
        f.patch.pop("p_w")
        lack.append(f)
    syn.write_fe_run(str(tmp_path / "fe"), lack, run.history, json.load(open(run.json_path)))
    capsys.readouterr()
    assert run_prepare(tmp_path / "fe" / run.name, tmp_path / "prep", table_dir) == 1
    err = capsys.readouterr().err
    assert "p_w" in err and "misses required contract items" in err
    doc = json.load(open(tmp_path / "prep" / PREP_NAME / "prepare.json"))
    assert doc["passed"] is False and doc["contract"]["missing"] == ["p_w"]
    assert doc["checks"]["contract_items"] == {**doc["checks"]["contract_items"], "value": 1, "passed": False}
    assert len(doc["frames"]) == 1                                # stopped at the first frame


def test_missing_items_of_the_run_are_listed_and_exit_1(tmp_path, table_dir, capsys):
    shutil.copytree(FIXTURE, tmp_path / "fe")
    p = tmp_path / "fe" / (syn.SPHERE_RUN_NAME + ".json")
    doc = json.loads(p.read_text())
    del doc["settings"]["macro_step_s"], doc["settings"]["frames_every"]
    p.write_text(json.dumps(doc))
    capsys.readouterr()
    assert run_prepare(tmp_path / "fe" / syn.SPHERE_RUN_NAME, tmp_path / "prep", table_dir) == 1
    err = capsys.readouterr().err
    assert "settings.macro_step_s" in err and "settings.frames_every" in err   # every missing item, in one message


def test_bad_arguments_and_missing_inputs_exit_2(tmp_path, table_dir):
    assert run_prepare(tmp_path / "nowhere", tmp_path / "prep", table_dir) == 2
    for frames_arg in ("2:1", "0:7", "a:b", "1:2:3"):
        assert run_prepare(RUN_DIR, tmp_path / "prep", table_dir, "--frames", frames_arg) == 2, frames_arg
    assert run_prepare(RUN_DIR, tmp_path / "prep", table_dir, "--every", "0") == 2
    assert run_prepare(RUN_DIR, tmp_path / "prep", str(tmp_path / "no_table")) == 2
    assert main(["prepare"]) == 2                                  # --fe-run is required
    # the run's own material is no finite-element material: without --material-table the material is missing
    assert main(["prepare", "--fe-run", RUN_DIR, "--outdir", str(tmp_path / "prep"), "--quiet"]) == 2
    assert not os.path.exists(tmp_path / "prep" / PREP_NAME / "prepare.json")


def test_command_line_exit_code_in_a_subprocess(tmp_path):
    r = subprocess.run([sys.executable, "-m", "spheral_frag", "prepare", "--fe-run", str(tmp_path / "nowhere"),
                        "--outdir", str(tmp_path)], cwd=REPO_ROOT, capture_output=True, text=True, timeout=120)
    assert r.returncode == 2 and "no finite-element run directory" in r.stderr


def test_frame_selection_and_name(tmp_path, table_dir):
    assert run_prepare(RUN_DIR, tmp_path, table_dir, "--frames", "0:2", "--every", "2", "--no-thickness") == 0
    name = naming.prepare_name(syn.SPHERE_RUN_NAME, 0, 2, every=2)
    doc = json.load(open(tmp_path / name / "prepare.json"))
    assert [f["k"] for f in doc["frames"]] == [0, 2] and doc["selection"] == {"k0": 0, "k1": 2, "every": 2,
                                                                              "n_frames": 2}
    assert "thickness" not in doc["frames"][0] and doc["depth_normals"] is None
    fr = frames.load_prepared_frame(str(tmp_path / name / "frames" / "frame_00002.npz"))
    assert "thickness" not in fr.patch and "slurry_depth" not in fr.patch


# ------------------------------------------------------------------------------------------ decision 9: normals
def test_depths_march_along_derived_normals_when_the_frame_carries_them(tmp_path, prepared, table_dir):
    d0, doc0 = prepared
    syn.with_derived_normals(FIXTURE, str(tmp_path / "fe"))
    assert run_prepare(tmp_path / "fe" / syn.SPHERE_RUN_NAME, tmp_path / "prep", table_dir) == 0
    d = str(tmp_path / "prep" / PREP_NAME)
    doc = json.load(open(os.path.join(d, "prepare.json")))
    jsonschema.validate(doc, PREPARE_SCHEMA)
    assert doc["depth_normals"] == "derived" and "n_derived" not in doc["contract"]["absent_optional"]
    for e, e0 in zip(doc["frames"], doc0["frames"]):
        assert e["depths"]["normals"] == "derived"
        # the facet-normal diagnostic is exactly the depth prepare computes without derived normals
        assert e["depths"]["facet"] == {k: e0["depths"][k] for k in e["depths"]["facet"]}
        assert e["n_derived"]["max_norm_error"] < 1e-15
        # radial against facet: none on the intact sphere; on frame 2 three of the staircase's side walls (measured)
        assert e["n_derived"]["n_against_facet"] == (3 if e["k"] == 2 else 0)
    for k in (1, 2):
        fr = frames.load_prepared_frame(os.path.join(d, "frames", frames.prepared_frame_name(k)))
        fr0 = frames.load_prepared_frame(os.path.join(d0, "frames", frames.prepared_frame_name(k)))
        assert fr.patch["n_derived"].shape == (len(fr.faces), 3)
        assert fr.patch["slurry_depth_facet"].tobytes() == fr0.patch["slurry_depth"].tobytes()
        assert fr.patch["liquid_depth_facet"].tobytes() == fr0.patch["liquid_depth"].tobytes()
        # the radial march: NaN (pointing out of the body) exactly where the radial direction is against the facet
        # normal -- three staircase side walls on frame 2, none on the intact frame 1 -- and different from the facet
        # march where the facets are tilted
        against = np.einsum("ij,ij->i", fr.patch["n_derived"], fr.patch["normal"]) <= 0.0
        assert np.array_equal(~np.isfinite(fr.patch["slurry_depth"]), against)
        assert np.count_nonzero(against) == (3 if k == 2 else 0)
        assert not np.array_equal(fr.patch["slurry_depth"], fr0.patch["slurry_depth"])
    # on the staircase frame the wrong-way radial rays are counted, and fewer rays leave the body than along facets
    e2 = doc["frames"][2]["depths"]
    assert e2["n_wrong_way"] == 3 and doc["frames"][2]["depths"]["facet"]["n_wrong_way"] == 0
    assert e2["n_left_body"] <= e2["facet"]["n_left_body"]
    # the drag reads the table at the radial normals' inclination; on the intact frame 1 they are the facets' within
    # the faceting, so the two drags agree to the fixture's faceting error (plan Task 7: below 1 %)
    for e in doc["frames"][1:]:
        assert e["drag"]["normals"] == "derived" and math.isfinite(e["drag"]["smooth_N"])
    f1 = doc["frames"][1]["drag"]
    assert abs(f1["smooth_N"] / f1["table_lee_N"] - 1.0) < 1e-2
