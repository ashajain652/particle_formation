"""The shared package of the combined refactor (stage 0 + S): the old module paths are the shared modules themselves,
the shared package imports nothing from either architecture, and the data paths still resolve after the move one
directory down. Later tasks extend MOVED."""
import ast
import glob
import importlib
import os
import subprocess

import pytest

from helpers import PROJECT_ROOT, PY

MOVED = ["trajectory", "earth", "atmosphere", "fap", "aero", "gas", "constants", "heating", "material", "sesam_io",
         "compare", "viz", "rheology", "dispersion", "surface_flow", "spray", "girin_case"]
ROOT_CONSTANTS = ("__version__", "PROJECT_ROOT", "REPO_ROOT", "DATA_DIR", "OUTPUT_DIR", "DRAMA_DIR", "DRAMA_DATA_DIR")


@pytest.mark.parametrize("name", MOVED)
def test_old_path_is_the_shared_module(name):
    old = importlib.import_module("reentry_model." + name)
    new = importlib.import_module("reentry_model.shared." + name)
    assert old is new
    assert getattr(importlib.import_module("reentry_model"), name) is new


@pytest.mark.parametrize("name", MOVED)
def test_import_order_does_not_matter(name):
    code = ("import reentry_model.shared.{n} as new\n"
            "import reentry_model.{n} as old\n"
            "from reentry_model import {n} as attr\n"
            "print(old is new and attr is new)\n").format(n=name)
    out = subprocess.run([PY, "-c", code], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout.strip()
    assert out == "True"


def test_data_paths_resolve():
    from reentry_model.shared import aero, fap, material, sesam_io
    for path in (material.DEFAULT_MATERIAL, aero.DEFAULT_ATDB, fap.DEFAULT_FAP_DAY, sesam_io.MELTING_REFERENCE_DIR):
        assert os.path.exists(path), path


def test_shared_modules_import_nothing_outside_shared():
    """`from .. import <root constant>` is the only way out of the shared package."""
    violations = []
    for path in sorted(glob.glob(os.path.join(PROJECT_ROOT, "reentry_model", "shared", "*.py"))):
        with open(path) as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            bad = False
            if isinstance(node, ast.ImportFrom):
                if node.level == 0:
                    mod = node.module or ""
                    bad = mod.startswith("reentry_model") and not (mod == "reentry_model.shared"
                                                                   or mod.startswith("reentry_model.shared."))
                elif node.level == 2:
                    bad = not (node.module is None and all(a.name in ROOT_CONSTANTS for a in node.names))
                elif node.level >= 3:
                    bad = True
            elif isinstance(node, ast.Import):
                bad = any(a.name.startswith("reentry_model") and not (a.name == "reentry_model.shared"
                                                                     or a.name.startswith("reentry_model.shared."))
                          for a in node.names)
            if bad:
                violations.append((os.path.relpath(path, PROJECT_ROOT), node.lineno))
    assert violations == []


def test_pickles_naming_old_paths_load():
    """A state pickled before the refactor names reentry_model.spray; it must load as the shared class."""
    import pickle
    from reentry_model.shared import material, spray
    mat = material.Material.from_drama_json("AA7075_scheil")
    model = spray.SprayModel(mat.liquid, k_r=1.0, k_t=1.0, we_critical=3.08, rt_spray=True, spray_gates="2026-11-10",
                             slip_band=True)
    blob = pickle.dumps(model, protocol=0)
    assert b"reentry_model.shared.spray" in blob
    old = blob.replace(b"reentry_model.shared.spray", b"reentry_model.spray")
    assert b"reentry_model.shared.spray" not in old
    back = pickle.loads(old)
    assert type(back) is spray.SprayModel
    # DispersionTable has no __eq__, so equal state means the same bytes when pickled again
    assert pickle.dumps(back, protocol=0) == blob


def test_coupled_reexports_the_droplet_writer():
    import inspect
    from reentry_model import coupled
    from reentry_model.shared import droplets
    assert coupled.write_particles is droplets.write_particles
    assert coupled.SUMMARY_COLUMNS is droplets.SUMMARY_COLUMNS
    assert coupled.SKIN_SUMMARY_COLUMNS is droplets.SKIN_SUMMARY_COLUMNS
    src = inspect.getsource(droplets)
    for bad in ("skin", "column"):
        assert not [ln for ln in src.splitlines()
                    if ln.lstrip().startswith(("import ", "from ")) and bad in ln.split("#")[0]]


def test_body_and_thermal_reexport_the_shared_names():
    import pickle
    import numpy as np
    from reentry_model import body, thermal
    from reentry_model.shared import body as sbody, thermal as sthermal
    for name in ("Body", "sphere_mass", "ConstantBody"):
        assert getattr(body, name) is getattr(sbody, name)
    for name in ("SIGMA_SB", "MissingBackend", "StepResult", "ThermalSolver"):
        assert getattr(thermal, name) is getattr(sthermal, name)
    r = sthermal.StepResult(T=np.array([1.0, 2.0]), Q_conv=3.0, Q_rad=4.0, iterations=5, Q_extra=6.0,
                            Q_dropped=7.0, Q_interface=8.0)
    back = pickle.loads(pickle.dumps(r))
    assert type(back) is sthermal.StepResult
    assert np.array_equal(back.T, r.T) and (back.Q_conv, back.Q_interface) == (3.0, 8.0)
    # the shared modules import nothing of the package outside shared
    import inspect
    for mod in (sbody, sthermal):
        assert not [ln for ln in inspect.getsource(mod).splitlines()
                    if ln.lstrip().startswith(("import reentry_model", "from reentry_model", "from .."))]


def test_cli_common_reexports_and_paths():
    import reentry_model
    from reentry_model import cli
    from reentry_model.shared import cli_common
    for name in ("DEFAULT_HEADING_DEG", "DEFAULT_LAT_DEG", "DEFAULT_LON_DEG", "DEFAULT_EPOCH", "DEFAULT_MATERIAL_DENSITY",
                 "ATMOSPHERES", "WINDS", "parse_epoch", "make_atmosphere", "git_commit", "provenance", "DEFAULT_OUTDIR"):
        assert getattr(cli, name) is getattr(cli_common, name), name
    assert cli_common.DEFAULT_OUTDIR == reentry_model.OUTPUT_DIR


def test_spheral_entry_points_resolve():
    """Spheral/spheral_frag/fe.py loads reentry_model.cli.make_atmosphere and reentry_model.material by name."""
    from reentry_model import cli, material
    assert callable(cli.make_atmosphere)
    assert hasattr(material.Material, "from_drama_json")
    mat = material.Material.from_drama_json(os.path.join(PROJECT_ROOT, "data", "materials", "AA7075_scheil.json"))
    for attr in ("_T_h", "_h_nodes", "_latent_slope"):
        assert hasattr(mat, attr), attr
