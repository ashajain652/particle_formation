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


SHARED_DIR = os.path.join(PROJECT_ROOT, "reentry_model", "shared")
SHARED_MODULES = sorted(MOVED + ["droplets", "body", "thermal", "cli_common"])   # the 21 modules besides __init__.py


def _outside_shared(node, own_name, shared_names):
    """True when the import statement `node`, in the shared module `own_name`, reaches outside reentry_model.shared."""
    def in_shared(dotted):
        return dotted == "reentry_model.shared" or dotted.startswith("reentry_model.shared.")
    if isinstance(node, ast.Import):
        return any(a.name.split(".")[0] == "reentry_model" and not in_shared(a.name) for a in node.names)
    if node.level == 0:
        mod = node.module or ""
        if mod.split(".")[0] != "reentry_model":
            return False
        if mod == "reentry_model":    # `from reentry_model import x` names the package root
            return True
        return not in_shared(mod)
    if node.level == 1:               # `from . import x` / `from .x import y`: x must be a shared module
        if node.module is None:
            return not all(a.name in shared_names for a in node.names)
        return node.module.split(".")[0] not in shared_names
    if node.level == 2:               # `from .. import CONSTANT` or `from ..shared[.x] import y`
        if node.module is None:
            return not all(a.name in ROOT_CONSTANTS or a.name == "shared" for a in node.names)
        return node.module.split(".")[0] != "shared"
    return True


def test_shared_modules_import_nothing_outside_shared_statically():
    """The only ways out of reentry_model.shared are `from .. import <package-root constant>`; relative imports of
    level 1 (`from . import x`, `from .x import y`) must name shared modules, absolute ones must name reentry_model.shared."""
    paths = sorted(glob.glob(os.path.join(SHARED_DIR, "*.py")))
    names = [os.path.basename(p)[:-3] for p in paths]
    assert sorted(n for n in names if n != "__init__") == SHARED_MODULES and "__init__" in names and len(paths) == 22
    shared_names = set(names)
    violations, scanned_imports = [], 0
    for path, name in zip(paths, names):
        with open(path) as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                scanned_imports += 1
                if _outside_shared(node, name, shared_names):
                    violations.append((os.path.relpath(path, PROJECT_ROOT), node.lineno))
    assert scanned_imports > 0
    assert violations == []


def test_the_static_scan_catches_each_kind_of_escape():
    """The scan itself: each of these statements, in a shared module, must be flagged; the allowed ones must not."""
    shared = set(SHARED_MODULES)

    def flagged(src):
        return _outside_shared(ast.parse(src).body[0], "aero", shared)
    for bad in ("from .skin import x", "from . import skin", "from .. import skin", "from ..skin import x",
                "from ... import x", "from reentry_model import skin", "from reentry_model.skin import x",
                "import reentry_model.skin", "from reentry_model import cli"):
        assert flagged(bad), bad
    for good in ("from . import material", "from .material import x", "from .. import PROJECT_ROOT, DATA_DIR",
                 "from reentry_model.shared import material", "from reentry_model.shared.material import x",
                 "import reentry_model.shared.material", "import numpy", "from scipy import optimize"):
        assert not flagged(good), good


def test_shared_package_imports_neither_architecture():
    """Importing every shared module, in a fresh interpreter, loads none of the skyscrapers modules."""
    code = ("import importlib, sys\n"
            "for n in {mods!r}:\n"
            "    importlib.import_module('reentry_model.shared.' + n)\n"
            "banned = {banned!r}\n"
            "bad = sorted(m for m in sys.modules if m in banned or (m.startswith('reentry_model.thermal.') "
            "and m.endswith('_backend')))\n"
            "assert any(m.startswith('reentry_model.shared.') for m in sys.modules)\n"
            "print(bad)\n").format(
        mods=SHARED_MODULES,
        banned={"reentry_model." + n for n in ("skin", "skinbody", "columns", "movemesh", "mesh", "spreading", "film",
                                               "restart", "coupled", "cli", "body", "skyscrapers")}
        | {"reentry_model.skyscrapers." + n for n in ("skin", "skinbody", "columns", "movemesh", "mesh", "spreading",
                                                      "restart", "coupled", "body")})
    out = subprocess.run([PY, "-c", code], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout.strip()
    assert out == "[]", out


def _reentry_imports(path):
    """[(module, attribute-or-None)] for every reentry_model import in a Python file (importlib.import_module with a
    literal or `_PKG + ".x"` name included, as Spheral's fe.py loads it)."""
    with open(path) as fh:
        tree = ast.parse(fh.read())
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found += [(a.name, None) for a in node.names if a.name.split(".")[0] == "reentry_model"]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and (node.module or "").split(".")[0] == "reentry_model":
            found += [(node.module, a.name) for a in node.names]
        elif isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "import_module" and node.args:
            arg = node.args[0]
            if isinstance(arg, ast.BinOp) and isinstance(arg.op, ast.Add) and isinstance(arg.left, ast.Name) \
                    and arg.left.id == "_PKG" and isinstance(arg.right, ast.Constant):
                found.append(("reentry_model" + arg.right.value, None))
            elif isinstance(arg, ast.Constant) and str(arg.value).split(".")[0] == "reentry_model":
                found.append((arg.value, None))
    return found


def test_every_reentry_model_import_in_tests_and_analysis_resolves():
    files = sorted(glob.glob(os.path.join(PROJECT_ROOT, "tests", "*.py")) +
                   glob.glob(os.path.join(PROJECT_ROOT, "analysis", "**", "*.py"), recursive=True))
    fe = os.path.join(os.path.dirname(PROJECT_ROOT), "Spheral", "spheral_frag", "fe.py")
    assert os.path.isfile(fe)
    files.append(fe)
    seen, unresolved = set(), []
    for path in files:
        for mod, attr in _reentry_imports(path):
            if (mod, attr) in seen:
                continue
            seen.add((mod, attr))
            try:
                m = importlib.import_module(mod)
                if attr is not None and not hasattr(m, attr):
                    importlib.import_module(mod + "." + attr)       # `from package import submodule`
            except ModuleNotFoundError as exc:
                if exc.name in ("dolfinx", "ufl", "basix", "mpi4py", "petsc4py"):
                    continue                                         # the optional FEniCSx backend, not installed
                unresolved.append((os.path.relpath(path, PROJECT_ROOT), mod, attr, repr(exc)))
            except Exception as exc:
                unresolved.append((os.path.relpath(path, PROJECT_ROOT), mod, attr, repr(exc)))
    assert len(seen) > 50, len(seen)
    assert unresolved == []


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
