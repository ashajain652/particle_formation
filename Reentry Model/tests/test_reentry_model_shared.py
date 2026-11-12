"""The shared package of the combined refactor (stage 0 + S): the old module paths are the shared modules themselves,
the shared package imports nothing from either architecture, and the data paths still resolve after the move one
directory down. Later tasks extend MOVED."""
import ast
import glob
import importlib
import os
import subprocess

import numpy as np
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
SHARED_MODULES = sorted(MOVED + ["droplets", "body", "thermal", "cli_common", "lubrication", "nose"])   # the 23 modules
                                                       # besides __init__.py (nose: level-set stage 3 Task 1)


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
    assert sorted(n for n in names if n != "__init__") == SHARED_MODULES and "__init__" in names and len(paths) == 24
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
                                                      "restart", "coupled", "body")}
        | {"reentry_model.skyscrapers.elements", "reentry_model.skyscrapers.elements.film"})
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


def test_moved_film_and_nose_names_are_the_shared_objects():
    """Level-set stage 3 Task 1: the physics the level set needs from the skyscrapers modules lives in shared, and the old
    modules re-export the shared objects themselves."""
    from reentry_model import body, film, skinbody
    from reentry_model.shared import lubrication, nose, spray
    assert film.lubrication is lubrication.lubrication
    for name in ("fit_sphere", "NOSE_CAP_ANGLE", "NOSE_CAP_FACTOR"):
        assert getattr(body, name) is getattr(nose, name), name
    assert skinbody.SLURRY_SIGMA is spray.SLURRY_SIGMA
    assert skinbody.region_extent is spray.region_extent
    assert skinbody.LEEWARD_TESTS is nose.LEEWARD_TESTS
    assert (nose.NOSE_CAP_ANGLE, nose.NOSE_CAP_FACTOR, spray.SLURRY_SIGMA) == (30.0, 1.67, 0.80)
    assert nose.LEEWARD_TESTS == ("facing", "centroid")


# The values below were recorded on the unmodified tree (commit 2291ded) by calling the old code itself:
# SkinBody.nose_fit's property (skinbody.py 872-889) on a stand-in with `_nose_raw` = (fitted, bound, points, concave),
# `size_feedback` and a mesh radius of 0.05 m; SkinBody.load_angles (824-833) on a stand-in with the given theta and
# windward corners; and SkinBody._leeward's plane test (810-816) on the given points, centroid and flight direction.
NOSE_CASES = [  # (fitted, bound, points, concave) -> used under "current"; "initial" gave 0.05 for every one
    ((0.04, 0.06, 7, False), 0.04), ((0.09, 0.06, 7, False), 0.06), ((0.03, 0.06, 7, True), 0.06),
    ((float("nan"), 0.06, 9, True), 0.06), ((float("nan"), 0.06, 3, False), 0.05), ((0.04, float("nan"), 6, False), 0.04),
    ((0.03, float("nan"), 6, True), 0.05), ((0.04, 0.0, 6, False), 0.04)]
LOAD_THETA = [0.0, 0.3, 1.2, 1.5707963267948966, 1.9, 2.8, 1.0]
LOAD_WINDWARD = [True, False, True, False, False, True, False]
LOAD_CENTROID = [0.0, 1.5707973267948965, 1.2, 1.5707973267948965, 1.9, 2.8, 1.5707973267948965]
AHEAD_POINTS = [[0.01, 0.0, 0.0], [-0.02, 0.01, 0.0], [0.0, 0.5, 0.3], [0.003, -0.01, 0.02], [-1e-9, 0.0, 0.0]]
AHEAD_CENTRE = [0.0, 0.1, -0.2]
AHEAD = {(1.0, 0.0, 0.0): [True, False, False, True, False], (0.6, 0.8, 0.0): [False, False, True, False, False]}


def test_shared_nose_rules_match_the_old_inline_rules():
    import numpy as np
    from reentry_model.shared import nose
    for (fitted, bound, n, concave), used in NOSE_CASES:
        assert nose.nose_radius(fitted, bound, concave, n, 0.05, "current") == used, (fitted, bound, n, concave)
        assert nose.nose_radius(fitted, bound, concave, n, 0.05, "initial") == 0.05
    got = nose.load_angles(np.array(LOAD_THETA), np.array(LOAD_WINDWARD))
    assert got.tolist() == LOAD_CENTROID                                   # bit for bit
    for v, expected in AHEAD.items():
        assert nose.ahead_of(np.array(AHEAD_POINTS), np.array(AHEAD_CENTRE), np.array(v)).tolist() == expected


def test_shared_slurry_properties_are_the_skins_old_arithmetic():
    """The slurry's properties as skinbody.slurry_modes built them at 2291ded (lines 528-532 and 554): Li's viscosity with
    the bridge at the cell's solid fraction and shear rate, the given mixture density, SLURRY_SIGMA."""
    import numpy as np
    from reentry_model.shared import rheology, spray
    from reentry_model.shared.material import Material
    mat = Material.from_drama_json(os.path.join(PROJECT_ROOT, "data", "materials", "AA7075_scheil.json"))
    T = np.array([880.0, 895.5, 900.0, 905.0, 912.0])
    f_s = 1.0 - mat.liquid_fraction(T)
    gamma = np.array([0.0, 10.0, 1e3, 1e5, 1e7])
    rho = np.array([2700.0, 2650.0, 2600.0, 2550.0, 2500.0])
    f_bridge = 1.0 - float(mat.liquid_fraction(np.array([rheology.BRIDGE_LOW]))[0])     # old line 528
    old_mu = rheology.viscosity(T, f_s, gamma, mat.liquid.mu, f_bridge)                   # old line 530
    props = spray.slurry_properties(mat, T, f_s, gamma, rho)
    assert np.array_equal(props.mu, old_mu) and np.array_equal(props.rho, rho) and props.sigma == spray.SLURRY_SIGMA


def test_size_feedback_names_is_the_shared_tuple():
    from reentry_model import body
    from reentry_model.shared import nose
    assert body.SIZE_FEEDBACK_NAMES is nose.SIZE_FEEDBACK_NAMES and nose.SIZE_FEEDBACK_NAMES == ("current", "initial")


def _cap(R, sign, n=12):
    a, b = np.linspace(0, 0.5, n), np.linspace(0, 2 * np.pi, n, endpoint=False)
    A, B = (m.ravel() for m in np.meshgrid(a, b))
    return np.column_stack([R * np.cos(A) * sign, R * np.sin(A) * np.cos(B), R * np.sin(A) * np.sin(B)])


def test_fit_cap_reproduces_the_old_inline_decision():
    """Old values recorded from SkinBody._fit_nose's inline code at commit 799375d, flight direction +x."""
    from reentry_model.shared import nose
    g = np.linspace(-0.01, 0.01, 5)
    Y, Z = np.meshgrid(g, g)
    cases = {"convex": (_cap(0.03, 1) - [0.03, 0, 0], (0.02999999999999999, False)),
             "concave": (_cap(0.03, -1) + [0.03, 0, 0], (0.02999999999999999, True)),
             "flat": (np.column_stack([np.zeros(25), Y.ravel(), Z.ravel()]), (float("nan"), True)),
             "three": (_cap(0.03, 1)[:3], (float("nan"), False))}
    for name, (pts, (r, concave)) in cases.items():
        fitted, c = nose.fit_cap(pts, np.array([1.0, 0.0, 0.0]))
        assert c is concave or c == concave, name
        assert (np.isnan(fitted) and np.isnan(r)) or fitted == r, name


def test_moved_spray_plan_pieces_are_the_shared_objects():
    """Level-set stage 3 Task 8 (fix round 1): the spray plan's architecture-neutral pieces live in shared/spray.py and
    skinbody keeps its old names as the shared objects themselves."""
    from reentry_model import skinbody
    from reentry_model.shared import spray
    assert skinbody._rows is spray.rows_of
    assert skinbody._no_slurry is spray.no_slurry
    assert skinbody._no_step_spray is spray.no_step_spray
    assert set(spray.no_slurry(2)) == set(spray.NO_SLURRY_KEYS) | {"branch"}
    assert spray.no_slurry(2)["branch"].tolist() == [-1, -1]


def _spray_plan_case():
    """The skyscrapers' Girin plan with the slurry, on fixed skins and a fixed gas state (the scratch case the old values
    below were recorded on): ten columns of 10 um cells -- liquid on rigid metal, thin liquid, liquid on slurry, slurry
    over liquid, a crust, rigid alone, a warm slurry cell, liquid over slurry, liquid over deep warm slurry, deep warm
    slurry -- every cell moving; shock-layer gas with Girin's closure, tau 200 Pa, G -3e4 Pa/m, deceleration 30 m/s2."""
    from types import SimpleNamespace
    from reentry_model import skin as sky, skinbody
    from reentry_model.shared import spray, surface_flow as sf
    from reentry_model.shared.material import LiquidProperties, Material
    mat = Material.from_drama_json("AA7075_scheil")
    L, S, R = 950.0, 902.0, 880.0
    temps = [[L] * 30 + [R] * 10, [L] * 4 + [R] * 10, [L] * 6 + [S] * 20 + [R] * 10, [S] * 3 + [L] * 5 + [R] * 10,
             [R] * 3 + [L] * 5 + [R] * 10, [R] * 12, [L] * 2 + [907.0] + [S] * 25 + [R] * 5, [L] * 3 + [S] * 40 + [R] * 5,
             [L] * 2 + [907.8] * 300 + [R] * 5, [907.5] * 200 + [R] * 5]
    n = len(temps)
    f = lambda v: np.full(n, v, float)
    flow = SimpleNamespace(closure=np.full(n, sf.CLOSURE_GIRIN), branch=sf.BRANCH_SHOCK_LAYER, u_eff=f(2000.0),
                           u_e=f(2000.0), rho_e=f(7e-4), mu_e=f(1.6e-4), delta_a=f(7e-3), mach_e=f(1.0), tau=f(200.0),
                           G=f(-3e4), kn_local=f(1e-3), deceleration=30.0)
    state = SimpleNamespace(freestream=SimpleNamespace(rho=8e-5), V=7200.0, ma=24.0, h=71e3)
    area = np.full(n, 1e-6) * np.linspace(1, 2, n)
    sk = sky.Skins(mat, area.copy())
    count = np.array([len(t) for t in temps])
    K = count.max()
    T = np.array([t + [300.0] * (K - len(t)) for t in temps])
    valid = np.arange(K)[None, :] < count[:, None]
    m = np.where(valid, mat.rho * area[:, None] * 1e-5, 0.0)
    U = np.where(valid[..., None], np.random.default_rng(4).normal(0, 0.3, (n, K, 3)), 0.0)
    sk._install(m, m * mat.enthalpy(T), count, U=U)
    theta = np.linspace(0.05, 0.9, n)
    points = np.column_stack([np.arange(n) * 1e-3, np.sin(theta) * 1e-2, np.cos(theta) * 1e-2])
    model = spray.SprayModel(LiquidProperties(2400.0, 1.3e-3, 0.8), spray_gates="2026-11-10")
    plan = skinbody.girin_plan(sk, flow, state, model, theta, area, np.column_stack([np.arange(n - 1), np.arange(1, n)]),
                               0.05, 0.5, slurry=True, points=points, centre=np.array([0.0, 0.0, -0.02]),
                               v_hat=np.array([0.0, 0.0, -1.0]), nose_radius=0.04)
    rng = np.random.default_rng(1)
    plan.mass[:] = np.where(plan.rate > 0, rng.uniform(1e-9, 1e-8, n), 0.0)
    plan.mass_slurry[:] = np.where(plan.rate_slurry > 0, rng.uniform(1e-9, 1e-8, n), 0.0)
    for k in ("enthalpy", "enthalpy_slurry", "solid", "solid_slurry"):
        getattr(plan, k)[:] = rng.uniform(0, 1e-3, n) * (plan.mass + plan.mass_slurry > 0)
    plan.momentum[:] = rng.normal(0, 1e-8, (n, 3))
    plan.momentum_slurry[:] = rng.normal(0, 1e-8, (n, 3))
    return plan


def test_moved_spray_plan_pieces_reproduce_the_old_values():
    """Old values recorded at commit e7a6539 by calling skinbody's own code (girin_plan with slurry_modes inline,
    SprayPlan.count, histogram, release_angle, source_block and SkinBody._step_columns) on `_spray_plan_case`; the moved
    code gives every value recorded there to the last bit (checked on the whole record in a scratch script, 2026-11-12)."""
    from reentry_model import skinbody
    plan = _spray_plan_case()
    assert plan.rate_slurry[8:].tolist() == [5798.850104260878, 4199.991695125593]
    rows = [2, 8, 9]
    old = {"mu": [1.8045651910830642, 0.0024874216047943098, 0.006583524814691101], "gamma": [108.62127427485834, 367.0, 367.0],
           "delta_m": [0.022997419083111988, 0.0002889863669755167, 0.0005524317337923679],
           "depth": [0.0002218463108255324, 0.0035089344045234287, 0.0023324387010949447]}
    for key, values in old.items():
        assert plan.slurry[key][rows].tolist() == values, key
    assert np.isnan(plan.slurry["r"][2]) and plan.slurry["r"][[8, 9]].tolist() == [0.00025144144518461127,
                                                                                0.00048066016042328256]
    plan.count(2400.0)
    assert (plan.number.sum(), plan.number_slurry.sum()) == (1.3896693658723713, 0.020678677301138745)
    hn, hm = plan.histogram(np.vstack([plan.mass, plan.mass_slurry]))
    assert (hn.sum(), hm.sum()) == (1.41034804317351, 4.2728380394538384e-08) and np.flatnonzero(hn).tolist() == [20, 22, 24, 26]
    assert plan.release_angle()[[0, 9]].tolist() == [0.9548412207142292, 24.471238492884595]
    block = plan.source_block(1.5)
    sums = block.sum(axis=0)
    assert block.shape == (8, 27)
    assert sums[[3, 9, 11, 12, 15, 19, 24, 25, 26]].tolist() == [244.78030247533502, 0.010864095740776177, 3291.286395003719,
                                                                 0.0016554610177059712, 21.79251679995043, 17.16381983157551,
                                                                 1157896.875534645, 705378.7740446577, 113.96312274432282]
    assert skinbody.SkinBody._step_columns(plan, block) == {
        "r_median_um": 163.24331099104404, "r_max_um": 480.66016042328255, "theta_cr_deg": 2.8647889756541165,
        "spraying_area_m2": 1.0666666666666666e-05, "release_angle_register_mean_deg": 12.956645121391393,
        "release_angle_register_min_deg": 0.9548412207142292, "release_angle_register_max_deg": 24.471238492884595}
