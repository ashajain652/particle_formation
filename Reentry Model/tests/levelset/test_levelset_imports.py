"""The level set uses only the shared physics (Asha, 2026-11-11, level-set stage 3 Task 1): every reentry_model import in
reentry_model/levelset is the shared package, the level-set package itself or a package-root constant, never one of the
skins' modules, and the package defines no physical constant of its own."""
import ast
import glob
import os

from helpers import PROJECT_ROOT

LEVELSET_DIR = os.path.join(PROJECT_ROOT, "reentry_model", "levelset")
ROOT_CONSTANTS = ("__version__", "PROJECT_ROOT", "REPO_ROOT", "DATA_DIR", "OUTPUT_DIR", "DRAMA_DIR", "DRAMA_DATA_DIR")
SKINS = ("skin", "skinbody", "columns", "film", "spreading", "movemesh", "mesh", "restart", "coupled", "cli", "body",
         "thermal")


def _paths():
    paths = sorted(glob.glob(os.path.join(LEVELSET_DIR, "*.py")))
    assert len(paths) >= 10
    return paths


def _allowed(dotted):
    return any(dotted == p or dotted.startswith(p + ".") for p in ("reentry_model.shared", "reentry_model.levelset"))


def _escapes(node, levelset_names):
    """True when the import statement `node`, in a level-set module, reaches a reentry_model module outside the shared
    package and the level-set package (package-root constants excepted)."""
    if isinstance(node, ast.Import):
        return any(a.name.split(".")[0] == "reentry_model" and not _allowed(a.name) for a in node.names)
    if node.level == 0:
        mod = node.module or ""
        if mod.split(".")[0] != "reentry_model":
            return False
        if mod == "reentry_model":
            return not all(a.name in ROOT_CONSTANTS or a.name == "shared" for a in node.names)
        return not _allowed(mod)
    if node.level == 1:                   # inside levelset: `from . import x`, `from .x import y`
        if node.module is None:
            return not all(a.name in levelset_names for a in node.names)
        return node.module.split(".")[0] not in levelset_names
    if node.level == 2:                   # `from .. import CONSTANT` or `from ..shared[.x] import y`
        if node.module is None:
            return not all(a.name in ROOT_CONSTANTS or a.name == "shared" for a in node.names)
        return node.module.split(".")[0] != "shared"
    return True


def test_levelset_imports_only_shared_physics():
    paths = _paths()
    names = {os.path.basename(p)[:-3] for p in paths}
    bad, scanned = [], 0
    for path in paths:
        with open(path) as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                scanned += 1
                if _escapes(node, names):
                    bad.append((os.path.basename(path), node.lineno))
    assert scanned > 0
    assert bad == []
    # the scan itself flags each of the skins' modules, by every spelling, and passes the allowed ones
    for skin in SKINS:
        for src in ("from reentry_model.{0} import x", "import reentry_model.{0}", "from reentry_model import {0}",
                    "from .. import {0}", "from ..{0} import x", "from .{0} import x"):
            if src.startswith("from .{") and skin in names:
                continue                  # `.skin` inside levelset/ is the level set's own skin.py (stage 3 Task 2)
            assert _escapes(ast.parse(src.format(skin)).body[0], names), src.format(skin)
    assert not _escapes(ast.parse("from .skin import Skins").body[0], names)
    for src in ("from reentry_model.shared.thermal import SIGMA_SB", "from reentry_model.shared import spray",
                "from ..shared.material import Material", "from .. import DATA_DIR", "from .cutcell import patches",
                "from reentry_model.levelset.transfer import locate", "import numpy"):
        assert not _escapes(ast.parse(src).body[0], names), src


SHARED_CONSTANTS = (5.670374419e-8, 895.1, 910.0, 293.0)   # SIGMA_SB, AA7075's T_rigid (rounded), its T_feed, T_REF


def test_levelset_defines_no_physical_constants():
    """No module-level assignment in levelset/ gives a name the value of a shared physical constant."""
    found = []
    for path in _paths():
        with open(path) as fh:
            tree = ast.parse(fh.read())
        for node in tree.body:
            if isinstance(node, (ast.Assign, ast.AnnAssign)) and isinstance(node.value, ast.Constant) \
                    and isinstance(node.value.value, (int, float)) and not isinstance(node.value.value, bool):
                if float(node.value.value) in SHARED_CONSTANTS:
                    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                    found.append((os.path.basename(path), [ast.unparse(t) for t in targets]))
    assert found == []
