"""The skyscrapers subpackage: the old module paths are the moved modules themselves (same object, same pickles),
and shared/ and levelset/ never import it. Later tasks extend MOVED."""
import ast
import glob
import importlib
import os
import pickle
import subprocess

import pytest

from helpers import PROJECT_ROOT, PY

MOVED = ["body", "coupled", "columns", "mesh", "movemesh", "restart", "skin", "skinbody", "spreading"]


@pytest.mark.parametrize("name", MOVED)
def test_old_path_is_the_moved_module(name):
    old = importlib.import_module("reentry_model." + name)
    new = importlib.import_module("reentry_model.skyscrapers." + name)
    assert old is new
    assert getattr(importlib.import_module("reentry_model"), name) is new


@pytest.mark.parametrize("name", MOVED)
def test_import_order_does_not_matter(name):
    code = ("import reentry_model.skyscrapers.{n} as new\n"
            "import reentry_model.{n} as old\n"
            "from reentry_model import {n} as attr\n"
            "print(old is new and attr is new)\n").format(n=name)
    out = subprocess.run([PY, "-c", code], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout.strip()
    assert out == "True"


def test_old_pickle_path_loads_the_moved_class():
    """A state pickled before the move names the class as reentry_model.body.MeltSettings; it must still load."""
    from reentry_model.skyscrapers.body import MeltSettings
    data = pickle.dumps(MeltSettings(runoff=False), protocol=0)
    assert b"reentry_model.skyscrapers.body" in data
    old = data.replace(b"reentry_model.skyscrapers.", b"reentry_model.")
    assert b"reentry_model.skyscrapers" not in old
    loaded = pickle.loads(old)
    assert type(loaded) is MeltSettings
    assert loaded.runoff is False


def _resolved_imports(path, package):
    """Dotted names every import in the file resolves to; `package` is the file's own package, dotted."""
    with open(path) as fh:
        tree = ast.parse(fh.read())
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.split(".")[:len(package.split(".")) - (node.level - 1)]
                base = ".".join(base + ([node.module] if node.module else []))
            else:
                base = node.module or ""
            found.append(base)
            found += [base + "." + a.name for a in node.names]
    return found


def test_shared_and_levelset_do_not_import_skyscrapers():
    scanned = 0
    bad = []
    for sub in ("shared", "levelset"):
        for path in sorted(glob.glob(os.path.join(PROJECT_ROOT, "reentry_model", sub, "*.py"))):
            for dotted in _resolved_imports(path, "reentry_model." + sub):
                scanned += 1
                if dotted == "reentry_model.skyscrapers" or dotted.startswith("reentry_model.skyscrapers."):
                    bad.append((os.path.relpath(path, PROJECT_ROOT), dotted))
    assert scanned > 0
    assert bad == []


def test_the_scan_resolves_relative_imports(tmp_path):
    f = tmp_path / "x.py"
    f.write_text("from .. import skyscrapers\nfrom ..skyscrapers.body import B\nfrom . import y\n")
    got = _resolved_imports(str(f), "reentry_model.shared")
    assert "reentry_model.skyscrapers" in got and "reentry_model.skyscrapers.body" in got
    assert "reentry_model.shared.y" in got


def test_film_is_split_into_lubrication_and_runoff():
    """Task 3: the old film path is the element model's film module; the lubrication formulas live in shared."""
    old = importlib.import_module("reentry_model.film")
    new = importlib.import_module("reentry_model.skyscrapers.elements.film")
    lub = importlib.import_module("reentry_model.shared.lubrication")
    assert old is new
    assert old.lubrication is lub.lubrication
    assert old.deep_flux is lub.deep_flux
    assert old.Runoff.__module__ == "reentry_model.skyscrapers.elements.film"
    assert lub.lubrication.__module__ == "reentry_model.shared.lubrication"


def test_shared_lubrication_imports_only_numpy():
    path = os.path.join(PROJECT_ROOT, "reentry_model", "shared", "lubrication.py")
    assert os.path.exists(path)
    with open(path) as fh:
        tree = ast.parse(fh.read())
    mods = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            mods.append(("." * node.level) + (node.module or ""))
    assert mods == ["numpy"]
