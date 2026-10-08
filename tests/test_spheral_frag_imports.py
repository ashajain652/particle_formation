"""The import boundary of spheral_frag's core (M1 plan, Global constraints and Task 1; review focus 4).

Every core module imports only the standard library and numpy at module level, so that the runner can import the
core under Spheral's Python; only `fe.py` names the finite-element package. The modules are listed from the source
tree, so modules added by later tasks are covered without editing this file. `runner/`, `m0/` and `container/` run
under Spheral's Python only and are not part of the core."""
import os
import re
import subprocess
import sys

import pytest

from helpers import REPO_ROOT

PACKAGE = os.path.join(REPO_ROOT, "spheral_frag")
SPHERAL_ONLY_DIRS = {"runner", "m0", "container"}
BLOCKED = ("scipy", "pyvista", "vtk", "vtkmodules", "gmsh", "cantera", "pymsis", "reentry_model")
FE_PACKAGE = "reentry_model"


def core_modules():
    """Dotted names of every core module (packages under spheral_frag/ with an __init__.py, Spheral-only dirs out)."""
    mods = []
    for dirpath, dirnames, filenames in os.walk(PACKAGE):
        rel = os.path.relpath(dirpath, REPO_ROOT)
        dirnames[:] = sorted(d for d in dirnames if d not in SPHERAL_ONLY_DIRS and not d.startswith((".", "__"))
                             and os.path.exists(os.path.join(dirpath, d, "__init__.py")))
        for fn in sorted(filenames):
            if not fn.endswith(".py") or fn == "__main__.py":
                continue
            parts = rel.split(os.sep) + ([] if fn == "__init__.py" else [fn[:-3]])
            mods.append(".".join(parts))
    return mods


def all_python_sources():
    for dirpath, dirnames, filenames in os.walk(PACKAGE):
        dirnames[:] = [d for d in dirnames if not d.startswith((".", "__"))]
        for fn in filenames:
            if fn.endswith(".py"):
                yield os.path.join(dirpath, fn)


PROBE = r"""
import importlib, importlib.abc, sys
BLOCKED = {blocked!r}
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in BLOCKED:
            raise ImportError("blocked by the import-boundary test: " + name)
        return None
sys.meta_path.insert(0, Block())
sys.path.insert(0, {root!r})
before = set(sys.modules)
importlib.import_module({module!r})
allowed = set(sys.stdlib_module_names) | {{"numpy", "spheral_frag"}}
extra = sorted({{m.split(".")[0] for m in set(sys.modules) - before}} - allowed)
extra = [m for m in extra if not m.startswith("_")]
if extra:
    raise SystemExit("imported outside the standard library and numpy: " + ", ".join(extra))
"""


def test_core_modules_listed():
    mods = core_modules()
    assert {"spheral_frag", "spheral_frag.contract", "spheral_frag.naming"} <= set(mods)
    assert "spheral_frag.__main__" not in mods and not any(".m0" in m or ".runner" in m for m in mods)


@pytest.mark.parametrize("module", core_modules())
def test_imports_with_numpy_and_the_standard_library_only(module):
    code = PROBE.format(blocked=BLOCKED, root=REPO_ROOT, module=module)
    r = subprocess.run([sys.executable, "-I", "-c", code], cwd=REPO_ROOT, capture_output=True, text=True,
                       timeout=60)
    assert r.returncode == 0, r.stderr + r.stdout


def test_only_fe_names_the_finite_element_package():
    pattern = re.compile(r"\b" + FE_PACKAGE + r"\b")
    offenders = [os.path.relpath(p, REPO_ROOT) for p in all_python_sources()
                 if os.path.basename(p) != "fe.py" and pattern.search(open(p, encoding="utf-8").read())]
    assert offenders == []
    fe = [p for p in all_python_sources() if os.path.basename(p) == "fe.py"]
    assert all(os.path.dirname(p) == PACKAGE for p in fe), fe     # spheral_frag/fe.py, not another fe.py


def test_main_stub_exits_2():
    r = subprocess.run([sys.executable, "-m", "spheral_frag"], cwd=REPO_ROOT, capture_output=True, text=True,
                       timeout=60)
    assert r.returncode == 2 and "usage" in r.stderr
