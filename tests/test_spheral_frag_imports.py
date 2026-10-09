"""The import boundary of spheral_frag's core (M1 plan, Global constraints and Task 1; review focus 4).

Every core module imports only the standard library and numpy at module level, so that the runner can import the
core under Spheral's Python; only `fe.py` names the finite-element package. The modules are listed from the source
tree, so modules added by later tasks are covered without editing this file. `runner/`, `m0/` and `container/` run
under Spheral's Python only and are not part of the core.

The runner's side of the boundary (Spheral M2 plan, Global constraints and Task 1), checked by reading the sources'
import statements (every one, at module level or inside a function), so that drama_env checks it without Spheral:
the core (and its command line) never imports `spheral_frag.runner`; `runner` modules import only the standard
library, numpy, scipy, mpi4py, Spheral's own modules (SPHERAL_MODULES), each other and the core's modules other than
the prepare side (RUNNER_BLOCKED_CORE)."""
import ast
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
RUNNER = "spheral_frag.runner"
RUNNER_THIRD_PARTY = {"numpy", "scipy", "mpi4py"}
# Spheral's Python modules a runner module may name: its compiled packages (Spheral, Spheral3d, SpheralRZ, ...,
# SpheralController) by prefix, and the helper modules listed here. Add one when a runner module needs it.
SPHERAL_PREFIX = "Spheral"
SPHERAL_MODULES = {"mpi", "PeanoHilbertDistributeNodes", "DistributeNodes", "distributeNodesGeneric",
                   "NodeGeneratorBase", "GenerateNodeDistribution2d", "GenerateNodeDistribution3d"}
# The prepare side and the command line: the runner reads the body stage's files, never these modules.
RUNNER_BLOCKED_CORE = {"spheral_frag.fe", "spheral_frag.prepare", "spheral_frag.__main__"}


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


# ---------------------------------------------------------------------------------------------------------------
# The runner's boundary (by reading the import statements)
# ---------------------------------------------------------------------------------------------------------------

def dotted_name(path):
    """spheral_frag/runner/nodes.py -> spheral_frag.runner.nodes; a package's __init__.py -> the package."""
    parts = os.path.relpath(path, REPO_ROOT)[:-3].split(os.sep)
    return ".".join(parts[:-1] if parts[-1] == "__init__" else parts)


def imported_modules(path):
    """Every module an import statement in the file names, relative imports resolved to dotted names."""
    tree = ast.parse(open(path, encoding="utf-8").read(), filename=path)
    name = dotted_name(path)
    package = name if os.path.basename(path) == "__init__.py" else name.rpartition(".")[0]
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.split(".")
                base = base[:len(base) - node.level + 1]
                prefix = ".".join(base)
                if node.module:
                    out.add(prefix + "." + node.module)
                else:
                    out.update(prefix + "." + a.name for a in node.names)
            else:
                out.add(node.module)
    return out


def runner_sources():
    root = os.path.join(PACKAGE, "runner")
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if not d.startswith((".", "__"))]
        for fn in sorted(filenames):
            if fn.endswith(".py"):
                yield os.path.join(dirpath, fn)


def core_sources():
    for p in all_python_sources():
        top = os.path.relpath(p, PACKAGE).split(os.sep)[0]
        if top not in SPHERAL_ONLY_DIRS:
            yield p


def runner_import_allowed(mod, core):
    top = mod.split(".")[0]
    if top == "spheral_frag":
        if mod == RUNNER or mod.startswith(RUNNER + "."):
            return True
        return mod not in RUNNER_BLOCKED_CORE and mod.rsplit(".", 1)[0] not in RUNNER_BLOCKED_CORE and (
            mod in core or mod.rsplit(".", 1)[0] in core)       # a module, or a name imported from one
    return (top in sys.stdlib_module_names or top in RUNNER_THIRD_PARTY or top.startswith(SPHERAL_PREFIX)
            or top in SPHERAL_MODULES)


def test_runner_sources_exist():
    names = {dotted_name(p) for p in runner_sources()}
    assert {RUNNER, RUNNER + ".nodes", RUNNER + ".checkpoint"} <= names


def test_core_never_imports_the_runner():
    offenders = {os.path.relpath(p, REPO_ROOT): sorted(m for m in imported_modules(p)
                                                       if m == RUNNER or m.startswith(RUNNER + "."))
                 for p in core_sources()}
    assert {k: v for k, v in offenders.items() if v} == {}


def test_runner_imports_only_what_spheral_python_has():
    core = set(core_modules())
    offenders = {os.path.relpath(p, REPO_ROOT): sorted(m for m in imported_modules(p)
                                                       if not runner_import_allowed(m, core))
                 for p in runner_sources()}
    assert {k: v for k, v in offenders.items() if v} == {}


def test_runner_rule_catches_violations():
    core = set(core_modules())
    for bad in ("reentry_model", "spheral_frag.fe", "spheral_frag.prepare", "spheral_frag.m0.bench", "pyvista",
                "gmsh", "spheral_frag.prepare.prepare_run"):
        assert not runner_import_allowed(bad, core), bad
    for ok in ("numpy", "scipy.spatial", "mpi4py", "Spheral3d", "SpheralRZ", "SpheralController", "mpi",
               "PeanoHilbertDistributeNodes", "spheral_frag.material", "spheral_frag.geometry.TetLocator",
               "spheral_frag.runner.nodes", "math", "__future__"):
        assert runner_import_allowed(ok, core), ok
