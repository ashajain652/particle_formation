"""Spheral's Python evaluates the material table without reentry_model (M1 plan Task 3; spec §8.1).

Through the M0 launcher, Spheral's Python (numpy 1.26) imports `spheral_frag.material`, loads tables written by
drama_env -- AA7075_nomelt from the repository's own package, and the committed AA7075_scheil fixture -- and evaluates
enthalpy, its inverse, the liquid fraction, c_p, the free density and the viscosity at 10,000 temperatures; the
results equal drama_env's bitwise. Without `material.interp` they would not: Spheral's numpy does not fuse
np.interp's last multiply-add and drama_env's does (6 of these 10,000 liquid fractions and 38 c_p values differed by
one unit in the last place before it, measured 2026-10-07)."""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from conftest import spheral_launcher
from helpers import FIXTURES, REPO_ROOT

from spheral_frag import material

SCRIPT = os.path.join(FIXTURES, "spheral_frag", "make_material_table.py")
CODE = r"""
import sys, numpy as np
from spheral_frag import material
assert "reentry_model" not in sys.modules
for name, path in (("nomelt", {nomelt!r}), ("scheil", {scheil!r})):
    t = material.MaterialTable.load(path)
    T = np.load({temps!r})["T"]
    h = t.enthalpy(T)
    np.savez({out!r}.format(name), h=h, T_back=t.temperature(h), f_l=t.liquid_fraction(T), mu=t.viscosity(T),
             cp=t.cp(T), rho=t.rho_free(T), exact=np.array(t.interp_exact), numpy=np.array([int(x) for x in np.__version__.split(".")[:2]]))
assert "reentry_model" not in sys.modules
"""


@pytest.fixture
def workdir():
    base = os.path.join(REPO_ROOT, "spheral_output", "tests")       # inside the repository: the container sees it
    os.makedirs(base, exist_ok=True)
    d = tempfile.mkdtemp(prefix="material_", dir=base)
    yield d
    shutil.rmtree(d, ignore_errors=True)


@pytest.mark.spheral
def test_spheral_python_reproduces_the_tables_bitwise(workdir):
    r = subprocess.run([sys.executable, SCRIPT, "--material", "AA7075_nomelt", "--out", workdir, "--samples", "10"],
                       cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr
    rel = lambda p: os.path.relpath(p, REPO_ROOT)
    tables = {"nomelt": os.path.join(workdir, "material_AA7075_nomelt.npz"),
              "scheil": os.path.join(FIXTURES, "spheral_frag", "material", "material_AA7075_scheil.npz")}
    T = np.random.default_rng(7).uniform(250.0, 1500.0, 10_000)
    np.savez(os.path.join(workdir, "temps.npz"), T=T)
    code = CODE.format(nomelt=rel(tables["nomelt"]), scheil=rel(tables["scheil"]), temps=rel(os.path.join(workdir, "temps.npz")),
                       out=rel(os.path.join(workdir, "spheral_{}.npz")))
    r = subprocess.run([spheral_launcher(), "-c", code], cwd=REPO_ROOT, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr + r.stdout
    for name, path in tables.items():
        t = material.MaterialTable.load(path)
        with np.load(os.path.join(workdir, "spheral_{}.npz".format(name))) as z:
            assert tuple(z["numpy"]) < (2, 0)                        # really the other numpy
            h = t.enthalpy(T)
            assert np.array_equal(z["h"], h)
            assert np.array_equal(z["T_back"], t.temperature(h))
            assert np.array_equal(z["f_l"], t.liquid_fraction(T))
            assert np.array_equal(z["mu"], t.viscosity(T))
            assert np.array_equal(z["cp"], t.cp(T)) and np.array_equal(z["rho"], t.rho_free(T))
            assert bool(z["exact"])
