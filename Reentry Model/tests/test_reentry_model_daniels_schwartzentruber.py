"""Daniels & Schwartzentruber (2026) sphere pressure and shear, ported to the shared package: it must reproduce the
reference implementation's distributions and parameters (fixtures in tests/fixtures/gas_loads/)."""
import csv
import os

import numpy as np
import pytest

from reentry_model.shared import daniels_schwartzentruber as ds

# The fixtures are printed to six significant figures (relative rounding up to 5e-6), so they are matched to 1e-5;
# test_matches_reference_module checks the port against the reference implementation itself to rounding error.
FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "gas_loads")
CASES = [("0.002", 1.007e-4), ("0.01", 2.014e-5), ("0.05", 3.987e-6), ("0.25", 8.057e-7)]


def _cond(rho):
    return ds.condition(rho, 200.0, 7108.0, 0.3048, 1500.0)


def _rows(name):
    with open(os.path.join(FIX, name)) as f:
        return list(csv.DictReader(f))


@pytest.mark.parametrize("label,rho", CASES)
def test_reproduces_reference_distributions(label, rho):
    rows = [r for r in _rows("ds_model.csv") if r["case_kn"] == label]
    assert len(rows) == 181
    theta = np.array([float(r["theta_deg"]) for r in rows])
    cp, cf = ds.coefficients(np.radians(theta), _cond(rho))
    np.testing.assert_allclose(cp, [float(r["Cp"]) for r in rows], rtol=1e-5, atol=1e-9)
    np.testing.assert_allclose(cf, [float(r["Ctau"]) for r in rows], rtol=1e-5, atol=1e-9)


@pytest.mark.parametrize("label,rho", CASES)
def test_reproduces_reference_parameters(label, rho):
    row = next(r for r in _rows("ds_params.csv") if r["case_kn"] == label)
    P = ds.parameters(_cond(rho))
    checked = 0
    for k, v in row.items():
        if k in ("case_kn", "Cd_HB_DSMC", "Cd_model", "Cd_continuum_only"):
            continue
        if v == "":
            continue
        if v == "nan":
            continue
        assert P[k] == pytest.approx(float(v), rel=1e-5), k
        checked += 1
    assert checked >= 28


@pytest.mark.parametrize("rho,cd", [(c[1], d) for c, d in zip(CASES, [1.0702, 1.236, 1.534, 1.796])])
def test_drag(rho, cd):
    assert ds.drag_coefficient(_cond(rho)) == pytest.approx(cd, rel=2e-3)


def test_separation_angle_capped():
    c = ds.condition(8.057e-7, 200.0, 7108.0, 0.3048, 1500.0)
    assert ds.parameters(c)["theta_s"] <= 180.0
    assert ds.parameters(c)["theta_s"] == 180.0


@pytest.mark.parametrize("label,rho", CASES)
def test_matches_reference_module(label, rho):
    import importlib.util
    path = os.path.join(os.path.dirname(FIX), "..", "..", "analysis", "gas_loads", "reference", "ds_model.py")
    if not os.path.exists(path):
        pytest.skip("reference implementation not present")
    spec = importlib.util.spec_from_file_location("ds_model_reference", path)
    ref = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ref)
    P, _ = ref.hb_parameters(rho)
    theta = np.arange(0.0, 181.0)
    cf_ref, cp_ref = ref.distribution(theta, P)
    cp, cf = ds.coefficients(np.radians(theta), _cond(rho))
    np.testing.assert_allclose(cp, cp_ref, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(cf, cf_ref, rtol=1e-12, atol=1e-14)
