"""spheral_frag.material and fe: the material table against the finite-element material (spec §8.1, §8.3; M1 plan
Task 3, check "Material tables").

Thresholds are a priori (the plan's "to round-off"): |h_core - (h_FE - h_FE(300))| <= 1e-12 max|h| + 1e-9 J/kg,
f_l and c_p bitwise, temperature(enthalpy(T)) within 1e-9 K, h_full_liquid within 1e-12 J/kg.

Three sources of finite-element values, so that the tests run on every machine:
- the committed fixture `tests/fixtures/spheral_frag/material/` (AA7075_scheil from the reconstructed Step 3
  prototype, with 5,543 of the FE's own values; the generating command is in its JSON header and in
  make_material_table.py) -- always;
- the repository's own `reentry_model` (Step 2's AA7075_nomelt), 100,000 samples in a subprocess -- always;
- Step 3's package at $SPHERAL_FRAG_FE_PACKAGE (AA7075_scheil and AA7075_range), 100,000 samples in a subprocess --
  skipped with a message when unset or absent (prototype/ is git-ignored).
The FE values come from a subprocess because one process holds one reentry_model (fe.import_fe_package)."""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from helpers import FIXTURES, REPO_ROOT

from spheral_frag import contract, fe, frames, material, record

FIX = os.path.join(FIXTURES, "spheral_frag", "material")
SCHEIL = os.path.join(FIX, "material_AA7075_scheil.npz")
SCRIPT = os.path.join(FIXTURES, "spheral_frag", "make_material_table.py")
ENV = fe.FE_PACKAGE_ENV


def fe_package_or_skip():
    path = os.environ.get(ENV)
    if not path:
        pytest.skip("{} is not set: set it to Step 3's package directory (e.g. prototype/proto3) to compare against "
                    "AA7075_scheil and AA7075_range; the committed table fixture covers the core".format(ENV))
    path = os.path.join(REPO_ROOT, path) if not os.path.isabs(path) else path
    if not os.path.isfile(os.path.join(path, "reentry_model", "__init__.py")):
        pytest.skip("{}={} holds no reentry_model package (prototype/ is git-ignored)".format(ENV, path))
    return path


def generate(tmp_path, name, package=None, samples=100_000):
    """Run make_material_table.py in a subprocess; (table, FE samples, the script's own comparison)."""
    cmd = [sys.executable, SCRIPT, "--material", name, "--out", str(tmp_path), "--samples", str(samples)]
    if package:
        cmd += ["--fe-package", package]
    r = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr + r.stdout
    table = material.MaterialTable.load(os.path.join(tmp_path, "material_{}.npz".format(name)))
    with np.load(os.path.join(tmp_path, "material_{}_fe_samples.npz".format(name))) as z:
        s = {k: z[k] for k in z.files}
    return table, s, json.loads(r.stdout.strip().splitlines()[-1])


def assert_round_off(table, s, n_min):
    T = s["T"]
    assert len(T) >= n_min
    h = table.enthalpy(T)
    assert np.max(np.abs(h - s["h300"])) <= 1e-12 * np.max(np.abs(s["h300"])) + 1e-9
    assert np.array_equal(table.liquid_fraction(T), s["f_l"])                     # bitwise
    assert np.array_equal(table.cp(T), s["cp"])
    assert np.max(np.abs(table.temperature(h) - T)) <= 1e-9
    assert np.max(np.abs(table.temperature(s["h300"]) - s["T_from_h_fe"])) <= 1e-9
    hfl, want = table.h_full_liquid(), float(s["h_full_liquid"])
    assert hfl == want if not math.isfinite(want) else abs(hfl - want) <= 1e-12
    assert float(s["h_fe_300"]) == table.h_offset_300
    # off the grid: the samples hold the nodes +-1e-9 K and random temperatures between them
    assert np.count_nonzero(T != np.round(T)) > 0.9 * len(T) - 3 * len(table.h_T)


@pytest.fixture(scope="module")
def scheil():
    return material.MaterialTable.load(SCHEIL)


@pytest.fixture(scope="module")
def scheil_samples():
    with np.load(os.path.join(FIX, "material_AA7075_scheil_fe_samples.npz")) as z:
        return {k: z[k] for k in z.files}


# --- round-off against the finite element -----------------------------------------------------------------------

def test_committed_scheil_table_reproduces_the_fe_values(scheil, scheil_samples):
    assert scheil.name == "AA7075_scheil" and scheil.fl_kind == material.FL_INTERP
    assert_round_off(scheil, scheil_samples, 5000)
    assert np.array_equal(scheil.enthalpy(scheil_samples["T"]), scheil_samples["h300"])   # in fact bitwise


def test_main_package_nomelt_to_round_off(tmp_path):
    table, s, report = generate(tmp_path, "AA7075_nomelt")
    assert table.fl_kind == material.FL_NONE and not table.melts
    assert table.header["fe_package"]["file"] == os.path.join(REPO_ROOT, "reentry_model", "__init__.py")
    assert_round_off(table, s, 100_000)
    assert report["fl_bitwise"] and report["h_bitwise"]
    assert np.all(table.liquid_fraction(s["T"]) == 0.0) and math.isinf(table.h_full_liquid())
    assert np.all(np.isinf(table.viscosity([300.0, 900.0, 1400.0])))


@pytest.mark.parametrize("name, kind", [("AA7075_scheil", material.FL_INTERP), ("AA7075_range", material.FL_LINEAR)])
def test_step3_package_to_round_off(tmp_path, name, kind):
    pkg = fe_package_or_skip()
    table, s, report = generate(tmp_path, name, pkg)
    assert table.fl_kind == kind
    assert_round_off(table, s, 100_000)
    assert report["fl_bitwise"]


def test_step3_scheil_table_equals_the_committed_fixture(tmp_path, scheil):
    pkg = fe_package_or_skip()
    table, _, _ = generate(tmp_path, "AA7075_scheil", pkg, samples=10)
    with np.load(SCHEIL) as a, np.load(os.path.join(tmp_path, "material_AA7075_scheil.npz")) as b:
        assert a.files == b.files and all(np.array_equal(a[k], b[k]) for k in a.files)
    assert table.header["fe_package"]["sha256"] == scheil.header["fe_package"]["sha256"]


def test_one_finite_element_package_per_process(tmp_path):
    import reentry_model                                      # the repository's own, as the other tests load it
    fake = tmp_path / "other"
    (fake / "reentry_model").mkdir(parents=True)
    (fake / "reentry_model" / "__init__.py").write_text("")
    with pytest.raises(RuntimeError, match="already imported"):
        fe.import_fe_package(str(fake))
    with pytest.raises(FileNotFoundError):
        fe.import_fe_package(str(tmp_path / "nothing"))
    assert fe.import_fe_package(None) is reentry_model
    prov = fe.package_provenance(reentry_model)
    assert prov["path"] == REPO_ROOT and len(prov["sha256"]) == 64


# --- the Scheil material's anchors (decision 1) ------------------------------------------------------------------

def test_scheil_anchors(scheil):
    assert (scheil.T_solidus, scheil.T_liquidus, scheil.T_feed) == (748.0, 908.0, 910.0)
    assert (scheil.rho_solid, scheil.rho_liquid, scheil.mu_liquid, scheil.sigma_liquid) == (2813.0, 2400.0, 1.3e-3, 0.8)
    assert scheil.latent_heat == 390e3 and scheil.T_ref_fe == 293.0
    remainder = float(scheil.liquid_fraction(752.0))           # the eutectic remainder, top of the +-2 K ramp
    assert remainder == pytest.approx(0.0369, abs=5e-4)        # spec §2: 3.6 %
    assert float(scheil.liquid_fraction(750.0)) == pytest.approx(0.5 * remainder, rel=1e-12)
    assert float(scheil.liquid_fraction(748.0)) == 0.0 and float(scheil.liquid_fraction(908.0)) == 1.0
    assert scheil.T_half == pytest.approx(895.1, abs=0.05)
    assert float(scheil.liquid_fraction(scheil.T_half)) == pytest.approx(0.5, abs=1e-15)
    assert scheil.T_bridge == pytest.approx(906.4, abs=0.05)   # f_l = 0.9: the spec's 906.4 K
    assert abs(scheil.enthalpy(300.0)) == 0.0                   # the zero of energy at 300 K
    assert scheil.h_full_liquid() == pytest.approx(float(scheil.enthalpy(908.0)), abs=1e-9)


# --- viscosity (spec §8.3; decision 7) ---------------------------------------------------------------------------

def test_viscosity_rows(scheil):
    li = material.MaterialTable.li_viscosity
    assert float(li(0.1)) == pytest.approx(1.265, abs=0.005)
    assert float(li(0.5)) == pytest.approx(5.63, abs=0.01)
    assert float(li(0.3)) == pytest.approx(0.871 * math.exp(3.7311 * 0.3), rel=1e-15)
    assert float(scheil.viscosity(scheil.T_half)) == pytest.approx(5.63, abs=0.01)
    assert float(scheil.viscosity(scheil.T_bridge)) == pytest.approx(float(li(0.1)), rel=1e-12)
    assert np.all(scheil.viscosity([908.0, 950.0, 1500.0]) == 1.3e-3)
    assert np.all(np.isinf(scheil.viscosity([300.0, 748.0, 850.0, scheil.T_half - 1e-6])))
    for Tc in (scheil.T_bridge, 908.0):                       # continuous at both ends of the bridge
        lo, hi = scheil.viscosity([Tc - 1e-7, Tc + 1e-7])
        assert abs(math.log(hi / lo)) < 1e-5
    T = np.linspace(scheil.T_half, 908.0, 2001)
    eta = scheil.viscosity(T)
    assert np.all(np.diff(eta) <= 0.0)                        # thins monotonically to the liquid
    # the bridge is log-linear in T
    mid = 0.5 * (scheil.T_bridge + 908.0)
    assert float(scheil.viscosity(mid)) == pytest.approx(math.sqrt(float(li(0.1)) * 1.3e-3), rel=1e-12)


def test_viscosity_shear_rate_is_held_at_367(scheil):
    T = np.array([896.0, 900.0, 905.0, 907.0])
    at = scheil.viscosity(T, shear_rate=367.0)
    assert np.array_equal(scheil.viscosity(T, shear_rate=1e4), at)
    assert np.all(at < scheil.viscosity(T, shear_rate=0.0))
    f_s = 1.0 - scheil.liquid_fraction(900.0)
    assert float(scheil.viscosity(900.0, 100.0)) == pytest.approx((0.871 - 0.00849 * 100.0 ** 0.74924)
                                                                  * math.exp(3.7311 * f_s), rel=1e-14)
    assert scheil.viscosity(T, shear_rate=np.array([0.0, 10.0, 100.0, 1000.0])).shape == (4,)


# --- the 1 K table and the placeholders (decisions 2 and 3) ------------------------------------------------------

def test_one_kelvin_table(scheil):
    T = scheil.T_grid
    assert T[0] == 250.0 and T[-1] == 1500.0 and np.all(np.diff(T) == 1.0)
    assert np.array_equal(scheil.grid["h300"], scheil.enthalpy(T))
    assert np.array_equal(scheil.grid["f_l"], scheil.liquid_fraction(T))
    rho = scheil.rho_free(T)
    below = T < scheil.T_liquidus
    assert np.all(np.diff(rho[below]) < 0.0) and np.all(np.diff(rho) <= 0.0)
    assert np.all(rho[T >= scheil.T_liquidus] == 2400.0)       # decision 3: the FE liquid density
    alpha = 23.4e-6
    assert float(scheil.rho_free(700.0)) == pytest.approx(2813.0 / (1 + alpha * 407.0) ** 3, rel=1e-14)
    f = float(scheil.liquid_fraction(900.0))
    want = 1.0 / ((1 - f) * (1 + alpha * 607.0) ** 3 / 2813.0 + f / 2400.0)
    assert float(scheil.rho_free(900.0)) == pytest.approx(want, rel=1e-14)


def test_mechanical_placeholders(scheil):
    E0, nu = 71.7e9, 0.33
    assert np.all(scheil.young([300.0, 600.0, 748.0]) == E0)
    assert float(scheil.young(821.0)) == pytest.approx(E0 * (scheil.T_half - 821.0) / (scheil.T_half - 748.0), rel=1e-12)
    assert np.all(scheil.young([896.0, 908.0, 1200.0]) == 0.0) and np.all(scheil.shear_modulus([896.0, 1200.0]) == 0.0)
    assert np.all(scheil.poisson([300.0, 1000.0]) == nu)
    assert float(scheil.shear_modulus(500.0)) == pytest.approx(E0 / (2 * (1 + nu)), rel=1e-15)
    K = scheil.bulk_modulus(scheil.T_grid)
    assert np.all(K == K[0]) and K[0] == pytest.approx(E0 / (3 * (1 - 2 * nu)), rel=1e-15)   # held at the solidus
    assert set(scheil.provisional) == {"rho_free", "E", "K", "G"}
    mech = scheil.header["mechanical"]["entries"]
    assert any(e["status"] == "provisional" for e in mech.values())
    assert mech["volume_change_on_melting"]["used"] is False
    assert scheil.implied_melting_volume_change() == pytest.approx(0.1229, abs=1e-4)


def test_provisional_follows_the_input_file(scheil):
    mech = json.loads(json.dumps(scheil.header["mechanical"]))
    for e in mech["entries"].values():
        e["status"] = "approved"
    T = scheil.T_grid
    _, prov = material.mechanical_columns(T, scheil.liquid_fraction(T), scheil.T_solidus, scheil.T_half,
                                          scheil.rho_solid, scheil.rho_liquid, mech)
    assert prov == ()


# --- storage and interfaces ---------------------------------------------------------------------------------------

def test_write_load_round_trip_and_integrity(tmp_path, scheil):
    with np.load(SCHEIL) as z:
        arrays = {k: z[k] for k in z.files}
    out = tmp_path / "prep" / "material_table.npz"
    material.write_table(out, arrays, {"fe_material": "AA7075_scheil", "provisional": ["E"]})
    t = material.MaterialTable.load(tmp_path / "prep")          # a directory resolves to material_table.npz
    assert t.provisional == ("E",) and np.array_equal(t.h_nodes, scheil.h_nodes)
    with np.load(out, allow_pickle=False) as z:                 # plain numeric npz (numpy 1.26 reads it)
        assert all(z[k].dtype.kind in "fi" for k in z.files)
    with open(out, "ab") as fh:
        fh.write(b"\0")
    with pytest.raises(ValueError, match="sha256"):
        material.MaterialTable.load(out)
    with pytest.raises(ValueError, match="lacks"):
        material.write_table(tmp_path / "x.npz", {"h_T": arrays["h_T"]}, {})


def test_satisfies_the_record_protocol(scheil):
    assert isinstance(scheil, record.RecordMaterial)
    assert float(scheil.viscosity(950.0, shear_rate=0.0)) == 1.3e-3
    assert isinstance(scheil.h_full_liquid(), float)


def test_consistency_check_with_the_table(scheil):
    T = np.array([300.0, 750.0, 895.0, 906.0, 930.0])
    fr = contract.Frame(k=0, time_s=0.0, points=np.zeros((5, 3)), tets=np.zeros((0, 4), np.int64),
                        faces=np.zeros((0, 3), np.int64), node={"T": T, "f_l": scheil.liquid_fraction(T)},
                        tet={}, patch={})
    c = frames.consistency_checks(fr, scheil)
    assert c["fl_bitwise"] is True and c["max_abs_dfl"] == 0.0
    fr.node["f_l"] = fr.node["f_l"] + np.array([0, 0, 1e-16, 0, 0])
    assert frames.consistency_checks(fr, scheil)["fl_bitwise"] is False


def test_interp_reproduces_either_kind_of_numpy(scheil):
    """material.interp: np.interp as a fusing or a non-fusing numpy computes it (the probe is a measured case)."""
    from fractions import Fraction
    x, xp, fp, fused = material._FMA_PROBE
    slope, d = (fp[1] - fp[0]) / (xp[1] - xp[0]), x - xp[0]
    assert float(Fraction(slope) * Fraction(d) + Fraction(fp[0])) == fused       # one rounding: the fma
    assert float(material._interp_emulated(x, np.array(xp), np.array(fp), False)) == slope * d + fp[0] != fused
    assert scheil.interp_fma == scheil.header["interp_fma"]
    T = np.random.default_rng(3).uniform(250.0, 1500.0, 20_000)
    T = np.concatenate([T, scheil.fl_T, [100.0, 2000.0, scheil.fl_T[0], scheil.fl_T[-1]]])
    here = material.interp_fuses()
    assert np.array_equal(material.interp(T, scheil.fl_T, scheil.fl, here), np.interp(T, scheil.fl_T, scheil.fl))
    assert np.array_equal(material._interp_emulated(T, scheil.fl_T, scheil.fl, False),
                          material.interp(T, scheil.fl_T, scheil.fl, False))
    if material.QUAD:                                            # Linux aarch64: the fused form in quad precision
        assert float(material._interp_emulated(x, np.array(xp), np.array(fp), True)) == fused
    else:
        assert scheil.interp_exact == (scheil.interp_fma == here or not scheil.interp_fma)


# --- fe.py's two helpers for Task 9 (the repository's package; the synthetic sphere run) ---------------------------

SPHERE_RUN = os.path.join(FIXTURES, "spheral_frag", "sphere_run")


def test_air_temperature_is_the_runs_atmosphere():
    from reentry_model import atmosphere
    run = frames.read_fe_run(SPHERE_RUN)
    T_air = fe.air_temperature(run.run_json, run.history)
    us76 = atmosphere.US76TableAtmosphere()
    want = [us76.state(0.0, h * 1e3, 0.0, 0.0).T for h in run.history["altitude_km"]]
    assert run.run_json["settings"]["atmosphere"] == "us76" and np.array_equal(T_air, want)


def test_fe_heat_content():
    run = frames.read_fe_run(SPHERE_RUN)
    mat = fe.fe_material("AA7075_nomelt")
    assert fe.fe_heat_content(frames.read_frame(run, 0), mat) == 0.0        # 300 K uniform: the zero of energy
    fr = frames.read_frame(run, 1)
    p, tets = fr.points, fr.tets
    h = mat.enthalpy(fr.node["T"]) - mat.enthalpy(300.0)
    want = 0.0
    for e, tet in enumerate(tets):
        a, b, c, d = p[tet]
        V = abs(np.dot(np.cross(b - a, c - a), d - a)) / 6.0
        want += fr.tet["phi"][e] * mat.rho * V * h[tet].sum() / 4.0
    assert fe.fe_heat_content(fr, mat) == pytest.approx(want, rel=1e-12)
