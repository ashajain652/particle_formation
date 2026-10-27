import numpy as np
import pytest

from reentry_model.wake_column.ambient import AmbientColumn


@pytest.mark.parametrize("src", ["waccm", "msis_afgl"])
def test_columns_cover_45_to_100_km_and_close(src):
    for tod in ("day", "night"):
        col = AmbientColumn.load(src, tod)
        for z in (45.0, 60.0, 76.0, 90.0, 100.0):
            a = col.at(z)
            assert 130 < a.T_K < 320
            assert sum(a.n.values()) == pytest.approx(a.n_total_m3, rel=0.02)
            assert all(v >= 0 for v in a.n.values())


def test_day_night_atomic_oxygen_contrast_below_80_km():
    day, night = AmbientColumn.load("waccm", "day").at(70.0), AmbientColumn.load("waccm", "night").at(70.0)
    assert day.n["O"] / night.n["O"] > 100
    assert day.n["H2O"] == pytest.approx(night.n["H2O"], rel=0.3)


def test_waccm_and_msis_agree_on_total_density_within_factor_2():
    for z in (60.0, 76.0, 90.0):
        a, b = AmbientColumn.load("waccm", "night").at(z), AmbientColumn.load("msis_afgl", "night").at(z)
        assert 0.5 < a.n_total_m3 / b.n_total_m3 < 2.0


def test_outside_table_raises_with_altitude_in_message():   # Review Focus 2
    with pytest.raises(ValueError, match="44.0 km"):
        AmbientColumn.load("waccm", "day").at(44.0)


# ---- continuity of the tables --------------------------------------------------------------------------------------
# The largest change between adjacent 1 km levels in the WACCM tables is 0.64 decades (night atomic oxygen at 80 km,
# where WACCM's night depletion ends); the MSIS + AFGL tables stay below 0.5. 0.7 decades per km therefore passes
# every real profile with a small margin, but fails the splice cliffs a plain WACCM fill into the MSIS table produced
# (night O: 9 decades at 50-51 km, from 3e4 to 6e13 m^-3; H: 1.4 decades at 75-76 km), where log-linear interpolation would be meaningless.
MAX_LOG10_JUMP_PER_KM = 0.7


@pytest.mark.parametrize("src", ["waccm", "msis_afgl"])
@pytest.mark.parametrize("tod", ["day", "night"])
def test_no_jump_between_adjacent_levels(src, tod):
    from reentry_model.wake_column.species import AMBIENT_SPECIES
    col = AmbientColumn.load(src, tod)
    series = {"n_total": col.n_total_m3, "p": col.p_Pa, **{s: col.n[s] for s in AMBIENT_SPECIES}}
    for name, y in series.items():
        jump = np.abs(np.diff(np.log10(y)))
        assert jump.max() < MAX_LOG10_JUMP_PER_KM, (name, col.z_km[jump.argmax()], jump.max())


# ---- the WACCM conversion code, on a synthetic stand-in for the netCDF variables ------------------------------------
class _Var:
    """Stands in for a scipy netcdf variable: indexable, and a scalar has getValue()."""
    def __init__(self, a):
        self.a = np.asarray(a, dtype=float)

    def __getitem__(self, idx):
        return self.a[idx]

    def getValue(self):
        return float(self.a)


def _isothermal_file(T=250.0, ps=1.0e5, p0=1.0e5, zint_geop_m=None, q=None):
    """Isothermal atmosphere, pure sigma levels (hyai = 0, hybi = p/PS), constant mass mixing ratios; 1 lat x 1 lon.

    Interfaces at the given geopotential heights Z (top first in the file); the hypsometric integral of an isothermal
    layer is exact, so Z_i = H ln(PS/p_i) with H = Rd T/g.
    """
    from reentry_model.wake_column.analysis import extract_waccm_columns as E
    H = E.R_D * T / E.G0
    z_int = np.asarray(zint_geop_m if zint_geop_m is not None else np.arange(0.0, 130e3 + 1, 2e3))[::-1]   # top first
    hybi = np.exp(-z_int / H) * ps / p0
    hyai = np.zeros_like(hybi)
    hybm, hyam = 0.5 * (hybi[:-1] + hybi[1:]), np.zeros(len(hybi) - 1)
    nlev = len(hybm)
    v = {"hyai": _Var(hyai), "hybi": _Var(hybi), "hyam": _Var(hyam), "hybm": _Var(hybm), "P0": _Var(p0),
         "PS": _Var(np.full((1, 1, 1), ps)), "T": _Var(np.full((1, nlev, 1, 1), T))}
    q = q or {"O2": 0.2314, "O": 1e-4, "O3": 2e-6, "H": 1e-9, "H2": 4e-8, "Q": 5e-6, "CO": 1e-6, "CO2": 5e-4}
    for name, val in q.items():
        v[name] = _Var(np.full((1, nlev, 1, 1), val))
    return v, H, z_int, q, hybm * p0


def test_native_column_geometric_altitude_densities_and_n2_remainder():
    from reentry_model.wake_column.analysis import extract_waccm_columns as E
    from reentry_model.wake_column.species import load_species
    mass = {s: sp.mass_amu for s, sp in load_species().items()}
    v, H, z_int, q, p_mid = _isothermal_file()
    col = E.native_column(v, 0, 0, mass)
    # altitude: layer-midpoint geopotential height, converted to geometric altitude z = R Z/(R - Z)
    Z = 0.5 * (z_int[:-1] + z_int[1:])
    z_expected = E.R_EARTH * Z / (E.R_EARTH - Z) / 1e3
    assert col["z_km"] == pytest.approx(z_expected, rel=1e-12)
    assert col["z_km"][0] > Z[0] / 1e3 + 1.0            # near 130 km the geometric altitude is some km above Z
    assert col["p_Pa"] == pytest.approx(p_mid, rel=1e-12)
    n_tot = p_mid / (E.K_B * 250.0)
    assert col["n_total_m3"] == pytest.approx(n_tot, rel=1e-12)
    for s, var in E.WACCM_VAR.items():
        assert col[s] == pytest.approx(q[var] * 28.96 / mass[s] * n_tot, rel=1e-12)
    others = sum(col[s] for s in E.WACCM_VAR)
    assert col["N2"] == pytest.approx(n_tot - others, rel=1e-12)
    assert np.all(col["N2"] > 0)


def test_native_column_starts_from_surface_geopotential_when_phis_present():
    from reentry_model.wake_column.analysis import extract_waccm_columns as E
    from reentry_model.wake_column.species import load_species
    mass = {s: sp.mass_amu for s, sp in load_species().items()}
    v, H, z_int, q, p_mid = _isothermal_file()
    base = E.native_column(v, 0, 0, mass)["z_km"]
    v["PHIS"] = _Var(np.full((1, 1, 1), 1000.0 * E.G0))          # a 1000 m surface
    lifted = E.native_column(v, 0, 0, mass)["z_km"]
    assert np.all(lifted > base + 0.99)                             # 1 km higher geopotential, a little more geometric


def test_to_grid_interpolates_temperature_linearly_and_densities_log_linearly():
    from reentry_model.wake_column.analysis import extract_waccm_columns as E
    zn = np.array([43.3, 47.1, 52.9, 61.0, 77.7, 90.4, 103.2, 120.0])       # irregular native altitudes, km
    col = {"z_km": zn, "T_K": 200.0 + 1.5 * zn, "p_Pa": 1e5 * np.exp(-zn / 7.0), "n_total_m3": 1e25 * np.exp(-zn / 6.0),
           "O": 1e12 * np.exp(0.2 * zn)}
    grid = np.arange(45.0, 100.0 + 0.5, 1.0)
    out = E.to_grid(col, grid)
    assert out["z_km"] == pytest.approx(grid)
    assert out["T_K"] == pytest.approx(200.0 + 1.5 * grid, rel=1e-12)           # linear in z, so exact
    assert out["n_total_m3"] == pytest.approx(1e25 * np.exp(-grid / 6.0), rel=1e-12)   # exponential in z, so exact
    assert out["O"] == pytest.approx(1e12 * np.exp(0.2 * grid), rel=1e-12)
    assert out["p_Pa"] == pytest.approx(1e5 * np.exp(-grid / 7.0), rel=1e-12)
    with pytest.raises(ValueError, match="native levels cover"):
        E.to_grid({**col, "z_km": zn[1:-1] + 0.0}, np.arange(40.0, 100.5, 1.0))


def test_native_column_then_to_grid_recovers_the_isothermal_profile():
    """End to end on the isothermal stand-in: the gridded density follows exp(-Z/H) with Z the geopotential height of the
    geometric grid altitude. The only departures are the layer-midpoint pressure (a constant factor cosh(d/2) for a
    layer of thickness d = 2 km / H, 0.9 %) and the curvature of log-linear interpolation (far below 0.1 %)."""
    from reentry_model.wake_column.analysis import extract_waccm_columns as E
    from reentry_model.wake_column.species import load_species
    mass = {s: sp.mass_amu for s, sp in load_species().items()}
    v, H, z_int, q, p_mid = _isothermal_file()
    out = E.to_grid(E.native_column(v, 0, 0, mass), E.Z_GRID_KM)
    z_m = E.Z_GRID_KM * 1e3
    Z = E.R_EARTH * z_m / (E.R_EARTH + z_m)                   # geopotential height of each geometric grid altitude
    d = 2e3 / H
    p_expected = 1e5 * np.exp(-Z / H) * np.cosh(d / 2)
    assert out["p_Pa"] == pytest.approx(p_expected, rel=2e-3)
    assert out["T_K"] == pytest.approx(250.0)
    assert out["n_total_m3"] == pytest.approx(p_expected / (E.K_B * 250.0), rel=2e-3)
    assert out["O"] / out["n_total_m3"] == pytest.approx(q["O"] * 28.96 / mass["O"], rel=1e-3)
