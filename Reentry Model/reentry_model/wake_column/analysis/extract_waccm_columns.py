"""Extract the WACCM day and night ambient columns at 40.7 N into data/wake_column/ambient_waccm_40.7N.csv.

Usage:  python -m reentry_model.wake_column.analysis.extract_waccm_columns <path to the WACCM cam.i .nc file> [out.csv]

The 1.6 GB starting-condition file (64-bit-offset netCDF) lies outside the repository, so the committed CSV is the
product. Method (design section 4.3): latitude nearest 40.7 N; day = longitude 180 (local noon at 00 UT), night =
longitude 0 (local midnight); pressure p = hyam*P0 + hybm*PS; hypsometric integration of the layer interfaces
(Rd = 287.05 J/kg/K, constant g = 9.80665 m/s^2, layer temperature) gives the geopotential height Z (starting from the
surface geopotential PHIS/g, or from 0 when the file has no PHIS, as is the case for this file); the geometric
altitude is z = R Z/(R - Z) with R = 6371 km (z is about 1.6 km above Z at 100 km); number density of a species with
mass mixing ratio q is q * 28.96/M * p/(k_B T); N2 is the remainder of the air (n_total minus the other eight
species, so it also absorbs helium, argon, atomic nitrogen and nitric oxide). The native levels are interpolated to
a regular 1 km grid from 45 to 100 km (temperature linearly, pressure and every density log-linearly in altitude), so
that the table has exact end points for ambient.AmbientColumn.at.
"""
import os
import sys

import numpy as np
from scipy.io import netcdf_file

from reentry_model.wake_column import DATA
from reentry_model.wake_column.species import AMBIENT_SPECIES, load_species

K_B = 1.380649e-23      # J/K
R_D = 287.05            # J/kg/K, dry air
G0 = 9.80665            # m/s^2
M_AIR = 28.96           # amu, mean molecular mass of dry air used in the conversion
R_EARTH = 6371e3        # m, mean Earth radius for the geopotential-to-geometric conversion
TARGET_LAT = 40.7
Z_GRID_KM = np.arange(45.0, 100.0 + 0.5, 1.0)
# WACCM variable name for each ambient species whose mixing ratio is read (water vapour is "Q"; N2 is the remainder)
WACCM_VAR = {"O2": "O2", "O": "O", "O3": "O3", "H": "H", "H2": "H2", "H2O": "Q", "CO": "CO", "CO2": "CO2"}
COLUMNS = {"day": 180.0, "night": 0.0}
OUT_CSV = os.path.join(DATA, "ambient_waccm_40.7N.csv")
FIELDS = ("time_of_day", "z_km", "p_Pa", "T_K", "n_total_m3") + AMBIENT_SPECIES


def native_column(v, j, i, mass_amu):
    """Native-level profile {z_km, p_Pa, T_K, n_total_m3, species...} (top level first) at lat index j, lon index i."""
    hyam, hybm, hyai, hybi = (np.array(v[k][:], dtype=float) for k in ("hyam", "hybm", "hyai", "hybi"))
    p0 = float(v["P0"].getValue())
    ps = float(v["PS"][0, j, i])
    T = np.array(v["T"][0, :, j, i], dtype=float)
    p = hyam * p0 + hybm * ps
    p_int = hyai * p0 + hybi * ps
    z_int = np.zeros(len(p_int))                       # the bottom interface (last index) is the ground
    if "PHIS" in v:                                    # surface geopotential (m^2/s^2) -> surface height
        z_int[-1] = float(v["PHIS"][0, j, i]) / G0
    for k in range(len(p_int) - 2, -1, -1):
        z_int[k] = z_int[k + 1] + R_D * T[k] / G0 * np.log(p_int[k + 1] / p_int[k])
    z_geop = 0.5 * (z_int[:-1] + z_int[1:])            # geopotential height of each layer midpoint, m
    z_km = R_EARTH * z_geop / (R_EARTH - z_geop) / 1e3  # geometric altitude, km
    n_total = p / (K_B * T)
    out = {"z_km": z_km, "p_Pa": p, "T_K": T, "n_total_m3": n_total}
    for s, var in WACCM_VAR.items():
        q = np.array(v[var][0, :, j, i], dtype=float)
        out[s] = q * M_AIR / mass_amu[s] * n_total
    out["N2"] = n_total - sum(out[s] for s in WACCM_VAR)
    return out


def to_grid(col, z_grid_km):
    """Interpolate a native column to z_grid_km: T linear, everything else log-linear (linear if any value <= 0)."""
    order = np.argsort(col["z_km"])
    z = col["z_km"][order]
    if z_grid_km[0] < z[0] or z_grid_km[-1] > z[-1]:
        raise ValueError(f"native levels cover {z[0]:.1f}-{z[-1]:.1f} km, not the requested "
                         f"{z_grid_km[0]}-{z_grid_km[-1]} km")
    out = {"z_km": z_grid_km}
    for k, y in col.items():
        if k == "z_km":
            continue
        y = y[order]
        if k != "T_K" and np.all(y > 0):
            out[k] = np.exp(np.interp(z_grid_km, z, np.log(y)))
        else:
            out[k] = np.interp(z_grid_km, z, y)
    return out


def extract(path):
    """Return {"day": grid column, "night": grid column, "lat": latitude used} from the WACCM file at path."""
    mass = {s: sp.mass_amu for s, sp in load_species().items()}
    f = netcdf_file(path, "r", mmap=True)
    try:
        v = f.variables
        lat, lon = np.array(v["lat"][:], dtype=float), np.array(v["lon"][:], dtype=float)
        j = int(np.argmin(abs(lat - TARGET_LAT)))
        res = {"lat": float(lat[j])}
        for tod, lon0 in COLUMNS.items():
            i = int(np.argmin(abs(lon - lon0)))
            res[tod] = to_grid(native_column(v, j, i, mass), Z_GRID_KM)
            res[tod + "_lon"] = float(lon[i])
        return res
    finally:
        try:
            f.close()
        except RuntimeWarning:
            pass


def write_csv(res, path, source_name):
    with open(path, "w") as fh:
        fh.write(f"# WACCM ambient columns at {res['lat']:.1f} N (nearest grid latitude to {TARGET_LAT}), "
                 "extracted by reentry_model/wake_column/analysis/extract_waccm_columns.py\n")
        fh.write(f"# source file: {source_name} (full-chemistry WACCM starting condition, date 2001-01-01 00 UT; "
                 f"day = longitude {res['day_lon']:.0f} E (local noon), night = longitude {res['night_lon']:.0f} E "
                 "(local midnight))\n")
        fh.write("# pressure from hybrid coefficients; altitude = geometric altitude z = R Z/(R - Z), R = 6371 km, from the "
                 "geopotential height Z by hypsometric integration (Rd 287.05, g 9.80665) from the surface; this file has "
                 "no PHIS (surface geopotential), so the surface height is taken as 0 (terrain height at these two "
                 "columns is not applied, at most a few hundred metres);\n"
                 "# number density = mass mixing ratio * 28.96/M * p/(k_B T); N2 = n_total minus the other eight species;\n")
        fh.write("# native levels interpolated to a 1 km grid, 45-100 km (T linear, p and densities log-linear). "
                 "Units: z km, p Pa, T K, densities m^-3. Held fixed for 12 h in the model.\n")
        fh.write(",".join(FIELDS) + "\n")
        for tod in COLUMNS:
            c = res[tod]
            for k in range(len(c["z_km"])):
                row = [tod, f"{c['z_km'][k]:.1f}"] + [f"{c[name][k]:.6e}" for name in FIELDS[2:]]
                fh.write(",".join(row) + "\n")


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    res = extract(argv[1])
    out = argv[2] if len(argv) > 2 else OUT_CSV
    write_csv(res, out, os.path.basename(argv[1]))
    print(f"wrote {out}")
    for tod in COLUMNS:
        c = res[tod]
        k = int(np.argmin(abs(c["z_km"] - 71.0)))
        print(f"{tod} at {c['z_km'][k]:.0f} km: T {c['T_K'][k]:.1f} K, n_total {c['n_total_m3'][k]:.3e}, "
              f"O {c['O'][k]:.3e}, H {c['H'][k]:.3e}, H2O {c['H2O'][k]:.3e} m^-3")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
