"""Build the cross-check ambient table data/wake_column/ambient_msis_afgl_40.7N.csv from NRLMSIS 2.1 and AFGL.

Usage:  python -m reentry_model.wake_column.analysis.build_msis_afgl_columns [out.csv]

Needs data/wake_column/ambient_waccm_40.7N.csv (for H2 and CO) and data/wake_column/afgl_midlatitude_winter.csv.
Method (design section 4.3):
  * pymsis.calculate for 2001-01-01 00:00 at latitude 40.7 N, longitude 180 E (day, local noon) and 0 E (night, local
    midnight), on the 1 km grid from 45 to 100 km, with moderate fixed solar indices (F10.7 = 150 and its 81-day mean
    150 solar flux units, Ap = 4 in all seven slots) so that no space-weather file is downloaded.
  * N2, O2, O and H are the MSIS number densities; n_total is the sum of all MSIS number densities (N2, O2, O, He, H,
    Ar, N, anomalous O, NO); T is the MSIS temperature; p = n_total k_B T. The output indices are taken from
    pymsis.Variable, not hard-coded.
  * O3, H2O and CO2 are the AFGL midlatitude-winter volume mixing ratios times n_total (log-linear in altitude).
  * H2 and CO are copied from the WACCM column of the same time of day, because AFGL has no H2 and its CO is far too
    low a mixing ratio to represent the mesosphere (assumed; stated in the header).
  * pymsis returns NaN for a species the model does not supply at that altitude (atomic oxygen below about 51 km,
    atomic hydrogen below about 76 km, atomic nitrogen and nitric oxide lower down still). Missing N, anomalous O, NO
    count as zero in n_total (they are negligible there). Missing O or H is taken from the WACCM column of the same
    time of day (assumed; the filled altitude ranges are listed in the header), because zero would be wrong by orders
    of magnitude for the chemistry.
  The stored species therefore need not sum exactly to n_total: helium, argon, atomic nitrogen and nitric oxide are
  not ambient species of the network, and the other sources are added on top (a deviation below 2 percent is checked
  in the tests).
"""
import csv
import os
import sys

import numpy as np
import pymsis

from reentry_model.wake_column import DATA
from reentry_model.wake_column.ambient import AmbientColumn
from reentry_model.wake_column.analysis.extract_waccm_columns import COLUMNS, FIELDS, TARGET_LAT, Z_GRID_KM
from reentry_model.wake_column.species import AMBIENT_SPECIES

K_B = 1.380649e-23
AFGL_CSV = os.path.join(DATA, "afgl_midlatitude_winter.csv")
OUT_CSV = os.path.join(DATA, "ambient_msis_afgl_40.7N.csv")
DATE = np.datetime64("2001-01-01T00:00")
NUMBER_DENSITY_VARS = (pymsis.Variable.N2, pymsis.Variable.O2, pymsis.Variable.O, pymsis.Variable.HE,
                       pymsis.Variable.H, pymsis.Variable.AR, pymsis.Variable.N, pymsis.Variable.ANOMALOUS_O,
                       pymsis.Variable.NO)
MSIS_SPECIES = {"N2": pymsis.Variable.N2, "O2": pymsis.Variable.O2, "O": pymsis.Variable.O, "H": pymsis.Variable.H}
AFGL_VMR = {"O3": "O3_ppmv", "H2O": "H2O_ppmv", "CO2": "CO2_ppmv"}     # ambient species -> AFGL column (ppmv)
COPIED_FROM_WACCM = ("H2", "CO")


def read_afgl(path=AFGL_CSV):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def build(z_grid_km=Z_GRID_KM):
    """Return {"day": column, "night": column}, each a dict of arrays keyed like the CSV fields (after time_of_day)."""
    afgl = read_afgl()
    out = {}
    for tod, lon in COLUMNS.items():
        res = pymsis.calculate(DATE, lon, TARGET_LAT, z_grid_km, f107s=[150], f107as=[150], aps=[[4] * 7])
        res = np.asarray(res)[0, 0, 0]                       # (n_alt, n_variables)
        nd = {v: res[:, int(v)] for v in NUMBER_DENSITY_VARS}
        waccm = AmbientColumn.load("waccm", tod)
        filled = {}
        for s, var in MSIS_SPECIES.items():
            missing = np.isnan(nd[var])
            if missing.any():
                nd[var] = np.where(missing, [waccm.at(float(z)).n[s] for z in z_grid_km], nd[var])
                filled[s] = (float(z_grid_km[missing].min()), float(z_grid_km[missing].max()))
        n_total = sum(np.nan_to_num(v, nan=0.0) for v in nd.values())
        T = res[:, int(pymsis.Variable.TEMPERATURE)]
        col = {"z_km": np.asarray(z_grid_km, dtype=float), "T_K": T, "n_total_m3": n_total, "p_Pa": n_total * K_B * T}
        for s, var in MSIS_SPECIES.items():
            col[s] = nd[var]
        for s, name in AFGL_VMR.items():
            vmr = np.exp(np.interp(z_grid_km, afgl["z_km"], np.log(afgl[name]))) * 1e-6
            col[s] = vmr * n_total
        for s in COPIED_FROM_WACCM:
            col[s] = np.array([waccm.at(float(z)).n[s] for z in z_grid_km])
        for name in FIELDS[2:]:
            if not np.all(np.isfinite(col[name])):
                raise ValueError(f"non-finite {name} in the {tod} MSIS/AFGL column")
        col["filled"] = filled
        out[tod] = col
    return out


def write_csv(cols, path):
    with open(path, "w") as fh:
        fh.write("# Cross-check ambient columns at 40.7 N: NRLMSIS 2.1 (pymsis 0.13) for T, n_total, N2, O2, O, H; "
                 "AFGL midlatitude winter (Anderson et al. 1986, data/wake_column/afgl_midlatitude_winter.csv) "
                 "volume mixing ratios times n_total for O3, H2O, CO2;\n")
        fh.write("# H2 and CO copied from the WACCM column of the same time of day because AFGL lacks them (assumed). "
                 "Date 2001-01-01 00 UT; day = longitude 180 E (local noon), night = longitude 0 E (local midnight); "
                 "F10.7 = 150 and 81-day mean F10.7 = 150 solar flux units, Ap = 4 in all slots (fixed, moderate);\n")
        fh.write("# n_total = sum of all MSIS number densities (N2, O2, O, He, H, Ar, N, anomalous O, NO); "
                 "p = n_total k_B T; 1 km grid 45-100 km. Built by "
                 "reentry_model/wake_column/analysis/build_msis_afgl_columns.py. Units: z km, p Pa, T K, densities m^-3.\n")
        for tod in COLUMNS:
            for s, (lo, hi) in sorted(cols[tod]["filled"].items()):
                fh.write(f"# NRLMSIS returns no {s} at {lo:.0f}-{hi:.0f} km ({tod}); those values are taken from the "
                         "WACCM column of the same time of day (assumed).\n")
        fh.write(",".join(FIELDS) + "\n")
        for tod in COLUMNS:
            c = cols[tod]
            for k in range(len(c["z_km"])):
                row = [tod, f"{c['z_km'][k]:.1f}"] + [f"{c[name][k]:.6e}" for name in FIELDS[2:]]
                fh.write(",".join(row) + "\n")


def main(argv):
    out = argv[1] if len(argv) > 1 else OUT_CSV
    cols = build()
    write_csv(cols, out)
    print(f"wrote {out}")
    for tod in COLUMNS:
        c = cols[tod]
        closure = sum(c[s] for s in AMBIENT_SPECIES) / c["n_total_m3"]
        print(f"{tod}: sum(species)/n_total from {closure.min():.4f} to {closure.max():.4f}; "
              f"filled from WACCM: {c['filled']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
