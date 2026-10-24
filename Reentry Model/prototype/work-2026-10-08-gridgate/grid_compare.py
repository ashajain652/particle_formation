"""Stage 0 grid gate (sub-plan 18, Task 0): the droplet population and the masses of runs on different surface meshes.

    grid_compare.py <run dir> [<run dir> ...]

Each run dir holds the model's <name>.csv, <name>.json and <name>/particles.npz (a model run or a switched-step harness
run). Over the continuum window (from the step switch, 49.5 s, to the end), per run: the droplet count, the thick
branch's share of the sprayed mass (Girin's thick and thin branches only, as fact 93), the median radius by number and
by mass, the mass sprayed in the window; and from the history the sprayed mass and the body's mass at the last row,
and the run time. The definitions are those of the switched-step series' compare.py and series_summary (facts 91-95).
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

THICK, THIN = 0, 1
T_SWITCH = 49.5


def weighted_median(values, weights):
    o = np.argsort(values)
    c = np.cumsum(weights[o])
    return float(values[o][np.searchsorted(c, 0.5 * c[-1])]) if c.size and c[-1] > 0 else float("nan")


def summarise(run):
    csv = [p for p in glob.glob(os.path.join(run, "*.csv"))]
    js = [p for p in glob.glob(os.path.join(run, "*.json")) if not p.endswith("record.json")]
    npz = glob.glob(os.path.join(run, "*", "particles.npz"))
    hist = pd.read_csv(csv[0])
    doc = json.load(open(js[0])) if js else {}
    src = np.load(npz[0])
    t = src["time_s"]
    sel = t > T_SWITCH
    m, br, dm, r, n = src["dm_kg"][sel], src["branch"][sel], src["delta_m_m"][sel], src["r_m"][sel], src["dn"][sel]
    gir = np.isfinite(dm) & ((br == THICK) | (br == THIN))
    M = m.sum()
    res = doc.get("results", {})
    return {"run": os.path.basename(os.path.normpath(run)),
            "droplets": float(n.sum()),
            "thick_share": float(m[gir & (br == THICK)].sum() / max(M, 1e-30)),
            "r_median_number_um": weighted_median(r, n) * 1e6,
            "r_median_mass_um": weighted_median(r, m) * 1e6,
            "sprayed_window_g": float(M) * 1e3,
            "t_end_s": float(hist.time_s.iloc[-1]),
            "sprayed_g": float(hist.sprayed_mass_kg.iloc[-1]) * 1e3,
            "mass_end_g": float(hist.mass_kg.iloc[-1]) * 1e3,
            "runtime_h": float(res.get("runtime_s", float("nan"))) / 3600.0,
            "n_elements": doc.get("settings", {}).get("n_elements"),
            "n_patches": doc.get("settings", {}).get("n_patches")}


if __name__ == "__main__":
    rows = [summarise(r) for r in sys.argv[1:]]
    keys = list(rows[0])
    for k in keys:
        print("{:22s}".format(k) + "".join("{:>24}".format(("{:.6g}".format(row[k]) if isinstance(row[k], float) else str(row[k]))[:23]) for row in rows))
    json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "grid_compare.json"), "w"), indent=1)
