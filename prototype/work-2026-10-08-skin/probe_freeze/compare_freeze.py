"""Compare two skin runs' freeze-back and thickness diagnostics: the capped rule (278091b) and the energy-limited one.

    compare_freeze.py <capped run dir> <energy run dir>
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd


def load(run):
    csv = glob.glob(os.path.join(run, "*.csv"))[0]
    js = [p for p in glob.glob(os.path.join(run, "*.json")) if not p.endswith("freeze_probe.json")][0]
    return pd.read_csv(csv), json.load(open(js))


rows = {}
for label, run in (("capped", sys.argv[1]), ("energy", sys.argv[2])):
    h, doc = load(run)
    r = doc["results"]
    def col(name, default=np.nan):
        return h[name] if name in h else pd.Series([default] * len(h))
    rows[label] = {
        "t_end_s": float(h.time_s.iloc[-1]),
        "frozen_g (cumulative re-solidified)": float(h.frozen_mass_kg.iloc[-1]) * 1e3,
        "skin_frozen_g": float(col("skin_frozen_mass_kg").iloc[-1]) * 1e3,
        "skin_returned_g": float(col("skin_returned_mass_kg").iloc[-1]) * 1e3,
        "skin_full_events": float(col("skin_full_events").iloc[-1]),
        "skins above 2s (max over steps)": float(np.nanmax(col("skin_thick_count"))),
        "skins above 2s (mean over steps)": float(np.nanmean(col("skin_thick_count"))),
        "max skin thickness mm (max over steps)": float(np.nanmax(col("skin_thickness_max_mm"))),
        "sprayed_g": float(h.sprayed_mass_kg.iloc[-1]) * 1e3,
        "mass_end_g": float(h.mass_kg.iloc[-1]) * 1e3,
        "film_left_g": float(h.film_mass_kg.iloc[-1]) * 1e3,
        "skin_drawn_g": float(col("skin_drawn_mass_kg").iloc[-1]) * 1e3,
        "skin_shortfall_g (max over steps)": float(np.nanmax(col("skin_shortfall_kg"))) * 1e3,
        "interface_repeats": float(col("interface_repeats").iloc[-1]),
        "energy residual": float(r.get("melt_energy_balance_residual", np.nan)),
        "droplets (million)": float(r.get("n_released", np.nan)) / 1e6,
        "runtime_min": float(r.get("runtime_s", np.nan)) / 60.0,
    }
for k in rows["capped"]:
    print("{:42s}{:>16.6g}{:>16.6g}".format(k, rows["capped"][k], rows["energy"][k]))
