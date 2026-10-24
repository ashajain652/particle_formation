#!/usr/bin/env python3
"""Flight-level facts of one prototype run, measured the way facts 49-65 of the shared context state them.

    python3 flight_facts.py <run.json> [--t 120]

Prints, for the whole run and (with --t) up to time t: melt onset (time, altitude), end reason and time, landing or
demise state, final mass and its share of m0, sprayed mass, droplet count, the median radius by number and by mass of
all droplets released (from particles.npz, not the per-step medians), the largest droplet, the front-surface
Rayleigh-Taylor release, the thin branch's share of the sprayed mass, the re-solidified mass, the cascade and deep
masses, the time spraying ends, the enthalpy per kg sprayed, the energy balance, steps, Newton iterations, runtime."""
import argparse
import csv
import json
import os

import numpy as np


def weighted_median(x, w):
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return float(x[o][np.searchsorted(c, 0.5 * c[-1])])


def facts(doc, rows, p, t_cut=None):
    cols = {k: np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in rows]) for k in rows[0]}
    t = cols["time_s"]
    i = len(t) - 1 if t_cut is None else int(np.searchsorted(t, t_cut + 1e-9) - 1)
    sel = p["time_s"] <= t[i] + 1e-9
    r, dn, dm, br = p["r_m"][sel], p["dn"][sel], p["dm_kg"][sel], p["branch"][sel]
    rel = cols["released_mass_kg"]
    spraying = np.flatnonzero(rel[: i + 1] > 0.0)
    out = {
        "t_end_s": float(t[i]), "altitude_end_km": float(cols["altitude_km"][i]),
        "mass_end_kg": float(cols["mass_kg"][i]), "mass_end_frac_m0": float(cols["mass_kg"][i] / cols["mass_kg"][0]),
        "sprayed_kg": float(cols["sprayed_mass_kg"][i]), "n_droplets": float(dn.sum()),
        "r_median_number_um": weighted_median(r, dn) * 1e6, "r_median_mass_um": weighted_median(r, dm) * 1e6,
        "r_max_mm": float(r[dn > 0].max() * 1e3) if (dn > 0).any() else float("nan"),
        "rt_release_g": float(dm[br == 3].sum() * 1e3), "thin_share_of_sprayed": float(dm[br == 1].sum() / dm.sum()),
        "frozen_g": float(cols["frozen_mass_kg"][i] * 1e3),
        "cascade_g": float(cols["cascade_mass_kg"][i] * 1e3) if "cascade_mass_kg" in cols else None,
        "deep_runoff_g": float(cols["deep_runoff_mass_kg"][i] * 1e3) if "deep_runoff_mass_kg" in cols else None,
        "deep_surfaced_g": float(cols["deep_surfaced_mass_kg"][i] * 1e3) if "deep_surfaced_mass_kg" in cols else None,
        "deep_account_peak_g": float(np.nanmax(cols["deep_mass_kg"][: i + 1]) * 1e3) if "deep_mass_kg" in cols else None,
        "spraying_last_step_s": float(t[spraying[-1]]) if spraying.size else None,
        "thick_branch_median": float(np.nanmedian(np.where(cols["thick_branch_fraction"][: i + 1] > 0, cols["thick_branch_fraction"][: i + 1], np.nan))),
        "molten_depth_max_mm": float(np.nanmax(cols["molten_depth_max_mm"][: i + 1])),
        "cascade_active_steps": int((cols["cascade_passes"][: i + 1] > 0).sum()) if "cascade_passes" in cols else None,
        "cascade_passes_max": int(np.nanmax(cols["cascade_passes"][: i + 1])) if "cascade_passes" in cols else None,
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("json")
    ap.add_argument("--t", type=float, action="append", default=[])
    a = ap.parse_args()
    doc = json.load(open(a.json))
    res = doc["results"]
    rows = list(csv.DictReader(open(doc["files"]["csv"])))
    p = np.load(doc["files"]["particles"])
    head = {k: res.get(k) for k in ("end_reason", "final_time_s", "final_altitude_km", "impact_time_s", "melt_onset_time_s",
                                     "melt_onset_altitude_km", "spraying_onset_time_s", "demise_time_s", "demise_altitude_km",
                                     "initial_mass_kg", "final_mass_kg", "sprayed_mass_kg", "n_released", "frozen_mass_kg",
                                     "cascade_mass_kg", "cascade_passes_max", "cascade_capped_steps", "deep_runoff_mass_kg",
                                     "deep_surfaced_mass_kg", "deep_mass_kg", "n_source_rows", "n_macro_steps",
                                     "mean_newton_iterations", "energy_balance_residual", "melt_energy_balance_residual",
                                     "absorbed_heat_J", "removed_enthalpy_J", "runtime_s", "n_frames")}
    head["enthalpy_per_kg_sprayed_MJ"] = res["removed_enthalpy_J"] / res["removed_mass_kg"] / 1e6 if res.get("removed_mass_kg") else None
    head["seed"] = doc["settings"].get("seed")
    print(json.dumps({"run": doc["run_name"], "results": head, "whole": facts(doc, rows, p),
                      **{"to_{:g}s".format(t): facts(doc, rows, p, t) for t in a.t}}, indent=1))


if __name__ == "__main__":
    main()
