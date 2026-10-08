"""Throwaway summary of the 2026-10-07 time-step series (Scheil, rigid substrate on), in fact 75's measures: over the
continuum window 49.5-120 s from the source table -- droplet count, median radius by number and by mass, the thick
branch's share of the sprayed mass, each branch's mass and median radii, the front-surface Rayleigh-Taylor release, the
mass-weighted mean film at release -- and from the history the sprayed mass and the mass at 120 s, the deep runoff, the
re-solidification counter and the medians of the rigid-substrate columns over the window; plus the run's wall time."""
import glob, json, os, sys
import numpy as np, pandas as pd
NAMES = {0: "thick", 1: "thin", 2: "rarefied", 3: "RT", 4: "regime2"}

def wmed(x, w):
    if not len(x): return float("nan")
    o = np.argsort(x); x, w = x[o], w[o]; c = np.cumsum(w)
    return float(x[np.searchsorted(c, 0.5 * c[-1])])

def summary(run):
    rec = json.load(open(os.path.join(run, "record.json")))
    h = pd.read_csv(glob.glob(os.path.join(run, "*.csv"))[0])
    p = np.load(glob.glob(os.path.join(run, "*", "particles.npz"))[0])
    j = json.load(open(glob.glob(os.path.join(run, "*.json"))[-1] if not glob.glob(os.path.join(run, "model*.json")) else glob.glob(os.path.join(run, "model*.json"))[0]))["results"]
    t, dn, dm, r, br, b = p["time_s"], p["dn"], p["dm_kg"], p["r_m"], p["branch"], p["b_m"]
    w = (t > 49.5) & (dm > 0)
    out = {"droplets_M": dn[w].sum() / 1e6, "r_num_um": wmed(r[w], dn[w]) * 1e6, "r_mass_um": wmed(r[w], dm[w]) * 1e6,
           "thick_share_pct": 100 * dm[w & (br == 0)].sum() / dm[w].sum(), "film_at_release_mm": 1e3 * (dm[w] * b[w]).sum() / dm[w].sum(),
           "sprayed_kg": float(h["sprayed_mass_kg"].iloc[-1]), "mass_120_kg": float(h["mass_kg"].iloc[-1]),
           "deep_runoff_g": 1e3 * float(h["deep_runoff_mass_kg"].iloc[-1]), "frozen_g": 1e3 * float(h["frozen_mass_kg"].iloc[-1]),
           "rt_g_window": 1e3 * dm[w & (br == 3)].sum(), "rt_g_all": 1e3 * dm[(dm > 0) & (br == 3)].sum(),
           "sprayed_before_switch_g": 1e3 * dm[(t < 49.5) & (dm > 0)].sum(), "sprayed_at_49p5_g": 1e3 * dm[(t == 49.5) & (dm > 0)].sum(), "energy": j["energy_balance_residual"],
           "wall_h": rec.get("wall_s", float("nan")) / 3600, "exit": rec.get("exit_code")}
    for c in (0, 1):
        m = w & (br == c)
        out[NAMES[c] + "_g"] = 1e3 * dm[m].sum()
        out[NAMES[c] + "_r_num_um"] = wmed(r[m], dn[m]) * 1e6
        out[NAMES[c] + "_r_mass_um"] = wmed(r[m], dm[m]) * 1e6
        out[NAMES[c] + "_droplets_M"] = dn[m].sum() / 1e6
    hw = h[h["time_s"] > 49.5]
    for c in ("thick_branch_fraction", "slurry_thick_fraction", "rigid_thin_fraction", "nonrigid_depth_mean_mm", "delta_m_mean_um", "molten_depth_mean_mm"):
        v = hw[c].to_numpy(float); v = v[np.isfinite(v)]
        out["med_" + c] = float(np.median(v)) if v.size else float("nan")
    # thick share by sub-window, as fact 75 reported for the finest step
    for a, z in ((49.5, 60), (60, 80), (80, 100), (100, 120)):
        m = (t > a) & (t <= z) & (dm > 0)
        out["thick_pct_%g_%g" % (a, z)] = 100 * dm[m & (br == 0)].sum() / dm[m].sum() if m.any() else float("nan")
    return out

runs = sys.argv[1:]
rows = {os.path.basename(r.rstrip("/")): summary(r) for r in runs}
keys = list(next(iter(rows.values())).keys())
print("%-24s" % "", " ".join("%14s" % k for k in rows))
for k in keys:
    print("%-24s" % k, " ".join("%14.5g" % v[k] if isinstance(v[k], (int, float)) and v[k] is not None else "%14s" % v[k] for v in rows.values()))
