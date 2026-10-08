"""Compare dtswitch runs: how much of the melt runs off before it is sprayed, how far it gets, and what the droplets are.

    compare.py <run dir> [<run dir> ...]       (each holds record.json, the history CSV and particles.npz)

Per run, over the whole run and over the continuum window (from the switch time of the first switched run on):
  * runoff share: the film mass that arrived on another patch by runoff before the spray stage, as a share of the mass
    sprayed (net per patch and step, so mass passing through a patch within one step is not counted twice);
  * the net mass that ran out across each latitude circle psi (the cap psi < psi_k loses it to the rest of the body),
    as a share of the film available inside that cap;
  * where the mass is sprayed against where it was available, by latitude band;
  * the equatorial ring: film that ran onto it, and what it sprayed;
  * the deep transport's arrivals (liquid below the conjugate depth, never sprayed);
  * the droplet population and the thick/thin split from the source table (b against delta_m at release).
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

THICK, THIN = 0, 1


def load(run):
    rec = json.load(open(os.path.join(run, "record.json")))
    csv = [p for p in glob.glob(os.path.join(run, "*.csv"))]
    hist = pd.read_csv(csv[0]) if csv else None
    npz = glob.glob(os.path.join(run, "*", "particles.npz"))
    src = np.load(npz[0]) if npz else None
    return rec, hist, src


def window_sums(steps, t0, t1, key):
    rows = [s[key] for s in steps if key in s and t0 < s["t"] <= t1]
    return np.sum(np.array(rows, dtype=float), axis=0) if rows else None


def summarise(run, t0, t1, label):
    rec, hist, src = load(run)
    st = rec["steps"]
    edges = rec["edges_deg"]
    avail = window_sums(st, t0, t1, "avail")
    arrive = window_sums(st, t0, t1, "arrive")
    leave = window_sums(st, t0, t1, "leave")
    rel = window_sums(st, t0, t1, "rel")
    released = sum(s["released_kg"] for s in st if t0 < s["t"] <= t1)
    ring_arr = sum(s.get("ring_arrive_kg", 0.0) for s in st if t0 < s["t"] <= t1)
    ring_rel = sum(s.get("ring_rel_kg", 0.0) for s in st if t0 < s["t"] <= t1)
    deep = sum(s.get("deep_arrive_kg", 0.0) for s in st if t0 < s["t"] <= t1)
    nsteps = sum(1 for s in st if t0 < s["t"] <= t1)
    print(f"\n{label}: {os.path.basename(run)}  window {t0:g}-{t1:g} s, {nsteps} melt steps, switch at {rec.get('switched_at_s')}")
    if avail is None:
        print("   no steps in the window")
        return {}
    print(f"   sprayed {released * 1e3:8.2f} g; deep liquid moved (per-step net arrivals) {deep * 1e3:7.2f} g")
    print(f"   equatorial ring: film arriving by runoff {ring_arr * 1e3:7.2f} g, sprayed there {ring_rel * 1e3:7.2f} g"
          f" ({ring_rel / max(released, 1e-30):5.1%} of all sprayed)")
    # Step-independent measures. The per-patch net arrivals above are summed per step, so film that creeps across N
    # patches over N small steps counts N times while the same journey inside one large step counts once: they are not
    # comparable between step sizes. The net flux across a fixed latitude circle is (a conserved quantity crossing a
    # fixed line), and so is where the mass is sprayed against where it entered the film.
    net_loss = leave - arrive                          # per band: net film the band lost to runoff over the window
    produced = rel + net_loss                          # per band: film that entered the band by every other route
    net_out = np.cumsum(net_loss)                      # net mass leaving the cap psi < edge_{k+1} by runoff
    cap_prod = np.cumsum(produced)
    centre = 0.5 * (np.array(edges[:-1]) + np.array(edges[1:]))
    centre[-1] = 95.0                                  # the lee band: beyond the equator
    shift = (np.sum(rel * centre) - np.sum(produced * centre)) / max(rel.sum(), 1e-30)
    print("   latitude band [deg]  entered film g  sprayed g | net runoff out of the cap below the band's upper edge")
    for k in range(len(edges) - 1):
        print(f"   {edges[k]:5.0f}-{edges[k + 1]:<5.0f}      {produced[k] * 1e3:9.2f}  {rel[k] * 1e3:9.2f}"
              f" | {net_out[k] * 1e3:8.2f} g = {net_out[k] / max(cap_prod[k], 1e-30):6.1%} of the film that entered that cap")
    print(f"   mass-weighted shift between where film entered and where it was sprayed: {shift:+.2f} deg of latitude;"
          f" per-step net arrivals (not comparable between steps) {arrive.sum() * 1e3:.1f} g")
    out = {"released_g": released * 1e3, "arrivals_per_step_g": arrive.sum() * 1e3, "shift_deg": float(shift),
           "ring_rel_g": ring_rel * 1e3, "ring_arrive_g": ring_arr * 1e3, "deep_moved_g": deep * 1e3,
           "net_out_g": (net_out * 1e3).tolist(), "net_out_share": (net_out / np.maximum(cap_prod, 1e-30)).tolist(),
           "produced_g": (produced * 1e3).tolist(), "rel_g": (rel * 1e3).tolist()}
    if src is not None:
        t = src["time_s"]
        sel = (t > t0) & (t <= t1)
        m, br, b, dm, r, n = src["dm_kg"][sel], src["branch"][sel], src["b_m"][sel], src["delta_m_m"][sel], src["r_m"][sel], src["dn"][sel]
        M = m.sum()
        gir = np.isfinite(dm) & ((br == THICK) | (br == THIN))
        thick = gir & (br == THICK)
        o = np.argsort(r)
        cn = np.cumsum(n[o])
        r_med = float(r[o][np.searchsorted(cn, 0.5 * cn[-1])]) if cn.size else float("nan")
        print(f"   droplets {n.sum():.4g}, median radius by number {r_med * 1e6:6.1f} um; thick branch {m[thick].sum() / max(M, 1e-30):5.1%} of the mass"
              f" (film deeper than delta_m {m[thick & (b > dm)].sum() / max(M, 1e-30):5.1%}); mass-weighted median film at release"
              f" {np.average(b[gir], weights=m[gir]) * 1e3 if gir.any() else float('nan'):.3f} mm (mean)")
        out.update({"droplets": float(n.sum()), "r_median_um": r_med * 1e6, "thick_mass_share": float(m[thick].sum() / max(M, 1e-30))})
    if hist is not None:
        h = hist[(hist.time_s > t0) & (hist.time_s <= t1)]
        if len(h):
            print(f"   history: sprayed {hist.sprayed_mass_kg.iloc[-1]:.4f} kg at {hist.time_s.iloc[-1]:g} s, mass {hist.mass_kg.iloc[-1]:.4f} kg,"
                  f" re-solidified {hist.frozen_mass_kg.iloc[-1] * 1e3:.2f} g, model runoff counter {hist.runoff_mass_kg.iloc[-1]:.3f} kg")
    return out


if __name__ == "__main__":
    runs = sys.argv[1:]
    recs = [json.load(open(os.path.join(r, "record.json"))) for r in runs]
    t_sw = min([r["switched_at_s"] for r in recs if r.get("switched_at_s") is not None] or [49.5])
    t_end = max(s["t"] for r in recs for s in r["steps"])
    res = {}
    for r in runs:
        res[r] = {"all": summarise(r, 0.0, t_end, "whole run"), "cont": summarise(r, t_sw, t_end, "continuum window")}
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(runs[0])), "compare.json"), "w"), indent=1)
