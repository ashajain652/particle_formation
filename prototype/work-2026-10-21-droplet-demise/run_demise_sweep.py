"""Droplet demise over 15 log-spaced diameters (1-200 um) x 21 release altitudes (90-50 km).

Writes runs/demise_sweep/final_states.csv (one row per case, skipped cases included), the per-case histories, and,
with --sensitivity, final_states_mach.csv for heating_mode="mach" on three diameters.
"""
import argparse
import csv
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

import droplet_demise as dd

DIAMETERS = tuple(float("{:.3g}".format(v)) for v in np.logspace(0.0, np.log10(200.0), 15))
# 1 nm to 500 nm, 15 points evenly spaced in log, in micrometres
DIAMETERS_NANO = tuple(float("{:.3g}".format(v)) for v in np.logspace(np.log10(0.001), np.log10(0.5), 15))
ALTITUDES = tuple(range(90, 49, -2))
OUT = os.path.join(dd.HERE, "runs", "demise_sweep")
OUT_NANO = os.path.join(dd.HERE, "runs", "demise_sweep_nano")


def history_name(d, h, mode="regime", drag="singh"):
    return "demise_d{:g}um_h{:02d}km{}{}.csv".format(d, int(h), "" if mode == "regime" else "_" + mode,
                                                   "" if drag == "singh" else "_" + drag)


def worker(job):
    d, h, mode, write, out, drag, wind, wind_s = job
    t0 = time.time()
    try:
        s, rows = dd.run_case(d, h, heating_mode=mode, drag_law=drag, w_meso=wind, w_strat=wind_s)
    except Exception as exc:
        return dict(d0_um=d, h0_km=h, heating_mode=mode, skipped=False, end_reason="error: {!r}".format(exc))
    s["heating_mode"], s["drag_law"], s["wall_s"] = mode, drag, time.time() - t0
    if rows and write:
        keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[-1]).index(k) if k in rows[-1] else 999)
        with open(os.path.join(out, "histories", history_name(d, h, mode, drag)), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader(); w.writerows(rows)
    return s


def write(rows, name):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sensitivity", action="store_true")
    ap.add_argument("--nano", action="store_true", help="1-500 nm diameters into runs/demise_sweep_nano")
    ap.add_argument("--drag", choices=("singh", "blend"), default="singh", help="drag law (default: singh)")
    ap.add_argument("--wind", type=float, default=None,
                    help="mean mesospheric vertical wind [cm/s, + up]; also follows every droplet to the ground. "
                         "Output goes to <sweep folder>_w<value>cms_s<strat>mms")
    ap.add_argument("--wind-strat", type=float, default=0.0, help="mean stratospheric vertical wind [mm/s, + up]")
    a = ap.parse_args()
    if a.nano:
        OUT, DIAMETERS = OUT_NANO, DIAMETERS_NANO
    wind = None if a.wind is None else a.wind / 100.0
    if wind is not None:
        OUT = OUT + "_w{:+g}cms_s{:+g}mms".format(a.wind, a.wind_strat)
    os.makedirs(os.path.join(OUT, "histories"), exist_ok=True)
    wind_s = a.wind_strat / 1000.0
    jobs = [(d, h, "regime", True, OUT, a.drag, wind, wind_s) for d in DIAMETERS for h in ALTITUDES]
    if a.sensitivity:
        jobs += [(d, h, "mach", True, OUT, a.drag, wind, wind_s) for d in (1.0, 14.1, 200.0) for h in ALTITUDES]
    # largest, slowest cases first so the pool stays busy to the end
    jobs.sort(key=lambda j: -j[0])
    t0 = time.time()
    with Pool() as pool:
        out = pool.map(worker, jobs, chunksize=1)
    print("{} cases in {:.0f} s".format(len(out), time.time() - t0))
    for o in out:
        if str(o.get("end_reason", "")).startswith("error"):
            print(o["d0_um"], o["h0_km"], o["heating_mode"], o["end_reason"])
    key = lambda r: (r["d0_um"], -r["h0_km"])
    write(sorted([o for o in out if o["heating_mode"] == "regime"], key=key), "final_states{}.csv".format("" if a.drag == "singh" else "_" + a.drag))
    if a.sensitivity:
        write(sorted([o for o in out if o["heating_mode"] == "mach"], key=key), "final_states_mach{}.csv".format("" if a.drag == "singh" else "_" + a.drag))
    print(sorted({o["end_reason"] for o in out}))
