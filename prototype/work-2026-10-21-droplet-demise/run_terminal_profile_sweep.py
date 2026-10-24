"""Final states and local terminal-velocity profiles for 15 chosen diameters (0.01-250 um) x 21 release altitudes
(90-50 km), Singh drag, no vertical wind. Writes into runs/terminal_profile_sweep/:

  final_states.csv          one row per case: d0_um, h0_km, V_final_ms, h_final_km, r_final_um
  terminal_velocity_profile.csv
                            one row per case: d0_um, h0_km, d_final_um, then v_t_<h>m [m/s] for h = 0, 50, ...,
                            100000 m: the speed at which drag balances weight for the droplet's final mass and
                            diameter in the local air (US76 table, droplet at the air temperature) at the final
                            state's latitude and longitude. Columns above h_final_km are what the droplet would
                            have there; it only ever sinks through those below.
"""
import csv
import math
import os
import time
from multiprocessing import Pool

import numpy as np

import droplet_demise as dd

DIAMETERS = (0.01, 0.1, 2.0, 6.0, 10.0, 14.0, 18.0, 25.0, 35.0, 45.0, 50.0, 100.0, 150.0, 200.0, 250.0)
ALTITUDES = tuple(range(90, 49, -2))
H_GRID_M = np.arange(0, 100000 + 1, 50)
OUT = os.path.join(dd.HERE, "runs", "terminal_profile_sweep")


def worker(job):
    d0, h0 = job
    s, rows = dd.run_case(d0, h0)
    if rows is None:
        return s, None
    last = rows[-1]
    lat, lon = math.radians(last["lat_deg"]), math.radians(last["lon_deg"])
    m, d = last["mass_kg"], last["d_um"] * 1e-6
    sim = dd.DropletSimulator(dd.parent_release(dd.sphere_sweep.load_parent(dd.PARENT_JSON, dd.PARENT_OBJECT), h0), d0)
    atm = dd.US76TableAtmosphere()
    vt = []
    for h in H_GRID_M:
        fs = atm.state(0.0, float(h), lat, lon)
        vt.append(sim.terminal_velocity(fs, dd.earth.geodetic_to_ecef(float(h), lat, lon), fs.T, d, m))
    return s, vt


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    jobs = sorted(((d, h) for d in DIAMETERS for h in ALTITUDES), key=lambda j: -j[0])
    t0 = time.time()
    with Pool() as pool:
        res = pool.map(worker, jobs, chunksize=1)
    print("{} cases in {:.0f} s".format(len(res), time.time() - t0))
    res.sort(key=lambda r: (r[0]["d0_um"], -r[0]["h0_km"]))
    nan = float("nan")
    with open(os.path.join(OUT, "final_states.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["d0_um", "h0_km", "V_final_ms", "h_final_km", "r_final_um"])
        for s, _ in res:
            w.writerow([s["d0_um"], s["h0_km"], s.get("V_final_ms", nan), s.get("h_final_km", nan),
                        s.get("d_final_um", 2 * nan) / 2.0])
    with open(os.path.join(OUT, "terminal_velocity_profile.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["d0_um", "h0_km", "d_final_um"] + ["v_t_{:d}m".format(int(h)) for h in H_GRID_M])
        for s, vt in res:
            w.writerow([s["d0_um"], s["h0_km"], s.get("d_final_um", nan)] + (vt or [nan] * len(H_GRID_M)))
    print(sorted({s["end_reason"] for s, _ in res}))
