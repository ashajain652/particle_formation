"""When does each baseline droplet reach its terminal velocity, and what is that velocity?

The terminal velocity is a moving target: it is the local settling speed (particle_trajectory.terminal_velocity
at the droplet's current altitude) and keeps falling as the droplet sinks into denser air.  A droplet is said to
have reached it at the first sampled instant after which its speed stays within TOL of the local terminal
velocity for the rest of the history.  The history is scanned backwards from the ground, so only the
samples near the approach need a terminal-velocity solve.

Writes runs/droplet_sweep/terminal_approach.csv: d_um, h0_km, t_terminal_s, h_terminal_km, v_terminal_ms,
and the speed ratio V/v_t at the last sample still outside the band (a check that the band was crossed cleanly).
"""
import math
import os
import sys
from multiprocessing import Pool

import pandas as pd

import droplet_heating_sweep as dh
import particle_trajectory as pt
from reentry_model.atmosphere import US76TableAtmosphere

RUN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(dh.HERE, "runs", "droplet_sweep")
TOL = 0.01


def one(job):
    d, h = job
    atm = US76TableAtmosphere()
    prt = pt.Particle(d * 1e-6, dh.RHO_LIQUID)
    r = pd.read_csv(os.path.join(RUN, "histories", dh.history_name(d, h)))
    t, hh, V = r.t_s.values, r.h_m.values, r.V.values
    i_in = None
    ratio_out = float("nan")
    for i in range(len(r) - 1, -1, -1):
        vt = pt.terminal_velocity(prt, max(hh[i], 0.0), atm)
        if abs(V[i] / vt - 1.0) > TOL:
            ratio_out = V[i] / vt
            break
        i_in = i
    vt_in = pt.terminal_velocity(prt, max(hh[i_in], 0.0), atm)
    # interpolate the band crossing between the last sample outside it and the first inside, in log time
    t_cross, h_cross = t[i_in], hh[i_in]
    if i_in > 0 and t[i_in - 1] > 0:
        e_out, e_in = abs(ratio_out - 1.0), abs(V[i_in] / vt_in - 1.0)
        f = (e_out - TOL) / (e_out - e_in)
        t_cross = math.exp(math.log(t[i_in - 1]) + f * (math.log(t[i_in]) - math.log(t[i_in - 1])))
        h_cross = hh[i_in - 1] + f * (hh[i_in] - hh[i_in - 1])
    vt_cross = pt.terminal_velocity(prt, max(h_cross, 0.0), atm)
    return dict(d_um=d, h0_km=h, t_terminal_s=t_cross, h_terminal_km=h_cross / 1e3, v_terminal_ms=vt_cross,
                t_terminal_sample_s=t[i_in], sample_spacing=t[i_in] / t[i_in - 1] - 1.0 if i_in > 0 and t[i_in - 1] > 0 else 0.0,
                speed_over_vt_just_before=ratio_out, n_samples=len(r))


if __name__ == "__main__":
    jobs = [(d, h) for d in dh.DIAMETERS for h in dh.ALTITUDES]
    with Pool() as pool:
        out = pd.DataFrame(pool.map(one, jobs))
    out.to_csv(os.path.join(RUN, "terminal_approach.csv"), index=False)
    print(out.pivot(index="h0_km", columns="d_um", values="v_terminal_ms").iloc[::4].to_string(float_format=lambda v: "%.3g" % v))
    print(out.pivot(index="h0_km", columns="d_um", values="t_terminal_s").iloc[::4].to_string(float_format=lambda v: "%.3g" % v))
    print(out.pivot(index="h0_km", columns="d_um", values="h_terminal_km").iloc[::4].to_string(float_format=lambda v: "%.3g" % v))
    print("speed/v_t just before entering the band: min {:.3f} max {:.3f}".format(
        out.speed_over_vt_just_before.min(), out.speed_over_vt_just_before.max()))
