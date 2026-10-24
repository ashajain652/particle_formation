"""Checks for droplet_heating_sweep.py; each prints its measured value and PASS/FAIL."""
import math
import sys

import numpy as np
from scipy import integrate

import droplet_heating_sweep as dh
import particle_trajectory as pt
import sphere_sweep

ok = True


def report(name, value, passed):
    global ok
    ok &= bool(passed)
    print("{:<70s} {:>14s}  {}".format(name, value, "PASS" if passed else "FAIL"))


# 1. element fluxes against direct quadrature over the Maxwellian (unit rho, T)
T, rho = 250.0, 1.0
c = math.sqrt(2 * pt.R_AIR * T)          # 1/beta
for s, sn in ((0.0, 0.0), (1.0, 0.3), (2.0, -0.8), (5.0, 2.0)):
    U = s * c
    un = sn * c                           # flow component along the inward normal
    ut = math.sqrt(max(U * U - un * un, 0.0))
    f = lambda w, v, u: rho / (math.pi ** 1.5 * c ** 3) * math.exp(-((u - un) ** 2 + (v - ut) ** 2 + w * w) / c ** 2)
    L = 8 * c
    mass = integrate.tplquad(lambda w, v, u: u * f(w, v, u), 0, max(un, 0) + L, ut - L, ut + L, -L, L, epsabs=0, epsrel=1e-8)[0]
    etr = integrate.tplquad(lambda w, v, u: 0.5 * (u * u + v * v + w * w) * u * f(w, v, u), 0, max(un, 0) + L,
                            ut - L, ut + L, -L, L, epsabs=0, epsrel=1e-8)[0]
    m2, e2 = dh.fm_element_fluxes(rho, T, s, np.array([sn]))
    report("element mass flux, s={:g} sn={:g}: rel. error".format(s, sn), "{:.2e}".format(abs(m2[0] / mass - 1)), abs(m2[0] / mass - 1) < 1e-6)
    report("element energy flux, s={:g} sn={:g}: rel. error".format(s, sn), "{:.2e}".format(abs(e2[0] / etr - 1)), abs(e2[0] / etr - 1) < 1e-6)

# 2. sphere limits
d = 10e-6
for s in (10.0, 30.0):
    V = s * c
    r = dh.fm_sphere_heat(1e-5, T, V, 300.0, d) / dh.fm_high_speed_limit(1e-5, V, d)
    report("sphere heat / (1/2 rho V^3 A) at s={:g} (-> 1 + O(1/s^2))".format(s), "{:.5f}".format(r), abs(r - 1) < 3.0 / s ** 2)
q0 = dh.fm_sphere_heat(1e-5, T, 1e-6, T, d)
report("sphere heat, V -> 0, T_wall = T_air (must vanish) [W]", "{:.2e}".format(q0), abs(q0) < 1e-12)
# s = 0 conduction against the closed form: alpha (2 + zeta/2) R (T - Tw) * rho/(2 sqrt(pi) beta) * 4 pi r^2
q = dh.fm_sphere_heat(1e-5, T, 0.0, 300.0, d)
ref = (2 + 1) * pt.R_AIR * (T - 300.0) * 1e-5 * c / (2 * math.sqrt(math.pi)) * math.pi * d * d
report("sphere heat at rest, wall 50 K hotter: closed form ratio", "{:.6f}".format(q / ref), abs(q / ref - 1) < 1e-9)
# quadrature convergence at the high-s kink
mu_save, w_save = dh.MU, dh.MU_W
x, w = np.polynomial.legendre.leggauss(400)
q48 = dh.fm_sphere_heat(1e-5, T, 7000.0, 1500.0, d)
dh.MU = np.concatenate((0.5 * (x - 1), 0.5 * (x + 1))); dh.MU_W = np.concatenate((0.5 * w, 0.5 * w))
q400 = dh.fm_sphere_heat(1e-5, T, 7000.0, 1500.0, d)
dh.MU, dh.MU_W = mu_save, w_save
report("sphere heat 2x48 vs 2x400 Gauss points: rel. difference", "{:.1e}".format(abs(q48 / q400 - 1)), abs(q48 / q400 - 1) < 1e-10)

# 3. parent state keyed on altitude against the measured values
parent = sphere_sweep.load_parent(dh.PARENT_JSON, dh.PARENT_OBJECT)
for h, v, g in ((90, 7.83, -0.37), (80, 7.62, -0.79), (70, 6.91, -1.53), (60, 4.21, -3.57), (50, 1.57, -10.6)):
    st = dh.parent_state_at_altitude(parent, h)
    report("parent at {} km: V = {:.3f} km/s, gamma = {:.2f} deg (expected ~{} / ~{})".format(
        h, st["velocity_kms"], st["flight_path_deg"], v, g), "",
        abs(st["velocity_kms"] - v) < 0.1 * v and abs(st["flight_path_deg"] - g) < 0.1 * abs(g))
    # the table quotes the nearest parent row; the interpolated state must lie between the bracketing rows
    a_, b_ = next((a_, b_) for a_, b_ in zip(parent.branch, parent.branch[1:]) if a_["altitude"] >= h >= b_["altitude"])
    report("  lies between parent rows at {:.2f} and {:.2f} km".format(a_["altitude"], b_["altitude"]), "",
           min(a_["velocity"], b_["velocity"]) <= st["velocity_kms"] <= max(a_["velocity"], b_["velocity"])
           and min(a_["path"], b_["path"]) <= st["flight_path_deg"] <= max(a_["path"], b_["path"]))
    report("  Vx, Vy recombine to V", "{:.1e}".format(abs(math.hypot(st["vx"], st["vy"]) / (st["velocity_kms"] * 1e3) - 1)),
           abs(math.hypot(st["vx"], st["vy"]) / (st["velocity_kms"] * 1e3) - 1) < 1e-12)
# against the velocity-keyed original: invert at the interpolated velocity, recover the altitude
st = dh.parent_state_at_altitude(parent, 76.3)
back = sphere_sweep.interpolate_state(parent, st["velocity_kms"])
report("altitude-keyed vs sphere_sweep.interpolate_state round trip at 76.3 km [km]",
       "{:.1e}".format(abs(back["altitude_km"] - 76.3)), abs(back["altitude_km"] - 76.3) < 1e-3)

# 4. energy closure and the hybrid hand-off on two cases
for dd, hh in ((10, 80), (250, 50)):
    st = dh.parent_state_at_altitude(parent, hh)
    s_, rows = dh.run_case(dd, hh, st["vx"], st["vy"], history=True)
    p1 = [r for r in rows if r["phase"] == 1]
    bal = s_["E_in_Jkg"] - s_["E_rad_Jkg"] - s_["E_conv_out_Jkg"]
    dT = dh.CP_LIQUID * (s_["T_handoff_K"] - dh.T_RELEASE)
    report("{} um at {} km: energy closure (in - rad - out) vs c_p dT".format(dd, hh),
           "{:.1e}".format(abs(bal / dT - 1) if dT else abs(bal)), abs(bal - dT) < 1e-4 * max(abs(dT), 1.0) + 1.0)
    report("{} um at {} km: hybrid hand-off speed mismatch with terminal".format(dd, hh),
           "{:.1e}".format(s_["hybrid_switch_mismatch"]), not (s_["hybrid_switch_mismatch"] > 1e-6))

print("\nALL PASS" if ok else "\nSOME CHECKS FAILED")
sys.exit(0 if ok else 1)
