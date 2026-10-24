"""Limits of singh_drag and a map of its ratio to the current composed drag law (particle_trajectory)."""
import math
import sys

import numpy as np

import droplet_demise as dd  # noqa: F401  (puts the settling folder on the path)
import particle_trajectory as pt
import singh_drag as sg

ok = True


def report(name, value, passed):
    global ok
    ok &= bool(passed)
    print("{:<74s} {:>14s}  {}".format(name, value, "PASS" if passed else "FAIL"))


# limits
cd, _ = sg.drag_coefficient(1e-6, 1e-4, 1.0)        # creeping, continuum: Kn = 1.5e-2 * ... tiny
r = cd / (24.0 / 1e-4)
cd2, g2 = sg.drag_coefficient(1e-8, 1e-4, 1.0)
report("Stokes limit (Re 1e-4, Kn 1.5e-4): C_d / (24/Re)", "{:.4f}".format(cd2 / 2.4e5), abs(cd2 / 2.4e5 - 1) < 0.01)
cd, g = sg.drag_coefficient(10.0, 1e-3, 1.0)
report("free-molecule limit (M 10, Kn {:.1e}): C_d / Schaaf-Chambre".format(g["Kn"]), "{:.5f}".format(cd / g["cd_fm"]),
       abs(cd / g["cd_fm"] - 1) < 1e-3)
cd, g = sg.drag_coefficient(20.0, 1e7, 1.0)
report("hypersonic continuum (M 20, Re 1e7)", "{:.3f}".format(cd), abs(cd - 0.9) < 0.05)
for Re in (1.0, 10.0, 30.0):
    cd, g = sg.drag_coefficient(0.01, Re, 1.0)
    sn = 24 / Re * (1 + 0.15 * Re ** 0.687)
    report("low Mach, Re {:g}: Singh / Schiller-Naumann (Kn {:.1e})".format(Re, g["Kn"]), "{:.3f}".format(cd / sn),
           abs(cd / sn - 1) < 0.1)
for M in (0.999999, 1.000001):
    pass
lo, _ = sg.drag_coefficient(0.999999, 50.0, 1.0); hi, _ = sg.drag_coefficient(1.000001, 50.0, 1.0)
report("continuity at Mach 1 (Re 50)", "{:.1e}".format(abs(hi / lo - 1)), abs(hi / lo - 1) < 1e-4)

# map: ratio Singh / current law in air at 220 K, 2 Pa (about 77 km), particle at 220 K and at 1800 K
T, p = 220.0, 2.0
rho = p / (pt.R_AIR * T)
mu = pt.viscosity(T)
a = math.sqrt(pt.GAMMA * pt.R_AIR * T)
print("\nratio C_d(Singh) / C_d(current), air 220 K; rows: Knudsen number lambda/D, columns: speed ratio s")
ss = (1e-3, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 25.0)
print("   Kn \\ s " + "".join("{:>8g}".format(s) for s in ss))
lam = mu / rho * math.sqrt(math.pi / (2 * pt.R_AIR * T))     # lambda consistent with Singh's Kn
for Tw in (T, 1800.0):
    print("  wall {:.0f} K".format(Tw))
    for kn in (1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0):
        d = lam / kn
        row = []
        for s in ss:
            V = s * math.sqrt(2 * pt.R_AIR * T)
            M, Re = V / a, rho * V * d / mu
            c_s, _ = sg.drag_coefficient(M, Re, Tw / T)
            c_c, _ = pt.drag_coefficient(V, rho, T, p, d, Tw / T, "blend")
            row.append(c_s / c_c)
        print("  {:>7g} ".format(kn) + "".join("{:>8.3f}".format(x) for x in row))
print("\nALL PASS" if ok else "\nSOME CHECKS FAILED")
sys.exit(0 if ok else 1)
