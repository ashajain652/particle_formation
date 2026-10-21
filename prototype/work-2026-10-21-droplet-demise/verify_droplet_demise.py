"""Checks for droplet_demise.py; each prints its measured value and PASS/FAIL."""
import math
import sys

import numpy as np
import pandas as pd

import droplet_demise as dd

ok = True


def report(name, value, passed):
    global ok
    ok &= bool(passed)
    print("{:<78s} {:>12s}  {}".format(name, value, "PASS" if passed else "FAIL"))


parent = dd.sphere_sweep.load_parent(dd.PARENT_JSON, dd.PARENT_OBJECT)

# 1. vapour pressures hit their anchors
for el in dd.VOLATILE:
    r = dd.p_sat(el, el.T_boil) / 101325.0
    report("{} vapour pressure at its normal boiling point [atm]".format(el.name), "{:.4f}".format(r), abs(r - 1) < 0.01)

# 2. Hertz-Knudsen flux of pure Al at 1500 K against the textbook number p sqrt(M / 2 pi R T)
sim = dd.DropletSimulator(dd.parent_release(parent, 80.0), 10.0)
fs = sim.atmosphere.state(0.0, 80e3, 0.0, 0.0)
m_vol = np.array([sim.m0 * e.mass_fraction for e in dd.VOLATILE])
mdot, y = sim.evaporation(fs, 1.0, 1500.0, m_vol, 10e-6)
x_al = sim.mole_fractions(m_vol)[0]
ref = math.pi * 1e-10 * x_al * dd.p_sat(dd.VOLATILE[0], 1500.0) * math.sqrt(0.026982 / (2 * math.pi * dd.R_GAS * 1500.0))
report("Al evaporation at 1500 K, 80 km (rarefied): model / Hertz-Knudsen", "{:.6f}".format(mdot[0] / ref),
       abs(mdot[0] / ref - 1) < 0.01)
report("  Zn+Mg vapour pressure over ambient at release temperature, 80 km", "{:.0f}".format(
    sim.evaporation(fs, 1.0, dd.T_RELEASE, m_vol, 10e-6)[1]), sim.evaporation(fs, 1.0, dd.T_RELEASE, m_vol, 10e-6)[1] > 1)
# at sea-level pressure and below boiling, diffusion must throttle the kinetic rate
fs_sl = sim.atmosphere.state(0.0, 0.0, 0.0, 0.0)
md_sl, y_sl = sim.evaporation(fs_sl, 1.0, 1100.0, m_vol, 100e-6)
hk = math.pi * 1e-8 * sim.mole_fractions(m_vol)[1] * dd.p_sat(dd.VOLATILE[1], 1100.0) * math.sqrt(
    0.06538 / (2 * math.pi * dd.R_GAS * 1100.0))
report("Zn at 1100 K, sea level, 100 um (below boiling): diffusion-limited / kinetic", "{:.3g}".format(md_sl[1] / hk),
       y_sl < 1 and md_sl[1] / hk < 0.5)

# 3. Weber pre-check skips an unstable droplet
we_save = dd.WE_CRIT
dd.WE_CRIT = 1e-3
s, rows = dd.run_case(10.0, 80.0, parent)
dd.WE_CRIT = we_save
report("Weber pre-check with a critical value of 0.001 skips the case", s["end_reason"][:7], s["skipped"] and rows is None)

# 4. evaporation off: peak temperature against the 2026-10-14 free-molecule sweep at accommodation 0.8
s, rows = dd.run_case(10.0, 80.0, parent, alpha_evap=0.0)
old = pd.read_csv(dd.SETTLING + "/runs/droplet_sweep/summary_eps0.1_alpha0.8.csv")
o = old[(old.d_um == 10)]
t_old = float(np.interp(80.0, o.h0_km.values[::-1], o.T_peak_K.values[::-1]))
report("no evaporation, 10 um at 80 km: peak T vs flat-Earth free-molecule sweep", "{:.0f}/{:.0f} K".format(
    s["T_peak_K"], t_old), abs(s["T_peak_K"] / t_old - 1) < 0.05)
report("  mass conserved with evaporation off (lost fraction)", "{:.1e}".format(s["mass_lost_frac"]),
       abs(s["mass_lost_frac"]) < 1e-9)

# 5. the hand-over between PhysicsHeating and the exchange model, wherever a flight actually blended them,
#    as a fraction of the incident energy flux alpha 1/2 rho V^3 A (the relative gap is meaningless near recovery)
gaps = []
for d_, h_ in ((10.0, 80.0), (100.0, 60.0), (200.0, 52.0)):
    s_, rows_ = dd.run_case(d_, h_, parent)
    sim = dd.DropletSimulator(dd.parent_release(parent, h_), d_)
    for r in rows_:
        if 0.0 < r["w_physics_heating"] < 1.0:
            fs = sim.atmosphere.state(0.0, r["h_km"] * 1e3, 0.0, 0.0)
            dm = r["d_um"] * 1e-6
            scale = dd.ALPHA_HEAT * 0.5 * fs.rho * r["V_ms"] ** 3 * math.pi * dm * dm / 4.0
            gaps.append(abs(sim.q_physics(fs, r["V_ms"], r["T_K"], dm) - sim.q_exchange(fs, r["V_ms"], r["T_K"], dm)) / scale)
report("heat models where blended: max |PhysicsHeating - exchange| / incident flux ({} states)".format(len(gaps)),
       "{:.3f}".format(max(gaps)) if gaps else "none", (max(gaps) < 0.3) if gaps else True)
s, rows = dd.run_case(10.0, 80.0, parent)
report("final state: speed within 1 % of terminal, T within 1 K of equilibrium", "{:.4f} / {:.2f} K".format(
    s["V_final_ms"] / s["v_terminal_final_ms"], s["T_final_K"] - s["T_eq_final_K"]),
       s["end_reason"] == "steady" and abs(s["V_final_ms"] / s["v_terminal_final_ms"] - 1) < 0.01
       and abs(s["T_final_K"] - s["T_eq_final_K"]) < 1.0)
lost = sum(s["lost_{}_kg".format(e.name)] for e in dd.VOLATILE)
report("element losses add up to the total mass lost", "{:.1e}".format(abs(lost / s["mass_lost_kg"] - 1)),
       abs(lost / s["mass_lost_kg"] - 1) < 1e-9)

print("\nALL PASS" if ok else "\nSOME CHECKS FAILED")
sys.exit(0 if ok else 1)
