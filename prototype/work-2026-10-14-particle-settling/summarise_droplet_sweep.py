"""Post-process the droplet sweep: add rough vaporisation brackets to every summary, write the
combined matrix table, and draw the two report figures."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import droplet_heating_sweep as dh
from reentry_model.atmosphere import US76TableAtmosphere

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(dh.HERE, "runs", "droplet_sweep")
atm = US76TableAtmosphere()
frames = []
for f in sorted(os.listdir(OUT)):
    if not f.startswith("summary_eps") or f.endswith("_with_vapour.csv"):
        continue
    s = pd.read_csv(os.path.join(OUT, f))
    tsat = np.array([dh.al_saturation_temperature(atm.state(0, h * 1e3, 0, 0).p) for h in s.h0_km])
    s["T_sat_Al_release_K"] = tsat
    # lower bracket: only the sensible heat above saturation at the peak goes into vapour
    s["vap_frac_lower"] = np.clip(dh.CP_LIQUID * (s.T_peak_K - tsat) / dh.H_VAP_AL, 0, 1)
    # upper bracket: everything absorbed beyond heating to saturation goes into vapour, none radiated
    s["vap_frac_upper"] = np.clip((s.E_in_Jkg - dh.CP_LIQUID * (tsat - dh.T_RELEASE)) / dh.H_VAP_AL, 0, 1)
    s["flag_curvature_matters"] = s.curvature_drop_at_stop_m > 300.0
    s["physical_status"] = np.where(s.flag_T_above_Al_boiling_1atm, "UNPHYSICAL: above 1-atm Al boiling point",
                           np.where(s.flag_T_above_Al_saturation_local, "UNPHYSICAL: above local Al saturation",
                                    "no Al vaporisation expected"))
    s.to_csv(os.path.join(OUT, f.replace(".csv", "_with_vapour.csv")), index=False)
    frames.append(s)
allv = pd.concat(frames)
b = allv[(allv.eps == dh.EPS_DEFAULT) & (allv.alpha == dh.ALPHA_DEFAULT)]
cols = ["d_um", "h0_km", "V0_ms", "gamma0_deg", "T_peak_K", "dT_peak_K", "E_in_Jkg", "E_rad_Jkg",
        "frac_E_in_s_gt3", "peak_decel_g", "t_e_s", "x_e_km", "t_stop_s", "x_stop_km", "path_stop_km",
        "x_ground_km", "t_ground_days", "t_cool_below_liquidus_s", "Kn_min", "Kn_max", "Kn_energy_weighted",
        "frac_E_in_Kn_lt10", "s_release", "s_min", "blend_warning", "flag_T_above_Al_boiling_1atm",
        "flag_T_above_Al_saturation_local", "T_sat_Al_release_K", "vap_frac_lower", "vap_frac_upper",
        "flag_curvature_matters", "physical_status"]
b[cols].sort_values(["d_um", "h0_km"]).to_csv(os.path.join(OUT, "matrix_summary_baseline.csv"), index=False)

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for (eps, al), g in allv.groupby(["eps", "alpha"]):
    g = g[g.d_um == 10].sort_values("h0_km")
    ax[0].plot(g.h0_km, g.T_peak_K, "o-", ms=3, label="emissivity {:g}, accommodation {:g}".format(eps, al))
g = b[b.d_um == 10].sort_values("h0_km")
ax[0].plot(g.h0_km, g.T_sat_Al_release_K, "k--", label="Al saturation temperature at local pressure")
ax[0].axhline(dh.T_BOIL_AL, color="r", ls=":", label="Al boiling point at 1 atm (2743 K)")
ax[0].set_xlabel("release altitude [km]"); ax[0].set_ylabel("peak droplet temperature [K]")
ax[0].set_title("Peak temperature (10 um; every diameter within 1 %)"); ax[0].legend(fontsize=8)
for h in (90, 80, 70, 60, 50):
    r = pd.read_csv(os.path.join(OUT, "histories", dh.history_name(10, h)))
    ax[1].semilogx(r.t_s.clip(lower=1e-7), r.T_K, label="released at {} km".format(h))
ax[1].set_xlabel("time since release [s]"); ax[1].set_ylabel("droplet temperature [K]")
ax[1].set_title("10 um droplet temperature history (baseline)"); ax[1].legend(fontsize=8)
ax[1].set_xlim(1e-6, 1e2)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "droplet_peak_temperature.png"), dpi=140)
print(b.groupby("physical_status").size())
print(b.pivot(index="h0_km", columns="d_um", values="vap_frac_lower").iloc[::2, [0, 4, 8]].round(3))
print(b.pivot(index="h0_km", columns="d_um", values="vap_frac_upper").iloc[::2, [0, 4, 8]].round(3))
print("curvature flag", b.flag_curvature_matters.sum(), b[b.flag_curvature_matters][["d_um", "h0_km"]].values.tolist())
