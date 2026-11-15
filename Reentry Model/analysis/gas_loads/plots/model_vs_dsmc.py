"""Surface loads round the whole sphere (0-180 deg) from the re-entry model itself, on the interim gas loads, against
Holman & Boyd's (2011) DSMC and the reference calculation they were built from, at Holman & Boyd's conditions
(0.3048 m sphere, Mach 25 air at 7108 m/s, 200 K freestream, 1500 K wall; freestream Knudsen numbers 0.002, 0.01,
0.05, 0.25, which are Kn_body 0.0026, 0.013, 0.067 and 0.33 in the model's own hard-sphere definition).

The model's loads are what a default `--thermal fem --heating physics` flight builds in `cli.build_thermal`: one
`shared.gas_loads.GasLoads(gamma_pm=1.15)` (2-degree bins) shared by `SurfaceFlow(rarefied_shear="slip",
gamma_pm=1.15, kn_body_shock=0.01, gas_loads=...)` and `heating.SimilarityHeating(..., bridging="matting",
matting_n=1, accommodation=0.8)`, evaluated on the model's own `AeroState` (built as `trajectory.Simulator.aero_state`
builds it, as in `reference/model_loads.py`) at every degree with a uniform 1500 K wall. Below the shear switch
(Kn_body 0.05) the model's shear is the similarity shear through the model's slip factor and, on the merged branch
(Kn_body 0.01-0.05), the SESAM-table blend toward free-molecular shear; the reference line there is the continuum
similarity shear before them.

Lines: Holman & Boyd DSMC (blue; reference/hb_dsmc.csv, their pressure coefficient read as p/(0.5 rho U^2), design A1);
the reference calculation (black dotted): the joined pressure of reference/similarity_lee.csv, the shear and bridged
heating of reference/similarity_heat.csv below the switch and Daniels & Schwartzentruber's shear (reference/ds_model.csv)
above it; the model on the interim gas loads (green); and the model's classic loads, `--gas-loads classic`
(reference/model_loads.csv, made before the change with the same code the classic option keeps; green dash-dot).
Pressure is p / p_inf; C_tau = tau/(0.5 rho U^2), Cq = q/(0.5 rho U^3). Writes output/gas_loads/model_vs_dsmc.png
and prints the ratios quoted in docs/model_assumptions.md section 14.6, and the shear jump at the Kn_body 0.05 switch (design A4b).
"""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]                                     # Reentry Model/
sys.path.insert(0, str(REPO))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from reentry_model.shared import aero, atmosphere, gas_loads, heating, surface_flow, trajectory  # noqa: E402
from reentry_model.shared.constants import M_BAR_AIR  # noqa: E402

DATA = HERE.parent / "reference"
T_INF, V_INF, DIAMETER, T_WALL = 200.0, 7108.0, 0.3048, 1500.0
RADIUS = 0.5 * DIAMETER
R_GAS = 287.05
CASES = ((0.002, 1.007e-4), (0.01, 2.014e-5), (0.05, 3.987e-6), (0.25, 8.057e-7))     # label, rho_inf [kg/m3]
MATERIAL_DENSITY = 2813.0      # kg/m3, the CLI's default; only sets the drag deceleration, which enters G, not the loads
RHO_LIQUID = 2400.0            # kg/m3, nominal; only enters the film's inertial term G
Q_OVER_PINF = 0.5 * 1.4 * (V_INF / (1.4 * R_GAS * T_INF) ** 0.5) ** 2      # 0.5 gamma M^2 = 440
THETA_DEG = np.arange(0.0, 181.0, 1.0)
ROWS = (("p", "Pressure / freestream pressure", "log", (1e-4, 2e3)),
        ("Ctau", "Shear coefficient", "linear", None),
        ("Cq", "Heat-flux coefficient", "log", (1e-5, 1.5)))
STYLE = {"hb": ("Holman & Boyd 2011, DSMC (digitized)", dict(color="#2a78d6", ls="-", lw=2.2)),
         "ref": ("Reference calculation (analysis/gas_loads/reference)", dict(color="#0b0b0b", ls=":", lw=1.6)),
         "model": ("Re-entry model, interim gas loads (SurfaceFlow and SimilarityHeating on GasLoads)",
                   dict(color="#1baf7a", ls="-", lw=1.8)),
         "classic": ("Re-entry model, classic loads (--gas-loads classic)", dict(color="#1baf7a", ls="-.", lw=1.1,
                                                                                 alpha=0.6))}


def aero_state(rho):
    """The model's AeroState for this freestream, as `trajectory.Simulator.aero_state` builds it (reference/model_loads.py)."""
    p = rho * R_GAS * T_INF
    fs = atmosphere.Freestream(rho, T_INF, p, M_BAR_AIR, np.zeros(3))
    v_rel = np.array([V_INF, 0.0, 0.0])
    ma = aero.mach(V_INF, T_INF, M_BAR_AIR)
    kn = aero.knudsen(rho, M_BAR_AIR, DIAMETER)
    cd = aero.drag_coefficient(kn, ma, aero.SphereDragTables.from_json(), aero.bridging_by_name("sesam-table"))
    area = math.pi * RADIUS ** 2
    mass = MATERIAL_DENSITY * 4.0 / 3.0 * math.pi * RADIUS ** 3
    a_drag = -0.5 * rho * V_INF * v_rel * cd * area / mass
    return trajectory.AeroState(float("nan"), 0.0, 0.0, fs, v_rel, V_INF, kn, ma, cd, a_drag, 0.5 * rho * V_INF ** 2)


def models(provider):
    flow = surface_flow.SurfaceFlow(rarefied_shear="slip", gamma_pm=surface_flow.GAMMA_PM,
                                    kn_body_shock=surface_flow.KN_BODY_SHOCK, gas_loads=provider)
    heat = heating.SimilarityHeating(provider, bridging="matting", matting_n=1.0, accommodation=0.8)
    return flow, heat


def model_loads():
    """{label: {"p", "Ctau", "Cq": arrays on THETA_DEG}, plus the provider's meet and separation angles and Kn_body}."""
    provider = gas_loads.GasLoads(gamma_pm=surface_flow.GAMMA_PM)
    flow_model, heat_model = models(provider)
    theta, T_w = np.radians(THETA_DEG), np.full(THETA_DEG.size, T_WALL)
    out = {}
    for label, rho in CASES:
        state = aero_state(rho)
        flow = flow_model.evaluate(state, theta, RADIUS, RHO_LIQUID, T_wall=T_w)
        heat = heat_model.evaluate(state, theta, T_w, RADIUS, T_mean=T_WALL)
        L = provider.meridian(state, RADIUS, T_WALL)
        q2, q3 = 0.5 * rho * V_INF ** 2, 0.5 * rho * V_INF ** 3
        out[label] = {"p": flow.p_w / state.freestream.p, "Ctau": flow.tau / q2, "Cq": heat.q_conv / q3,
                      "kn_body": state.kn, "meet": L.meet_deg, "separation": L.separation_deg, "source": L.shear_source,
                      "q_stag": heat.q_stag, "branch": flow.branch_name}
    return out


def read(name, key, transform=lambda v: v):
    out = defaultdict(dict)
    with open(DATA / name) as f:
        for r in csv.DictReader(f):
            if r.get(key) not in ("", None):
                out[float(r["case_kn"])][float(r["theta_deg"])] = transform(float(r[key]))
    return out


def reference_curves(model):
    """{source: {label: {key: {theta: value}}}} for DSMC, the reference calculation and the classic loads."""
    hb = {"p": read("hb_dsmc.csv", "Cp", lambda cp: cp * Q_OVER_PINF), "Ctau": read("hb_dsmc.csv", "Ctau"),
          "Cq": read("hb_dsmc.csv", "Cq")}
    sim_tau, ds_tau = read("similarity_heat.csv", "Ctau_continuum"), read("ds_model.csv", "Ctau")
    ref = {"p": read("similarity_lee.csv", "p_over_pinf"), "Cq": read("similarity_heat.csv", "Cq"),
           "Ctau": {label: (ds_tau[label] if model[label]["kn_body"] >= gas_loads.KN_SHEAR_SWITCH else sim_tau[label])
                    for label, _ in CASES}}
    classic = {"p": read("model_loads.csv", "Cp", lambda cp: 1.0 + cp * Q_OVER_PINF), "Ctau": read("model_loads.csv", "Ctau"),
               "Cq": read("model_loads.csv", "Cq")}
    return {"hb": hb, "ref": ref, "classic": classic}


def plot(model, curves):
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#8a8984", "xtick.color": "#52514e", "ytick.color": "#52514e"})
    fig, axes = plt.subplots(3, 4, figsize=(13.5, 9.4), sharex=True)
    fig.patch.set_facecolor("#fcfcfb")
    for i, (key, ylabel, scale, ylim) in enumerate(ROWS):
        for j, (label, _) in enumerate(CASES):
            ax = axes[i, j]
            ax.set_facecolor("#fcfcfb")
            ax.axvline(90, color="#c3c2b7", lw=1.0)
            for name in ("hb", "classic", "ref", "model"):
                text, style = STYLE[name]
                if name == "model":
                    pts = list(zip(THETA_DEG, model[label][key]))
                else:
                    pts = sorted(curves[name][key][label].items())
                if scale == "log":
                    pts = [(t, v) for t, v in pts if v > 0]
                if pts:
                    ax.plot(*zip(*pts), label=text, **style)
            ax.set_yscale(scale)
            if ylim:
                ax.set_ylim(*ylim)
            else:
                ax.set_ylim(bottom=0)
            ax.set_xlim(0, 180)
            ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
            ax.grid(True, which="major", color="#e4e3df", lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == 0:
                ax.set_title(f"Freestream Kn {label:g} (Kn_body {model[label]['kn_body']:.2g})", fontsize=10)
            if j == 0:
                ax.set_ylabel(ylabel)
            if i == 2:
                ax.set_xlabel("Angle from stagnation point (deg)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("The re-entry model's interim gas loads round the sphere: 0.305 m sphere, Mach 25 air at 7.1 km/s, "
                 "1500 K wall (Holman & Boyd 2011)", fontsize=11, y=0.99)
    fig.text(0.5, 0.893, "Shear below the Kn_body 0.05 switch (first two columns): the model's similarity shear after its "
             "slip factor (and, at Kn_body 0.013, the SESAM-table blend);\nthe reference line is the continuum shear before "
             "them. Above the switch both are Daniels & Schwartzentruber's.", ha="center", fontsize=8.5, color="#52514e")
    fig.tight_layout(rect=(0, 0, 1, 0.885))
    out = REPO / "output" / "gas_loads" / "model_vs_dsmc.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print("saved", out)


def ratios(model, curves):
    """The model against the reference calculation (how closely the shared code reproduces it on 2-degree bins, after
    the slip factor and blend) and against DSMC (the known accuracy)."""
    for label, _ in CASES:
        m = model[label]
        print(f"Kn {label:<5} Kn_body {m['kn_body']:.4g} branch {m['branch']:<12} shear {m['source']:<24} "
              f"meet {m['meet']:g} deg separation {m['separation']:g} deg q_stag {m['q_stag']:.4g} W/m2")
        for key, angles in (("p", (0, 45, 90, 120, 150, 180)), ("Ctau", (10, 30, 45, 60, 75, 90, 100, 120, 140)),
                            ("Cq", (0, 45, 90, 120, 135, 150))):
            row_ref, row_hb = [], []
            for t in angles:
                v = m[key][int(t)]
                r, h = curves["ref"][key][label].get(float(t)), curves["hb"][key][label].get(float(t))
                row_ref.append(f"{t}:{v / r:.3f}" if r else f"{t}:-")
                row_hb.append(f"{t}:{v / h:.3g}" if h else f"{t}:-")
            print(f"    {key:<4} model/reference {' '.join(row_ref)}")
            print(f"    {key:<4} model/DSMC      {' '.join(row_hb)}")


def switch_jump(angles_deg=(30.0, 60.0, 90.0), kn_targets=(0.0495, 0.0505)):
    """The model's shear at the Kn_body 0.05 switch (design A4b): SurfaceFlow on a GasLoads at the Holman-Boyd freestream
    (200 K, 7108 m/s, 0.3048 m, 1500 K wall) with rho set so that Kn_body is just below and just above 0.05; the ratio
    above/below is the jump, not smoothed. For scale, the same two densities with the switch moved out of the way (one
    shear source at both) give the change the 2 % density step alone makes."""
    rho_ref = 3.987e-6
    kn_ref = aero.knudsen(rho_ref, M_BAR_AIR, DIAMETER)
    rhos = [rho_ref * kn_ref / k for k in kn_targets]
    theta = np.radians(np.asarray(angles_deg))
    T_w = np.full(theta.size, T_WALL)

    def taus(kn_switch):
        flow, _ = models(gas_loads.GasLoads(gamma_pm=surface_flow.GAMMA_PM, kn_switch=kn_switch))
        return [flow.evaluate(aero_state(rho), theta, RADIUS, RHO_LIQUID, T_wall=T_w).tau for rho in rhos]

    below, above = taus(gas_loads.KN_SHEAR_SWITCH)
    sim_lo, sim_hi = taus(1.0)                         # similarity shear (slip, blend) at both densities
    ds_lo, ds_hi = taus(0.0)                           # Daniels & Schwartzentruber's shear at both
    print("shear jump at the switch: Kn_body {} with rho {} kg/m3".format(
        [round(aero.knudsen(r, M_BAR_AIR, DIAMETER), 5) for r in rhos], [f"{r:.5g}" for r in rhos]))
    for k, t in enumerate(angles_deg):
        print(f"  {t:g} deg: tau below {below[k]:.4g} Pa (similarity), above {above[k]:.4g} Pa (Daniels & "
              f"Schwartzentruber), ratio above/below {above[k] / below[k]:.3f}; one source at both densities: "
              f"similarity {sim_hi[k] / sim_lo[k]:.4f}, Daniels & Schwartzentruber {ds_hi[k] / ds_lo[k]:.4f}")


def main():
    model = model_loads()
    curves = reference_curves(model)
    plot(model, curves)
    ratios(model, curves)
    switch_jump()


if __name__ == "__main__":
    main()
