"""Evaluate the re-entry model's current surface pressure, wall shear and convective heat flux around a sphere at the
four Holman and Boyd (2011) conditions, for comparison with their DSMC results.

The repository is imported read-only from its path (no file in it is changed; bytecode writing is off). The freestream
is injected by building the model's own `atmosphere.Freestream` and `trajectory.AeroState` objects directly, exactly as
`trajectory.Simulator.aero_state` builds them, so no altitude lookup and no monkeypatch is involved. The post-shock
state is the model's own Cantera equilibrium construction (inside SurfaceFlow and PhysicsHeating).

Models and options are the defaults a `reentry_model run --thermal fem` flight builds in `cli.build_thermal`:
  PhysicsHeating(stagnation="fay-riddell", bridging="matting", matting_n=1.0, accommodation=0.8, catalycity=1.0)
  SurfaceFlow(rarefied_shear="slip", gamma_pm=1.15, kn_body_shock=0.01)  (bridging = aero.SesamTable(), kn_body_fm = 10)

Outputs (next to this script):
  model_loads.csv              case_kn, theta_deg, Cp, Ctau, Cq (0..180 deg every 1 deg)
  model_loads_diagnostics.csv  per angle: dimensional p_w, tau_w, q_w, the shear sub-branch, closure, Kn_local, edge state
  model_loads_summary.json     per case: Knudsen and Mach numbers, branch, stagnation pressure and heating, DSMC peak
"""
import csv
import json
import math
import os
import sys

sys.dont_write_bytecode = True
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # Reentry Model/
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from reentry_model.shared import aero, atmosphere, heating, surface_flow, trajectory  # noqa: E402
from reentry_model.shared.constants import M_BAR_AIR  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# Holman and Boyd (2011) conditions
T_INF, V_INF, DIAMETER, T_WALL = 200.0, 7108.0, 0.3048, 1500.0
RADIUS = 0.5 * DIAMETER
R_GAS = 287.05
CASES = [(0.002, 1.007e-4), (0.01, 2.014e-5), (0.05, 3.987e-6), (0.25, 8.057e-7)]
DSMC_PEAK_HEATING = {0.002: 9.30e5, 0.01: 8.04e5, 0.05: 3.54e5, 0.25: 1.05e5}   # W/m2, Holman and Boyd's DSMC
MATERIAL_DENSITY = 2813.0      # kg/m3, the CLI's default (drama-AA7075); only sets the drag deceleration, which enters G
RHO_LIQUID = 2400.0            # kg/m3, nominal; only enters the film's inertial term G, not p_w, tau_w or q_w

THETA_DEG = np.arange(0.0, 180.0 + 0.5, 1.0)
THETA = np.radians(THETA_DEG)


def aero_state(rho):
    """The model's AeroState for this freestream, built the way trajectory.Simulator.aero_state builds it."""
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


def shear_branch(flow, kn_local, theta):
    """The rule SurfaceFlow.evaluate used for tau_w at each angle."""
    out = []
    for kl, th in zip(kn_local, theta):
        if th > 0.5 * math.pi:
            out.append("leeward: tau = 0, p = p_inf")
        elif flow.branch == surface_flow.BRANCH_FREE_MOLECULAR:
            out.append("free molecular (Schaaf-Chambre)")
        elif flow.branch == surface_flow.BRANCH_MERGED:
            out.append("merged: (1-f) tau_slip + f tau_fm")
        elif kl >= surface_flow.KN_LOCAL_SLIP:
            out.append("shock layer: free-molecular tau (Kn_local >= 0.1)")
        elif kl >= surface_flow.KN_LOCAL_CONTINUUM:
            out.append("shock layer: slip tau, Couette closure (0.01 <= Kn_local < 0.1)")
        else:
            out.append("shock layer: continuum tau/(1+Kn_local), Girin closure (Kn_local < 0.01)")
    return out


def main():
    flow_model = surface_flow.SurfaceFlow(rarefied_shear="slip", gamma_pm=surface_flow.GAMMA_PM,
                                          kn_body_shock=surface_flow.KN_BODY_SHOCK)
    heat_model = heating.PhysicsHeating(stagnation="fay-riddell", bridging="matting", matting_n=1.0, accommodation=0.8,
                                        catalycity=1.0)
    T_w = np.full(THETA.size, T_WALL)
    rows, diag, summary = [], [], {}
    for label, rho in CASES:
        state = aero_state(rho)
        fs = state.freestream
        q_inf = 0.5 * rho * V_INF ** 2
        flow = flow_model.evaluate(state, THETA, RADIUS, RHO_LIQUID, T_wall=T_w)
        heat = heat_model.evaluate(state, THETA, T_w, RADIUS, T_mean=T_WALL)
        Cp = (flow.p_w - fs.p) / q_inf
        Ctau = flow.tau / q_inf
        Cq = heat.q_conv / (0.5 * rho * V_INF ** 3)
        branches = shear_branch(flow, flow.kn_local, THETA)
        for k, th in enumerate(THETA_DEG):
            rows.append((label, th, Cp[k], Ctau[k], Cq[k]))
            diag.append({"case_kn": label, "theta_deg": th, "p_w_Pa": flow.p_w[k], "tau_w_Pa": flow.tau[k],
                         "tau_continuum_Pa": flow.tau_continuum[k], "tau_fm_Pa": flow.tau_fm[k], "q_w_Wm2": heat.q_conv[k],
                         "body_branch": flow.branch_name, "shear_rule": branches[k],
                         "closure": "girin" if flow.closure[k] == surface_flow.CLOSURE_GIRIN else "couette",
                         "kn_local": flow.kn_local[k], "u_e_ms": flow.u_e[k], "T_e_K": flow.T_e[k], "rho_e_kgm3": flow.rho_e[k],
                         "delta_a_m": flow.delta_a[k]})
        # stagnation fluxes and the post-shock state, recomputed with the model's own functions for the record
        q_c, q_fm, h_s = heat_model.stagnation_fluxes(state, T_WALL, RADIUS)
        stag = heat_model.air.stagnation(rho, T_INF, V_INF)
        wind = THETA <= 0.5 * math.pi
        floor_hit = int(np.sum(wind & (flow.p_w <= fs.p * (1 + 1e-12)) & (THETA > 0)))
        i90 = int(np.argmin(np.abs(THETA_DEG - 90.0)))
        summary[str(label)] = {
            "rho_inf_kgm3": rho, "p_inf_Pa": fs.p, "q_inf_Pa": q_inf, "T_inf_K": T_INF, "V_ms": V_INF,
            "kn_body_model_hard_sphere_d3.65A": state.kn, "mach_model_gamma1.4": state.ma, "cd_model": state.cd,
            "body_branch": flow.branch_name, "re_shock_Re2": flow.re_shock,
            "shear_bridging_f_sesam_table": flow_model.bridging(state.kn),
            "p_stag_Pa": flow.p_stag, "Cp_stag": (flow.p_stag - fs.p) / q_inf,
            "phi_sonic_deg": math.degrees(flow.phi_sonic) if np.isfinite(flow.phi_sonic) else None,
            "p_w_90deg_Pa": float(flow.p_w[i90]), "Cp_90deg": float(Cp[i90]),
            "p_floor_triggered_windward_angles": floor_hit,
            "post_shock_T_K": stag.T, "post_shock_rho_kgm3": stag.rho, "post_shock_h_Jkg": stag.h,
            "post_shock_hD_Jkg": stag.h_D, "post_shock_X": stag.X,
            "kn_local_stag": float(flow.kn_local[0]), "kn_local_45deg": float(flow.kn_local[45]),
            "kn_local_90deg": float(flow.kn_local[i90]),
            "q_stag_Wm2": float(heat.q_conv[0]), "q_stag_bridged_Wm2": heat.q_stag,
            "q_stag_continuum_fay_riddell_Wm2": heat.q_stag_c, "q_stag_free_molecular_Wm2": heat.q_stag_fm,
            "q_stag_check_c_fm": [q_c, q_fm], "heating_fm_weight_w": heat.blend,
            "dsmc_peak_heating_Wm2": DSMC_PEAK_HEATING[label],
            "model_over_dsmc": float(heat.q_conv[0]) / DSMC_PEAK_HEATING[label],
            "tau_max_Pa": float(flow.tau.max()), "theta_tau_max_deg": float(THETA_DEG[int(np.argmax(flow.tau))]),
            "shear_rules_used": sorted(set(branches)),
        }
    with open(os.path.join(HERE, "model_loads.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case_kn", "theta_deg", "Cp", "Ctau", "Cq"])
        for r in rows:
            w.writerow([r[0], "{:g}".format(r[1])] + ["{:.8e}".format(x) for x in r[2:]])
    with open(os.path.join(HERE, "model_loads_diagnostics.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(diag[0]))
        w.writeheader()
        for d in diag:
            w.writerow({k: ("{:.8e}".format(v) if isinstance(v, float) and k not in ("case_kn", "theta_deg") else
                            ("{:g}".format(v) if isinstance(v, float) else v)) for k, v in d.items()})
    with open(os.path.join(HERE, "model_loads_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    for label, s in summary.items():
        print("Kn {:>5}: Kn_body {:.4g} Ma {:.2f} branch {:<13} f {:.3f} p0 {:.4g} Pa Cp0 {:.4f} phi* {} "
              "q_stag {:.4g} (FR {:.4g}, FM {:.4g}, w {:.3f}) DSMC {:.3g} ratio {:.3f} Kn_loc(0,90) {:.2e},{:.2e} floor {}".format(
                  label, s["kn_body_model_hard_sphere_d3.65A"], s["mach_model_gamma1.4"], s["body_branch"],
                  s["shear_bridging_f_sesam_table"], s["p_stag_Pa"], s["Cp_stag"],
                  None if s["phi_sonic_deg"] is None else round(s["phi_sonic_deg"], 2), s["q_stag_Wm2"],
                  s["q_stag_continuum_fay_riddell_Wm2"], s["q_stag_free_molecular_Wm2"], s["heating_fm_weight_w"],
                  s["dsmc_peak_heating_Wm2"], s["model_over_dsmc"], s["kn_local_stag"], s["kn_local_90deg"],
                  s["p_floor_triggered_windward_angles"]))
        print("          rules:", s["shear_rules_used"], "Re2 {:.4g} tau_max {:.4g} Pa at {} deg".format(
            s["re_shock_Re2"], s["tau_max_Pa"], s["theta_tau_max_deg"]))


if __name__ == "__main__":
    main()
