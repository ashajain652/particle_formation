"""Sub-plan 18, Task 9, Step 2 (spec 12.2): the box against the fine box.

A box heated by a uniform 2 MW/m2 on its x = 0 face only, radiation off, Scheil material, 700 K at the start, 2.0 s.
With insulated sides this is one-dimensional. Liquid is removed as soon as it appears (after every step the script books
each patch's film and deep liquid into removed_mass / removed_enthalpy at film_enthalpy() and zeroes them), so the
recession is the removed mass per unit front area over the metal density.

  coarse box with skins : mesh.box_mesh(8e-3, 10e-3, 10e-3, 2e-3), surface_model "skin", presets dev (20 um) and
                          production (10 um), macro steps 0.05, 0.1, 0.25, 0.5 s, sub-step 0.010 s
  fine box, no skins    : mesh.box_mesh(8e-3, 2e-3, 2e-3, 0.1e-3), surface_model "elements", macro step 0.005 s

Usage (the model is imported from the frozen measurement copy, read-only):
    python box_ablation.py --list
    python box_ablation.py fine skin_dev_dt0.5 ...        # run the named runs, one JSON each under runs/
    python box_ablation.py --collect                      # runs/*.json -> results.json, summary.md, box_ablation.png
"""
import argparse
import glob
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
MEAS = os.path.normpath(os.path.join(HERE, "..", "..", "meas"))
sys.dont_write_bytecode = True
sys.path.insert(0, MEAS)

import numpy as np

from reentry_model import body, heating, material, mesh, skin, thermal

RUNS = os.path.join(HERE, "runs")
Q_FRONT = 2.0e6            # W/m2 on the x = 0 face (the spec's value; variants override it per run)
T_START = 700.0            # K
T_END = 2.0                # s
RULE_FRONT099 = "front099"   # the controller's ruling: q on the patches with n_x < -0.99, zero elsewhere
RULE_PROJECTED = "projected" # q proportional to max(0, -n_x) on every patch, scaled each step so the total power is exactly
                             # q times the cross-section: the same power on a front that is stepped and not a flat plane
SUBSTEP = 0.010            # s, the skins' sub-step
RHO_METAL = 2813.0         # kg/m3, the density the recession is read with
COMPARE_TIMES = (0.5, 1.0, 1.5, 2.0)
PRESETS = {"dev": skin.SKIN_PRESETS["dev"], "production": skin.SKIN_PRESETS["production"]}
MACRO = (0.05, 0.1, 0.25, 0.5)
COARSE = (8e-3, 10e-3, 10e-3, 2e-3)
FINE = (8e-3, 2e-3, 2e-3, 0.1e-3)
FINE_DT = 0.005
SIDE_TOL = 1e-9            # m: how near a box plane a patch centroid must be to lie on that side face


def run_specs(tag="spec", q=Q_FRONT, rule=RULE_FRONT099, t_end=T_END, fine_box=FINE):
    """The nine runs of one variant. tag "spec" is the controller's measurement (file names without a prefix)."""
    common = dict(tag=tag, q=q, rule=rule, t_end=t_end)
    specs = {"fine": dict(box=fine_box, model="elements", dt=FINE_DT, preset=None, **common)}
    for preset in ("dev", "production"):
        for dt in MACRO:
            specs["skin_{}_dt{:g}".format(preset, dt)] = dict(box=COARSE, model="skin", dt=dt, preset=preset, **common)
    for dt in MACRO:     # context, not asked for: the coarse box without skins, to show what the skins change
        specs["elem_dt{:g}".format(dt)] = dict(box=COARSE, model="elements", dt=dt, preset=None, **common)
    for hh in (1.0e-3, 0.5e-3):    # context, not asked for: the same skins on a finer 3D mesh (is the lag the 3D mesh's?)
        specs["diag_skin_dev_dt0.05_h{:g}mm".format(hh * 1e3)] = dict(box=(COARSE[0], COARSE[1], COARSE[2], hh), model="skin",
                                                                       dt=0.05, preset="dev", **common)
    return specs


def file_name(spec, name):
    return name if spec["tag"] == "spec" else "{}__{}".format(spec["tag"], name)


def classify(b, box):
    """Boolean masks over the current surface's patches: the front (n_x < -0.99), the interior front by the controller's
    ruling (the patch's owner element owns no patch on a side face, i.e. a y or z boundary of the box), and a stricter
    interior (patch centroid more than one coarse cell from every side plane)."""
    lx, ly, lz, h = box
    s = b.surface
    front = s.normals[:, 0] < -0.99
    c, n = s.centroids, s.normals
    on_side = (((np.abs(n[:, 1]) > 0.99) & ((np.abs(c[:, 1]) < SIDE_TOL) | (np.abs(c[:, 1] - ly) < SIDE_TOL)))
               | ((np.abs(n[:, 2]) > 0.99) & ((np.abs(c[:, 2]) < SIDE_TOL) | (np.abs(c[:, 2] - lz) < SIDE_TOL))))
    side_owner = np.zeros(b.mesh.n_elements, dtype=bool)
    side_owner[s.owner[on_side]] = True
    interior = front & ~side_owner[s.owner]
    lateral = (c[:, 1] > h * (1 - 1e-9)) & (c[:, 1] < ly - h * (1 - 1e-9)) & (c[:, 2] > h * (1 - 1e-9)) & (c[:, 2] < lz - h * (1 - 1e-9))
    strict = front & lateral
    return front, interior, strict, lateral


def make_body(spec):
    lx, ly, lz, h = spec["box"]
    msh = mesh.box_mesh(lx, ly, lz, h)
    mat = material.Material.from_drama_json("AA7075_scheil")
    sk = skin.SkinSettings(cell=PRESETS[spec["preset"]], substep=SUBSTEP) if spec["model"] == "skin" else None
    settings = body.MeltSettings(removal="girin", size_feedback="initial", surface_model=spec["model"], skin=sk)
    solver = thermal.thermal_solver("skfem", linear_solver="direct")
    mass = float(mat.rho * msh.volume())
    return body.MeltingBody(msh, mat, solver, mass, settings=settings, T0=T_START, emissivity=0.0), mat


def run_one(name, spec, log=print):
    """One run; returns the record (series and totals). Exceptions are caught by the caller."""
    started = time.perf_counter()
    b, mat = make_body(spec)
    lx, ly, lz, h = spec["box"]
    dt, q_flux, rule, t_end = spec["dt"], spec["q"], spec["rule"], spec["t_end"]
    n_steps = int(round(t_end / dt))
    A_front0 = ly * lz
    A_lat0 = (ly - 2.0 * h) * (lz - 2.0 * h)               # the lateral interior's projected area
    S = {k: [] for k in ("t", "r_front_mm", "r_interior_mm", "r_strict_mm", "r_lat_mm", "r_total_mm", "T_top_interior",
                          "T_top_front", "T_top_strict", "T_top_lat", "residual", "mismatch", "repeats", "shortfall",
                          "skin_wall", "heated_fraction", "n_active", "n_skins", "step_wall", "interior_area_mm2",
                          "thickness_min_um", "newton")}
    cum = dict(front=0.0, interior=0.0, strict=0.0)       # mm of recession, summed over the steps (mass per own area)
    removed_lat = 0.0
    worst = 0.0
    prev_repeats = 0
    r0 = b.removed_mass
    for k in range(n_steps):
        t = (k + 1) * dt
        n_x = b.surface.normals[:, 0]
        if rule == RULE_PROJECTED:
            w = np.maximum(0.0, -n_x)
            q = q_flux * A_front0 * w / float((w * b.surface.areas).sum())
        else:
            q = np.where(n_x < -0.99, q_flux, 0.0)
        loads = heating.HeatingResult(q, q_flux, 0.0, 0.0, 0.0)
        heated = float((q * b.surface.areas).sum() / (q_flux * A_front0))
        t0 = time.perf_counter()
        b.advance(t, dt, loads)                                   # no flight state: nothing sprays, all liquid is film
        step_wall = time.perf_counter() - t0
        if b.consumed or b.m_f.size != b.surface.n_patches:
            raise RuntimeError("the body was consumed at t = {:.3f} s (the metal ran out)".format(t))
        # -- remove the liquid the step produced, booking its mass and its enthalpy --
        front, interior, strict, lateral = classify(b, spec["box"])
        L = b.m_f + b.m_d
        dE = L * b.film_enthalpy()
        b.removed_mass += float(L.sum())
        b.removed_enthalpy += float(dE.sum())
        areas = b.surface.areas
        for key, sel in (("front", front), ("interior", interior), ("strict", strict)):
            A = float(areas[sel].sum())
            if A > 0.0:
                cum[key] += float(L[sel].sum()) / A / RHO_METAL * 1e3
        removed_lat += float(L[lateral].sum())
        b.m_f[:] = 0.0
        b.m_d[:] = 0.0
        b.solver.set_film_mass(b._film_nodal())
        # -- record --
        res = b.energy_balance_residual()
        worst = max(worst, abs(res))
        Tp = b.surface_temperature()
        hot = lateral & (b.surface.normals[:, 0] < -0.1)          # the lateral interior's heated side
        wmean = lambda sel: float((Tp[sel] * areas[sel]).sum() / areas[sel].sum()) if sel.any() else float("nan")
        S["t"].append(t)
        S["r_front_mm"].append(cum["front"])
        S["r_interior_mm"].append(cum["interior"])
        S["r_strict_mm"].append(cum["strict"])
        S["r_lat_mm"].append(removed_lat / (A_lat0 * RHO_METAL) * 1e3)
        S["r_total_mm"].append((b.removed_mass - r0) / (A_front0 * RHO_METAL) * 1e3)
        S["T_top_interior"].append(wmean(interior))
        S["T_top_front"].append(wmean(front))
        S["T_top_strict"].append(wmean(strict))
        S["T_top_lat"].append(wmean(hot))
        S["residual"].append(res)
        S["mismatch"].append(b.last_mismatch)
        S["repeats"].append(b.interface_repeats - prev_repeats)
        prev_repeats = b.interface_repeats
        S["shortfall"].append(b.last_skin_shortfall)
        S["skin_wall"].append(b.skin_wall)
        S["heated_fraction"].append(heated)
        S["n_active"].append(int(b.mesh.n_active))
        S["n_skins"].append(int(b.skins.n) if b.skins is not None else 0)
        S["step_wall"].append(step_wall)
        S["interior_area_mm2"].append(float(areas[interior].sum() * 1e6))
        S["thickness_min_um"].append(float(b.skins.thickness().min() * 1e6) if b.skins is not None and b.skins.n else float("nan"))
        S["newton"].append(int(b.iterations[-1]))
        if (k + 1) % max(1, n_steps // 20) == 0 or k == n_steps - 1:
            log("  {} t={:.3f} r_int={:.4f} r_front={:.4f} r_lat={:.4f} r_tot={:.4f} mm T_lat={:.1f} resid={:.2e} heated={:.3f} active={} wall={:.1f}s".format(
                name, t, cum["interior"], cum["front"], S["r_lat_mm"][-1], S["r_total_mm"][-1], S["T_top_lat"][-1], res, heated,
                b.mesh.n_active, time.perf_counter() - started))
    out = {"name": name, "tag": spec["tag"], "spec": {k: (list(v) if isinstance(v, tuple) else v) for k, v in spec.items()},
           "n_steps": n_steps, "series": S, "run_time_s": time.perf_counter() - started,
           "worst_energy_residual": worst, "final_residual": S["residual"][-1],
           "interface_repeats_total": int(b.interface_repeats), "skin_wall_s": float(b.skin_wall),
           "max_shortfall_kg": float(max(S["shortfall"])), "skins_created": int(b.skins_created),
           "mass_balance": float(b.mass(0.0) - (b.mass0 - b.removed_mass)),
           "removed_mass_kg": float(b.removed_mass), "mass0_kg": float(b.mass0),
           "skin_short_events": int(b.skin_short_events), "skin_full_events": int(b.skin_full_events),
           "n_elements": int(b.mesh.n_elements), "n_nodes": int(b.mesh.n_nodes), "status": "ok"}
    return out


def save(fname, record):
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, fname + ".json"), "w") as f:
        json.dump(record, f, indent=1)


# -- collection ---------------------------------------------------------------------------------------------------
# The variants. "spec" is the controller's measurement to the letter. "q6p" is the supplementary one the spec's own numbers
# call for (see summary.md): the spec's 2 MW/m2 for 2 s never melts the box, and the n_x < -0.99 rule drops the heating the
# moment the first layer of tetrahedra dies, so the variant raises the flux and heats every patch by its projection on the
# cross-section.
VARIANTS = {
    "spec": dict(q=Q_FRONT, rule=RULE_FRONT099, ref="r_front_mm", coarse="r_interior_mm", T_coarse="T_top_interior", T_fine="T_top_front"),
    "q6p": dict(q=6.0e6, rule=RULE_PROJECTED, ref="r_total_mm", coarse="r_lat_mm", T_coarse="T_top_lat", T_fine="T_top_lat"),
}
T_PLATEAU = 905.0   # K: the top temperature has reached the feed ramp (T_feed = 910 K) once it is at or above this
RECESSION_KEYS = ("r_interior_mm", "r_front_mm", "r_strict_mm", "r_lat_mm", "r_total_mm")


def at(series_t, series_v, times):
    return [float(np.interp(tt, series_t, series_v)) for tt in times]


def drawdown(T):
    """Largest fall of a series below its running maximum [K], and the number of slope sign changes above 0.5 K/step."""
    T = np.asarray(T, dtype=float)
    T = T[np.isfinite(T)]
    if T.size < 2:
        return 0.0, 0
    peak = np.maximum.accumulate(T)
    d = np.diff(T)
    big = d[np.abs(d) > 0.5]
    flips = int((np.sign(big[1:]) != np.sign(big[:-1])).sum()) if big.size > 1 else 0
    return float((peak - T).max()), flips


def load_variant(tag):
    recs = {}
    for path in sorted(glob.glob(os.path.join(RUNS, "*.json"))):
        with open(path) as f:
            r = json.load(f)
        if r.get("tag", r.get("spec", {}).get("tag", "spec")) != tag:
            continue
        recs[r["name"]] = r
    return recs


def collect_variant(tag, recs):
    v = VARIANTS[tag]
    if "fine" not in recs or recs["fine"]["status"] != "ok":
        return None
    fine = recs["fine"]["series"]
    ft = np.concatenate([[0.0], fine["t"]])
    fr = np.concatenate([[0.0], fine[v["ref"]]])
    fr_at = at(ft, fr, COMPARE_TIMES)
    out = {"tag": tag, "q_W_m2": v["q"], "rule": v["rule"], "ref_series": v["ref"], "coarse_series": v["coarse"],
           "compare_times_s": list(COMPARE_TIMES), "fine_recession_mm_at_times": fr_at,
           "fine_time_to_0.05mm_s": float(np.interp(0.05, fr, ft)) if fr[-1] >= 0.05 else None, "runs": {}}
    for name, r in recs.items():
        if r["status"] != "ok":
            out["runs"][name] = {"status": "failed", "error": r.get("error")}
            continue
        S = r["series"]
        t = np.concatenate([[0.0], S["t"]])
        o = {"status": "ok", "dt": r["spec"]["dt"], "preset": r["spec"]["preset"], "n_steps": r["n_steps"]}
        for key in RECESSION_KEYS:
            if key not in S:       # the spec-variant fine box was run by the first version of the script, before r_lat existed
                continue
            val = np.concatenate([[0.0], S[key]])
            rec = at(t, val, COMPARE_TIMES)
            rel = [(a - b) / b if b > 0.05 else None for a, b in zip(rec, fr_at)]    # only where the fine box has >= 0.05 mm
            o[key + "_at_times"] = rec
            o[key + "_rel_diff_at_times"] = rel
            used = [abs(x) for x in rel if x is not None]
            o[key + "_max_abs_rel_diff"] = max(used) if used else None
            tt = t[1:]
            fi = np.interp(tt, ft, fr)
            m = fi >= 0.05
            o[key + "_max_abs_rel_diff_curve"] = float(np.max(np.abs(val[1:][m] - fi[m]) / fi[m])) if m.any() else None
            o[key + "_max_abs_diff_curve_mm"] = float(np.max(np.abs(val[1:] - fi)))
        Tkey = v["T_fine"] if name == "fine" else v["T_coarse"]
        T = np.array(S[Tkey])
        dd, flips = drawdown(T)
        o["T_top_series"] = Tkey
        o["T_top_final_K"] = float(T[np.isfinite(T)][-1]) if np.isfinite(T).any() else None
        o["T_top_max_K"] = float(np.nanmax(T)) if np.isfinite(T).any() else None
        o["T_top_max_drawdown_K"] = dd
        o["T_top_slope_flips"] = flips
        o["T_top_nan_steps"] = int((~np.isfinite(T)).sum())
        tt_s = np.array(S["t"])
        o["T_top_at_times"] = [float(np.interp(x, tt_s, T)) if np.isfinite(T).all() else None for x in COMPARE_TIMES]
        above = np.flatnonzero(T >= T_PLATEAU)
        monotone_to_onset = bool((np.diff(T[: above[0] + 1] if above.size else T) >= -1e-9).all())
        o["T_top_monotone_before_onset"] = monotone_to_onset
        if above.size:
            Tp = T[above[0]:]
            o["T_top_onset_s"] = float(tt_s[above[0]])
            o["T_plateau_min_K"], o["T_plateau_max_K"] = float(Tp.min()), float(Tp.max())
            o["T_plateau_mean_K"], o["T_plateau_std_K"] = float(Tp.mean()), float(Tp.std())
            o["T_plateau_n"] = int(Tp.size)
        else:
            o["T_top_onset_s"] = None
        o["worst_energy_residual"] = r["worst_energy_residual"]
        o["final_energy_residual"] = r["final_residual"]
        o["mass_balance"] = r["mass_balance"]
        o["repeats_total"] = r["interface_repeats_total"]
        o["repeats_per_step"] = r["interface_repeats_total"] / r["n_steps"]
        o["max_shortfall_kg"] = r["max_shortfall_kg"]
        o["max_mismatch_J"] = float(max(S["mismatch"]))
        o["run_time_s"] = r["run_time_s"]
        o["skin_wall_s"] = r["skin_wall_s"]
        o["n_active_final"] = S["n_active"][-1]
        o["n_skins_final"] = S["n_skins"][-1]
        o["heated_fraction_min"] = float(min(S["heated_fraction"]))
        o["heated_fraction_mean"] = float(np.mean(S["heated_fraction"]))
        o["interior_area_mm2_min"] = float(min(S["interior_area_mm2"]))
        o["interior_area_mm2_max"] = float(max(S["interior_area_mm2"]))
        o["skin_short_events"] = r["skin_short_events"]
        o["skin_full_events"] = r["skin_full_events"]
        o["removed_mass_kg"] = r["removed_mass_kg"]
        out["runs"][name] = o
    return out


def collect():
    results, all_recs = {}, {}
    for tag in VARIANTS:
        recs = load_variant(tag)
        all_recs[tag] = recs
        if recs:
            res = collect_variant(tag, recs)
            if res is not None:
                results[tag] = res
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=1)
    return all_recs, results


def fmt(x, spec="{:.3f}"):
    return "n/a" if x is None else spec.format(x)


def verdict(tag, name, o, res):
    """A plain sentence per run against the criterion (recession within 5 % at every compared time; no oscillation)."""
    v = VARIANTS[tag]
    fine = res["runs"]["fine"]
    if name == "fine":
        return "Reference."
    key = v["coarse"]
    rel = o.get(key + "_rel_diff_at_times")
    if rel is None:
        return "No recession series."
    Tdesc = ""
    if o.get("T_top_onset_s") is not None:
        Tdesc = ("Top temperature rises {} to {:.0f} K at t = {:.2f} s, then sits between {:.0f} and {:.0f} K (standard deviation {:.1f} K; "
                 "the fine box: {:.1f} K).".format("monotonically" if o["T_top_monotone_before_onset"] else "non-monotonically", T_PLATEAU,
                 o["T_top_onset_s"], o["T_plateau_min_K"], o["T_plateau_max_K"], o["T_plateau_std_K"], fine.get("T_plateau_std_K", float("nan"))))
    else:
        Tdesc = "Top temperature rises {} to {:.1f} K at 2 s (the fine box: {:.1f} K).".format(
            "monotonically" if o["T_top_monotone_before_onset"] else "non-monotonically", o["T_top_final_K"], fine["T_top_final_K"])
    if all(x is None for x in rel):
        return "Not assessable: the fine box has removed under 0.05 mm at every compared time. " + Tdesc
    used = [(t, x) for t, x in zip(COMPARE_TIMES, rel) if x is not None]
    worst_t, worst = max(used, key=lambda tx: abs(tx[1]))
    good = [t for t, x in used if abs(x) <= 0.05]
    if len(good) == len(used):
        head = "Meets the 5 % recession criterion (largest difference {:+.1f} % at t = {:g} s).".format(100 * worst, worst_t)
    else:
        head = "Fails the 5 % recession criterion: largest difference {:+.1f} % at t = {:g} s; within 5 % only at {}.".format(
            100 * worst, worst_t, ", ".join("{:g} s".format(t) for t in good) if good else "none of the compared times")
    return head + " " + Tdesc


def write_summary(all_recs, results):
    lines = ["# Box against the fine box (sub-plan 18, Task 9, Step 2): tables", "",
             "Generated by `box_ablation.py --collect` from `runs/`. Recession is the removed mass per unit front area over 2813 kg/m3, in mm. "
             "Relative differences are (coarse - fine) / fine, the fine curve interpolated to t = 0.5, 1.0, 1.5 and 2.0 s, and are given "
             "only where the fine box has removed at least 0.05 mm (below that the quotient divides by nothing).", ""]
    for tag, res in results.items():
        v = VARIANTS[tag]
        lines.append("## Variant `{}`: {:.0e} W/m2, heating rule `{}`".format(tag, v["q"], v["rule"]))
        lines.append("")
        lines.append("Coarse recession read from `{}`, fine from `{}`. Fine box recession at the four times (mm): {}.".format(
            v["coarse"], v["ref"], ", ".join("{:.4f}".format(x) for x in res["fine_recession_mm_at_times"])))
        lines.append("")
        lines.append("| run | recession (mm) at 0.5 / 1.0 / 1.5 / 2.0 s | rel. diff at those times | largest abs rel. diff (4 times) | largest abs diff over the whole curve (mm) | worst energy residual | repeats per step | max shortfall (kg) | top-T drawdown (K) | run time (s) |")
        lines.append("|" + "---|" * 10)
        for name in sorted(res["runs"], key=lambda n: (n != "fine", n)):
            o = res["runs"][name]
            if o["status"] != "ok":
                lines.append("| {} | FAILED: {} | | | | | | | | |".format(name, o.get("error")))
                continue
            key = v["ref"] if name == "fine" else v["coarse"]
            rel = o[key + "_rel_diff_at_times"]
            lines.append("| {} | {} | {} | {} | {} | {:.1e} | {:.2f} | {:.1e} | {:.2f} | {:.0f} |".format(
                name, " / ".join("{:.4f}".format(x) for x in o[key + "_at_times"]),
                "reference" if name == "fine" else " / ".join("n/a" if x is None else "{:+.1f}%".format(100 * x) for x in rel),
                "-" if name == "fine" or o[key + "_max_abs_rel_diff"] is None else "{:.1f}%".format(100 * o[key + "_max_abs_rel_diff"]),
                "-" if name == "fine" else "{:.4f}".format(o[key + "_max_abs_diff_curve_mm"]),
                o["worst_energy_residual"], o["repeats_per_step"], o["max_shortfall_kg"], o["T_top_max_drawdown_K"], o["run_time_s"]))
        lines.append("")
        lines.append("Top temperature of the front patches (K) at 0.5 / 1.0 / 1.5 / 2.0 s, and what it does after reaching 905 K:")
        lines.append("")
        lines.append("| run | top T at the four times (K) | after 905 K: min / max / std (K) | cumulative interface repeats | skin wall time (s) | verdict |")
        lines.append("|" + "---|" * 6)
        for name in sorted(res["runs"], key=lambda n: (n != "fine", n)):
            o = res["runs"][name]
            if o["status"] != "ok":
                continue
            Tt = " / ".join("n/a" if x is None else "{:.1f}".format(x) for x in o["T_top_at_times"])
            pl = "never reached" if o.get("T_top_onset_s") is None else "{:.1f} / {:.1f} / {:.1f}".format(
                o["T_plateau_min_K"], o["T_plateau_max_K"], o["T_plateau_std_K"])
            lines.append("| {} | {} | {} | {} | {:.2f} | {} |".format(name, Tt, pl, o["repeats_total"], o["skin_wall_s"], verdict(tag, name, o, res)))
        lines.append("")
    with open(os.path.join(HERE, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


def plot(all_recs):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    tags = [t for t in VARIANTS if all_recs.get(t)]
    if not tags:
        return
    fig, ax = plt.subplots(len(tags), 2, figsize=(13, 5 * len(tags)), squeeze=False)
    for i, tag in enumerate(tags):
        v = VARIANTS[tag]
        for name, r in sorted(all_recs[tag].items()):
            if r["status"] != "ok":
                continue
            S = r["series"]
            t = [0.0] + S["t"]
            if name == "fine":
                ax[i][0].plot(t, [0.0] + S[v["ref"]], "k-", lw=2.5, label="fine box (elements, 0.005 s)")
                ax[i][1].plot(S["t"], S[v["T_fine"]], "k-", lw=2.5, label="fine box")
            else:
                ls = "-" if "dev" in name else "--"
                ax[i][0].plot(t, [0.0] + S[v["coarse"]], ls, lw=1.2, label=name)
                ax[i][1].plot(S["t"], S[v["T_coarse"]], ls, lw=1.2, label=name)
        ax[i][0].set_title("variant {} ({:.0e} W/m2, {})".format(tag, v["q"], v["rule"]))
        ax[i][0].set_xlabel("time (s)"); ax[i][0].set_ylabel("recession (mm)"); ax[i][0].legend(fontsize=7)
        ax[i][1].set_xlabel("time (s)"); ax[i][1].set_ylabel("mean top temperature (K)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "box_ablation.png"), dpi=130)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*", help="run names (see --list); 'all' for every run, 'coarse' for all but the fine box")
    ap.add_argument("--variant", default="spec", choices=sorted(VARIANTS))
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--collect", action="store_true")
    a = ap.parse_args()
    v = VARIANTS[a.variant]
    specs = run_specs(a.variant, v["q"], v["rule"])
    if a.list:
        print("\n".join(specs))
        return
    if a.collect:
        all_recs, results = collect()
        write_summary(all_recs, results)
        plot(all_recs)
        for tag, res in results.items():
            print(tag, {n: (o.get("r_interior_mm_max_abs_rel_diff") if tag == "spec" else o.get("r_lat_mm_max_abs_rel_diff"), o.get("worst_energy_residual")) for n, o in res["runs"].items()})
        return
    names = list(specs) if a.names == ["all"] else [n for n in specs if n != "fine" and not n.startswith("diag_")] if a.names == ["coarse"] else a.names
    for name in names:
        if name not in specs:
            raise SystemExit("unknown run " + name)
        print("== {} / {} ==".format(a.variant, name), flush=True)
        try:
            rec = run_one(name, specs[name], log=lambda s: print(s, flush=True))
        except Exception as exc:                                  # a failed run is recorded; the others continue
            traceback.print_exc()
            rec = {"name": name, "tag": a.variant, "status": "failed", "error": "{}: {}".format(type(exc).__name__, exc),
                   "spec": {k: (list(x) if isinstance(x, tuple) else x) for k, x in specs[name].items()}}
        save(file_name(specs[name], name), rec)


if __name__ == "__main__":
    main()
