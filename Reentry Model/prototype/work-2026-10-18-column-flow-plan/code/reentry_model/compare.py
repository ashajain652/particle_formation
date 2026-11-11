"""Model history vs a SESAM reference: the model sampled at the reference's own time stamps, error
metrics over the whole flight and the hypersonic phase (V_ref > 1 km/s), and six overlay plots."""
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "analysis"))
try:
    from plot_style import INK, MUTED, SECOND, apply_rcparams, strip_top_right_spines
except ImportError:                                    # analysis/ not present: fall back to plain matplotlib
    INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#898781"

    def apply_rcparams(plt, font_size=10.5):
        return None

    def strip_top_right_spines(ax):
        return None

MODEL_COLOR = "#3987e5"
REF_COLOR = INK
HYPERSONIC_MS = 1000.0
PLOT_NAMES = ("velocity_time.png", "altitude_time.png", "altitude_velocity.png", "angles_time.png",
              "ground_track.png", "regime_drag.png")
THERMAL_PLOT_NAMES = ("heating_time.png", "temperature_time.png", "integrated_heat.png")
MELT_PLOT_NAMES = ("mass_time.png", "mass_altitude.png", "mass_budget.png", "spraying_time.png", "closures_time.png",
                   "droplet_size_time.png", "size_distribution.png")
CONTINUUM_KN = 0.01
DEMISE_FRACTION = 0.01


def align(history, reference):
    t_mod = history.columns["time_s"]
    # Exclude reference points outside the model's time window to avoid extrapolation
    mask = (reference.time >= t_mod[0]) & (reference.time <= t_mod[-1])
    t = reference.time[mask]

    def at(col, scale=1.0):
        return np.interp(t, t_mod, history.columns[col]) * scale

    return {
        "t": t,
        "V_model": at("velocity_kms", 1e3), "V_ref": reference.velocity[mask],
        "h_model": at("altitude_km", 1e3), "h_ref": reference.altitude[mask],
        "gamma_model": at("flight_path_deg"), "gamma_ref": np.degrees(reference.flight_path[mask]),
        "heading_model": at("heading_deg"), "heading_ref": np.degrees(reference.heading[mask]),
        "lat_model": at("lat_deg"), "lat_ref": np.degrees(reference.lat[mask]),
        "lon_model": at("lon_deg"), "lon_ref": np.degrees(reference.lon[mask]),
        "kn_model": at("knudsen"), "kn_ref": reference.knudsen[mask],
        "cd_model": at("drag"), "cd_ref": reference.drag[mask],
    }


def _phase(a, mask):
    dV = a["V_model"][mask] - a["V_ref"][mask]
    dh = a["h_model"][mask] - a["h_ref"][mask]
    rel = np.abs(dV) / np.maximum(a["V_ref"][mask], 1.0)
    rms = lambda x: float(math.sqrt(np.mean(x * x))) if x.size else float("nan")
    return {"n_points": int(mask.sum()),
            "dV_max_ms": float(np.abs(dV).max()) if dV.size else float("nan"), "dV_rms_ms": rms(dV),
            "dV_rel_max": float(rel.max()) if rel.size else float("nan"), "dV_rel_rms": rms(rel),
            "dh_max_m": float(np.abs(dh).max()) if dh.size else float("nan"), "dh_rms_m": rms(dh)}


def _overlay(ax, x_ref, y_ref, x_model, y_model, ylabel, xlabel=None, title=None, legend=True, log_y=False):
    """Draw reference and model overlay lines with consistent styling."""
    plot_fn = ax.semilogy if log_y else ax.plot
    plot_fn(x_ref, y_ref, color=REF_COLOR, lw=1.6, label="SESAM")
    plot_fn(x_model, y_model, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title, color=SECOND, fontsize=10)
    if legend:
        ax.legend(frameon=False)
    strip_top_right_spines(ax)


def metrics(history, reference):
    a = align(history, reference)
    t_model_end, t_ref_end = float(history.columns["time_s"][-1]), float(reference.time[-1])
    v_model_final = float(history.columns["velocity_kms"][-1] * 1e3)
    v_ref_final = float(reference.velocity[-1])
    return {
        "n_points": int(a["t"].size),
        "model_end_time_s": t_model_end, "reference_end_time_s": t_ref_end,
        "d_end_time_s": t_model_end - t_ref_end, "d_end_time_rel": (t_model_end - t_ref_end) / t_ref_end,
        "final_velocity_model_ms": v_model_final,
        "final_velocity_reference_ms": v_ref_final,
        "d_final_velocity_ms": v_model_final - v_ref_final,
        "all": _phase(a, np.ones(a["t"].size, dtype=bool)),
        "hypersonic": _phase(a, a["V_ref"] > HYPERSONIC_MS),
    }


def _overlay_with_residual(a, key, ylabel, resid_label, scale, path, title):
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    _overlay(ax, a["t"], a[key + "_ref"] * scale, a["t"], a[key + "_model"] * scale, ylabel, title=title)
    rx.plot(a["t"], (a[key + "_model"] - a[key + "_ref"]), color=MODEL_COLOR, lw=1.0)
    rx.axhline(0.0, color=MUTED, lw=0.6)
    rx.set_ylabel(resid_label); rx.set_xlabel("time [s]")
    strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def plot_all(history, reference, outdir, title):
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    a = align(history, reference)
    paths = [os.path.join(outdir, n) for n in PLOT_NAMES]
    _overlay_with_residual(a, "V", "velocity [km/s]", "model - SESAM [m/s]", 1e-3, paths[0], title)
    _overlay_with_residual(a, "h", "altitude [km]", "model - SESAM [m]", 1e-3, paths[1], title)

    fig, ax = plt.subplots(figsize=(7, 5))
    _overlay(ax, a["V_ref"] / 1e3, a["h_ref"] / 1e3, a["V_model"] / 1e3, a["h_model"] / 1e3,
             ylabel="altitude [km]", xlabel="velocity [km/s]", title=title)
    fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)

    fig, (g, hd) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    _overlay(g, a["t"], a["gamma_ref"], a["t"], a["gamma_model"], ylabel="flight-path angle [deg]", title=title)
    _overlay(hd, a["t"], a["heading_ref"], a["t"], a["heading_model"], ylabel="heading [deg]", xlabel="time [s]", legend=False)
    fig.tight_layout(); fig.savefig(paths[3], dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    _overlay(ax, a["lon_ref"], a["lat_ref"], a["lon_model"], a["lat_model"],
             ylabel="latitude [deg]", xlabel="longitude [deg]", title=title)
    fig.tight_layout(); fig.savefig(paths[4], dpi=150); plt.close(fig)

    fig, (kx, cx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    _overlay(kx, a["t"], np.maximum(a["kn_ref"], 1e-7), a["t"], np.maximum(a["kn_model"], 1e-7),
             ylabel="Knudsen number", title=title, log_y=True)
    for thr in (10.0, 1.0, 0.1, 0.01):
        kx.axhline(thr, color=MUTED, lw=0.5, ls=":")
    _overlay(cx, a["t"], a["cd_ref"], a["t"], a["cd_model"],
             ylabel="C_D", xlabel="time [s]", legend=False)
    fig.tight_layout(); fig.savefig(paths[5], dpi=150); plt.close(fig)
    return paths


def has_thermal(history, reference=None):
    """True when the model history carries the coupled run's heat columns (and the reference has all three of SESAM's heat columns)."""
    ok = "convective_heat_W" in history.columns
    return ok and (reference is None or all(x is not None for x in (reference.convective_heat, reference.rad_cooling, reference.integrated_heat)))


def align_thermal(history, reference):
    """Model heat/temperature columns interpolated at the reference's time stamps (within the model's window)."""
    t_mod = history.columns["time_s"]
    mask = (reference.time >= t_mod[0]) & (reference.time <= t_mod[-1])
    t = reference.time[mask]
    at = lambda col: np.interp(t, t_mod, history.columns[col])
    return {
        "t": t, "V_ref": reference.velocity[mask], "kn_ref": reference.knudsen[mask],
        "Q_model": at("convective_heat_W"), "Q_ref": reference.convective_heat[mask],
        "H_model": at("integrated_heat_J"), "H_ref": reference.integrated_heat[mask],
        "R_model": -at("rad_cooling_W"), "R_ref": -reference.rad_cooling[mask],
        "T_model": at("temperature_K"), "T_ref": reference.temperature[mask],
        "Tmax_model": at("surface_T_max_K"), "Tmin_model": at("surface_T_min_K"),
    }


def _peak_relative(delta, ref, mask):
    """Errors normalised by the reference's peak over the mask (robust where the reference is small)."""
    peak = float(np.abs(ref[mask]).max()) if mask.any() else float("nan")
    d = np.abs(delta[mask]) / peak if mask.any() else np.array([])
    return {"max": float(d.max()) if d.size else float("nan"), "rms": float(math.sqrt(np.mean(d * d))) if d.size else float("nan"),
            "peak_reference": peak}


def thermal_metrics(history, reference):
    """Spec section 9 metrics. Convective and radiated power: max/rms error relative to SESAM's peak over the hypersonic
    phase (SESAM's hot-wall factor clamps its heating to zero late in the flight (facts note s.15), so a point-wise
    relative error next to the clamp is unbounded and the high-Kn start carries negligible absolute heat); point-wise relative error
    of Q_conv restricted to the continuum part (Kn_ref < 0.01 and Q_ref above 10 % of its peak); integrated heat at the end of the hypersonic phase and
    at the end; energy-equivalent temperature vs SESAM's lumped temperature over the whole flight."""
    a = align_thermal(history, reference)
    hyp = a["V_ref"] > HYPERSONIC_MS
    peak = float(np.abs(a["Q_ref"][hyp]).max()) if hyp.any() else 0.0
    # point-wise only where SESAM's heat is at least 10 % of its peak: its hot-wall factor clamps the heating to zero
    # late in the flight, and relative errors next to the clamp are unbounded by construction
    cont = hyp & (a["kn_ref"] < CONTINUUM_KN) & (a["Q_ref"] > 0.1 * peak)
    dQ = a["Q_model"] - a["Q_ref"]
    rel_cont = np.abs(dQ[cont]) / a["Q_ref"][cont]
    i_hyp = int(np.nonzero(hyp)[0][-1]) if hyp.any() else len(a["t"]) - 1
    dT = a["T_model"] - a["T_ref"]
    return {
        "n_points": int(a["t"].size), "n_continuum_points": int(cont.sum()),
        "Q_conv": {**_peak_relative(dQ, a["Q_ref"], hyp), "continuum_rel_max": float(rel_cont.max()) if rel_cont.size else float("nan"),
                   "continuum_rel_rms": float(math.sqrt(np.mean(rel_cont * rel_cont))) if rel_cont.size else float("nan")},
        "integrated_heat": {"rel_error_end_of_hypersonic": float(a["H_model"][i_hyp] / a["H_ref"][i_hyp] - 1.0) if a["H_ref"][i_hyp] else float("nan"),
                            "rel_error_end": float(a["H_model"][-1] / a["H_ref"][-1] - 1.0) if a["H_ref"][-1] else float("nan"),
                            "model_end_J": float(a["H_model"][-1]), "reference_end_J": float(a["H_ref"][-1]),
                            "ratio_end": float(a["H_model"][-1] / a["H_ref"][-1]) if a["H_ref"][-1] else float("nan")},
        "temperature": {"dT_max_K": float(np.abs(dT).max()), "dT_rel_max": float((np.abs(dT) / a["T_ref"]).max()),
                        "dT_rms_K": float(math.sqrt(np.mean(dT * dT))),
                        "model_peak_K": float(a["T_model"].max()), "reference_peak_K": float(a["T_ref"].max())},
        "radiated": _peak_relative(a["R_model"] - a["R_ref"], a["R_ref"], hyp),
    }


def plot_thermal(history, reference, outdir, title):
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    a = align_thermal(history, reference)
    paths = [os.path.join(outdir, n) for n in THERMAL_PLOT_NAMES]
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    _overlay(ax, a["t"], a["Q_ref"] / 1e3, a["t"], a["Q_model"] / 1e3, "convective power [kW]", title=title)
    rx.plot(a["t"], (a["Q_model"] - a["Q_ref"]) / 1e3, color=MODEL_COLOR, lw=1.0)
    rx.axhline(0.0, color=MUTED, lw=0.6)
    rx.set_ylabel("model - SESAM [kW]"); rx.set_xlabel("time [s]")
    strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(paths[0], dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    _overlay(ax, a["t"], a["T_ref"], a["t"], a["T_model"], "temperature [K]", xlabel="time [s]", title=title)
    ax.plot(a["t"], a["Tmax_model"], color=MODEL_COLOR, lw=0.8, ls=":", label="model surface max")
    ax.plot(a["t"], a["Tmin_model"], color=MUTED, lw=0.8, ls=":", label="model surface min")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(paths[1], dpi=150); plt.close(fig)

    fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    _overlay(ax, a["t"], a["H_ref"] / 1e3, a["t"], a["H_model"] / 1e3, "integrated convective heat [kJ]", title=title)
    with np.errstate(divide="ignore", invalid="ignore"):
        rx.plot(a["t"], 100.0 * (a["H_model"] / a["H_ref"] - 1.0), color=MODEL_COLOR, lw=1.0)
    rx.axhline(0.0, color=MUTED, lw=0.6)
    rx.set_ylabel("model / SESAM - 1 [%]"); rx.set_xlabel("time [s]"); rx.set_ylim(-20.0, 20.0)
    strip_top_right_spines(rx)
    fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)
    return paths


# ---------------------------------------------------------------------------------------------------------------
# Step 3: mass loss, spraying and size distributions (spec Step 3 sections 12-13)

def has_melt(history, reference=None):
    """True when the model history carries the melt columns (and, if a reference is given, its mass varies)."""
    ok = "sprayed_mass_kg" in history.columns and "film_mass_kg" in history.columns
    if reference is None:
        return ok
    return ok and reference.mass is not None and float(reference.mass.min()) < 0.999 * float(reference.mass[0])


def _first_time(t, mask):
    return float(t[mask][0]) if mask.any() else None


def _crossing_time(t, m, level):
    """Time at which the decreasing series m first falls below `level`, interpolated linearly between samples."""
    below = np.flatnonzero(m < level)
    if below.size == 0:
        return None
    i = int(below[0])
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (m[i - 1] - level) / (m[i - 1] - m[i]) * (t[i] - t[i - 1]))


def melt_metrics(history, reference):
    """Mass vs SESAM at the reference's time stamps (relative to the initial mass), melt-onset altitude and the time
    at which the mass falls below DEMISE_FRACTION of the initial (SESAM removes the last gram at its 'demise'; the model
    stops at the fraction, so the 1 % crossing is the like-for-like end time; for the model it is taken on the body's
    material, film excluded, since a static leeward film can stay attached)."""
    c = history.columns
    t_mod, m_mod = c["time_s"], c["mass_kg"]
    mask = (reference.time >= t_mod[0]) & (reference.time <= t_mod[-1])
    t = reference.time[mask]
    m_ref = reference.mass[mask]
    m_at = np.interp(t, t_mod, m_mod)
    m0 = float(reference.mass[0])
    onset_ref = _first_time(reference.time, reference.mass < 0.999999 * m0)
    onset_mod = _first_time(t_mod, (c["removed_mass_kg"] + c["film_mass_kg"]) > 0.0)
    h_ref = float(np.interp(onset_ref, reference.time, reference.altitude)) / 1e3 if onset_ref is not None else None
    h_mod = float(np.interp(onset_mod, t_mod, c["altitude_km"])) if onset_mod is not None else None
    end_ref = _crossing_time(reference.time, reference.mass, DEMISE_FRACTION * m0)
    body_mass = m_mod - c["film_mass_kg"]                                                          # the body's material (a leeward film may stay attached)
    end_mod = _crossing_time(t_mod, body_mass, DEMISE_FRACTION * float(body_mass[0]))
    return {
        "n_points": int(t.size), "initial_mass_model_kg": float(m_mod[0]), "initial_mass_reference_kg": m0,
        "mass": {"max_abs_kg": float(np.abs(m_at - m_ref).max()), "max_rel_m0": float(np.abs(m_at - m_ref).max() / m0),
                 "rms_rel_m0": float(math.sqrt(np.mean((m_at - m_ref) ** 2)) / m0)},
        "onset_time_model_s": onset_mod, "onset_time_reference_s": onset_ref,
        "onset_altitude_model_km": h_mod, "onset_altitude_reference_km": h_ref,
        "onset_altitude_diff_km": (h_mod - h_ref) if (h_mod is not None and h_ref is not None) else None,
        "demise_time_model_s": end_mod, "demise_time_reference_s": end_ref,
        "demise_time_diff_s": (end_mod - end_ref) if (end_mod is not None and end_ref is not None) else None,
        "demise_time_rel": ((end_mod - end_ref) / end_ref) if (end_mod is not None and end_ref is not None) else None,
        "final_mass_model_kg": float(m_mod[-1]), "sprayed_mass_kg": float(c["sprayed_mass_kg"][-1]), "film_mass_kg": float(c["film_mass_kg"][-1]),
    }


def plot_melt(history, outdir, title, reference=None, size_distribution_csv=None):
    """The Step 3 plots (MELT_PLOT_NAMES); with a melting reference the mass plots carry the SESAM overlay and a
    residual panel. Returns the paths."""
    import csv
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    c = history.columns
    t, h, m = c["time_s"], c["altitude_km"], c["mass_kg"]
    paths = [os.path.join(outdir, n) for n in MELT_PLOT_NAMES]
    # 1. mass vs time (+ residual)
    if reference is not None:
        mask = (reference.time >= t[0]) & (reference.time <= t[-1])
        tr, mr = reference.time[mask], reference.mass[mask]
        fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
        _overlay(ax, tr, mr, t, m, "mass [kg]", title=title)
        rx.plot(tr, 100.0 * (np.interp(tr, t, m) - mr) / reference.mass[0], color=MODEL_COLOR, lw=1.0)
        rx.axhline(0.0, color=MUTED, lw=0.6)
        rx.set_ylabel("model - SESAM [% of m0]"); rx.set_xlabel("time [s]")
        strip_top_right_spines(rx)
    else:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(t, m, color=MODEL_COLOR, lw=1.2, label="model")
        ax.set_ylabel("mass [kg]"); ax.set_xlabel("time [s]"); ax.set_title(title, color=SECOND, fontsize=10)
        strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[0], dpi=150); plt.close(fig)
    # 2. mass vs altitude
    fig, ax = plt.subplots(figsize=(7, 5))
    if reference is not None:
        ax.plot(reference.mass, reference.altitude / 1e3, color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(m, h, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_xlabel("mass [kg]"); ax.set_ylabel("altitude [km]"); ax.set_title(title, color=SECOND, fontsize=10); ax.legend(frameon=False)
    strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[1], dpi=150); plt.close(fig)
    # 3. mass budget
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, c["sprayed_mass_kg"], color=MODEL_COLOR, lw=1.2, label="sprayed (cumulative)")
    ax.plot(t, c["removed_mass_kg"] - c["sprayed_mass_kg"], color=SECOND, lw=1.0, ls=":", label="removed instantly (cumulative)")
    ax.plot(t, c["film_mass_kg"], color=MUTED, lw=1.2, label="film")
    ax.plot(t, m, color=REF_COLOR, lw=1.2, ls="--", label="remaining (body + film)")
    ax.set_xlabel("time [s]"); ax.set_ylabel("mass [kg]"); ax.set_title(title, color=SECOND, fontsize=10); ax.legend(frameon=False)
    strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)
    # 4. theta_cr and spraying area
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax.plot(t, c["theta_cr_deg"], color=MODEL_COLOR, lw=1.0); ax.set_ylabel("theta_cr [deg]"); ax.set_title(title, color=SECOND, fontsize=10)
    bx.plot(t, c["spraying_area_m2"] * 1e4, color=MODEL_COLOR, lw=1.0); bx.set_ylabel("spraying area [cm2]"); bx.set_xlabel("time [s]")
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[3], dpi=150); plt.close(fig)
    # 5. closure fractions and the two Knudsen numbers
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.stackplot(t, c["closure_fraction_girin"], c["closure_fraction_couette_slip"], c["closure_fraction_couette_fm"],
                 labels=["Girin closure (Kn_local < 0.01)", "Couette, slip-corrected shear", "Couette, free-molecular shear"],
                 colors=[INK, MODEL_COLOR, MUTED], alpha=0.8)
    ax.set_xlabel("time [s]"); ax.set_ylabel("windward area fraction"); ax.set_ylim(0.0, 1.0); ax.legend(frameon=False, loc="upper left", fontsize=8)
    bx2 = ax.twinx()
    bx2.semilogy(t, c["kn_body"], color=SECOND, lw=1.0, ls="--")
    bx2.semilogy(t, c["kn_local_stag"], color=SECOND, lw=1.0, ls=":")
    bx2.axhline(0.01, color=MUTED, lw=0.6)
    bx2.set_ylabel("Kn_body (dashed), Kn_local at the nose (dotted)", color=SECOND, fontsize=8)
    ax.set_title(title, color=SECOND, fontsize=10); strip_top_right_spines(ax)
    fig.tight_layout(); fig.savefig(paths[4], dpi=150); plt.close(fig)
    # 6. droplet size and film thickness
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax.semilogy(t, c["r_median_um"], color=MODEL_COLOR, lw=1.0, label="median r")
    ax.semilogy(t, c["r_max_um"], color=MUTED, lw=0.8, ls=":", label="max r")
    ax.set_ylabel("droplet radius [um]"); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    bx.semilogy(t, np.maximum(c["film_thickness_mean_mm"], 1e-6), color=MODEL_COLOR, lw=1.0, label="mean film thickness")
    bx.semilogy(t, np.maximum(c["film_thickness_max_mm"], 1e-6), color=MUTED, lw=0.8, ls=":", label="max")
    bx.set_ylabel("film thickness [mm]"); bx.set_xlabel("time [s]"); bx.legend(frameon=False)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[5], dpi=150); plt.close(fig)
    # 7. size distributions (from size_distribution.csv when present next to the run, else skipped with an empty axes)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.5))
    path = size_distribution_csv or os.path.join(outdir, "size_distribution.csv")
    if os.path.isfile(path):
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        windows = sorted({(float(r["window_start_s"]), float(r["window_end_s"])) for r in rows})
        for lo, hi in windows:
            sel = [r for r in rows if float(r["window_start_s"]) == lo and float(r["window_end_s"]) == hi]
            r_mid = np.sqrt(np.array([float(r["r_lo_m"]) for r in sel]) * np.array([float(r["r_hi_m"]) for r in sel])) * 1e6
            dn, dm = np.array([float(r["dn"]) for r in sel]), np.array([float(r["dm_kg"]) for r in sel])
            total = (lo, hi) == windows[-1] and lo == 0.0
            kw = dict(color=INK, lw=1.6, label="flight") if total else dict(color=MODEL_COLOR, lw=0.7, alpha=0.5)
            if dn.sum() > 0:
                ax.loglog(r_mid, np.maximum(dn, 1e-300), **kw)
                bx.loglog(r_mid, np.maximum(dm, 1e-300), **kw)
        ax.set_ylim(bottom=max(ax.get_ylim()[0], 1e-1)); bx.set_ylim(bottom=max(bx.get_ylim()[0], 1e-12))
    ax.set_xlabel("droplet radius [um]"); ax.set_ylabel("dn per bin"); bx.set_xlabel("droplet radius [um]"); bx.set_ylabel("dM per bin [kg]")
    ax.set_title(title, color=SECOND, fontsize=10)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False)
    strip_top_right_spines(ax); strip_top_right_spines(bx)
    fig.tight_layout(); fig.savefig(paths[6], dpi=150); plt.close(fig)
    return paths
