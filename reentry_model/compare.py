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

MODEL_COLOR, REF_COLOR = "#3987e5", "#0b0b0b"
HYPERSONIC_MS = 1000.0
PLOT_NAMES = ("velocity_time.png", "altitude_time.png", "altitude_velocity.png", "angles_time.png",
              "ground_track.png", "regime_drag.png")


def align(history, reference):
    t_mod = history.columns["time_s"]
    mask = reference.time <= t_mod[-1]
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


def metrics(history, reference):
    a = align(history, reference)
    t_model_end, t_ref_end = float(history.columns["time_s"][-1]), float(reference.time[-1])
    return {
        "n_points": int(a["t"].size),
        "model_end_time_s": t_model_end, "reference_end_time_s": t_ref_end,
        "d_end_time_s": t_model_end - t_ref_end, "d_end_time_rel": (t_model_end - t_ref_end) / t_ref_end,
        "final_velocity_model_ms": float(history.columns["velocity_kms"][-1] * 1e3),
        "final_velocity_reference_ms": float(reference.velocity[-1]),
        "all": _phase(a, np.ones(a["t"].size, dtype=bool)),
        "hypersonic": _phase(a, a["V_ref"] > HYPERSONIC_MS),
    }


def _overlay_with_residual(a, key, ylabel, resid_label, scale, path, title):
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(a["t"], a[key + "_ref"] * scale, color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(a["t"], a[key + "_model"] * scale, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_ylabel(ylabel); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    rx.plot(a["t"], (a[key + "_model"] - a[key + "_ref"]), color=MODEL_COLOR, lw=1.0)
    rx.axhline(0.0, color=MUTED, lw=0.6)
    rx.set_ylabel(resid_label); rx.set_xlabel("time [s]")
    for x in (ax, rx):
        strip_top_right_spines(x)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def plot_all(history, reference, outdir, title):
    os.makedirs(outdir, exist_ok=True)
    apply_rcparams(plt)
    a = align(history, reference)
    paths = [os.path.join(outdir, n) for n in PLOT_NAMES]
    _overlay_with_residual(a, "V", "velocity [km/s]", "model - SESAM [m/s]", 1e-3, paths[0], title)
    _overlay_with_residual(a, "h", "altitude [km]", "model - SESAM [m]", 1e-3, paths[1], title)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(a["V_ref"] / 1e3, a["h_ref"] / 1e3, color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(a["V_model"] / 1e3, a["h_model"] / 1e3, color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_xlabel("velocity [km/s]"); ax.set_ylabel("altitude [km]"); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    strip_top_right_spines(ax); fig.tight_layout(); fig.savefig(paths[2], dpi=150); plt.close(fig)

    fig, (g, hd) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for ax, key, lab in ((g, "gamma", "flight-path angle [deg]"), (hd, "heading", "heading [deg]")):
        ax.plot(a["t"], a[key + "_ref"], color=REF_COLOR, lw=1.6, label="SESAM")
        ax.plot(a["t"], a[key + "_model"], color=MODEL_COLOR, lw=1.2, ls="--", label="model")
        ax.set_ylabel(lab); strip_top_right_spines(ax)
    g.legend(frameon=False); g.set_title(title, color=SECOND, fontsize=10); hd.set_xlabel("time [s]")
    fig.tight_layout(); fig.savefig(paths[3], dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(a["lon_ref"], a["lat_ref"], color=REF_COLOR, lw=1.6, label="SESAM")
    ax.plot(a["lon_model"], a["lat_model"], color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    ax.set_xlabel("longitude [deg]"); ax.set_ylabel("latitude [deg]"); ax.legend(frameon=False); ax.set_title(title, color=SECOND, fontsize=10)
    strip_top_right_spines(ax); fig.tight_layout(); fig.savefig(paths[4], dpi=150); plt.close(fig)

    fig, (kx, cx) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    kx.semilogy(a["t"], np.maximum(a["kn_ref"], 1e-7), color=REF_COLOR, lw=1.6, label="SESAM")
    kx.semilogy(a["t"], np.maximum(a["kn_model"], 1e-7), color=MODEL_COLOR, lw=1.2, ls="--", label="model")
    for thr in (10.0, 1.0, 0.1, 0.01):
        kx.axhline(thr, color=MUTED, lw=0.5, ls=":")
    kx.set_ylabel("Knudsen number"); kx.legend(frameon=False); kx.set_title(title, color=SECOND, fontsize=10)
    cx.plot(a["t"], a["cd_ref"], color=REF_COLOR, lw=1.6); cx.plot(a["t"], a["cd_model"], color=MODEL_COLOR, lw=1.2, ls="--")
    cx.set_ylabel("C_D"); cx.set_xlabel("time [s]")
    for x in (kx, cx):
        strip_top_right_spines(x)
    fig.tight_layout(); fig.savefig(paths[5], dpi=150); plt.close(fig)
    return paths
