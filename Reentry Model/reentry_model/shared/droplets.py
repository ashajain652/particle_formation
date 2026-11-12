"""The droplet outputs of a run: particles.npz, particles_summary.csv and size_distribution.csv."""
import os

import numpy as np

# particles_summary.csv, one row a macro step: the element model's history columns
SUMMARY_COLUMNS = ["time_s", "altitude_km", "velocity_kms", "released_mass_kg", "n_released", "r_median_um", "r_max_um",
                   "theta_cr_deg", "spraying_area_m2", "film_mass_kg", "film_thickness_mean_mm"]
# a skin flight's (sub-plan 22 Task 4): (the summary's column, the skin history's column it is read from) -- the same
# droplet columns, the step's released mass from released_step_mass_kg, since the skin history's released_mass_kg is
# cumulative, and where the element model has its film, which the skins do not, the skins' liquid and flowing depths
SKIN_SUMMARY_COLUMNS = [("time_s", "time_s"), ("altitude_km", "altitude_km"), ("velocity_kms", "velocity_kms"),
                        ("released_mass_kg", "released_step_mass_kg"), ("n_released", "n_released"),
                        ("r_median_um", "r_median_um"), ("r_max_um", "r_max_um"), ("theta_cr_deg", "theta_cr_deg"),
                        ("spraying_area_m2", "spraying_area_m2"), ("liquid_depth_mean_mm", "liquid_depth_mean_mm"),
                        ("flowing_depth_mean_mm", "flowing_depth_mean_mm")]


def write_particles(run_dir, body, history, window=10.0):
    """particles.npz (the source table, SOURCE_COLUMNS as arrays), particles_summary.csv (per macro step) and
    size_distribution.csv (dn, dM per log bin per `window`-second window and in total). Returns the paths.

    The table is the body's `source_rows`, the element model's rows, the layered skins' blocks of them (sub-plan 22 Task
    2) or their `spray.SourceSpool` on disk (Task 4), read column by column (`spray.source_columns`), so that a flight's
    records need not be in memory at once: the numbers written are the same either way. The summary's columns are the
    element model's history columns (SUMMARY_COLUMNS); a skin flight's are SKIN_SUMMARY_COLUMNS, read from its history's
    droplet columns (Task 4); a history with neither writes no summary, and the paths returned then have none."""
    import csv
    from . import spray
    os.makedirs(run_dir, exist_ok=True)
    cols = spray.source_columns(body.source_rows)
    npz = os.path.join(run_dir, "particles.npz")
    np.savez_compressed(npz, **cols)
    c = history.columns
    summary = os.path.join(run_dir, "particles_summary.csv")
    if set(SUMMARY_COLUMNS) <= set(c):
        header, keys = SUMMARY_COLUMNS, SUMMARY_COLUMNS
    elif {k for _, k in SKIN_SUMMARY_COLUMNS} <= set(c):
        header, keys = [h for h, _ in SKIN_SUMMARY_COLUMNS], [k for _, k in SKIN_SUMMARY_COLUMNS]
    else:
        summary = None
    if summary is not None:
        with open(summary, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            for i in range(len(history)):
                w.writerow(["{:.9g}".format(c[k][i]) for k in keys])
    dist = os.path.join(run_dir, "size_distribution.csv")
    edges = spray.BIN_EDGES
    t, r, dn, dm = cols["time_s"], cols["r_m"], cols["dn"], cols["dm_kg"]
    t_end = float(c["time_s"][-1])
    starts = np.arange(0.0, t_end + 1e-9, window)
    with open(dist, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["window_start_s", "window_end_s", "r_lo_m", "r_hi_m", "dn", "dm_kg"])
        for t0 in list(starts) + [None]:
            if t0 is None:
                sel = np.ones(len(t), dtype=bool)
                lo, hi = 0.0, t_end
            else:
                lo, hi = t0, t0 + window
                sel = (t >= lo) & (t < hi)
            n_hist, m_hist = spray.histogram(r[sel], dn[sel], dm[sel]) if sel.any() else (np.zeros(spray.N_BINS), np.zeros(spray.N_BINS))
            for j in range(spray.N_BINS):
                w.writerow(["{:.9g}".format(lo), "{:.9g}".format(hi), "{:.9g}".format(edges[j]), "{:.9g}".format(edges[j + 1]),
                            "{:.9g}".format(n_hist[j]), "{:.9g}".format(m_hist[j])])
    out = {"particles": npz, "particles_summary": summary, "size_distribution": dist}
    return {k: v for k, v in out.items() if v is not None}
