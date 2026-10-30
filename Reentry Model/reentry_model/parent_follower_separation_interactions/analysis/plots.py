"""The plots (spec 2 goal 6 and spec 5, Task 10), read back from the files of the sweep and the tables:

    "$PY" -m reentry_model.parent_follower_separation_interactions.analysis.plots [--out DIR] [--theta DEG] [--bound low|high]
        [--trajectory-wake-flow MODE]

`--out` is the folder `sweep` and `tables` wrote to (default `parent_follower_separation_interactions.OUTPUT`); the
pictures go to `<out>/plots/`. Nothing is integrated here: the plots read summary.csv, runs/<name>.csv|json,
plume_size.csv, ballistic_coefficient.csv and affix_regimes.csv. Prints how many runs summary.csv holds and which
plots it wrote or skipped.

  plume_size.png             r_b (low bound solid, high dashed) against distance behind the parent s and against the time
                             since it passed (s / U), one line per release altitude with the parent's plume regime in the
                             legend; Taylor's fit limit 1200 sqrt(C_D A) and Park's 30 km limit marked on each line.
  margin_<h>km.png           for each release altitude: the distance from the parent's track (thick) and the plume radius
                             r_b (thin) against time since release, for every integrated case at that altitude (far-wake,
                             escape and surface releases), one panel per wake_flow, for one bound (`--bound`). Colour is
                             the diameter, line style the release kind, and the legend lists the thetas of each group.
                             A far-wake integration is shared by all its thetas, so it is drawn once.
  trajectories_<h>km.png     altitude against downrange of the parent (the twin) and of each integrated case (left), and
                             altitude against the follower's lag behind the parent (right), for one wake_flow
                             (`--trajectory-wake-flow`, default off) and one bound.
  residence_map.png          over diameter and altitude at one release angle (`--theta`, or the nearest in the sweep),
                             one panel per wake_flow and bound. A far-wake release (it starts inside the plume) shows the
                             residence time; an escape or surface release (it starts outside) shows the total time spent
                             in the plume. The release kind is hatched and each cell's special outcome is a letter.
  affix_zones.png            the AFFIx zone over diameter and release angle with theta_bif drawn; one panel per altitude
                             when the validity flags differ between altitudes, otherwise one panel.
  shock_layer_paths_<h>km.png  for each altitude with rarefied followers: their paths in the parent's frame (x forward, y
                             from the axis), the body, the base plane and the shock curve; left near the body, right
                             over the first second (or to the shock-layer exit if later).
"""
import argparse
import csv
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Rectangle

from reentry_model import DRAMA_DIR, earth
from reentry_model.parent_follower_separation_interactions import OUTPUT
from reentry_model.parent_follower_separation_interactions.adapters import PARENT_FRONTAL_AREA
from reentry_model.parent_follower_separation_interactions.wake_gas import WAKE_FLOW_MODES

sys.path.insert(0, os.path.join(DRAMA_DIR, "analysis"))           # plot_style lives with the sweep plots
try:
    from plot_style import INK, MUTED, SECOND, SEQ_BLUE, apply_rcparams, strip_top_right_spines
except ImportError:                                    # DRAMA Sphere Experiments/analysis/ not present: plain matplotlib
    INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#898781"
    SEQ_BLUE = plt.get_cmap("Blues")

    def apply_rcparams(plt, font_size=10.5):
        return None

    def strip_top_right_spines(ax):
        return None

INTEGRATED = ("far_wake", "escape", "surface")            # release kinds whose follower was integrated
RESIDENCE_KINDS = ("far_wake",)                           # start inside the plume: shown by residence_time (the others by total_inside)
KIND_STYLE = {"far_wake": "-", "escape": "-.", "surface": ":"}
KIND_HATCH = {"far_wake": "", "escape": "///", "surface": "xxx"}
ZONE_COLOR = {"far_wake": "#256abf", "shock_band": "#e69f00", "escape": "#2a9d6f"}
EXIT_MARKER = {"through_shock": "o", "past_base_plane": "s", "hit_parent": "X", "inside_at_end": "^"}
# letters printed in the residence map's cells (reason or censored flag), explained in its legend
CELL_FLAGS = {"censored": "c", "reached_parent": "r", "ahead_of_parent": "a", "hit_parent": "h",
              "rarefied_no_law": "n", "too_slow_no_law": "n", "shock_band": "b", "near_wake_trapped": "t"}
CELL_FLAG_TEXT = {"c": "censored (still inside when the integration ended)", "r": "reached_parent (fell behind, caught up)",
                  "a": "ahead_of_parent (never fell behind)", "h": "hit_parent (struck the parent)",
                  "n": "no plume law (rarefied or too slow parent)", "b": "shock_band (not integrated)",
                  "t": "near_wake_trapped (not integrated)"}
# one- or two-letter codes printed in the AFFIx zone cells for the validity flags that vary
AFFIX_FLAG_CODES = {"flag_not_hypersonic": "H", "flag_r_below_fit": "F", "flag_r_above_rstar": "R",
                    "flag_rarefied_follower": "f", "flag_rarefied_parent": "P", "flag_unequal_density": "U"}
VALUE_FLOOR_S = 1e-3          # s, the lowest value of the residence map's log colour scale (a zero is drawn at this floor)
SHOCK_VIEW_RADII = (-3.0, 2.2, 3.0)       # shock-layer panel: x from -3 R to 2.2 R, y from 0 to 3 R
PATH_WINDOW_S = 1.0           # s, a shock-layer path is drawn over this time after release, or to its exit if later
MARKERS = ("o", "s", "^", "D", "v", "P")      # one per release angle, to tell the lines of one diameter and kind apart
MAX_POINTS = 300              # points of a history drawn in the trajectory plot


# ------------------------------------------------------------------------------------------------------ reading

def read_rows(path):
    """Rows (dicts of text) of a CSV file written by `sweep.write_csv`, which may start with '#' comment lines."""
    with open(path) as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def num(text):
    """A float from CSV text; NaN for a blank."""
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def is_true(text):
    return str(text).strip() == "True"


def load_summary(out):
    """summary.csv as rows with the numeric fields as floats; None when the file is missing."""
    path = os.path.join(out, "summary.csv")
    if not os.path.exists(path):
        return None
    rows = read_rows(path)
    for r in rows:
        for k in ("diameter_mm", "altitude_km", "theta_deg", "residence_time", "first_interval", "t_first_entry",
                  "total_inside", "t_shock_layer_exit", "x_exit_parent_frame", "y_exit_parent_frame"):
            if k in r:
                r[k] = num(r[k])
        r["censored"] = is_true(r.get("censored"))
    return rows


def load_history(out, name):
    """{column: array} of runs/<name>.csv (the numeric columns), or None (with a printed note) when the file is missing
    or empty."""
    path = os.path.join(out, "runs", name + ".csv")
    if not os.path.exists(path):
        print(f"note: runs/{name}.csv not found; that case is left out of the plot")
        return None
    rows = read_rows(path)
    if not rows:
        print(f"note: runs/{name}.csv has no rows; that case is left out of the plot")
        return None
    return {k: np.array([num(r[k]) for r in rows]) for k in rows[0] if k != "regime"}


def load_case_json(out, name):
    with open(os.path.join(out, "runs", name + ".json")) as fh:
        return json.load(fh)


def integrated(rows, h=None, bound=None):
    """The rows of integrated cases (release kind far-wake, escape or surface), at altitude h and bound if given."""
    return [r for r in rows if r.get("release_kind") in INTEGRATED and (h is None or r["altitude_km"] == h)
            and (bound is None or r["bound"] == bound)]


def pick(values, preferred):
    """`preferred` if it is among `values`, else the first of them (None when there are none)."""
    values = list(values)
    return preferred if preferred in values else (values[0] if values else None)


def diameter_colors(diameters):
    """{diameter: colour}, a perceptually ordered ramp from the smallest to the largest."""
    ds = sorted(set(diameters))
    cols = plt.cm.viridis(np.linspace(0.05, 0.88, max(len(ds), 2)))
    return {d: cols[i] for i, d in enumerate(ds)}


def groups_by_kind(rows):
    """{(diameter, release kind): [rows]}; a far-wake integration is shared by every theta of its (diameter, bound,
    wake_flow), so only the first row of each far-wake (diameter, bound, wake_flow) is kept."""
    groups, seen = {}, set()
    for r in sorted(rows, key=lambda r: (r["diameter_mm"], r["theta_deg"])):
        if r["release_kind"] == "far_wake":
            key = (r["diameter_mm"], r["bound"], r["wake_flow"])
            if key in seen:
                continue
            seen.add(key)
        groups.setdefault((r["diameter_mm"], r["release_kind"]), []).append(r)
    return groups


def thetas_text(rows):
    return "/".join(f"{t:g}" for t in sorted({r["theta_deg"] for r in rows}))


def save(fig, out, name):
    folder = os.path.join(out, "plots")
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, name)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")
    return path


def new_axes(*args, **kwargs):
    apply_rcparams(plt)
    fig, axes = plt.subplots(*args, **kwargs)
    for ax in np.atleast_1d(axes).ravel():
        strip_top_right_spines(ax)
    return fig, axes


# --------------------------------------------------------------------------------------------------- plume_size

def _loglog_at(x, xs, ys):
    """y at x by linear interpolation of log y against log x, NaN outside the data."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    ok = (xs > 0) & (ys > 0)
    if ok.sum() < 2 or not (xs[ok].min() <= x <= xs[ok].max()):
        return math.nan
    return float(math.exp(np.interp(math.log(x), np.log(xs[ok]), np.log(ys[ok]))))


def plot_plume_size(out):
    path = os.path.join(out, "plume_size.csv")
    if not os.path.exists(path):
        print("plume_size.png skipped: plume_size.csv not found (run tables first)")
        return None
    rows = read_rows(path)
    cda = {}
    if os.path.exists(os.path.join(out, "ballistic_coefficient.csv")):
        cda = {num(r["altitude_km"]): num(r["parent_cda_m2"]) for r in read_rows(os.path.join(out, "ballistic_coefficient.csv"))}
    altitudes = sorted({num(r["altitude_km"]) for r in rows}, reverse=True)
    colors = dict(zip(altitudes, plt.cm.viridis(np.linspace(0.05, 0.88, max(len(altitudes), 2)))))
    fig, (ax_s, ax_t) = new_axes(1, 2, figsize=(12.5, 5.0))
    for h in altitudes:
        rs = [r for r in rows if num(r["altitude_km"]) == h]
        regime = rs[0]["regime"]
        s = np.array([num(r["s_m"]) for r in rs])
        t = np.array([num(r["t_since_passage_s"]) for r in rs])
        lo = np.array([num(r["r_b_low_m"]) for r in rs])
        hi = np.array([num(r["r_b_high_m"]) for r in rs])
        label = f"{h:g} km, {regime}" + (" (no plume law)" if not np.isfinite(lo).any() else "")
        ax_s.plot([], [], color=colors[h], label=label)
        for ax, x in ((ax_s, s), (ax_t, t)):
            ax.plot(x, lo, "-o", ms=3, color=colors[h])
            ax.plot(x, hi, "--", color=colors[h])
        if not np.isfinite(lo).any():
            continue
        speed = num(rs[0]["parent_speed_ms"])
        if h in cda and math.isfinite(cda[h]):
            s_tay = 1200.0 * math.sqrt(cda[h])               # Taylor's fit limit x <= 1200 sqrt(C_D A)
            r_tay = _loglog_at(s_tay, s, lo)
            if math.isfinite(r_tay):
                ax_s.plot([s_tay], [r_tay], "v", color=colors[h], ms=8, mec=INK)
                ax_t.plot([s_tay / speed], [r_tay], "v", color=colors[h], ms=8, mec=INK)
        r_park = _loglog_at(30e3, s, lo)
        if math.isfinite(r_park):
            ax_t.plot([30e3 / speed], [r_park], "D", color=colors[h], ms=6, mec=INK)
    ax_s.axvline(30e3, color=MUTED, lw=1, ls=":")
    for ax, xl in ((ax_s, "distance behind the parent, s [m]"), (ax_t, "time since the parent passed, s / U [s]")):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(xl); ax.set_ylabel("plume radius r_b [m]")
        ax.grid(alpha=0.3, which="both", lw=0.4)
    handles, labels = ax_s.get_legend_handles_labels()
    handles += [Line2D([], [], color=INK, ls="-"), Line2D([], [], color=INK, ls="--"),
                Line2D([], [], color=INK, marker="v", ls="", mec=INK),
                Line2D([], [], color=INK, marker="D", ls="", mec=INK)]
    labels += ["low bound", "high bound", "Taylor fit limit, 1200 (C_D A)^(1/2)", "Park limit, 30 km (dotted line, left)"]
    ax_s.legend(handles, labels, fontsize=8, loc="upper left", frameon=False)
    fig.suptitle("Plume radius of the parent against distance and against time since passage", color=INK)
    return save(fig, out, "plume_size.png")


# ------------------------------------------------------------------------------------------- margin, trajectories

def plot_margin(out, rows, h, bound):
    """margin_<h>km.png; None (with a note) when the altitude has no integrated case."""
    cases = integrated(rows, h, bound)
    if not cases:
        print(f"margin_{h:g}km.png skipped: no integrated case at {h:g} km (bound {bound})")
        return None
    flows = [f for f in WAKE_FLOW_MODES if any(r["wake_flow"] == f for r in cases)]
    colors = diameter_colors(r["diameter_mm"] for r in cases)
    thetas = sorted({r["theta_deg"] for r in cases})
    fig, axes = new_axes(1, len(flows), figsize=(6.2 * len(flows) + 2.6, 5.0), squeeze=False)
    used = {}
    for ax, flow in zip(axes[0], flows):
        groups = groups_by_kind([r for r in cases if r["wake_flow"] == flow])
        for (d, kind), rs in groups.items():
            used[(d, kind)] = rs
            for r in rs:
                hist = load_history(out, r["name"])
                if hist is None:
                    continue
                t = hist["t_since_release"]
                ok = t > 0
                mark = MARKERS[thetas.index(r["theta_deg"]) % len(MARKERS)]      # one marker per release angle
                ax.plot(t[ok], hist["distance_m"][ok], KIND_STYLE[kind], color=colors[d], lw=1.6, marker=mark,
                        ms=4, markevery=max(int(ok.sum()) // 12, 1))
                ax.plot(t[ok], np.where(hist["r_b"][ok] > 0, hist["r_b"][ok], np.nan), KIND_STYLE[kind],
                        color=colors[d], lw=0.7, alpha=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("time since release [s]"); ax.set_title(f"wake_flow {flow}", color=SECOND)
        ax.grid(alpha=0.3, which="both", lw=0.4)
    axes[0][0].set_ylabel("distance to the track (thick) and plume radius r_b (thin) [m]")
    handles = [Line2D([], [], color=colors[d], lw=2) for d in sorted(colors)]
    labels = [f"{d:g} mm" for d in sorted(colors)]
    for kind in INTEGRATED:
        ks = [(d, rs) for (d, k), rs in used.items() if k == kind]
        if ks:
            handles.append(Line2D([], [], color=INK, ls=KIND_STYLE[kind]))
            labels.append(f"{kind} (theta {thetas_text([r for _, rs in ks for r in rs])})")
    for i, th in enumerate(thetas):
        handles.append(Line2D([], [], color=INK, marker=MARKERS[i % len(MARKERS)], ms=4, ls=""))
        labels.append(f"marker: theta {th:g} degrees")
    axes[0][-1].legend(handles, labels, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.suptitle(f"Distance to the parent's track against plume radius, {h:g} km release, {bound} bound "
                 f"(a far-wake run is shared by all its release angles)", color=INK)
    return save(fig, out, f"margin_{h:g}km.png")


def altitude_downrange(hist):
    """(altitude [km], downrange [km] from the parent's first position, lag behind the parent [km]) of a case's
    follower, and the same two for the parent: ((alt, down), lag, (alt_p, down_p)). Sampled at most MAX_POINTS times."""
    n = len(hist["t"])
    idx = np.unique(np.linspace(0, n - 1, min(n, MAX_POINTS)).astype(int))
    lat0 = lon0 = None
    cols = {}
    for who in ("follower", "parent"):
        alt, down = [], []
        for i in idx:
            h, lat, lon = earth.ecef_to_geodetic([hist[f"{who}_x"][i], hist[f"{who}_y"][i], hist[f"{who}_z"][i]])
            if lat0 is None:
                lat0, lon0 = lat, lon
            alt.append(h / 1e3)
            down.append(earth.great_circle_distance(lat0, lon0, lat, lon) / 1e3)
        cols[who] = (np.array(alt), np.array(down))
    return cols["follower"], hist["lag_m"][idx] / 1e3, cols["parent"]


def plot_trajectories(out, rows, h, bound, wake_flow):
    """trajectories_<h>km.png; None (with a note) when the altitude has no integrated case."""
    cases = integrated(rows, h, bound)
    flow = pick([f for f in WAKE_FLOW_MODES if any(r["wake_flow"] == f for r in cases)], wake_flow)
    cases = [r for r in cases if r["wake_flow"] == flow]
    if not cases:
        print(f"trajectories_{h:g}km.png skipped: no integrated case at {h:g} km (bound {bound})")
        return None
    colors = diameter_colors(r["diameter_mm"] for r in cases)
    fig, (ax, ax_lag) = new_axes(1, 2, figsize=(13.5, 5.2))
    used, parent = {}, None                       # the parent's own line: the longest history among the cases
    for (d, kind), rs in groups_by_kind(cases).items():
        used[(d, kind)] = rs
        for r in rs:
            hist = load_history(out, r["name"])
            if hist is None:
                continue
            (alt, down), lag, (alt_p, down_p) = altitude_downrange(hist)
            if parent is None or hist["t"][-1] > parent[0]:
                parent = (hist["t"][-1], down_p, alt_p)
            ax.plot(down, alt, KIND_STYLE[kind], color=colors[d], lw=1.3)
            ax_lag.plot(lag, alt, KIND_STYLE[kind], color=colors[d], lw=1.3)
    if parent is not None:
        ax.plot(parent[1], parent[2], color=MUTED, lw=4, alpha=0.5, zorder=1)
    ax.set_xlabel("downrange from the release point [km]"); ax.set_ylabel("altitude [km]")
    ax_lag.set_xlabel("lag behind the parent along its track [km] (negative: ahead of it)")
    ax_lag.set_ylabel("altitude [km]"); ax_lag.axvline(0.0, color=MUTED, lw=0.8)
    for a in (ax, ax_lag):
        a.grid(alpha=0.3, lw=0.4)
    handles = [Line2D([], [], color=MUTED, lw=4, alpha=0.5)] + [Line2D([], [], color=colors[d], lw=2) for d in sorted(colors)]
    labels = ["parent (the twin)"] + [f"{d:g} mm" for d in sorted(colors)]
    for kind in INTEGRATED:
        ks = [rs for (d, k), rs in used.items() if k == kind]
        if ks:
            handles.append(Line2D([], [], color=INK, ls=KIND_STYLE[kind]))
            labels.append(f"{kind} (theta {thetas_text([r for rs in ks for r in rs])})")
    ax_lag.legend(handles, labels, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.suptitle(f"Trajectories after release at {h:g} km, wake_flow {flow}, {bound} bound", color=INK)
    return save(fig, out, f"trajectories_{h:g}km.png")


# --------------------------------------------------------------------------------------------- residence map

def map_cell(r):
    """The residence map's cell of a summary row: {"value": seconds or NaN (not integrated), "kind", "flags": [letters]}.
    A far-wake release starts inside the plume and shows its residence time. An escape or surface release starts outside
    it, so its residence time is 0 and its first interval can begin long after release; it shows the total time spent in
    the plume instead (0 when it never entered)."""
    kind = r.get("release_kind") or ""
    flags = []
    if r.get("censored"):
        flags.append(CELL_FLAGS["censored"])
    if r.get("reason") in CELL_FLAGS:
        flags.append(CELL_FLAGS[r["reason"]])
    if kind in RESIDENCE_KINDS:
        value = r["residence_time"]
    elif kind in INTEGRATED:
        value = 0.0 if math.isnan(r["total_inside"]) else r["total_inside"]
    else:
        value = math.nan
    return {"value": value, "kind": kind, "flags": sorted(set(flags))}


def nearest_theta(rows, theta):
    """The release angle present in the rows that is nearest `theta`."""
    return min({r["theta_deg"] for r in rows}, key=lambda t: abs(t - theta))


def plot_residence_map(out, rows, theta):
    if not rows:
        print("residence_map.png skipped: summary.csv has no rows")
        return None
    use = nearest_theta(rows, theta)
    if use != theta:
        print(f"residence_map: no case at theta {theta:g} degrees; using the nearest, {use:g}")
    rows = [r for r in rows if r["theta_deg"] == use]
    diameters = sorted({r["diameter_mm"] for r in rows})
    altitudes = sorted({r["altitude_km"] for r in rows})
    flows = [f for f in WAKE_FLOW_MODES if any(r["wake_flow"] == f for r in rows)]
    bounds = [b for b in ("low", "high") if any(r["bound"] == b for r in rows)]
    cells = {(r["wake_flow"], r["bound"], r["altitude_km"], r["diameter_mm"]): map_cell(r) for r in rows}
    values = [c["value"] for c in cells.values() if math.isfinite(c["value"])]
    norm = LogNorm(vmin=VALUE_FLOOR_S, vmax=max(max(values, default=1.0), 10.0 * VALUE_FLOOR_S))
    fig, axes = new_axes(len(flows), len(bounds), figsize=(max(7.5, 3.0 + 0.8 * len(diameters) * len(bounds)),
                                                           2.2 + 0.55 * len(altitudes) * len(flows) + 1.6), squeeze=False)
    for i, flow in enumerate(flows):
        for j, bound in enumerate(bounds):
            ax = axes[i][j]
            for yi, h in enumerate(altitudes):
                for xi, d in enumerate(diameters):
                    c = cells.get((flow, bound, h, d))
                    if c is None:
                        ax.add_patch(Rectangle((xi - 0.5, yi - 0.5), 1, 1, fc="none", ec=MUTED, lw=0.3))
                        continue
                    nan = not math.isfinite(c["value"])
                    face = "#d9d8d2" if nan else SEQ_BLUE(norm(max(c["value"], VALUE_FLOOR_S)))
                    ax.add_patch(Rectangle((xi - 0.5, yi - 0.5), 1, 1, fc=face, ec="white", lw=1.0,
                                           hatch=KIND_HATCH.get(c["kind"], "")))
                    text = "" if nan else f"{c['value']:.2g}"
                    dark = (not nan) and norm(max(c["value"], VALUE_FLOOR_S)) > 0.55
                    ax.text(xi, yi, text, ha="center", va="center", fontsize=8, color="white" if dark else INK,
                            bbox=None if nan else dict(boxstyle="round,pad=0.1", fc=face, ec="none", alpha=0.75))
                    if c["flags"]:
                        ax.text(xi - 0.45, yi + 0.42, "".join(c["flags"]), ha="left", va="top", fontsize=7,
                                color=INK, fontweight="bold")
            ax.set_xlim(-0.5, len(diameters) - 0.5); ax.set_ylim(-0.5, len(altitudes) - 0.5)
            ax.set_xticks(range(len(diameters))); ax.set_xticklabels([f"{d:g}" for d in diameters])
            ax.set_yticks(range(len(altitudes))); ax.set_yticklabels([f"{h:g}" for h in altitudes])
            ax.set_title(f"wake_flow {flow}, {bound} bound", color=SECOND, fontsize=10)
            if i == len(flows) - 1:
                ax.set_xlabel("sphere diameter [mm]")
            if j == 0:
                ax.set_ylabel("release altitude [km]")
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=SEQ_BLUE), ax=axes.ravel().tolist(), shrink=0.7, pad=0.02)
    cbar.set_label("time in the plume [s], log scale\n(a zero is drawn at the floor)", fontsize=9)
    kinds_present = [k for k in INTEGRATED if any(c["kind"] == k for c in cells.values())]
    handles = [Patch(fc="#e9e8e2", ec=MUTED, hatch=KIND_HATCH[k] or None, label=lbl) for k, lbl in (
        ("far_wake", "far_wake release: residence time (time to the first exit)"),
        ("escape", "escape release: total time inside the plume"),
        ("surface", "surface release (rarefied follower): total time inside the plume")) if k in kinds_present]
    handles.append(Patch(fc="#d9d8d2", ec=MUTED, label="not integrated"))
    present = sorted({f for c in cells.values() for f in c["flags"]})
    handles += [Patch(fc="none", ec="none", label=f"{f}: {CELL_FLAG_TEXT[f]}") for f in present]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8, frameon=False, bbox_to_anchor=(0.45, -0.02 - 0.03 * len(handles)))
    fig.suptitle(f"Time spent in the parent's plume, release angle {use:g} degrees\n"
                 f"(far-wake releases: residence time; escape and surface releases: total time inside)", color=INK, fontsize=11, y=1.06)
    return save(fig, out, "residence_map.png")


# -------------------------------------------------------------------------------------------------- affix zones

def affix_flag_names(rows):
    """The yes/no validity flags of affix_regimes.csv (the `flag_` columns that hold True or False)."""
    return [k for k in rows[0] if k.startswith("flag_") and all(r[k] in ("True", "False") for r in rows)]


def affix_panels(rows):
    """[(title, rows)]: one panel when the zone and the validity flags of every (diameter, theta) are the same at every
    altitude, else one per altitude (highest first)."""
    flags = affix_flag_names(rows)
    by_alt = {}
    for r in rows:
        by_alt.setdefault(num(r["altitude_km"]), {})[(num(r["diameter_mm"]), num(r["theta_deg"]))] = (
            r["zone"], tuple(r[k] for k in flags))
    altitudes = sorted(by_alt, reverse=True)
    if all(by_alt[a] == by_alt[altitudes[0]] for a in altitudes):
        return [("all altitudes", [r for r in rows if num(r["altitude_km"]) == altitudes[0]])]
    return [(f"{a:g} km", [r for r in rows if num(r["altitude_km"]) == a]) for a in altitudes]


def plot_affix_zones(out):
    path = os.path.join(out, "affix_regimes.csv")
    if not os.path.exists(path):
        print("affix_zones.png skipped: affix_regimes.csv not found (run tables first)")
        return None
    rows = read_rows(path)
    if not rows:
        print("affix_zones.png skipped: affix_regimes.csv has a header and no rows")
        return None
    flags = affix_flag_names(rows)
    varying = [k for k in flags if len({r[k] for r in rows}) > 1]
    panels = affix_panels(rows)
    diameters = sorted({num(r["diameter_mm"]) for r in rows})
    thetas = sorted({num(r["theta_deg"]) for r in rows})
    edges = [thetas[0] - 0.5 * (thetas[1] - thetas[0]) if len(thetas) > 1 else thetas[0] - 5.0]
    edges += [0.5 * (a + b) for a, b in zip(thetas, thetas[1:])]
    edges.append(thetas[-1] + 0.5 * (thetas[-1] - thetas[-2]) if len(thetas) > 1 else thetas[0] + 5.0)
    ncol = min(len(panels), 3)
    nrow = math.ceil(len(panels) / ncol)
    fig, axes = new_axes(nrow, ncol, figsize=(4.6 * ncol + 1.8, 0.42 * len(diameters) * nrow + 2.4 * nrow), squeeze=False,
                         sharex=True, sharey=True)
    for ax in axes.ravel()[len(panels):]:
        ax.set_visible(False)
    letters = {k: AFFIX_FLAG_CODES.get(k, k[len("flag_"):][:2]) for k in varying}
    for ax, (title, rs) in zip(axes.ravel(), panels):
        for r in rs:
            xi, yi = thetas.index(num(r["theta_deg"])), diameters.index(num(r["diameter_mm"]))
            ax.add_patch(Rectangle((edges[xi], yi - 0.5), edges[xi + 1] - edges[xi], 1, fc=ZONE_COLOR.get(r["zone"], MUTED),
                                   ec="white", lw=1.0))
            mark = "".join(letters[k] + " " for k in varying if r[k] == "True").strip()
            if mark:
                ax.text(0.5 * (edges[xi] + edges[xi + 1]), yi, mark, ha="center", va="center", fontsize=6.5, color="white")
        bif = {num(r["diameter_mm"]): num(r["theta_bif"]) for r in rs}
        ys = [diameters.index(d) for d in diameters if d in bif]
        for d in diameters:                                    # theta_bif of each diameter as a bar across its row
            if d in bif:
                yi = diameters.index(d)
                ax.plot([bif[d], bif[d]], [yi - 0.5, yi + 0.5], color=INK, lw=2.2)
        ax.plot([bif[d] for d in diameters if d in bif], ys, color=INK, lw=0.8, ls="--")
        ax.set_xlim(edges[0], edges[-1]); ax.set_ylim(-0.5, len(diameters) - 0.5)
        ax.set_yticks(range(len(diameters))); ax.set_yticklabels([f"{d:g}" for d in diameters])
        ax.set_xticks(thetas); ax.set_xticklabels([f"{t:g}" for t in thetas])
        ax.set_title(title, color=SECOND, fontsize=10)
    for ax in axes[-1]:
        ax.set_xlabel("release angle theta [degrees]")
    for ax in axes[:, 0]:
        ax.set_ylabel("sphere diameter [mm]")
    zones = [z for z in ZONE_COLOR if any(r["zone"] == z for r in rows)] + sorted(
        {r["zone"] for r in rows} - set(ZONE_COLOR))
    handles = [Patch(fc=ZONE_COLOR.get(z, MUTED), label=z) for z in zones]
    handles.append(Line2D([], [], color=INK, lw=2.2, label="theta_bif"))
    handles += [Patch(fc="none", ec="none", label=f"{letters[k]}: {k[len('flag_'):]}") for k in varying]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, frameon=False)
    note = ("validity flags differ between altitudes: one panel per altitude" if len(panels) > 1
            else "validity flags are the same at every altitude")
    fig.suptitle(f"AFFIx release zone over diameter and release angle ({note})", color=INK, fontsize=11)
    return save(fig, out, "affix_zones.png")


# ---------------------------------------------------------------------------------------- shock-layer paths

def parent_frame_path(hist):
    """(x, y) [m] of a case's follower in the parent's frame at the history's times: x along the parent's direction of
    flight from its centroid (positive forward), y the distance from the flight axis. The direction of flight is the
    direction of the parent's position's time derivative (central differences of the history's parent positions)."""
    t = hist["t"]
    P = np.column_stack([hist["parent_x"], hist["parent_y"], hist["parent_z"]])
    F = np.column_stack([hist["follower_x"], hist["follower_y"], hist["follower_z"]])
    v = np.gradient(P, t, axis=0)
    e_x = v / np.linalg.norm(v, axis=1)[:, None]
    d = F - P
    x = np.einsum("ij,ij->i", d, e_x)
    y = np.linalg.norm(d - x[:, None] * e_x, axis=1)
    return x, y


def rebuild_shock_curve(out, name, h_km):
    """(R [m], shock shape or None, note): the body radius from the case's frontal area and, rebuilt with
    `ShockLayerField.build` from the parent's free stream at release and the case's own shock-layer settings, the
    field's shock curve x_s(y). The shape is None, with the reason in the note, when the rebuild fails."""
    R = 0.5 * math.sqrt(4.0 * PARENT_FRONTAL_AREA / math.pi)      # the default parent, until the case's own is read
    try:
        cfg = load_case_json(out, name)["config"]["settings"]
        R = 0.5 * math.sqrt(4.0 * cfg["frontal_area"] / math.pi)
        from reentry_model.parent_follower_separation_interactions import scenarios
        from reentry_model.parent_follower_separation_interactions.analysis.sweep import get_twin
        from reentry_model.parent_follower_separation_interactions.shock_layer import (
            FreeStream, ShockLayerField, ShockLayerSettings)
        rel = scenarios._release(get_twin(cfg["frontal_area"], cfg["beta_scale"]), h_km)
        e = rel.e
        fs = FreeStream(rho=e.rho, T=e.T, p=e.p, U=e.U, mfp=e.mfp, m_bar=rel.m_bar)
        field = ShockLayerField.build(fs, 0.5 * e.d, e.r0, ShockLayerSettings(**cfg["shock"]), air=scenarios._air())
        return 0.5 * e.d, field.shape, ""
    except Exception as exc:                                  # the plot is still useful without the curve
        return R, None, f"shock curve not drawn ({type(exc).__name__}: {exc})"


def plot_shock_layer_paths(out, rows, h, bound, wake_flow):
    """shock_layer_paths_<h>km.png; None (with a note) when the altitude has no rarefied (surface-release) follower."""
    cases = [r for r in integrated(rows, h, bound) if r["release_kind"] == "surface"]
    flow = pick([f for f in WAKE_FLOW_MODES if any(r["wake_flow"] == f for r in cases)], wake_flow)
    cases = sorted((r for r in cases if r["wake_flow"] == flow), key=lambda r: (r["diameter_mm"], r["theta_deg"]))
    if not cases:
        print(f"shock_layer_paths_{h:g}km.png skipped: no rarefied follower at {h:g} km (bound {bound})")
        return None
    R, shape, note = rebuild_shock_curve(out, cases[0]["name"], h)
    thetas = sorted({r["theta_deg"] for r in cases})
    styles = ["-", "--", "-.", ":", (0, (5, 1, 1, 1)), (0, (1, 3))]
    colors = diameter_colors(r["diameter_mm"] for r in cases)
    fig, (ax_near, ax_far) = new_axes(1, 2, figsize=(13.5, 5.4))
    kinds = set()
    for r in cases:
        hist = load_history(out, r["name"])
        if hist is None:
            continue
        x, y = parent_frame_path(hist)
        t_exit = r["t_shock_layer_exit"]
        keep = hist["t_since_release"] <= max(PATH_WINDOW_S, t_exit if math.isfinite(t_exit) else 0.0)
        style = styles[thetas.index(r["theta_deg"]) % len(styles)]
        for ax in (ax_near, ax_far):
            ax.plot(x[keep], y[keep], style, color=colors[r["diameter_mm"]], lw=1.2)
            ax.plot([x[0]], [y[0]], ".", color=colors[r["diameter_mm"]], ms=4)
            if r["shock_layer_exit"] in EXIT_MARKER and math.isfinite(r["x_exit_parent_frame"]):
                ax.plot([r["x_exit_parent_frame"]], [r["y_exit_parent_frame"]], EXIT_MARKER[r["shock_layer_exit"]],
                        color=colors[r["diameter_mm"]], mec=INK, ms=7, zorder=5)
        kinds.add(r["shock_layer_exit"])
    for ax in (ax_near, ax_far):
        ax.add_patch(Circle((0.0, 0.0), R, fc="#d9d8d2", ec=INK, lw=1.0, zorder=3))
        ax.axvline(0.0, color=MUTED, lw=0.9, ls="--", zorder=2)
        ax.set_xlabel("x, along the parent's flight direction from its centroid [m]")
    ax_near.set_ylabel("y, distance from the flight axis [m]")
    if shape is not None:
        ys = np.linspace(0.0, SHOCK_VIEW_RADII[2] * R, 200)
        xs = np.array([shape.x_s(v) for v in ys])
        ax_near.plot(xs, ys, color="#d55e00", lw=1.6, zorder=4)
        ax_far.plot(xs, ys, color="#d55e00", lw=1.6, zorder=4)
    ax_near.set_xlim(SHOCK_VIEW_RADII[0] * R, SHOCK_VIEW_RADII[1] * R)
    ax_near.set_ylim(0.0, SHOCK_VIEW_RADII[2] * R)
    ax_near.set_aspect("equal", adjustable="box")
    ax_near.set_title("near the body", color=SECOND, fontsize=10)
    ax_far.set_ylim(bottom=0.0)
    ax_far.set_aspect("equal", adjustable="box")
    ax_far.set_title("over the first second after release (or to the shock-layer exit)", color=SECOND, fontsize=10)
    handles = [Line2D([], [], color=colors[d], lw=2) for d in sorted(colors)]
    labels = [f"{d:g} mm" for d in sorted(colors)]
    handles += [Line2D([], [], color=INK, ls=styles[i % len(styles)]) for i in range(len(thetas))]
    labels += [f"theta {t:g} degrees" for t in thetas]
    handles += [Line2D([], [], color="#d9d8d2", marker=EXIT_MARKER[k], ls="", mec=INK, ms=7) for k in EXIT_MARKER if k in kinds]
    labels += [f"exit: {k}" for k in EXIT_MARKER if k in kinds]
    handles += [Patch(fc="#d9d8d2", ec=INK)] + ([Line2D([], [], color="#d55e00", lw=1.6)] if shape is not None else [])
    labels += [f"parent (radius {R:.4g} m)"] + (["shock curve of the field"] if shape is not None else [])
    ax_far.legend(handles, labels, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.suptitle(f"Rarefied followers' paths in the parent's frame, {h:g} km release, wake_flow {flow}, {bound} bound"
                 + (f"\n{note}" if note else ""), color=INK, fontsize=11)
    return save(fig, out, f"shock_layer_paths_{h:g}km.png")


# ------------------------------------------------------------------------------------------------------ command

def _parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], prog="plots")
    p.add_argument("--out", default=OUTPUT, help="folder the sweep and the tables wrote to; plots go to <out>/plots")
    p.add_argument("--theta", type=float, default=120.0, help="release angle of the residence map, degrees")
    p.add_argument("--bound", choices=("low", "high"), default="low", help="bound of the margin, trajectory and path plots")
    p.add_argument("--trajectory-wake-flow", choices=WAKE_FLOW_MODES, default="off", dest="trajectory_wake_flow",
                   help="wake_flow of the trajectory and path plots (the nearest present is used if absent)")
    return p


def main(argv=None):
    args = _parser().parse_args(argv)
    out = os.fspath(args.out)
    rows = load_summary(out)
    if rows is None:
        print(f"no summary.csv under {out}: run sweep first")
        return 1
    print(f"{len(rows)} runs found in {os.path.join(out, 'summary.csv')}")
    failures = []

    def run(label, fn, *a):
        """One plot; an exception is reported and the other plots are still written."""
        try:
            fn(*a)
        except Exception as exc:
            failures.append(label)
            print(f"FAILED {label}: {type(exc).__name__}: {exc}")

    run("plume_size.png", plot_plume_size, out)
    altitudes = sorted({r["altitude_km"] for r in rows}, reverse=True)
    bound = pick(sorted({r["bound"] for r in rows}), args.bound)
    for h in altitudes:
        run(f"margin_{h:g}km.png", plot_margin, out, rows, h, bound)
        run(f"trajectories_{h:g}km.png", plot_trajectories, out, rows, h, bound, args.trajectory_wake_flow)
        run(f"shock_layer_paths_{h:g}km.png", plot_shock_layer_paths, out, rows, h, bound, args.trajectory_wake_flow)
    run("residence_map.png", plot_residence_map, out, rows, args.theta)
    run("affix_zones.png", plot_affix_zones, out)
    if failures:
        print(f"{len(failures)} plot(s) failed: {', '.join(failures)}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
