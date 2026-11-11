"""The opening of the plume video: where the wake column sits in a droplet's re-entry. Four scenes zoom in from the
parent's trajectory to the droplet, the aluminium leaving it, and the atoms slowing to the speed of the air, then back
out to the column that the box model starts from.

Everything drawn comes from the models: the parent's DRAMA trajectory, the droplet run (`run_case`), and Stage A's own
atom histories (`slow_down(n_tracks=...)`). `plume_video` calls `compute` and `render`.
"""
import math
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
import numpy as np

from reentry_model.droplet_reentry import parent_trajectory
from reentry_model.droplet_reentry.model import run_case
from reentry_model.wake_column.ambient import AmbientColumn
from reentry_model.wake_column.reactions import load_reactions
from reentry_model.wake_column.settings import Settings
from reentry_model.wake_column.slowing import TRACKED, calibrate, slow_down
from reentry_model.wake_column.species import load_species
from reentry_model.wake_column.trail import trail_from_rows

FPS = 30
SCENES = (("parent", 405), ("zoom", 225), ("droplet", 360), ("atoms", 585), ("column", 675))
N_FRAMES = sum(n for _, n in SCENES)
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8b8a85", "#e4e3df"
STAGE_ON, STAGE_OFF = "#7f9fc9", "#e9eef5"
ATOM_COLOUR = {"Al": "#2a78d6", "AlO": "#eb6834", "OAlO": "#1baf7a", "OAlO2": "#eda100", "AlCO3": "#e87ba4",
               "AlOH": "#008300", "AlOH2": "#6250d6"}
SETTLED_WINDOW_S = 24e-6             # the atoms scene shows this much real time
N_TRACKS = 500
T_COLOUR_K = (900.0, 2100.0)         # droplet temperatures spanned by the colour bar, K
SHADE = (0.25, 0.90)                 # the part of the inferno colour map used for them
PARENT_SIZE = 3600                   # scatter area of the satellite icon, points squared


def satellite_marker():
    """A satellite icon as one matplotlib marker: a body, a strut each side, and a solar panel of three separated
    cells each side, with a small dish on top, tilted. It is a symbol only and is not to scale."""
    def box(x0, y0, x1, y1):
        return [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    boxes = [box(-0.32, -0.32, 0.32, 0.32), box(-0.48, -0.04, -0.32, 0.04), box(0.32, -0.04, 0.48, 0.04),
             box(-0.06, 0.32, 0.06, 0.46)]
    for k in range(3):
        x0 = 0.48 + k * 0.32
        boxes += [box(x0, -0.3, x0 + 0.27, 0.3), box(-x0 - 0.27, -0.3, -x0, 0.3)]
    verts, codes = [], []
    for b in boxes:
        verts += b
        codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
    dish = Path.arc(200, 340).transformed(Affine2D().scale(0.34).translate(0, 0.72))   # the bowl opens upward
    verts += list(dish.vertices) + [dish.vertices[0]]
    codes += [Path.MOVETO] + [Path.CURVE4 if c == Path.CURVE4 else Path.LINETO for c in dish.codes[1:]] + [Path.CLOSEPOLY]
    return Path(verts, codes).transformed(Affine2D().rotate_deg(-25))


SATELLITE = satellite_marker()


def tapered_trail(d, main):
    """The trail at the moment the droplet has laid its length L: at each distance xi behind the droplet, how long ago
    that piece was laid, its radius and its gas composition, read from the box run at that age. The box itself treats
    every piece as one age; this is what its clock implies for the real, tapered trail."""
    q = d["droplet"]
    t_end = float(np.interp(d["L"], q["x"], q["t"]))
    lost = q["mAl"][0] - q["mAl"]
    x_first = float(np.interp(0.001 * lost[-1], lost, q["x"]))        # where aluminium loss really starts
    xi = np.logspace(-3, np.log10(d["L"] - x_first), 500)
    age = t_end - np.interp(d["L"] - xi, q["x"], q["t"])
    lt = np.log10(main["t"])
    la = np.log10(np.maximum(age, main["t"][0]))
    r = np.interp(la, lt, main["r"])
    names = list(main["conc"])
    shares = np.array([np.interp(la, lt, main["frac"][a]) for a in names]).T
    return dict(t_end=t_end, xi=xi, age=age, r=r, shares=shares / shares.sum(axis=1, keepdims=True), names=names)


def compute(case, main):
    """Everything the opening draws, for one case; `main` is the box run's data from plume_video.compute."""
    settings = Settings()
    parent = parent_trajectory.load_parent()
    final, rows = run_case(case.d0_um, case.h0_km, parent=parent)
    trail = trail_from_rows(rows, case.d0_um, case.h0_km, settings.n_emission_samples)
    species = load_species()
    column = AmbientColumn.load(settings.ambient, case.time_of_day)
    channels = calibrate(load_reactions(stage=("hot", "both")), species)
    hot = slow_down(trail.samples, column, settings, species, channels)
    # one emission sample near the middle of the loss, followed atom by atom
    mid = replace(trail.samples[len(trail.samples) // 2], weight=1.0)
    traced = slow_down((mid,), column, replace(settings, n_test_particles=N_TRACKS), species, channels, n_tracks=N_TRACKS)

    br = parent.branch
    p = {k: np.array([r[key] for r in br]) for k, key in (("t", "time"), ("x", "downrange"), ("h", "altitude"),
                                                          ("V", "velocity"))}
    d = {k: np.array([r[key] for r in rows]) for k, key in (("t", "t_s"), ("h", "h_km"), ("V", "V_ms"), ("T", "T_K"),
                                                            ("d", "d_um"), ("mAl", "m_Al_kg"), ("mZn", "m_Zn_kg"),
                                                            ("mMg", "m_Mg_kg"), ("Kn", "Kn"))}
    d["x"] = np.concatenate([[0.0], np.cumsum(0.5 * (d["V"][1:] + d["V"][:-1]) * np.diff(d["t"]))])  # m along the path
    release = float(np.interp(-case.h0_km, -p["h"], p["t"]))
    x_release = float(np.interp(release, p["t"], p["x"]))
    t_loss = float(np.interp(0.999 * (d["mAl"][0] - d["mAl"][-1]), d["mAl"][0] - d["mAl"], d["t"]))
    out = dict(case=case, parent=p, droplet=d, release=release, x_release=x_release, t_loss=t_loss,
               L=trail.L_m, r0=hot.r0_m, fractions=dict(hot.fractions), collisions=hot.mean_collisions,
               stop_time=hot.mean_stop_time_s, N_Al=trail.N_Al, T_amb=column.at(trail.z_mean_km).T_K,
               mid_V=mid.V_ms, mid_T=mid.T_d_K, tracks=traced.tracks, mass_al=species["Al"].mass_amu,
               masses={s: species[s].mass_amu for s in TRACKED})
    out["taper"] = tapered_trail(out, main)
    return out


class Fixed:
    """Random but fixed choices shared by all frames."""

    def __init__(self, d):
        rng = np.random.default_rng(20261105)
        self.emit = rng.uniform(0.0, SETTLED_WINDOW_S, len(d["tracks"]))
        self.theta = rng.uniform(0, 2 * np.pi, 160)
        self.phase = rng.random(160)
        # dots along the tapered trail, log-spaced behind the droplet so that some stay in view at every zoom
        tp = d["taper"]
        n = 3200
        self.dot_xi = 10 ** rng.uniform(math.log10(0.01), math.log10(tp["xi"][-1]), n)
        self.dot_y = rng.normal(0.0, 1.0 / math.sqrt(2.0), n) * np.interp(np.log10(self.dot_xi), np.log10(tp["xi"]), tp["r"])
        cum = np.cumsum(tp["shares"], axis=1)
        q = rng.random(n)
        idx = np.searchsorted(tp["xi"], self.dot_xi).clip(0, len(tp["xi"]) - 1)
        self.dot_col = [ATOM_COLOUR[tp["names"][int(np.searchsorted(cum[j], u))]] for j, u in zip(idx, q)]
        self.scene_of, self.local, self.length = [], [], []
        for name, n in SCENES:
            self.scene_of += [name] * n
            self.local += list(range(n))
            self.length += [n] * n


def smooth(x):
    x = min(max(x, 0.0), 1.0)
    return x * x * (3.0 - 2.0 * x)


def fmt_t(t):
    if t < 1e-3:
        return f"{t * 1e6:.1f} µs"
    if t < 1:
        return f"{t * 1e3:.0f} ms"
    if t < 120:
        return f"{t:.2f} s"
    return f"{t / 60:.1f} min"


def ribbon(fig, active, d):
    ax = fig.add_axes([0.012, 0.865, 0.976, 0.085])
    ax.set_xlim(0, 4); ax.set_ylim(0, 1); ax.axis("off")
    dur = d["release"] - d["parent"]["t"][0]
    stages = (("1  Parent re-entry", f"DRAMA trajectory, {dur / 60:.0f} minutes from 111 to {d['case'].h0_km:g} km"),
              ("2  Droplet heats and vaporises", f"droplet model, Al lost within about {d['t_loss']:.1f} s"),
              ("3  Atoms slow to the air's speed", f"Stage A Monte Carlo, about {d['stop_time'] * 1e6:.0f} µs per atom"),
              ("4  Wake column", "box model, 12 hours of chemistry and spreading"))
    for i, (title, note) in enumerate(stages):
        on = i == active
        ax.add_patch(FancyBboxPatch((i + 0.012, 0.06), 0.976, 0.88, boxstyle="round,pad=0,rounding_size=0.03",
                                    mutation_aspect=1 / 4 * 4, fc=STAGE_ON if on else STAGE_OFF, ec="none"))
        ax.text(i + 0.03, 0.66, title, fontsize=14, fontweight="bold" if on else "normal", color=INK if on else INK2,
                va="center")
        ax.text(i + 0.03, 0.28, note, fontsize=10.5, color=INK if on else MUTED, va="center")


def caption(fig, text):
    fig.text(0.5, 0.045, text, ha="center", va="center", fontsize=15, color=INK)


def frame(d, fx, i, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": INK})
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor=SURFACE)
    c = d["case"]
    fig.text(0.012, 0.975, f"How the wake column fits into a {c.d0_um:g} µm droplet's re-entry", fontsize=19,
             fontweight="bold", va="center")
    scene, k, n = fx.scene_of[i], fx.local[i], fx.length[i]
    u = k / (n - 1)
    {"parent": scene_parent, "zoom": scene_zoom, "droplet": scene_droplet, "atoms": scene_atoms,
     "column": scene_column}[scene](fig, d, fx, u, i)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    return i


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=10.5)
    ax.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)


def big_axes(fig):
    return fig.add_axes([0.06, 0.13, 0.60, 0.70])


def scene_parent(fig, d, fx, u, i):
    ribbon(fig, 0, d)
    p = d["parent"]
    t = p["t"][0] + (d["release"] - p["t"][0]) * smooth(u)
    ax = big_axes(fig)
    style(ax)
    ax.plot(p["x"], p["h"], color=GRID, lw=2)
    m = p["t"] <= t
    x, h, V = np.interp(t, p["t"], p["x"]), np.interp(t, p["t"], p["h"]), np.interp(t, p["t"], p["V"])
    ax.plot(p["x"][m], p["h"][m], color="#2a78d6", lw=2.5)
    ax.scatter([x], [h], s=PARENT_SIZE, marker=SATELLITE, c=INK, linewidths=0, zorder=5)
    ax.annotate("parent", xy=(x, h), xytext=(x + 90, h + 8), fontsize=12.5)
    ax.set_xlim(p["x"][0] - 50, d["x_release"] + 1500); ax.set_ylim(0, 125)
    ax.set_xlabel("distance along the ground track (km)", fontsize=12, color=INK2)
    ax.set_ylabel("altitude (km)", fontsize=12, color=INK2)
    ax.axhline(d["case"].h0_km, color=MUTED, ls=":", lw=1.2)
    ax.text(p["x"][0], d["case"].h0_km + 1.5, f"droplet release at {d['case'].h0_km:g} km", fontsize=11, color=INK2)
    side = fig.add_axes([0.70, 0.13, 0.285, 0.70]); side.axis("off")
    side.set_xlim(0, 1); side.set_ylim(0, 1)
    for row, (name, val) in enumerate((("Altitude", f"{h:.0f} km"), ("Speed", f"{V:.2f} km/s"),
                                       ("Time from 111 km", fmt_t(t - p["t"][0])))):
        side.text(0, 0.85 - row * 0.15, name, fontsize=14, color=INK2)
        side.text(1, 0.85 - row * 0.15, val, fontsize=20, fontweight="bold", ha="right")
    caption(fig, "A re-entering satellite falls through the atmosphere. Molten droplets are released from it "
                 "at the parent's speed.")


def scene_zoom(fig, d, fx, u, i):
    ribbon(fig, 0 if u < 0.5 else 1, d)
    p = d["parent"]
    w0, w1 = 3500.0, 6.0                     # km across
    w = w0 * (w1 / w0) ** smooth(u)
    ax = fig.add_axes([0.06, 0.13, 0.92, 0.70])
    style(ax)
    xc = d["x_release"]
    ax.plot(p["x"], p["h"], color="#2a78d6", lw=2.5)
    ax.scatter([xc], [d["case"].h0_km], s=PARENT_SIZE, marker=SATELLITE, c=INK, linewidths=0, zorder=5)
    ax.set_xlim(xc - 0.75 * w, xc + 0.25 * w)
    hw = max(0.5 * w / 35.0, 0.5)
    ax.set_ylim(max(d["case"].h0_km - 1.1 * hw, 0), min(d["case"].h0_km + 1.1 * hw, 125))
    ax.set_xlabel("distance along the ground track (km)", fontsize=12, color=INK2)
    ax.set_ylabel("altitude (km)", fontsize=12, color=INK2)
    ax.text(0.5, 0.92, f"Zooming in: {w:,.0f} km across" if w >= 10 else f"Zooming in: {w:.1f} km across",
            transform=ax.transAxes, ha="center", fontsize=15, color=INK2)
    caption(fig, f"Zooming in to the release point at {d['case'].h0_km:g} km, where one droplet leaves the parent "
                 f"at {d['droplet']['V'][0] / 1e3:.1f} km/s.")


def droplet_at(d, t):
    q = d["droplet"]
    return {k: float(np.interp(t, q["t"], q[k])) for k in ("x", "h", "V", "T", "d", "mAl", "mZn", "mMg", "Kn")}


def temperature_shade(T):
    """The droplet's temperature [K] as a position on the inferno colour map, limited to the bar's range."""
    heat = (T - T_COLOUR_K[0]) / (T_COLOUR_K[1] - T_COLOUR_K[0])
    return SHADE[0] + (SHADE[1] - SHADE[0]) * min(max(heat, 0.0), 1.0)


def temperature_colourbar(fig, T):
    """A bar of the colours the droplet takes, with a line at its temperature T [K] now."""
    cb = fig.add_axes([0.022, 0.15, 0.012, 0.36])
    cb.imshow(np.linspace(*SHADE, 256)[:, None], origin="lower", aspect="auto", cmap="inferno", vmin=0, vmax=1,
              extent=(0, 1, *T_COLOUR_K))
    cb.axhline(min(max(T, T_COLOUR_K[0]), T_COLOUR_K[1]), color="white", lw=2.5)
    cb.axhline(min(max(T, T_COLOUR_K[0]), T_COLOUR_K[1]), color=INK, lw=1)
    cb.set_xticks([]); cb.yaxis.tick_right()
    cb.set_yticks(np.arange(T_COLOUR_K[0], T_COLOUR_K[1] + 1, 200))
    cb.tick_params(colors=INK2, labelsize=10.5)
    cb.set_ylabel("droplet temperature (K)", fontsize=11.5, color=INK2, labelpad=4)
    cb.yaxis.set_label_position("left")


def scene_droplet(fig, d, fx, u, i):
    ribbon(fig, 1, d)
    q = d["droplet"]
    tend = max(d["t_loss"] * 1.15, 0.2)
    t = tend * u
    now = droplet_at(d, t)
    al0, zn0, mg0 = q["mAl"][0], q["mZn"][0], q["mMg"][0]
    al_total = al0 - q["mAl"][-1]
    # the droplet's path: the trail it lays, on stretched axes
    ax = fig.add_axes([0.06, 0.60, 0.60, 0.23])
    style(ax)
    tt = np.linspace(0, t, 300) if t > 0 else np.array([0.0])
    px, ph = np.interp(tt, q["t"], q["x"]) / 1e3, (np.interp(tt, q["t"], q["h"]) - q["h"][0]) * 1e3
    loss = np.gradient(np.interp(tt, q["t"], q["mAl"])) if len(tt) > 1 else np.zeros(1)
    ax.scatter(px, ph, s=8 + 60 * np.abs(loss) / max(np.abs(loss).max(), 1e-30), c=ATOM_COLOUR["Al"], alpha=0.55,
               linewidths=0)
    ax.scatter([now["x"] / 1e3], [(now["h"] - q["h"][0]) * 1e3], s=70, c=INK, zorder=5)
    xmax = float(np.interp(tend, q["t"], q["x"])) / 1e3
    ax.set_xlim(-0.05, xmax * 1.05)
    ymin = float(np.interp(tend, q["t"], q["h"]) - q["h"][0]) * 1e3
    ax.set_ylim(ymin * 1.15 - 2, 4)
    ax.set_xlabel("distance flown since release (km)", fontsize=11.5, color=INK2)
    ax.set_ylabel("height change (m)", fontsize=11.5, color=INK2)
    if now["x"] >= d["L"]:
        ax.annotate("", xy=(d["L"] / 1e3, ymin * 0.75), xytext=(0, ymin * 0.75),
                    arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.4))
        ax.text(d["L"] / 2e3, ymin * 0.70, f"trail length {d['L'] / 1e3:.2f} km: 90 % of the Al is lost here",
                ha="center", fontsize=11.5, color=INK2)
    # the droplet itself, with atoms leaving it
    sp = fig.add_axes([0.06, 0.12, 0.27, 0.42]); sp.set_xlim(-1.2, 1.2); sp.set_ylim(-1.2, 1.2); sp.set_aspect("equal")
    sp.axis("off")
    rad = 0.5 * now["d"] / d["case"].d0_um
    colour = plt.get_cmap("inferno")(temperature_shade(now["T"]))
    sp.add_patch(Circle((0, 0), rad, fc=colour, ec=INK, lw=1))
    rate = np.interp(t, q["t"], -np.gradient(q["mAl"], q["t"]))
    peak = np.max(-np.gradient(q["mAl"], q["t"]))
    n = int(160 * min(max(rate / peak, 0), 1) ** 0.5)
    for j in range(n):
        a = fx.theta[j]
        r = rad + 0.72 * ((fx.phase[j] + i * 0.03) % 1.0)
        sp.scatter([r * math.cos(a)], [r * math.sin(a)], s=14, c=ATOM_COLOUR["Al"], linewidths=0)
    sp.text(0, -1.15, f"droplet {now['d']:.0f} µm across, {now['T']:.0f} K", ha="center", fontsize=13)
    sp.text(0, 1.1, "Al atoms leaving the surface", ha="center", fontsize=12, color=INK2)
    temperature_colourbar(fig, now["T"])
    # the droplet's state
    for pos, (y, lab, col, marker) in enumerate(((q["T"], "droplet temperature (K)", "#eb6834", now["T"]),
                                                 (q["V"] / 1e3, "speed (km/s)", "#2a78d6", now["V"] / 1e3),
                                                 (100 * (al0 - q["mAl"]) / al0, "aluminium vaporised (%)", "#6250d6",
                                                  100 * (al0 - now["mAl"]) / al0))):
        a2 = fig.add_axes([0.725, 0.645 - pos * 0.255, 0.26, 0.185])
        style(a2)
        a2.plot(q["t"], y, color=col, lw=2)
        a2.scatter([t], [marker], s=50, c=col, zorder=5)
        a2.set_xlim(0, tend); a2.set_ylim(0, max(y[q["t"] <= tend]) * 1.12)
        a2.set_title(lab, fontsize=12, color=INK2, loc="left", pad=4)
    fig.axes[-1].set_xlabel("time since release (s)", fontsize=11.5, color=INK2)
    sp2 = fig.add_axes([0.345, 0.12, 0.33, 0.42]); sp2.axis("off"); sp2.set_xlim(0, 1); sp2.set_ylim(0, 1)
    rows = (("Time since release", fmt_t(t)), ("Aluminium lost", f"{100 * (al0 - now['mAl']) / al0:.0f} %"),
            ("Zinc lost", f"{100 * (zn0 - now['mZn']) / zn0:.0f} %"), ("Magnesium lost", f"{100 * (mg0 - now['mMg']) / mg0:.0f} %"),
            ("Gas mean free path / droplet size", f"{now['Kn']:.1f}"))
    for row, (name, val) in enumerate(rows):
        sp2.text(0.02, 0.86 - row * 0.17, name, fontsize=13, color=INK2)
        sp2.text(0.98, 0.86 - row * 0.17, val, fontsize=16, fontweight="bold", ha="right")
    caption(fig, f"Air heats the droplet to about {q['T'].max():,.0f} K in a fraction of a second and aluminium "
                 "evaporates from its surface as it flies, laying a trail of vapour.")


def track_states(d, fx, clock):
    """Position (m) in the air's frame, and species, of each traced atom at time `clock` since the window began."""
    xs, ys, cols = [], [], []
    V = d["mid_V"]
    for tr, te in zip(d["tracks"], fx.emit):
        tau = clock - te
        if tau < 0:
            continue
        j = int(np.searchsorted(tr[:, 0], tau)) - 1
        j = min(max(j, 0), len(tr) - 1)
        if j == len(tr) - 1:
            x, y = tr[j, 1], tr[j, 2]
        else:
            f = (tau - tr[j, 0]) / (tr[j + 1, 0] - tr[j, 0])
            x, y = tr[j, 1] + f * (tr[j + 1, 1] - tr[j, 1]), tr[j, 2] + f * (tr[j + 1, 2] - tr[j, 2])
        xs.append(V * te + x); ys.append(y)
        cols.append(ATOM_COLOUR[TRACKED[int(tr[j, 4])]])
    return np.array(xs), np.array(ys), cols


def kinetic_temperature(d, tt):
    """Mean kinetic energy of the traced atoms against time since emission, as a temperature (K), and the share of
    atoms still moving faster than the stop criterion."""
    k_b = 1.380649e-23
    out = np.zeros(len(tt))
    for q, tau in enumerate(tt):
        e = [tr[min(max(int(np.searchsorted(tr[:, 0], tau)) - 1, 0), len(tr) - 1), 5] for tr in d["tracks"]]
        out[q] = 2.0 * float(np.mean(e)) / (3.0 * k_b)
    return out


def scene_atoms(fig, d, fx, u, i):
    ribbon(fig, 2, d)
    clock = SETTLED_WINDOW_S * u
    V = d["mid_V"]
    ax = fig.add_axes([0.06, 0.26, 0.60, 0.57])
    ax.set_xlim(-10, V * SETTLED_WINDOW_S * 1e3 * 1.02); ax.set_ylim(-28, 28); ax.set_aspect("equal", adjustable="box")
    style(ax)
    xs, ys, cols = track_states(d, fx, clock)
    if len(xs):
        ax.scatter(xs * 1e3, ys * 1e3, s=11, c=cols, linewidths=0, alpha=0.9)
    x_drop = V * clock * 1e3
    ax.add_patch(Circle((x_drop, 0), 0.5 * d["case"].d0_um * 1e-3, fc=INK, ec=INK, lw=0.5, zorder=6))
    ax.annotate(f"droplet, to scale: {d['case'].d0_um:g} µm at {V / 1e3:.1f} km/s", xy=(x_drop, 0), xytext=(x_drop - 2, 22),
                fontsize=11.5, ha="right" if x_drop > 60 else "left",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    for s in (1, -1):
        ax.axhline(s * d["r0"] * 1e3, color=INK2, ls="--", lw=1.1)
    ax.text(ax.get_xlim()[0] + 1, d["r0"] * 1e3 + 1, f"r.m.s. radius where atoms stop: {d['r0'] * 1e3:.1f} mm", fontsize=11,
            color=INK2)
    ax.set_xlabel("distance along the droplet's path (mm), in the frame of the still air", fontsize=11.5, color=INK2)
    ax.set_ylabel("across the path (mm)", fontsize=11.5, color=INK2)
    ax.text(0.99, 0.95, f"{clock * 1e6:.1f} µs", transform=ax.transAxes, ha="right", fontsize=22, fontweight="bold")
    # how fast the atoms cool
    a2 = fig.add_axes([0.725, 0.50, 0.26, 0.33]); style(a2)
    tt = np.linspace(0, SETTLED_WINDOW_S, 60)
    T = kinetic_temperature(d, tt)
    a2.semilogy(tt * 1e6, T, color="#eb6834", lw=2.2)
    a2.axhline(d["T_amb"], color=INK2, ls=":", lw=1.3)
    a2.text(0, d["T_amb"] * 1.15, f"air, {d['T_amb']:.0f} K", ha="left", fontsize=10.5, color=INK2)
    a2.axvline(min(clock, SETTLED_WINDOW_S) * 1e6, color=INK, lw=1.2)
    a2.set_title("energy of an average atom, as a temperature (K)", fontsize=12, color=INK2, loc="left", pad=4)
    a2.set_xlabel("time since the atom left the droplet (µs)", fontsize=11.5, color=INK2)
    # what they become
    a3 = fig.add_axes([0.725, 0.13, 0.26, 0.26]); a3.axis("off"); a3.set_xlim(0, 1); a3.set_ylim(0, 1)
    a3.text(0, 0.95, "Where they end up", fontsize=13, fontweight="bold", va="top")
    for row, (s, f) in enumerate(d["fractions"].items()):
        y = 0.68 - row * 0.2
        a3.scatter([0.03], [y], s=90, c=ATOM_COLOUR[s], linewidths=0)
        a3.text(0.08, y, s, fontsize=14, va="center")
        a3.text(0.98, y, f"{f:.1%}", fontsize=15, fontweight="bold", ha="right", va="center")
    a3.text(0, 0.06, f"about {d['collisions']:.0f} collisions and {d['stop_time'] * 1e6:.0f} µs per atom", fontsize=12,
            color=INK2)
    a4 = fig.add_axes([0.06, 0.135, 0.60, 0.1]); a4.axis("off"); a4.set_xlim(0, 1); a4.set_ylim(0, 1)
    x = 0.0
    for s, col in ATOM_COLOUR.items():
        a4.scatter([x + 0.012], [0.75], s=80, c=col, linewidths=0)
        a4.text(x + 0.03, 0.75, s, fontsize=12.5, va="center")
        x += 0.12
    a4.text(0.0, 0.2, "Each dot is a real test atom from the Monte Carlo.", fontsize=12, color=INK2, va="center")
    caption(fig, "Each atom keeps the droplet's speed until collisions with air slow it down. About half its energy "
                 "goes in each collision, and some Al reacts with O₂ to give AlO.")


def fmt_len(x):
    if x < 1:
        return f"{x * 100:.3g} cm" if x >= 0.01 else f"{x * 1e3:.3g} mm"
    return f"{x:.3g} m" if x < 1000 else f"{x / 1e3:.3g} km"


def nice_length(x):
    """The largest of 1, 2 or 5 times a power of ten that does not exceed x."""
    p = 10.0 ** math.floor(math.log10(x))
    return max(m * p for m in (1, 2, 5) if m * p <= x)


def scene_column(fig, d, fx, u, i):
    ribbon(fig, 3, d)
    tp = d["taper"]
    a, b, c = smooth(u / 0.36), smooth((u - 0.38) / 0.22), smooth((u - 0.68) / 0.22)
    span = float(tp["xi"][-1])
    ax = fig.add_axes([0.05, 0.145, 0.90, 0.675]); ax.axis("off")
    aspect = 0.675 * 10.8 / (0.90 * 19.2)
    w = 0.14 * (1.14 * span / 0.14) ** a                   # metres across the view
    half = 0.5 * w * aspect
    rmax = float(tp["r"].max())
    full = 0.72 * 0.5 * (1.14 * span) * aspect / rmax      # the stretch that makes the wide end fill the view
    stretch = full ** b
    ax.set_xlim(-0.93 * w, 0.07 * w); ax.set_ylim(-half / stretch, half / stretch)
    fade = 1.0 - c
    if fade > 0:
        ax.fill_between(-tp["xi"], -tp["r"], tp["r"], color="#7f9fc9", alpha=0.35 * fade, linewidth=0)
        ax.plot(-tp["xi"], tp["r"], color=INK2, lw=1.1, alpha=fade)
        ax.plot(-tp["xi"], -tp["r"], color=INK2, lw=1.1, alpha=fade)
        ax.scatter(-fx.dot_xi, fx.dot_y, s=13, c=fx.dot_col, linewidths=0, alpha=0.9 * fade)
        ax.scatter([0], [0], s=46, c=INK, alpha=fade, zorder=6)
        ax.text(0.003 * w, half / stretch * 0.1, "droplet", fontsize=12, alpha=fade)
        # the ruler: its length, not the trail's width, is what changes on screen
        y = -0.86 * half / stretch
        bar = nice_length(0.22 * w)
        x1 = -0.88 * w
        ax.plot([x1, x1 + bar], [y, y], color=INK, lw=3, alpha=fade, solid_capstyle="butt")
        for xr in (x1, x1 + bar):
            ax.plot([xr, xr], [y - 0.025 * half / stretch, y + 0.025 * half / stretch], color=INK, lw=2, alpha=fade)
        ax.text(x1 + bar / 2, y + 0.05 * half / stretch, fmt_len(bar), ha="center", fontsize=15, fontweight="bold",
                alpha=fade)
    if b < 0.5:
        ax.text(0.5, 0.93, "The camera is pulling back. The trail is not shrinking: only the view is widening.",
                transform=ax.transAxes, ha="center", fontsize=17, fontweight="bold", alpha=(1 - 2 * b) * (0.0 + min(u / 0.03, 1)))
        ax.text(0.5, 0.86, f"view {fmt_len(w)} across, trail about {2 * d['r0'] * 1e3:.0f} mm wide at the droplet, drawn "
                "at true proportions", transform=ax.transAxes, ha="center", fontsize=14, color=INK2, alpha=1 - 2 * b)
    if b > 0 and fade > 0:
        ax.text(0.5, 0.93, f"Same trail, with its width stretched {stretch:,.0f} times so the shape can be seen",
                transform=ax.transAxes, ha="center", fontsize=17, fontweight="bold", alpha=b * fade)
        old_r, old_age = float(tp["r"][-1]), float(tp["age"][-1])
        ax.text(0.07, 0.17, f"oldest end: laid {old_age:.2f} s ago\nradius about {fmt_len(old_r)}", transform=ax.transAxes,
                fontsize=14, alpha=b * fade)
        ax.text(0.93, 0.70, f"newest end: just laid\nradius {d['r0'] * 1e3:.0f} mm", transform=ax.transAxes, ha="right",
                fontsize=14, alpha=b * fade)
        ax.text(0.5, 0.875, "Older vapour has had longer to spread, and its colour shows it has reacted further. Width "
                "grows with the square root of age, so the edges curve.", transform=ax.transAxes, ha="center", fontsize=13,
                color=INK2, alpha=b * fade)
    s = smooth((u - 0.74) / 0.2)
    if s > 0:
        ax2 = fig.add_axes([0.05, 0.145, 0.90, 0.675]); ax2.axis("off"); ax2.set_xlim(0, 1); ax2.set_ylim(-1, 1)
        r0, L = d["r0"], d["L"]
        h = 0.32
        col = "#7f9fc9"
        ax2.add_patch(Rectangle((0.14, -h), 0.72, 2 * h, fc=col, ec=INK2, lw=1.4, alpha=0.45 * s, zorder=3))
        ax2.add_patch(Ellipse((0.86, 0), 0.035, 2 * h, fc=col, ec=INK2, lw=1.4, alpha=0.7 * s, zorder=3))
        ax2.add_patch(Ellipse((0.14, 0), 0.035, 2 * h, fc=SURFACE, ec=INK2, lw=1.4, alpha=s, zorder=4))
        rng = np.random.default_rng(3)
        n = 300
        xs, ys = rng.uniform(0.155, 0.85, n), rng.uniform(-0.93, 0.93, n) * h
        names = list(d["fractions"]); pr = np.array([d["fractions"][k] for k in names]); pr = pr / pr.sum()
        ks = rng.choice(len(names), n, p=pr)
        ax2.scatter(xs, ys, s=22, c=[ATOM_COLOUR[names[k]] for k in ks], alpha=s, linewidths=0, zorder=5)
        ax2.text(0.5, 0.74, "The model's simplification: one column, all the same age", ha="center", fontsize=19,
                 fontweight="bold", alpha=s)
        ax2.text(0.5, 0.52, f"The taper is ignored: the whole trail starts at t = 0 with one radius. Age differences of "
                 f"{tp['age'][-1]:.1f} s matter less and less as the plume ages.", ha="center", fontsize=13.5, color=INK2, alpha=s)
        info = (f"radius {r0 * 1e3:.1f} mm     length {L / 1e3:.2f} km     air temperature {d['T_amb']:.0f} K     "
                f"{d['fractions']['Al']:.0%} Al, {d['fractions']['AlO']:.0%} AlO, {d['fractions']['OAlO']:.1%} OAlO")
        ax2.text(0.5, -0.56, info, ha="center", fontsize=15, alpha=s)
        ax2.text(0.5, -0.76, "The wake column model's clock starts here: t = 0", ha="center", fontsize=19,
                 fontweight="bold", color="#1f5fa8", alpha=s)
    caption(fig, "The slowed atoms form a very thin, very long trail. The box model then follows it for 12 hours.")
