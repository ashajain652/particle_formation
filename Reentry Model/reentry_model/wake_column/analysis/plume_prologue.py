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
SCENES = (("parent", 270), ("zoom", 150), ("droplet", 240), ("atoms", 390), ("column", 240))
N_FRAMES = sum(n for _, n in SCENES)
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8b8a85", "#e4e3df"
STAGE_ON, STAGE_OFF = "#7f9fc9", "#e9eef5"
ATOM_COLOUR = {"Al": "#2a78d6", "AlO": "#eb6834", "OAlO": "#1baf7a"}
SETTLED_WINDOW_S = 24e-6             # the atoms scene shows this much real time
N_TRACKS = 500


def compute(case):
    """Everything the opening draws, for one case."""
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
    return dict(case=case, parent=p, droplet=d, release=release, x_release=x_release, t_loss=t_loss,
                L=trail.L_m, r0=hot.r0_m, fractions=dict(hot.fractions), collisions=hot.mean_collisions,
                stop_time=hot.mean_stop_time_s, N_Al=trail.N_Al, T_amb=column.at(trail.z_mean_km).T_K,
                mid_V=mid.V_ms, mid_T=mid.T_d_K, tracks=traced.tracks, mass_al=species["Al"].mass_amu,
                masses={s: species[s].mass_amu for s in TRACKED})


class Fixed:
    """Random but fixed choices shared by all frames."""

    def __init__(self, d):
        rng = np.random.default_rng(20261105)
        self.emit = rng.uniform(0.0, SETTLED_WINDOW_S, len(d["tracks"]))
        self.theta = rng.uniform(0, 2 * np.pi, 160)
        self.phase = rng.random(160)
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
    ax.scatter([x], [h], s=150, marker="s", c=INK, zorder=5)
    ax.annotate("parent", xy=(x, h), xytext=(x + 40, h + 7), fontsize=12.5)
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
    ax.scatter([xc], [d["case"].h0_km], s=150, marker="s", c=INK, zorder=5)
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
    heat = (now["T"] - 900) / 1000
    colour = plt.get_cmap("inferno")(0.25 + 0.65 * min(max(heat, 0), 1))
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


def scene_column(fig, d, fx, u, i):
    ribbon(fig, 2 if u < 0.45 else 3, d)
    ax = fig.add_axes([0.06, 0.18, 0.88, 0.62]); ax.axis("off")
    r0, L = d["r0"], d["L"]
    ax.set_xlim(0, 1); ax.set_ylim(-1, 1)
    z = smooth(u / 0.55)                                   # 0: a 14 cm view of the trail; 1: the whole trail
    view = 0.14 * (1.3 * L / 0.14) ** z
    x0 = 0.05 + 0.0 * z
    length = min(L / view, 1.0) * 0.9 if z > 0 else 0.9
    # the trail seen at this zoom; thickness drawn at true scale until it is under a pixel
    half = max(4.5 * r0 / view, 0.003)                     # 1 y unit is about 335 px and the view spans 1520 px
    s = smooth((u - 0.5) / 0.35)
    ax.add_patch(Rectangle((x0, -half), max(length, 0.02), 2 * half, fc=ATOM_COLOUR["AlO"], ec="none",
                           alpha=0.75 * (1 - s)))
    w = 1.0 - smooth((u - 0.5) / 0.3)
    ax.plot([0.05, 0.95], [-0.72, -0.72], color=INK2, lw=2, alpha=max(w, 0.0))
    ax.text(0.5, -0.82, f"view is {view:,.2f} m across; the trail is {2 * r0 * 1e3:.0f} mm wide and {L / 1e3:.2f} km long",
            ha="center", fontsize=14, color=INK2, alpha=max(w, 0.0))
    if s > 0:
        # the box model's column, the first thing the wake column section draws
        h = 0.34
        col = "#7f9fc9"
        ax.add_patch(Rectangle((0.14, -h), 0.72, 2 * h, fc=col, ec=INK2, lw=1.4, alpha=0.45 * s, zorder=3))
        ax.add_patch(Ellipse((0.86, 0), 0.035, 2 * h, fc=col, ec=INK2, lw=1.4, alpha=0.7 * s, zorder=3))
        ax.add_patch(Ellipse((0.14, 0), 0.035, 2 * h, fc=SURFACE, ec=INK2, lw=1.4, alpha=s, zorder=4))
        rng = np.random.default_rng(3)
        n = 300
        xs, ys = rng.uniform(0.155, 0.85, n), rng.uniform(-0.93, 0.93, n) * h
        names = list(d["fractions"]); pr = np.array([d["fractions"][k] for k in names]); pr = pr / pr.sum()
        ks = rng.choice(len(names), n, p=pr)
        ax.scatter(xs, ys, s=22, c=[ATOM_COLOUR[names[k]] for k in ks], alpha=s, linewidths=0, zorder=5)
        ax.text(0.5, 0.62, "One piece of the trail, treated as a well-mixed column", ha="center", fontsize=19,
                fontweight="bold", alpha=s)
        info = (f"radius {r0 * 1e3:.1f} mm     length {L / 1e3:.2f} km     air temperature {d['T_amb']:.0f} K     "
                f"{d['fractions']['Al']:.0%} Al, {d['fractions']['AlO']:.0%} AlO, {d['fractions']['OAlO']:.1%} OAlO")
        ax.text(0.5, -0.6, info, ha="center", fontsize=15, alpha=s)
        ax.text(0.5, -0.78, "The wake column model's clock starts here: t = 0", ha="center", fontsize=19,
                fontweight="bold", color="#1f5fa8", alpha=s)
    caption(fig, "Zooming back out: the slowed atoms form a thin trail. Each piece of it starts its own clock, and "
                 "the box model follows one such piece for 12 hours.")
