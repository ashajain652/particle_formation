"""Four-panel figure for the single-sphere trajectory."""
import csv
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter

import particle_trajectory as pt
from reentry_model.atmosphere import US76TableAtmosphere

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
BLUE, ORANGE = "#2a78d6", "#eb6834"

src = sys.argv[1] if len(sys.argv) > 1 else "trajectory_full.csv"
rows = list(csv.DictReader(open(src)))
for r in rows:
    for k, v in r.items():
        if k != "regime":
            r[k] = float(v)
t = np.array([r["t_s"] for r in rows])
V = np.array([r["V"] for r in rows])
h = np.array([r["h_m"] for r in rows]) / 1e3
cd = np.array([r["C_D"] for r in rows])
kn = np.array([r["Kn"] for r in rows])
sr = np.array([r["speed_ratio"] for r in rows])
tp = np.maximum(t, 1e-6)                       # for the log time axis

atm = US76TableAtmosphere()
part = pt.Particle(10e-6, 2813.0)
vt = np.array([pt.terminal_velocity(part, max(x, 0.0) * 1e3, atm) for x in h])

fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4), facecolor=SURFACE)
for ax in axes.ravel():
    ax.set_facecolor(SURFACE)
    ax.grid(True, which="major", color=GRID, lw=0.7, zorder=0)
    ax.grid(True, which="minor", color=GRID, lw=0.35, alpha=0.6, zorder=0)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        ax.spines[s_].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9, which="both")
    ax.yaxis.set_minor_formatter(NullFormatter())

# (a) speed against the local terminal velocity
ax = axes[0, 0]
ax.set_xscale("log"); ax.set_yscale("log")
ax.plot(tp, V, color=BLUE, lw=2.0, zorder=3, label="speed")
ax.plot(tp, vt, color=ORANGE, lw=2.0, ls=(0, (5, 2)), zorder=3,
        label="terminal velocity at that altitude")
ax.set_xlabel("time since release  [s]", color=INK2, fontsize=10)
ax.set_ylabel("speed  [m/s]", color=INK2, fontsize=10)
ax.set_title("(a)  Speed collapses onto terminal velocity in seconds",
             color=INK, fontsize=10.5, loc="left", pad=8)
leg = ax.legend(frameon=False, fontsize=8.5, loc="lower left", labelcolor=INK2, handlelength=2.2)

# (b) altitude
ax = axes[0, 1]
ax.set_xscale("log")
ax.plot(tp, h, color=BLUE, lw=2.0, zorder=3)
ax.set_xlabel("time since release  [s]", color=INK2, fontsize=10)
ax.set_ylabel("altitude  [km]", color=INK2, fontsize=10)
ax.set_title("(b)  Altitude: 34 days to the ground", color=INK, fontsize=10.5, loc="left", pad=8)
for sec, lab in ((3600, "1 h"), (86400, "1 d"), (86400 * 10, "10 d")):
    ax.axvline(sec, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=2)
    ax.annotate(lab, xy=(sec, 78), xytext=(3, 0), textcoords="offset points",
                color=INK2, fontsize=8, ha="left")

# (c) drag coefficient against speed ratio, with the two limits
ax = axes[1, 0]
ax.set_xscale("log"); ax.set_yscale("log")
ss = np.logspace(-5, 1.6, 400)
ax.plot(ss, [pt.cd_free_molecular(x) for x in ss], color=ORANGE, lw=1.6, ls=(0, (5, 2)),
        zorder=2, label="Schaaf–Chambre free-molecule limit")
ax.plot(sr, cd, color=BLUE, lw=2.2, zorder=3, label="trajectory")
ax.set_xlabel("molecular speed ratio  $s = V/\\sqrt{2RT}$", color=INK2, fontsize=10)
ax.set_ylabel("drag coefficient  $C_D$", color=INK2, fontsize=10)
ax.set_title("(c)  $C_D$ rises from 2.07 to 43 000 as the sphere slows",
             color=INK, fontsize=10.5, loc="left", pad=8)
ax.legend(frameon=False, fontsize=8.5, loc="upper right", labelcolor=INK2, handlelength=2.2)

# (d) where the flight sits on the rarefaction map
ax = axes[1, 1]
ax.set_xscale("log"); ax.set_yscale("log")
sc = ax.scatter(kn, sr, c=np.maximum(t, 1e-3), cmap="viridis", norm=LogNorm(), s=7, zorder=3)
cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label("time since release  [s]", color=INK2, fontsize=9)
cb.ax.tick_params(colors=INK2, labelsize=8)
cb.outline.set_edgecolor(GRID)
ax.set_xlim(2e3, 5e-3)
ax.set_ylim(5e-6, 6e1)
# the uncalibrated corner needs BOTH a moderate Knudsen number AND a speed ratio above the
# blend floor, so it is a rectangle, not a vertical band
ax.add_patch(plt.Rectangle((20.0, pt.S_LO), 0.1 - 20.0, 6e1 - pt.S_LO,
                           color=ORANGE, alpha=0.14, lw=0, zorder=1))
ax.axhline(pt.S_LO, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=2)
ax.annotate("blend floor $s$ = {:g}".format(pt.S_LO), xy=(1.5e3, pt.S_LO),
            xytext=(0, 5), textcoords="offset points", color=INK2, fontsize=8, ha="left")
ax.annotate("uncalibrated\ncorner", xy=(1.4, 4.0),
            color="#9c4418", fontsize=8.5, ha="center", va="center")
ax.set_xlabel("Knudsen number  $Kn = 2\\lambda/d$", color=INK2, fontsize=10)
ax.set_ylabel("molecular speed ratio  $s$", color=INK2, fontsize=10)
ax.set_title("(d)  The flight misses the uncalibrated corner", color=INK, fontsize=10.5,
             loc="left", pad=8)

fig.suptitle("10 μm AA7075 sphere released at 80 km with $V_x$ = 6.5 km/s, "
             "$V_y$ = +0.5 km/s", color=INK, fontsize=12.5, x=0.008, ha="left", y=0.985)
fig.tight_layout(rect=(0, 0, 1, 0.955))
fig.savefig("particle_trajectory.png", dpi=200, facecolor=SURFACE)
print("wrote particle_trajectory.png")
