"""One 16:9 slide figure: the pieces of the droplet drag law and where each one applies.

The model's drag (drag_law="singh", the default) is the single correlation of Singh et al. (2021), singh_drag.py:

    C_D = C_D,cont(M, Re) f(Kn, Wr) (1 - w) + C_D,fm(s, T_p/T) w,      w = Br^1.8 / (1 + Br^1.8)

Left, the Knudsen-Mach plane coloured by which piece sets C_D there, with the paths of three droplets of the demise
sweep: solid while braking (the saved history), dashed while settling at terminal speed from the end of braking to the
ground (computed here with DropletSimulator.terminal_velocity, droplet at the air temperature, no wind); right, one card per piece. The pieces are told apart by the free-molecule weight w and the rarefaction
factor f: w > 0.9 free-molecule; w < 0.1 and f > 0.9 continuum (Stokes below Re 1, Abraham's boundary-layer law
below Mach 1, the bow-shock law above); w < 0.1 and f < 0.9 continuum with slip; otherwise the bridge. Wall at the
air temperature.

  output/droplet_reentry_output/demise_sweep/figures/13_drag_law_diagram.{png,pdf,svg}
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from reentry_model.droplet_reentry import model as dd
from reentry_model.droplet_reentry import particle_trajectory as pt
from reentry_model.droplet_reentry import singh_drag as sg
from reentry_model.droplet_reentry.analysis import demise_sweep as rs

FIG = os.path.join(rs.OUT, "figures")
os.makedirs(FIG, exist_ok=True)
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
# region fills (light tints) and their card accents, in the order of REGIONS
FILLS = ["#f3e3c3", "#dce9fb", "#f7d6c8", "#b4d1f5", "#d3eee0", "#e8d8f4"]
ACCENTS = ["#b07a00", "#5598e7", "#d4572a", "#1c5cab", "#16895f", "#8a4cc0"]
PATHS = [((1.0, 90), "#c2410c"), ((20.6, 70), "#0b0b0b"), ((200.0, 50), "#0d6b4a")]
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "font.size": 12, "mathtext.fontset": "cm", "axes.spines.top": False, "axes.spines.right": False})

REGIONS = [   # (map label, card title, equation, when it applies)
    ("Stokes", "Stokes (slow viscous flow)",
     r"$C_D \to 24/Re$",
     "Re < 1, Mach < 1, Kn < about 0.03: slow settling low in the atmosphere"),
    ("Abraham", "Abraham boundary layer",
     r"$C_D = C_0\,\Theta(M)\,(1 + \delta_0/\sqrt{Re})^2$",
     "Re > 1, Mach < 1, Kn < about 0.03: subsonic, thin viscous layer"),
    ("bow shock", "Supersonic continuum (bow shock)",
     r"$C_D = C_1\,(1 - \alpha U_s/U) + C_{D,Abraham}$ (post-shock)",
     "Mach > 1, Kn < about 0.03: curved bow shock; C_D tends to 0.9"),
    ("continuum\n+ slip", "Continuum with slip correction",
     r"$C_D = C_{D,cont}\; f(Kn, Wr)$",
     "Kn about 0.03 to 0.2: air slips at the surface (Cunningham), less drag"),
    ("bridge", "Transitional bridge",
     r"$C_D = (1-w)\,C_{D,cont} f + w\,C_{D,fm}$",
     "Kn about 0.1 to 2: weighted mix of the two limits, weight w set by Br"),
    ("free-molecule", "Free-molecule sphere",
     r"$C_D = C_{D,fm}(s,\;T_p/T)$",
     "Kn > about 2: molecules hit one by one (Schaaf & Chambre); no shock"),
]


def classify(kn, ma):
    """Index into REGIONS of the piece that sets C_D at Knudsen number lambda/D kn and Mach number ma."""
    re = (ma / kn) * math.sqrt(math.pi * sg.GAMMA / 2.0)
    _, g = sg.drag_coefficient(ma, re, 1.0)
    if g["w_fm"] > 0.9:
        return 5
    if g["w_fm"] > 0.1:
        return 4
    if g["f_kn"] < 0.9:
        return 3
    if ma >= 1.0:
        return 2
    return 0 if re < 1.0 else 1


def knudsen(ma, re):
    """lambda / D as singh_drag defines it (Eq. 2), so the paths sit on the map's axis."""
    return (ma / re) * math.sqrt(math.pi * sg.GAMMA / 2.0)


def settling_path(sim, last, step_m=100.0):
    """(Kn, Mach) of a droplet sinking at terminal speed from the end of braking to the ground."""
    lat, lon = math.radians(last.lat_deg), math.radians(last.lon_deg)
    d, m = last.d_um * 1e-6, last.mass_kg
    kn, ma = [], []
    for h in np.arange(last.h_km * 1e3, -1.0, -step_m):
        fs = sim.atmosphere.state(0.0, h, lat, lon)
        v = sim.terminal_velocity(fs, dd.earth.geodetic_to_ecef(h, lat, lon), fs.T, d, m)
        ma.append(v / math.sqrt(sg.GAMMA * pt.R_AIR * fs.T))
        kn.append(knudsen(ma[-1], fs.rho * v * d / pt.viscosity(fs.T)))
    return np.array(kn), np.array(ma)


KN = np.logspace(-4, 5, 280)
MA = np.logspace(-8, np.log10(30.0), 340)
Z = np.array([[classify(k, m) for k in KN] for m in MA])

fig = plt.figure(figsize=(13.333, 7.5))
fig.text(0.04, 0.94, "Droplet drag: one correlation, six regimes", fontsize=22, weight="bold", color=INK)
fig.text(0.04, 0.895, "Singh et al. (2021) sphere drag, used from release to the ground; "
         "the Knudsen and Mach numbers decide which piece sets the drag", fontsize=13, color=INK2)

ax = fig.add_axes([0.065, 0.21, 0.43, 0.64])
ax.pcolormesh(KN, MA, Z, cmap=ListedColormap(FILLS), vmin=-0.5, vmax=5.5, shading="auto", rasterized=True)
re1 = np.sqrt(2.0 / (math.pi * sg.GAMMA)) * KN          # Re = 1 line: Ma = Kn sqrt(2/(pi gamma))
ax.plot(KN[KN < 0.05], re1[KN < 0.05], color=INK2, lw=0.8, ls=":")
ax.text(3e-4, 3e-4 * 0.67 * 1.8, "Re = 1", rotation=25, fontsize=9.5, color=INK2)
ax.axhline(1.0, color=INK2, lw=0.8, ls=":")
ax.text(1.3e-4, 1.25, "Mach 1", fontsize=9.5, color=INK2)
LABELS = ["Stokes", "Abraham", "bow shock", "slip", "bridge", "free-molecule"]
for i, (lab, (k, m)) in enumerate(zip(LABELS, [(2e-4, 1e-7), (2e-4, 0.06), (2e-4, 6.0), (0.095, 1e-6),
                                                (0.75, 1e-6), (150.0, 0.2)])):
    ax.text(k, m, lab, fontsize=11.5, weight="bold", color=INK, va="center", ha="center" if i in (3, 4) else "left",
            rotation=90 if i in (3, 4) else 0)
parent = dd.parent_trajectory.load_parent()
for (d, h), c in PATHS:
    p = os.path.join(rs.OUT, "histories", rs.history_name(d, h))
    if not os.path.exists(p):
        continue
    r = pd.read_csv(p)
    ma = r.Ma.values
    kn = knudsen(ma, r.rho_air.values * r.V_ms.values * r.d_um.values * 1e-6 / [pt.viscosity(t) for t in r.T_air])
    kn_s, ma_s = settling_path(dd.DropletSimulator(dd.parent_release(parent, float(h)), d), r.iloc[-1])
    for x, y, ls in ((kn, ma, "-"), (kn_s, ma_s, (0, (4, 2)))):
        ax.plot(x, y, color=SURFACE, lw=4.5, solid_capstyle="round", dash_capstyle="round")
        ax.plot(x, y, color=c, lw=2.0, ls=ls)
    ax.plot([], [], color=c, lw=2.0, label="{:g} µm released at {} km".format(d, h))
    ax.plot(kn[0], ma[0], "o", ms=8, color=c, mec=SURFACE, mew=1.5)
    ax.plot(kn_s[-1], ma_s[-1], "s", ms=7, color=c, mec=SURFACE, mew=1.5)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(KN[0], KN[-1])
ax.set_ylim(MA[0], MA[-1])
ax.set_xlabel("Knudsen number  (air mean free path / droplet diameter)")
ax.set_ylabel("Mach number")
ax.legend(title="Droplet paths: solid braking, dashed settling (dot release, square ground)", loc="upper left", bbox_to_anchor=(-0.12, -0.12), ncol=3,
          frameon=False, fontsize=10, title_fontsize=10, columnspacing=1.0, handlelength=1.5)
ax.legend_._legend_box.align = "left"

# ---- right: one card per piece
C = fig.add_axes([0.55, 0.06, 0.42, 0.8])
C.axis("off")
C.set_xlim(0, 1)
C.set_ylim(0, 1)
step = 1.0 / len(REGIONS)
for i, (_, title, eq, when) in enumerate(REGIONS):
    y0 = 1.0 - (i + 1) * step
    C.add_patch(plt.Rectangle((0.0, y0 + 0.012), 1.0, step - 0.024, facecolor=FILLS[i], edgecolor="none"))
    C.add_patch(plt.Rectangle((0.0, y0 + 0.012), 0.012, step - 0.024, facecolor=ACCENTS[i], edgecolor="none"))
    C.text(0.035, y0 + step - 0.03, title, fontsize=12, weight="bold", va="top", color=INK)
    C.text(0.035, y0 + step * 0.47, eq, fontsize=12, va="center", color=INK)
    C.text(0.035, y0 + 0.032, when, fontsize=10, va="bottom", color=INK2)

fig.text(0.04, 0.012, "Re Reynolds number;  s speed ratio V/sqrt(2RT);  Br bridging parameter (about "
         "Kn at low Mach);  w free-molecule weight;  f slip factor;  T_p/T droplet over air temperature "
         "(map drawn at 1).",
         fontsize=9.5, color=MUTED)
for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(FIG, "13_drag_law_diagram." + ext), dpi=200)
print(os.path.join(FIG, "13_drag_law_diagram.png"))
