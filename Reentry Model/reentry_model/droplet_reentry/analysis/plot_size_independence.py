"""One 16:9 slide figure: why the share of a droplet lost to vapour does not depend on its initial size.

Left, the scaling argument: drag, heating and evaporation all act on the surface (area ~ d^2) while the droplet's
momentum, heat capacity and mass scale with its volume (~ d^3), so every rate carries the same factor 6/(rho_p d).
Taking the speed V as the independent variable instead of time cancels that factor, and with it every trace of d.
Right, the evidence from output/droplet_reentry_output/demise_sweep/final_states.csv: mass lost against initial
diameter at five release altitudes.

  output/droplet_reentry_output/demise_sweep/figures/11_size_independence_slide.{png,pdf,svg}
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from reentry_model.droplet_reentry.analysis import demise_sweep as rs

RUN = sys.argv[1] if len(sys.argv) > 1 else rs.OUT
FIG = os.path.join(RUN, "figures")
os.makedirs(FIG, exist_ok=True)

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
ACCENT, BOX = "#2a78d6", "#f0efec"
ORDINAL = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]
ALTITUDES = [50, 62, 70, 78, 90]          # km, ordered low to high; ramp light to dark
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0, "font.size": 13, "mathtext.fontset": "cm",
    "axes.axisbelow": True})

f = pd.read_csv(os.path.join(RUN, "final_states.csv"))
f = f[f.skipped.astype(str) != "True"]

fig = plt.figure(figsize=(13.333, 7.5))
fig.text(0.04, 0.93, "The share of a droplet that vaporises does not depend on its size",
         fontsize=22, weight="bold", color=INK)
fig.text(0.04, 0.885, "Every surface process scales the same way with diameter, so the diameter cancels "
         "once speed, not time, is the clock", fontsize=14, color=INK2)

# ---- left: the argument
L = fig.add_axes([0.04, 0.05, 0.47, 0.80])
L.axis("off")
L.set_xlim(0, 1)
L.set_ylim(0, 1)
L.text(0, 0.99, "1   Surface over volume", fontsize=15, weight="bold", color=ACCENT, va="top")
L.text(0.04, 0.905, r"Area $\pi d^2$, mass $\rho_p \pi d^3/6$:  every rate per unit mass carries "
       r"$\dfrac{6}{\rho_p\, d}$", fontsize=14, va="center")
L.text(0, 0.83, "2   The three rates", fontsize=15, weight="bold", color=ACCENT, va="top")
eqs = [(r"$\dfrac{dV}{dt} = -\dfrac{6}{\rho_p d}\;\tau(\rho_a, V)$", "drag"),
       (r"$\dfrac{dT}{dt} = \dfrac{6}{\rho_p d\, c}\;\left[q_{in}(\rho_a, V, T) - q_{rad}(T) - L\,J(T)\right]$",
        "heating"),
       (r"$\dfrac{d\ln m}{dt} = -\dfrac{6}{\rho_p d}\;J(T)$", "evaporation")]
for k, (eq, name) in enumerate(eqs):
    yy = 0.735 - 0.09 * k
    L.text(0.04, yy, eq, fontsize=15, va="center")
    L.text(1.0, yy, name, fontsize=12, color=INK2, va="center", ha="right")
L.text(0, 0.505, "3   Use speed as the clock: divide by dV/dt", fontsize=15, weight="bold", color=ACCENT, va="top")
L.add_patch(plt.Rectangle((0.02, 0.335), 0.96, 0.115, facecolor=BOX, edgecolor="none"))
L.text(0.05, 0.3925, r"$\dfrac{dT}{dV} = -\dfrac{q_{in} - q_{rad} - L J}{c\;\tau}$", fontsize=16, va="center")
L.text(0.58, 0.3925, r"$\dfrac{d\ln m}{dV} = \dfrac{J}{\tau}$", fontsize=16, va="center")
L.text(0.04, 0.305, r"No $d$ is left, and no $\rho_p$.  The same start ($V_0$, $T_0$) gives the same" "\n"
       r"$T(V)$ and the same $m/m_0$ for every size.", fontsize=13.5, va="top", linespacing=1.5)
L.text(0.04, 0.198, r"Size only stretches the clock: braking time $\propto \rho_p\, d$.",
       fontsize=13.5, va="top", color=INK2)
L.text(0, 0.135, "Holds while", fontsize=13, weight="bold", color=INK2, va="top")
L.text(0.04, 0.09, "free-molecule flow (Kn ≫ 1): drag and heat coefficients do not depend on d;\n"
       "one uniform droplet temperature;  braking at almost constant altitude.",
       fontsize=11.5, color=INK2, va="top", linespacing=1.4)
L.text(0.0, -0.005, r"$\tau$ drag per unit surface area,  $q$ heat fluxes,  $J$ evaporation mass flux,  "
       r"$L$ latent heat,  $c$ specific heat,  $\rho_a$ air density", fontsize=10.5, color=MUTED, va="top")

# ---- right: the evidence
R = fig.add_axes([0.6, 0.2, 0.3, 0.57])
for h, col in zip(ALTITUDES, ORDINAL):
    s = f[f.h0_km == h].sort_values("d0_um")
    R.plot(s.d0_um, 100 * s.mass_lost_frac, color=col, marker="o", markersize=5,
           markeredgecolor=SURFACE, markeredgewidth=1.0)
    R.text(s.d0_um.iloc[-1] * 1.25, 100 * s.mass_lost_frac.iloc[-1], "{} km".format(h),
           color=INK, fontsize=12, va="center")
R.set_xscale("log")
R.set_xticks([1, 10, 100])
R.set_xticklabels(["1", "10", "100"])
R.set_ylim(0, 50)
R.set_xlabel(r"initial droplet diameter  [$\mu$m]")
R.set_ylabel("mass lost to vapour  [%]")
R.set_title("Mass lost, 1 to 200 µm: flat", fontsize=14, weight="bold", loc="left", pad=12)
R.text(1.03, 1.0, "release\naltitude", transform=R.transAxes, fontsize=11, color=INK2, ha="left", va="bottom")
s62 = f[(f.h0_km == 62)].set_index("d0_um").mass_lost_frac
R.annotate("Kn ≈ 3: leaving free-molecule\nflow, heating drops", xy=(200, 100 * s62.iloc[-1]),
           xytext=(9, 10), fontsize=11, color=INK2,
           arrowprops=dict(arrowstyle="-", color=MUTED, lw=1.0))
fig.text(0.6, 0.06, "Droplet demise model, AA7075 released along the parent's DRAMA trajectory\n"
         "output/droplet_reentry_output/demise_sweep/final_states.csv", fontsize=9.5, color=MUTED)

for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(FIG, "11_size_independence_slide." + ext), dpi=200)
print(os.path.join(FIG, "11_size_independence_slide.png"))
