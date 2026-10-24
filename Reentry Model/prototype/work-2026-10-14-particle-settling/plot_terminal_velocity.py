"""Figure for the terminal-velocity sweep: settling speed and 1 km settling time vs diameter."""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

SURFACE = "#fcfcfb"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#dcdbd6"
# sequential blue ramp, ordinal steps 250..700: altitude is an ordered magnitude, so one hue
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
        "#256abf", "#1c5cab", "#104281", "#0d366b"]

rows = list(csv.DictReader(open("terminal_velocity.csv")))
for r in rows:
    for k, v in r.items():
        if k != "regime":
            r[k] = float(v)
alts = sorted({r["h_km"] for r in rows})
assert len(alts) == len(RAMP), (len(alts), len(RAMP))

fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.9), facecolor=SURFACE)
for ax in axes:
    ax.set_facecolor(SURFACE)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.grid(True, which="major", color=GRID, lw=0.7, zorder=0)
    ax.grid(True, which="minor", color=GRID, lw=0.35, alpha=0.6, zorder=0)
    ax.yaxis.set_minor_formatter(NullFormatter())
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9, which="both")
    ax.set_xlabel("sphere diameter  [μm]", color=INK2, fontsize=10)

for c, a in zip(RAMP, alts):
    sub = sorted([r for r in rows if r["h_km"] == a], key=lambda r: r["d_um"])
    d = [r["d_um"] for r in sub]
    axes[0].plot(d, [r["v_t"] for r in sub], color=c, lw=2.0,
                 label="{:.0f} km".format(a), zorder=3)
    axes[1].plot(d, [r["fall_1km_s"] / 86400.0 for r in sub], color=c, lw=2.0,
                 label="{:.0f} km".format(a), zorder=3)

axes[0].set_ylabel("terminal velocity  [m/s]", color=INK2, fontsize=10)
axes[0].set_title("(a)  Gravitational settling velocity", color=INK, fontsize=11,
                  loc="left", pad=10)
axes[1].set_ylabel("time to settle 1 km  [days]", color=INK2, fontsize=10)
axes[1].set_title("(b)  Time to settle 1 km at that velocity", color=INK, fontsize=11,
                  loc="left", pad=10)
axes[1].axhline(1.0 / 24.0, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=2)
axes[1].annotate("1 hour", xy=(95.0, 1.0 / 24.0), xytext=(0, 4), textcoords="offset points",
                 color=INK2, fontsize=8.5, va="bottom", ha="right")
axes[1].axhline(365.25, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=2)
axes[1].annotate("1 year", xy=(95.0, 365.25), xytext=(0, 4), textcoords="offset points",
                 color=INK2, fontsize=8.5, va="bottom", ha="right")

leg = axes[0].legend(title="altitude", frameon=False, fontsize=8.5, ncol=2,
                     loc="lower right", labelcolor=INK2, handlelength=1.6,
                     borderaxespad=0.8)
leg.get_title().set_color(INK2); leg.get_title().set_fontsize(9)

fig.suptitle("Terminal velocity of AA7075 spheres (2813 kg/m³) in US76 air\n"
             "Stokes drag with Cunningham slip correction and the Schiller–Naumann "
             "finite-Reynolds factor",
             color=INK, fontsize=12, x=0.008, ha="left", y=0.985)
fig.tight_layout(rect=(0, 0, 1, 0.945))
fig.savefig("terminal_velocity.png", dpi=200, facecolor=SURFACE)
print("wrote terminal_velocity.png")
