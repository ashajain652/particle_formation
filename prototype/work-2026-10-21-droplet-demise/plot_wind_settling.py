"""Figures for the settling phase under a mean mesospheric vertical wind (run_demise_sweep.py --wind W [--nano]).

  8_time_to_ground_w<W>cms.png   per wind case and size range: time from release to the ground, every 4 km of
                                 release altitude; droplets held aloft are coloured separately and labelled with
                                 the altitude where they hover
  9_wind_effect.png              time to the ground against initial diameter (1 nm - 200 um together) at three
                                 release altitudes, for the three wind cases
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.patches import Patch

import run_demise_sweep as rs

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
SEQ = LinearSegmentedColormap.from_list("blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
HELD = "#eb6834"                        # held aloft (categorical slot 2)
# second sequential ramp, for the altitude where a held droplet hovers (light = low, dark = high)
ORANGE = LinearSegmentedColormap.from_list("orange", ["#fde3d6", "#f8b99a", "#f28d60", "#eb6834", "#c4501f", "#8f3a16"])
HOVER_KM = (15.0, 100.0)                # one fixed hover-altitude scale for every figure, so they can be compared
CASES = [(5.0, 0.3), (0.0, 0.0), (-5.0, -0.3)]          # (mesospheric cm/s, stratospheric mm/s), + up
WIND_COLOURS = {CASES[0]: "#2a78d6", CASES[1]: "#1baf7a", CASES[2]: "#eb6834"}
WIND_LABEL = {CASES[0]: "upward: +5 cm/s mesosphere, +0.3 mm/s stratosphere", CASES[1]: "no vertical wind",
              CASES[2]: "downward: -5 cm/s mesosphere, -0.3 mm/s stratosphere"}
TARGETS = (("t_ground_s", "ground", "the ground"), ("t_tropopause_s", "tropopause", "the tropopause (15 km)"))
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.2, "font.size": 13, "axes.titlesize": 15,
    "axes.titleweight": "bold", "legend.frameon": False, "axes.axisbelow": True})
YEAR = 365.25 * 86400.0


def run_dir(base, case):
    return "{}_w{:+g}cms_s{:+g}mms".format(base, case[0], case[1])


def load(base, w):
    p = os.path.join(run_dir(base, w), "final_states.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def tlabel(t):
    if not np.isfinite(t):
        return ""
    if t < 3600:
        return "{:.2g} s".format(t)
    if t < 86400:
        return "{:.2g} h".format(t / 3600)
    if t < YEAR:
        return "{:.3g} d".format(t / 86400)
    y = t / YEAR
    return "{:,.0f} yr".format(float("{:.3g}".format(y))) if y >= 1000 else "{:.3g} yr".format(y)


# ---------------------------------------------------------------- 8. time to the ground, per wind and size range
for base, unit, scale in ((rs.OUT, "micrometres", 1.0), (rs.OUT_NANO, "nanometres", 1000.0)):
    for w, (tcol, tname, tlong) in [(c, t) for c in CASES for t in TARGETS]:
        f = load(base, w)
        if f is None:
            continue
        D = sorted(f.d0_um.unique())
        H5 = [h for h in sorted(f.h0_km.unique()) if (h - 50) % 4 == 0]
        sub = f[f.h0_km.isin(H5)]
        piv = lambda c: sub.pivot(index="h0_km", columns="d0_um", values=c).reindex(index=H5, columns=D).values
        t = piv(tcol).astype(float)
        held = ~np.isfinite(t)
        hs = piv("h_suspended_km").astype(float)
        in_years = bool(np.isfinite(t).any() and np.nanmax(np.where(np.isfinite(t), t, np.nan)) > YEAR)
        unit_s, unit_name = (YEAR, "years") if in_years else (86400.0, "days")
        tdays = np.ma.masked_where(held | ~np.isfinite(t), t / unit_s)
        cmap = SEQ.copy(); cmap.set_bad(HELD)
        norm = LogNorm(tdays.min(), tdays.max()) if tdays.count() else LogNorm(1e-3, 1)
        fig, a = plt.subplots(figsize=(30, 17))
        im = a.imshow(tdays, origin="lower", aspect="auto", cmap=cmap, norm=norm)
        hheld = np.ma.masked_where(~held | ~np.isfinite(hs), hs)
        hnorm = Normalize(*HOVER_KM) if hheld.count() else None
        im_h = a.imshow(hheld, origin="lower", aspect="auto", cmap=ORANGE, norm=hnorm) if hheld.count() else None
        for j in range(len(D) + 1):
            a.axvline(j - 0.5, color=SURFACE, lw=2)
        for i in range(len(H5) + 1):
            a.axhline(i - 0.5, color=SURFACE, lw=2)
        for i in range(len(H5)):
            for j in range(len(D)):
                if held[i, j]:
                    txt = "held aloft\nat {:.1f} km".format(hs[i, j])
                    col = "#ffffff" if (hnorm is not None and np.isfinite(hs[i, j]) and hnorm(hs[i, j]) > 0.55) else INK
                else:
                    txt, col = tlabel(t[i, j]), ("#ffffff" if norm(tdays[i, j]) > 0.5 else INK)
                a.text(j, i, txt, ha="center", va="center", fontsize=16, linespacing=1.2, color=col)
        a.set_xticks(range(len(D))); a.set_xticklabels(["{:.3g}".format(d * scale) for d in D], fontsize=18)
        a.set_yticks(range(len(H5))); a.set_yticklabels([str(h) for h in H5], fontsize=18)
        a.set_xlabel("initial droplet diameter [{}]".format(unit), fontsize=20, labelpad=10)
        a.set_ylabel("release (starting) altitude [km]", fontsize=20, labelpad=10); a.grid(False)
        a.set_title("Time from release to {}; vertical wind {:+g} cm/s (mesosphere), {:+g} mm/s (stratosphere)".format(
            tlong, w[0], w[1]),
                    loc="left", fontsize=22, pad=14)
        if tdays.count():                  # time scale only where some droplets reach the target
            cb = fig.colorbar(im, ax=a, shrink=0.8, pad=0.012); cb.outline.set_visible(False)
            cb.set_label("time to {} [{}]".format(tlong, unit_name), color=INK2, fontsize=18)
            cb.ax.tick_params(labelsize=16)
        if im_h is not None:               # altitude scale for the droplets held aloft
            cbh = fig.colorbar(im_h, ax=a, shrink=0.8, pad=0.012); cbh.outline.set_visible(False)
            cbh.set_label("altitude where a held droplet hovers [km]", color=INK2, fontsize=18)
            cbh.ax.tick_params(labelsize=16)
        n_held = int((~np.isfinite(f[tcol].astype(float))).sum())
        fig.legend(handles=[Patch(fc=HELD, label="held aloft (orange scale = hover altitude): the updraft exceeds the "
                                                 "terminal velocity ({} of {} cases on the full 2 km grid)".format(n_held, len(f)))],
                   loc="lower left", bbox_to_anchor=(0.01, -0.02), fontsize=16)
        fig.text(0.01, -0.045, "Wind profile: mesospheric value 70-90 km (tapered to 0 at 100 km), linear to the stratospheric "
                 "value at 50 km, stratospheric value 15-50 km, none below 13 km. Bounding case. Units: s, h, d (days), yr.",
                 color=INK2, fontsize=16)
        fig.tight_layout(rect=(0, 0.02, 1, 1))
        out = os.path.join(base, "figures")
        os.makedirs(out, exist_ok=True)
        name = "8_time_to_{}_w{:+g}cms_s{:+g}mms.png".format(tname, w[0], w[1])
        fig.savefig(os.path.join(out, name), dpi=110, bbox_inches="tight")
        plt.close(fig)
        print("wrote", base.split("/")[-1], name)

# ---------------------------------------------------------------- 9. wind effect across all sizes
fig, axs = plt.subplots(2, 3, figsize=(24, 15), sharey=True, sharex=True)
for row, (tcol, tname, tlong) in enumerate(TARGETS):
    for a, h0 in zip(axs[row], (90, 70, 50)):
        for w in CASES:
            parts = [load(b, w) for b in (rs.OUT_NANO, rs.OUT)]
            parts = [p for p in parts if p is not None]
            if not parts:
                continue
            f = pd.concat(parts).sort_values("d0_um")
            g = f[f.h0_km == h0]
            t = g[tcol].astype(float)
            ok = g[np.isfinite(t)]
            a.loglog(ok.d0_um * 1000.0, ok[tcol] / YEAR, color=WIND_COLOURS[w], label=WIND_LABEL[w])
            held = g[~np.isfinite(t)]
            if len(held):
                a.plot(held.d0_um * 1000.0, np.full(len(held), 3e5), "^", color=WIND_COLOURS[w], ms=9)
        a.axhline(3e5, color=MUTED, lw=0.8, ls=":")
        a.text(1.1, 4.5e5, "held aloft (never reaches it)", fontsize=12, color=INK2)
        for yr, lab in ((1 / 365.25, "1 day"), (1.0, "1 year"), (1000.0, "1000 years")):
            a.axhline(yr, color=GRID, lw=1.0)
            a.text(2.5e5, yr * 1.3, lab, fontsize=11, color=INK2, ha="right")
        a.axvline(1000.0, color=MUTED, lw=0.8, ls="--")
        a.set_ylim(1e-5, 2e6)
        a.set_title("To {}, released at {} km".format(tlong, h0), loc="left")
        if row == 1:
            a.set_xlabel("initial droplet diameter [nanometres]")
    axs[row][0].set_ylabel("time from release [years]")
axs[0][0].legend(loc="lower left", fontsize=12)
fig.text(0.01, -0.01, "Nanometre sweep (1-500 nm) and micron sweep (1-200 um) together; the dashed line marks 1 um. "
         "Bounding wind cases; see droplet_demise.MeanVerticalWind for the layer profile.", color=INK2, fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(rs.OUT, "figures", "9_wind_effect.png"), dpi=120, bbox_inches="tight")
print("wrote 9_wind_effect.png")
