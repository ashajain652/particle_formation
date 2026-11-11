"""Video of one wake-column run: the plume drawn as a column growing over time, the processes spreading it, the gas
concentrations, the reactions dominating at each moment, the nanoparticles as they form, and gravitational settling.
An opening (plume_prologue.py) first zooms in from the parent's re-entry to the droplet and the atoms leaving it.

    "$PY" -m reentry_model.wake_column.analysis.plume_video --diameter 200 --altitude 70 --time-of-day day
    "$PY" -m reentry_model.wake_column.analysis.plume_video ... --frames 0,600 --opening-frames 100    # sample PNGs only

Everything shown comes from the model's own state and rate terms: `compute` runs the chain of `run.run` (trail, Stage A,
Box) on a fine time grid and asks the Box for its event rates, so nothing is reimplemented. What is schematic is the
drawing: the radius axis is logarithmic, the dots show shares and not counts, and the falling particles are exaggerated
(the true fall distance is printed beside them). Writes to OUTPUT/video/<run name>/.
"""
import argparse
import math
import multiprocessing
import os
import pickle
import shutil
import subprocess

import imageio_ffmpeg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
import numpy as np
from scipy.integrate import cumulative_trapezoid

from reentry_model.wake_column import OUTPUT
from reentry_model.wake_column.analysis import plume_prologue
from reentry_model.wake_column.ambient import AmbientColumn
from reentry_model.wake_column.box import Box
from reentry_model.wake_column.particles import formation_pathways
from reentry_model.wake_column.reactions import load_reactions
from reentry_model.wake_column.settings import Case, Settings, run_name
from reentry_model.wake_column.slowing import calibrate, slow_down
from reentry_model.wake_column.species import AL_SPECIES, load_species
from reentry_model.wake_column.trail import trail_for
from reentry_model.wake_column.transport import D_molecular

FPS = 30
SECONDS = 105                     # of time-lapse, logarithmic in model time; holds at both ends are added on top
HOLD_START_S, HOLD_END_S = 3.0, 5.0
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8b8a85", "#e4e3df"
SPECIES_COLOUR = dict(zip(AL_SPECIES, ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6")))
NAME = {"Al": "Al", "AlO": "AlO", "OAlO": "OAlO", "OAlO2": "OAlO₂", "AlCO3": "AlCO₃", "AlOH": "AlOH", "AlOH2": "Al(OH)₂",
        "O2": "O₂", "O3": "O₃", "H2": "H₂", "H2O": "H₂O", "CO2": "CO₂"}
TYPES = (("Al2O3", "Al₂O₃", "no hydrogen", "#2a78d6"), ("Orich", "Al₂O₄ to Al₂O₆", "oxygen-rich", "#eb6834"),
         ("AlOH2pair", "Al₂O₂(OH)₄", "from Al(OH)₂ + Al(OH)₂", "#1baf7a"),
         ("Hother", "other with H", "e.g. Al₂O₃H₃, Al₂O₅H₂", "#eda100"))
SMOKE_COLOUR, EDDY_COLOUR, SETTLE_COLOUR = "#8a7a66", "#0e8f8f", "#7a4a1e"
N_DOTS, N_GLYPHS, N_SMOKE = 460, 150, 26
RING_FRAMES = 50                       # a new nanoparticle or newly coated smoke particle keeps its fading ring this long


# ---- the data --------------------------------------------------------------------------------------------------------

def label(name):
    return NAME.get(name, name)


def reaction_label(rx):
    def side(items):
        return " + ".join(label(p) for p in items)
    return f"{side(rx.reactants)} → {side(rx.products)}" if rx.form != "photolysis" else \
        f"{side(rx.reactants)} + sunlight → {side(rx.products)}"


def pair_type(species, a, b):
    """The type of the particle formed by molecules a and b (oxygen is topped up to Al2O3 by the particle's uptake)."""
    o = species[a].deliver_O + species[b].deliver_O
    h = species[a].deliver_H + species[b].deliver_H
    if h == 0:
        return "Al2O3" if o <= 3 else "Orich"
    return "AlOH2pair" if a == b == "AlOH2" else "Hother"


def compute(case, settings):
    """The run and every per-time quantity the frames need."""
    species = load_species()
    trail = trail_for(case.d0_um, case.h0_km, settings.n_emission_samples)
    column = AmbientColumn.load(settings.ambient, case.time_of_day)
    amb = column.at(trail.z_mean_km)
    hot = slow_down(trail.samples, column, settings, species, calibrate(load_reactions(stage=("hot", "both")), species))
    box = Box(trail, hot, amb, case.time_of_day, settings, species=species)
    sol = box.solve()
    lay = box.layout

    keep = np.where(sol.t > 0.0)[0]
    t = sol.t[keep]
    paths = formation_pathways()
    labels = [f"{reaction_label(rx)}" for rx in box.reactions]
    labels += [f"{label(p.a)} + {label(p.b)} → new particle" for p in paths]
    labels += [f"{label(a)} sticks to smoke" for a in AL_SPECIES]
    n_col = len(labels)
    labels += [f"{label(a)} + particle → growth" for a in AL_SPECIES]
    labels += ["O, O₂, O₃ + particle → oxygen uptake", "particle + particle → coagulation"]
    rates = np.zeros((len(t), len(labels)))
    out = {k: np.zeros(len(t)) for k in ("r", "V", "dr2", "dlnV", "N", "O_Al", "H_Al", "d_nm", "lam", "smoke_uptake_ref")}
    conc = {a: np.zeros(len(t)) for a in AL_SPECIES}
    frac = {a: np.zeros(len(t)) for a in tuple(AL_SPECIES) + ("particles", "smoke", "settled")}
    for j, i in enumerate(keep):
        x = sol.Y[i]
        gas, r2 = x[lay.gas], float(x[lay.r2][0])
        V, dV, dr2 = box._geometry(r2)
        r = math.sqrt(r2)
        rates[j, :n_col] = box._network_rates(gas, V)
        p = box._pstate(x, V, r)
        N, al = p["N"], p["N"] * box.bins.n_atoms
        rates[j, n_col:n_col + len(AL_SPECIES)] = p["rho_ks"] @ N
        rates[j, n_col + len(AL_SPECIES)] = float((N * p["Z"] * p["gval"]).sum())
        rates[j, n_col + len(AL_SPECIES) + 1] = 0.5 * float((N * p["KN"]).sum())
        b = sol.bins(i)
        w = al / al.sum() if al.sum() > 0 else al
        out["r"][j], out["V"][j], out["dr2"][j], out["dlnV"][j] = r, V, dr2, dV / V
        out["N"][j] = N.sum()
        out["O_Al"][j], out["H_Al"][j] = float(w @ b["O_per_Al"]), float(w @ b["H_per_Al"])
        out["d_nm"][j] = float(w @ b["d_m"]) * 1e9
        out["lam"][j] = float(w @ p["lam"])
        bud = sol.al_budget(i)
        for a in frac:
            frac[a][j] = bud[a]
        for a in AL_SPECIES:
            conc[a][j] = sol.concentration(a, i)

    # particles formed so far, by type, from the box's own formation rates
    n_rx = len(box.reactions)
    formed = {k: np.zeros(len(t)) for k, *_ in TYPES}
    form_rate = {k: np.zeros(len(t)) for k, *_ in TYPES}
    for q, p in enumerate(paths):
        kind = pair_type(species, p.a, p.b)
        form_rate[kind] += rates[:, n_rx + q]
    for k in formed:
        formed[k] = cumulative_trapezoid(form_rate[k], t, initial=0.0) + form_rate[k][0] * t[0] / 2
    fall = cumulative_trapezoid(out["lam"] * math.pi * out["r"] / 2.0, t, initial=0.0)   # lam = 2 v_s / (pi r)
    return dict(case=case, t=t, labels=labels, rates=rates, conc=conc, frac=frac, formed=formed, form_rate=form_rate,
                fall=fall, D_mol=D_molecular(species["AlO"], amb), K_eddy=box.K_eddy, L=box.L, r0=box.r0,
                N_Al=box.N_Al, z_km=amb.z_km, T_K=amb.T_K, hot=dict(hot.fractions),
                m_Al=trail.m_Al_kg, **out)


def frame_times(d):
    n = int(SECONDS * FPS)
    return np.logspace(np.log10(d["t"][0]), np.log10(d["t"][-1]), n)


def lerp(d, frame_t, y, log=True):
    """y at time frame_t, log-log between the outputs (zeros stay near zero)."""
    lt = np.log10(d["t"])
    if log:
        return 10 ** np.interp(np.log10(frame_t), lt, np.log10(np.maximum(y, 1e-300)))
    return np.interp(np.log10(frame_t), lt, y)


# ---- the frames ------------------------------------------------------------------------------------------------------

def fmt_time(t):
    if t < 1e-3:
        return f"{t * 1e6:.3g} µs"
    if t < 1:
        return f"{t * 1e3:.3g} ms"
    if t < 120:
        return f"{t:.3g} s"
    if t < 7200:
        return f"{t / 60:.3g} min"
    return f"{t / 3600:.3g} h"


def fmt_len(x):
    if x < 1e-6:
        return f"{x * 1e9:.3g} nm"
    if x < 1e-3:
        return f"{x * 1e6:.3g} µm"
    if x < 1:
        return f"{x * 1e3:.3g} mm"
    if x < 1e3:
        return f"{x:.3g} m"
    return f"{x / 1e3:.3g} km"


def count(x):
    return "under 1" if x < 1 else sci(x)


def sub(x):
    return f"{x:.2f}".translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))


def sci(x, unit=""):
    if x <= 0:
        return f"0{unit}"
    e = int(math.floor(math.log10(x)))
    if -2 <= e <= 2:
        return f"{x:.3g}{unit}"
    sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{x / 10 ** e:.1f}×10{sup}{unit}"


class Static:
    """Random but fixed positions, so a dot keeps its place in the column from frame to frame."""

    def __init__(self, d):
        rng = np.random.default_rng(20261105)
        self.u, self.v, self.q = rng.uniform(0.075, 0.925, N_DOTS), rng.uniform(-0.96, 0.96, N_DOTS), rng.random(N_DOTS)
        total = sum(d["formed"][k][-1] for k, *_ in TYPES)
        self.per_glyph = total / N_GLYPHS
        self.glyphs = []                       # (x, v, type index, appearance time)
        for ti, (k, *_ ) in enumerate(TYPES):
            final = d["formed"][k][-1]
            n = max(int(round(final / self.per_glyph)), 3)
            levels = (np.arange(n) + 0.5) / n * final
            tt = 10 ** np.interp(levels, d["formed"][k], np.log10(d["t"]))
            for a in tt:
                self.glyphs.append((rng.uniform(0.08, 0.92), rng.uniform(-0.9, 0.9), ti, a))
        self.smoke = np.column_stack([rng.uniform(0.05, 0.95, N_SMOKE), rng.uniform(-1.12, 1.12, N_SMOKE)])
        self.h_min = math.log10(d["r0"] * 0.35)
        self.h_max = math.log10(d["r"][-1] * 1.12)
        # each smoke particle inside the column picks up aluminium once: when the plume's uptake passes an evenly spaced
        # share of its final value (from the uptake that first shows, 1e-6) and the column has reached it
        uptake = np.maximum.accumulate(d["frac"]["smoke"])
        lt = np.log10(d["t"])
        heights = (np.log10(d["r"]) - self.h_min) / (self.h_max - self.h_min)
        taken = np.argsort(np.argsort(rng.random(N_SMOKE)))
        levels = np.geomspace(1e-6, max(uptake[-1], 1e-6), N_SMOKE)[taken]
        arrive = np.interp(np.abs(self.smoke[:, 1]), np.maximum.accumulate(heights), lt, right=np.inf)
        self.smoke_t = 10 ** np.maximum(np.interp(np.log10(levels), np.log10(np.maximum(uptake, 1e-300)), lt, right=np.inf), arrive)

    def height(self, r):
        return (math.log10(r) - self.h_min) / (self.h_max - self.h_min)


def draw_column(ax, d, st, t, dt_flash, f):
    r = lerp(d, t, d["r"])
    h = st.height(r)
    gas_c = sum(lerp(d, t, d["conc"][a]) for a in AL_SPECIES)
    alpha = 0.10 + 0.34 * np.clip((math.log10(gas_c) - 6) / 10, 0, 1)
    ax.set_xlim(-0.075, 1.0); ax.set_ylim(-1.3, 1.26); ax.axis("off")
    # the log radius scale, on both sides of the axis
    for dec in range(-4, 6):
        x = 10.0 ** dec
        if d["r0"] * 0.35 < x < d["r"][-1] * 1.12:
            y = st.height(x)
            lab = fmt_len(x)
            for s in (1, -1):
                ax.plot([0.012, 0.024], [s * y, s * y], color=MUTED, lw=1)
                ax.plot([0.024, 0.985], [s * y, s * y], color=GRID, lw=0.6, zorder=0)
            ax.text(0.008, y, lab, ha="right", va="center", fontsize=10, color=INK2)
            ax.text(0.008, -y, lab, ha="right", va="center", fontsize=10, color=INK2)
    ax.text(-0.068, 0, "column radius (logarithmic scale)", rotation=90, ha="center", va="center", fontsize=11, color=INK2)
    # the column: a body and two end caps
    body = "#7f9fc9"
    ax.add_patch(Rectangle((0.06, -h), 0.88, 2 * h, fc=body, ec=INK2, lw=1.2, alpha=alpha, zorder=2))
    capw = 0.012 + 0.05 * h
    ax.add_patch(Ellipse((0.94, 0), capw, 2 * h, fc=body, ec=INK2, lw=1.2, alpha=min(alpha * 1.5, 1), zorder=2))
    ax.add_patch(Ellipse((0.06, 0), capw, 2 * h, fc=SURFACE, ec=INK2, lw=1.2, zorder=2))
    # background smoke: fixed in the air; the plume sweeps over it and aluminium collects on it
    ax.scatter(st.smoke[:, 0], st.smoke[:, 1], s=46, c=SMOKE_COLOUR, alpha=0.75, zorder=3, linewidths=0)
    with np.errstate(divide="ignore"):
        age = np.log10(t / st.smoke_t)
    for (x, y), a in zip(st.smoke[(age >= 0) & (age < dt_flash)], age[(age >= 0) & (age < dt_flash)]):
        grow = 1 - a / dt_flash
        ax.scatter([x], [y], s=46 + 170 * grow, facecolors="none", edgecolors=SPECIES_COLOUR["AlOH2"],
                   linewidths=1.8, alpha=grow, zorder=3)
    # gas aluminium: each dot one fixed share of the aluminium, coloured by the molecule that share is in now
    shares = np.array([lerp(d, t, d["frac"][a]) for a in AL_SPECIES])
    edges = np.cumsum(shares)
    idx = np.searchsorted(edges, st.q)
    live = idx < len(AL_SPECIES)
    cols = np.array([SPECIES_COLOUR[AL_SPECIES[i]] for i in idx[live]])
    ax.scatter(st.u[live], st.v[live] * h, s=20, c=cols, linewidths=0, alpha=0.95, zorder=4)
    # new nanoparticles: one glyph per st.per_glyph particles, flashing as it forms
    for x, v, ti, ta in st.glyphs:
        if ta <= t:
            c = TYPES[ti][3]
            ax.scatter([x], [v * h], s=70, marker="h", c=c, edgecolors=INK, linewidths=0.6, zorder=6)
            age = math.log10(t / ta)
            if age < dt_flash:
                grow = 1 - age / dt_flash
                ax.scatter([x], [v * h], s=70 + 520 * grow, marker="o", facecolors="none", edgecolors=c,
                           linewidths=2.2 * grow + 0.3, alpha=grow, zorder=7)
    # what is spreading the column
    total = d["D_mol"] + eddy_now(d, t)
    share_mol, share_eddy = d["D_mol"] / total, eddy_now(d, t) / total
    for x in np.linspace(0.14, 0.86, 7):
        for s in (1, -1):
            if share_mol > 0.02:
                ax.annotate("", xy=(x, s * (h + 0.07 * share_mol + 0.02)), xytext=(x, s * max(h - 0.02, 0)),
                            arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.7, alpha=0.3 + 0.7 * share_mol), zorder=8)
            ax.annotate("", xy=(x + 0.04, s * (h + 0.015)), xytext=(x + 0.04, s * (h + 0.095)),
                        arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2, alpha=0.8), zorder=8)
    if share_eddy > 0.01:
        for x in np.linspace(0.1, 0.9, 6):
            for s in (1, -1):
                a = np.linspace(0, 4 * np.pi, 70) + f * 0.22 * s
                rad = np.linspace(0.002, 0.024, 70)
                ax.plot(x + rad * np.cos(a), s * h + 1.75 * rad * np.sin(a), color=EDDY_COLOUR, lw=1.8,
                        alpha=0.25 + 0.75 * share_eddy, zorder=8)
    # settling: particles leave through the lower side (exaggerated; the true distance is in the settling panel)
    if lerp(d, t, d["N"]) * d["lam"][-1] > 0 and sum(lerp(d, t, d["formed"][k]) for k, *_ in TYPES) > 0:
        for i, x in enumerate(np.linspace(0.12, 0.88, 9)):
            phase = ((f * 0.012 + i * 0.37) % 1.0)
            ax.scatter([x], [-h - 0.015 - phase * 0.13], s=26, marker="h", c=SETTLE_COLOUR, alpha=(1 - phase) * 0.9,
                       linewidths=0, zorder=8)
    ax.text(0.5, 1.2, f"Now {sci(sum(lerp(d, t, d['conc'][a]) for a in AL_SPECIES) * 1e-6)} Al molecules per cm³ of plume",
            ha="center", va="center", fontsize=11, color=INK2)
    # the column's length, as a ruler; the box's volume adds the spreading along the axis once r nears L
    y = -1.235
    ax.annotate("", xy=(0.94, y), xytext=(0.06, y), arrowprops=dict(arrowstyle="|-|", color=INK, lw=1.6, mutation_scale=7))
    effective = lerp(d, t, d["V"]) / (math.pi * r * r)
    text = f"column length {fmt_len(d['L'])}"
    if effective > 1.01 * d["L"]:
        text += f" ({fmt_len(effective)} with spreading along the axis)"
    ax.text(0.5, y, f"  {text}, not drawn to the radius scale  ", ha="center", va="center", fontsize=11.5, color=INK,
            bbox=dict(fc=SURFACE, ec="none", pad=1.5), zorder=9)


def eddy_now(d, t):
    return max(lerp(d, t, d["dr2"], log=False) / 4.0 - d["D_mol"], 0.0)


def draw_key(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    items = [("arrow", INK, "molecular diffusion"), ("swirl", EDDY_COLOUR, "eddy diffusion"),
             ("arrow", MUTED, "air entering"), ("glyph", SETTLE_COLOUR, "settling (exaggerated)"),
             ("smoke", SMOKE_COLOUR, "meteoric smoke"), ("dot", SPECIES_COLOUR["AlO"], "gas aluminium")]
    for x, (kind, c, text) in zip((0.01, 0.175, 0.315, 0.445, 0.665, 0.825), items):
        if kind == "arrow":
            ax.annotate("", xy=(x + 0.025, 0.5), xytext=(x, 0.5), arrowprops=dict(arrowstyle="-|>", color=c, lw=2))
        elif kind == "swirl":
            a = np.linspace(0, 4 * np.pi, 60); rad = np.linspace(0.001, 0.011, 60)
            ax.plot(x + 0.012 + rad * np.cos(a), 0.5 + 16 * rad * np.sin(a), color=c, lw=1.6)
        else:
            ax.scatter([x + 0.012], [0.5], s=20 if kind == "dot" else 46, c=c, marker="h" if kind == "glyph" else "o",
                       linewidths=0)
        ax.text(x + 0.032, 0.5, text, fontsize=10, va="center", color=INK2)


def draw_timeline(ax, d, t):
    ax.set_xscale("log"); ax.set_xlim(d["t"][0], d["t"][-1]); ax.set_ylim(0, 1)
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_yticks([])
    marks = ((1e-6, "1 µs"), (1e-3, "1 ms"), (1, "1 s"), (60, "1 min"), (3600, "1 h"), (d["t"][-1], "12 h"))
    ax.set_xticks([m[0] for m in marks]); ax.set_xticklabels([m[1] for m in marks], fontsize=11, color=INK2)
    ax.minorticks_off()
    ax.axvspan(d["t"][0], t, color="#7f9fc9", alpha=0.55)
    ax.axvline(t, color=INK, lw=2)


def draw_concentrations(ax, d, t):
    names = list(AL_SPECIES)
    rows = [(label(a), lerp(d, t, d["conc"][a]) * 1e-6, lerp(d, t, d["frac"][a]), SPECIES_COLOUR[a]) for a in names]
    V = lerp(d, t, d["V"])
    n = sum(lerp(d, t, d["formed"][k]) for k, *_ in TYPES)
    rows.append(("Nanoparticles", lerp(d, t, d["N"]) / V * 1e-6, lerp(d, t, d["frac"]["particles"]), INK))
    rows.append(("Al on smoke", lerp(d, t, d["frac"]["smoke"]) * d["N_Al"] / V * 1e-6, lerp(d, t, d["frac"]["smoke"]),
                 SMOKE_COLOUR))
    rows = rows[::-1]
    y = np.arange(len(rows))
    lo, hi = 1e-4, 3e11
    for i, (name, c, fr, col) in enumerate(rows):
        ax.barh(i, max(c, lo) - lo, left=lo, height=0.68, color=col, linewidth=0)
        ax.text(hi * 1.3, i, f"{sci(c)}", va="center", ha="left", fontsize=10.5, color=INK)
        ax.text(hi * 1.3e4, i, f"{fr * 100:.3g} %" if fr >= 1e-4 else "under 0.01 %", va="center", ha="left",
                fontsize=10.5, color=INK2)
    ax.set_xscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=11, color=INK)
    ax.set_xticks([1e-3, 1e0, 1e3, 1e6, 1e9])
    ax.tick_params(axis="x", labelsize=9.5, colors=INK2); ax.tick_params(axis="y", length=0)
    ax.minorticks_off(); ax.grid(True, axis="x", color=GRID, lw=0.7); ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title("Concentrations in the plume now", loc="left", fontsize=13, fontweight="bold", color=INK, pad=8)
    ax.text(hi * 1.3, len(rows) - 0.25, "per cm³", fontsize=9.5, color=MUTED, ha="left", va="bottom")
    ax.text(hi * 1.3e4, len(rows) - 0.25, "share of Al", fontsize=9.5, color=MUTED, ha="left", va="bottom")


def draw_reactions(ax, d, t, top=6):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    lt, tt = np.log10(t), np.log10(d["t"])
    i = min(int(np.searchsorted(tt, lt)), len(tt) - 1)
    i0 = max(i - 1, 0)
    w = 0.0 if i == i0 else (lt - tt[i0]) / (tt[i] - tt[i0])
    rate = 10 ** ((1 - w) * np.log10(np.maximum(d["rates"][i0], 1e-300)) + w * np.log10(np.maximum(d["rates"][i], 1e-300)))
    order = np.argsort(rate)[::-1][:top]
    total = rate.sum()
    ax.text(0, 1.0, "Dominant reactions now", fontsize=13, fontweight="bold", color=INK, va="top")
    ax.text(1, 0.985, "events per second in the whole plume", fontsize=9.5, color=MUTED, va="top", ha="right")
    for k, j in enumerate(order):
        y = 0.83 - k * 0.15
        share = rate[j] / total
        ax.add_patch(Rectangle((0, y - 0.058), share, 0.116, fc="#dce6f2", ec="none", zorder=0))
        ax.text(0.008, y, f"{k + 1}", fontsize=11, color=MUTED, va="center")
        ax.text(0.04, y, d["labels"][j], fontsize=12, color=INK, va="center", fontweight="bold" if k == 0 else "normal")
        ax.text(0.83, y, sci(rate[j]), fontsize=11, color=INK2, va="center", ha="right")
        ax.text(1.0, y, f"{share * 100:.3g} %", fontsize=11, color=INK2, va="center", ha="right")


def draw_processes(ax, d, st, t):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    eddy = eddy_now(d, t)
    total = d["D_mol"] + eddy
    dlnV = lerp(d, t, d["dlnV"])
    lam = lerp(d, t, np.maximum(d["lam"], 1e-30))
    formed_any = sum(lerp(d, t, d["formed"][k]) for k, *_ in TYPES) > 1
    r = lerp(d, t, d["r"])
    fall = lerp(d, t, np.maximum(d["fall"], 1e-30)) if formed_any else 0.0
    rows = [("Molecular diffusion", d["D_mol"] / total, f"{d['D_mol']:.3g} m²/s, {d['D_mol'] / total:.0%} of the spreading"),
            ("Eddy diffusion", eddy / total, f"{eddy:.3g} m²/s, {eddy / total:.0%} of the spreading"),
            ("Air entering the plume", 1.0, f"volume doubles in {fmt_time(math.log(2) / dlnV)}"),
            ("Particle settling", 1.0 if formed_any else 0.0,
             f"loses particles at {sci(lam)} per second" if formed_any else "no particles yet")]
    ax.text(0, 1.0, "Processes active now", fontsize=13, fontweight="bold", color=INK, va="top")
    for k, (name, level, text) in enumerate(rows):
        y = 0.79 - k * 0.18
        ax.scatter([0.015], [y], s=150, c=[(0.11, 0.66, 0.48, 0.12 + 0.88 * level)], edgecolors=MUTED, linewidths=0.7)
        ax.text(0.045, y, name, fontsize=12, color=INK, va="center")
        ax.text(0.40, y, text, fontsize=11, color=INK2, va="center")
    y = 0.05
    ax.text(0.0, y, f"Settling to true scale: a particle has fallen {fmt_len(fall) if fall >= 1e-9 else 'under 1 nm'}; "
            f"the column radius is {fmt_len(r)}; {sci(lerp(d, t, d['frac']['settled']) * 100)} % of the Al has settled out.",
            fontsize=11, color=SETTLE_COLOUR, va="center")


def draw_particles(ax, d, st, t):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    N = lerp(d, t, d["N"])
    ax.text(0, 1.0, "Nanoparticles forming now", fontsize=13, fontweight="bold", color=INK, va="top")
    ax.text(1, 0.985, f"each hexagon in the column is {sci(st.per_glyph)} particles", fontsize=9.5, color=MUTED, va="top", ha="right")
    ax.text(0.47, 0.82, "formed so far", fontsize=9.5, color=MUTED, ha="right", va="center")
    ax.text(0.72, 0.82, "per second now", fontsize=9.5, color=MUTED, ha="right", va="center")
    total_rate = sum(lerp(d, t, d["form_rate"][k]) for k, *_ in TYPES)
    for k, (key, name, note, col) in enumerate(TYPES):
        y = 0.69 - k * 0.135
        rate = lerp(d, t, d["form_rate"][key])
        ax.scatter([0.015], [y], s=120, marker="h", c=col, edgecolors=INK, linewidths=0.6)
        ax.text(0.045, y, name, fontsize=12, color=INK, va="center")
        ax.text(0.47, y, count(lerp(d, t, d["formed"][key])), fontsize=11, color=INK2, va="center", ha="right")
        ax.text(0.72, y, count(rate), fontsize=11, color=INK2, va="center", ha="right")
        ax.text(0.74, y, f"{rate / total_rate:.0%}" if total_rate > 0 else "", fontsize=10.5, color=MUTED, va="center")
    O, H = lerp(d, t, d["O_Al"], log=False), lerp(d, t, d["H_Al"], log=False)
    dn = lerp(d, t, d["d_nm"], log=False)
    ax.text(0, 0.07, f"In the plume: {sci(N)} particles, about {dn:.2f} nm across, average Al₂O{sub(2 * O)}H{sub(2 * H)}",
            fontsize=11, color=INK, va="center")


_cache = {}


def render(args):
    data_path, i, out = args
    if data_path not in _cache:
        with open(data_path, "rb") as f:
            d = pickle.load(f)
        _cache[data_path] = d, Static(d), frame_times(d)
    d, st, times = _cache[data_path]
    t = times[i]
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": INK})
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor=SURFACE)
    c = d["case"]
    fig.text(0.012, 0.965, f"Aluminium plume from a {c.d0_um:g} µm droplet released at {c.h0_km:g} km, "
             f"{'daytime' if c.time_of_day == 'day' else 'night-time'} air", fontsize=19, fontweight="bold", va="center")
    fig.text(0.988, 0.965, f"t = {fmt_time(t)}", fontsize=26, fontweight="bold", va="center", ha="right")
    fig.text(0.012, 0.930, f"radius {fmt_len(lerp(d, t, d['r']))}     volume {sci(lerp(d, t, d['V']), ' m³')}     "
             f"{sci(d['N_Al'])} aluminium atoms vaporised ({sci(d['m_Al'] * 1e6)} mg)     air {d['T_K']:.0f} K     "
             f"model output on a schematic drawing", fontsize=12, color=INK2, va="center")
    draw_timeline(fig.add_axes([0.012, 0.885, 0.976, 0.022]), d, t)
    ring_decades = RING_FRAMES * math.log10(times[-1] / times[0]) / (len(times) - 1)     # frames are evenly spaced in log time
    draw_column(fig.add_axes([0.012, 0.335, 0.60, 0.515]), d, st, t, ring_decades, i)
    draw_key(fig.add_axes([0.03, 0.305, 0.59, 0.03]))
    draw_concentrations(fig.add_axes([0.70, 0.525, 0.175, 0.31]), d, t)
    draw_reactions(fig.add_axes([0.625, 0.235, 0.363, 0.25]), d, t)
    draw_processes(fig.add_axes([0.012, 0.022, 0.60, 0.27]), d, st, t)
    draw_particles(fig.add_axes([0.625, 0.022, 0.363, 0.205]), d, st, t)
    if i == 0:
        h = d["hot"]
        fig.text(0.32, 0.755, f"Start: the slowing-down stage is over. The vapour is {h['Al']:.0%} Al, {h['AlO']:.0%} AlO, "
                 f"{h['OAlO']:.1%} OAlO\nin a column of radius {fmt_len(d['r0'])}", ha="center", va="center", fontsize=14,
                 bbox=dict(boxstyle="round,pad=0.5", fc=SURFACE, ec=MUTED))
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    return i


def render_opening(args):
    data_path, i, out = args
    if data_path not in _cache:
        with open(data_path, "rb") as f:
            d = pickle.load(f)
        _cache[data_path] = d, plume_prologue.Fixed(d)
    d, fx = _cache[data_path]
    return plume_prologue.frame(d, fx, i, out)


def cached(path, make):
    if not os.path.exists(path):
        data = make()
        with open(path, "wb") as f:
            pickle.dump(data, f)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--diameter", type=float, default=200.0)
    ap.add_argument("--altitude", type=float, default=70.0)
    ap.add_argument("--time-of-day", default="day", choices=("day", "night"))
    ap.add_argument("--points-per-decade", type=int, default=40)
    ap.add_argument("--workers", type=int, default=max((os.cpu_count() or 2) - 1, 1))
    ap.add_argument("--frames", default=None, help="comma-separated frame numbers: write only these PNGs")
    ap.add_argument("--opening-frames", default=None, help="the same for the opening scenes")
    ap.add_argument("--out", default=os.path.join(OUTPUT, "video"))
    a = ap.parse_args(argv)
    case = Case(a.diameter, a.altitude, a.time_of_day)
    settings = Settings(points_per_decade=a.points_per_decade)
    folder = os.path.join(a.out, run_name(case, Settings()))
    frames_dir = os.path.join(folder, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    data_path = cached(os.path.join(folder, "data.pkl"), lambda: compute(case, settings))
    with open(data_path, "rb") as f:
        main_data = pickle.load(f)
    n = len(frame_times(main_data))
    opening_path = cached(os.path.join(folder, "opening.pkl"), lambda: plume_prologue.compute(case, main_data))
    def job(prefix, path, i):
        return path, i, os.path.join(frames_dir, f"{prefix}{i:05d}.png")
    picked = a.frames is not None or a.opening_frames is not None
    main_frames = [int(x) for x in (a.frames or "").split(",") if x] if picked else list(range(n))
    opening_frames = ([int(x) for x in (a.opening_frames or "").split(",") if x] if picked
                      else list(range(plume_prologue.N_FRAMES)))
    with multiprocessing.get_context("spawn").Pool(a.workers) as pool:
        for fn, jobs in ((render_opening, [job("p", opening_path, i) for i in opening_frames]),
                         (render, [job("f", data_path, i) for i in main_frames])):
            for k, _ in enumerate(pool.imap_unordered(fn, jobs, chunksize=4)):
                if (k + 1) % 200 == 0:
                    print(f"{k + 1}/{len(jobs)} frames", flush=True)
    if picked:
        print("wrote frames in", frames_dir)
        return 0
    video = os.path.join(folder, "plume.mp4")
    fade = 1.0
    offset = plume_prologue.N_FRAMES / FPS - fade
    graph = (f"[1:v]tpad=start_mode=clone:start_duration={HOLD_START_S}:stop_mode=clone:stop_duration={HOLD_END_S}[main];"
             f"[0:v]tpad=stop_mode=clone:stop_duration={fade}[open];"
             f"[open][main]xfade=transition=fade:duration={fade}:offset={offset + fade}[v]")
    subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error",
                    "-framerate", str(FPS), "-i", os.path.join(frames_dir, "p%05d.png"),
                    "-framerate", str(FPS), "-i", os.path.join(frames_dir, "f%05d.png"),
                    "-filter_complex", graph, "-map", "[v]", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                    "-preset", "slow", video], check=True)
    shutil.rmtree(frames_dir)
    print("wrote", video)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
