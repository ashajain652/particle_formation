"""The reaction network of the wake-column model as one figure, drawn from the tables the model itself reads.

    "$PY" -m reentry_model.wake_column.analysis.plot_reaction_network [--out DIR]
    "$PY" -m reentry_model.wake_column.analysis.plot_reaction_network --flux [--d0 14.1 --h0 76 --tod day] [--out DIR]

Without --flux, writes <out>/plots/reaction_network.png and .pdf, the structure alone. With --flux, runs one case and
writes reaction_flux_<case>.png, .pdf and .csv: arrow thickness is then the aluminium carried by that reaction over the
run, as a share of the inventory, on a logarithmic scale. A share can exceed 1, because an atom going round a cycle
(AlO to OAlO and back, say) is counted each time. The flux run stops if the time integrals fail to close each species'
aluminium budget to BALANCE_TOL. One arrow is one row of data/wake_column/reactions.csv with stage box
or both, so the figure cannot drift from the table: the partners, products and reaction ids come from the row, and the
script stops if a row has no entry in ROUTES (a new reaction must be given a place on the drawing). A row whose
species has no node stops it too.

ROUTES holds drawing choices only (curvature of the arrow and where its label sits). The species placement in NODES
is chosen to keep crossings few; the network cannot be drawn without any, because its six species other than Al are
joined by 13 pairs and a planar drawing of six nodes holds at most 12.
"""
import argparse
import csv
import dataclasses
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402
import numpy as np                                                   # noqa: E402
from matplotlib.lines import Line2D                                  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch       # noqa: E402

from reentry_model.wake_column import OUTPUT                         # noqa: E402
from reentry_model.wake_column.ambient import AmbientColumn          # noqa: E402
from reentry_model.wake_column.box import Box                        # noqa: E402
from reentry_model.wake_column.particles import formation_pathways   # noqa: E402
from reentry_model.wake_column.reactions import load_reactions       # noqa: E402
from reentry_model.wake_column.settings import Case, Settings        # noqa: E402
from reentry_model.wake_column.slowing import calibrate, slow_down   # noqa: E402
from reentry_model.wake_column.species import AL_SPECIES, LABEL, load_species  # noqa: E402
from reentry_model.wake_column.trail import trail_for                # noqa: E402

INK, INK2, HAIR, PANEL = "#1f1f1e", "#5f5e58", "#c9c8c0", "#f6f5f1"
# what the ambient partner of a reaction is made of (dataviz categorical slots 1-4, validated on white)
FAMILY_COLOUR = {"oxygen": "#2a78d6", "carbon": "#eb6834", "hydrogen": "#1baf7a", "light": "#eda100"}
FAMILY_NAME = {"oxygen": "O, O$_2$, O$_3$", "carbon": "CO, CO$_2$", "hydrogen": "H, H$_2$, H$_2$O",
               "light": "sunlight (day only)"}

FLUX_POINTS_PER_DECADE = 40                     # output times for the flux integrals (the default 10 is for plots)
BALANCE_TOL = 1e-3                              # largest allowed budget residual, as a share of the inventory
FLUX_FLOOR, FLUX_CEIL = 1e-4, 300.0             # shares drawn at the thinnest and thickest line
LW_MIN, LW_MAX = 0.5, 7.0
FS_LABEL, FS_TEXT, FS_HEAD, FS_NODE, FS_TITLE = 12, 12, 13.5, 17, 20      # sized for a 16:9 slide
PANEL_X1 = 19.0                                 # right edge of the gas-chemistry panel
HOT_START = ("Al", "AlO", "OAlO")               # slowing.TRACKED: what Stage A hands to Stage B
NODES = {"Al": (1.5, 8.7), "AlO": (5.8, 5.0), "OAlO": (11.4, 8.8), "OAlO2": (17.6, 8.8), "AlCO3": (14.5, 5.0),
         "AlOH": (11.4, 1.2), "AlOH2": (17.6, 1.2)}
NODE_NOTE = {"Al": "atom", "AlO": "monoxide", "OAlO": "dioxide", "OAlO2": "trioxide", "AlCO3": "carbonate",
             "AlOH": "hydroxide", "AlOH2": "dihydroxide"}

# rid: (arc3 curvature, position of the label along the arrow 0-1, label offset x, label offset y)
ROUTES = {
    "R1": (-0.42, 0.5, 0.15, 0.0), "G26": (-0.14, 0.5, 0.0, 0.0),
    "R5": (-0.14, 0.5, 0.0, 0.0), "G33": (-0.42, 0.5, -0.1, 0.0),
    "G42": (-0.42, 0.56, -0.2, 0.2), "G38": (-0.24, 0.56, -0.2, 0.2), "R2": (-0.07, 0.56, -0.2, 0.25),
    "R3": (-0.08, 0.5, 0.25, -0.1), "R4": (-0.26, 0.62, 0.35, -0.1),
    "R6": (-0.52, 0.34, 0.0, 0.0),
    "R7": (0.0, 0.32, 0.0, 0.0),
    "R8": (0.0, 0.5, 0.0, 0.0),
    "R9": (0.0, 0.58, -0.8, 0.0),
    "R10": (-0.27, 0.8, 0.0, 0.0),
    "R11": (0.0, 0.5, 0.0, 0.0),
    "R12": (0.0, 0.5, 0.0, 0.0),
    "R13": (0.0, 0.4, 0.0, 0.0), "R14": (0.0, 0.7, 0.0, 0.0), "R15": (0.0, 0.3, 0.0, 0.0),
    "R16": (0.0, 0.55, 0.0, 0.0),
    "R17": (-0.3, 0.5, -0.45, 0.0), "GH2": (0.0, 0.5, 0.25, 0.0),
}


SOURCES_NOTE = "R: Plane et al. 2021; G: Gobrecht et al. 2022 (GH2: Mangan et al. 2021). M: any air molecule."
STRUCTURE_NOTE = "Label: id: ambient reactant → ambient product. Arrows show which reactions exist, not how fast.\n" \
                 + SOURCES_NOTE
FLUX_NOTE = "Label: id: ambient reactant → ambient product | Al atoms carried per Al atom released (above 1: atoms " \
            "cycle). One case only.\n" + SOURCES_NOTE


def parts(rx):
    """(aluminium-bearing reactant, aluminium-bearing product, ambient reactants, ambient products) of a reaction;
    the third body M is left in the ambient reactants and dropped from the products."""
    al_in = [s for s in rx.reactants if s in AL_SPECIES]
    al_out = [s for s in rx.products if s in AL_SPECIES]
    if len(al_in) != 1 or len(al_out) != 1:
        raise ValueError(f"{rx.rid}: expected one aluminium-bearing reactant and product, got {al_in} and {al_out}")
    return (al_in[0], al_out[0], [s for s in rx.reactants if s not in AL_SPECIES],
            [s for s in rx.products if s not in AL_SPECIES and s != "M"])


def family(rx, partners):
    """What the arrow is coloured by: the ambient reactant that is not the third body."""
    if rx.form == "photolysis":
        return "light"
    name = next(s for s in partners if s != "M")
    return "hydrogen" if "H" in name else "carbon" if "C" in name else "oxygen"


def edge_label(rx):
    _, _, ins, outs = parts(rx)
    left = "$h\\nu$" if rx.form == "photolysis" else " + ".join(LABEL.get(s, s) for s in ins)
    right = f" → {' + '.join(LABEL[s] for s in outs)}" if outs else ""
    return f"{rx.rid}: {left}{right}"


def line_width(share):
    """Line width for a share of the inventory: logarithmic between FLUX_FLOOR and FLUX_CEIL, clipped at both."""
    if share <= FLUX_FLOOR:
        return LW_MIN
    x = min(math.log10(share / FLUX_FLOOR) / math.log10(FLUX_CEIL / FLUX_FLOOR), 1.0)
    return LW_MIN + x * (LW_MAX - LW_MIN)


def share_text(share):
    """A share of the inventory as short label text."""
    if share == 0.0:
        return "0"
    return f"{share:.0e}" if share < 1e-3 else f"{share:.2g}" if share < 10.0 else f"{share:.0f}"


def reaction_flux(case, settings=None):
    """Run one case and integrate its aluminium flows in time. Returns a dict: `reaction` {rid: share of the inventory
    carried over the run}, `start` {species: share Stage A hands over}, `formation`, `growth` and `smoke` (shares that
    left the gas each way), the end-of-run `budget`, `residual` {species: budget error as a share of the inventory} and
    the ambient `z_km` and `T_K`. Reads the Box's own rates, as plume_video does, so no rate is reimplemented."""
    settings = settings or dataclasses.replace(Settings(), points_per_decade=FLUX_POINTS_PER_DECADE)
    species = load_species()
    trail = trail_for(case.d0_um, case.h0_km, settings.n_emission_samples)
    column = AmbientColumn.load(settings.ambient, case.time_of_day)
    amb = column.at(trail.z_mean_km)
    hot = slow_down(trail.samples, column, settings, species, calibrate(load_reactions(stage=("hot", "both")), species))
    box = Box(trail, hot, amb, case.time_of_day, settings, species=species)
    sol = box.solve()
    lay, paths = box.layout, formation_pathways()
    n_rx, n_al = len(box.reactions), len(AL_SPECIES)
    net, growth = np.zeros((len(sol.t), n_rx + len(paths) + n_al)), np.zeros((len(sol.t), n_al))
    for i, x in enumerate(sol.Y):
        r2 = float(x[lay.r2][0])
        V = box._geometry(r2)[0]
        net[i] = box._network_rates(x[lay.gas], V)
        p = box._pstate(x, V, math.sqrt(r2))
        growth[i] = p["rho_ks"] @ p["N"]
    done = np.trapezoid(net, sol.t, axis=0) / box.N_Al                  # events (one Al atom each) per inventory atom
    grown = np.trapezoid(growth, sol.t, axis=0) / box.N_Al
    formed = done[n_rx:n_rx + len(paths)]
    smoke = done[n_rx + len(paths):]
    residual = {}
    for a, name in enumerate(AL_SPECIES):
        change = hot.fractions.get(name, 0.0) - sol.count(name, len(sol.t) - 1) / box.N_Al
        for j, rx in enumerate(box.reactions):
            al_in, al_out = parts(rx)[:2]
            change += done[j] * ((al_out == name) - (al_in == name))
        change -= sum(formed[q] * ((p.a == name) + (p.b == name)) for q, p in enumerate(paths))
        residual[name] = float(change - smoke[a] - grown[a])
    worst = max(abs(v) for v in residual.values())
    if worst > BALANCE_TOL:
        raise RuntimeError(f"the flux integrals miss the aluminium budget by {worst:.3g} of the inventory: {residual}")
    return {"reaction": {rx.rid: float(done[j]) for j, rx in enumerate(box.reactions)},
            "start": {s: float(hot.fractions.get(s, 0.0)) for s in HOT_START},
            "formation": float(2.0 * formed.sum()), "growth": float(grown.sum()), "smoke": float(smoke.sum()),
            "budget": sol.al_budget(len(sol.t) - 1), "residual": residual, "z_km": float(amb.z_km),
            "T_K": float(amb.T_K)}


def write_flux_csv(path, flux, rx_list):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "kind", "share_of_inventory"])
        for rx in rx_list:
            w.writerow([rx.rid, "reaction", flux["reaction"].get(rx.rid, 0.0)])
        for k, v in flux["start"].items():
            w.writerow([k, "handed over by Stage A", v])
        for k in ("formation", "growth", "smoke"):
            w.writerow([k, "left the gas", flux[k]])
        for k, v in flux["budget"].items():
            w.writerow([k, "at end of run", v])


def draw_reaction(ax, rx, patches, flux=None):
    if rx.rid not in ROUTES:
        raise ValueError(f"{rx.rid} has no entry in ROUTES: give the new reaction a place on the drawing")
    a, b, ins, _ = parts(rx)
    rad, t, dx, dy = ROUTES[rx.rid]
    (x0, y0), (x1, y1) = NODES[a], NODES[b]
    share = flux["reaction"].get(rx.rid, 0.0) if flux else None      # R17 is not in the network at night
    lw = line_width(share) if flux else 1.6
    colour = FAMILY_COLOUR[family(rx, ins)]
    style = (0, (5, 2.5)) if "M" in ins else (0, (1.5, 1.5)) if rx.form == "photolysis" else "solid"
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                                 mutation_scale=9 + 2.2 * lw if flux else 13, lw=lw, color=colour, linestyle=style,
                                 shrinkA=1, shrinkB=1,
                                 patchA=patches[a], patchB=patches[b], zorder=2, joinstyle="miter"))
    cx, cy = 0.5 * (x0 + x1) + rad * (y1 - y0), 0.5 * (y0 + y1) - rad * (x1 - x0)          # the arc's control point
    lx = (1 - t) ** 2 * x0 + 2 * t * (1 - t) * cx + t ** 2 * x1 + dx
    ly = (1 - t) ** 2 * y0 + 2 * t * (1 - t) * cy + t ** 2 * y1 + dy
    text = f"{edge_label(rx)}  |  {share_text(share)}" if flux else edge_label(rx)
    ax.text(lx, ly, text, ha="center", va="center", fontsize=FS_LABEL, color=INK, zorder=3,
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=colour, lw=0.7))


def box(ax, x, y, w, h, fc="white", ec=HAIR, lw=1.0, z=1):
    return ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.22", fc=fc, ec=ec, lw=lw,
                                       zorder=z))


def arrow(ax, p, q, rad=0.0, lw=1.6):
    ax.add_patch(FancyArrowPatch(p, q, connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>", lw=lw, color=INK2,
                                 mutation_scale=max(14, 9 + 2.2 * lw), shrinkA=0, shrinkB=0, zorder=2))


def draw_species(ax):
    """The species boxes, returned by name so that the arrows can stop at their edges."""
    patches = {}
    for name, (x, y) in NODES.items():
        hot = name in HOT_START
        patches[name] = box(ax, x - 1.05, y - 0.6, 2.1, 1.2, fc="#e4e3da" if hot else "white", ec=INK, lw=1.3, z=4)
        ax.text(x, y + 0.14, LABEL[name], ha="center", va="center", fontsize=FS_NODE, fontweight="bold", color=INK,
                zorder=5)
        ax.text(x, y - 0.34, NODE_NOTE[name], ha="center", va="center", fontsize=FS_LABEL - 1, color=INK2, zorder=5)
    return patches


def draw_hot_start(ax, hot_rids, patches, flux=None):
    y0, h = 11.0, 1.4
    box(ax, 0.0, y0, PANEL_X1, h, fc="white", ec=INK)
    ax.text(0.3, y0 + h - 0.32, "Stage A: hot start", fontsize=FS_HEAD, fontweight="bold", color=INK, va="center")
    ax.text(0.3, y0 + 0.42, "Aluminium vapour slows in the air, reacting in energetic collisions by "
            + ", ".join(hot_rids) + " and D1 (AlO + M → Al + O + M).\n"
            "It hands Stage B its aluminium as the three shaded species only.", fontsize=FS_TEXT, color=INK2,
            va="center", linespacing=1.25)
    for name in HOT_START:
        x, y = NODES[name]
        share = flux["start"][name] if flux else None
        lw = line_width(share) if flux else 1.6
        ax.add_patch(FancyArrowPatch((x, y0), (x, y), arrowstyle="-|>", mutation_scale=max(14, 9 + 2.2 * lw), lw=lw,
                                     color=INK2, patchB=patches[name], shrinkA=0, shrinkB=1, zorder=2))
        if flux:
            ax.text(x + 0.2, y0 - 0.2, share_text(share), ha="left", va="center", fontsize=FS_LABEL, color=INK)


def draw_condensation(ax, n_pathways, flux=None):
    n_al = len(AL_SPECIES)
    X, W = PANEL_X1 + 1.3, 6.2
    top, bottom = 12.4, 4.7
    box(ax, X, bottom, W, top - bottom, fc="white", ec=INK)
    ax.text(X + W / 2, top - 0.4, "Plume particles", ha="center", va="center", fontsize=FS_HEAD + 1, fontweight="bold",
            color=INK)
    ax.text(X + W / 2, top - 0.85, "size bins counted in Al atoms", ha="center", va="center", fontsize=FS_TEXT,
            color=INK2)
    y = top - 1.45
    for title, text in (
            ("Dimer formation", f"any two of the {n_al} species collide:\n{n_pathways} pathways, each a new\n"
                                "particle of two Al atoms"),
            ("Growth", f"any of the {n_al} species sticks,\nbringing its Al, O and H"),
            ("Oxygen uptake", "O, O$_2$, O$_3$ each add one O atom,\nslowing to stop at O/Al = 1.5"),
            ("Coagulation", "particle + particle")):
        ax.text(X + 0.3, y, title, fontsize=FS_HEAD - 0.5, fontweight="bold", color=INK, va="center")
        lines = text.count("\n") + 1
        ax.text(X + 0.3, y - 0.27, text, fontsize=FS_TEXT, color=INK2, va="top", linespacing=1.2)
        y -= 0.72 + 0.4 * lines
    box(ax, X, 2.25, W, 0.85, fc="white", ec=INK)
    ax.text(X + W / 2, 2.67, "Settled out of the plume", ha="center", va="center", fontsize=FS_HEAD, fontweight="bold",
            color=INK)
    box(ax, X, -0.2, W, 1.5, fc="white", ec=INK)
    ax.text(X + W / 2, 0.92, "Background meteoric smoke", ha="center", va="center", fontsize=FS_HEAD,
            fontweight="bold", color=INK)
    ax.text(X + W / 2, 0.3, f"takes up any of the {n_al} species;\nfixed, not evolved", ha="center", va="center",
            fontsize=FS_TEXT, color=INK2)
    if flux:                                            # the three flows out of the gas phase, and settling
        to_particles, to_smoke = flux["formation"] + flux["growth"], flux["smoke"]
        settled = flux["budget"]["settled"]
        for x, y, v, ha in ((PANEL_X1 + 0.65, 8.0, to_particles, "center"), (PANEL_X1 + 0.65, 1.1, to_smoke, "center"),
                            (X + W / 2 + 0.3, 3.7, settled, "left")):
            ax.text(x, y, share_text(v), ha=ha, va="center", fontsize=FS_LABEL, color=INK)
    arrow(ax, (PANEL_X1, 7.5), (X, 7.5), lw=line_width(to_particles) if flux else 1.6)
    arrow(ax, (PANEL_X1, 0.55), (X, 0.55), lw=line_width(to_smoke) if flux else 1.6)
    arrow(ax, (X + W / 2, bottom), (X + W / 2, 3.1), lw=line_width(settled) if flux else 1.6)


def _legend(ax, handles, title, x, ncol):
    return ax.legend(handles=handles, title=title, loc="upper left", ncol=ncol, frameon=False, handlelength=1.6,
                     columnspacing=1.1, handletextpad=0.5, bbox_to_anchor=(x, 0.109), fontsize=FS_TEXT,
                     title_fontsize=FS_TEXT, alignment="left", borderaxespad=0.0, borderpad=0.0)


def draw_legend(ax, flux=None):
    names = dict(FAMILY_NAME, light="sunlight")
    ax.add_artist(_legend(ax, [Line2D([], [], color=FAMILY_COLOUR[f], lw=2.6, label=names[f]) for f in FAMILY_COLOUR],
                          "Colour: ambient gas reacting", 0.01, 4))
    leg = _legend(ax, [Line2D([], [], color=INK2, lw=2.2, label="two-body"),
                       Line2D([], [], color=INK2, lw=2.2, ls=(0, (4, 2)), label="three-body"),
                       Line2D([], [], color=INK2, lw=2.2, ls=(0, (1.2, 1.2)), label="photolysis")],
                  "Line: kind of reaction", 0.40, 3)
    if flux:
        ax.add_artist(leg)
        _legend(ax, [Line2D([], [], color=INK2, lw=line_width(v), label=t)
                     for v, t in ((1e-4, "1e-4 or less"), (1e-2, "0.01"), (1.0, "1"), (100.0, "100"))],
                "Thickness: Al carried per Al atom released", 0.68, 4)


def plot_network(path_stem, flux=None, title="Wake-column reaction network", note=None):
    box_rx = load_reactions()
    hot_rids = [r.rid for r in box_rx if r.stage == "both"]
    pathways = formation_pathways()
    fig, ax = plt.subplots(figsize=(16.0, 9.0))
    ax.set_xlim(-0.3, 26.7)
    ax.set_ylim(-2.1, 13.09)
    ax.axis("off")
    ax.set_aspect("equal")
    box(ax, 0.0, -0.2, PANEL_X1, 10.9, fc=PANEL, ec="none", z=0)
    ax.text(0.3, 1.05, "Stage B: gas chemistry", fontsize=FS_HEAD, fontweight="bold", color=INK, va="center")
    ax.text(0.3, 0.3, f"{len(box_rx)} reactions among\n{len(AL_SPECIES)} aluminium species", fontsize=FS_TEXT,
            color=INK2, va="center", linespacing=1.2)
    patches = draw_species(ax)
    for rx in box_rx:
        draw_reaction(ax, rx, patches, flux)
    draw_hot_start(ax, hot_rids, patches, flux)
    draw_condensation(ax, len(pathways), flux)
    ax.text(0.0, 12.8, title, fontsize=FS_TITLE, fontweight="bold", color=INK, va="center")
    draw_legend(ax, flux)
    ax.text(0.0, -2.05, note or STRUCTURE_NOTE, fontsize=FS_LABEL - 1, color=INK2, va="bottom")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.995, bottom=0.005)
    for ext in ("png", "pdf"):
        fig.savefig(f"{path_stem}.{ext}", dpi=200, facecolor="white")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", default=OUTPUT, help="folder that holds plots/")
    parser.add_argument("--flux", action="store_true", help="run one case and draw thickness by aluminium carried")
    parser.add_argument("--d0", type=float, default=14.1, help="droplet diameter, micrometres")
    parser.add_argument("--h0", type=float, default=76.0, help="release altitude, km")
    parser.add_argument("--tod", choices=("day", "night"), default="day", help="time of day")
    args = parser.parse_args(argv)
    plots = os.path.join(args.out, "plots")
    os.makedirs(plots, exist_ok=True)
    if not args.flux:
        stem = os.path.join(plots, "reaction_network")
        plot_network(stem)
        print(f"wrote {stem}.png and .pdf")
        return
    flux = reaction_flux(Case(args.d0, args.h0, args.tod))
    stem = os.path.join(plots, f"reaction_flux_d{args.d0:g}um_h{args.h0:g}km_{args.tod}")
    title = f"Wake-column reaction flux: {args.d0:g} µm droplet, {args.h0:g} km, {args.tod}, 12 hours"
    plot_network(stem, flux, title, FLUX_NOTE)
    write_flux_csv(f"{stem}.csv", flux, load_reactions())
    print(f"wrote {stem}.png, .pdf and .csv")
    print(f"ambient {flux['z_km']:.1f} km, {flux['T_K']:.1f} K; largest budget residual "
          f"{max(abs(v) for v in flux['residual'].values()):.2e}; budget {flux['budget']}")


if __name__ == "__main__":
    sys.exit(main())
