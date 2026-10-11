"""Experiment X2 (Asha, 2026-10-09): prism-layer gridding, heating on smoothed normals, individual element recession.

Runs the box harness of analysis/box/box_ablation.py (imported, not edited: its run_specs, run_one, classify and
collect_variant) on the experiment copy experiments/x2/code, with
  - the box meshed as prism layers (mesh.layered_box_mesh: 2 mm triangles, layers of t = 0.5, 1.0, 2.0 mm, each prism
    split into three tetrahedra that die one by one -- nothing in the model kills a whole prism);
  - the heating computed on the smoothed surface's normals and carried onto the actual faces: per step
    n_s = MeltingBody.smoothed_normals_x2(), a patch's weight w = max(0, n . n_s) max(0, -n_s . x), its flux
    q A0 w / sum(w A), so that the total power is exactly q times the cross-section, as box_ablation's projected rule;
  - the switch cover_exposed_v2 of the copy (an exposed face that faces away from the smoothed surface gets no skin by
    the solidus test while it faces away).

Two deviations from the brief, both measured and reported (summary.md, the report):
  D1  smoothed_normals_x2 keeps each patch's own orientation (mesh.smoothed_normals(orient=False)). The model's
      smoothed_normals flips every smoothed normal onto the patch's own, so with it n . n_s is never negative: the
      brief's max(0, n . n_s) would do nothing and cover_exposed_v2 could never act.
  D2  the box's walls (the four side planes and the back plane) are held fixed by the smoothing (frozen=...). Without
      that the smoothing rounds the box's edges and the power normalisation sends 12-21 % of the power into the side
      walls, leaving the lateral interior 11.5-15.6 % short of 6 MW/m2 on the intact box (`--probe` prints it); with
      the walls held, an intact flat front gets exactly the projected rule, which is the brief's own premise.
  The heating as the brief literally wrote it (neither deviation) is run as context for the element model ("literal").

Variant q6p (6 MW/m2, 2 s, Scheil, 700 K), dev preset (20 um cells), sub-step 10 ms, direct linear solver.

    python x2_box.py --list
    python x2_box.py --probe                          # the heating rules on the intact boxes
    python x2_box.py skin__t0.5__dt0.05 ...           # run the named runs, one JSON each under runs/
    python x2_box.py --collect                        # runs/*.json -> results.json, summary.md, x2_box.png
"""
import glob
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "code")
BOX = os.path.normpath(os.path.join(HERE, "..", "..", "analysis", "box"))
X1_RUNS = os.path.normpath(os.path.join(HERE, "..", "x1", "runs"))
RUNS = os.path.join(HERE, "runs")
sys.dont_write_bytecode = True
sys.path.insert(0, CODE)
import reentry_model  # noqa: E402  the experiment copy, imported before box_ablation puts meas/ on the path

sys.path.insert(0, BOX)
import box_ablation as bx  # noqa: E402

import numpy as np  # noqa: E402

from reentry_model import body, heating, material, mesh, skin, thermal  # noqa: E402

for _m in (reentry_model, bx.body, bx.skin, bx.mesh, bx.thermal, bx.material, bx.heating, body, mesh):
    assert os.path.realpath(_m.__file__).startswith(os.path.realpath(CODE)), _m.__file__

TAG = "q6p"
V = bx.VARIANTS[TAG]
BOXDIM = bx.COARSE                       # (8e-3, 10e-3, 10e-3, 2e-3): 8 mm deep, 10 x 10 mm front, 2 mm triangles
LAYERS = (0.5e-3, 1.0e-3, 2.0e-3)
MACRO = (0.05, 0.5)
MODELS = {"elem": ("elements", False), "skin": ("skin", False), "skin_v2": ("skin", True)}
FEED_ZERO = 0.0                          # kg: a step whose lateral liquid is at or below this fed nothing
SLOW_RATE = 1.1e-3                       # m/s: X1's softer stall measure (a quarter of the fine box's 1-2 s mean rate)
ONSET_MM = 0.05
AWAY_NX = -0.1                           # X1's convention: n_x >= -0.1 does not face the heat
WALL_TOL = 1e-12
COLOURS = {"elem": "tab:blue", "skin": "tab:orange", "skin_v2": "tab:green", "elem_literal": "tab:gray", "skin_nowall": "tab:red"}


def tname(t):
    return "t{:g}".format(t * 1e3)


def specs():
    out = {}
    common = dict(tag=TAG, q=V["q"], rule=bx.RULE_PROJECTED, t_end=bx.T_END, box=BOXDIM)
    for t in LAYERS:
        for mname, (model, v2) in MODELS.items():
            for dt in MACRO:
                out["{}__{}__dt{:g}".format(mname, tname(t), dt)] = dict(
                    common, model=model, dt=dt, preset="dev" if model == "skin" else None, layer_t=t, heating="smoothed",
                    cover_exposed_v2=v2, group=mname)
        out["elem_literal__{}__dt0.05".format(tname(t))] = dict(       # context: the brief's rule without D1 and D2
            common, model="elements", dt=0.05, preset=None, layer_t=t, heating="literal", cover_exposed_v2=False,
            group="elem_literal")
        for dt in MACRO:      # diagnostic beyond the brief: the box's wall planes (adiabatic symmetry planes of the slab it
            # stands for) left out of the skins' solidus test, to separate the wall artefact from the front's mechanism
            out["skin_nowall__{}__dt{:g}".format(tname(t), dt)] = dict(
                common, model="skin", dt=dt, preset="dev", layer_t=t, heating="smoothed", cover_exposed_v2=False,
                nowall_skins=True, group="skin_nowall")
    # the reproduction check: the copy, switches off, on X1's Kuhn box with box_ablation's projected rule
    base = bx.run_specs(TAG, V["q"], V["rule"])["skin_dev_dt0.05"]
    out["kuhn_base__dt0.05"] = dict(base, layer_t=None, heating="projected", cover_exposed_v2=False, group="kuhn_base",
                                    x1="base__dt0.05")
    return out


def wall_nodes(msh, box):
    """Nodes on the box's four side planes and its back plane x = lx: the frame D2 holds fixed in the smoothing."""
    lx, ly, lz, _ = box
    p = msh.points
    return np.flatnonzero((np.abs(p[:, 1]) < WALL_TOL) | (np.abs(p[:, 1] - ly) < WALL_TOL) | (np.abs(p[:, 2]) < WALL_TOL)
                          | (np.abs(p[:, 2] - lz) < WALL_TOL) | (np.abs(p[:, 0] - lx) < WALL_TOL))


def wall_patches(s, box):
    lx, ly, lz, _ = box
    c, n = s.centroids, s.normals
    return (((np.abs(n[:, 1]) > 0.99) & ((np.abs(c[:, 1]) < 1e-9) | (np.abs(c[:, 1] - ly) < 1e-9)))
            | ((np.abs(n[:, 2]) > 0.99) & ((np.abs(c[:, 2]) < 1e-9) | (np.abs(c[:, 2] - lz) < 1e-9)))
            | ((n[:, 0] > 0.99) & (np.abs(c[:, 0] - lx) < 1e-9)))


def smoothed_weights(s, ns):
    """The brief's weight per unit of patch area, w = max(0, n . n_s) max(0, -n_s . x), and n . n_s."""
    d = np.einsum("ij,ij->i", s.normals, ns)
    return np.maximum(0.0, d) * np.maximum(0.0, -ns[:, 0]), d


def heating_flux(b, spec):
    """Per-patch flux [W/m2] of the run's heating rule on the current surface, and the smoothed normals it used."""
    s = b.surface
    lx, ly, lz, h = spec["box"]
    A0 = ly * lz
    if spec["heating"] == "projected":
        w, ns = np.maximum(0.0, -s.normals[:, 0]), None
    else:
        ns = b.smoothed_normals_x2() if spec["heating"] == "smoothed" else s.smoothed_normals()
        w, _ = smoothed_weights(s, ns)
    return spec["q"] * A0 * w / float((w * s.areas).sum()), ns


_last = {}


def make_body(spec):
    """box_ablation's body (its own make_body on the Kuhn box; the same construction on the layered box), the switch set,
    D2's frame given to the body, and its advance wrapped: the heating replaced by the run's rule, and the X2 diagnostics
    recorded after every step (before the harness removes the liquid). Recording only reads."""
    lx, ly, lz, h = spec["box"]
    if spec.get("layer_t"):
        msh = mesh.layered_box_mesh(lx, ly, lz, h, spec["layer_t"])
        mat = material.Material.from_drama_json("AA7075_scheil")
        sk = skin.SkinSettings(cell=bx.PRESETS[spec["preset"]], substep=bx.SUBSTEP) if spec["model"] == "skin" else None
        settings = body.MeltSettings(removal="girin", size_feedback="initial", surface_model=spec["model"], skin=sk)
        solver = thermal.thermal_solver("skfem", linear_solver="direct")
        b = body.MeltingBody(msh, mat, solver, float(mat.rho * msh.volume()), settings=settings, T0=bx.T_START, emissivity=0.0)
    else:
        b, mat = _orig_make_body(spec)
        msh = b.mesh
    b.settings.cover_exposed_v2 = bool(spec.get("cover_exposed_v2", False))
    if spec["heating"] == "smoothed":
        b.smooth_frozen_nodes = wall_nodes(msh, spec["box"])
    if spec.get("nowall_skins"):                    # the diagnostic: no solidus skins on the box's wall planes
        b.x2_solidus_excluded = lambda s: wall_patches(s, spec["box"])
    ce = msh.points[msh.tets].mean(axis=1)
    lat_el = (ce[:, 1] > h) & (ce[:, 1] < ly - h) & (ce[:, 2] > h) & (ce[:, 2] < lz - h)
    A_lat0 = (ly - 2.0 * h) * (lz - 2.0 * h)
    keys = ("t", "q_lat_proj_MWm2", "q_lat_own_MWm2", "P_lat_share", "P_lat_share_projected", "P_wall_share", "ns_dev_lat_max_deg",
            "ns_dev_lat_mean_deg", "tilt_lat_power_deg", "tilt_lat_max_deg", "folded_n", "folded_P_share_literal",
            "fed_lat_kg", "fed_total_kg", "n_active", "n_dead_lat", "skins_total", "skins_lat", "skins_lat_heated",
            "skins_lat_away_nx", "skins_lat_away_nx_exposed", "skins_lat_away_back", "skins_lat_away_ns",
            "skins_lat_exposed_noheat", "skins_created", "cover_v2_withheld", "uncovered_lat_now", "exposed_away_ns_now",
            "exposed_away_ns_literal_now")
    rec = {k: [] for k in keys}
    rec["uncovered_first"] = {}          # face id -> [t, n_x, n_y, n_z, area] when first seen on the lateral interior
    inner = b.advance

    def advance(t, dt, loads, state=None):
        s = b.surface
        if spec["heating"] != "projected":
            q, ns = heating_flux(b, spec)
            loads = heating.HeatingResult(q, spec["q"], 0.0, 0.0, 0.0)
        else:
            q, ns = loads.q_conv, None
        _, _, _, lateral = bx.classify(b, spec["box"])
        P = q * s.areas
        Ptot = spec["q"] * ly * lz
        rec["t"].append(t)
        rec["q_lat_proj_MWm2"].append(float(P[lateral].sum() / A_lat0 / 1e6))
        heated_lat = lateral & (q > 0.0)
        rec["q_lat_own_MWm2"].append(float(P[heated_lat].sum() / s.areas[heated_lat].sum() / 1e6) if heated_lat.any() else 0.0)
        rec["P_lat_share"].append(float(P[lateral].sum() / Ptot))
        wp = np.maximum(0.0, -s.normals[:, 0])
        rec["P_lat_share_projected"].append(float((wp * s.areas)[lateral].sum() / (wp * s.areas).sum()))
        rec["P_wall_share"].append(float(P[wall_patches(s, spec["box"])].sum() / Ptot))
        tilt = np.degrees(np.arccos(np.clip(-s.normals[:, 0], -1.0, 1.0)))
        rec["tilt_lat_power_deg"].append(float((tilt * P)[heated_lat].sum() / P[heated_lat].sum()) if heated_lat.any() else float("nan"))
        rec["tilt_lat_max_deg"].append(float(tilt[heated_lat].max()) if heated_lat.any() else float("nan"))
        if ns is not None:
            dev = np.degrees(np.arccos(np.clip(-ns[:, 0], -1.0, 1.0)))
            rec["ns_dev_lat_max_deg"].append(float(dev[heated_lat].max()) if heated_lat.any() else float("nan"))
            rec["ns_dev_lat_mean_deg"].append(float((dev * P)[heated_lat].sum() / P[heated_lat].sum()) if heated_lat.any() else float("nan"))
            d = np.einsum("ij,ij->i", s.normals, ns)
            folded = d <= 0.0
            rec["folded_n"].append(int(folded.sum()))
            nsl = np.where(folded[:, None], -ns, ns)               # what the model's own flip would have read
            wl, _ = smoothed_weights(s, nsl)
            rec["folded_P_share_literal"].append(float((wl * s.areas)[folded].sum() / (wl * s.areas).sum()))
        else:
            for k in ("ns_dev_lat_max_deg", "ns_dev_lat_mean_deg", "folded_P_share_literal"):
                rec[k].append(float("nan"))
            rec["folded_n"].append(0)
        inner(t, dt, loads, state)
        # -- after the step (the liquid not yet removed) --
        s = b.surface
        _, _, _, lateral = bx.classify(b, spec["box"])
        rec["fed_lat_kg"].append(float((b.m_f + b.m_d)[lateral].sum()))
        rec["fed_total_kg"].append(float((b.m_f + b.m_d).sum()))
        rec["n_active"].append(int(b.mesh.n_active))
        rec["n_dead_lat"].append(int((lat_el & ~b.mesh.active).sum()))
        exp_ids = np.fromiter(b.exposed_faces, dtype=np.int64, count=len(b.exposed_faces))
        exposed = np.isin(s.face_ids, exp_ids)
        ns2 = b.smoothed_normals_x2()
        w2, d2 = smoothed_weights(s, ns2)
        d_lit = np.einsum("ij,ij->i", s.normals, s.smoothed_normals(frozen=b.smooth_frozen_nodes))
        back = s.centroids[:, 0] > spec["box"][0] * (1.0 - 1e-9)
        for k in np.flatnonzero(lateral & exposed):
            f = int(s.face_ids[k])
            if f not in rec["uncovered_first"]:
                rec["uncovered_first"][f] = [t, float(s.normals[k, 0]), float(s.normals[k, 1]), float(s.normals[k, 2]), float(s.areas[k])]
        rec["uncovered_lat_now"].append(int((lateral & exposed).sum()))
        rec["exposed_away_ns_now"].append(int((exposed & (d2 <= 0.0)).sum()))
        rec["exposed_away_ns_literal_now"].append(int((exposed & (d_lit <= 0.0)).sum()))
        if b.skins is not None and b.skins.n:
            on = b.skin_of_patch >= 0
            nx = s.normals[:, 0]
            lat_on = lateral & on
            rec["skins_lat"].append(int(lat_on.sum()))
            rec["skins_lat_heated"].append(int((lat_on & (nx < AWAY_NX)).sum()))
            rec["skins_lat_away_nx"].append(int((lat_on & (nx >= AWAY_NX)).sum()))
            rec["skins_lat_away_nx_exposed"].append(int((lat_on & (nx >= AWAY_NX) & exposed).sum()))
            rec["skins_lat_away_back"].append(int((lat_on & (nx >= AWAY_NX) & back).sum()))
            rec["skins_lat_away_ns"].append(int((lat_on & (d2 <= 0.0)).sum()))
            rec["skins_lat_exposed_noheat"].append(int((lat_on & exposed & (w2 <= 0.0)).sum()))
        else:
            for k in ("skins_lat", "skins_lat_heated", "skins_lat_away_nx", "skins_lat_away_nx_exposed", "skins_lat_away_back",
                      "skins_lat_away_ns", "skins_lat_exposed_noheat"):
                rec[k].append(0)
        rec["skins_total"].append(int(b.skins.n) if b.skins is not None else 0)
        rec["skins_created"].append(int(b.skins_created))
        rec["cover_v2_withheld"].append(int(b.cover_v2_withheld))

    b.advance = advance
    _last["body"], _last["rec"] = b, rec
    return b, mat


_orig_make_body = bx.make_body
bx.make_body = make_body     # run_one looks make_body up in its module at call time


def run(name, spec):
    print("== X2 {} ==".format(name), flush=True)
    try:
        r = bx.run_one(name, spec, log=lambda s: print(s, flush=True))
        b, x = _last["body"], _last["rec"]
        r["series"]["heated_fraction"] = [1.0] * len(r["series"]["t"])   # every rule here delivers exactly q A0
        r["x2"] = x
        r["mass_balance_rel"] = abs(r["mass_balance"]) / r["mass0_kg"]
        r["cover_v2_withheld"] = int(b.cover_v2_withheld)
        r["exposed_faces_total"] = len(b.exposed_faces)
        r["code"] = os.path.dirname(reentry_model.__file__)
        r["status"] = "ok"
    except Exception as exc:
        import traceback
        traceback.print_exc()
        r = {"name": name, "tag": TAG, "status": "failed", "error": "{}: {}".format(type(exc).__name__, exc)}
    r["group"] = spec["group"]
    r["x2_spec"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in spec.items()}
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, name + ".json"), "w") as f:
        json.dump(r, f, indent=1)


# -- the heating rules on the intact boxes ---------------------------------------------------------------------------
def probe():
    lx, ly, lz, h = BOXDIM
    q, A0 = V["q"], ly * lz
    A_lat0 = (ly - 2 * h) * (lz - 2 * h)
    out = {}
    for t in (None,) + LAYERS:
        msh = mesh.box_mesh(lx, ly, lz, h) if t is None else mesh.layered_box_mesh(lx, ly, lz, h, t)
        s = msh.surface()
        c = s.centroids
        lat = (c[:, 1] > h * (1 - 1e-9)) & (c[:, 1] < ly - h * (1 - 1e-9)) & (c[:, 2] > h * (1 - 1e-9)) & (c[:, 2] < lz - h * (1 - 1e-9))
        front = s.normals[:, 0] < -0.99
        walls = wall_patches(s, BOXDIM)
        row = {}
        for rule, ns in (("literal", s.smoothed_normals()), ("smoothed", s.smoothed_normals(orient=False, frozen=wall_nodes(msh, BOXDIM)))):
            w, _ = smoothed_weights(s, ns)
            qq = q * A0 * w / (w * s.areas).sum()
            P = qq * s.areas
            dev = np.degrees(np.arccos(np.clip(-ns[:, 0], -1.0, 1.0)))
            row[rule] = {"q_lat_proj_MWm2": float(P[lat].sum() / A_lat0 / 1e6), "P_wall_share": float(P[walls].sum() / (q * A0)),
                         "ns_dev_front_max_deg": float(dev[front].max()), "ns_dev_lat_max_deg": float(dev[front & lat].max()),
                         "q_front_min_MWm2": float(qq[front].min() / 1e6), "q_front_max_MWm2": float(qq[front].max() / 1e6)}
        out["kuhn" if t is None else tname(t)] = row
        print("kuhn" if t is None else tname(t), json.dumps(row))
    return out


# -- collection -----------------------------------------------------------------------------------------------------
def stall(t, fed, dt):
    """X1's stall: the longest interval with no lateral feed after the first lateral feed (duration, from, to) [s]."""
    fed = np.asarray(fed)
    on = np.flatnonzero(fed > FEED_ZERO)
    if not on.size:
        return None, None, None
    best, start = (0, None), None
    for k in range(on[0] + 1, len(fed)):
        if fed[k] <= FEED_ZERO:
            start = k if start is None else start
            if k - start + 1 > best[0]:
                best = (k - start + 1, start)
        else:
            start = None
    if best[0] == 0:
        return 0.0, None, None
    n, k0 = best
    return n * dt, t[k0] - dt, t[k0 + n - 1]


def slow(t, fed, dt, A_lat):
    """X1's softer measure: the longest interval after the first lateral feed with the lateral rate below SLOW_RATE."""
    rate = np.asarray(fed) / (A_lat * bx.RHO_METAL) / dt
    on = np.flatnonzero(np.asarray(fed) > FEED_ZERO)
    if not on.size:
        return None, None, None
    best, start = (0, None), None
    for k in range(on[0] + 1, len(rate)):
        if rate[k] < SLOW_RATE:
            start = k if start is None else start
            if k - start + 1 > best[0]:
                best = (k - start + 1, start)
        else:
            start = None
    if best[0] == 0:
        return 0.0, None, None
    n, k0 = best
    return n * dt, t[k0] - dt, t[k0 + n - 1]


def onset(t, r, level=ONSET_MM):
    """The first time the recession reaches `level` mm, linear between the steps (0 at t = 0)."""
    tt = np.concatenate([[0.0], t])
    rr = np.concatenate([[0.0], r])
    k = np.flatnonzero(rr >= level)
    if not k.size:
        return None
    k = int(k[0])
    return float(tt[k - 1] + (level - rr[k - 1]) / (rr[k] - rr[k - 1]) * (tt[k] - tt[k - 1]))


def at_time(t, v, x):
    i = int(np.argmin(np.abs(np.asarray(t) - x)))
    return v[i] if abs(t[i] - x) < 1e-9 else None


def tilt_stats(first):
    """The tilt of the faces uncovered on the lateral interior, each counted once (when first seen): the angle [deg]
    between the heat-facing ones' (n_x < -0.1) outward normal and -x, and how many were vertical or faced back."""
    if not first:
        return {"n_uncovered": 0}
    a = np.array(list(first.values()))
    nx = a[:, 1]
    heat = nx < AWAY_NX
    ang = np.degrees(np.arccos(np.clip(-nx[heat], -1.0, 1.0)))
    bins = [0, 5, 15, 25, 35, 50, 80]
    hist = np.histogram(ang, bins=bins)[0].tolist() if ang.size else [0] * (len(bins) - 1)
    return {"n_uncovered": int(len(a)), "n_heat_facing": int(heat.sum()), "n_vertical": int((np.abs(nx) <= -AWAY_NX).sum()),
            "n_facing_back": int((nx > -AWAY_NX).sum()), "angle_bins_deg": bins, "angle_hist": hist,
            "angle_p10_deg": float(np.percentile(ang, 10)) if ang.size else None,
            "angle_median_deg": float(np.median(ang)) if ang.size else None,
            "angle_p90_deg": float(np.percentile(ang, 90)) if ang.size else None,
            "angle_max_deg": float(ang.max()) if ang.size else None,
            "distinct_angles_deg": sorted({round(float(x), 1) for x in ang})}


def deaths_before_feed(x, n_elements):
    """Elements dead (body-wide, and in the lateral column) at the end of the last step before the first step that
    produced any liquid (body-wide) / any lateral liquid."""
    ft, fl = np.asarray(x["fed_total_kg"]), np.asarray(x["fed_lat_kg"])
    na, nl = np.asarray(x["n_active"]), np.asarray(x["n_dead_lat"])
    k = np.flatnonzero(ft > 0.0)
    kl = np.flatnonzero(fl > 0.0)
    body_dead = None if not k.size else (0 if k[0] == 0 else int(n_elements - na[k[0] - 1]))
    lat_dead = None if not kl.size else (0 if kl[0] == 0 else int(nl[kl[0] - 1]))
    return body_dead, lat_dead, (x["t"][k[0]] if k.size else None), (x["t"][kl[0]] if kl.size else None)


def collect():
    with open(os.path.join(bx.RUNS, "q6p__fine.json")) as f:
        fine = json.load(f)
    recs = {}
    for path in sorted(glob.glob(os.path.join(RUNS, "*.json"))):
        with open(path) as f:
            r = json.load(f)
        recs[r["name"]] = r
    ok = {n: r for n, r in recs.items() if r.get("status") == "ok"}
    res = bx.collect_variant(TAG, dict(ok, fine=fine))
    sp = specs()
    fs = fine["series"]
    out = {"fine_recession_mm_at_times": res["fine_recession_mm_at_times"], "compare_times_s": res["compare_times_s"],
           "fine_onset_s": onset(fs["t"], fs["r_total_mm"]),
           "criterion": "lateral recession within 5 % of the fine box's at 0.5, 1.0, 1.5 and 2.0 s (where the fine box has >= 0.05 mm)",
           "runs": {}}
    for name, r in recs.items():
        if r.get("status") != "ok":
            out["runs"][name] = {"status": "failed", "error": r.get("error"), "group": r.get("group")}
            continue
        o = res["runs"][name]
        S, x = r["series"], r["x2"]
        spec = r["x2_spec"]
        dt = spec["dt"]
        lx, ly, lz, h = spec["box"]
        A_lat0 = (ly - 2 * h) * (lz - 2 * h)
        d, a, b_ = stall(S["t"], x["fed_lat_kg"], dt)
        ds, as_, bs_ = slow(S["t"], x["fed_lat_kg"], dt, A_lat0)
        dbf, dbf_lat, t_feed, t_feed_lat = deaths_before_feed(x, r["n_elements"])
        on_s = onset(S["t"], S["r_lat_mm"])
        k_on = 0 if on_s is None else int(np.searchsorted(S["t"], on_s))
        tp = np.asarray(x["tilt_lat_power_deg"], dtype=float)
        e = {"status": "ok", "group": r["group"], "layer_t_mm": None if not spec.get("layer_t") else spec["layer_t"] * 1e3,
             "dt": dt, "heating": spec["heating"], "cover_exposed_v2": spec.get("cover_exposed_v2", False),
             "r_lat_mm_at_times": o["r_lat_mm_at_times"], "rel_diff_at_times": o["r_lat_mm_rel_diff_at_times"],
             "max_abs_rel_diff_4times": o["r_lat_mm_max_abs_rel_diff"], "max_abs_diff_curve_mm": o["r_lat_mm_max_abs_diff_curve_mm"],
             "meets_5pct": all(abs(v) <= 0.05 for v in o["r_lat_mm_rel_diff_at_times"] if v is not None),
             "onset_s": on_s, "stall_s": d, "stall_from_s": a, "stall_to_s": b_, "slow_s": ds, "slow_from_s": as_, "slow_to_s": bs_,
             "first_feed_s": t_feed, "first_lateral_feed_s": t_feed_lat,
             "deaths_before_any_feed": dbf, "lateral_deaths_before_lateral_feed": dbf_lat,
             "tilt": tilt_stats(x["uncovered_first"]),
             "tilt_lat_power_after_onset_mean_deg": float(np.nanmean(tp[k_on:])) if tp[k_on:].size and np.isfinite(tp[k_on:]).any() else None,
             "tilt_lat_power_max_deg": float(np.nanmax(tp)) if np.isfinite(tp).any() else None,
             "q_lat_proj_MWm2_mean": float(np.mean(x["q_lat_proj_MWm2"])), "q_lat_proj_MWm2_min": float(np.min(x["q_lat_proj_MWm2"])),
             "q_lat_proj_MWm2_max": float(np.max(x["q_lat_proj_MWm2"])),
             "q_lat_own_MWm2_min": float(np.min(x["q_lat_own_MWm2"])),
             "P_wall_share_max": float(np.max(x["P_wall_share"])),
             "ns_dev_lat_max_deg": float(np.nanmax(x["ns_dev_lat_max_deg"])) if np.isfinite(x["ns_dev_lat_max_deg"]).any() else None,
             "ns_dev_lat_mean_deg_mean": float(np.nanmean(x["ns_dev_lat_mean_deg"])) if np.isfinite(x["ns_dev_lat_mean_deg"]).any() else None,
             "folded_patch_steps": int(np.sum(x["folded_n"])),
             "folded_P_share_literal_max": float(np.nanmax(x["folded_P_share_literal"])) if np.isfinite(x["folded_P_share_literal"]).any() else None,
             "exposed_away_ns_max": int(max(x["exposed_away_ns_now"])), "exposed_away_ns_literal_max": int(max(x["exposed_away_ns_literal_now"])),
             "cover_v2_withheld": r["cover_v2_withheld"],
             "skins_lat_away_nx_at_1s": at_time(x["t"], x["skins_lat_away_nx"], 1.0),
             "skins_lat_away_nx_at_2s": at_time(x["t"], x["skins_lat_away_nx"], 2.0),
             "skins_lat_away_nx_exposed_at_1s": at_time(x["t"], x["skins_lat_away_nx_exposed"], 1.0),
             "skins_lat_away_nx_exposed_at_2s": at_time(x["t"], x["skins_lat_away_nx_exposed"], 2.0),
             "skins_lat_away_back_at_1s": at_time(x["t"], x["skins_lat_away_back"], 1.0),
             "skins_lat_away_back_at_2s": at_time(x["t"], x["skins_lat_away_back"], 2.0),
             "skins_lat_away_nx_exposed_max": int(max(x["skins_lat_away_nx_exposed"])),
             "skins_lat_away_ns_max": int(max(x["skins_lat_away_ns"])),
             "skins_lat_exposed_noheat_max": int(max(x["skins_lat_exposed_noheat"])),
             "skins_lat_exposed_noheat_at_1s": at_time(x["t"], x["skins_lat_exposed_noheat"], 1.0),
             "skins_lat_heated_at_1s": at_time(x["t"], x["skins_lat_heated"], 1.0),
             "skins_lat_heated_at_2s": at_time(x["t"], x["skins_lat_heated"], 2.0),
             "worst_energy_residual": r["worst_energy_residual"], "mass_balance_rel": r["mass_balance_rel"],
             "books_close": bool(abs(r["worst_energy_residual"]) < 1e-8 and r["mass_balance_rel"] < 1e-12),
             "repeats_per_step": r["interface_repeats_total"] / r["n_steps"], "max_shortfall_kg": r["max_shortfall_kg"],
             "skins_created": r["skins_created"], "n_active_final": S["n_active"][-1], "n_elements": r["n_elements"],
             "run_time_s": r["run_time_s"], "skin_wall_s": r["skin_wall_s"],
             "T_top_lat_at_times": o["T_top_at_times"], "T_top_max_K": o.get("T_top_max_K"), "T_plateau_std_K": o.get("T_plateau_std_K")}
        if spec.get("x1"):
            q = json.load(open(os.path.join(X1_RUNS, spec["x1"] + ".json")))["series"]
            e["reproduces_x1"] = {
                "file": "experiments/x1/runs/" + spec["x1"] + ".json",
                "max_abs_diff_r_lat_mm": float(np.max(np.abs(np.array(S["r_lat_mm"]) - np.array(q["r_lat_mm"])))),
                "max_abs_diff_r_total_mm": float(np.max(np.abs(np.array(S["r_total_mm"]) - np.array(q["r_total_mm"])))),
                "bit_identical": all(S[k] == q[k] for k in ("r_lat_mm", "r_total_mm", "n_skins", "n_active", "T_top_lat", "residual"))}
        out["runs"][name] = e
    # each skin run against the element model on the same mesh and macro step
    for name, e in out["runs"].items():
        if e["status"] != "ok" or e["group"] not in ("skin", "skin_v2", "skin_nowall"):
            continue
        ref = "elem__{}__dt{:g}".format(tname(e["layer_t_mm"] * 1e-3), e["dt"])
        if ref in recs and recs[ref].get("status") == "ok":
            a = np.array(recs[name]["series"]["r_lat_mm"])
            bb = np.array(recs[ref]["series"]["r_lat_mm"])
            er = out["runs"][ref]
            e["vs_elem"] = {"elem_run": ref, "elem_max_abs_rel_diff_4times": er["max_abs_rel_diff_4times"],
                            "elem_max_abs_diff_curve_mm": er["max_abs_diff_curve_mm"],
                            "max_abs_diff_to_elem_curve_mm": float(np.max(np.abs(a - bb))),
                            "closer_to_fine_at_times": [None if (x is None or y is None) else abs(x) < abs(y)
                                                        for x, y in zip(e["rel_diff_at_times"], er["rel_diff_at_times"])]}
        if e["group"] == "skin_v2":
            twin = name.replace("skin_v2__", "skin__")
            if twin in recs and recs[twin].get("status") == "ok":
                e["identical_to_skin_run"] = recs[name]["series"]["r_lat_mm"] == recs[twin]["series"]["r_lat_mm"]
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out, recs, fine


def pct(x):
    return "n/a" if x is None else "{:+.1f} %".format(100 * x)


def verdict(e, out):
    rel = [(t, v) for t, v in zip(bx.COMPARE_TIMES, e["rel_diff_at_times"]) if v is not None]
    worst_t, worst = max(rel, key=lambda tv: abs(tv[1]))
    good = [t for t, v in rel if abs(v) <= 0.05]
    if len(good) == len(rel):
        s = "Meets the 5 % criterion: its largest difference from the fine box is {} at {:g} s.".format(pct(worst), worst_t)
    else:
        s = "Fails the 5 % criterion: its largest difference from the fine box is {} at {:g} s, and it is within 5 % only at {}.".format(
            pct(worst), worst_t, ", ".join("{:g} s".format(t) for t in good) if good else "none of the compared times")
    s += " Its recession reaches 0.05 mm at {} (the fine box at {:.2f} s).".format(
        "never" if e["onset_s"] is None else "{:.2f} s".format(e["onset_s"]), out["fine_onset_s"])
    ve = e.get("vs_elem")
    if ve:
        better = ve["elem_max_abs_rel_diff_4times"] is not None and e["max_abs_rel_diff_4times"] < ve["elem_max_abs_rel_diff_4times"]
        s += (" Against the element model on the same mesh and step ({} largest difference {}), it is {}; the two curves differ by at "
              "most {:.3f} mm.".format(ve["elem_run"], "n/a" if ve["elem_max_abs_rel_diff_4times"] is None else "{:.1f} %".format(100 * ve["elem_max_abs_rel_diff_4times"]),
                                       "closer to the fine box" if better else "further from the fine box", ve["max_abs_diff_to_elem_curve_mm"]))
    if e.get("identical_to_skin_run"):
        s += " Its recession series is bit for bit that of the same run without cover_exposed_v2: the switch never changed a creation that mattered."
    return s


def order(out):
    names = []
    for t in LAYERS:
        for dt in MACRO:
            names += ["{}__{}__dt{:g}".format(m, tname(t), dt) for m in MODELS]
        names.append("elem_literal__{}__dt0.05".format(tname(t)))
        names += ["skin_nowall__{}__dt{:g}".format(tname(t), dt) for dt in MACRO]
    names.append("kuhn_base__dt0.05")
    return [n for n in names if n in out["runs"]]


def write_summary(out):
    L = ["# Experiment X2: prism-layer box, heating on smoothed normals, individual element recession (variant q6p)", "",
         "Generated by `x2_box.py --collect` from `runs/`. Lateral recession (`r_lat_mm`: the liquid removed from the lateral interior, "
         "the central 6 x 6 mm of the 10 x 10 mm front, over its projected area and 2813 kg/m3) against the fine box's `r_total_mm` "
         "(analysis/box/runs/q6p__fine.json), which is {} mm at 0.5 / 1.0 / 1.5 / 2.0 s and first reaches 0.05 mm at {:.3f} s. "
         "Relative differences are (coarse - fine) / fine where the fine box has at least 0.05 mm, so not at 0.5 s. "
         "Heating `smoothed` is the brief's rule with deviations D1 and D2 (see x2_box.py); `literal` is the brief's rule as written.".format(
             " / ".join("{:.4f}".format(v) for v in out["fine_recession_mm_at_times"]), out["fine_onset_s"]), ""]
    L += ["| run | recession at 0.5 / 1.0 / 1.5 / 2.0 s (mm) | rel. diff at 1.0 / 1.5 / 2.0 s | largest abs rel. diff | largest abs diff over the curve (mm) | "
          "onset, 0.05 mm (s) | stall: no lateral feed (s) | slow: below 1.1 mm/s (s) | deaths before any feed (body / lateral column) | "
          "lateral skins with n_x >= -0.1 at 1 / 2 s [on exposed faces; back face] | max lateral skins on exposed faces with no heat | "
          "max lateral skins with n.n_s <= 0 | cover_v2 withheld | worst energy residual | mass books (rel.) | repeats per step | max shortfall (kg) | "
          "skins created | active elements at end | run time (s) |", "|" + "---|" * 20]
    for n in order(out):
        e = out["runs"][n]
        if e["status"] != "ok":
            L.append("| {} | FAILED: {} |".format(n, e.get("error")) + " |" * 18)
            continue
        st = "{:.2f} ({:.2f}-{:.2f})".format(e["stall_s"], e["stall_from_s"], e["stall_to_s"]) if e["stall_from_s"] is not None else (
            "none" if e["stall_s"] == 0.0 else "no feed")
        sl = "{:.2f} ({:.2f}-{:.2f})".format(e["slow_s"], e["slow_from_s"], e["slow_to_s"]) if e["slow_from_s"] is not None else (
            "none" if e["slow_s"] == 0.0 else "no feed")
        L.append("| {} | {} | {} | {} | {:.4f} | {} | {} | {} | {} / {} | {} / {} [{} / {}; {} / {}] | {} | {} | {} | {:.1e} | {:.1e} | {:.2f} | {:.1e} | {} | {} / {} | {:.0f} |".format(
            n, " / ".join("{:.4f}".format(v) for v in e["r_lat_mm_at_times"]), " / ".join(pct(v) for v in e["rel_diff_at_times"][1:]),
            "n/a" if e["max_abs_rel_diff_4times"] is None else "{:.1f} %".format(100 * e["max_abs_rel_diff_4times"]),
            e["max_abs_diff_curve_mm"], "never" if e["onset_s"] is None else "{:.2f}".format(e["onset_s"]), st, sl,
            e["deaths_before_any_feed"], e["lateral_deaths_before_lateral_feed"],
            e["skins_lat_away_nx_at_1s"], e["skins_lat_away_nx_at_2s"], e["skins_lat_away_nx_exposed_at_1s"], e["skins_lat_away_nx_exposed_at_2s"],
            e["skins_lat_away_back_at_1s"], e["skins_lat_away_back_at_2s"], e["skins_lat_exposed_noheat_max"], e["skins_lat_away_ns_max"],
            e["cover_v2_withheld"], e["worst_energy_residual"], e["mass_balance_rel"], e["repeats_per_step"], e["max_shortfall_kg"],
            e["skins_created"], e["n_active_final"], e["n_elements"], e["run_time_s"]))
    L += ["", "## The heating and the uncovered faces", "",
          "| run | lateral flux per projected area: mean / min / max (MW/m2) | lowest lateral flux per own area (MW/m2) | power to the walls, max | "
          "smoothed normal from -x on the heated lateral interior: max / power-weighted mean (deg) | folded patch-steps (n.n_s <= 0) | their power under the model's own flip, max share | "
          "uncovered lateral faces: heat-facing / vertical / facing back | tilt of the heat-facing ones: p10 / median / p90 / max (deg), distinct values | "
          "power-weighted tilt of the heated lateral faces after onset: mean / max (deg) |", "|" + "---|" * 10]
    for n in order(out):
        e = out["runs"][n]
        if e["status"] != "ok":
            continue
        tl = e["tilt"]
        tilt = "none uncovered" if not tl.get("n_uncovered") else "{} / {} / {}".format(tl["n_heat_facing"], tl["n_vertical"], tl["n_facing_back"])
        ang = "-" if not tl.get("n_heat_facing") else "{:.1f} / {:.1f} / {:.1f} / {:.1f}; {}".format(
            tl["angle_p10_deg"], tl["angle_median_deg"], tl["angle_p90_deg"], tl["angle_max_deg"], ", ".join("{:g}".format(a) for a in tl["distinct_angles_deg"][:8]))
        L.append("| {} | {:.3f} / {:.3f} / {:.3f} | {:.3f} | {:.4f} | {} | {} | {} | {} | {} | {} |".format(
            n, e["q_lat_proj_MWm2_mean"], e["q_lat_proj_MWm2_min"], e["q_lat_proj_MWm2_max"], e["q_lat_own_MWm2_min"], e["P_wall_share_max"],
            "-" if e["ns_dev_lat_max_deg"] is None else "{:.2f} / {:.2f}".format(e["ns_dev_lat_max_deg"], e["ns_dev_lat_mean_deg_mean"]),
            e["folded_patch_steps"], "-" if e["folded_P_share_literal_max"] is None else "{:.1e}".format(e["folded_P_share_literal_max"]),
            tilt, ang, "-" if e["tilt_lat_power_after_onset_mean_deg"] is None else "{:.1f} / {:.1f}".format(e["tilt_lat_power_after_onset_mean_deg"], e["tilt_lat_power_max_deg"])))
    L += ["", "## Verdicts per run (5 % criterion; against the element model on the same mesh and step)", ""]
    for n in order(out):
        e = out["runs"][n]
        L.append("- `{}`: {}".format(n, verdict(e, out) if e["status"] == "ok" else "The run failed: {}.".format(e.get("error"))))
    rep = out["runs"].get("kuhn_base__dt0.05", {}).get("reproduces_x1")
    if rep:
        L += ["", "## Does the copy with its switches off reproduce X1's base run?", "",
              "`kuhn_base__dt0.05` against `{}`: the largest difference in the lateral recession is {:.3g} mm and in the total recession {:.3g} mm; "
              "the recession, skin-count, active-element, top-temperature and residual series are {}bit for bit the same.".format(
                  rep["file"], rep["max_abs_diff_r_lat_mm"], rep["max_abs_diff_r_total_mm"], "" if rep["bit_identical"] else "NOT ")]
    with open(os.path.join(HERE, "summary.md"), "w") as f:
        f.write("\n".join(L) + "\n")


def plot(recs, fine):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5), squeeze=False)
    ft = [0.0] + fine["series"]["t"]
    fr = [0.0] + fine["series"]["r_total_mm"]
    for j, t in enumerate(LAYERS):
        for i, dt in enumerate(MACRO):
            a = ax[i][j]
            a.plot(ft, fr, "k-", lw=2.5, label="fine box")
            names = ["{}__{}__dt{:g}".format(m, tname(t), dt) for m in MODELS] + ["skin_nowall__{}__dt{:g}".format(tname(t), dt)]
            if dt == 0.05:
                names.append("elem_literal__{}__dt0.05".format(tname(t)))
            for n in names:
                r = recs.get(n)
                if r is None or r.get("status") != "ok":
                    continue
                g = r["group"]
                a.plot([0.0] + r["series"]["t"], [0.0] + r["series"]["r_lat_mm"], lw=1.3, label=g, color=COLOURS[g],
                       ls="--" if g == "skin_v2" else "-")
            a.set_title("layers {:g} mm, macro step {:g} s".format(t * 1e3, dt))
            a.set_xlabel("time (s)")
            a.set_ylabel("lateral recession (mm)")
            a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "x2_box.png"), dpi=120)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--collect", action="store_true")
    a = ap.parse_args()
    sp = specs()
    if a.list:
        print("\n".join(sp))
        return
    if a.probe:
        with open(os.path.join(HERE, "probe.json"), "w") as f:
            json.dump(probe(), f, indent=1)
        return
    if a.collect:
        out, recs, fine = collect()
        write_summary(out)
        plot(recs, fine)
        return
    for name in a.names:
        if name not in sp:
            raise SystemExit("unknown run " + name)
        run(name, sp[name])


if __name__ == "__main__":
    main()
