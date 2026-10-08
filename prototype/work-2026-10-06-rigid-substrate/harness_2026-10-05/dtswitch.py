"""A run whose macro step shortens once the sphere enters the continuum regime, instrumented to show where the melt is
made and where it is sprayed. No package file is edited: the harness wraps methods from the outside, like
casc/measure3.py, so a baseline and a switched run are measured by the same code.

    dtswitch.py <copy dir> <record.json> <kn switch> <dt fine> <frame interval s> -- <cli run arguments>

The step: the CLI's --dt until the body Knudsen number (the trajectory's AeroState.kn, lambda_inf / D, the model's
own `knudsen` column) first falls below <kn switch>; from the next step on, <dt fine>, latched. A switch of 0 never
fires (the baseline). Frames are written every <frame interval> seconds of flight time rather than every n steps, so
runs with different steps have frames at the same times.

Recorded per melt step, on the surface the spray stage sees (deaths come later in the step), by the polar angle psi of
each patch's centroid about the flight axis through the mass centre (0 at the front, 90 at the equator):
  * avail   -- the film on the patch when the spray stage begins, before the runoff transport: the step's melt feed,
               the deep liquid that surfaced, and whatever the previous step left (cascade and death hand-overs);
  * arrive / leave -- the runoff transport's net gain / net loss of film per patch over the step;
  * rel     -- the mass the spray stage released from the patch;
  * the same arrivals and releases on the equatorial ring (windward patches bordering a leeward one and lying at least
    0.9 of the transverse radius from the axis, the definition of casc/measure3.py);
  * the deep transport's net arrivals (liquid below the conjugate depth moving along the surface, never sprayed).
"""
import json
import os
import sys
import time

copy_dir, out_json = sys.argv[1], sys.argv[2]
KN_SWITCH, DT_FINE, FRAME_DT = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
argv = sys.argv[sys.argv.index("--") + 1:]
sys.path.insert(0, copy_dir)
import numpy as np  # noqa: E402

from reentry_model import body as body_mod, cli, coupled, film  # noqa: E402

EDGES = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 85.0, 90.0, 180.0]
REC = {"argv": argv, "copy": copy_dir, "kn_switch": KN_SWITCH, "dt_fine": DT_FINE, "frame_dt": FRAME_DT,
       "edges_deg": EDGES, "steps": [], "switched_at_s": None}
STATE = {"run": None, "body": None, "phase": None, "film_io": None, "deep_moved": None, "deep_bins": None,
         "next_frame": 0.0, "t0": time.perf_counter()}


def psi_of(b):
    x = b.surface.centroids - b.mass_centre
    r = np.linalg.norm(x, axis=1)
    c = np.clip((x @ b.v_hat) / np.where(r > 0.0, r, 1.0), -1.0, 1.0)
    return np.degrees(np.arccos(c)), x


def ring_of(b, x):
    s = b.surface
    _, i, j = s.edges()
    w = b.windward
    mixed = w[i] != w[j]
    row = np.zeros(s.n_patches, dtype=bool)
    row[i[mixed & w[i]]] = True
    row[j[mixed & w[j]]] = True
    rad = np.sqrt(np.maximum((x * x).sum(axis=1) - (x @ b.v_hat) ** 2, 0.0))
    return row & (rad >= 0.9 * b.transverse_radius)


def binned(psi, values):
    k = np.clip(np.digitize(psi, EDGES) - 1, 0, len(EDGES) - 2)
    return np.bincount(k, values, len(EDGES) - 1).tolist()


# -- the coupled run: capture it, and write frames by flight time --------------------------------------------------------
orig_init = coupled.CoupledRun.__init__


def init(self, *a, **k):
    orig_init(self, *a, **k)
    STATE["run"] = self


def frame(self, step, t, loads):
    s = self.settings
    if not s.frames_every or t + 1e-9 < STATE["next_frame"]:
        return
    k = len(s.frames)
    coupled.write_vtk_frame(s.output_dir, k, self.body, loads)
    s.frames.append((t, k))
    while STATE["next_frame"] <= t + 1e-9:
        STATE["next_frame"] += FRAME_DT


# -- the step switch -----------------------------------------------------------------------------------------------------
orig_advance = body_mod.MeltingBody.advance


def advance(self, t, dt, loads, state=None):
    STATE["body"] = self
    orig_advance(self, t, dt, loads, state)
    if REC["steps"] and REC["steps"][-1].get("t") == t:
        REC["steps"][-1]["wall_s"] = time.perf_counter() - STATE["t0"]
    run = STATE["run"]
    if (REC["switched_at_s"] is None and KN_SWITCH > 0.0 and state is not None and run is not None
            and np.isfinite(state.kn) and state.kn < KN_SWITCH):
        run.settings.dt = DT_FINE
        REC["switched_at_s"] = float(t)
        REC["switched_kn"] = float(state.kn)
        print("dtswitch: Kn {:.5f} < {} at t = {:.3f} s -> dt = {} s".format(state.kn, KN_SWITCH, t, DT_FINE), flush=True)


# -- the instrumentation -------------------------------------------------------------------------------------------------
orig_transport = film.Runoff.transport


def transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps=4):
    m_in = np.array(m_f, dtype=float)
    out = orig_transport(self, m_f, q_of_thickness, t_hat, rho_l, areas, dt, substeps)
    d = out[0] - m_in
    if STATE["phase"] == "film":
        STATE["film_io"] = (m_in, out[0].copy(), out[1])
    elif STATE["phase"] == "deep":
        STATE["deep_moved"] = (STATE["deep_moved"] or 0.0) + float(np.maximum(d, 0.0).sum())
        b = STATE["body"]
        if b is not None and len(d) == b.surface.n_patches:
            psi, _ = psi_of(b)
            prev = STATE["deep_bins"] or [0.0] * (len(EDGES) - 1)
            STATE["deep_bins"] = [p + q for p, q in zip(prev, binned(psi, np.maximum(d, 0.0)))]
    return out


orig_deep = getattr(body_mod.MeltingBody, "_deep_runoff", None)


def deep_runoff(self, *a, **k):
    STATE["phase"] = "deep"
    try:
        return orig_deep(self, *a, **k)
    finally:
        STATE["phase"] = None


orig_fs = body_mod.MeltingBody._film_and_spray


def film_and_spray(self, t, dt, state, h_p, flow=None, delta_m=None):
    STATE["body"] = self
    m_entry = self.m_f.copy()
    STATE["phase"], STATE["film_io"] = "film", None
    try:
        released = orig_fs(self, t, dt, state, h_p, flow, delta_m)
    finally:
        STATE["phase"] = None
    n = self.surface.n_patches
    rec = {"t": float(t), "dt": float(dt), "released_kg": float(released)}
    if n:
        psi, x = psi_of(self)
        if STATE["film_io"] is not None and len(STATE["film_io"][0]) == n:
            m0, m1, n_sub = STATE["film_io"]
        else:
            m0, m1, n_sub = m_entry, m_entry, 0
        d = m1 - m0
        res = self.last_spray
        rel = np.asarray(res.dm, dtype=float) if (res is not None and len(res.dm) == n) else np.zeros(n)
        ring = ring_of(self, x)
        rec.update({"avail": binned(psi, m0), "arrive": binned(psi, np.maximum(d, 0.0)),
                    "leave": binned(psi, np.maximum(-d, 0.0)), "rel": binned(psi, rel),
                    "ring_arrive_kg": float(np.maximum(d, 0.0)[ring].sum()), "ring_rel_kg": float(rel[ring].sum()),
                    "ring_avail_kg": float(m0[ring].sum()), "n_ring": int(ring.sum()), "n_sub": int(n_sub),
                    "film_end_kg": float(self.m_f.sum()), "net_d_kg": float(d.sum())})
    rec["deep_arrive_kg"] = float(STATE["deep_moved"] or 0.0)
    rec["deep_arrive_bins"] = STATE["deep_bins"]
    STATE["deep_moved"], STATE["deep_bins"] = None, None
    REC["steps"].append(rec)
    return released


coupled.CoupledRun.__init__ = init
coupled.CoupledRun._frame = frame
body_mod.MeltingBody.advance = advance
film.Runoff.transport = transport
if orig_deep is not None:
    body_mod.MeltingBody._deep_runoff = deep_runoff
body_mod.MeltingBody._film_and_spray = film_and_spray

started = time.perf_counter()
code = cli.main(argv)
REC["exit_code"] = code
REC["wall_s"] = time.perf_counter() - started
b = STATE["body"]
if b is not None:
    REC["final"] = {"sprayed_kg": b.sprayed_mass, "runoff_kg": b.runoff_mass, "mass_kg": b.mass(0.0),
                    "film_kg": float(np.sum(b.m_f)), "balance": b.energy_balance_residual(), "frozen_kg": b.frozen_mass,
                    "n_released": b.n_released, "deep_runoff_kg": float(getattr(b, "deep_runoff_mass", 0.0)),
                    "cascade_kg": float(getattr(b, "cascade_mass", 0.0))}
os.makedirs(os.path.dirname(os.path.abspath(out_json)), exist_ok=True)
with open(out_json, "w") as fh:
    json.dump(REC, fh)
print("dtswitch: exit {} in {:.0f} s -> {}".format(code, REC["wall_s"], out_json), flush=True)
sys.exit(code)
