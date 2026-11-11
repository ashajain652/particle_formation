"""Why do the skins run short of metal in the box (q6p variant, dev preset, 0.05 s)? Phase 1 evidence, no fix.

Reruns box_ablation's q6p__skin_dev_dt0.05 step by step and records, per step, over the skins on the lateral interior
(patches at least one coarse cell from every side plane, any orientation):
  * counts: skins, short (< s/2), empty (< 1 um); deaths; skins created (and how many on front-facing patches, n_x < -0.1);
    mass drawn, fed (liquid removed) and handed over;
  * for the short skins, the reason their element cannot supply them:
      dead_wait  -- the owner element is at or below the death level (it has nothing above PHI_DEATH) but is still alive
      shared     -- the owner holds metal above PHI_DEATH, but it also owns patches without a skin, whose area share of the
                    owner's metal is reserved for them
      exhausted  -- the owner has given everything above PHI_DEATH to its skins this step (avail ~ 0, no skinless patches)
      other
  * for skins created this step: on newly exposed patches (not on the surface the step before) or on old ones.

    python diag_stall.py [run name, default skin_dev_dt0.05] [variant, default q6p]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import box_ablation as bx  # noqa: E402  (imports the model from the frozen measurement copy)
from reentry_model import body, heating  # noqa: E402

name = sys.argv[1] if len(sys.argv) > 1 else "skin_dev_dt0.05"
tag = sys.argv[2] if len(sys.argv) > 2 else "q6p"
v = bx.VARIANTS[tag]
spec = bx.run_specs(tag=tag, q=v["q"], rule=v["rule"])[name]
b, mat = bx.make_body(spec)
lx, ly, lz, h = spec["box"]
dt, q_flux = spec["dt"], spec["q"]
A_front0 = ly * lz
s_target = b.skins.settings.thickness if b.skins is not None else float("nan")
rows = []
prev_faces = set(b.surface.face_ids.tolist())
prev_created = prev_draw = prev_hand = 0.0
prev_active = b.mesh.n_active
cum_lat = 0.0
A_lat0 = (ly - 2 * h) * (lz - 2 * h)
for k in range(int(round(spec["t_end"] / dt))):
    t = (k + 1) * dt
    n_x = b.surface.normals[:, 0]
    w = np.maximum(0.0, -n_x)
    q = q_flux * A_front0 * w / float((w * b.surface.areas).sum())
    b.advance(t, dt, heating.HeatingResult(q, q_flux, 0.0, 0.0, 0.0))
    s = b.surface
    front, interior, strict, lateral = bx.classify(b, spec["box"])
    L = b.m_f + b.m_d
    fed_lat = float(L[lateral].sum())
    cum_lat += fed_lat / (A_lat0 * bx.RHO_METAL) * 1e3
    b.removed_mass += float(L.sum())
    b.removed_enthalpy += float((L * b.film_enthalpy()).sum())
    b.m_f[:] = 0.0
    b.m_d[:] = 0.0
    b.solver.set_film_mass(b._film_nodal())
    sk = b.skins
    rec = {"t": round(t, 4), "r_lat_mm": cum_lat, "deaths": int(prev_active - b.mesh.n_active)}
    prev_active = b.mesh.n_active
    faces_now = set(s.face_ids.tolist())
    new_faces = faces_now - prev_faces
    prev_faces = faces_now
    rec["created"] = int(b.skins_created - prev_created)
    prev_created = b.skins_created
    rec["drawn_g"] = (b.skin_drawn_mass - prev_draw) * 1e3
    prev_draw = b.skin_drawn_mass
    rec["handed_g"] = (b.skin_handover_mass - prev_hand) * 1e3
    prev_hand = b.skin_handover_mass
    rec["fed_lat_g"] = fed_lat * 1e3
    if sk is not None and sk.n:
        p = b.patch_of_skin
        th = sk.thickness()
        lat = lateral[p]
        short = lat & (th < 0.5 * s_target)
        e = s.owner[p]
        n_patches_of = np.bincount(s.owner, minlength=b.mesh.n_elements)
        n_skins_of = np.bincount(e, minlength=b.mesh.n_elements)
        skinless = n_patches_of[e] > n_skins_of[e]
        avail = b._skin_available(p)
        dead_level = body.PHI_DEATH * (1.0 + body.DEATH_TOLERANCE)
        dead_wait = short & (b.phi[e] <= dead_level)
        shared = short & ~dead_wait & skinless
        exhausted = short & ~dead_wait & ~skinless & (avail < 1e-12)
        other = short & ~dead_wait & ~shared & ~exhausted
        rec.update({"skins_lat": int(lat.sum()), "short_lat": int(short.sum()),
                    "empty_lat": int((lat & (th < 1e-6)).sum()),
                    "why_dead_wait": int(dead_wait.sum()), "why_shared": int(shared.sum()),
                    "why_exhausted": int(exhausted.sum()), "why_other": int(other.sum()),
                    "owner_phi_short_median": float(np.median(b.phi[e[short]])) if short.any() else float("nan"),
                    "skins_total": int(sk.n),
                    "skins_on_new_faces": int(np.isin(sk.face_id, list(new_faces)).sum()) if new_faces else 0,
                    "new_faces": len(new_faces),
                    "new_faces_front_facing": int((np.isin(s.face_ids, list(new_faces)) & (s.normals[:, 0] < -0.1)).sum()) if new_faces else 0,
                    "skins_front_facing_lat": int((lat & (s.normals[p, 0] < -0.1)).sum()),
                    "skins_not_front_facing_lat": int((lat & (s.normals[p, 0] >= -0.1)).sum()),
                    "thick_median_um_lat": float(np.median(th[lat]) * 1e6) if lat.any() else float("nan"),
                    "shortfall_g": b.last_skin_shortfall * 1e3})
    rows.append(rec)
    if k % 2 == 1:
        print(json.dumps({kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in rec.items()}), flush=True)
json.dump(rows, open(os.path.join(HERE, "runs", "diag_stall__{}__{}.json".format(tag, name)), "w"), indent=1)
