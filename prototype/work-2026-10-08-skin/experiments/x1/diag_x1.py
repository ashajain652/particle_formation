"""Experiment X1, why the switches do what they do (diagnostic, no change to the model). Reruns one x1_box run with
extra per-step records:

  * the heated lateral interior (patches with n_x < -0.1 at least one coarse cell from every side plane): how many
    patches, how many carry a skin, and how many bare ones are above the solidus and in the exposed set (the faces
    cover_exposed refuses a created skin);
  * skins removed empty at a hand-over (`_hand_over_skins` keeps only skins with mass);
  * the draw's energy booking at the start of the melt step, sum of (m_drawn h_e - H_drawn) over elements [J] -- what
    the elements' nodes are debited because the drawn metal entered the skins at h(T_base) rather than at the element's
    own mean enthalpy -- split into patch-owning and buried elements; and the pending (deferred) load after the step.

    python diag_x1.py both__dt0.05 [more run names]     -> runs/diag/diag__<name>.json and a short print
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x1_box as xb  # noqa: E402  (puts the experiment copy first and patches box_ablation's make_body)

import numpy as np  # noqa: E402

bx = xb.bx


def instrument(b, spec, rows):
    mat = b.material
    melt_step, hand = b.melt_step, b._hand_over_skins
    state = {"empty": 0}

    def melt_step_spy(t, dt, st):
        T = b.solver.temperature()
        h_e = mat.enthalpy(T)[b.mesh.tets].mean(axis=1)
        book = b._draw_m * h_e - b._draw_H
        owner = b.owner_area > 0.0
        state["draw_owner_J"] = float(book[owner].sum())
        state["draw_buried_J"] = float(book[~owner].sum())
        state["drawn_buried_kg"] = float(b._draw_m[~owner].sum())
        return melt_step(t, dt, st)

    def hand_spy2(old_surface):
        sk = b.skins
        n_lost = int((b.patch_of_face[sk.face_id] < 0).sum())
        n_empty_kept = int(((b.patch_of_face[sk.face_id] >= 0) & (sk.m.sum(axis=1) <= 0.0)).sum())
        state["empty"] += n_empty_kept
        hand(old_surface)

    b.melt_step = melt_step_spy
    b._hand_over_skins = hand_spy2
    inner = b.advance

    def advance(t, dt, loads, st=None):
        state.update(draw_owner_J=0.0, draw_buried_J=0.0, drawn_buried_kg=0.0)
        heat_in = float((loads.q_conv * b.surface.areas).sum() * dt)
        inner(t, dt, loads, st)
        _, _, _, lateral = bx.classify(b, spec["box"])
        s = b.surface
        heated = lateral & (s.normals[:, 0] < -0.1)
        Tf = s.facet_mean(b.solver.temperature())
        exp = np.isin(s.face_ids, np.fromiter(b.exposed_faces, dtype=np.int64, count=len(b.exposed_faces)))
        bare = b.skin_of_patch < 0
        rows.append({"t": round(t, 4), "heated_patches": int(heated.sum()), "heated_with_skin": int((heated & ~bare).sum()),
                     "heated_bare_hot_exposed": int((heated & bare & exp & (Tf >= mat.T_solidus)).sum()),
                     "heated_bare_hot": int((heated & bare & (Tf >= mat.T_solidus)).sum()),
                     "skins_removed_empty_cum": state["empty"],
                     "draw_booking_owner_J": state["draw_owner_J"], "draw_booking_buried_J": state["draw_buried_J"],
                     "drawn_from_buried_g": state["drawn_buried_kg"] * 1e3,
                     "pending_load_J": float(b.pending_load.sum()),
                     "heat_in_J": heat_in})

    b.advance = advance


def main():
    for name in sys.argv[1:]:
        spec = xb.specs()[name]
        rows = []
        orig = bx.make_body

        def make_body(sp):
            b, mat = orig(sp)
            instrument(b, sp, rows)
            return b, mat

        bx.make_body = make_body
        try:
            r = bx.run_one(name, spec, log=lambda s: None)
        finally:
            bx.make_body = orig
        out = {"name": name, "rows": rows, "r_lat_mm": r["series"]["r_lat_mm"], "worst_energy_residual": r["worst_energy_residual"]}
        with open(os.path.join(xb.RUNS, "diag", "diag__{}.json".format(name)), "w") as f:
            json.dump(out, f, indent=1)
        print("==", name, "(r_lat at 2 s = {:.4f} mm; the same as the measured run: {})".format(
            r["series"]["r_lat_mm"][-1], r["series"]["r_lat_mm"] == json.load(open(os.path.join(xb.RUNS, name + ".json")))["series"]["r_lat_mm"]))
        for w in rows:
            print("  t={t:.2f} heated {heated_with_skin}/{heated_patches} with skin, bare & hot & exposed {heated_bare_hot_exposed}, "
                  "removed empty (cum) {skins_removed_empty_cum}, draw booking owner {draw_booking_owner_J:+.2f} J buried "
                  "{draw_booking_buried_J:+.2f} J ({drawn_from_buried_g:.4f} g), pending {pending_load_J:+.2f} J".format(**w))


if __name__ == "__main__":
    main()
