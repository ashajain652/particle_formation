"""How often does the skins' freeze-back cap bind, and what would an energy-limited freeze allow instead?

Runs the 100 mm flight with skins through cli.main, instrumenting MeltingBody._freeze_into_skins (the code is not
changed). Per call it records, over the skin patches whose liquid wants to refreeze:
  * want   -- what the rule asks to freeze: (1 - feed_fraction(T0)) of the liquid on the patch, bounded by the liquid;
  * capped -- what the cap (MAX_THICKNESS_FACTOR x the target thickness) lets freeze;
  * energy -- the largest x <= want for which the top node (top cell + liquid), frozen at constant energy, stays at or
              below the bottom of the feed ramp (T_feed - 2 MELT_RAMP): freezing releases latent heat into the node,
              so only part of a thick pile can freeze before the node is back in the melting range;
and where the cap binds (angle from the stagnation point, liquid depth, skin thickness).

    python probe_freeze_cap.py OUTDIR T_MAX [extra cli args...]
"""
import json
import os
import sys

import numpy as np

CODE = "/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/work-2026-10-08-skin/code"
sys.path.insert(0, CODE)
os.chdir(CODE)

from reentry_model import body, cli, material, skin  # noqa: E402

OUT, T_MAX = sys.argv[1], sys.argv[2]
records = []
_melt_step = body.MeltingBody.melt_step
_freeze = body.MeltingBody._freeze_into_skins


def melt_step(self, t, dt, state):
    self._probe_t = t
    return _melt_step(self, t, dt, state)


def energy_limited(self, p, r, want):
    """Largest x in [0, want] per patch with the frozen node at or below the bottom of the feed ramp."""
    mat, sk = self.material, self.skins
    L = (self.m_f + self.m_d)[p]
    m0 = sk.m[r, 0]
    E = sk.H[r, 0] + L * mat.enthalpy_liquid(sk.T[r, 0])
    T_low = mat.T_feed - 2.0 * material.MELT_RAMP

    def T_after(x):
        mm, ll = m0 + x, L - x
        tot = mm + ll
        return mat.temperature_from_enthalpy_mixed(E / tot, np.where(tot > 0, ll / tot, 0.0))

    ok_all = T_after(want) <= T_low
    lo, hi = np.zeros_like(want), want.copy()
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        good = T_after(mid) <= T_low
        lo, hi = np.where(good, mid, lo), np.where(good, hi, mid)
    return np.where(ok_all, want, lo)


def freeze(self, solid):
    sk = self.skins
    p = np.flatnonzero(solid > 0.0)
    if p.size:
        r = self.skin_of_patch[p]
        L = (self.m_f + self.m_d)[p]
        want = np.minimum(solid[p], L)
        room = np.maximum(self.material.rho * sk.area[r] * skin.MAX_THICKNESS_FACTOR * sk.settings.thickness
                          - sk.m[r].sum(axis=1), 0.0)
        capped = np.minimum(want, room)
        energy = energy_limited(self, p, r, want)
        binds = want > room * (1.0 + 1e-12)
        depth_liquid = L / (self.liquid.rho * self.surface.areas[p])
        thick = sk.thickness()[r]
        theta = np.degrees(self.theta[p])
        records.append({
            "t": float(getattr(self, "_probe_t", np.nan)), "n_skins": int(sk.n), "n_freezing": int(p.size),
            "n_binding": int(binds.sum()), "want_g": float(want.sum() * 1e3), "capped_g": float(capped.sum() * 1e3),
            "held_by_cap_g": float((want - capped).sum() * 1e3), "energy_g": float(energy.sum() * 1e3),
            "n_energy_limits": int((energy < want * (1 - 1e-9)).sum()),
            "n_energy_below_cap": int((energy < capped * (1 - 1e-9)).sum()),
            "binding_theta_deg": theta[binds].round(1).tolist(),
            "binding_liquid_depth_mm": (depth_liquid[binds] * 1e3).round(3).tolist(),
            "binding_skin_thickness_mm": (thick[binds] * 1e3).round(3).tolist(),
            "binding_want_over_room": (want[binds] / np.maximum(room[binds], 1e-30)).round(2).tolist(),
            "freezing_theta_deg_median": float(np.median(theta)),
            "skin_thickness_max_mm": float(sk.thickness().max() * 1e3),
            "frozen_share_of_liquid_on_patches": float(want.sum() / max(L.sum(), 1e-30)),
        })
    return _freeze(self, solid)


body.MeltingBody.melt_step = melt_step
body.MeltingBody._freeze_into_skins = freeze

argv = ["run", "--diameter", "100", "--velocity", "7.5", "--altitude", "77.500133", "--flight-path-angle", "-0.959331",
        "--atmosphere", "us76", "--thermal", "fem", "--melt", "on", "--heating", "physics", "--t-max", T_MAX,
        "--surface-model", "skin", "--outdir", OUT] + sys.argv[3:]
rc = cli.main(argv)
json.dump(records, open(os.path.join(OUT, "freeze_probe.json"), "w"))
print("exit", rc, "calls", len(records))
