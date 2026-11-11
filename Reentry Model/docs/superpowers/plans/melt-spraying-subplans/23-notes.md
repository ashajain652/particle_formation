# Sub-plan 23 (skyscrapers architecture) — carried-forward notes

Sub-plan 23 has not been written yet. This file holds requirements to fold into it when it is.

## 1. Release angle (Register et al. 2020 definition) in the output history files (Asha, 2026-10-29)

**Requirement.** Every droplet release written by the skyscrapers model carries its release angle as defined by
Register et al. (2020), and the per-macro-step history file reports it.

**Definition (Register, P. J., et al. 2020, "Interactions between asteroid fragments during atmospheric entry",
Icarus 337, 113468, section 4.2).** The angle θ is between the trajectory axis of the larger body F1 and the line
joining the centroids of F1 and the released fragment F2.
- **Larger body F1:** the melting sphere.
- **Released fragment F2:** the droplet, at its release position.
- **Measured from the downstream direction:** θ = 0 degrees is directly behind the body and θ = 180 degrees is
  directly ahead (the stagnation point). The paper's zones fix this sense: "near wake" for small θ, "in front" for
  large θ.
- **Converting from the model's angle.** The skyscrapers model measures a surface position's angle from the
  stagnation point (θ_stag). The conversion is θ_Register = 180 degrees - θ_stag. Use the actual body-centroid to
  release-point vector against the body's velocity (air-relative) at release, so that a deformed or receded body
  still gets the right angle.

**Output.**
- **Source rows** (the per-column, per-mode, per-macro-step release records): a column
  `release_angle_register_deg`.
- **History file** (one row per macro step): the mass-weighted mean of the step's release angles, plus their
  minimum and maximum, as `release_angle_register_mean_deg`, `..._min_deg` and `..._max_deg`. Columns are empty
  when nothing is released.

**Optional, if cheap.** The droplet-to-body diameter ratio r = D2/D1 is known at release, so each row can also
record Register's AFFIx final zone z_f (Eq. 5) for that r and θ.
- **The zones:** 3 = trapped in the near wake, 2 = far wake, 1 = carried along the shock, 0 = independent.
- **The two thresholds:**
  - θ_inn(r) = 122 r + 47 degrees, for r <= 0.353;
  - θ_bif(r) = 100 exp(-1.618 r) + 40 degrees.
- **The shock band** has a half-width of 1 to 3 degrees.
- **History summaries:** the step's mass fractions in each zone.
- **Caveat:** AFFIx is fitted for 1/16 <= r <= 1, in continuum hypersonic flow, for equal-density spheres. Droplets
  have r far below 1/16 and are often released in rarefied flow, so z_f is an extrapolation; label it so in the
  column description.

**Why.** The wake-residence work (spec `docs/superpowers/specs/2026-10-29-sphere-wake-residence-design.md` on
branch `wake-residence`; literature notes in
`Particle Wake Evolution/citations/fragment_shock_separation.md`) found that the release angle decides whether a
released body is trapped in its parent's near wake, drifts back in the far wake, rides the shock or escapes. That
matters more than any sideways release speed. Recording the angle now lets the droplet outputs be fed to the
overlap engine (a droplet in its sphere's plume) without re-running flights.

**Tests to add with it.**
- A release at the stagnation point gives 180 degrees.
- A release at the equator relative to the velocity gives 90 degrees.
- The history's mass-weighted mean equals the mean recomputed from the step's source rows.
