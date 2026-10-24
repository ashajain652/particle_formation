# Handoff prompt

Copy everything below the line into a new session.

---

I want to run a two-axis sweep of molten AA7075 droplets released into the wake of a
descending satellite, and report how much each droplet heats up while it decelerates.

## Start by reading these

All paths are relative to the repo root
`/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution`.

- `prototype/work-2026-10-14-particle-settling/particle_trajectory.py` — the single-sphere
  integrator built in a previous session. **Read this first and reuse it; do not rewrite it.**
  Its module docstring explains the composed drag law in full.
- `prototype/work-2026-10-14-particle-settling/terminal_velocity.py` — the steady-state
  settling-speed sweep it was validated against.
- `CLAUDE.md` — repo conventions. Note the interpreter is
  `/Users/ashajain/miniforge3/envs/drama_env/bin/python` (no virtualenv, no `pip install -e`).

These two scripts are working and verified. They are in `prototype/` because that is where
this repo keeps throwaway working copies; nothing is committed yet.

## What the existing script already does

`particle_trajectory.py` integrates one sphere from an initial velocity vector and altitude
down to the ground, in a flat-Earth vertical plane, with the drag coefficient re-evaluated at
every instant. Its drag law spans the whole rarefaction map, which is the hard part:

- **Low-speed branch**: Stokes with the Schiller–Naumann finite-Reynolds factor and the
  Cunningham–Millikan–Davies slip correction. Valid at any Knudsen number, low speed ratio.
- **High-speed branch**: Schaaf & Chambre free-molecule sphere (exact at any speed ratio)
  bridged into the compressible continuum sphere via `reentry_model.aero.SesamErf`.
- The two are blended on the **molecular speed ratio** `s = V/sqrt(2RT)` over 0.05 to 0.5,
  **not** on Knudsen number. This matters: the additive Knudsen bridge in
  `reentry_model/aero.py` is a hypersonic construct and pushes drag the *wrong way* at low
  speed (measured error a factor of 1.94 at Kn 0.1). Do not "simplify" back to it.
- Stiffness is handled by `integrate_hybrid`: full momentum equation for the ballistic phase,
  then the reduced equation `dh/dt = -v_t(h)` once the sphere provably sits on its terminal
  velocity. The switch is verified, not assumed. Without this a 2 um sphere will not finish.
- `check_blend_domain` reports whether a run entered the one uncalibrated corner (Knudsen
  0.1–20 together with speed ratio above 0.05). **Honour this warning in the results** — the
  100 um case at these conditions does trigger it.

Verification already done, all passing: Schaaf–Chambre limits (2.0 at high speed ratio,
4.1906 as the low-speed constant); series/full agreement 1.7e-7 across the switch; the two
branches agree to 5–8 % in the blend window; terminal velocities reproduce `terminal_velocity.py`
to the printed digit; and the hybrid integrator reproduces the full stiff integration exactly
(2.93553e6 s both ways for the 10 um check case).

## What I want built

### 1. The parent state, keyed on altitude

The release velocity vector comes from the parent satellite descent already used by
`sphere_sweep.py`: the file `Generic_Satellite Reentry/output/dmf_output.json`, object
`Main Body`. Reuse `sphere_sweep.load_parent` — do not re-parse the JSON by hand.

`sphere_sweep.interpolate_state` keys on **velocity**, but I need to key on **altitude**.
Write an altitude-keyed companion on the descending branch (`parent.branch`), which has 316
rows and 83 of them between 45 and 95 km, so the band is well sampled.

Convert to the script's planar frame: `Vx = V cos(gamma)`, `Vy = V sin(gamma)`, with `gamma`
the flight path angle in degrees from the `path` column (negative is downward).

Measured parent values for sanity-checking your interpolation:

| altitude | speed | flight path angle |
|---|---|---|
| 90 km | about 7.83 km/s | about -0.37 deg |
| 80 km | about 7.62 km/s | about -0.79 deg |
| 70 km | about 6.91 km/s | about -1.53 deg |
| 60 km | about 4.21 km/s | about -3.57 deg |
| 50 km | about 1.57 km/s | about -10.6 deg |

The speed is nearly flat across the top of the ladder and then collapses over the bottom
third, so the 50 km cases start roughly five times slower than the 90 km cases. Expect the
two ends of the sweep to behave completely differently.

### 2. A thermal model

The droplets are **already molten** — do not model melting, and do not re-run any melting
logic. Use molten AA7075 properties and say explicitly in the report which values you used
and where they came from. Suggested starting values, but check them:

- density about 2400 kg/m3 (the solid value 2813 is in
  `reentry_model/data/materials/AA7075_nomelt.json`; molten is lower)
- liquid specific heat about 1180 J/(kg K)
- emissivity 0.4 is the value in that same material file, but it is a solid-surface value and
  molten aluminium is considerably lower — treat this as an uncertainty and test sensitivity

Add a lumped-capacitance energy balance alongside the momentum equation. The Biot number for
micron spheres is tiny so isothermal is safe, but state that rather than assuming it silently.

- Convective input: free-molecular. In the high-speed-ratio limit the absorbed power is
  `alpha * 0.5 * rho * V^3 * A_frontal`. If you use that limit, **verify** what fraction of
  the total absorbed energy arrives while the speed ratio is above about 3, and report it; if
  a meaningful fraction arrives below that, implement the general Schaaf & Chambre sphere heat
  transfer instead.
- Radiative loss: `epsilon * sigma * (T^4 - T_bg^4)` over the full sphere area.
- Accommodation coefficient: state what you assume.

**Do not model vaporisation**, but flag every case whose temperature crosses the boiling point
of aluminium (about 2740 K), because a prior rough estimate suggested a 10 um sphere released
at 80 km absorbs roughly 9.9e6 J/kg, which is close to the vaporisation enthalpy. I expect a
large share of this matrix to be unphysical for that reason, and I want that stated plainly
rather than hidden behind a temperature number.

### 3. The sweep

- diameters: 2, 6, 10, 14, 18, 50, 100, 150, 250 microns
- release altitudes: 90 km down to 50 km in 2 km steps (21 altitudes)
- 189 cases total

For each case report at minimum: peak temperature and total temperature rise, peak
deceleration, the time and distance over which the droplet loses its release speed, the
downrange distance, the time to reach the ground, and the regime diagnostics
(Knudsen and speed ratio ranges, and the blend-domain warning).

Write the per-case history to CSV and a summary table across the matrix. One run of the
existing script takes a second or two, so the full matrix is cheap.

## How I want it reported

Follow the reporting rules in my global CLAUDE.md: full sentences and complete paragraphs,
every number explained with its units and why it matters, terms spelled out on first use.
Say plainly what was verified, what was assumed, and what remains unknown.

Be aware that long physics write-ups in the previous session repeatedly tripped a safety
classifier and were truncated mid-response. If that happens to you, keep individual messages
shorter and deliver the detail in files rather than in chat.
