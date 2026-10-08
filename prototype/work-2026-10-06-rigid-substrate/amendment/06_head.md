## Amendment of 2026-10-06 — what the film sees of the rigid substrate (no code change)

> Part of the rigid-substrate amendment (sub-plan 09's amendment of this date holds the design; facts 78–87 in
> `00-shared-context.md` the measurements).

Nothing in `film.py` changes. With `--rigid-substrate on` the body hands `film.lubrication` the regime layer — the film,
the deep account and the non-rigid depth down to the first point at most half liquid — as `b_layer`, in both of its
calls, instead of the film, the deep account and the fully molten depth. Consequences, checked rather than assumed:

- **The branch flag and the surface velocity follow the spray's regime test.** On a patch the slurry makes thick,
  `lubrication` returns Girin's surface velocity V_s = τ δ_m / μ_l, which the spray's thick branch needs for its Weber
  number and growth time; one patch never carries two surface velocities.
- **The film on such a patch moves as the top b of the conjugate layer**, by the 2026-10-05 formula V_s b (1 − b/(2 δ_m))
  plus the pressure-driven G b³/(3 μ_l), instead of as a Couette film over a rigid wall, τ b²/(2 μ_l). The ratio of the
  two shear-driven fluxes is 2 δ_m / b − 1: 59 for a 10 µm film under a 300 µm conjugate depth. This is the treatment the
  film has had over fully molten elements since 2026-10-05; the slurry's own viscosity (700 to 4 000 times the liquid's)
  is not represented, as the deep liquid's is not.
- **The emptying rate stays bounded** by the same argument as on 2026-10-05: the shear part of the flux is at most V_s b,
  so the rate is at most V_s ℓ / A however thin the film. None of the flights of facts 82–84 logged an overflow
  warning or a singular matrix, and every one closed its energy balance to round-off; the rates themselves were not
  probed as fact 74's were.

