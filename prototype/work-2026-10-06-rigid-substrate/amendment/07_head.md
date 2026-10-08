## Amendment of 2026-10-06 — the regime test reads the non-rigid layer, and off the closure a film on slurry takes no shear mode

> Part of the rigid-substrate amendment (sub-plan 09's amendment of this date holds the design; facts 78–87 in
> `00-shared-context.md` the measurements and Asha's decisions). The code below is the tested code, as a diff against
> the throwaway copy described there.

**The change: two optional arguments to `SprayModel.evaluate`.**

- `regime_layer` replaces `b_layer` in the regime test only: `deep = has_film & isfinite(delta_m) & (regime > delta_m)`,
  with `regime` the non-rigid layer the body passes — the film, the deep account and everything more than half liquid
  beneath the wall. Girin's first stage asks whether the conjugate layer can form in the melt that is present, and with
  a slurry base it can: the patch is deep melt, and the second stage (kinematic viscosities) then gives the thick branch
  for liquid aluminium, as before. On a patch the test makes deep the shear depth of the Weber number is δ_m (Girin's
  We_s = ρ_l V_s² δ_m / Σ); everywhere else it is `min(delta_m, layer)` as before — the expression is written so that
  without `regime_layer` it is the old one bit for bit.
- `on_slurry` flags the patches whose film rests on slurry. Where δ_m is not finite such a film is `held`: it takes
  neither the thin nor the rarefied mode, since the thin mode needs a rigid wall and no deep-melt mode exists without
  Girin's closure (decision (3) of fact 78's amendment).

**What does not change.** The Rayleigh–Taylor criteria — the reported one and the applied front-surface mode — keep
`b_layer`, the liquid layer; the applied mode is not covered by the hold, so a held film that collects on the cap can
still be shed by it (fact 85; fact 87 (a) asks whether it should). The wave-fits test, the critical Weber number gate, the
release rates and the droplet caps are untouched. With neither argument given every result is the old one bit for bit
(the flights of fact 81).

**Tests.** Three, appended to `tests/test_reentry_model_spray.py` (diff below): a 20 µm film whose slurry ends at half δ_m
stays thin with its Weber number on the film, one whose slurry reaches three times δ_m takes the thick branch with
Girin's Weber number on δ_m, and without `regime_layer` both are thin; the Rayleigh–Taylor flag reads `b_layer` at a
deceleration where the regime layer would pass its depth criterion; and off the closure a film on slurry takes no shear
mode under the Couette closure and in the free-molecular branch, while one on a rigid wall keeps the thin or rarefied mode.

