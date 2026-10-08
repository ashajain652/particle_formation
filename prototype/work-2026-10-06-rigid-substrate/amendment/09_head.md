## Amendment of 2026-10-06 — the thin branch needs a rigid substrate: the non-rigid depth and the regime layer

> Asha's request and decisions of 2026-10-06 (facts 78–87 in `00-shared-context.md` carry the measurements and the
> decisions in full). Measured in a throwaway copy of `prototype/proto3/` with the 24 diff blocks of the amendments of
> 2026-10-02, 2026-10-03 and both of 2026-10-05 applied in date order — every one exactly — and the copy verified before
> any change to reproduce fact 73's flights to the last digit. The code below is the tested code, as a diff against that
> copy. Sub-plan 07's amendment of this date holds the spray's half of the change, 06's what the film sees of it, 10's
> the history columns and the frame field, 13's the flag and 02's `Material.T_rigid`.

**Why.** The spray's regime test (fact 28(a), sub-plan 07) sends a patch to Girin's thin branch — his dominant ablation,
Girin & Kopyt's 1994 mode — when the liquid under it is shallower than the conjugate depth δ_m, on the premise that "the
rigid core still stabilises the disturbances". On an alloy with a 158 K melting range the material under the liquid is
mush, and mush more than half liquid is a slurry that flows (Li et al. 2014), not a rigid core: a film on it is deep melt
and belongs on the thick branch. Coherent mush, less than half liquid, carries load (Chen et al. 2016) and is a rigid
substrate. Fact 54 showed what lies under the film on the flights: the element directly beneath it always has a node in
the mushy range, at a median 897 K — 93 % liquid on the linear law.

**The change, in `body.py`.**

1. `MeltSettings.rigid_substrate: bool = True`, set by `--rigid-substrate on|off` (sub-plan 13).
2. `MeltingBody.nonrigid_depth(T=None, max_crossings=NONRIGID_MAX_CROSSINGS)` — the depth of non-rigid material under each
   patch (fact 79): a line from the patch centre along the inward normal through the elements it actually crosses, the
   P1 field's crossing of `Material.T_rigid` found exactly inside an element, stopping at the first rigid point, the
   active boundary or 64 elements, each element's stretch weighed by φ_e. `NONRIGID_MAX_CROSSINGS = 64`.
3. In `_film_and_spray`, with the switch on: the regime layer `b + nonrigid + deep` replaces `b + molten + deep` in both
   calls of `film.lubrication` (the runoff's flux and the spray's surface velocity and branch flag) and is passed to
   `spray.evaluate` as `regime_layer`; `on_slurry = nonrigid > 0` is passed as well, and the film mass on the windward
   patches it holds off Girin's closure is recorded. `layer = b + molten + deep` still goes to the spray as `b_layer`
   (the Rayleigh–Taylor criteria) and `molten` still defines the wave-fits region.
4. Four entries in `last_melt` and `melt_stats()`: `nonrigid_depth_mean_mm` (mean over the wet windward patches),
   `slurry_thick_fraction` (share of the wet windward patches under Girin's closure made thick only by the slurry),
   `rigid_thin_fraction` (share made thin by the non-rigid depth where the liquid layer said thick — added during the
   measurement, to report fact 58's count being removed) and `slurry_held_mass_kg` (the film held off the closure); all
   nan with the switch off. `last_nonrigid` keeps the step's depth for the surface frame (sub-plan 10).

With the switch off nothing of this runs and every flight is bit-identical to the copy before the change (fact 81).

**What does not change.** The feed, the feed gate of fact 28(b), the deep runoff (it still drains `_molten_chain`, the
fully molten chain), the molten cascade, freeze-back, deaths and every enthalpy booking: the change moves no mass and no
heat itself, it only decides which branch a patch takes. The energy balance stays exact on every flight (facts 82–84).

**Alternatives considered, and why they were not taken.**

- *Any slurry beneath the film triggers the thick branch.* Offered and declined (decision (2)): a slurry skin shallower
  than δ_m over coherent mush leaves rigid material inside the shear's reach, which is Girin's dominant-ablation case.
- *Judge the substrate element by element* (the owner's mean liquid fraction, or `molten_depth`'s march with a T_rigid
  test). Whole elements are 0.67–1.16 mm deep under the patches, all deeper than δ_m (fact 58), so every patch with a
  slurry owner would be thick however thin its slurry: the element-size artefact of fact 58 again. The line through the
  P1 field reads the depth inside the element.
- *Keep the liquid layer for `film.lubrication` and give only the spray Girin's surface velocity.* A patch would carry two
  surface velocities, the film flowing as over a rigid wall while spraying as deep melt; the regime layer is read by both
  instead, and on a patch made thick by slurry the film moves as the top b of the conjugate layer, as it has over fully
  molten elements since 2026-10-05 (sub-plan 06's amendment of this date).
- *Count the slurry in the Rayleigh–Taylor criteria.* They describe a liquid pool, the slurry is 700 to 4 000 times more
  viscous than the liquid (0.87–5.6 Pa s against 1.3 mPa s), and counting it would widen the applied front-surface mode
  over the nose; not asked for.
- *Keep the thin mode off Girin's closure.* Offered and declined (decision (3)). Facts 83 and 85 measure what the strict
  rule does; fact 87 (a) holds the follow-up.

**Two tests changed and one added.** `test_film_spraying_death_and_balances` starts the body at 880 K, which is 82 %
liquid on the linear law, so its held film piles up and the front-surface Rayleigh–Taylor mode sheds droplets above the
size histogram's 10 mm top edge (fact 85); its histogram check now compares the bins with the released droplets inside
their range and requires the rest to lie above it and within R/4. The 2026-10-05 two-backend test
(`test_the_thick_film_flux_matches_the_skfem_backend`) runs with the rule off: its interior at 850 K makes the whole body
slurry, and the backends then part at the spray floor (fact 86); it stays the record of the runoff-flux agreement, and
`test_the_rigid_substrate_matches_the_skfem_backend` checks the backends with the rule on, the interior at 820 K.

