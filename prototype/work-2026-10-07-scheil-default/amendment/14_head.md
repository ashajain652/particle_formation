## Amendment of 2026-10-07 — the drivers on Scheil's curve with the rigid substrate, and the switched-step series repeated

> Asha's decisions of 2026-10-07 (facts 88–96 in `00-shared-context.md`; the rigid substrate is the amendment of
> 2026-10-06, facts 78–87). The code below is the tested code, as a diff against the copy of the rigid-substrate
> amendment. These changes stack on this sub-plan's amendments above, which still apply except where an item below
> replaces their wording.

1. **The drivers run the new default.** `melt_verification.py`'s `physics` mode and `melt_sensitivity.py`'s `base` run
   name `--material AA7075_scheil` (they named `AA7075_range`), so every physics row and every sensitivity row is now on
   Scheil's curve with the rigid substrate on. The bookkeeping device (`AA7075`, `--removal instant`, `--runoff off`) is
   untouched: with instant removal there is no film, so the rule never acts, and its rows and thresholds stand. The
   resolved mode (`AA7075`, Girin removal, D₀ and R₀) transports a film and, under Girin's closure, can now take the
   thick branch over slurry — on `AA7075` the slurry band is the 2 K above its 850 K half-liquid point, so the effect is
   expected to be small, which was not measured — and its two rows are to be re-measured with the melting references of
   Task 11.
2. **Two new sensitivity variants.** `norigid` (`--rigid-substrate off`: the regime test reads the fully molten depth,
   as before 2026-10-06) and `range` (`--material AA7075_range`: the linear law, `base`'s material before 2026-10-07).
   Between them they carry the two changes the table's earlier rows were measured without, so each can be read against
   `base` and `seed1` at the default step.
3. **What moves, measured rather than assumed** (facts 82–84 and 89–95). At the default step on the 100 mm flight to 120
   s the rule raises the sprayed mass by 1.2 % on the linear law and 0.6 % on Scheil's, lowers the droplet count by 32 %
   and 26 %, and raises the median radius by number by 4 % and 3 %; on the 50 mm flight, which never has Girin's
   closure, two-thirds of the sprayed mass leaves by the front-surface Rayleigh–Taylor mode and the median radius by
   mass is 7.4 times larger. Every physics and sensitivity row is therefore to be re-run; none was here.
4. **The switched-step series, repeated with the rule on** (facts 91–95), with the 2026-10-05 harness unchanged
   (`dtswitch.py` and `compare.py`, kept with the throwaway copy; fact 77 (b) still asks whether the switch should
   become a model option). It replaces item 4 of the 2026-10-05 amendment's measured values as the series for the
   default model; the 2026-10-05 series stays the record of the linear law with the rule off.
5. **Quote every difference against the scatter at its own step** (item 6 of 2026-10-05, unchanged): this series ran
   seed 1 at 0.05 and 0.0125 s again, and fact 92 gives the spreads.

