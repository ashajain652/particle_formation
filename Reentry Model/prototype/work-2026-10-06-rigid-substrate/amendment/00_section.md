## Amendment of 2026-10-06 — the thin branch needs a rigid substrate (facts 78–87)

Asha's request of 2026-10-06: Girin's thin branch — his dominant ablation, the Girin & Kopyt (1994) mode, the case where
"the rigid core still stabilises the disturbances" — is valid only where the film rests on a rigid surface. A film on
slurry is deep melt and takes the thick branch; a film on coherent (semi-solid) mush or on solid AA7075 keeps the thin
one. Her decisions while the design was settled, the same day: (1) the boundary between coherent mush and slurry is
50 % liquid, the boundary Step 4 and the large-fragment design already use (Chen et al. 2016's semi-solid strength law
ends there and Li et al. 2014's slurry viscosity data begin there); (2) the regime test's melt layer runs from the top of
the film down to the first point at that boundary, and the thick branch fires where that layer is deeper than δ_m — a
slurry skin shallower than δ_m over rigid material keeps the thin branch (chosen over "any slurry beneath the film
triggers the thick branch"); (3) where there is no δ_m (off Girin's closure) a film whose base is slurry does not spray
by the shear modes at all — the strict reading, chosen over keeping the thin mode there; (4) the design as presented and
approved: lubrication's branch flag and surface velocity read the same layer as the spray's regime test, the
Rayleigh–Taylor criteria and the wave-fits region keep the liquid layer, the front-surface Rayleigh–Taylor mode is not
covered by (3), `--rigid-substrate on|off` defaults to on, and the change is recorded in history columns and a frame
field. It amends sub-plans 02 (`Material.T_rigid`), 06 (no code change), 07 (`spray.evaluate`), 09
(`MeltingBody.nonrigid_depth` and the melt step; one melting test and one two-backend test changed, one two-backend
test added), 10 (four history columns and a surface-frame field) and 13 (the flag). Sub-plans 14 and 15 are not amended
(fact 87 (d)). Facts 1–77 stand except where these say otherwise. **With the rule on, fact 28(a)'s decision that the
contiguous liquid depth decides the branch is replaced by fact 79's non-rigid depth, so fact 58's whole-element count
and double count no longer enter the branch test; decision 77(e) — whether the merged branch should spray — is answered
for films on slurry by (3), and facts 83 and 85 say what that does.** Measured in a throwaway copy
(`prototype/work-2026-10-06-rigid-substrate/`, git-ignored) of `prototype/proto3/` — unchanged since 2026-10-03, its 114
files plus six cached meshes — with the 24 diff blocks of the amendments of 2026-10-02 (9), 2026-10-03 (8), 2026-10-05
seeding (3) and 2026-10-05 runoff flux (4) applied in date order, every one exactly (no fuzz, no offset), and sub-plan
02's material code of 2026-09-28 (the Scheil variant and `AA7075-empiricaldata`) added for the Scheil runs; the
generator in sub-plan 02 rewrote `AA7075.json` and `AA7075_range.json` byte-identical to the copy's (`cmp`). Before any
change the copy reproduced fact 73's flights to the last digit (fact 81), so everything facts 49, 54, 62 and 69 list
about the copy holds: `AA7075_range` unless stated, US76, physics heating, the dense band, `PHI_DEATH = 0.05`, the deep
runoff and the molten cascade on, every run seeded with the default 12345. The six measurement runs were made from the
amended copy, never edited while they ran, all six at once, on mains power under `caffeinate -i`, the power source
logged at the start and end of each. (The baseline 100 mm run's log records battery power at its end — a short
interruption of the mains; its results match fact 73 to the last digit, as a seeded run's must.)

78. **The half-liquid temperature.** `Material.T_rigid` is the temperature at which the material is
    `RIGID_LIQUID_FRACTION = 0.5` liquid, found by bisection on the material's own liquid-fraction law (monotonic):
    829.0 K for `AA7075_range` (750 + 0.5 × 158), 895.10 K for `AA7075_scheil` (Scheil's formula gives 895.11 K; its
    1 K table, 895.10 K), 850.0 K for the single-temperature `AA7075` (the middle of its ±2 K ramp), and infinite for a
    material that does not melt. The slurry band between T_rigid and the liquidus is therefore 79 K wide on the linear
    law and 12.9 K on Scheil's, and that difference is the one the flights depend on most (facts 82 and 84).

79. **The non-rigid depth, and what reads it.** `MeltingBody.nonrigid_depth` follows a line from each patch's centre
    along its inward normal through the tetrahedra it actually crosses, leaving each by the face its barycentric
    coordinates reach first. The P1 temperature is linear along the line inside an element, so the point where it falls
    to T_rigid is found exactly rather than counted in whole elements. Only material continuous with the wall counts: the
    line stops at the first point at or below T_rigid, at the active mesh's boundary, or after
    `NONRIGID_MAX_CROSSINGS = 64` elements. Each element's stretch counts φ_e of its length, because what the element has
    fed to the film is already in the film account added on top — fact 58's double count, removed by an assumption:
    that the fed part came evenly from along the line (fact 87 (f)). Tested exact to 1e-9 relative on a linear field at
    0.4 mm (inside the second prism layer) and 3 mm (several elements down), halved to 1e-12 by φ = 0.5, and stopped at
    the first rigid point when the field is hot again deeper. The regime layer is the film depth plus the deep account's
    plus the non-rigid depth. It replaces `molten_depth` in four places only: `film.lubrication`'s branch flag and its
    surface velocity (both calls, so the runoff's flux and the spray's Weber number see the same branch), the spray's
    regime test, and the shear depth of Girin's Weber number on the patches the test makes deep (δ_m there, as on every
    thick patch). `molten_depth` still feeds the Rayleigh–Taylor criteria, the molten region the wave-fits test measures
    against, the deep runoff's chain, the molten cascade and the `molten_depth_*` columns. Off Girin's closure a film
    whose wall is above T_rigid (non-rigid depth > 0, `on_slurry`) takes neither the thin nor the rarefied mode; the
    front-surface Rayleigh–Taylor mode is not affected (decision (4)). With the switch off every one of these reads the
    liquid layer as before, bit for bit (fact 81).

80. **Tests.** Added: `test_the_rigid_temperature_is_where_half_the_material_is_liquid` (02); three spray tests — a thin
    film on slurry deeper than δ_m takes the thick branch with Girin's Weber number on δ_m while one whose slurry ends
    within δ_m stays thin, the Rayleigh–Taylor criteria keep the liquid layer, and off the closure a film on slurry takes
    no shear mode in the Couette and free-molecular branches alike (07); five melting tests — the three depth tests of
    fact 79, a 20 µm film over a millimetre of slurry at 70 km thick on every wet windward patch with the rule on and
    thin with it off, and at 30 s (no closure) the film held over slurry and sprayed over a rigid wall (09); and a
    two-backend test with the rule on (09). Extended: the CLI's run name, columns, settings and bad-argument tests (13)
    and the coupled run's frame fields (10). Changed, both because their scenarios make the whole body slurry under the
    linear law: (a) `test_film_spraying_death_and_balances` starts the body at 880 K, 82 % liquid, so its held film
    piles up and the front-surface Rayleigh–Taylor mode sheds droplets above the size histogram's 10 mm top edge (fact
    85); its check that the histogram holds every droplet released now compares the histogram with the released droplets
    inside its range and requires the rest to lie above it and below R/4 — the histogram drops out-of-range radii by
    design; (b) the 2026-10-05 two-backend test runs with the rule off, because with its interior at 850 K the two
    backends stop agreeing at the spray floor (fact 86) and it is the record of the runoff-flux agreement. Unit tier
    (`not drama and not reference`): 251 passed, 1 skipped, and the same 2 failures and 3 errors as the copy before the
    change (236 passed) — all five need Task 11's melting SESAM references, which the prototype does not have. The 15
    tests added are these 9, sub-plan 02's five Scheil and empirical-data tests, and one new bad-argument case. FEniCSx
    tier: 12 passed (11 before).

81. **Reproduction, and the switch off bit for bit.** The copy before the change reproduced fact 73: the 100 mm flight
    to 120 s sprays 1.0585 kg, leaves 0.4136 kg, writes 197 075 source-table rows, moves 90.9 g by the deep runoff and
    158.0 g by the cascade, re-solidifies 6.62 g and closes its energy balance to −8.5e-11; the 50 mm flight demises at
    199.0 s and 70.74 km with 399 history rows and 21 010 source rows. With `--rigid-substrate off` the amended code
    reproduces both bit for bit: all 241 and 399 rows of all 86 history columns the runs share, every result field but
    the run time, and all 22 source-table columns.

82. **The 100 mm flight to 120 s (`AA7075_range`), the rule on against off**, each change read against fact 65's
    scatter between four seeds: sprayed mass 1.0585 to 1.0708 kg (+1.2 %, scatter 0.13 %); the body at 120 s 0.4136 to
    0.4012 kg (−3.0 %, scatter 0.33 %); droplets 17.89 to 12.10 million (−32 %, scatter 2.0 %); median radius by number
    179.7 to 186.5 µm (+3.8 %, scatter 0.29 %) and by mass 190.6 to 199.5 µm (+4.6 %, scatter 0.16 %); re-solidified
    6.62 to 4.27 g (−36 %, scatter 1.4 %); largest droplet 5.3 to 5.2 mm, none above 10 mm. By branch, as shares of the
    sprayed mass: thick 85.1 to 88.2 % (901 to 944 g), thin 10.4 to 2.5 % (110 to 27 g), front-surface Rayleigh–Taylor
    4.5 to 9.4 % (47.4 to 100.1 g). The thin branch's droplets fall from 6.61 to 0.92 million, which is 5.69 million of
    the 5.79 million fewer droplets. Before Girin's closure (to 49.5 s) the flight sprays 87.3 g, all by the thin mode,
    with the rule off, and 74.8 g with it on — 16.4 g by the thin mode on walls below T_rigid (the thin mode's sprayed
    mass sits at a mass-weighted 70° from the stagnation point, against 52° with the rule off) and 58.5 g by the
    Rayleigh–Taylor mode; spraying begins at 28.5 s instead of 25.5 s, and the film held from the shear modes peaks at
    16.0 g at 48.5 s, one step before the closure begins (the film on the body peaks at 15.5 g against 6.4 g). Under the
    closure, as medians over its 142 steps: the thick share of the wet windward patches rises from 27.4 to 62.8 %; within
    the run with the rule on, 23.3 % of those patches are thick only because of slurry and 11.6 % are thin by the
    non-rigid depth where the whole-element liquid layer said thick (fact 58's count, removed); and the non-rigid depth
    under them is a median 8.1 mm, against a molten depth of 1.15 mm and a δ_m of 290 µm. Medians of fractions taken
    within one run do not add up to the difference between two runs, whose wet patches differ. The energy balance stays
    exact (−1.2e-10 of the absorbed heat), Newton takes 3.00 iterations per step (3.02), the peak surface temperature is
    975.8 K (981.1 K), the deep runoff moves 91.1 g (90.9 g) and the cascade feeds 156.9 g (158.0 g). The run took
    990 s against 973 s (+1.7 %), but all six runs ran at once, so the cost is not measured cleanly (fact 76).

83. **The 50 mm flight, which never has Girin's closure.** Decision (3) acts on the whole flight. Sprayed mass 0.1770 to
    0.1797 kg (+1.6 %); demise at 199.0 s and 70.74 km against 199.5 s and 70.58 km; droplets 1.917 to 0.555 million
    (−71 %); median radius by number 173.3 to 173.7 µm, but by mass 227 to 1 682 µm, 7.4 times larger; largest droplet
    1.8 to 3.9 mm. The front-surface Rayleigh–Taylor mode releases 0.3 g with the rule off and 119.4 g with it on, 66.4 %
    of the sprayed mass; the thin mode 176.6 g and 60.3 g, the latter from film that reached walls below T_rigid (a
    mass-weighted 56° from the stagnation point, against 47°). The held film peaks at 10.3 g at 191 s; spraying begins
    at 176.5 s instead of 174.5 s; re-solidified 3.7 mg against 12.2 mg; the energy balance 6.0e-11.

84. **The 100 mm flight with `AA7075_scheil`, the rule on against off.** Compare within the material, not across: Scheil
    also differs from the linear range in its latent heat (390 against 400 kJ/kg) and liquid surface tension (0.80
    against 0.86 N/m; sub-plan 02, 2026-09-28). Sprayed mass 1.1724 to 1.1789 kg (+0.6 %); the body at 120 s 0.2996 to
    0.2931 kg (−2.2 %); droplets 18.10 to 13.40 million (−26 %); median radius by number 181.9 to 186.9 µm (+2.7 %) and
    by mass 192.3 to 201.6 µm (+4.8 %); thick 85.7 to 88.4 % of the sprayed mass, thin 9.3 to 2.3 % (109 to 27 g),
    Rayleigh–Taylor 5.0 to 9.3 % (59.2 to 109.4 g); before the closure 93.8 g by the thin mode against 16.0 g thin and
    64.3 g Rayleigh–Taylor; the held film peaks at 17.3 g at 49.0 s. Under the closure (medians over 142 steps) the thick
    share of the wet windward patches rises from 32.2 to 55.9 %, but only 6.2 % are thick only because of slurry
    (23.3 % on the linear law) and 14.4 % are thin by the non-rigid depth where the liquid layer said thick; the
    non-rigid depth is a median 1.45 mm, 5.6 times thinner than the linear law's 8.1 mm, as the 12.9 K band against
    79 K predicts. So the slurry test itself matters four times less on Scheil's curve, while the no-closure rule acts
    the same. Energy balance 1.1e-10, Newton 3.51 iterations per step (3.51), re-solidified 5.23 to 4.55 g.

85. **Where the held film goes: the front-surface Rayleigh–Taylor mode.** Off the closure a film on slurry is held from
    the shear modes. Some of it runs to cooler walls, where the thin mode takes it (facts 82 and 83); the rest collects on
    the cap until its depth passes the Rayleigh–Taylor depth criterion, and the applied front-surface mode sheds it within
    one growth time, as droplets of its wavelength capped by the film mass on the patch. On the flights none exceeds
    5.2 mm. In `test_film_spraying_death_and_balances` (the body at 880 K throughout, 8 s at 71 km, no closure) the mode
    sheds 5 droplets above 10 mm carrying 57.8 g, 5.2 % of the 1.118 kg sprayed (a second, unseeded run of the same
    setting: 6 droplets, 10.1–10.9 mm, 6.1 %), from facets whose film is 18–34 mm deep by mass per area — blob facets
    where the lubrication picture has already failed (fact 27); with the rule off the same setting sprays 45.5 million
    droplets, the largest 3.2 mm, against 9.0 million. The masses barely move because the release is supply-limited
    (fact 9): what the rule changes is which mechanism releases the film, and so the droplet count and sizes.

86. **The two backends, and a threshold the rule sharpens.** In the 2026-10-05 two-backend test's setting — the 6 mm pool
    at 960 K over an interior at 850 K, ten steps of 0.05 s on the coarse mesh — the whole body is slurry under the
    linear law, and the non-rigid line runs 99.95 mm, through the sphere. With the rule on the two backends make the same
    branch decision on every patch for five steps and then differ on 8 to 29 of about 4 000; every difference is a film
    of 0.01–1 µm at the spray floor `B_MIN`, which Girin's thick branch releases whole (its rate does not depend on the
    film's depth) on one side and not the other. After ten steps the runoff differs by 1.6e-4, the sprayed mass by
    4.4e-7 and the temperatures by up to 0.29 K, against about 1e-8 and 2e-4 K with the rule off. With the interior at
    820 K, below T_rigid, the rule on agrees as tightly as ever: identical branch decisions at every step, mass 6.7e-10,
    sprayed mass 1.0e-8, runoff 1.2e-8, film 6.1e-8, temperatures within 1.8e-4 K — the new two-backend test. The
    threshold is the model's, not the backends', and the rule sharpens it because it puts more micron films on thick
    patches. How much it moves a whole flight was not measured: the flights of facts 82–84 ran on one backend, and the
    seed scatter of fact 65 was measured with the rule off.

87. **Not done, and for Asha to decide.** (a) **The front-surface Rayleigh–Taylor mode over slurry off the closure —
    decide this first.** Decision (3) holds the film only from the shear modes, and the Rayleigh–Taylor mode then sheds it
    as millimetre droplets (facts 83 and 85): two-thirds of the 50 mm flight's sprayed mass and a mass median radius 7.4
    times larger. Options: (1) keep it — the mode is a different mechanism, a deceleration-driven instability of a deep
    pool, and film piling up on slurry is such a pool; (2) hold the Rayleigh–Taylor mode as well over slurry off the
    closure — nothing would then spray before 49.5 s on the 100 mm flight and almost nothing on the 50 mm flight, whose
    film (already 10–17 g held) would ride the body until the closure, the deaths hand it on, it freezes, or the body
    demises with it; (3) revisit decision (3) and keep the thin mode there. Recommendation: (1), recorded as an
    assumption, because the mode's depth criterion and growth time are its own and it is the only mechanism this model
    has for a pool too deep for the thin mode; but the droplet sizes it gives on blob facets inherit fact 27's limit, so
    the 50 mm droplet population should be quoted only with that caveat until fact 27's blob issue is addressed.
    (b) **The material for quoted branch splits.** The slurry test's own effect depends on the liquid-fraction law: 23 %
    of wet windward patches thick only because of slurry on the linear law, 6 % on Scheil's (facts 82 and 84). Sub-plan
    02 already records the decision to make `AA7075_scheil` the melting default; until that is done, branch splits
    should be quoted with their material. (c) **The time-step series with the rule on** (fact 77 (a)). The non-rigid
    depth removes the whole-element count fact 75 found the step changing; the film's one-step supply remains in the
    regime layer. Recommendation: repeat the five-run series (about five hours) after (a). (d) **Sub-plans 14 and 15 are
    not amended**: the sensitivity script has no `rigidsubstrate-off` variant, and the README and
    `docs/model_assumptions.md` entries for the rule, `T_rigid` and the four columns are not written. (e) **The 64-element
    cap** is not known to bind on the flights; where the whole body is above T_rigid the line runs through it (99.95 mm,
    fact 86), which leaves the branch unchanged — anything deeper than δ_m is thick — but truncates
    `nonrigid_depth_mean_mm`. Not measured. (f) **The φ weighting assumes the fed part of an element came evenly from
    along the line**; in a wall-owning element it is in fact the hottest part, at the wall. Not measured.

