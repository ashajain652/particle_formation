## Amendment of 2026-10-07 — the Rayleigh–Taylor mode stays, Scheil's curve is the melting default, a freeze-back round-off fixed, and the time-step series with the rule on (facts 88–96)

Asha's decisions of 2026-10-07 on fact 87: (a) the front-surface Rayleigh–Taylor mode is **not** held back over slurry
off Girin's closure — option (1), so the film held from the shear modes there is still shed by it (facts 83 and 85
stand, and their caveat on the sizes it gives on blob facets with them); (b) `AA7075_scheil` becomes the melting
default; (c) the switched-step series is repeated with the rule on, on Scheil's curve, the new default (chosen over
repeating it on the linear range, which would have isolated the rule's effect from the material's); (d) sub-plans 14 and
15 are brought up to date. The series found a round-off in the deep runoff's freeze-back (fact 89), which is fixed here.
It amends sub-plans 02 (no code change), 09 (the freeze-back), 13 (the default), 14 (the drivers and the series) and 15
(the README, the assumptions and the spec, for this amendment and the one of 2026-10-06). Facts 1–87 stand except where
these say otherwise; **fact 87 (a), (b), (c) and (d) are answered here, and fact 77 (a)'s question — does the droplet
population converge in the step — is answered for the rule on by facts 92–94.** Measured in a throwaway copy
(`prototype/work-2026-10-07-scheil-default/`, git-ignored) made from the copy of fact 78's amendment, which that
amendment's ten diff blocks rebuild exactly from the earlier copy (checked file by file); the measurement copy was
frozen with a checksum manifest before the runs, re-frozen once for the fix of fact 89, and never edited while a run was
in flight; every run was made on mains power under `caffeinate -i`, the power source logged at the start and the end.

88. **What changed, in date order.** (1) The default (sub-plan 13): `--melt on` without `--material` selects
    `AA7075_scheil` — Scheil's curve between 750 and 908 K, 390 kJ/kg, Σ 0.80 N/m — instead of `AA7075_range`; the four
    other packaged materials stay selectable by name, and the help text lists all five. (2) The drivers (sub-plan 14):
    `melt_verification.py`'s physics mode and `melt_sensitivity.py`'s base run use `AA7075_scheil`, and the sensitivity
    table gains `range` (the linear law) and `norigid` (`--rigid-substrate off`). (3) The freeze-back (sub-plan 09, fact
    89). The CLI test of the melting run expects the Scheil default (it failed on the material first, as it should); one
    melting test is added (fact 89). Unit tier 252 passed, 1 skipped, with only Task 11's five known reference failures;
    FEniCSx tier 12 passed.

89. **A negative film from the freeze-back's round-off, and its fix.** The series' 0.00625 s run, first launched on the
    copy with only (1) and (2), stopped at its first fine step, 49.50625 s, with the thermal solver's input check "film
    mass must be one finite non-negative value per node" (exit code 2, fact 77 (c)); the six other runs, then about
    twenty minutes in, were stopped and the seven relaunched after the fix. A probe that re-ran the seeded flight and
    checked the film and deep accounts after every stage of the melt step found the first bad value right after
    `_freeze_back`: a film of −1.29e-25 kg on one windward patch 51° from the stagnation point, whose deep liquid the
    freeze-back had just taken whole; nothing was non-finite, and the non-rigid depth, `film.lubrication`'s outputs and
    the runoff transport were clean. The cause is the 2026-10-02 amendment's arithmetic: freeze-back takes the deep
    liquid first and debits the film by the capped amount less the deep part, m_f − ((m_f + m_d) − m_d), which in
    floating point is negative by a rounding unit for about half of all pairs whose deep liquid dwarfs the film (99 860
    of 200 000 random pairs), and the spray stage's clean-up of films below 1e-30 kg runs before freeze-back. The fix
    takes each account's part directly — from_deep = min(m_d, wanted), from_film = min(m_f, wanted − from_deep) — so no
    subtraction can round below zero; amounts where the cap does not bind are those of before. A unit test with one such
    pair (2.34e-14 kg of film over 1.88e-9 kg of deep liquid) failed before the fix and passes after it. **None of the
    2026-10-06 flights is changed by it:** re-run on the fixed copy, the 100 mm flight to 120 s and the 50 mm flight on
    the linear range with the rule on and off, the 100 mm flight on Scheil's curve with the rule off, and (as the
    series' 0.5 s run, fact 91) with the rule on are bit-identical to the runs of facts 81–84 in all 90 history columns,
    every result field but the run time, and all 22 source-table columns. So facts 81–86 stand as written.

90. **The default changes no explicit run.** The re-runs of fact 89 name `--material AA7075_range` or `AA7075_scheil`
    explicitly and reproduce the runs made under the old default bit for bit; the series' 0.5 s run, which names no
    material, reproduces fact 84's explicit `AA7075_scheil` run bit for bit, so the default selects exactly that
    material. The run name does not carry the material (fact 96 (c)).

91. **The series.** The 2026-10-05 harness, unchanged (`dtswitch.py`, recovered from that session's scratchpad and kept
    with the copy; fact 69): the 100 mm flight to 120 s at 0.5 s throughout, and switched at the body Knudsen number's
    first value below 0.01 — 49.5 s, as before — to 0.05, 0.025, 0.0125 and 0.00625 s, latched; seed 1 beside the
    default 12345 at 0.05 and 0.0125 s, for the scatter between seeds at those steps (fact 75 measured it the same way);
    the same command line as fact 75's runs but without `--material`, so `AA7075_scheil`, the rule on and every other
    setting at its default. Before 49.5 s the runs with the default seed are identical: each sprays 80.3 g, 16.0 g by
    the thin mode where the wall is below T_rigid and 64.3 g by the front-surface Rayleigh–Taylor mode (fact 84), and
    the step at 49.5 s, the first under Girin's closure and still 0.5 s long, 18.1 g (the seed-1 runs 80.1 and 18.2 g).
    Values below are over the continuum window 49.5–120 s unless stated, as fact 75's are, and each is read against the
    scatter between seeds measured at its own step.
