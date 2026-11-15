# Decisions taken on Asha's behalf while building the low-Mach pressure and the solve cadence

These are the decisions Claude made without asking while the plan `2026-11-14-gas-loads-low-mach-and-cadence.md` was
built (branch `interim-gas-loads`, 2026-11-14 to 2026-11-15), in the order they were made. The plan carries out Asha's
three decisions of 2026-11-14: Daniels and Schwartzentruber's pressure law at low Mach numbers, the aerodynamics solved
only once per 0.5 s macro step rather than at every sub-step, and the Step 2 reference test on the new default loads.
Each entry says what was decided, why, and what it would cost if it is wrong. Any of them can be reversed.

"The provider" is `GasLoads`, the shared object in `reentry_model/shared/gas_loads.py` that computes the pressure, shear
and heating round the meridian; "the solve" is its Lees local-similarity boundary-layer march; "D&S" is Daniels and
Schwartzentruber's 2026 sphere model; "the joined law" is the earlier pressure (modified Newtonian on the windward side,
Prandtl-Meyer expansion past it, continued into the D&S lee); "the hold" is the provider serving one solution for a
stretch of flight time; "the skip" is the heating making no solve when the wall cannot be heated.

## Before the build

1. **"Low Mach" means below Mach 3.5, blended to the joined law by Mach 5, and the D&S pressure is anchored to the
   model's own stagnation pressure (P1).** The joined law's step at 90 degrees starts below Mach 3.39, so the D&S law is
   used alone below Mach 3.5 and blended linearly in freestream Mach number (the perfect-gas Mach number of the D&S
   condition) to the joined law at Mach 5, above which nothing changes. An abrupt switch would have changed the windward
   pressure by up to 45 % in one step between 60 and 90 degrees. The anchoring scales the D&S pressure rise so that its
   stagnation value equals the model's unreacted stagnation pressure, because the solve takes the nose velocity gradient
   from the pressure fall between the nose and the 4-degree bin, and an unanchored offset of 0.3 to 2 % at the nose
   would have distorted the stagnation heat flux. The final review found that below about Mach 0.5 the anchoring also
   carries a separate offset (entry 13). The cost if wrong is small: the two blend edges are constants.
2. **"Every 0.5 macrostep, not at every sub-step" was read as once per `--dt` of flight time (P2).** The loads were
   already computed once per macro step; the sub-steps that multiplied the solves are the continuum steps of 0.0125 s
   after the switch to continuum flow, 40 of them per 0.5 s. The provider now holds one solution for one `--dt`
   (default 0.5 s). The final reviewer listed this as Asha's reading to confirm. The cost if wrong is one argument,
   the hold length.
3. **The skip was built (P3).** Asha chose "Yes, start building" rather than "Build without the skip". The cost if
   wrong is one task to revert.
4. **A statement about D&S's validation was corrected (P4).** The plan first said D&S had validated their model against
   simulations at Mach 0.37, 2, 3, 4, 7 and 10. The research note records shear comparisons at those Mach numbers; the
   Mach numbers of their pressure comparisons were not read, because the paper could not be read as text here. The plan
   was corrected (fa6acbb), and this correction was reported to Asha. There is no cost for the code.

## Task 1: the low-Mach pressure

5. **The pressure-floor test now guards the anchored law (P5).** The test that checks the pressure is floored before
   Cantera's expansion was moved from the joined law to the anchored D&S law at its Mach 2 state, because at Mach 2
   that is the pressure that now reaches Cantera; its assertions are unchanged. The final reviewer found the floor is
   never reached below Mach 5 (the lowest blended pressure over 7,812 states is 0.0023 of the freestream pressure). The
   cost if wrong is that the floor is no longer tested on the joined law, which is used only above Mach 3.5, where no
   non-positive pressure has been seen.

## Task 2: the hold

6. **One resume test got its own flight (P6).** Under the hold, the test that a resumed flight repeats the original to
   the last bit no longer had droplets after its save point, so its precondition failed. It was given its own flight
   and save time (the flight to 0.2125 s, saved at 0.1875 s) rather than moving a shared fixture that about 15 resume
   tests rely on, or relaxing the check. The final review notes that the test passes on a single droplet record. The
   cost if wrong is one test on a different flight.
7. **A last step shorter than `--dt` is served from the window before it (P7).** That is what "at most one solve per
   `--dt` of flight" means. The measured effect on a cut-short last step is 0.05 % in heat flux, and it is documented.
   The cost if wrong is that the final row of a flight whose last step is cut short can differ slightly from before.
8. **Held fields are the solution's (P8).** Within a window the held solution's wall enthalpy and nose radius are those
   of the solve at the window's start, used as the existing code uses them; section 14.1 says so. There is no cost
   beyond what is documented.
9. **A resumed run records the hold it flew with (P9).** A resume writes `settings.gas_loads_hold_s` from the loaded
   provider, as it already writes `settings.gas_loads`, with a test. The cost if wrong is one key in the results file.
10. **The hold serves only calls strictly later than the held time (P10).** The plan's interface let a call at the same
    flight time as the solve be served from the hold. On the element model the heating calls the provider at the wall
    temperature before the thermal step and the surface flow calls it after, at the same flight time; the task review
    found the flow was being served the heating's solution, which would have changed melting element flights even when
    stepping at `--dt`. The window was amended so a call at the held time goes to the key cache instead; the element
    model again solves heating and flow separately (two solves per 0.5 s after the switch, against 80 before the hold),
    the skins solve once (against 40 before), and a flight stepping at `--dt` throughout gives the same rows as before.
    The plan carries an amendment line saying so. The cost if wrong is one extra solve per window on the element model.

## Task 3: the skip

11. **Three points from the skip's implementation (P11).** A timing test's failure (5.63 s against a 5 s bound) was
    judged to be machine load (the load average was between about 40 and 160, from other sessions' runs) and was
    re-checked in Task 4's unit tier, where it passed (0.68 s when run alone). The restart test module had to finish
    before review, and did (73 passed). The documents give the measured stagnation wall temperatures at which the skip
    starts (1878.1 K at 168 s on the 100 mm flight, 1960.6 K at 283 s on the 50 mm flight) instead of a rounded figure.
    There is no cost.

## Task 4: the reference flights

12. **Main was rebuilt to finish the classic comparison (P12).** The check that `--gas-loads classic` still reproduces
    main 4e857e1 first reported six extra "named file is missing" lines, because the earlier scratch copy of main had
    been deleted and its mesh files could not be compared. Main 4e857e1 was rebuilt from git at the same path, its four
    cases were run again, and the comparison then gave exactly the four expected lines (the new `settings.gas_loads`
    key three times, and the absolute path of the fap file). The cost was about five minutes.

## The final review

13. **The low-Mach offset is documented, not removed (P13).** The final reviewer found that the provider's stagnation
    pressure comes from Cantera's air (nitrogen and oxygen without argon), whose freestream pressure is 0.396 % above
    the atmosphere's pressure that the provider uses as the freestream value. At hypersonic speeds this is negligible;
    at subsonic speeds the dynamic pressure is small, so the offset becomes a large share of the pressure rise, and the
    anchoring scales the whole D&S distribution by it. The pressure difference from freestream comes out 1.173 times
    D&S's at Mach 0.18, 1.061 at Mach 0.3, 1.021 at Mach 0.5 and 1.002 at Mach 3. Every flight to the ground ends in
    this range (the Step 2 flights land at Mach 0.24 and 0.17), but no result measured so far is affected: the Step 2
    flights receive no heat below Mach 6.2 and have no surface flow, and the melting flights so far stop near Mach 3.5.
    The offset is older than this branch. It is now documented in sections 14.1 and 14.9 and in the provider's
    docstring. Removing it is Asha's decision, because the clean fix (one freestream pressure consistent with the
    provider's air) would also move results at and above Mach 5 by about 1 part in 100,000, which this branch's rule
    that those results stay bit-identical forbids. The cost if wrong is low-Mach pressure differences 2 to 17 % high
    below Mach 0.5 until she decides.
14. **This document was written after the fix wave (P14).** The code, tests and documents cited rulings by number while
    the ledger that held them is not part of the repository. Citations now state the rule and point here. There is no
    cost.
15. **The final reviewer's "declined to judge" items were left outside this branch (P15).** They are: the neglect of
    convective cooling of a wall hotter than the gas (older than this branch; the skip keeps it), the blend edges and
    the linear blend (entry 1), D&S's windward pressure below Mach 3.5 although the design kept the model's windward law
    because D&S's fails near the shoulder (Asha's decision 1 overrides it), the offset's effect outside the anchored law
    (the classic loads must stay identical to main), the anchoring in rarefied flow below Mach 5 (no re-entry trajectory
    is both that rarefied and that slow), the difference of about 7 parts in 10 million between the library and command
    paths of the Step 2 flights (not separated, negligible), the hold's accuracy on melting flights (quantified only on
    non-melting histories), Part B of the design, and history before this plan. Three of them come to Asha as questions:
    confirming the reading in entry 2, whether to remove the offset in entry 13, and one melting flight through Mach 5
    to 3.5 before any melting result below Mach 5 is quoted (no melting flight has yet crossed that range on the new
    pressure). There is no cost for this merge.
