# Decisions taken on Asha's behalf while building the fitted-sphere windward loads

These are the decisions Claude made without asking while the plan `2026-11-13-fitted-sphere-loads.md` was built (branch
`fitted-sphere-loads`, 2026-11-13 to 2026-11-15), in the order they were made. Each says what was decided, why, and what
it would cost if it is wrong. Any of them can be reversed. "The sphere rule" is `--windward-loads sphere`, the default:
one sphere fitted to the exposed windward corners, each exposed corner loaded at the angle of the ray from the fitted
centre. "Normals" is `--windward-loads normals`, the earlier rule, each corner loaded at the angle of its own surface
normal. "The provider" is `GasLoads`, the interim gas loads that landed on main as 70deee0
(`reentry_model/shared/gas_loads.py`); "classic" is `--gas-loads classic`, the loads before it. Section numbers refer to
`docs/model_assumptions.md` after the merge, where the windward loads are section 15.

## Tasks 2 to 6: the rule, its push and its documents

1. **The command line was left alone in Task 2.** The command-line tests that the new default broke were listed and
   fixed in Task 4, where `--windward-loads` was added. The cost, had it been wrong, was a branch whose command-line
   tests were failing between Tasks 2 and 4; nothing was merged in that state.
2. **Frames written under normals also carry the per-corner load angle.** The design adds the per-corner values to the
   frames, and the byte-identity check against main ignores frame files. If this is wrong, a normals run's frame files
   differ from main's by one extra array.
3. **The fitted radius is recorded as not-a-number when no valid sphere exists**, a concave fit included, as the
   design's `SphereFit` says; the separate `nose_cap_concave` flag records that the fit was concave. If this is wrong,
   the history under the sphere rule loses a concave cap's diagnostic radius.
4. **One spraying test was pinned to normals.** `test_a_profile_skin_body_sprays_and_keeps_its_books` failed under the
   sphere rule on a refused mesh move at the rear (radius ratio 0.0439 against the limit 0.05). The same move was already
   refused before Task 3 and passed only because the rebuilt mesh happened to be better, so the thin rear margin predates
   the rule; the test's subject is the spraying books, not the load rule. If this is wrong, a real spraying flight on the
   sphere default fails on a refused rear move. The resumed flights of 2026-11-14 (15.4) checked this on both saved
   bodies and found no extra refusals under the sphere rule.
5. **The unit-tier test count was left for Task 7**, after the merged run, so the count was stale on the branch until
   then (it is now 1110 in `Reentry Model/CLAUDE.md`).
6. **The final review's declined-to-judge list was set aside as the reviewer proposed**: the rule's physics, the push
   choice, the merge, the test count, the frames under normals, the plots, leeward heating past 90 degrees, the
   near-liquidus probe script, the reference test and ParaView series across a switch. Each is Asha's decision, an
   earlier decision, existing behaviour or cosmetic. If this is wrong, one of them surfaces after the merge.
7. **The push along the fitted meridian is normalised after it is projected onto the corner's tangent plane**, as the
   design specifies and as the normals rule already does. On a steep rim this keeps the full push where a shorter
   projection would weaken it; the reviewer estimated the push at most about 1.37 times too strong on the 50 mm body and
   1.06 times on the 100 mm body (15.3).
8. **The push direction is set per part of a step, and its magnitude per macro step**, as under normals. If this is
   wrong, a corner that changes shelter in the middle of a step mixes the two for that one step.
9. **The normals arm of the resumed comparison flew today's run defaults**, not the saved flight's own facing and
   initial settings, because that is the rule a new flight on main flies today. If this is wrong, the comparison in 15.4
   is against a normals variant the saved state never flew.
10. **The resumed flights' results were written into the documents by Claude directly** (15.4 and 15.7, numbers checked
    against the flights' history files), because the section said "not known until flown" and every claim must carry its
    measurement. The cost, had it been wrong, was a paragraph without a separate review; the merge review read it since.

## Task 7: the merge with the interim gas loads

11. **What each corner reads on the provider (ruling 1 of the merge).** Under the sphere rule on the provider, a corner
    the ray test shelters reads the wake angle, 180 degrees, whichever way it faces; a corner behind the mass-centroid
    plane keeps its own angle where it faces away from the flow (past 90 degrees) and reads 180 degrees where it faces
    the flow; an exposed corner keeps its ray angle, even past 90 degrees. On classic the rule is the branch's, bit for
    bit. The plan said every corner behind the plane reads 180 degrees, but it was written before the interim gas loads'
    own ruling that corners facing away keep their lee loads (R15), and following R15 makes the two windward rules agree
    on the rear of the body. If this is wrong, under the sphere rule rear corners facing away carry the provider's lee
    heating and shear where 180 degrees would give them none; it is a one-line change.
12. **The gas-loads section keeps the number 14 and the windward loads become section 15**, because the gas loads landed
    on main first. Cross-references were renumbered in the living documents and code comments, not in the historical
    specs and plans. The cost, had it been wrong, is renumbering only.
13. **Four minor findings of the merge review were fixed in the same round as the two important ones** rather than being
    deferred, because the merge was the last task and no later review would have looked at them. They added one test
    (the sphere rule's corners evaluated through the surface flow and the heating on the provider), strengthened another,
    and corrected the wording of 14.7, 15.5 and the README. The cost was one short test.

The implementer's own smaller decisions during the merge (for example that "behind the centroid's plane" is the centroid
test's own mask, and which of main's tests were pinned to normals) are listed with their costs in the merge report,
which is not kept in the repository; the tests and section 15.5 record their outcome.
