# Skyscrapers architecture: open to-do items (2026-11-12)

The open items for the layered-skin (skyscrapers) architecture after the stage 0+S refactor and the skyscrapers
subpackage landed on main (655d1a2). Built from the post-refactor open-items investigation (reports in the git-ignored
`output/open_items/`) and updated with what happened since.

## Waiting on Asha's decision

1. **The fixes branch.** Seven fixes on `skyscrapers-fixes` (ce40852..e4d84f0), unit tier 856 passed and 1 skipped,
   melting results unchanged. Choose: fast-forward main and push, pull request first, or keep local.
2. **The composite viscosity law** replacing the "no data" bridge (906.4-908 K). Draft:
   `reports/near_liquidus_viscosity_composite.py` and `.png` in the main checkout; sources in
   `research_notes/Near liquidus aluminium alloy viscosity/fulltext_*.md`. Open choices: Li's per-rate fits or his
   averaged law; structural stress below 10 % solid as f_s squared or cubed; the real liquid under the high-shear limit;
   the present law kept selectable by name; the cell's history choosing its branch. Then design document, plan, build,
   and the first test (resume the d050 state at 178.5 s under both laws).
3. **Is only the material windward of the mass centroid heated?** A Mac Studio test shows the rear material heating.
   Check: which leeward rule is in force (`--leeward-test facing|centroid`) and that it does what it says; whether the
   rear heating is conduction through the body, radiation, or flux wrongly applied to rear-facing columns; whether the
   sideways heat moves heat rearward. Start from the Mac Studio run's rows and check one step by hand.
   Evidence: [figures/d050__fixes__slipband-on_cross-sections.png](figures/d050__fixes__slipband-on_cross-sections.png),
   the Mac Studio flight d050__fixes__slipband-on (cross-section through the flight direction, 171-189 s, stopped at
   189.075 s). By 179-183 s the 3D core's rear (leeward) half is at about 780-820 K, against about 700 K at 171 s; the
   skins grow and fold on the rear and sides (183-185.5 s); from 187.5 s a liquid layer appears on the rear surface,
   about x = -55 mm at 188.5 s, and the column tops enclose 133-434 cm3 against 63-65 cm3 of metal, so the geometry is
   non-physical by then (see item 4).
4. **A geometry stop for skin flights**: end the flight, with a message, once the core passes the body's centre. The
   running 50 mm Mac Studio flight has none; judge it only up to that point (earlier 50 mm flights went non-physical at
   about 185.5 s).
5. **Tell the dispersion session about the new law.** The kinked-profile build (branch `kinked-dispersion`) places its
   kink with the bridge band, which the new law removes. Its own pending decisions (`FLOOR_RE`, `KINK_DROP`, fluid-like
   columns) belong here too.

## Investigations and analysis

6. Analyse the 100 mm Mac Studio flights; move their results branch (`skin-fixes-flights-macstudio`) onto main's layout.
7. ~~Test the jump-fault rule (thick graded bottom cells against the 5 K tolerance) on the faults that prompted it
   (a re-run of about 10 h).~~ Done 2026-11-14 (fact 198). The re-run counted 30,848 faults under the old rule and
   25,620 under Task 5k's; they were a step of about 10 K across thick skins' 2-4 mm bottom cells, which a single-column
   test showed puts the rigid front 17-64 % too far. Fixed by fine cells at a skin's base (`--skin-base-cell`, 0.2 mm,
   within about 1 % on that test; `docs/model_assumptions.md` section 16); the jump-fault test and Task 5k's rule stay.
   Flown 2026-11-15 on the Mac Studio (runbook `2026-11-14-skin-base-cells-flights-macstudio.md`, results on the branch
   `skin-base-cells-flights-macstudio`): the 50 mm flight from 178.5 to 186.5 s counted 25,364 faults under Task 5k's
   rule with the base cell off and 4 with it at 0.2 mm; without it the thickest skin grew to 46 mm against 31 mm in the
   last second and took 12 g more metal into the skins. The extra wall time (14.0 against 10.9 h) was the continuum
   step's switch, item 23, not the cells.
8. Facts 163(b) and 168(a)-(e) from `we-critical-308`: the thin branch's threshold, the spray within 10 degrees of
   stagnation, the momentum-balance record, the droplet population, the subsonic patches.
9. Why the 75 deposit columns sit at about 0.7 of their steady state; whether the melt-speed result holds at other times
   and on the 100 mm sphere.
10. The two-valued slurry viscosity band (top speed changes by a median 11 % with the branch); the composite law's
    history rule would settle it.
11. Transient shear-thickening after a sudden rise in shear rate (Koke and Modigell); every law in the model is steady.
12. Whether re-frozen (dendritic) deposit needs its own viscosity; the skins do not track where a cell's metal came from.
13. The oxide skin, which may resist stripping independently of the bulk viscosity.

## Sub-plan 22's remaining work

14. Base flights: 100 mm to 120 s and 50 mm, on the Mac Studio by runbook.
15. Convergence repeats and the acceptance criterion.
16. The comparison with the element model, and the cost and memory of a flight.
17. Task 6: documentation.

## Small tidy-ups

18. A test that resumes a spraying element-model state.
19. The flight launcher returns 1, not 2, on a bad argument, and resolves relative output folders from the caller.
20. `CLAUDE.md` still quotes the port's 702-test count (now 856).
21. Left alone: the Spheral CRLF fixture test (outside this architecture) and the duplicate facts 101-107 (kept
    unnumbered so references to either set work).

## Added later on 2026-11-12

22. **Improve how the skyscrapers represent the metal's volume** (Asha, 2026-11-12). Each column keeps its base's dual
    area at every height (a cell's thickness is its mass over its density and that area) and stands along its base
    normal, so tall columns on a curved base fan apart and their tops enclose space that holds no metal. On the Mac
    Studio flight `d050__fixes__slipband-on` the volume inside the tops was 1.03 times the metal's at 179 s, 1.17 at
    183 s, 1.42 at 185.5 s and about 2 at 187.5 s (figure `figures/d050__fixes__slipband-on_cross-sections.png`); on the
    2026-11-10 comparison flights the deepest rear column, 66.5 mm on a base 19.3 mm from the centre, would be about 24 mm
    deep as a spherical shell of the same metal. The geometry stop (item 4, `--geometry-stop`, limit 1.25) only ends
    flights before this; it does not fix it. Directions to weigh: columns whose cross-section grows with height away from
    a convex base, metal moved sideways where neighbouring columns diverge, and the level-set architecture's two tracked
    interfaces (`docs/superpowers/specs/levelset/`), which have no columns.

## Added 2026-11-15

23. **The continuum step latches on a knife edge** (Asha, 2026-11-15: record it; nothing changed). A skins flight
    switches from `--dt` (0.5 s) to `--dt-continuum` (0.0125 s) at the first step that ends with any skin's top cell
    above the liquidus, 908 K (`CoupledSimulation._switch` with `SkinBody.liquid_present`), or at the Knudsen switch if
    that comes first. On the base-cell flights of 2026-11-15 (item 7) the hottest surface sat within half a kelvin of
    the liquidus from 179 to 185 s, so the switch fell at 179.5 s on the flight without the base cell and at 182.0 s on
    the one with it. The second flew 2.5 s more on 0.5 s steps, each split into 10 to 13 parts, at about 4.6 times the
    skins' wall time per second of flight, which is the whole of its 3.1 h of extra wall time; and two flights meant to
    differ in one setting differed in their step for 2.5 s. The results do not look wrong for it (the two agreed within
    about 1 % by 183 s), but comparisons between flights are muddied and the cost is set by a fraction of a kelvin.
    Directions to weigh: latch with a margin (a top cell some kelvin from the liquidus, or a liquid depth above a
    floor), latch on what the spraying first needs rather than on any liquid, and measure whether 0.5 s steps split
    into parts save anything at all once the skins are thick (on this flight they cost more per second of flight than
    the fine step).

24. **A mesh move refused on a fresh mesh, from the base cells and the fitted-sphere loads together** (Asha,
    2026-11-15: "Go with option 1; test the split by branching from this branch"; branch `refused-move-split`, off
    `skin-base-cells` at 6663f4a). Merged with main at aa0525f, the base cells made
    `test_a_spraying_body_in_the_slip_band_keeps_the_three_books` stop in its second step: a move that took one element
    of four surface corners from a radius ratio of 0.25 to 0.0486 on the freshly rebuilt mesh (the floor is 0.05), a
    5.0 mm pool's base handed back 0.36 mm beside thin skins drawn in by up to 0.44 mm. Neither change alone fails it.
    *Done:* such a move is now made in pieces before it stops the flight (`docs/model_assumptions.md` section 17); in
    halves the test passes, and so does the whole unit tier (1127 of 1128 tests, one skipped). *Still open:* no flight has split a move, so whether real flights on the base cells and the
    sphere rule refuse moves like this is not known, and `MovingMesh.split_moves` is not yet in the history, so a
    flight would not show that it split; the pieces of a move share one part of the macro step, with no conduction
    between them.
