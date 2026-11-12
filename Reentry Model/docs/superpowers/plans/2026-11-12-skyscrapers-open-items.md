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
4. **A geometry stop for skin flights**: end the flight, with a message, once the core passes the body's centre. The
   running 50 mm Mac Studio flight has none; judge it only up to that point (earlier 50 mm flights went non-physical at
   about 185.5 s).
5. **Tell the dispersion session about the new law.** The kinked-profile build (branch `kinked-dispersion`) places its
   kink with the bridge band, which the new law removes. Its own pending decisions (`FLOOR_RE`, `KINK_DROP`, fluid-like
   columns) belong here too.

## Investigations and analysis

6. Analyse the 100 mm Mac Studio flights; move their results branch (`skin-fixes-flights-macstudio`) onto main's layout.
7. Test the jump-fault rule (thick graded bottom cells against the 5 K tolerance) on the faults that prompted it
   (a re-run of about 10 h).
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
