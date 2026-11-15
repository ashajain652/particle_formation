# Decisions taken on Asha's behalf while building the interim gas loads

These are the decisions Claude made without asking while the plan `2026-11-13-interim-gas-loads.md` was built (branch
`interim-gas-loads`, 2026-11-13 to 2026-11-14), in the order they were made. Each says what was decided, why, and what
it would cost if it is wrong. Any of them can be reversed. "The provider" is `GasLoads`, the shared object in
`reentry_model/shared/gas_loads.py` that computes the pressure, shear and heating round the meridian once a step; "the
similarity solve" is the Lees local-similarity boundary-layer march in `reentry_model/shared/similarity.py`; "D&S" is
Daniels and Schwartzentruber's 2026 rarefied sphere model; "Kn_body" is the body Knudsen number, the freestream mean free
path divided by the sphere's diameter.

## Before the build

1. **The surface flow does not import the provider at module level.** `surface_flow.py` receives the provider as an
   argument and only calls its `meridian(...)` method, because the provider itself imports the surface flow's pressure
   functions and a module-level import would be circular. If this is wrong, the cost is an import-order failure, which
   the shared-package import test would catch at once.
2. **The classic comparison needs no renaming.** The plan expected the classic run's output files to need renaming
   before they could be compared with main's, but `same_results.py` runs every case under its own case name, so the new
   `_loads-classic` suffix never reaches those files. If this is wrong, the comparison reports missing files and they
   are renamed as the plan says.
3. **Code is located by name, not by the plan's line numbers.** The plan's line references date from an older main
   (08cca55); the branch starts from 4e857e1, after the geometry-stop and viscosity merges moved lines. The only cost
   is search time.
4. **The default-bin tolerances were to be decided on measurement.** If the plan's estimated tolerances for the 2-degree
   bins (3 % in the loads, 2 degrees in separation) failed by a small margin, the implementer would report the measured
   differences and Claude would choose. They passed (1.3 % in shear, 0.37 % in heating, 2 degrees in separation), so
   this decision was never exercised.

## Tasks 1 to 4: D&S model, frozen gas, similarity solve, provider

5. **Comparisons against the D&S fixture files use a relative tolerance of 1e-5.** Those files were printed to six
   significant figures, so a tighter check would fail on rounding alone. The port is also compared with the reference
   script `analysis/gas_loads/reference/ds_model.py` to 1e-12, which would expose any real error that the looser check
   could hide. The cost if wrong is none beyond that.
6. **A failed stagnation-point solve fails the flight's step instead of being swallowed.** If the similarity solve
   fails at the nose even after a fresh start, the provider re-raises the error with the freestream state in its
   message, so the flight stops at that step with its state saved and can be resumed, rather than flying on wrong loads.
   The cost if wrong is a flight that stops on a rare state; none has been met.
7. **The boundary-layer solver is capped at 20,000 mesh nodes.** One solve at the separation bin of the Knudsen 0.01
   reference state ran to 200,000 nodes and took 11 seconds, only to fail, which is read as separation anyway. The
   largest solve that converges uses 1,866 nodes, so the cap is more than ten times that. The cost if wrong is that a
   solve needing more nodes would be called separated early; the reference fixtures would catch it.
8. **About 2.2 to 2.9 seconds per meridian solve was accepted as meeting the plan's soft 2-second target.** These were
   measured on a loaded machine at the Knudsen 0.01 and 0.25 states, and the test's bound of 5 seconds holds. The cost
   if wrong is slower flights; Task 8 then measured the real cost per step (see the slowdown decision below).

## Task 5: the surface flow

9. **In thin shock-layer flow, the switch to free-molecular shear applies on the windward side only.** On the lee the
   shear is the provider's continuum shear multiplied by the slip factor. The brief's literal rule would have set the
   lee shear abruptly to zero at a patch threshold, which defeats the purpose of having shear round the whole body. The
   cost if wrong is lee shear that is somewhat high in thin shock-layer flow.
10. **The windward slip quantities stay as before; only the lee's come from the provider's pressure.** On the windward
    side the local Knudsen number, slip factor and closure choice are still computed from the classic equilibrium wall
    pressure, so that a pressure change of about 2 % cannot flip a closure at its 0.01 threshold. The cost if wrong is
    that the windward slip factor uses a pressure about 2 % different from the provider's.

## Task 6: the heating

11. **The provider gained a small public method for the wall enthalpy.** The hot-wall factor on each patch needs the
    frozen-air enthalpy at that patch's temperature, so `GasLoads.wall_enthalpy(T)` was added rather than reaching into
    private attributes or using the equilibrium-air table. The cost is a one-method addition outside Task 6's file list.
12. **Separation one bin later than the reference at Knudsen 0.01 was accepted.** The first bin where the pressure rises
    again sits at a point where the solver can go either way; the reference failed there and the provider converges.
    The loads there are about 0.1 % of the nose's. The cost if wrong is separation 2 degrees late in some states.
13. **The free-molecular heating test runs at a density of 1e-10 kg/m3.** The continuum heating falls only as the square
    root of density, so the bridge to free-molecular heating reaches 99 % weight only far into thin air. The cost is
    none.

## Task 7: wiring into the flights

14. **Sheltered corners read 180 degrees only on runs that use the provider.** Under `--gas-loads classic`, and on states
    saved before this change, they keep the old angle of just past 90 degrees, so classic stays identical to main. This
    was amended by the next decision.
15. **Only the corners the centroid test shelters are sent to 180 degrees.** The brief's literal rule would have sent
    every corner of the rear hemisphere to 180 degrees (1,129 of 2,254 on the test sphere, 974 of them in attached flow
    between 90 and 150 degrees), discarding the provider's lee loads on every default skin flight. Your decision 2a
    named "the corners the centroid test shelters", which are the 3 corners previously set just past 90 degrees, so
    only those go to 180 degrees; corners that face away keep their true angle. You asked to revisit this rule together
    with the fitted-sphere loading bug.
16. **Six smaller interface choices.** `--gas-loads` is accepted and ignored when there is no thermal model; every
    thermal run records `settings.gas_loads` and a classic run carries `_loads-classic` in its name; a similarity run
    keeps its stagnation and catalycity settings in the results and lists them under `heating_unused`;
    `--bridging-heat sesam-table` with the similarity loads exits with an error naming `--gas-loads classic`; a resume has
    no `--gas-loads` option and keeps the loads it was saved with; the element model and the non-melting body also get
    the provider by default. Each is small and reversible.
17. **Two mechanism tests now name `--gas-loads classic`.** The 71 km continuum-step latch test and the 65 km spray-supply
    failure test were built on states tuned to trigger those mechanisms under the old loads; on the new loads the 71 km
    state has no nose liquid at 0.1 s and spraying starts at 0.0875 s instead of 0.05 s. Their assertions are unchanged.
    The cost is that those two mechanisms are not exercised on the new default loads until a test is re-tuned for them.

## Task 8 and the final review

18. **The classic comparison with main was accepted as passed.** `same_results.py` reports two differences: the new
    `settings.gas_loads` key, which is expected, and the absolute path of the fap file, which names the same file from a
    different checkout. With those two keys made equal, all 21 files are identical. The documentation says so.
19. **One fix round after the final review, limited to documentation, tests and two guards.** It documents the pressure
    step at 90 degrees below about Mach 4 as a known limitation (the fix itself is your decision), corrects a comment and
    seven documentation details, adds five tests, makes a resume keep its saved loads label, raises a clear error when
    the edge enthalpy exceeds the frozen-air table's 20,000 K top, and floors the joined pressure at 1e-20 of the
    freestream pressure before Cantera's expansion. The cost if wrong is a slightly larger change.
20. **An explicit non-default `--stagnation` or `--catalycity` on the similarity loads is refused.** It exits with an
    error naming `--gas-loads classic`, because silently ignoring, for example, `--catalycity 0` would mislead. This
    amends part of decision 16. The cost is one refused combination, reversible.
21. **The final reviewer's "declined to judge" items were left outside this branch.** They are Part B's mass removal,
    the D&S lee accuracy and the zero heating behind separation, the slip closures, the SESAM blend below the switch, the
    isothermal wall, the missing convective cooling of a hot wall (older than this branch), and decision 12. The cost
    for this merge is none.
22. **Three small points from the last check were left as they are.** The provider's error message now reads "solution
    failed" and names its cause; the README credits a 16 % rise in convective heat, 6.2 % of it behind 90 degrees, to the
    element flight, which is where section 14.7 measured them; and section 14.1 explains the shear dip at Mach 3 by the
    pressure gradient parameter falling to zero, an explanation that no test checks. That section will be rewritten
    whichever fix you choose for the step at 90 degrees. The cost is one unverified explanatory sentence until then.
