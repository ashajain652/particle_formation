# Decisions taken on Asha's behalf while building stage 3 of the level set (plan 03b, one skin depth)

These are the decisions Claude made without asking while plan
`docs/superpowers/plans/levelset/03b-stage3-one-skin-depth.md` was built on the branch `levelset-stage3` (2026-11-14 to
2026-11-16). The plan builds design section 21 of
`docs/superpowers/specs/levelset/2026-10-11-level-set-deforming-body-design.md`, the amendment of 2026-11-13: one skin
depth everywhere, and the gas surface evolved. The decisions were made by the session that coordinated the build (the
controller), on the reports of the agents that built and reviewed each task, and were recorded in the plan's working
ledger, which was kept outside the repository and is deleted once the branch's final review is clean; this document
replaces it as the record. Each decision is numbered, keeps its id from the ledger (P1 to P5 for the pre-flight scan,
T1-a onward for the tasks), and says what was decided, why, and what it would cost if it is wrong. Any of them can be
reversed. The rulings made while the plan was written (planning rulings 1 to 11) are in the plan itself and are not
repeated here.

Terms used throughout. φ (phi) is the gas surface, the level set that is evolved; β (beta) is the skin base, the zero
contour of φ + H, derived from φ; H is the skin depth, `SkinSettings.depth`, 1.2 mm by default. A marker is a point on φ
carrying one skin; its column is the part of the layer between φ and β that it owns; the division (`division.divide`)
gives every column its cross-section by depth, its shared boundaries and its group of β's patches. A step is a 0.5 s
macro step; a part is one of the shorter pieces a step is split into, each ending in an acceptance (φ moved, the
division rebuilt, every column filled). The cases (`analysis/levelset/stage3_cases.py`) are the plates (one column of 8
mm of metal: cold, 700 K at 6 MW/m², and hot, 880 K at 2 MW/m²), the q6p box (8 by 10 by 10 mm, 700 K, 6 MW/m² on its 10
by 10 face), the shoulder (a segment of the 50 mm sphere about its 45 degree station) and the cost cases; the test cube
is the 10 mm cube of the body's tests, and the spraying cube is the `sprayed` fixture of
`tests/levelset/test_levelset_skinbody.py`, a 10 mm cube with 2 mm rounded edges. The record is
`analysis/levelset/stage3_recorded.json`. A strict expected failure is a test kept with its tolerance and marked
`xfail(strict=True)`. D1, D2 and D3 are the three diagnoses made in scratch after Task 7 (D1 the box's uneven face, D2
the spraying cube's hang and the shoulder's move-sum stop, both before Task 7b; D3 the box's mass books, before Task
7d). Section numbers of `docs/model_assumptions.md` are those after the merge of main at 9993a77, where stage 3 is
section 18. A measurement labelled S1 to S20 is described, with what was run, on which commit and how often, under
"Where the scratch measurements came from" below; the scratch folders themselves will not survive. A review finding or a
task report that the code, the tests or the record cites by name is listed, numbered F1 onward, under "Review findings
and reports that the code cites".

## Before Task 1: the pre-flight scan

1. **P1.** Task 2's sharpness tests take as their reference the ideal shape's exact signed distance evaluated on the
   same mesh nodes (`inside_volumes(mesh, sdf_ideal(nodes))`: the sharp box shrunk by H, the box after an equal
   recession, the rounded offset `sdf_box − d`, the eroded L `sdf_L + d`), with the tolerance kept at a tenth of the
   rounding volume and the analytic volume reported beside it. Why: the box is offset 0.13 mm from the grid, so φ's
   piecewise-linear zero contour cuts every edge with a grid error as large as the tolerance, and the deposit test's
   expected value had left out the corners' term, (8 − 4π/3)d³ = 0.48 mm³ at d = 0.5 mm, 79 per cent of its tolerance.
   With the grid's own reference the grid's cut at the corners lies on both sides, so the test judges sharp against
   rounded, which is planning ruling 11's intent; design 21.4 says sharpness below the grid is sub-grid by design. If
   wrong: a move that rounds edges at sub-grid scale passes; the analytic numbers reported beside it would show it.
2. **P2.** Task 4's reference for a temperature step relaxing along a row of columns is the error-function solution
   summed over the images of the row's insulated ends, which is exact on the finite row, with the tolerance kept at 1
   per cent of the step. Why: the row is 40 columns of 0.25 mm, 10 mm with insulated ends, and at a diffusion length of
   2 mm the plain error function misses the insulated end by erfc(1.25)/2 times 100 K, about 3.9 K, beyond the
   tolerance. If wrong: nothing; the reference is exact.
3. **P3.** When Task 6 deletes `Sheet.of`, it keeps Task 4's
   `test_depth_faces_with_constant_lengths_transport_as_before` by building the reference sheet from its constructor
   with the old face values on the flat sheet (or recorded transport values), rather than deleting the test. Why: the
   test pins that depth faces with constant lengths reduce to the old transport, and Task 4 wrote it on `Sheet.of`,
   which Task 6 deletes. If wrong: one test.
4. **P4.** The re-partition by depth (`combine(by="depth")` in `repartition.py`) is Task 6's, where its test is; Task 5
   leaves `repartition.py` re-partitioning by mass and then fills. Why: `repartition.py` is only in Task 6's file list,
   and design 21.8 item 4b's test is Task 6's. If wrong: between Tasks 5 and 6 a re-partition did not keep each depth's
   metal at its depth; Task 6 restored it.
5. **P5.** Task 10 records the spraying table path's records and released masses at its starting commit, before
   `spraying.py` is edited, as a committed fixture under `tests/levelset/fixtures/`, and its test compares to it byte
   for byte. Why: "unchanged" needs a recorded "before". If wrong: one fixture file.

## Task 1: the division of the layer

6. **T1-a.** The division's measure is the coarea weight: each level triangle's area divided by the length of φ's
   gradient in its cell. Why: design 21.14 leaves the integration to the builder, and the raw level areas fall short by
   0.7 to 5.2 per cent on edged bodies, because the piecewise-linear |∇φ| is not one in the edge cells; the weight
   changes nothing on flat faces and spheres. If wrong: one weight in `divide`.
7. **T1-b.** `DIVISION_TOL` stays 1e-3, and the body records the closure factor and never raises on it (planning ruling
   1). The number of levels, `DIVISION_LEVELS`, is the smallest of 24, 48 and 96 for which the box's factor is within
   1e-3 of one at both 0.7 mm and 0.175 mm cells, provided the division then takes at most four times the 24-level time
   at 0.7 mm; otherwise 24. Why: design 21.8 item 5 asks for a tolerance from the division's integration, recorded, and
   the noise at 24 levels on the body's meshes was 2e-3 to 1.2e-2 before Task 1's fix round. If wrong: a per-piece
   scaling of the areas by up to about 1 per cent on body meshes. As measured first (S1), 96 levels qualified (at 0.7 mm
   the factor's departure was −2.04e-3, +3.32e-3 and +5.60e-4 at 24, 48 and 96 levels, the division taking 2.9, 5.6 and
   11.4 s; at 0.175 mm −6.55e-4, −8.14e-5 and −1.55e-5, taking 32.0, 53.9 and 91.7 s; 96 levels took 3.9 times the
   24-level time); after the fix round (number 14) the same rule gave 24 levels again, because 48 levels then departed
   by 3.3e-3 at 0.7 mm.
8. **T1-c.** The feet of the level triangles on φ are found along the recovered (averaged) gradient, through an added
   keyword `ownership(..., measure=None)`. Why: gradients taken cell by cell open gaps at convex edges and lose the
   boundary across the edge that the design needs. If wrong: one keyword.
9. **T1-d.** Metal thinner than 2H is detected through the mesh's nodes: boundaries of zero length between columns whose
   normals oppose are flagged. It reports metal thinner than 2H as the mesh resolves it, to within about one cell where
   the cells are coarser than β's interior, which is the geometry the body uses. Why: where two sides of a thin sheet
   meet on its mid-plane the contact is parallel to the levels and has no length on any level. If wrong: a stop on metal
   between 2H and 2H plus one cell thick; recorded in the documents (section 18.7).
10. **T1-e.** The division's area is the column's volume per unit depth (the coarea), which the skins fill; the gas
    loads and the heat keep the markers' owned areas of φ (Task 6). Task 1 measured, on the box, the largest relative
    difference between the division's top area and the owned area. Why: the two measure different things, and the fluxes
    act on the surface. If wrong: a few per cent between the area the heat comes in on and the cross-section filled, at
    edge columns. Measured (S1): up to 10.4 per cent at edge columns of the box at 0.7 mm cells (95th percentile 8.6 per
    cent, median 0.20 per cent).
11. **T1-f.** The box test's column selection classifies a marker's face by its normal's largest component. Why: the
    test's selection of columns had failed, not the division. The ledger records no cost for this ruling; it changes
    which columns one test inspects.
12. **T1-g.** `Skins.set_profile(levels, area)` sets the cross-section profile only (the cells' volumes, the fill and
    the interfaces); `Skins.area`, the area the top fluxes act on (gas heat, radiation, spraying), stays the markers'
    owned areas that the body sets, not the division's top area. Why: the coarea top area is the piecewise-linear volume
    per unit depth and exceeds the surface's area by up to 10.4 per cent at edge columns, and `Skins.area` multiplies
    the top flux and feeds the spraying. If wrong: at edges a top cell's flux area and its profile area differ (the
    books still close).
13. **T1-h.** "Area strictly decreasing with depth" is judged band by band, each band one mesh cell deep, and the base
    below the top. Why: level by level the coarea of an edge column jumps by 12 to 44 per cent where level surfaces
    cross planes of nodes, an artefact of the piecewise-linear mesh. If wrong: nothing in Task 1. It was parked for
    Tasks 3 and 4, to smooth the profile over one-cell bands if the jagged profile disturbed conduction or the fill;
    ruling T3-c closed it without smoothing.
14. **T1-i.** The same-piece rule: a level triangle (and its foot) in one piece of the body may be owned only by a
    marker on that piece, so every column lies in one piece and the closure factor is a true measure per piece. Why:
    design 21.5 says ownership never crosses gas, and two pieces are by definition separated by gas; the review found a
    column spanning two pieces counted wholly in the larger, which shifted the factors by about 2e-3 at a 0.7 mm gap. If
    wrong: one restriction in `divide`. It was carried as the same rule into Task 2's groups and Task 6's owned areas.
    (In the same fix round, cells too thin in depth for the levels were given their exact volume; that change followed
    the review's other finding and was not a ruling.)

## Task 2: the skin base, and φ moved by exact volumes

15. **T2-a.** `move_by_volumes` returns φ reset to a signed distance (redistanced within the band at each round) and
    measures the miss on that φ; the correction of each piece's volume is the body's (Task 5); and the groups' swept
    volumes must sum to φ's total change. Why: planning ruling 2 left it unclear who resets φ. If wrong: where one call
    sits.
16. **T2-b.** The move's final reset is `redistance(keep=True)`: the nodes of the cells the zero contour cuts keep their
    values, so "signed distances after every reset" means except at those nodes. Task 5 uses it wherever it resets φ,
    the remesh included; `skin_base` was to keep the plain reset (overruled for β by T5-c). Why: planning ruling 2
    assumed the reset keeps the zero contour, and the plain piecewise-linear reset does not: it eroded the box's convex
    edges by 0.26 mm³ a call and raised a converged move's miss from 8.7e-7 to 7.4e-3 (S2); keeping the cut cells is the
    standard interface-preserving reinitialisation. If wrong: the kept nodes are not exact distances (within about 0.4
    of a cell near convex edges, estimated), so the feet and recovered gradients there carry that error.
17. **T2-c.** Each cell cut by a move is split exactly among the markers owning its points' feet, not given whole to one
    group. Why: whole cells disagreed with the division's columns by −21 to +14 per cent, which would make faces recede
    unevenly; the exact split agrees to −0.73 and +0.52 per cent (S2). If wrong: compute time, 4.1 s per box move.
18. **T2-d.** The offsets are found by Newton's method with the exact Jacobian and step halving; `MOVE_ROUNDS` bounds
    the Newton iterations within a round while `Moved.rounds` counts the rounds of a long move. Why: the groups are
    coupled through shared cells, and a secant per group left a miss of 1e-2 after 8 rounds (S2). If wrong: nothing to
    results; naming only.
19. **T2-e.** The test of extra loss along an edge uses 0.2 mm, since the brief gave no depth. Parked for Task 5: a
    group asked for more than its region holds returns a miss instead of meeting it (the box's corners exhaust within
    about 1 mm), so the body must decide where unmet volume goes and keep its band larger than each round's move plus
    the cut cells. The ledger records no cost for this ruling; Task 5 took the unmet volume into the per-piece
    correction and its `move_miss_max`.
20. **T2-f.** Newton's method starts from each group's own volume profile at eight offsets (`_first_offsets`) rather
    than from the change in volume over the top area. Why: from the latter it reached only 4.9e-6 after 8 iterations,
    because edge columns narrow with depth. If wrong: about eight extra coverings per move.
21. **T2-g.** `MOVE_FLOOR`, 1e-3 of the move's largest change, is the floor of the relative miss. Why: a group asked for
    nearly nothing has no meaningful relative miss. If wrong: a nearly empty group's miss is judged on the move's scale.
22. **T2-h.** The six box move tests stay in the unit tier, unmarked. Why: the plan's rule about box runs means body
    runs, and these are 2 to 6 s geometry checks of one function. If wrong: about 30 s more on the unit tier.

## Task 3: skins on one depth, with measured cross-sections

23. **T3-a.** With a profile set, the skins' divide, cut, combine and remove scale, sum or drop the profile's rows as
    they do the area. Why: the brief is silent on profiled columns that split and merge, and parts split and merge
    between acceptances. If wrong: a split column's volume would not match its share until the next `set_profile`.
24. **T3-b.** The old stub checks in `SkinSettings` stayed beside the new depth checks until Task 5 deleted the stub.
    Why: the brief replaced the two-cell check's argument, and the stub still drove the hand-down and the minimum mass.
    If wrong: one redundant check, for two tasks.
25. **T3-c.** The measured profile is used unsmoothed, closing T1-h. Why: the trigger (a factor of two between
    neighbouring conductances) was a proxy for the failures it guards against (open books, thin cells, sign changes,
    oscillating temperatures), none of which occurred; on the box's real division 28 of 5,406 neighbouring conductance
    pairs differ by more than two, at corner columns 1.01 to 1.10 mm deep (2.84 at most where not pinched, 22.3 in a
    pinched column), the regrid's grading alone gives 1.5 on flat columns, and the prescribed smoothing left the largest
    at 2.85 (S3). If wrong: corner columns near H conduct through interface areas carrying level-by-level noise (195
    reversals over 132 columns), not seen in their temperatures.
26. **T3-d.** Below an overfull column's end the cross-section at H is continued (the mean section for a pinched
    column), as an interim rule. Why: columns are exactly full at every acceptance, so overfill exists only within a
    step until the base exchange returns it. If wrong: within a step, swollen corner columns place their excess deeper
    than it is (2.26 mm in Task 3's stress run). Parked for Task 5 and closed by T5-f: the deepest cell below a column's
    end within a part was 0.152 mm, under a mesh cell.

## Task 4: sideways heat and runoff on depth-resolved faces

27. **T4-a.** Task 4 made `runoff.landing_conduct`'s skins carry the profile; the thickness step, the re-partition's
    stretch, the step end and the body's saved state were left to Tasks 5 and 6. Why: `landing_conduct` conducts the
    metal Task 4's transport lands, and the others are in Tasks 5 and 6's files. If wrong: landed metal conducted on the
    old curvature path until Task 5.
28. **T4-b.** The sideways solve is the builder's choice (design 21.14) but must be implicit, in the same solve as the
    conduction through the depth, with exact books; one conduction call on the box was measured with and without the
    sideways faces, for the cost of the level set's to-do item 5. If wrong: a slower solve than needed. Built as one
    Newton iteration over all cells, its Jacobian solved by GMRES preconditioned by each column's tridiagonal solve
    (`SIDE_RTOL` 1e-12, 18 to 45 iterations; a sparse direct solve was 75 times slower). Measured once: a box heat
    sub-step took 0.011 s without and 0.144 s with sideways heat in development (13 times), 0.014 and 0.228 s in
    production (16 times), about 9 microseconds per cell. (Task 12 re-measured it like for like, S20.)
29. **T4-c.** The trial pass's derivative of a column's base heat with sideways heat on is the one for a common shift of
    every base temperature. Why: an exact derivative per column needs one solve per column, and the mismatch at the
    interface is booked. If wrong: the linearisation of the coupling between the volume and the skins is approximate
    across neighbours. The mismatch with sideways heat on, measured in Task 6, was 7.26e-4 of the heat in on the box.
30. **T4-d.** The thickness step's predicted fluxes use the depth-resolved boundary lengths (the areas untouched). Why:
    the implicit thickness term must predict the transport that follows it. If wrong: nothing to the books.
31. **T4-e.** Metal that lands flows on with the profile's top area, as planning ruling 5's A_i(z_k) says. If wrong: up
    to 10 per cent in the landing rates of edge columns.

## Task 5: the body moves φ and keeps every column full

32. **T5-a.** The thickness step's height per unit of top volume uses the profile's top area, consistent with the part
    ratio's; the step end's depth from volume over area becomes the inversion of the profile; and the body's saved part
    state carries the profile. Why: the division's top cross-section is the column's volume per unit depth at the
    surface. If wrong: the thickness term's heights are off by up to 10 per cent at edge columns.
33. **T5-b (superseded by T5-l and T5-m).** No change to Task 1's division for the bumps that grew on edges under
    repeated recession, yet. Why: their growth per 0.1 mm step fell, 0.81, 2.14, 1.68, 0.84, 0.20 and −0.11 mm³,
    saturating at an equivalent edge radius of about 0.47 mm, under the 0.7 mm edge cell. A measurement was set before
    the review: twelve steps on the same cube, accepted if the last four steps' growth summed to at most 0.3 mm³ in
    magnitude, after which Task 7 would repeat it on the box's own mesh; otherwise the Task 1 change (building the
    division on the move's exact cell split) would come into Task 5's fix round. If wrong: Task 7's cases carry
    numerical edge rounding of up to an edge cell, which reads as melt rounding. The criterion failed (S4: the last four
    steps grew +1.159, +0.761, +0.156 and −0.186 mm³, 2.26 mm³ in magnitude), which led to T5-l.
34. **T5-c.** `skin_base` resets φ + H with `keep=True`, overruling T2-b's plain reset for β. Why: the plain reset left
    β short of φ + H by 1.6e-3 of the body's volume, so the volume's metal and the layer would not tile and the 1e-12
    mass books would fail. If wrong: β's node values in cut cells are not exact distances.
35. **T5-d.** Drawn metal capped at the rigid temperature was accepted as interim, with a measurement before the review:
    rerun the case that needed the cap without it, and remove the cap if its part count stays within twice the capped
    run's (restoring design 21.3), or else keep it as an open departure reported to Asha. Why: stage 3's full-depth draw
    had the same cap, and the enthalpy it leaves is booked as heat owed to the volume, so the books close. If wrong:
    drawn metal enters the skin colder than the volume, which keeps superheat at the base. On the fixed code the
    uncapped step took 10 parts against 9 (S5), so the cap was removed (718cf7e).
36. **T5-e.** The markers' areas are given afresh on the moved φ at every acceptance (`_own`). Why: with fixed areas the
    edges receded 1.9 mm within a step and the part ratio (a part's largest move over its allowance) reached 2.8e12; the
    part ratio's top area must be the current one. If wrong: one pass over the areas per acceptance. Carried into Task
    6: the owned areas of T1-e are recomputed at every acceptance too.
37. **T5-f.** Two fixes to the skins' regrid were accepted, and T3-d closed: the deepest cell below a column's end
    within a part was 0.152 mm, under a mesh cell. Why: the fixes cure a measured cycle of merges and splits, and an
    excess placed 26 m below its column's end. If wrong: the excess of narrowing columns sits shallower than Task 3's
    rule put it.
38. **T5-g.** `surfaces.outer_surface` is kept. Why: `analysis/levelset/measure_geometry.py` imports it and was outside
    Task 5's files. If wrong (cost): one function kept; the final review may delete it with that script's one-line
    change.
39. **T5-h.** The spraying fixture's stop on the thin-metal check (two tests marked as expected failures at run time)
    was accepted as interim and carried into Task 6: with the runoff on the division's faces, the fixture's 0.05 s step
    must run without the stop and the marks go; if it still stopped, Task 6 was to report where the folds form and
    whether halving the part length halves the rise, and stop for a ruling. Why: the body still transported on the old
    gas-surface faces while the whole 1.2 mm skin now flowed, and Task 6 is where runoff crosses edges on depth faces.
    If wrong: if the folds are real piling with no spray, the fixture must change, not the code. (Task 6 found capillary
    folding at the sharp edges, rulings T6-a and T6-d; the diagnostic run between Tasks 5 and 6, S7, showed that the
    swelling of drawn metal does not cause it.)
40. **T5-i.** The carry that the step end's own settle and fill leave (the volume change of molten metal returned to the
    volume, which holds its metal at the solid's density) waits for the next step's first move of φ; planning ruling 4's
    "before the closure check" is read as "the closure check bounds the carry left". The carry was to be measured on the
    heated cube and the test to assert twice the measured values; if the carry left exceeded 1e-3 of the body's volume,
    one more carry move, division and fill would follow the step end's fill instead. Why: that extra round costs a
    division (about 8.7 s, 6 per cent of a melting step), and the carry is a density difference of drawn metal, so φ
    lags the metal's true volume by it for one step. If wrong: φ lags by the step end's carry for one step. Measured
    (S6): at most 8.1e-12 m³ left (8e-6 of the body) and at most 1.42e-3 of a column; the test asserts 1.6e-11 m³ and
    2.8e-3.
41. **T5-j.** `division.py` (Task 1's module) could be edited in Task 5's fix round: `divide` returns the ownership of
    its base level (the marker owning each level-H triangle, with the plain-area measure too), so ownership is computed
    once per division and the body's grouping of β's patches is built from it; a test asserts each β patch triangle's
    owner is the division's level-H owner exactly, and the group test regains its checks or states why they cannot hold.
    Why: the review found that the body's `_beta_owners` copied the division's ownership rule rather than using it, so
    two copies of one rule could drift. The ledger records no cost for this ruling. Built: 3,900 patches identical to
    the division's owners, and the group test's checks restored (0.952 to 1.001).
42. **T5-k.** Each column with a base draws at its own base temperature, the area mean of its patch group's temperatures
    on β, computed once per fill; columns with no base keep drawing from their neighbours. T5-d's comparison was rerun
    after this fix. Why: the review found the draw's temperature read at x − H n, outside β for all 58 edge and corner
    columns of the test cube (median 0.39 mm outside, at most 0.55 mm), where the solve holds ghost or stale values,
    which may have driven the refill runaway the cap was hiding. The ledger records no cost for this ruling. Measured
    once: a corner column 0.549 mm outside β drew at 855.50 K where its base point read 865.56 K.
43. **T5-l.** Before the Task 1 change T5-b named, diagnose. Why: the 0.7 per cent mismatch per column between the
    division and the move can explain at most about 0.08 mm³ a step (0.7 micrometres over a 1 mm band along 120 mm of
    edge), against 0.8 to 2.2 mm³ measured, so that remedy might not cure it. The diagnostic: (A) the twelve steps with
    each group asked its exact share of an ideal uniform 0.1 mm offset; (B) each kind of column's ratio of the
    division's ask to that share; if (A) stayed within 1.0 mm³ and (B)'s edge or corner ratios were off by more than 5
    per cent, the remedy would be the division's volumes made to agree with the exact split (designed and reported
    first); if (A) also grew, bisect the body's step. If wrong: one diagnostic run, under 45 minutes. Result (S4): see
    T5-m.
44. **T5-m.** Recession is judged against the exact recession (φ0 + kD) on a fully solid cube (or with swelling off),
    not against the shrunk cube sampled on the nodes; no change to Task 1's division now. Why: the diagnostic (S4)
    showed that the node-sampled measure swings by up to 2.2 mm³ on its own, that at 850 K the growth came from the
    density carry (drawn metal swelling about 2.3 per cent), and that at 700 K, fully solid, exact asks track the exact
    recession to 0.017 mm³ while the division's asks depart by 2.6 mm³ over 1.2 mm with no net change of volume, an
    equivalent edge rounding of about 0.32 mm on 0.7 mm edge cells whose growth slowed (2.0 mm³ by step 5, 0.6 mm³ over
    steps 6 to 12). Carried into Task 7: repeat the measure on the box's own mesh and make the Task 1 change if after
    1.2 mm the departure exceeds one edge cell's rounding volume, (1 − π/4) times the cell squared times the edge
    length, or the last four steps grow it by more than 10 per cent. If wrong: a shape error of up to half an edge cell
    through Task 7's cases. (Task 7's measurement: T7-c.)
45. **T5-n.** The body follows design 21.6 as written: the volume holds its metal at the solid's density, so mushy or
    molten metal swells only when drawn into a skin. Raised with Asha as a design question and not changed without her
    word. Why: the specification binds. If wrong: φ sits low by the volume's unrealised swelling (up to about 15 per
    cent of the molten depth below β; liquid 2,400 against 2,813 kg/m³) and rises unevenly as metal is drawn.

## Task 6: the body on the division's faces

46. **T6-a.** No change to the capillary push in Task 6; the spraying fixture's expected failures stay with the
    diagnosed reason, and sharp-edge capillarity is raised with Asha (design section 17 item 8, the curvature from the
    level set, is hers to decide, and design 21.11 defers the thin case). Two scratch measurements were made for her
    (S8): eight parts on the same saved cube, and the same step on a cube with edges rounded to 2 mm, surface tension
    on. Why: the specification binds the curvature's source and defers the thin case, and convergence with the part
    count separates a numerical defect from physics. Task 6 had found the fold to be the capillary push from φ's
    recovered curvature at the cube's sharp edges (median 2.5e5 N/m³ against a pressure gradient of 3.8e4 N/m³), raising
    a ridge 1.0 mm high at one or two parts and 1.6 mm at four, with the transport converged; with surface tension off
    the face receded evenly (0.65 mm) and stopped on one pair at the receding corner's rim. If wrong: spraying through
    the body stayed unverified into Task 7. Result: at eight parts the cube still stops, the ridge +1.53 mm and the
    trough −1.04 mm against +1.61 and −0.99 mm at four, so the fold is converged in the part count; the rounded cube
    completes its step without the stop.
47. **T6-b.** Planning ruling 5 amended: the sideways conductance between two cells divides by the distance between the
    two columns' centroids over their shared span of depth (per level, from the division), keeping the harmonic mean of
    the conductivities and the shared length over the overlap. The box-edge test compares the body's column cell at the
    probe's depth with the reference's mean over that column's part of the layer at that depth, `EDGE_RTOL` 5 per cent
    unchanged; its strict expected failure is removed if it passes and kept with both numbers if not. Why: across the
    edge the centroids are 0.56 of the markers' distance apart at 0.5 mm depth (0.47 at 0.9 mm), so the conductance was
    about 1.8 times too weak there, and 0.60 to 1.84 times along the face; and the body's resolution is the column,
    across which the reference varies by 8 per cent (30.71 against 28.37 K). If wrong: every non-uniform sideways
    conductance changes, and the test is blind to variation within a column of up to 8 per cent. Result: the probe's
    column rose 26.30 K against the reference's 28.37 K column mean (7.3 per cent under; 9.3 per cent before), so the
    strict expected failure stays; the box's interface mismatch went from 7.85e-4 to 7.26e-4.
48. **T6-c.** `Markers.merge` combines its pair of skins by depth (`combine(by="depth")`), as design 21.5 requires when
    markers merge and planning ruling 6 does for a release; a test checks mass, enthalpy and momentum exactly and each
    level's metal at its level to 1e-12. Why: the design binds; merging by mass had moved depths by 0.63 H. If wrong:
    one call.
49. **T6-d (provisional; Asha told).** The spraying fixture's cube gets 2 mm rounded edges (more than H), since the
    fixture tests spraying and the books, not sharp-edge capillarity, and the rounded cube's step completes (S8). Made
    in Task 7. Sharp-edge capillarity (cap or smooth the curvature, or leave it to designs 21.10 and 21.11) stays Asha's
    question. If wrong: one fixture builder, and the spraying tests judged on a rounded cube.
50. **T6-e (amended by number 63).** The box-edge heat, 7.3 per cent under against 5, is parked with its strict expected
    failure and both numbers; its cause was not yet known, and it was carried to Task 7, whose box records its edges,
    and to the final review. Why: changing the conductances further without a cause would be a guess. If wrong: edge
    temperatures 7 per cent low near the box's edges until diagnosed.

## Task 7 and diagnoses D1 and D2: the verification cases on the new geometry

51. **T7-a.** Task 7 is reviewed as delivered (its brief allowed a failing case with a strict expected failure and the
    measured reason), and the box's uneven recession is diagnosed before Task 8 by a scratch-only diagnosis, D1, whose
    fix becomes Task 7b with its own review. The box had stopped on the thin-metal check at the end of its second step
    (1.0 s, 22 pairs), its face centre 0.907 mm against the reference's 0.789 mm (15 per cent ahead), with a trough 1.2
    to 1.4 mm deep 1 to 2 mm inside the heated face's edges and lips of up to 0.41 mm at them. Why: the box is design
    21.8's central check and the reason for the new geometry, and a numerical error in the face carried into Tasks 9 to
    11 would compound; Task 8's gate was likely shut anyway. If wrong: one diagnosis's time, if the cause needed Asha's
    design decision anyway.
52. **T7-b.** A hang is never acceptable. The spraying cube's part-count test hung (`split_part` asked for about 4.3e9
    sub-parts, for a rounded-edge column with a top of 2.4e-16 m² holding 1.35 mg), and the shoulder stopped in its
    third step on `MOVE_SUM_RTOL` (its groups swept 6.5e-12 m³ less than φ's change, 1.6e-3 of the move). Both are
    diagnosed together in scratch, D2, since both show a column's geometry disagreeing with its contents within a step
    at a curved edge, and the fix joins Task 7b; whatever the cause, `split_part` must stop with a RuntimeError naming
    the column past a named cap rather than hang. Why: a flight on the Mac Studio could hang for hours. If wrong:
    nothing beyond the diagnosis.
53. **T7-c.** T5-m's Task 1 change is not triggered by the box's own mesh (S10). The departure from the exact recession
    after 1.2 mm, 0.784 mm³, exceeds the criterion at the measured 0.175 mm edge cells by 7 per cent (the starting box's
    edges) to 43 per cent (the receded box's), but it equals an edge radius of 0.21 mm, about one edge cell, which this
    mesh already reads for a sharp edge (0.07 to 0.19 mm at the start), and it stopped growing (+0.4 per cent over the
    last four steps, against 10). T5-m guards against rounding that accumulates beyond a cell, and this does not. If
    wrong: the box's edges are rounded by about one cell (0.2 mm) through a run. Passed to the final review.
54. **T7-d.** The hot plate's base criterion stays out of its judgment, and stays recorded. Why: the brief judges the
    hot plate on its recession and non-rigid depth only; the criterion served design 20.6's full-depth skin, which
    design 21.13 retires, and with one depth the base lies in non-rigid metal whenever that is deeper than H. If wrong:
    an assertion Task 11 could restore.
55. **T7-e.** The cases' books check the closure factor against `DIVISION_TOL` (an added `abs(layer_factor_max) <=
    DIVISION_TOL`), in Task 7's fix round. Why: planning ruling 1 says a factor more than `DIVISION_TOL` from one fails
    the step's check, "raised by tests", and the cases' books did not check it; the plates (1e-14) and the box (5.9e-4)
    still pass, and the shoulder was already a strict expected failure. If wrong: one line.
56. **T7-f.** The speck fix (D2's option 1, S13): the nodes of a speck that a move lifts out of the metal are left out
    of the volume correction's claims, and after each correction a check raises if a piece made wholly of cells that
    held no metal reappears; if that did not clear both replays, the speck would be lifted again after the correction
    (option 2, tried in scratch). Why: D2 found the root cause numerical and inside design 21.6: the speck rule lifts a
    parted speck's nodes to 1 nm outside, the correction then shifts nodes within two rings by more than that (0.24
    micrometres on the cube, 3.1 on the shoulder) and the speck returns. If wrong: one more call to `levelset.pieces`
    per move.
57. **T7-g.** Clear stops: the fill raises, naming the marker, position, piece and volumes, when a column's miss stays
    above `FILL_RTOL` after its rounds; the division raises when a piece of φ has layer but no column, or its factor
    departs from one by more than `DIVISION_GROSS` = 0.1; and a split is bounded by `PARTS_MAX` = 1000 sub-parts,
    raising with the worst column named. Why: a hang or a silent skip is never acceptable, and every completed recorded
    step had filled to 4e-16 and closed within 4e-3, so none of these fires on a good step. If wrong: a stop where a
    looser rule would have run on.
58. **T7-h.** `DIVISION_TOL` stays 1e-3, and the shoulder is built 0.13 mm off its mesh's planes of nodes, as the box
    and the test cube are. Why: its closure departure of 3.67e-3 is an artefact of flat faces lying on node planes (451
    nodes at φ = 0), about a hundred times smaller off the grid (S13). If wrong: the shoulder's geometry moves 0.13 mm,
    which no judged quantity depends on.
59. **T7-i (amended by T7b-a).** Heat booking (D1's remedy 3a): each marker's top flux is taken under
    `Markers.owned_patches` and taken again on the moved surface at every part, from the same loads state (no new gas
    evaluation; flux and runoff forcing), in flights and in the case script. Why: design 21.5 (one partition of φ for
    the areas and the fluxes) and 21.6 (areas afresh on the moved surface); on the box it brought the heat in from +5.5
    per cent a step over the beam's power to +0.5 per cent (S12). If wrong (cost): a pass over the patches and a patch
    mean per part, a few seconds.
60. **T7-j.** The box's failure is 2 mm markers not resolving a sharp heated edge (D1, S12), which design 21.10 plans
    after 21.9; the plan's order stands and the box keeps its strict expected failure with this cause, and Asha is asked
    whether to bring 21.10 forward. Why: the specification orders 21.10 after 21.9, and no later task builds on the box
    passing. If wrong: sharp heated edges stay unverified until 21.10.
61. **T7-k.** A Task 7c after 7b diagnoses and fixes the thin-metal test so that it flags only metal thinner than 2H
    (design 21.11). Why: D1 and D2 found pairs flagged between markers 6.8 to 10.1 mm and 2.8 to 8.7 mm apart through
    solid metal; the stop must mean what 21.11 says, and it blocked the spraying cube after the speck fix and would
    block the box after 21.10. If wrong (cost): one more diagnosed task.
62. **T7-l.** The parted island (D2: after the speck fix the shoulder stopped on a parted piece holding volume metal,
    3.27e-10 m³ in 96 cells, whose correction could not converge) stops clearly after Task 7b (T7-g), and its remedy
    waits with 21.10 for Asha, since what parts it, neighbouring groups at an edge receding very differently, is what 2
    mm markers do at a sharp edge (inferred in D2, consistent with D1). If wrong: the shoulder stays a strict expected
    failure until then.
63. **(No id; amends T6-e.)** T6-e stays parked, its reason now naming D1's cause: the box-edge probe sits in a column
    that straddles the edge (14.4 per cent low at the point with markers 2 mm apart and 6.0 per cent with them 1 mm
    apart; the column mean 7.29 per cent low, and 6.05 per cent at 1 mm; S12, S14). Updated in Task 7b. The ledger
    records no cost for this ruling; T6-e's stands.
64. **T7-m (closed by T11-a).** The mesh offset is carried to Task 11, which judges the cold plate at a second offset
    (0.27 mm) and records both; reported to Asha. Why: D1 measured once that the plate shifted 0.27 mm receded 0.760 mm
    against 0.800 mm at 1.0 s (S12), and the plate's pass used 4.4 of its 5 per cent. If wrong (cost): one plate run.

## Task 7b: the fixes from D1 and D2

65. **T7b-a (amends T7-i).** Flights take the top flux again at every part under the part's ownership but hold the
    runoff's forcing at the step's start (option A1), and the code says why. Why: T7-i's argument is the heat book (flux
    times area must equal the beam's power on the surface as it is), which the flux alone satisfies; no book or
    diagnosis calls for taking the forcing again, and, measured once, doing so folded a sheet of liquid under the
    spraying cube's bottom face and stopped its first step on the thin-metal check, while the held forcing completed as
    in Task 7. Markers up to 2 mm from φ within a step would own a re-taken forcing, so the fold joins the marker-drift
    question for Asha. If wrong: within a step the runoff is driven by the step start's forcing, as before Task 7b, and
    a real thin sheet under a leeward face could be masked until the forcing is taken again.
66. **T7b-b.** The speck rule extends to the remesh (option B2): before the remesh's volume correction, a piece of φ on
    the new mesh other than the largest, with no marker on it and none of the volume's metal, is lifted out of the
    metal, its metal given to the largest piece, and its nodes held fixed in the correction, by the criterion the moves
    use; the division's stop stays for any piece that still has layer and no column after the settle. Why: the box's
    step-end remesh made a 2.3e-16 m³ markerless speck in the rim, the artefact the speck rule exists for, which should
    neither stop a run nor pass silently. If wrong: a markerless piece made by the remesh is removed rather than
    reported. With it the controller accepted that the check for a returning speck compares against φ before the
    correction, not before the move (D2 wrote the same; the earlier comparison could never see the speck the correction
    re-makes).
67. **T7b-c.** `FILL_RTOL` stays the fill stop's tolerance (option C1), with one condition: if the shoulder's fill stop
    recurred, and its miss was below `FILL_RTOL` of the body's whole layer volume, the stop would gain that absolute
    floor (firing only when the miss exceeds `FILL_RTOL` of both the column and the layer), with tests that such a
    sub-floor miss passes while the existing fill-stop test still raises. Why: the stop exists to catch D2's miss of
    8.7e10 times a column, not 1e-19 m³ in a 3.4e-13 m³ column, and one measurement decides. If wrong: one more shoulder
    run, about 15 minutes. It recurred: the miss was 5.244e-25 m³, 1.539e-12 of the 3.4e-13 m³ column and about 6e5
    times below the layer floor of 3.3e-19 m³ (S14), so the floor was added (bfb6344).
68. **T7b-d.** The shoulder's new stop (a third step of 1012 parts that ended when the sideways solve's GMRES did not
    converge, 6.05e-12 against 1e-12; S14) becomes `SHOULDER_FAILS`' reason and joins the design 21.10 question for Asha
    with D1's 1 mm box, which also split its parts small and stopped on GMRES; it is not diagnosed in stage 3. Why: no
    later task builds on the shoulder passing, and both GMRES stops arose where parts became very small at a sharp edge
    (inferred, not measured). If wrong: a solver weakness unrelated to edges stays unexamined until 21.10.

## Task 7c and diagnosis D3: the thin-metal stop

69. **T7c-a.** The levels' opposed boundaries are no longer candidates for the thin-metal stop; only edges found by the
    node pass are, and only when both the metal along the normals and the sum to the tangent planes confirm the metal
    thinner than 2H at the edge's midpoint. Why: design 21.11's stop means metal thinner than 2H or an edge sharper than
    60 degrees; the levels' points fooled both measures on the box (rim slivers 0.15 to 0.3 mm thick whose undersides
    face back into the body and went to side-face markers 10 mm apart) and on D2's cube (a rounded lip whose normals
    differ by 120 to 124 degrees, over at least 2.49 mm of metal), while the node pass found every real case (the slab,
    the wedge, the fold; S15). If wrong: a thin site that only the levels resolve goes unflagged. The review checked it
    on the merits at the body's 0.7 mm cells: a 1.4 mm plate still gives 82 pairs, a 2.0 mm plate 66, and wedges of 45
    and 55 degrees 9 and 12 (S15).
70. **T7c-b.** The box's mass books, which missed by 1.19e-11 on its second step (`MASS_TOL` 1e-12, growing about
    1.4e-12 per part), are a conservation defect in code Task 7c did not change but unmasked; diagnosis D3 runs in
    scratch beside Task 7c's review, and a Task 7d fixes it inside design 21. Why: the books must close to `MASS_TOL` on
    every completed step, and about 1.4e-12 a part would reach about 1.4e-9 over the shoulder's 1012 parts. If wrong
    (cost): one diagnosis and one task.
71. **T7c-c.** Three findings at the box's edge join the design 21.10 question for Asha and are not fixed in stage 3:
    the box's third step stopping, after 56 parts, on a markerless piece of the rim (4.1e-14 m³, beside the −y edge,
    with one accepted part at the 12-pass cap and an interface mismatch 29 times its tolerance); markerless slivers that
    are not flagged; and the same-side rule giving facets that face back into the body to markers far away (the
    boundaries kept; the exchange likely small, not measured). Why: all three arise at the unresolved sharp edge D1
    diagnosed. If wrong: the ownership defect's exchange across 10 mm stays unmeasured until 21.10.
72. **T7c-d.** Stage 3 does not stop on a sharp concave notch of gas in thick metal; Asha is told it as a reading of
    design 21.11 to confirm. Why: 21.11's "an edge sharper than about 60 degrees" sits beside "columns across an
    ordinary edge meet on the halfway plane", and a notch's columns meet on its bisector with length at every level,
    which the division handles as any shared boundary, and the metal there is not thin; T7-k and the old rule (which
    flagged notches only sporadically) agree. If wrong: a sharp groove runs on without a stop; no recorded case forms
    one.

## Task 7d: the closure's sum and the regrid's ties

73. **T7d-a.** The closure's layer volume per piece is summed with `math.fsum` (D3's option A); the columns' `bincount`
    over a few hundred columns stays, and the other sums over cell volumes in `levelset.py` are listed with whether they
    feed a target or a book, changed only if a measurement shows them leaking. Why: D3 (S16) found that the division
    summed each piece's layer with `np.bincount`, a running sum that over the box's 83,000 to 146,000 equal 0.175 mm
    layer cells rounds every addition the same way and comes out about 1.1e-18 m³ (2.3e-12 of the layer) too large; the
    fill filled to it while the books measured the layer exactly, so each acceptance that moved φ added about 1.3e-12 of
    the mass; with `fsum` (in scratch) the step's residual went from +1.187e-11 to −2.3e-15. Design 21.5 says the
    columns tile the layer exactly, `fsum` is exact whatever the number of cells (which grows with 21.10), and it costs
    under 1 per cent of a division. If wrong: nothing to the physics; the closure factor moves by about 2e-12 of itself.
74. **T7d-b.** Task 7d is sent again, on a stronger model, with a trace first, after the fix made
    `test_a_repartition_by_depth_moves_no_metal_through_the_base_and_keeps_each_depth` fail (a cell boundary moved by
    8.3e-3 of H and a level's metal by 1.6e-2, against 1e-12). The trace finds the comparison in `Skins.regrid` that
    flips, its margin before and after the fix and why the fixture sits on it. Remedy by what it shows: (i) if the
    regrid is not stable under round-off (the cells it makes sit on its own thresholds), make it stable inside `skin.py`
    with a stated basis, such as a margin tied to the thresholds' own scale, with a test that a column perturbed by
    1e-12 of its volume regrids to the same cells; (ii) if only the fixture lands on a threshold by coincidence, move
    the fixture off the tie without changing what it checks. Never loosen `DEPTH_SHARE` or `LEVEL_SHARE`, and never mark
    the test as an expected failure. Why: a change of 1e-13 must not move a cell by 8e-3 of H; the code or the fixture
    is fragile, and the trace says which. If wrong: (i) changes every case's cells at round-off level, which the
    re-gathered record measures anyway. The trace (S17) gave verdict (i): the fill's equal fine cells always put a cell
    exactly on the merge threshold, and the remedy is `REGRID_RTOL` = 1e-9.
75. **T7d-c (amended after Task 7d's review).** The box's slowed third step (not ended after 115 minutes, where before
    it stopped after 56 parts on a markerless piece) and the shoulder's new stop (after 23 parts on `MOVE_SUM_RTOL`,
    where before it ran 1012 parts to the GMRES stop) join the design 21.10 question for Asha, stated as not
    reproducible features of the method, since both stop points moved under a change of about 1e-12 of a column; not
    diagnosed in stage 3. Why: both arise at the unresolved sharp edges D1 diagnosed (the shoulder's message names cells
    no marker may own, the family of T7c-c's markerless pieces), and the smooth plates' recession is unchanged to 1e-14
    under the same change. If wrong: a sensitivity in the move's groups or the part splitting that also affects smooth
    bodies stays unexamined until 21.10. Amended: when it goes to Asha the shoulder's stop is described as a defect in
    φ's move at the sharp edge (its marker groups swept 4.27e-10 m³ where φ's volume changed by 2.67e-10 m³, a 60 per
    cent miss, far beyond round-off), not only as a sensitivity to round-off.
76. **T7d-d.** The box's record entry to 1.0 s is accepted with its note, and no Mac Studio run is made now; Task 11
    runs the box with `--steps 2`, on the laptop or by runbook on the Mac Studio, its third step judged as a strict
    expected failure, and the reference tier on the laptop runs with `-k "not q6p_box"`. Why: the full run's length is
    unknown and would buy a stop point that T7d-c says is sensitive to round-off. If wrong: the box's strict expected
    failure on 73b5a9a is reasoned from its assertions, not measured.
77. **T7d-e (provisional).** The regrid's margin counts a tie as beyond the threshold, so every cell kept lies inside
    the band; accepted for the reviewer to judge on the merits (the review approved). Why: it keeps the existing band
    assertions and moves no judged number. If wrong: a change in cell layout at round-off level in every case, already
    measured by the re-gathered record.

## Task 8: the merge of main (aa0525f)

78. **T8-a.** The skyscrapers reach main's slurry viscosity through `law=` (the shared `slurry_properties` calls the
    same `rheology.law_viscosity` with the same arguments as `skin.viscosity_of`) rather than through a viscosity
    callable. Why: the plan's interface names `law` and `rate_cap` only, and one formula path for both architectures
    keeps the shared module the single place the law lives (Asha, 2026-11-12). If wrong: nothing while the byte-identity
    check against main compares identical (it did, 19 files); a callable can be added if the level set ever needs a
    viscosity the law does not give.
79. **T8-b.** The moved-spray-plan test (`test_moved_spray_plan_pieces_reproduce_the_old_values`) pins its case to the
    bridge law rather than being re-recorded under the composite law. Why: the test proves the move kept the arithmetic
    to the bit, its record is the bridge law's (under the new composite default `rate_slurry[8]` was 7017.63 against the
    recorded 5798.85), and main pinned its own bridge-valued tests the same way (14 of them, 69aba28). If wrong: that
    file has no composite coverage of the shared `slurry_modes`; Task 9's tests cover the composite law on the level
    set's side.
80. **T8-c.** On the review's finding that the brief routed the skyscrapers' slurry viscosity through
    `skin_mod.viscosity_of`: T8-a stands, and a fix round adds a parity test that `spray.slurry_properties(...).mu` is
    array-equal to `skin.viscosity_of(...)` under both laws, and rewords `skinbody.slurry_modes`' docstring. Why: the
    brief's goal (the slurry mode follows `--rheology` exactly as main does) is met and the numbers agree; the test
    turns agreement by reading into agreement by assertion. If wrong: one more test.

## Task 9: the rheology law

81. **T9-a.** The spraying fixture stays pinned to the bridge law in Task 9; the composite law's thin-metal stop on the
    rounded cube (one pair near (3.7, 4.0, −3.2) mm) is carried into Task 11 as a case to diagnose (a real thin site or
    a spurious pair) and onto Asha's list. Why: the two tests check the spraying's books and part count, not the law,
    and the plan gives the cases' re-judging under the composite law to Task 11. If wrong: until Task 11 (and in the
    event, until a remedy for T11-c), no unit test sprays under the default law, and a default-law defect could hide
    behind the pin.

## Task 10: the kinked-profile spraying

82. **T10-a.** The profile path's blanking of the slurry rows' release depth and growth rate is written again in
    `levelset/spraying.py`, with a comment naming the skyscrapers' lines, and pinned by the parity test's equal records;
    moving it into the shared package with the skyscrapers calling it is a follow-up on main. Why: that move needs a
    skyscrapers edit this plan forbade. If wrong: a few duplicated lines until main unifies them.
83. **T10-b.** P5's fixture is recorded on the spraying module's own test skins (committed), and the check on the
    `sprayed` body is a one-off in Task 10's report (its table-path records and released masses byte-identical to two
    runs at the base commit). Why: P5 pins `spraying.py`'s table path, which the module's test skins exercise, while a
    stored record of a whole body would break on any change elsewhere in the body. If wrong: a table-path change that
    only a whole body reaches goes unpinned.
84. **T10-c.** "Equal records" in the profile path's parity tests with the skyscrapers means equal within 1e-9 relative;
    the ledger records the solid fraction's comparison at 1e-8 (the two architectures differ by 3.4e-9 at 907.8 K, each
    from its own property table) and the droplet count's difference as 5.3e-10 against 1e-9. Why: each architecture
    keeps its own property table (skins are never shared), so bitwise equality is not available. If wrong: a difference
    below 1e-9 between the architectures goes unseen.

## Task 11: the cases re-judged on main's physics

85. **T11-a (closes T7-m).** The cold plate's recession does not depend on where its skin base cuts the mesh's cells, to
    within one 20 micrometre skin cell, and the raised plate's failure is its 8.27 mm slab, recorded as measured with
    the 8.27 mm references beside it (S18). Why: the face cannot move among the node planes on this mesh box without
    changing the plate's depth, and the attribution was measured, not inferred: against references for an 8.27 mm slab
    the offset moves the recession by at most 0.9 percentage points. If wrong: a dependence on the cut smaller than one
    skin cell goes unseen; it would be reopened with a mesh box whose planes move instead of the face.
86. **T11-b.** The shoulder's new stop is design 21.11's thin case (S18: metal 1.39 mm thick through the cap's marker,
    under 2H = 2.4 mm; a fin 0.6 to 0.7 mm wide at an edge where the pairs' normals turn 142 and 146 degrees); its
    strict expected failure stays on it, and the shoulder cannot be judged at 2.0 s until the thin case is designed.
    Whether the composite law's runoff cutting a trough 1.2 mm deep in 0.5 s is physical goes to Asha with the 21.10 and
    21.11 questions. Why: the re-run under the bridge law isolates the law as the cause, not whether the result is
    plausible. If wrong: the shoulder waits on the thin case's design when the real defect is the runoff under the
    composite law.
87. **T11-c.** The rounded cube's stop under the composite law is spurious (S18): Task 7c's node-pass mechanism landing
    on its threshold (the pair shares no length at any level, its metal is 3.8 and 8.3 mm thick through the two markers,
    and the thin test's measures read 2.18 and 2.40 mm against 2.4). `division.py` is unchanged, the bridge-law pin
    stays (T9-a), and the fixture's re-worded comment is accepted (its old words, "inferred, not isolated" and "the
    controller's to decide", were stale). The remedy (a node-pass pair must also share a boundary at some level, or be
    neighbours among the markers) is a design decision for Asha. If wrong: no unit test sprays under the default law
    until a remedy lands.
88. **T11-e.** The runs that overlapped one another (scratch runs, the unit tier, and the reference tier beside the
    shoulder and the sphere) are accepted. Why: every judged number is deterministic, and the reference tier reproduced
    the record. If wrong: the wall times in the record overstate the code's cost.
89. **T11-d.** The profile spraying path resolving no column on the cube's first spraying step is the shared rule, not a
    gap in the level set (S19): `melt_profile.reducible` needs the top cell to move, and the cells' speeds at a spraying
    step's start are those the step before left. None of 135 top cells moved before the first spraying step (the melting
    step has no runoff), so every column took the table's path; before a second 0.05 s spraying step all 133 top cells
    moved (the fastest at 1.34 m/s), 49 were reducible and 43 resolved (35 floor-like, 8 fluid-like, the wall at a
    median depth of 0.039 mm), and the 69.9 mg sprayed was all by the gas-driven mode, with books of −1.2e-14 (mass) and
    9.7e-17 (enthalpy) over 10 parts. So the sphere's first step would resolve nothing either, and a measurement of the
    profile path on it needs at least two steps. If wrong: none found; the check is one run on one body.

## Task 12 and the step before the final review: the documents and the merge

90. **T12-a.** The brief's "model_assumptions section 12 (stage 3 as amended)" predated the merge of main, which put
    main's sections 12 to 15 (dispersion, rheology, gas loads, windward loads) at those numbers, and the branch never
    had a stage 3 section; Task 12 adds stage 3 as a new section 16, after main's 15, and sections 10 and 11 point to
    it. Why: renumbering main's sections would conflict with main and break their cross-references. If wrong: one
    heading renumbered later. It was: main's own sections 16 and 17 arrived with the merge of 9993a77 (number 93), and
    stage 3 became section 18.
91. **T12-b.** The record's `cost_rim` and `cost_skip_dev` (stamped e129c95, before design section 21) are measured
    again on this branch's code by the cases' own commands and inserted into the record by key (not by `record`, which
    would drop the record's added blocks), so that design 21.12 cites the current design; the sideways heat's cost is
    measured again like for like in scratch, on one box state, the same skins' heat sub-step with and without the side
    faces (S20), as Task 7's review required; and the sphere's cost is cited only as Task 11's lower bound. Why: the rim
    and skip runs take minutes and the sphere hours. If wrong: about ten minutes of laptop time.
92. **T12-c.** Before the branch is finished, the rulings and the provenance of the scratch measurements the documents
    cite are committed beside the plan (this document; the ledger named it `03b-stage3-rulings.md`), so that every
    citation resolves in the repository; the ledger is deleted when the final review is clean. If wrong: one more
    document to keep in step.
93. **T12-d.** origin/main had moved from aa0525f to 9993a77 (the skyscrapers' skin base cells, 07d45c3 to 6663f4a, and
    the split of a mesh move refused on a fresh mesh, be67361 and 9993a77: the skyscrapers, the command line, their
    tests and documents; nothing in `shared/` or `levelset/`), and main's model_assumptions has its own sections 16 and
    17. After Task 12's review, origin/main is merged into this branch (the house pattern: main merged into the branch
    before the branch goes to main; local, not pushed), the level set's section is renumbered 18 with every reference to
    it, and the level set's unit and reference tiers (`-k "not q6p_box"`) and the shared and command-line tests are run
    again; the final review then covers merge-base(origin/main)..HEAD, the tree that would go to main. If wrong: one
    merge commit to drop if Asha keeps the branch apart. Done as merge 4b38a7d.

## The final review (2026-11-16 to 17)

These were taken after this document was first written (f1c32bd), on the final whole-branch review of
merge-base(origin/main)..HEAD and its one fix wave (1447754, 01c2341, 5bcd629). The review found nothing critical and
three important points (I1 the strict expected failures accepting any failure, I2 documents that turn false on merge,
I3 the sideways solve's hard stop near round-off); its fix wave's scoped re-review found every finding addressed.

94. **PF-a.** The level set's assumptions section is renumbered from 16 to 18 inside the merge commit 4b38a7d, not as a
    commit of its own. Why: the merge's one conflict, in model_assumptions, could not be resolved without choosing the
    section's number, and a tree with two sections 16 would have been wrong; the merge message says so. If wrong: a
    reader of the merge sees the level set's documents renumbered beside main's changes; no code is involved.
95. **FR-a.** The review's M1, making `law` a required argument of the shared `slurry_properties` and `slurry_modes`
    (which default to the bridge law while `rheology.DEFAULT_LAW` is composite), is not done on this branch. Why: a
    required argument breaks any caller on main or on another session's branch that omits it, and the plan fixed the
    default; both architectures pass `law=` explicitly today. If wrong: a future caller that omits `law` flies the
    bridge law silently until the change is made on main (an open question below).
96. **FR-b.** The documents I2 named describe the tree they are in, in sentences true both on this branch and on main
    after a merge (what the level set builds, not where it lives); none says "merged" or names a merge commit before
    Asha decides; todo item 1 is marked built by plan 03b, its merge commit to be added when it lands; the old Mac
    Studio cost runbook is marked superseded and no new one is written. If wrong: one line per document to adjust at
    the merge.
97. **FR-c.** The review's I3 is fixed in the same wave as an inexact Newton step: `SIDE_RTOL` (1e-12) stays the sideways
    GMRES solve's target, a solve that stalls short of it is accepted when its relative residual is at most
    `SIDE_ACCEPT_RTOL` (1e-9) and counted per step (`StepBook.side_solves_at_floor`), and it raises above that, as
    before. The fix first checked in `conduct` that the books cannot depend on the solve's residual (the sideways
    exchange is antisymmetric, the base flux is the residual, and the outer Newton tolerance decides convergence). Why:
    the two recorded stops (D1's box at 1 mm markers, 2.1e-11, and the shoulder's 1012-part third step, 6.05e-12) were
    numerical artefacts that the next planned work, 21.10's finer markers, would meet first. Neither run was repeated.
    If wrong: an accepted solve shifts heat between neighbouring columns by up to 1e-9 of the solve's right-hand side
    within one Newton iteration.
98. **FR-d.** The review's M9, two imports of `skyscrapers/skinbody.py` left unused by Task 8's conflict resolution, is
    fixed only where the branch made the import unused (checked against 9993a77) and nothing uses it by attribute; the
    review's M13 (the step end's carry spread by volume, bounded at 8e-6 of the body) needs no change. If wrong: one
    line more or less in the skyscrapers' diff against main.
99. **FR-e.** The box's strict expected failure (`test_q6p_box`, kept off the laptop by T7d-d) is pinned to its recorded
    failure, the face-centre recession at 1.0 s outside its tolerance, and its judging was checked on the record's
    entry and on stubs, not by running it. The recorded failure at 1.0 s decides the outcome, so a crash or stop in the
    third step or later still counts as the expected failure while the 1.0 s recession stays out of tolerance. If
    wrong: such a later crash of the box, in its full-length run on the Mac Studio, is reported as the expected failure.
100. **FR-f.** The fix wave's three minor leftovers, found by its re-review, are left: the "books close" assertion in the
    test of an accepted solve holds by construction (the test's substantive checks are the count and the temperatures
    within `NEWTON_TOL` of the converged run's); a sentence in `skin.py`'s module docstring wraps raggedly; and the box's
    judging decides on the 1.0 s row first (FR-e). If wrong: a weaker test of the accepted solve than its name suggests.

## Where the scratch measurements came from

Each paragraph says what was run, on which commit, how many times, and what it gave. "Once" means one run. The scratch
scripts and logs named here lived in the build's scratch folder and will not survive; the results are stated in full
here. Loads are one-minute load averages on the laptop's 8 cores. Measurements made by committed tests or recorded in
`analysis/levelset/stage3_recorded.json` are not repeated here, since those can be re-run.

- **S1. Task 1's division measurements** (cited for T1-b and T1-e, model_assumptions 18.2, design 21.12 and the
  docstring of `division.py`, as "Task 1's report" or "Task 1's fix round"). In the first round (code committed as
  c6d5c79), the division was run on the stage 3 box with markers seeded 2 mm apart at 0.7 and 0.175 mm cells with 24, 48
  and 96 levels, once each, timing `divide` (the factors and times are in number 7); the same box at 0.7 mm gave the
  difference between the division's top area and the owned areas (up to 10.4 per cent, 95th percentile 8.6, median
  0.20), and the level-by-level jumps of 12 to 44 per cent in edge columns' coarea (T1-h), once. In the fix round
  (51be5e0) the closure factors were measured once each: the review's 6 mm cube 1 + 6.3e-4; the box 1 − 3.7e-4 at 0.7 mm
  (the division taking 1.6 s) and 1 − 3.6e-4 at 0.175 mm (17.6 s); the slab exact; a sphere 1 − 2.8e-5; two cubes 1 +
  6.0e-5 and 1 + 7.3e-4. The flat cells' closure shortfall of 2.13 per cent on a 6 mm cube was the review's measurement
  on c6d5c79.
- **S2. Task 2's moves** (cited for T2-b to T2-d and model_assumptions 18.3, as "plan 03b Task 2's report, section 4").
  The four sharpness tests' values in 18.3 (0.279 mm³ against 1.889 mm³, the edge at 0.68 against 0.20 mm, the deposit
  and the L) are the committed tests' printed values on 770259e, each run once. The comparisons behind the rulings were
  made once each in a prototype in scratch before the code went into the package: whole cells against the exact split
  (−21 to +14 per cent against −0.73 to +0.52 per cent of the division's columns), the secant (a miss of 1e-2 after 8
  rounds) and the plain reset (0.26 mm³ of edge eroded a call; a converged miss raised from 8.7e-7 to 7.4e-3). At that
  commit the plain reset's β held 321.830 mm³ against the grid reference's 322.423 mm³ (0.593 mm³ between); after T5-c
  it holds 322.500 mm³ (S20).
- **S3. Task 3's profile on the box** (cited for T3-c and model_assumptions 18.2). The skins were given the box's real
  division (24 levels, the division taking 8.0 to 8.3 s, factor 1 − 5.3e-4) on the code committed as a146971, once: the
  fill exact to 2.2e-16, 28 of 5,406 neighbouring conductance pairs differing by more than a factor of two, 195
  reversals of the profile over 132 columns, the largest ratio 2.85 after smoothing over one mesh cell; and a stress run
  of 30 melting sub-steps, in which every column overfilled by swelling by up to 14.5 per cent and the base exchange
  returned the excess to 3.2e-16.
- **S4. T5-b and T5-l, edge rounding under repeated recession** (cited for T5-b, T5-l, T5-m, model_assumptions 18.3 and
  18.4, the latter as "plan 03b Task 5's report, section 12"). T5-b: twelve cold 0.1 mm recession steps of the test cube
  (each taking the top 0.1 mm of every column, then one cold 0.05 s step), on 3d376c0, once, in scratch (`bumps12.py`),
  with about six cores busy with other work: the gap to the shrunk cube sampled on the nodes grew from 0.236 to 8.145
  mm³ in bursts about 0.7 mm apart, the last four steps +1.159, +0.761, +0.156 and −0.186 mm³. T5-l: on 718cf7e, in
  scratch, about 40 minutes, every run once: (A) the same twelve steps with exact asks at 850 K (where the Scheil alloy
  is 13.5 per cent liquid) still grew +7.75 mm³; the measure itself swung by up to 2.2 mm³ with a 0.7 mm period under an
  exact uniform offset; the bisection named the fill's density carry, drawn metal swelling about 2.3 per cent (2,749
  against the volume's 2,813 kg/m³) and moving φ out by 6.9 mm³ over twelve steps (a uniform 15 micrometres plus 3.8 to
  4.1 mm³ of shape); at 700 K (solid, no carry) exact asks tracked the exact recession to 0.017 mm³ while the division's
  asks departed by 0.11 mm³ after one step, 2.0 by step 5 and 2.6 by step 12; (B) the ratios of ask to exact share were
  0.9999 to 1.0036 on faces (worst +9.2 per cent), 0.9993 to 1.0016 on edges (worst +9.6 per cent) and 0.83 to 0.999 at
  corners (worst −20 per cent), the misplacement summing to at most 0.44 mm³ a step.
- **S5. T5-d, the cap on the draw's temperature** (cited for T5-d and model_assumptions 18.4). The test cube at 6 MW/m²
  on its +x face with the liquid fed off, run from one start with and without the cap, once each: first on 3d376c0
  (superseded), then on the fix round's code before it was committed as 718cf7e: one 0.5 s step took 9 parts capped and
  10 uncapped, two 0.25 s steps 3 and 3, then 8 and 8; mass residuals at most 3.4e-13, enthalpy 1.2e-16. Every 0.5 s run
  and the second 0.25 s step still stopped on the thin-metal check, with or without the cap.
- **S6. T5-i, the density carry** (cited for T5-i and model_assumptions 18.4). The heated cube (2 MW/m² on its +x face,
  three 0.5 s steps), on the fix round's code (committed as 718cf7e), once before and once after the cap's removal: the
  carry left at the three step ends was −8.10e-12, −6.13e-12 and −1.36e-12 m³ (8.0e-6, 6.0e-6 and 1.3e-6 of the cube's
  volume of about 1,017 mm³), and the largest row's share of its column 1.29e-3, 1.42e-3 and 1.00e-3; the run after the
  cap's removal gave the same largest values.
- **S7. T5-o, swelling and the thin-metal stop** (behind T5-h and T5-n). From a detached scratch worktree at 718cf7e,
  once each, with `SkinSettings(swell=False)` against the committed settings: the spraying fixture's step still failed
  without swelling (in the runoff's transport solve, residual 7.5e-13 against 1e-13, after 5 parts; its face rising to
  +1.87 mm against +1.24 mm with swelling, where the committed run stopped on the thin-metal check with 8 pairs), and
  the 6 MW/m² 0.5 s step stopped on 7 pairs without swelling against 20 with it (highest point +0.079 against +0.145
  mm). Mass residuals at most 4.9e-13, enthalpy 1.7e-16. Conclusion: the folds come from the runoff, not the swelling.
- **S8. T6-a and Task 6's fold** (cited for T6-a, T6-d and model_assumptions 18.6). On Task 6's first-round code
  (0134911), from one saved melted sharp cube, the 0.05 s continuum step was run once in each way: one part (3 opposed
  pairs at the stop, ridge about 1.0 mm), two parts (the same), four parts (ridge 1.6 mm), and with surface tension off
  (even recession, median 0.65 mm, one pair at a receding corner's rim); the push terms' medians on the face were
  capillary 2.5e5, pressure gradient 3.8e4, deceleration 126 and thickness 18 N/m³. Then, in a scratch worktree of
  0134911, once each: (i) eight parts, still stopping on 1 pair, ridge +1.53 mm and deepest corner trough −1.04 mm
  against +1.61 and −0.99 at four; (ii) a cube with 2 mm rounded edges, melted the same way, surface tension on, about
  100 s: no stop, 5 parts, mass residual −1.7e-13, the face from +1.18 to −1.23 mm (median −0.75), medians capillary
  1.18e5 and pressure gradient 2.78e5 N/m³.
- **S9. Task 6's box edge and sideways heat** (cited for T6-b, T6-e and model_assumptions 18.5). The box-edge reference
  test (marked `reference`), once on Task 6's first-round code (0134911) and once after T6-b's change (df971a2): the
  probe's column rose 25.73 K and then 26.30 K against the reference's 28.37 K column mean (30.71 K at the probe point);
  the centroid ratios (0.56 across the edge at 0.5 mm depth, 0.47 at 0.9 mm, 0.60 to 1.84 along a face) from the
  division of the stage 3 box, once; the lateral step test's largest miss 0.079 K of a 100 K step; the box's interface
  mismatch with sideways heat 7.26e-4 (7.85e-4 before).
- **S10. T7-c, exact recession on the box's own mesh** (cited for T7-c, model_assumptions 18.3 and design 21.8 item 4).
  On Task 7's code (the report names no commit; Task 7 committed 721efce and 71ef69d), in scratch (`t7/t5m.py`, and
  `t7/t5m_split.py` with the remesh split out, which reproduced every number), once: the box as the case builds it
  (0.175 mm cells at every edge), at 700 K, twelve steps each taking the top 0.1 mm of every column and one cold 0.05 s
  step, measured against φ0 + kD on the rebuilt nodes. The departure after 1.2 mm was 0.784 mm³, having risen to 0.71 by
  step 7 and then swung between 0.65 and 0.79 (+0.003 mm³ over the last four steps); the criterion at 0.175 mm cells was
  0.547 mm³ over the receded box's 83.2 mm of edge and 0.736 over the starting box's 112 mm; the remesh changed the
  departure by at most 0.0013 mm³; the closure factor departed by 1.14e-3 at step 5.
- **S11. Task 7's box with surface tension zero** (behind T7-a). On Task 7's code, the box's second step from its kept
  starting state, once as the case runs it and once with the liquid's surface tension zero: both stopped on the same
  edges, the face the same within 0.04 mm.
- **S12. Diagnosis D1, the box's uneven face** (cited for T7-i, T7-j, T7-m, T6-e and model_assumptions 18.5, 18.9 and
  18.11). On 71ef69d (a commit made meanwhile, ed30160, changed only the cases' profiler and its test, which the runs
  did not use), every run once, one at a time, on mains power, about 50 minutes in all, re-running the box's second step
  from Task 7's kept state at 0.5 s (which reproduced Task 7's step exactly: centre 0.907 mm, 22 pairs). Found: 20
  side-face markers own the heated face's outer strip (32.7 of 100 mm²), their columns straddling the edge with the
  strip 24 to 38 per cent of their tops, which peak at 906.4 K against T_feed 910 K and remove nothing; the unspent heat
  raises the inner columns' base temperatures 2.5 to 3.4 K above the plate's (up to 10 K under the outer face columns),
  so the inner columns recede 0.915 to 0.938 mm against the plate's 0.800. The plate's column with 9 markers 2 mm apart
  and no edge recedes 0.798 to 0.800 mm at five of six points read and 0.826 mm at the sixth, as the plate. Heat taken
  again under one ownership each part: heat in +0.5 per cent of the beam's instead of +5.5, no trough, centre 0.888 mm,
  4 pairs. Markers 1 mm apart (499 against 132): rim, lips and trough gone and the centre within about 4 per cent of the
  plate, but the step stopped at 0.979 s on the sideways GMRES (2.1e-11 against 1e-12), its parts split to 5.2 ms. The
  box-edge probe: 14.4 per cent low at the point at 2 mm spacing and 6.0 per cent at 1 mm; the column mean 7.29 and 6.05
  per cent. The cold plate shifted 0.27 mm receded 0.760 against 0.800 mm at 1.0 s (later attributed to its deeper slab,
  S18). Ruled out: the runoff, sideways heat (off is worse), β's cut of the volume's cells, the first step, a tighter
  same-side rule.
- **S13. Diagnosis D2, the hang and the shoulder's move-sum stop** (cited for T7-b, T7-f to T7-h, T7-l, T7-k and
  model_assumptions 18.2, 18.4 and 18.11). On 71ef69d (ed30160 as for D1), from Task 7's kept states, every run once,
  about 50 minutes in all, each symptom reproduced to Task 7's digits. Found: the volume correction shifts nodes near a
  lifted speck by more than its 1 nm lift (0.24 micrometres on the cube, 3.1 on the shoulder), so the speck returns; on
  the cube it lay nearer than the body to two markers that had drifted 1.67 and 2.08 mm into the gas, which the division
  put on it (top 2.4e-16 m²), so the fill skipped them and the split asked for 4.3e9 parts; on the shoulder a returned
  markerless speck was swallowed at move 22, ten cells changing 6.507e-12 m³ in no group. Lifting again after the
  correction (in scratch) cleared both, after which the cube stopped on thin-metal pairs 2.8 to 8.7 mm apart through
  solid metal and the shoulder on a parted island (3.27e-10 m³, 96 cells). The shoulder's closure at 24, 48 and 96
  levels: 3.67e-3, 2.24e-3, 1.13e-3 as built (451 nodes at φ = 0 on its faces), 3.7e-5 and 4.5e-5 at 24 and 96 levels
  shifted 0.13 mm; the rounded cube at most 2e-5. Within one 50 ms step φ moved up to 2.08 mm (1.7 H) from markers that
  stay where the step began.
- **S14. Task 7b's re-measurements** (cited for T7b-c, T7b-d, the amended T6-e and model_assumptions 18.2, 18.4, 18.5).
  On Task 7b's working tree (ed30160 with its uncommitted fixes), D1's probe script re-measured the box-edge test, once
  each: 26.300 K against 28.367 K (7.29 per cent low) with markers 2 mm apart, 28.870 K against 30.730 K (6.05 per cent)
  at 1 mm, equal to D1's. The shoulder built 0.13 mm off the grid closed within 4.3e-4. The shoulder's fill-stop miss,
  5.244e-25 m³ in a 3.4e-13 m³ column, was read from its stopped run under T7b-a, once, before bfb6344. The shoulder's
  1012-part third step, which stopped on the sideways GMRES (6.05e-12 against 1e-12) after 4230 s with the division 54
  per cent of it, came from the reference run on bfb6344, once, sharing the laptop with a second shoulder run.
- **S15. Task 7c's diagnosis and review checks** (cited for T7c-a, the docstring of `division.py` and model_assumptions
  18.7 and 18.11). On fc548e4, in scratch (`t7c/repro.py`, `t7c/diag.py`), once each on the stopped bodies: the markers
  stood on φ at the step end to 1.4e-14 m; the box's four pairs came from the levels on rim slivers 0.15 to 0.3 mm thick
  whose undersides faced side-face markers 10 mm away; the cube's two from a rounded lip with normals 120 to 124 degrees
  apart and at least 2.49 mm of metal, reached through 1.2 mm mesh edges that jump a column. The review's check at the
  body's 0.7 mm cells, on the fixed code (92f4424), once: opposed pairs on a 1.4 mm plate 85 before and 82 after, a 2.0
  mm plate 85 and 66, a 45 degree wedge 14 and 9, a 55 degree wedge 15 and 12.
- **S16. Diagnosis D3, the box's mass books** (cited for T7c-b, T7d-a and model_assumptions 18.2). On 47dedb9, with
  ed30160 and fc548e4 checked out in detached scratch worktrees for comparison, about 27 minutes, every run once: the
  box's second step's miss is identical part by part on all three commits (its cause dates from Task 1's c6d5c79); it
  follows the number of full fine layer cells; `math.fsum` per piece took the step's residual from +1.187e-11 to
  −2.3e-15 over the same 8 parts; the proposed test was red on 47dedb9 (a 6 mm cube 1.0e-13 of the layer, the box
  fixture 3.9e-14, against a bound of 1e-14) and green with `fsum`.
- **S17. Task 7d's trace and ties probe** (cited for T7d-b, model_assumptions 18.4 and 18.10 and design 21.12). On
  47dedb9 with the `fsum` line added, in scratch, once each: the trace (`t7d/trace.py`, with the committed and the fixed
  closure) found the flipping comparison, the merge test on a 20 micrometre cell at 0.6 mm depth exactly half its 40
  micrometre target, its margin +2.44e-15 of the threshold after the test's regrid and −3.44e-15 in one column after the
  settle with `fsum`; with the committed closure the settle's regrid already changed 16 of 109 columns, none of them
  checked by the test. The ties probe (`t7d/ties.py`) filled single 1.2 mm columns, flat and narrowing, at both presets:
  every column filled to exactly H had a cell on the merge threshold to within 8.7e-15, and scaling its metal by 1 +
  1e-12 changed its cells from 43 to 42 (development) and 72 to 71 (production). The repartition test's numbers in 18.4
  (depths within 4.5e-16 of H, metal 1.8e-16, books −2.2e-16 and 2.2e-16) are the committed test's on 73b5a9a, once.
- **S18. Task 11's runs** (cited for T11-a to T11-c, the docstring clause of `stage3_cases.py` on T11-a, and
  model_assumptions 18.7, 18.9 to 18.11, the latter as "Task 11's report, section 4.4"). All on e57945f, once each, with
  load averages of about 6 to 28, higher while they overlapped: the box and the shoulder for two steps each under
  `SkinSettings(rheology="bridge")`, which reproduced Task 7d's numbers apart from wall times; `analysis/ablation_1d.py`
  (the references' own code, 1 micrometre cells, 0.5 ms steps) for 8 mm, reproducing both committed references exactly,
  and for 8.27 mm, receding 0.755, 2.623 and 5.044 mm with removal at T_feed and 0.803, 2.683 and 5.120 mm at the
  liquidus; the shoulder's stopped step re-run from its kept state with probes, twice (metal 1.39 mm through the cap's
  marker along its normal; the node-pass edges' measures 0.84 to 1.77 mm and 1.39 to 2.06 mm; a map of φ on a 0.1 mm
  grid across the edge showing a trough about 1.2 mm deep, a fin 0.6 to 0.7 mm wide and a run-over about 0.4 mm thick);
  and the rounded cube under each law (the composite law's pair 4.17 mm apart, 3.83 and 8.30 mm of metal through its
  markers, zero shared length at all 25 levels, measures 2.18 and 2.40 mm; under the bridge law the step completed, 5
  parts, 42.7 mg sprayed). The 50 mm sphere (`cost --what sphere --preset dev --steps 2`) built in 2137 s and its first
  step had not finished after 39 minutes, once.
- **S19. T11-d, the cube through two spraying steps on the profile path** (cited for T11-d, model_assumptions 18.8,
  design 21.8 item 11 and parity ledger row 12). The controller's check (`t11/ctl/cube_two_sprays.py`), on 216296b,
  under the bridge law, once, on 2026-11-15: the rounded cube melted 0.25 s with no runoff, then two 0.05 s spraying
  steps on the profile path, counting before each step the skins whose top cell moves and those reducible; the numbers
  are in number 89.
- **S20. Task 12's probes** (cited for T12-b and model_assumptions 18.1, 18.2, 18.6 and 18.10, design 21.8 items 1, 5, 8
  and 10 and 21.12, as "Task 12's scratch probe" or "Task 12's run"). On 216296b, with load averages of 3.4 to 3.7:
  `side_heat.py`, once (23:36 to 23:38 on 2026-11-15), built the box as the case does (132 markers, 327,182 cells) and
  advanced its first step with the case's own function, wrapping `Skins.conduct`: at the middle heat sub-step of the
  first real pass (132 columns of 42 cells, 403 shared boundaries, cells at 762 to 899 K) the same call took a median of
  0.0278 s with side faces and 0.0024 s without (five each, interleaved, the cells restored before each), 11.6 times;
  over the step's 600 heat sub-steps with side faces, 20.3 s against 1.6 s for the same calls without, 18.8 s or 54 per
  cent of the step's 34.9 s (median ratio 12.4, 10th to 90th percentile 10.2 to 16.4); the step's 34.9 s is the quiet
  re-timing of the box's first step. The item tests were run once with their printed output (31 passed in 65 s): the
  skin-base test's 322.500 mm³ against 322.423 mm³ and 0.077 mm³ between; the box's test-mesh factor 1 − 3.30e-4 and a
  sphere's 1 − 2.84e-5; the cold box's φ moved 0.0 m; the runoff test's 1.818e-7 kg. The re-measured `cost_rim` and
  `cost_skip_dev` are in the record (load averages about 3.1 to 3.4).

## Review findings and reports that the code cites

The package's docstrings, the tests' docstrings and reasons, and the record cite some findings of the task reviews and
some task reports by name ("Task 1's report", "stage 3 Task 11's review, M6"). Those reviews and reports are working
files kept outside the repository (`.superpowers/sdd/03-stage3-skins/` for plan 03, the first stage 3 plan, and
`.superpowers/sdd/03b-stage3-one-skin-depth/` for plan 03b), and they are deleted with the ledger. Each cited item is
listed here once, numbered F1 onward, with what it found; a citation in the code is followed by "cited finding F<n>".
"Stage 3 Task N" without a plan means plan 03's task, as the code wrote it before plan 03b; "plan 03b Task N" is this
plan's. The numbers are the reviews' and reports' own, each measured once unless said otherwise.

*Plan 03 (the first stage 3 plan, two skin modes).*

- **F1. Plan 03 Task 9, fix round 2, the ruling on its concern 1.** A marker's base law was concentrated on its nearest
  patches of β, and the pending heat on one patch; the fix spreads both over a group of patches whose area is within
  `GROUP_AREA_RANGE` (0.5 to 2.0) of the skin's base area, which raised the full-size cube's coldest volume temperature
  after its first step from 545.8 K to 830.8 K. Cited by
  `test_a_marker_whose_nearest_patches_are_few_gets_a_grown_group_and_the_law_stays_exact`.
- **F2. Plan 03 Task 11's review, C1 (critical).** The skins' areas (the tessellated, globally scaled Voronoi cells) and
  the ownership the loads and the base groups used were two partitions of φ: the q6p box's markers took 0.9606 of a
  beam's projected heat at the start and 0.6398 after the first melting step by their areas, against exactly 1 by the
  areas they owned, the ratio of the two running from 0.20 to 1.62. Fixed in the fix round 1 by one partition, each
  marker's area the area of φ it owns (`markers.ownership`). Cited in `loads.py`, `markers.py` and
  `test_levelset_loads.py`.
- **F3. Plan 03 Task 11's review, I2.** The shoulder's 4 mm freeze-reach criterion (the implementer's own) failed on
  liquid in transit 3 K above the rigid temperature, not on absorbed base metal as `SHOULDER_FAILS` then said (deposit
  20.05 mg, hot metal 4.30 mm past the heat's end). The controller's ruling on it judges the deposit at 2.0 s instead:
  at least `DEPOSIT_RIGID_SHARE` (0.9) of it rigid and at most the skyscrapers' 3.0 per cent past 5 mm, the reach
  reported. Cited by `stage3_cases.py` and `test_levelset_stage3_cases.py` ("the controller's ruling on the review's
  I2").
- **F4. Plan 03 Task 11's review, M1.** The volume solve was not reproducible run to run, because the algebraic
  multigrid preconditioner draws a random start vector from numpy's global state; it is now drawn under `AMG_SEED` and
  the state put back. Cited in `cutthermal.py` and `test_levelset_cutthermal.py`.
- **F5. Plan 03 Task 11's review, M3.** The box's heat-coverage unit test could not see C1 (F2): on its fresh box the
  markers took 1.013 of the projected heat inside a band of 0.95 to 1.053; with one partition it asserts exactness
  (`COVERAGE_EXACT`). Cited in `test_levelset_loads.py`.
- **F6. Plan 03 Task 11's review, M5.** The committed record mixed runs from before and after the fix round with nothing
  to tell them apart; each case's summary now carries the commit it ran on and whether the code differed from it
  (`code_stamp`). Cited in `stage3_cases.py` and `test_levelset_stage3_cases.py`.
- **F7. Plan 03 Task 11's review, M6.** A failing case left nothing to diagnose (the states it named were gone, so the
  box's "volume correction of piece 1 did not converge" at 2.0 s could not be diagnosed); with `STAGE3_STATE_DIR` set
  the case script now pickles the body before each step and keeps it when the step raises, and `correct_volumes` reports
  the piece's target, bracket, iterations and claimed cells when it raises. Cited in `levelset.py`, `stage3_cases.py`,
  `test_levelset_levelset.py` and `test_levelset_stage3_cases.py`.
- **F8. Plan 03 Task 11's review, M7.** The seed for an island of new metal was evaluated at every node outside the old
  metal within the fine margin on every base move, though only island nodes read it (base moves were 48 per cent of
  the box's step to 1.5 s); it is now evaluated only at the nodes the extension cannot reach. Cited in `cutthermal.py`.
- **F9. Plan 03 Task 11's fix, item 1.** A skin's base group is the base area beneath its own cell on φ, by the
  ownership the loads use; before, by nearest base point, the box's face skins' groups were 0.55 to 0.82 of their base
  areas. `GROUP_AREA_RANGE` stays as a guard, and the most groups it replaced in a part is recorded
  (`groups_guarded_max`). Cited in `skinbody.py` and `stage3_cases.py`.
- **F10. Plan 03 Task 11's fix, item 2.** The gas loads are each marker's area-weighted mean over the patches of φ it
  owns, not its own normal's, by which the box's markers caught 74 per cent of the projected heat. Cited in `loads.py`
  and `test_levelset_loads.py`.
- **F11. Plan 03 Task 11's fix, item 3.** An island of new metal that a base move makes apart from the old metal (on the
  box at 1.5 s, 24 nodes about 1.1 mm from it) takes the handed metal's temperatures, where before the move raised;
  the nodes so filled are counted (`island_nodes`). Cited in `skinbody.py`'s `StepBook` and `stage3_cases.py`.
- **F12. Plan 03 Task 11's fix, item 4.** The hot plate's volume beneath the skin was judged at most 1 K
  (`HOT_BASE_TOL_K`) above the rigid temperature, the skyscrapers' criterion (it stood 4.7 to 5.2 K above before);
  plan 03b reports it and no longer judges it, because the one skin depth puts the melting front below the base. Cited
  in `stage3_cases.py` and `test_levelset_stage3_cases.py`.
- **F13. Plan 03 Task 11's fix round 1, the review's concern (b).** The box's "volume correction of piece 1 did not
  converge" at 2.0 s, left undiagnosed for want of the state (F7); the fix round made the correction converge in one
  case that can produce it, a speck of metal that a base move makes beside the body, which must grow some 1e5 times and
  whose first Newton step overshot into saturation. Cited in `test_levelset_levelset.py`.
- **F14. Plan 03 Task 11's re-review, minor item 4.** The code stamp's difference check covered only `reentry_model` and
  `analysis/levelset`, so a changed material table or reference value would not mark a run; it now covers `data`
  (`STAMP_PATHS`). Cited in `stage3_cases.py` and `test_levelset_stage3_cases.py`.

*Plan 03b (this plan).*

- **F15. Task 1's report.** With each level triangle's area divided by |∇φ| (ruling T1-a) the columns' shortfall against
  the layer at sharp edges is gone, and what remains of the closure's departure is the levels' sampling of the cells
  along the edges (at 24, 48 and 96 levels the box's factor was −2.04e-3, +3.32e-3 and +5.60e-4 from one at 0.7 mm
  cells, before the fix round below); the level surfaces' own areas fall smoothly but for rises of 4e-4 of the top near
  the base in two columns. Cited in `division.py` and `test_levelset_division.py`.
- **F16. Task 1's review, Important 1, "the review's measurement".** P1 cells along an edge whose depth is constant are
  cut by no level: on a 6 mm cube so placed, 474 cells constant to 0.1 µm held 2.0 per cent of the layer and the levels
  alone left the columns 2.1 per cent short (96 levels), whatever the number of levels. The fix round gives such a
  cell's exact volume whole to the column owning its centroid's foot. Cited in `division.py` and
  `test_levelset_division.py` ("the review's Important 1").
- **F17. Task 1's review, Important 2, and its case for the fix round.** A column spanning two pieces of metal was
  counted wholly in its majority piece: of two 6 mm cubes 0.7 mm apart across gas, markers of each owned 0.39 per cent
  of the other's surface and the pieces' factors moved about 2e-3 apart. Ruling T1-i divides ownership piece by piece.
  Cited in `division.py` ("the review's measurement") and `test_levelset_division.py`.
- **F18. Task 2's report.** The plain distance reset moves the zero contour where it bends within a cell, inward at β's
  convex edges: on the stage 3 box β held 0.67 mm³ less than φ + depth. Cited in `surfaces.py`.
- **F19. Task 3's review.** Two deferred minors cited: no test ran the momentum step under a measured profile (one now
  does), and combining two columns by depth moves the keeper's cells by up to 4.9 per cent where swelling metals at
  different temperatures mix (their volumes do not add), which the fill after a combine undoes. Cited in `skin.py` and
  `test_levelset_skin.py`.
- **F20. Task 4's review.** Over long trial passes the sideways solve's worst residual at 10 µm cells was 9.9e-13, only
  a few times under `SIDE_RTOL` (1e-12), and a solve short of it raised with no fallback; the reviewer found that the
  books close whatever the solve's tolerance. The final review's I3 (F33) acted on it. Cited in `skin.py`.
- **F21. Task 6's diagnosis of the spraying fixture.** On the sharp cube the capillary push of φ's recovered curvature
  at the edges (some 2.5e5 N/m³ on the face markers) folded the liquid face within the step and stopped it on
  `THIN_METAL`; hence the spraying fixture's 2 mm rounded edges (ruling T6-d) and the edge-runoff test's surface tension
  off. Cited in `test_levelset_skinbody.py`.
- **F22. The reviews' note on `Markers.release` and Task 6's minor item.** A released marker's keeper was chosen by
  position alone, ignoring the same-side and same-piece rules (Task 5's review, carried into Task 6 as optional; now
  applied); and, from Task 6's review, the book's `area_repartitioned_share` was renamed `volume_repartitioned_share`
  in Task 7, a book pickled under the old name loading under the new one. Cited in `test_levelset_markers.py` and
  `test_levelset_skinbody.py`.
- **F23. Task 7's review, I1, and the controller's carry 5.** The solver watch did not measure the interface mismatch on
  `INTERFACE_TOLERANCE`'s own scale (sum |mismatch| over sum |E_top| + |E_base| per accepted part): the recorded share
  was the last part's over the step's heat in, so the 1.46e-3 then reported compared two different quantities; the
  watch's `interface_mismatch_max` is now the per-part measure. Cited in `test_levelset_stage3_cases.py`.
- **F24. Task 7's review, M2.** `BOX_FAILS` cited numbers measured only in scratch (the stopped step's face and edges);
  `run_steps`' `failed` now records the case's own measures of a stopped body. Cited in `stage3_cases.py` and
  `test_levelset_stage3_cases.py` ("Task 7's minor M2").
- **F25. Task 7's review, M5.** The case tests did not apply planning ruling 1's closure check; `_books` now asserts
  each step's `layer_factor_max` within `DIVISION_TOL`. Cited in `test_levelset_stage3_cases.py`.
- **F26. The controller's carry 4 into Task 7**, as the code cites it: the division's own cost line in each step's row
  (every division, the tie rule's second nearest-marker query and the sideways conductances' centroid distances:
  `division_cost`). Cited in `stage3_cases.py`.
- **F27. Task 7c's diagnosis.** The spurious thin-metal stops: on the box's second step all four marked pairs came from
  the levels, slivers of metal 0.15 to 0.3 mm thick off the receded face's rim given to markers on opposite side faces,
  10 mm apart through solid metal but holding only 0.24 mm along their normals; hence the two measures of the thin
  test, at the node pass's edges alone (ruling T7-k). Cited in `division.py`.
- **F28. Task 7d's report, round 2** (the record's `box_dev.record_note`, "task-7d-report.md, round 2"). On Task 7d's
  code (73b5a9a) the box's third step had not ended after 115 minutes when its run was stopped at the two-hour limit,
  and the shoulder's third step stopped after 23 parts on the move's volume-sum check; the record's box entry is
  therefore a run to 1.0 s (ruling T7d-d).
- **F29. Plan 03b Task 11's report** (the record's shoulder entry, "task-11-report.md"). Under the composite law the
  shoulder's runoff cuts a trough about 1.2 mm deep in the heated cap and leaves a fin 0.6 to 0.7 mm wide, metal 1.39
  mm through the cap's marker, a real thin site (ruling T11-b); and the cold plate raised 0.27 mm (8.27 mm of metal) is
  recorded, its slower recession that of the deeper slab (ruling T11-a). Cited in `test_levelset_stage3_cases.py`.
- **F30. Plan 03b Task 11's review, M1.** `record()` rebuilt the whole record from the output files, which would erase
  the blocks added by hand (`against_task7`, `against_offset_0`, `record_note`) and Task 12's cost entries; the final
  review's M5 (F35) made it merge. Cited in `stage3_cases.py` and `test_levelset_stage3_cases.py`.
- **F31. Plan 03b Task 11's review, M5.** `--offset-mm` was silently ignored for the box and the shoulder; it is now
  refused for any case but the plates (F35). Cited in `stage3_cases.py` and `test_levelset_stage3_cases.py`.

*The final review of plan 03b (`final-review.md`, 9993a77..f1c32bd).*

- **F32. I1.** The three strict expected failures (the box, the shoulder, the box edge's heat) passed for any failure,
  a crash included, because the case script catches every exception into its failure record; each is now marked
  `raises=AssertionError` and asserts only its recorded failure, anything else raising `NotTheRecordedFailure`. Cited
  in `tests/levelset/cases.py`, `test_levelset_stage3_cases.py` and `test_levelset_skinbody.py`.
- **F33. I3.** The sideways solve's 1e-12 tolerance sat at round-off as a hard stop, where the shoulder's 1012-part
  third step under the bridge law (6.05e-12) and the box at 1 mm markers (2.1e-11) stopped; a stall at most
  `SIDE_ACCEPT_RTOL` (1e-9) is now accepted and counted. Cited in `skin.py`, `skinbody.py`, `stage3_cases.py` and
  `test_levelset_skin.py`.
- **F34. M6 and M7.** The `PARTS_MAX` stop did not mark its part as restored, so a body with no rigid metal left
  reported the part count instead of `NO_CORE`; and `Skins._append_bottom` lacked the trial-pass guard of the method it
  replaced. Cited in `skinbody.py`, `skin.py` and `test_levelset_skinbody.py`.
- **F35. M5.** The cases script's two hazards, F30 and F31, fixed. Cited in `stage3_cases.py` and
  `test_levelset_stage3_cases.py`.
- **F36. M11 and M12.** The heated cube's first step was judged against the face's areas after its third step; and three
  behaviours had no direct test (the thin test's rejecting half, the profile spraying path with reordered skin rows,
  `ownership` given one of its two normals). Cited in `test_levelset_skinbody.py`, `test_levelset_division.py`,
  `test_levelset_spraying.py` and `markers.py`.

## Open questions for Asha

These are not decided. They are the ones model_assumptions section 18.11 states, with the same numbers, followed by
those the ledger names that section 18.11 does not.

- **The order of design 21.10 (finer markers near edges) against 21.9 (the volume melt flow), and what the box's rim and
  the shoulder need** (T7-j, T7-l, T7b-d, T7c-c, T7d-c). The plan keeps the design's order, 21.9 first. The box's
  unreceding strip and the box-edge test's straddling column were diagnosed as 2 mm markers not resolving a sharp heated
  edge (D1); the markerless rim pieces, the parted island and the third steps' stops arise at the same edges and are
  attributed to that cause by inference. The box's third-step stop point moved under a change of about 1e-12 of a
  column, and the shoulder's third-step stop is a defect in φ's move at the sharp edge (its groups swept 4.27e-10 m³
  where φ changed by 2.67e-10 m³, a 60 per cent miss). At 1 mm the box's face came within about 4 per cent of the plate
  but its step stopped on a sideways GMRES solve (2.1e-11 against 1e-12), and the shoulder's 1012-part third step
  stopped the same way (6.05e-12); since the final review (number 97) a stall at or below 1e-9 is accepted and counted,
  so neither solve would stop its run now, and neither run was repeated. On the old design, markers 0.4 mm apart took
  1528 s for the box's first step against
  31.4 s and 25 parts against 1 (design 21.1), so 21.10's two prerequisites come first: runoff's cost growing no faster
  than the number of markers, and the part limit untied from the marker spacing.
- **The thin case of design 21.11** (T11-b). The shoulder stops at 1.0 s on a real thin site under the composite law
  (metal 1.39 mm through the cap's marker, under 2H = 2.4 mm), so it cannot be judged at 2.0 s until the thin case is
  designed. Whether the composite law's runoff cutting a trough about 1.2 mm deep in the heated cap in 0.5 s, and
  running metal 0.4 mm thick over the edge, is physical, or a sign that the law's viscosity near the bridge is too low
  for the skins' runoff, is not judged by these runs.
- **The shoulder's fin may be the box's edge-ownership failure** (Task 11's review). The cut-face marker's skin holds no
  non-rigid metal (0.0 mm at 858.9 K) while the fin moved about 0.1 mm and the cap inboard lost about 1.2 mm, which is
  D1's mechanism, so 21.10's finer markers are a candidate remedy beside 21.11. And the box's face is much less even
  under the composite law than under the bridge law: its markers' recessions spread over 0.405 mm against 0.110 mm, the
  largest tilt of a face marker rose from 3.4 to 17.3 degrees, and inner points receded as little as 0.485 mm, 39 per
  cent behind the reference, which the face-centre criterion hides.
- **The node pass still flags a spurious pair at its threshold** (T11-c). On the rounded spraying cube under the
  composite law the stop fires on a pair whose columns share no boundary at any level and whose metal is 3.83 and 8.30
  mm thick through the markers; the thin test's measures read 2.18 and 2.40 mm against 2.4. Until a remedy is chosen
  (for example, a node-pass pair must also share a boundary at some level, or be neighbours among the markers), the
  spraying fixture stays on the bridge law (T9-a) and no unit test sprays under the default law.
- **Markers drift within a step** (D2; T7b-a). Markers stood up to 2.08 mm (1.7 H) off φ within one 50 ms step, and they
  own the re-taken top flux within the step; they are back on φ at the step's end (1.4e-14 m, S15). Whether they should
  settle within a step is open, and with it whether the runoff's forcing can be taken again each part.
- **Concave notches** (T7c-d). A sharp concave notch of gas in thick metal does not stop the run; a reading of 21.11 to
  confirm.
- **Swelling** (T5-n). The volume holds its metal at the solid's density, so φ sits low by the volume's unrealised
  swelling (up to about 15 per cent of the molten depth below β, liquid 2,400 against 2,813 kg/m³) and rises unevenly as
  metal is drawn.
- **Sharp-edge capillarity** (T6-a, T6-d). φ's recovered curvature at a sharp edge pushes a ridge of liquid 1.5 mm high;
  whether to cap or smooth the curvature, or leave it to 21.10 and 21.11, is design section 17 item 8's question. The
  spraying fixture's 2 mm rounded edges are provisional on it.
- **The skyscrapers' regrid ties** (Task 7d's review). The skyscrapers' regrid keeps the plain comparisons; if its fill
  also makes equal fine cells, round-off decides the same tie cells there. Not checked.
- **Also marked for Asha to check in design section 21:** the default depth of 1.2 mm, the superseded hand-down rule,
  the volume's one density and the base as the spraying floor (model_assumptions 18.1, 18.4, 18.8).

Named in the ledger but not listed in section 18.11:

- **"Floor-like" in two senses** (Task 10's review). The spraying module's docstring calls the skin's base design 21.7's
  "floor-like wall at depth H", while the records give a base wall `floor_like_wall` 0 (the rigid front's kind), so a
  skin flowing to its base counts as a rigid-front column in `floor_like_share`. The docstring is to be reworded, and
  design 21.7's line is marked for Asha to check (model_assumptions 18.8 states it as a wording point).
- **The level set's loads and end** (Task 12's report). The level set reads neither main's interim gas-loads provider
  nor the fitted-sphere windward loads, and it has no geometry stop of its own (parity ledger row 10; model_assumptions
  18.12 lists both as not yet verified). Whether it should, and when, is Asha's.
- **The sphere's cost and a whole body on the profile path** (T11-d, Task 11's review). The 50 mm sphere's step cost is
  known only as a lower bound (its first step unfinished after 39 minutes on a loaded laptop), and the profile spraying
  path has not run on a whole body over two or more flowing steps. The open item is the cost sphere through at least two
  steps on the Mac Studio, by runbook, when Asha chooses.
- **The shared slurry functions' default law** (number 95, the final review's M1). `shared/spray.py`'s
  `slurry_properties` and `slurry_modes` default to `law="bridge"` while `rheology.DEFAULT_LAW` is `"composite"`. Both
  architectures pass `law=` explicitly, so no result is affected today, but a future caller that omits it flies the
  bridge law silently. Making `law` required is a small change on main, for Asha to schedule with whoever owns the
  shared spray model.
