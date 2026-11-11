# Spreading the re-frozen deposit: design draft

**Status: built, not yet flown (2026-11-06). Built with subagents as Asha decided; checked on one 50 mm flight to 177.0 s, see "As built and checked" at the end.** First drafted the same day with nothing built. Asha's idea, in answer to the Task 5 skin flights stopping: "Could we
smooth the re-frozen mass by spreading the re-frozen mass around to neighboring cells to keep the angles reasonable?"
She chose to run snapshot flights and draft this on paper together, with no model code until the snapshots are in. The
snapshots came in on 2026-11-06 and answer the open questions below; they change one rule of the proposal and show a second
problem spreading does not reach. Measurements: `prototype/work-2026-10-28-skin-spraying/measure/task5/snapshots/findings.md`.

## Why

Both skin flights of 2026-11-05 stopped because a surface corner's move could no longer be absorbed, even by a fresh
mesh: the 100 mm sphere at 33.0 s of flight (fold 138 degrees) and the 50 mm sphere at 179.0 s (fold 144 degrees).
Spec §12 puts the limit of the 2 mm surface at about 75-95 degrees.

Most of the metal that enters the skins re-freezes and goes back to the body: on the 100 mm sphere 66.3 g entered and
64.1 g was handed down (97 %), on the 50 mm sphere 39.8 g and 27.7 g (70 %). Those are running totals: the deposit melts
again as the ring it forms moves outward, so the metal standing above the starting surface at the stop is about 26 g and
7.8 g. The fold is flat until the metal flows, then climbs with the deposit: 11 degrees at 2.4 g handed down, 66 at 8 g,
109 at 19 g, 137 at 64 g.

The snapshots place it. The folded edges follow the outer edge, or toe, of a ring of re-frozen metal that travels outward
ahead of the flow at about 4 mm/s on both spheres, and none sits on the receded nose. What stops a flight is the 3D mesh's
radius-ratio floor of 0.05. The lowest ratio fell to 0.098 at the 100 mm stop, in tetrahedra at the same place as the worst
folded edges, while the surface triangles kept their shape, so the fold is the likely cause; the refused move itself is
not recorded.

Spraying alone is unlikely to cure it: the element model sprayed 4.4 g over the same stretch of the 100 mm flight,
against 64 g re-frozen. Measurements: `prototype/work-2026-10-28-skin-spraying/measure/task5/stop_check/findings.md`.

## What the model does today

- A corner's move is one number: `inward = drawn/rho + depth*area - handed/rho` (`SkinBody._move_base`). Drawing and
  absorbing move it in, handing down moves it out.
- `surface_fold` is measured on the **base surface**, the 3D model's corners, not on the top of the metal. The snapshots
  show the top folding at least as much: 178 and 179 degrees at the stops, against 137 and 130 on the base.
- Existing smoothing (`movemesh.tangential_smoothing`) slides corners along the surface. It moves no metal and cannot
  flatten a deposit.
- A hand-down is already limited two ways (spec §6.6): metal more than 5 K warmer than the base node waits in the skin,
  and a part hands down at most the move allowance (`MOVE_MAX`, a quarter of the mean edge).
- A rebuild keeps the surface triangulation, so it cannot undo a fold (spec §12).

## The proposal

One new pass in `SkinBody._move_base`, after `hand_down` and before the corners move: where the base as it would stand
after this move folds past a limit angle, part of the metal a corner is handing down goes to its neighbours' corners
instead. This is the angle-of-repose relaxation used for loose material on terrain, applied to the deposit.

1. Find the edges whose neighbouring base triangles would fold past the limit.
2. Move handed mass, with its enthalpy, from the proud corner to its lower neighbours on the column graph, a few passes
   a part at most, so metal travels at most a ring or two of 2 mm neighbours in one part.
3. Corners move, and each column is reconciled, by the spread hand-down `handed_s` in place of `handed`.

**Books.** The skins have already given up `handed`, so nothing in them changes. The spread keeps the sums exactly:
`sum(handed_s) = sum(handed)` and likewise for the enthalpy. `h_outward` becomes the mean enthalpy each corner receives.

**Rules.**

- Only handed metal moves. Draws and absorption stay where they are.
- A receiving corner respects the move allowance, and a donor keeps what balances its draw. The 5 K fault tolerance does
  not apply at a receiving corner: the lower neighbours of a toe are the cold corners ahead of it, across folded edges the
  3D temperature jumps by a median 62 K on the 100 mm sphere and 31 K on the 50 mm, and about 90 % of the edges exceed 5 K, so
  the rule would refuse most of the spreading it is for. See "Decision: the 5 K rule".
- What cannot be placed stays where it was, and is counted, as `held_mass` counts what the tolerance holds.
- A limit that does not bind changes nothing to the last bit, as with the cap of 2026-10-27.

**Starting limit: 60 degrees** (my proposal), below spec §12's 75-95 and above the 37 degrees the 50 mm base shows before
any melting. Set by a command-line option, `--spread-limit`, off by default until the flights have tested it, so that every
existing run stays as it is.

## What it costs physically

It moves solid metal sideways with no physical mechanism. The justification is resolution: a 2 mm triangle cannot carry
a ridge narrower than itself, so the excess is put where the grid can hold it. The total deposited is unchanged. Its
shape at the 2 mm scale is then set by the limit we choose, and the sensitivity to that limit has to be measured and
reported. If the snapshots show a real rivulet rather than a grid product, this is a stated approximation, not a fix.

It changes the droplet result only through the shape the metal flows over. It does nothing to get metal off the body,
so the question of how the skins spray before the Knudsen number reaches 0.01 stays open and separate.

## Decision: the 5 K rule

**Decided by Asha, 2026-11-06: spread metal lands hot.** Hot metal freezing on cold metal is what happens at a real toe. A
move's enthalpy correction is shared over all the moving nodes (`MovingMesh.move`), so booked through `h_outward` the heat of
metal landing tens of kelvin above its node would go to the whole zone rather than to that corner. So the metal crosses the
base at the receiving node's enthalpy, and the rest of its heat goes to that node through `SkinBody.pending`, the deferred
per-corner heat the 3D solve takes as a load in the next part and which the energy books already count. The landing is not a
jump fault; its largest size is recorded.

## What spreading does not reach

Metal gathers at the tips of lobes that grow on the front, and spreading moves only metal being handed down. On the 50 mm
sphere at 179.0 s four pools 95-105 degrees from the hottest point hold 12.4 g in 101 corners deeper than 5 mm, up to 22.4 mm
on a 25 mm radius, and the deepest rose from 4.5 mm in the last 2 s. The element flight had sprayed 3.7 g by then. The
100 mm sphere shows the same start, 3.7 mm deep at about eleven lobe tips. With a toe held below the limit, these pools are
the likely next stop, and they need the spray outlet or a flow that can leave them.

## The snapshots' answers

1. **Where is the fold?** At the toe of the travelling ring of re-frozen metal, all the way round. On the 100 mm sphere at
   the stop, 347 edges exceed 60 degrees, all between 50.6 and 66.5 degrees from the hottest point. None is on the receded
   nose.
2. **Ridge or sawtooth?** A ridge. A corner's roughness correlates positively with its neighbours' (0.25-0.73), and the top
   follows the base, so spreading handed metal is the right shape of remedy and no hand-down margin is needed.
3. **How wide is the ridge?** About 9 mm at half height on the 100 mm sphere, five triangles, so its body is resolved. Its toe
   is not: at the worst edge the base steps about 1.4 mm across a triangle about 1.5 mm wide. On the 50 mm sphere the ridge
   itself is one to three triangles wide. No flight has yet tested this against a finer grid.
4. **Why did the 50 mm fold rise from 7 to 37 degrees by 19.5 s of flight**, long before any melting? *Partly answered
   from the snapshots, 2026-11-06:* it happens in the first 2-5 s, on 14-15 edges about 90 degrees from the hottest point,
   and then stays. Corners there moved up to 1.22 mm along the surface and 0.47 mm radially by 5 s (median 0.07 mm along
   the surface over the sphere), so it looks like the mesh relaxing under the tangential smoothing at the start, not a
   deposit. Why at that ring is not traced. Measured once, on this flight. Whatever its cause, it starts the 50 mm base
   37 degrees into any limit we choose.

## Alternatives

- **A hand-down limited by its neighbours' depth** (spec §12's option 3): the excess waits in the skin instead of moving
  sideways. Simpler, but the columns then grow: the 50 mm flight's were already 18 mm deep, on a 25 mm radius.
- **Local refinement** (spec §13): the physically honest answer, and the largest piece of work.
- **A margin on the hand-down**: was the cure if question 2 had come back as a sawtooth. It did not, so this is dropped.

## When built

Tests first, through the plan as usual. Mass and enthalpy exact to round-off with spreading on; bit-identical when the
limit does not bind and when off; a one-column ridge flattened below the limit on a test surface; move allowance respected
at the receiving corner; hot landing's heat reaching that corner; a bad limit exits 2. Then spec §6.6 and §12 amended and a fact added.

The check on the real flights: resume from the states the snapshot flights saved at 22.0 s and 170.5 s, which skips the
hours before the metal flows (this needs the resume command of Task 5b). Success is a fold held below the limit past
33.0 s and 179.0 s with the books closed, repeated at a second limit to measure the sensitivity. On the laptop those windows
took about 4.9 h and 2.8 h. A flight that then stops on the pools is a result about the spray outlet, not a failure of
spreading.

**Cost:** about one to two days to build and test (my estimate), plus the resume command, then the flights.

## As built and checked (2026-11-06)

The decisions above stand and were not changed in the build: the pass in `SkinBody._move_base`, 60 degrees as the starting
limit, the share of 0.25 and four passes a part, only handed metal moving, and Asha's decision that spread metal lands hot.
Built in sub-plan 22 Tasks 5b (the resume command) and 5c (the spreading), with `docs/superpowers/plans/2026-11-06-deposit-spreading-build.md`;
recorded in spec 2026-10-11 §6.6, §10.1, §10.2, §12 and §16 and in facts 191 and 192 of the sub-plans' shared context.

- **`--spread-limit` is off by default** until flights have tested it. The draft Asha read left the default unstated ("set by a command-line option, with `off` available so that a flight can compare"); the controller set it off when writing her decision of 2026-11-06 into the design, because no flight has tested the pass and every existing run should stay as it is, and told her so; it stands unless she changes it.
- **The check**, on the 50 mm state saved at 170.5 s, resumed once with a limit of 60 degrees, one mesh: to 172.0 s all 84
  shared columns equal the snapshot flight's at full precision. To 177.0 s the largest fold is 64.86 degrees against 116.75 in
  the snapshot flight, 5 edges past 60 degrees at 177.0 s against 108, 1.704 g spread against 13.63 g handed down, the largest
  landing 57.4 K above its node, mass closed to 4.5e-16, and no fault against 2. The fold stays a little above the limit, as a
  limit on an estimate would.
- **The numbers describe the code just before the review's `Spread.landed` fix** (corners that took hot metal and gave more
  onward sent their heat across the base and were left out of the largest landing). The fix was not flown.
- **After 176.0 s the two flights are not like with like**: the continuum switch latched at 176.0 s with spreading and 176.5 s
  without. To 176.0 s, 11 macro steps each, the lowest 3D radius ratio is 0.0756 against 0.0678 and the deepest skin 2.22 mm
  against 3.04 mm. A radius ratio near 0.05 is where a rebuild fires by design, so on its own it says little about the mesh;
  to 175.5 s it is 5 rebuilds against 8 and a lowest ratio of 0.0945 against 0.0787 (one flight each, identical to 174 s).
- **Untested**: whether a flight gets past the old stops, 33.0 s on the 100 mm sphere and 179.04 s on the 50 mm sphere, since
  the resumed flights have not been flown; the sensitivity to the limit; and the pools at the lobe tips, which spreading does
  not reach.
