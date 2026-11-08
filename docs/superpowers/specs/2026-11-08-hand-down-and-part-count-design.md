# Re-frozen metal kept in the skins, or the part count set by melting: design

**Status: approved by Asha 2026-11-08, to be built with subagents (sub-plan 22 Tasks 5d, 5e and, added the same day, 5f), then flown and
compared. Nothing built yet.**

## Why

The two spreading flights of 2026-11-07 (`--spread-limit 60`, resumed from the snapshot states, rows in
`prototype/work-2026-10-28-skin-spraying-plan/code/reentry_model_output/task5_spreading`) both got past the points where
the earlier skin flights stopped: the 50 mm sphere passed 179.04 s after about 2.1 h of running and reached 179.66 s, its
base fold at most 68.0 degrees against 144 when the earlier flight stopped; the 100 mm sphere passed 33.0 s after about
4.1 h and reached 36.18 s, fold at most 86.5 degrees against 138. The mesh stayed above its radius-ratio floor, mass closed
to 1e-15 kg, and the rebuilds' energy residual was 3.7e-4 and 3.9e-4 of the heat absorbed (4.9e-4 and 4.4e-4 in the
earlier flights).

Both then stalled. The number of parts a 12.5 ms macro step runs as grew by a quarter every step: on the 50 mm sphere
from 9 at the old stop to 9 765, on the 100 mm sphere from 1 to 1 537, the last steps taking 6.7 h and 3.75 h, while the
metal handed down stayed about 0.2-0.3 g a step. Asha stopped them on 2026-11-08 at 16:55.

The cause, inferred from the code and not instrumented: `SkinBody.advance` starts the next macro step with
`ceil(MOVE_HINT_MARGIN * ratio)` parts, where `ratio` is the largest over this step's parts of the part's move ratio
times the number of such parts in the macro step (`dt_macro / dt`). The move ratio includes the hand-down. Since the
hand-down cap of 2026-10-27 (spec 6.6, 7.1) a part hands down at most `HAND_DOWN_SHARE` of the allowance, so the
hand-down can never split a part; but a capped hand-down, which by the cap's own reasoning does not shrink with the
part ("a held cell released whole"), still reads in the hint as every part needing the full allowance, and the next step
gets 25 % more parts. On the 50 mm sphere each of the last eleven counts is exactly the previous times 1.25 rounded up,
which needs a part's move within one part in 10 000 of the allowance, where the capped hand-down sits. The earlier 50 mm
flight was already climbing (134 parts at its stop), so this predates spreading.

Asha then asked why metal is handed down at all, and observed that re-frozen metal need not pile up: it can heat again
and flow on to other skins, or be sprayed. The earlier flights agree: the 100 mm sphere handed down 64.1 g in all but only
about 26 g stood above the starting surface at its stop, the 50 mm sphere 27.7 g and 7.8 g; most of what froze melted
again, today by a round trip through the 3D model (handed down, then absorbed back).

## Decision (Asha, 2026-11-08)

Build both and compare them on flights:

1. **Re-frozen metal stays in the skins**: `--hand-down on|off`, default `on`. With `off` no metal is handed down at
   all: every rigid cell stays in its skin, the base moves only by draws and absorption, and the deposit heats, flows and
   sprays where it is. **Spreading is off with it** (Asha: spreading only works with the hand-down): the command refuses
   `--hand-down off` with a spread limit, given or carried by a resumed state.
2. **The part count set by melting**: `--part-count-from moves|melting`, default `moves` (today's rule). With `melting`
   the next macro step's part count is set by the draws alone: the largest over the step's parts of the part's largest
   draw, `drawn / (rho * area)`, over the allowance, times `dt_macro / dt`, times `MOVE_HINT_MARGIN`, rounded up. The test
   that splits a part is unchanged (it still reads `|drawn - handed|`): a capped hand-down is below the allowance, so only
   a draw can split a part, and the count loses nothing it was there to predict. **Flown with `--spread-limit 60`** (Asha).

3. **Sideways heat between cells that touch** (Asha's idea, 2026-11-08, first written in the side session's draft
   `docs/superpowers/specs/2026-11-08-refrozen-metal-in-skins-design.md` section 3; Asha then asked for the heat to be
   paired correctly on the comparison flights): `--sideways-heat depth|height`, default `depth` (today's pairing by depth
   below the top, spec 6.3). With `height` each cell's range is its height above its column's base, the neighbour's shifted
   by the offset between the two bases along the mean of the two columns' directions, and cells exchange heat by the
   overlap of their ranges times the shared face's conduction weight. Where the tops differ the part of a taller column
   above its neighbour's top is a side wall facing the gas: it exchanges no heat with the neighbour and is not heated or
   cooled from the gas (stated). Tracking every cell's Cartesian position would give the same answer at 2 mm spacing at
   more cost and is not chosen (the side session's call). The flow's sideways transport is unchanged. **On all four
   flights** (controller's reading of Asha's "I want the heat paired correctly"), and once the comparison decides the
   part-count fix, the flight is re-run with it (Asha).

**Defaults** stay today's (controller's call, in the description Asha approved): every earlier run stays reproducible
until the comparison decides, as with `--spread-limit`.

**A stall guard** (controller's call, in the description Asha approved): a macro step that would start with more than
`PARTS_MAX = 1000` parts stops the run with a `RuntimeError` naming the count, before it runs. Ordinary steps have needed
one to ten parts. It applies to every run, so it is the one change without a new option; only a run that would have
stalled is affected.

## What each costs, expected (controller's reading, not tested)

- **Hand-down off.** The base no longer moves outward, and every fold that stopped a flight sat at the toe of the
  re-frozen ring, where the hand-down pushed the base out; the receded nose had none. Costs: sideways conduction pairs
  cells at the same depth below their tops (spec 6.3), so across a toe whose top steps about 1.4 mm in one triangle it pairs
  cells at different heights; the fold moves to the top of the metal (178 and 179 degrees in the earlier flights, not a
  mesh, but read by the flow's push and Girin's spraying); a column is a straight stack on a fixed dual area, so deep
  columns on a curved body misplace volume (the 50 mm sphere's pools were already 22-45 mm deep on a 25 mm radius). Cells
  are cheap: about 90 for a 10 mm column against 59 at the minimum (spec 6.2, estimated).
- **Part count from melting.** Held metal drains at one allowance per part however many parts the step has; on the
  spreading flights it was never deeper than 1.48 mm (50 mm) and 1.12 mm (100 mm), about three parts. If draws also come in
  lumps the count can still grow; the check and the flights will show it.

## The comparison

Four flights on the laptop, on mains, resumed from the snapshot states at 170.5 s (50 mm) and 22.0 s (100 mm), before
the metal flows, to the old end times 190 s and 45 s, a frame every macro step (the top's fold is measured from frames),
states saved at 178.5 s and 32.5 s, slurry spray off as in the spreading flights:
`--hand-down off`, and `--spread-limit 60 --part-count-from melting`, on each sphere, all four with
`--sideways-heat height`.

Compared: how far each gets and its wall time; the fold of the base (`surface_fold_deg`) and of the top (frames); the
deepest skin; the books; the heat absorbed and the mass lost, so that a difference in the physics, not only in the
running, is seen.
