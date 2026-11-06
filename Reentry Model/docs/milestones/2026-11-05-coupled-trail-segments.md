# Fully coupled trail segments (back-up, 2026-11-05)

Asha asked on 2026-11-05 to keep this as a back-up. Nothing here is built. The wake column's segmented trail
(`reentry_model/wake_column/segments.py`) solves each segment alone and adds spreading along the trail afterwards
(`along_trail`). If that shortcut fails a check, this is the route that replaces it.

## Why the shortcut is used

Spreading along the trail only matters once the radius nears the segment length: at 70 km, 184 m segments and a
200 um droplet, its effect relative to spreading across the trail is about 1 % at 70 s and 100 % at 1,550 s. Late
chemistry does not depend on concentration: with four times less and four times more aluminium (70 km, one case) the
ratio of Al(OH)2 to AlOH at 12 hours is 2.89 in all three runs, and smoke uptake changes only with the gas left. Plume
particles do scale with concentration, but they form early. So spreading can be applied afterwards in closed form.

## When to pick this up

- A case where late gas chemistry does depend on concentration. Check the shares of gas species at 12 hours for
  changed aluminium, at the altitudes where plume particles are largest.
- A process added later that is second order in the aluminium and acts when the segments overlap, for example slow
  coagulation becoming important, or aluminium consuming enough of an ambient species to matter.
- Spreading that differs along and across the trail, for example from wind shear, if the kernel in `along_trail` can
  no longer be written in closed form.

## What it would take

One stiff system holding every segment, about 1,000 unknowns each, so about 24,000 for 24 segments. Each segment would
gain exchange with its neighbours, a flux between segments i and i+1 of D (c_i - c_i+1) A / l for every aluminium-bearing
species and particle bin, with l the segment length and A the cross-section, and segments would join the system at
their laying times. The Jacobian is block tridiagonal, so the solver needs sparse linear algebra, which `Box.solve`
does not use today. Segments could no longer run on separate processors. The cost has not been measured.

A cheaper variant runs segments alone for about the first minute, when exchange is under 1 %, and couples them after.
It needs the same solver work.

## The check it would give

One coupled run with a few segments would test `along_trail` directly. Asha set this aside with the coupled option
on 2026-11-05, so the shortcut's only checks so far are the closed-form tests and the concentration test above.
