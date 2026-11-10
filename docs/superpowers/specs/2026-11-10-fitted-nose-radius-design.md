# The windward loads from a sphere fitted to the skins' tops: design

**Status: decided by Asha on 2026-11-10; to be built as sub-plan 22 Task 5g.**

## Why

`SkinBody.nose_radius()` returns the mesh's original radius (`mesh.params["radius_m"]`: 50 mm on the 100 mm sphere,
25 mm on the 50 mm sphere) whatever the shape. That radius feeds the windward loads:

- the stagnation-point heating (`CoupledRun.loads_at`, through `heating.evaluate`: Fay-Riddell, Sutton-Graves or DKR, each
  scaling as one over the square root of the radius);
- the surface flow (`SkinBody` at about line 790, through `SurfaceFlow.evaluate`): the edge state along the surface at
  s = R theta, so the pressure gradient that pushes the melt (dp/ds, one over the radius), the boundary-layer thickness
  that sets the shear, the shock Reynolds number and the wall Knudsen number, and through them the state Girin's spraying
  reads.

The angles of the top surface's corners to the flight direction already follow the shape (`_refresh_geometry`); only the
radius is fixed. The element model has fitted its windward cap since 2026-09-21 (`MeltingBody._fit_nose`, with
`--size-feedback current`); the skins never did, and `cli.py` refuses `--size-feedback` with `--surface-model skin`.

Measured from the comparison flights' frames (2026-11-10; the tops of the columns whose normal is within 30 degrees of
the direction of motion, an unweighted least-squares sphere, `body.fit_sphere`): on the 100 mm sphere the cap's radius is
49.7 mm at 22.5 s and about 95-107 mm from about 30 s to the end, on both flights; on the 50 mm sphere 24.7 mm at 171.0 s,
83 mm at 176.0 s, 45 mm at 181.0 s, 34 mm at 184.0 s and 24-38 mm later on an irregular shape. At twice the radius the
stagnation heating is about 30 % lower and the pressure gradient half of what the flights applied.

## Decision (Asha, 2026-11-10)

1. **Fit the gas-molten interface**: the sphere is fitted to the tops of the skin columns (the skyscrapers), not to the
   base (the rigid surface).
2. **The bound is a comparison, not a limit**: the element model's upper bound, `NOSE_CAP_FACTOR` (1.67) times the
   windward transverse radius (a flat face heats like a sphere of 1.67 times its radius, Boison and Curtiss 1959), is
   computed and reported and is not applied. The fitted radius is reported every macro step.
3. **On by default**: size feedback is on by default for skin runs.

## The controller's calls (2026-11-10, for Asha to check)

- **Which tops**: the corners of the top surface (`self.surface`) whose angle to the direction of motion is at most
  `NOSE_CAP_ANGLE` (30 degrees, the element model's), because the stagnation formulas need the curvature at the
  stagnation point; the whole windward half would average in the shoulders. Unweighted, through the corner points, as
  the element model fits through its patch centroids.
- **The windward transverse radius** for the reported bound: the largest distance from the flight axis (through the 3D
  mesh's mean point, along the direction of motion) of the top corners at most 90 degrees from the direction of motion,
  so that a leeward pool standing off the back does not enter it.
- **No fit**: fewer than 4 cap corners, or a fit that is not finite and positive: the original radius is used and the
  fitted radius is reported as not-a-number.
- **A concave cap** (the fitted centre on the flow side of the cap's points, so no sphere of positive curvature describes
  it): the radius used is the reported bound, the flat-face equivalent, and the row records the cap as concave. Not seen
  in the frames; recorded so that it is not silent.
- **When**: the fit is refreshed with the top surface in `_refresh_geometry`, so the loads of each macro step read the
  shape at its start, as the angles do.
- **The option**: `--size-feedback initial|current` with `--surface-model skin`, default `current` on `run`; on `resume`
  the default is the saved run's, and states saved before this task (which have no setting) resume as `initial`, as they
  were flown. `initial` reproduces every earlier output to the last bit, apart from the added columns.
- **Columns** added to the skins' history: `nose_radius_mm` (the radius the loads used), `fitted_nose_radius_mm` (the
  fit, not-a-number where there is none), `nose_bound_mm` (1.67 times the windward transverse radius, reported only),
  `nose_cap_points` (the corners fitted) and `nose_cap_concave` (1 where the cap was concave).
- **Not changed**: the droplet cap's body radius (`transverse_radius()`), and the body Knudsen number's length (the
  element model's size feedback also moves that to the equivalent diameter of the remaining mass; the skins lose mass
  only by spraying).

## Further decisions (Asha, 2026-11-10, later the same day): sub-plan 22 Tasks 5h-5j

The comparison flights' spray records (`<flight>/source_rows`, 186.5-187.3 s on the 50 mm sphere) showed: no release on
Girin's thick branch (the body Knudsen number was about 0.039, the merged regime, where `surface_flow` gives every patch
the Couette closure and `spray.melt_layer` leaves delta_m undefined); Girin and Kopyt's thin branch releasing in that
regime, on layers 37-63 mm deep, reading an edge state from a shock-layer construction the model calls invalid there;
the Rayleigh-Taylor droplets (12-17 mm) set by the cap of a quarter of the body's transverse radius, itself inflated to
57-64 mm by misplaced column tops; and releases from top corners up to 82 mm behind the original centre whose normals
faced forward (median 12-20 degrees off the flight direction), which the angle-based windward test (normal within 90
degrees) counted as windward although they sit behind the body.

1. **Defaults** (Task 5h): for melting skins, `--hand-down off` and `--part-count-from melting` become the defaults of
   `run` (the approach the comparison chose).
2. **The Rayleigh-Taylor cap** (Task 5h): the quarter-of-the-body-radius cap is removed from the Rayleigh-Taylor branch;
   the crest limit (a droplet holds no more liquid than one crest of its wave on the depth available) is kept.
3. **The thin branch's Knudsen gate** (Task 5i): the thin branch sprays only in the continuum regime, where Girin's
   closure applies (today's gate on the thick branch: a distinct shock, body Knudsen number below 0.01, and a wall
   Knudsen number below 0.01). A new flag lets the thin and thick branches spray in the slip regime, Knudsen numbers from
   0.01 to 0.1; a test case is flown with it on.
4. **The leeward test by the centroid** (Task 5j): no Girin spraying behind the plane through the body's new centroid
   (the mass-weighted centre of the rigid and the molten metal, core and skins) at right angles to the flight
   direction, because there is no shear or pressure on the leeward side; a corner is windward only if it lies ahead of
   that plane **and** its normal faces the flow (Asha's yes to the controller's question 2); the Rayleigh-Taylor mode is
   excluded behind the plane too (Asha, question 3). The windward sphere fit takes its points from the windward corners
   so defined.

### The controller's calls on what Asha left open (2026-11-10, for her to check)

- **The fit region** (question 1, unanswered): the fit takes the corners that are windward by item 4 and whose normals are
  within 30 degrees of facing the flow, not the whole windward half: on the 100 mm sphere the whole half ahead of the
  centroid fits 50.9 mm at 32.5 s and 75.5 mm at 43.0 s, where the 30-degree cap fits 104.0 and 100.1 mm, because the
  flattening is confined to the nose, and the stagnation heating and the push on the melt near the nose are set by the
  curvature there.
- **Which Knudsen numbers** (question 5, unanswered): the continuum gate is today's Girin closure (shock-layer branch and
  wall Knudsen number below 0.01); with the slip flag on, Girin's closure, and with it both branches, applies where the
  body Knudsen number and the patch's wall Knudsen number are both below 0.1, the edge state then read from the
  shock-layer construction that the merged branch otherwise treats as invalid (the flag accepts that knowingly). The
  rarefied branch (Girin and Kopyt on the freestream, free-molecular body) is a thin-branch variant outside the continuum
  and is gated with it.
- **The test case** (question 6, unanswered): the 100 mm sphere resumed from its snapshot at 22.0 s to 45 s with all the
  fixes, once with the slip flag off and once on (its body Knudsen number runs from 0.019 to 0.012 over the window, all
  in the slip band).
- **The leeward test reaches every windward consumer**: the stagnation heating's angular law, the surface flow's shear and
  pressure gradient, and the spray plan read one effective angle per corner, the corner's own where it is windward by
  item 4 and past 90 degrees where it is not; the geometry (`theta`) itself is unchanged.
- **Options and reproducibility**: each change is an option of `run` and `resume` in the existing pattern, with the new
  behaviour the default of `run` and `resume` defaulting to the saved run's, so that states saved before resume as they
  were flown: `--leeward-test facing|centroid` (default `centroid`), `--slip-spray off|on` (default `off`), and the
  Rayleigh-Taylor cap and the thin branch's gate behind `--spray-gates 2026-11-10|before` (default `2026-11-10`). The
  element model's spraying is not changed.

### Amended by Asha later on 2026-11-10: the slip band follows the plan of 2026-11-06

Item 3 and the controller's call on which Knudsen numbers are replaced. Asha accepted the calls on questions 1, 4, 5 and
6, then pointed to a note a side session wrote, `docs/superpowers/plans/2026-11-06-slip-band-spraying.md` (never built),
whose six decisions of 2026-11-06 she confirmed: the band is the merged branch with a body Knudsen number from 0.01 to
0.1 (`KN_BODY_SLIP`); Girin's conjugate depth in the band comes from a momentum balance on the shear integrated from the
nose, delta_m^3 = 3 mu^2 I / (rho tau^2 r) with I the integral of r tau ds, not from Eq. (2) (so no gas-state or
boundary-layer construction is read for it; the depth jumps by about 2 to 2.5 times at 0.01, accepted); the flowing
depth may spray; `--slurry-spray on` lets the slurry spray in the band; the shear is the bridged one and the thin branch
keeps its edge-state wavelength (flagged); the element model gets the band too. The option is that plan's
`--slip-band on|off`, default `off`. Asha's additions of 2026-11-10: a patch is in the band only if its own wall
Knudsen number is below 0.1 as well; and the thin branch (with the rarefied branch) sprays only where a conjugate depth
exists, Girin's closure below 0.01 or the band's depth with `--slip-band on`, under `--spray-gates 2026-11-10`, in both
surface models. With Task 5j the band's shear integral reads the windward corners by the centroid test, through the
effective angle the surface flow receives.

### Amended by Asha later on 2026-11-10: the bound applies on the flat side (option B)

Decision 2 is replaced on its upper side. Task 5g's review found that with no bound a nearly flat cap reads discontinuously:
a 30 mm cap sagging by one micrometre fits 450 m read as convex (about 1 % of the original sphere's stagnation heating)
and takes the bound, 83.8 mm, read as concave (about 77 %); and on the 100 mm flights the cap fitted 95-107 mm from about
30 s against a bound of about 84 mm. Asha chose option B: the radius used is the fitted radius when the cap is convex and
no flatter than the flat-face value, and the flat-face value (`NOSE_CAP_FACTOR`, 1.67, times the windward transverse
radius; Boison and Curtiss 1959) when the convex fit is flatter, when the cap is concave, and when it is flat to rounding.
There is no lower bound. The fitted radius and the bound are still reported every step. The windward transverse radius is
taken over the columns whose base normal faces the flow (review finding 1).
