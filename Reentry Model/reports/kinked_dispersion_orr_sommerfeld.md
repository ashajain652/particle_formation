# The viscous check of the floor test: what an Orr–Sommerfeld solver says about FLOOR_RE and KINK_DROP

Task 10 of the kinked-profile dispersion plan (the spec's stage 6), 2026-11-12, branch `kinked-dispersion`. Analysis
only: no model code was changed, and neither constant was changed.

## 1. What you asked for, and whether it is done

You asked for an offline viscous check of the two constants that the new profile path uses to decide where the melt's
gradient instability (Girin's mode, driven by the kink in the speed profile) can grow; this was your decision 4 of
2026-11-11. The first constant is `FLOOR_RE`, the threshold of the floor Reynolds number below which the slurry under a
kink is treated as a rigid floor; it is 1.0 today. The floor Reynolds number is
Re_f = |u_top − u_kink| d / (1.225 ν_below): the speed jump between the surface and the kink, times the kink's depth d,
over the slurry's kinematic viscosity just below the kink times Girin's fastest wavenumber-depth product 1.225. It
compares how fast the wave carries momentum along with how fast the slurry's viscosity spreads it. The second constant is
`KINK_DROP`, the largest ratio of the shear rate below a breakpoint to the shear rate above it (the slope ratio) for that
breakpoint to count as a kink, a candidate floor, at all; it is 0.2 today, a five-fold drop. Both were labelled "assumed
until stage 6". The check was to be an Orr–Sommerfeld solver (the linear stability equation of a viscous shear flow) for
the whole melt column, validated against three published results, run on sampled flight columns and on synthetic columns
across a floor Reynolds number of 1, and used to recommend values, with the decision left to you.

It is done, with four concerns, listed here so that they are not buried. First, Hooper and Boyd's formula as our
research notes transcribe it has its viscosity ratio upside down; the solver matches the paper's own definition to
within 0.02 to 4.7 percent, and the notes' version is wrong by a factor of the viscosity ratio cubed (section 3.2).
Second, the two "earlier" saved states the brief named (d050 at 170.5 s and d100 at 22.0 s) have no flowing top cells,
so the flight sample covers two flight times, d050 at 178.5 s and d100 at 32.5 s, each from three flights saved at that
time (section 3.4). Third, the flight columns with their speeds scaled up were solved at four wavenumbers each, not by
a full search, because one solve of such a column took up to three minutes on the then heavily loaded machine; in eight
of eighteen the fastest of the four lay at an end of the set, and a later search beyond the end changed no verdict and
moved the mean errors by at most 0.013 (section 3.4). Fifteen of the eighteen are now checked at doubled resolution, and
three are not (section 3.5). Fourth, Albert and
Charru (2000), which the brief asked me to read first if it could be obtained, is not available as an open copy; one web
search found only bibliographic records (European Journal of Mechanics B 19, 229–252,
doi 10.1016/S0997-7546(00)00109-6), so nothing in the validations was changed by it.

## 2. What was found, in plain language

**The slurry does not switch from fluid to floor at one value.** As the slurry under a kink becomes more viscous, the
growth of the gradient instability does not drop to zero at a floor Reynolds number of 1. It falls smoothly. It is half
of the inviscid value (the value the broken-line solver gives when the slurry moves freely with the wave) at a floor
Reynolds number of about 3.5, a quarter at about 1, and a tenth at about 0.3. In the idealized two-layer columns this
curve is the same whether the top layer's own Reynolds number is 30 or 10,000, whether the slurry is 3 or 15 times as
deep as the top layer, and whether the Weber number is 30 or 1,000. That collapse is the main result: the floor Reynolds
number is the right variable, as the spec's scaling argument said, but a threshold of 1 sits on the floor side of the
transition, not at its middle, and no single threshold can represent a transition spread over two decades.

**In the flight columns the same curve holds, and the profile path's viscous gate, not the floor test, decides.** With
their speeds scaled until the gradient mode exists, the viscous growth of eighteen sampled flight columns follows the
same curve to within about 0.14 of the inviscid growth on average. All eighteen are fluid-like under the present
`FLOOR_RE`, yet the profile path predicts no gradient-mode growth in twelve of them, because its viscous gate (growth must
exceed twice the layer's effective kinematic viscosity times the wavenumber squared) uses the effective viscosity of the
whole layer down to the rigid front, which the slurry dominates. In those twelve the viscous solver keeps between 5 and
74 percent of the inviscid growth, e-folding times of microseconds to milliseconds.

**Raising `FLOOR_RE` alone would not help the flight columns.** At 3.5 the wall moves up to the bridge's base in about
half the scaled d050 columns, but the bridge above it is itself curved, so its inviscid growth is still about 95 percent
of the whole profile's, and four columns would then spray at nearly the full inviscid growth where the viscous solver
says a quarter to two fifths. Replacing the switch (and the gate, for this mode) by a smooth factor, the floor Reynolds
number over itself plus 3.5, matches the viscous solver far better than any threshold, on both the idealized and the
flight columns.

**`KINK_DROP` changes no verdict under the present gate at the present `FLOOR_RE`, but it decides which kink sets the
floor Reynolds number, and 0.2 is the best of the values tried.** Where the shear stress is continuous across a kink, the
slope ratio is the inverse of the viscosity ratio. With a nearly inviscid top layer (Re_1, the top layer's Reynolds
number, of 30 or more) every floor-like kink in the idealized columns has a slope ratio below 0.04, so `KINK_DROP` cannot
matter there; with a viscous top layer (Re_1 of 3 and 10) floor-like interfaces have slope ratios from 0.05 to 0.22, and
the verdicts stay the same only because the present viscous gate removes every one of those columns. In the flight
columns `KINK_DROP` moves the walls a great deal (0.5 would move 46 of the 69 sampled walls closer to the surface), and it
changes which breakpoint counts as the first kink, whose floor Reynolds number the smooth factor below would use: at 0.5
that kink changes in 58 of the 69 columns at their own speeds and in 14 of the 18 scaled ones, and the smooth factor's
agreement with the viscous solver falls from 18 to 13 of 18. So the evidence supports keeping 0.2.

**At their own speeds the flight columns do not grow by Girin's mode, as the profile path says, but the viscous solver
finds a slower instability of another kind.** Surface tension holds the gradient mode back at the flights' speeds
(the Weber number of the layer above the kink has a median of 0.0095, against a threshold of 3), and the profile path
sprays from no sampled column. The viscous solver nevertheless finds a growing wave in all 69 columns: a short wave that
lives inside the bridge's viscosity gradient, growing at a median of 6 to 15 per second in the d050 states and about 200
per second in the d100 states. It does not come from the cell-to-cell viscosity steps, but it sits on the skins'
cell-scale speed profile, whose kinks at the cell centres are themselves unstable without viscosity, and whether it can
deform the surface enough to spray is not known (section 3.4).

## 3. The evidence

### 3.1 The solver

The solver is `analysis/kinked_dispersion/orr_sommerfeld.py`; its entry point is
`growth_rate(k, y, U, mu, rho, sigma, top="free")`, which returns the complex wave speed c (metres per second) of the
most unstable wave of wavenumber k (per metre), whose growth rate is k times the imaginary part of c (per second). The
column is the skins' own: its flowing cells from the top, each with its own viscosity and density, the bottom face of
the last flowing cell being a no-slip wall, as in the skins' momentum step. The base speed inside each cell is straight
in each half-cell. The speed at a face between two cells is found as the skins find it, by requiring the shear stress to
be the same on both sides of the face; the surface speed is extrapolated from the top two cell centres, as the
straight-segment fit does. The kink this leaves at each cell's centre is where the skins' finite volumes lump the push
(the pressure gradient and the deceleration) that bends the profile.

At a face where the viscosity or density changes, the disturbance must keep velocity and both stresses continuous across
the displaced face. At a kink inside uniform fluid the solver uses the Orr–Sommerfeld equation itself. At the free
surface it uses zero tangential stress, the kinematic condition and the normal-stress condition with a surface tension
of 0.80 newtons per metre; with no viscosity these reduce exactly to the broken-line solver's conditions. Gravity, the
deceleration and the gas's own disturbance are left out, as in Girin's theory and the broken-line solver, and the
disturbance sees each cell's viscosity at its base shear rate, as if the melt were Newtonian. Each cell is discretized
with Chebyshev polynomials, with the stream function and its vorticity both as unknowns, so that no fourth derivative is
formed: an earlier version that formed it lost up to 0.4 percent of Girin's growth to round-off at a Reynolds number of
100,000, which the second-order form removed (a development observation, not recorded in any committed script). The eigenvalues are found by a shift-and-invert transform about a point
above the unstable region, which keeps them accurate although the problem's matrices are badly conditioned (inverting
the matrix directly lost 1 percent of a growth and reversed one validation's sign, again a development observation that
was not recorded); on long flight columns the same
transform is solved by ARPACK, an iterative method, with a check that it reached far enough from the shift. Spurious
eigenvalues are excluded twice: any growth faster than twice the column's speed range is discarded, and the chosen
eigenvalue's eigenfunction must be resolved, which means that its last two Chebyshev coefficients in every cell are below
a thousandth of its largest.

The surface condition contains one modelling choice that the papers do not settle. When the surface tilts, the gas's
base shear stress can be taken to act as a fluid's stress (then the tilt adds nothing, as at every internal face; this is
the default) or as a pull along the tangent only (then the normal stress gains a term of relative size one over the
Reynolds number). Switching between them changes the gradient mode's growth by less than 1 percent (0.156 against 0.157
velocity jumps per layer depth on the idealized column tried) and the flight columns' own-speed modes by less than one
part in a million, and changes one slow mode of a uniform film with no kink by a factor of five (section 3.3), so the
choice does not affect the recommendations. These numbers, and every later one marked so, are recorded in `checks.json`
by `floor_sweep checks`.

### 3.2 Validations

Every validation in the table below is in `tests/test_kinked_dispersion_orr_sommerfeld.py`, and every value was the
same at doubled resolution to at least six digits. The tests assert the tolerances (5 percent, 0.2 percent for the
density-jump rows and one part in a million for Lamb's wave); the values printed in the table come from one-off runs of
the tests' own functions on 2026-11-12, which are recorded nowhere else but can be reproduced from them. The table compares the solver with each reference on the same
quantity; it shows that the solver reproduces each one well inside the 5 percent tolerance, the largest difference being
one of Hooper and Boyd's cases, where their asymptotic formula itself is least accurate.

| Check | Reference value | Solver | Difference |
|---|---|---|---|
| Damped capillary wave on deep liquid at rest (Lamb's exact relation), wave speed | 0.957637 − 0.153018 i | the same to 13 digits | 1.2 × 10⁻¹³ relative |
| Girin's profile, Reynolds number 10⁴, growth at the broken line's fastest wavenumber (surface speeds per layer depth) | 0.24690 | 0.24685 | −0.020 % |
| Girin's profile, Reynolds number 10⁵, the same | 0.24690 | 0.24689 | −0.0026 % |
| Girin's profile, fastest wave over all wavenumbers, Reynolds number 10⁴: growth, and its wavenumber times the layer depth | 0.24690 at 1.2259 | 0.24685 at 1.2263 | −0.020 % and +0.036 % |
| Girin's profile, fastest wave over all wavenumbers, Reynolds number 10⁵ | 0.24690 at 1.2259 | 0.24689 at 1.2259 | −0.0026 % and 0 |
| Hooper and Boyd, short-wave parameter 5, viscosity ratio 0.5 (growth in units of the lower fluid's shear rate) | 1.6667 × 10⁻³ | 1.6610 × 10⁻³ | −0.34 % |
| Hooper and Boyd, parameter 5, ratio 2 | 1.3333 × 10⁻² | 1.2708 × 10⁻² | −4.69 % |
| Hooper and Boyd, parameter 10, ratio 0.5 | 4.1667 × 10⁻⁴ | 4.1658 × 10⁻⁴ | −0.022 % |
| Hooper and Boyd, parameter 10, ratio 2 | 3.3333 × 10⁻³ | 3.3219 × 10⁻³ | −0.34 % |
| Charru and Hinch, thin lower layer the more viscous: sign of the long-wave growth, and Yih's wave speed | unstable; 0.1817 | unstable; 0.1840 | +1.25 % |
| Thin lower layer the less viscous | stable; −0.1606 | stable; −0.1604 | −0.14 % |
| Thin upper layer the more viscous | unstable; −0.02271 | unstable; −0.02271 | below 0.01 % |
| Thin upper layer the less viscous | stable; 0.08032 | stable; 0.08031 | −0.01 % |
| Girin's profile with the liquid below the kink 1.5 times as dense, Reynolds number 10⁵ (broken line with its per-segment densities) | 0.21678 | 0.21667 | −0.050 % |
| The same with the liquid below 0.7 times as dense | 0.26812 | 0.26806 | −0.021 % |

Some of these numbers need their context. On Girin's profile (a straight shear layer over liquid at rest, here with a
wall twelve layer depths down and a Weber number of 10,000, compared with the broken-line solver on the same geometry)
the viscous growth differs from the inviscid one by two parts in ten thousand at a Reynolds number (surface speed times
layer depth over kinematic viscosity) of ten thousand, and by about eight times less at a hundred thousand, close to the
tenfold fall that a correction proportional to one over the Reynolds number would give; the brief's criterion, "within 5 percent, with a smaller error at
10⁵", is met with a wide margin. For Hooper and Boyd the reference is their Eq. 15: the growth rate equals the lower
fluid's shear rate times m(1 − m)²/[2(1 + m)] divided by the square of the short-wave parameter (the wavenumber times the
square root of the lower fluid's kinematic viscosity over its shear rate), with m the lower fluid's viscosity over the
upper fluid's. The research notes (`fulltext_viscous_papers_read.md`, lines 34 and 37, which the brief quotes) define m
the other way round, as the upper fluid's over the lower's; the paper's tangential-stress condition (its Eq. 5c) and its
sentence "suppose m < 1, i.e. the lower fluid is the less viscous" both fix it as I use it. With the notes' definition the
formula is wrong by a factor of m cubed, 8 for m of 0.5. The remaining errors fall as the fourth power of the short-wave
parameter, which is the truncation of their asymptotic series, not a solver error: the case at parameter 5 and ratio 2 is
the mirror image of a case at parameter 2.5 and ratio 0.5. Charru and Hinch's rule, "stable if the thinner layer is the
less viscous, and unstable otherwise", holds in all four arrangements, and the wave speed matches Yih's leading-order
formula (their Appendix C), measured relative to the interface's speed in units of the lower layer's shear rate times its
depth. The last two rows test the interface conditions used where neighbouring densities differ by more than 1 percent:
the density changes the inviscid growth by −12 and +9 percent, so the agreement is not trivial, and at a Reynolds number
of 10⁴ the gaps are 0.29 and 0.13 percent, so they shrink with viscosity as they should. A further test checks that the
dense shift-and-invert solve equals scipy's QZ algorithm on a Charru and Hinch column, and that ARPACK equals the dense
solve on a d050 fixture column of 50 cells, large enough for ARPACK to run (the test asserts which solve ran), each to a
billionth of the speed scale. In the first version of this report that ARPACK comparison silently ran the dense solve,
because ARPACK is tried only above 360 rows linear in the wave speed; every flight column is above that. The reviewer of
this task also compared ARPACK with the dense solve on two full flight columns, which agreed to 4 × 10⁻¹¹ of the top
speed.

### 3.3 Synthetic two-layer columns

The synthetic columns isolate the floor question. Each has a top layer of depth d over a slurry down to a wall at
depth b, under a free surface, with the shear stress continuous at the kink (a two-layer Couette column). In the top
layer's units (its depth, its velocity jump and its density all 1) a column is fixed by four numbers: the top layer's
Reynolds number Re_1 (its velocity jump times its depth over its kinematic viscosity), the slurry's viscosity over the
top layer's, M, the depth ratio b over d, and the Weber number (density times velocity jump squared times depth over
surface tension). The slope ratio at the kink is then 1 over M, and the floor Reynolds number is Re_1 over 1.225 M; a test
pins that this is what the profile path's own `dispersion.floor_reynolds` computes on the same column. The main grid has
b = 15 d, as in the d050 thick columns (830 micrometres of flowing depth under a 55-micrometre bridge), a Weber number of
1,000, and Re_1 of 3, 10, 30, 100, 1,000 and 10,000, each with ten values of M spanning three decades so that the floor
Reynolds number runs from about 80 down to about 0.08; a shallow slurry (b = 3 d) and a Weber number of 30 were run at
Re_1 = 1,000 as sensitivities. All 80 cases took 222 seconds on four processes on the then heavily loaded machine;
regenerating them with the committed code on a quiet machine took 50 seconds and reproduced every value exactly.

![Viscous over inviscid growth against the floor Reynolds number and the slope ratio](kinked_dispersion_orr_sommerfeld/synthetic_growth_ratio.png)

The left panel shows the main result: plotted against the floor Reynolds number, the viscous growth over the inviscid
growth falls along one curve for every top layer with Re_1 of 30 or more and for both sensitivities; the right panel
shows that against the slope ratio the same data do not collapse at all. The one-parameter curve Re_f / (Re_f + c),
fitted over floor Reynolds numbers from 0.3 to 80, gives c = 3.5, with a root-mean-square difference of 0.045 and a
largest difference of 0.10 over 28 columns. The table below gives, for each series, the floor Reynolds number at which
the viscous growth falls to half the inviscid, and the share kept at two lower values, each interpolated between the
grid's columns linearly in the logarithms of both quantities (interpolating the share itself linearly in the logarithm
of Re_f would give, for example, 0.14 instead of 0.13 for Re_1 = 30 at 0.3); it shows that the half-way point lies
between 2.9 and 5.7 in every series with a nearly inviscid top layer, and that a fifth to a third of the growth survives
at the present threshold of 1.

| Series | Re_f where half the growth is left | Share left at Re_f = 1 | Share left at Re_f = 0.3 |
|---|---|---|---|
| Re_1 = 30 | 2.9 | 0.32 | 0.13 |
| Re_1 = 100 | 3.8 | 0.25 | 0.10 |
| Re_1 = 1,000 | 4.6 | 0.23 | 0.16 |
| Re_1 = 10,000 | 3.3 | 0.25 | 0.13 |
| Re_1 = 1,000, slurry 3 d deep | 5.7 | 0.21 | 0.15 |
| Re_1 = 1,000, Weber number 30 | 4.9 | 0.20 | 0.14 |

Two features of the curve need explaining. Below a floor Reynolds number of about 0.3, the series with a moderately
viscous top layer (Re_1 of 3 to 100) keep falling in proportion to the floor Reynolds number, at about half of it (for
example to 0.013 of the inviscid growth at a floor Reynolds number of 0.025), which is what one expects if the wave
survives only through the slurry's small motion. The series with Re_1 of 1,000 and 10,000 instead level off at 0.13 and
0.09 of the inviscid growth: there the slurry is effectively rigid, and what remains is a slow viscous mode of the
straight top film itself, which the inviscid theory calls stable (Ghosh's result for a straight film on a wall); it grows
at 0.021 to 0.035 velocity jumps per layer depth, a tenth of Girin's mode, and the surface-condition choice of
section 3.1 changes it by about 7 percent (`checks.json`). A third, separate mode appears in a uniform film with no kink at all, at 0.013
times the shear rate; it depends on the surface-condition choice by a factor of five, so I do not rely on it.

The binary floor test can be scored against the viscous solver on these columns, calling the gradient mode "alive" where
it keeps at least half its inviscid growth. On the 57 main-grid columns with a positive inviscid growth, and with the gate
as the profile path applies it, the test with `FLOOR_RE` = 1 agrees with the viscous verdict on 52, with 3 on 54, with
10 on 51 and with 0.3 on 48; its mean error in growth, as a share of the inviscid growth, is 0.21, 0.20, 0.21 and 0.25.
No threshold does much better than 0.2, because the transition it replaces is spread over two decades. The smooth factor
Re_f / (Re_f + 3.5) applied to the inviscid growth, with no gate, has a mean error of 0.05 on the same columns.

`KINK_DROP`, from 0.05 to 1, changes none of these numbers at `FLOOR_RE` of 0.1 to 3, but it does elsewhere. At
`FLOOR_RE` = 10 the test agrees on 52 columns with `KINK_DROP` 0.05 against 51 with 0.1 to 1 (mean errors 0.207 and
0.212), and at 30 on 50, 49 and 47. Without the gate it matters at every threshold: at `FLOOR_RE` = 1 the test agrees on
43, 45, 45, 46 and 46 columns for `KINK_DROP` of 0.05, 0.1, 0.2, 0.5 and 1. The reason is in the columns with a viscous
top layer: with Re_1 of 3 and 10 there are floor-like interfaces with slope ratios of 0.05 to 0.22 (for example 0.215 at
Re_1 = 3 and M = 4.6, where Re_f is 0.53), which a small `KINK_DROP` does not test as kinks. With the gate in place this is
hidden, because the gate removes all of those columns: among the columns with a positive inviscid growth it removes the
mode in nine of nine with Re_1 of 3 and of 10, eight of nine with 30 and six of ten with 100, although in the columns it
removes the viscous solver keeps up to 18, 40, 60 and 49 percent of the inviscid growth in those four series.

### 3.4 Flight columns

The 69 flight columns come from six saved states. From d050 at 178.5 s: the Task 4 fixture's 36 columns (24 of the
spec's thick set, 6 plug columns and 6 others, from the hand-down-off flight) and 8 and 7 more, sampled the same way, from
the parts-melting comparison flight and the spread-60 flight saved at the same time. From d100 at 32.5 s: 6 columns from
each of the three flights saved then (hand-down off, parts-melting and spread-60), all "other" columns, because at that
time no column falls in the thick set or is a plug.
The snapshot states the brief named, d050 at 170.5 s and d100 at 22.0 s, have no flowing top cells. The sampled columns
are committed (`flight_columns.npz`, 0.27 MB), so the sweep needs no saved state. Every column was run at its own speeds,
and every fourth (18 columns) also with all its speeds scaled up until the Weber number of the layer above its first kink
reached 30, a factor of 14 to 930, as a what-if for the faster shear that hard deceleration would bring; the cells'
viscosities were held fixed although the slurry would thin under faster shear, which would raise the floor Reynolds
number further.

**At their own speeds.** Sixty-two of the 69 columns have a kink; their floor Reynolds numbers have a median of 0.050
(the middle 80 percent between 0.024 and 1.15) and their slope ratios a median of 0.045, so 55 of the 69 are floor-like
under the defaults. The Weber number of the layer above the kink has a median of 0.0095 and a largest value of 1.9,
below the threshold of 3 at which the inviscid gradient mode first exists, so the profile path sprays from none of the
columns (the inviscid growth on the whole profile is positive in 7 columns, and the cut and the gate remove it in all
seven). A useful number for the future is the floor Reynolds number each column would have when its speeds have risen
enough to reach the Weber threshold of 3: it has a median of 0.70 to 2.5 in the d050 states and 7.4 to 8.2 in the d100
states. By the curve of section 3.3, a d050 column would keep a sixth to two fifths of the inviscid growth when the mode
first appears, so the floor matters there; a d100 column would keep about seven tenths, so the floor hardly matters.

The viscous solver nevertheless finds a growing wave in every one of the 69 columns at their own speeds. The table
below compares the growth by state; it shows that the d100 columns grow 13 to 38 times faster than the d050 ones, at growth
rates that would e-fold within about five milliseconds.

| State | Columns | Viscous growth, median (per second) | Middle half (per second) | Wavenumber times kink depth, median | Wave speed over surface speed, median |
|---|---|---|---|---|---|
| d050 fixture, 178.5 s | 36 | 8.6 | 2.7 to 21 | 8.8 | 0.62 |
| d050 parts-melting, 178.5 s | 8 | 15 | 6.8 to 95 | 4.3 | 0.63 |
| d050 spread-60, 178.5 s | 7 | 5.7 | 0.37 to 11 | 6.7 | 0.48 |
| d100 hand-down off, 32.5 s | 6 | 200 | 155 to 333 | 7.2 | 0.29 |
| d100 parts-melting, 32.5 s | 6 | 218 | 180 to 263 | 7.1 | 0.32 |
| d100 spread-60, 32.5 s | 6 | 204 | 143 to 268 | 7.4 | 0.29 |

These are short waves, with a wavelength of about 0.9 of the kink's depth (wavenumber times depth about 7), travelling at
a speed found inside the bridge, not at the surface. Five checks, all recorded in `checks.json`, say what they are.
They do not come from the stepped viscosity: when every cell is split into two and then four with the viscosity
interpolated between neighbours, so that the steps are four times smaller, the growth of the three fastest d100 columns
falls by about a fifth and converges (from 688 to 565 to 531, from 664 to 579 to 557, and from 486 to 423 to 409 per
second). In that test the base speed was kept exactly, as the same straight pieces through the original faces and cell
centres, so the kinks in the speed at the original cell centres stayed where they were; the test isolates the viscosity
steps, not those kinks. They do not depend on the surface condition (the same growth with either closure to one part in
a million) or on surface tension (doubling it changes 688 per second to 688). The fastest of all, in the d050 fixture
(889 per second), is different: a long wave at the band's end, faster than the surface, which doubling the surface
tension slows to 677 per second. The fifth check, measured once on five columns (first by the reviewer, then recorded),
runs the inviscid broken-line solver on each column's own cell-scale profile, with a breakpoint at every face and at
every kinked cell centre: those profiles are unstable without viscosity, at a quarter to half of the viscous fastest
wavenumber (1,710, 3,425, 215, 1,401 and 5.8 per second, against viscous growths of 688, 664, 2.6, 6.4 and 0.37), but
neutral at the viscous fastest wavenumber in four of the five (the fifth gives 1.3 per second there). That supports
reading the fastest viscous modes as viscous in origin; it also shows that the cell-scale profile, unlike the fitted
one, is unstable in Rayleigh's sense, so the kinks the skins leave at their cell centres are not innocent. So within
the solver's assumptions, and on the skins' cell-scale speed profile, these are genuine viscous instabilities of the
bridge's viscosity gradient, of the family Hooper and Boyd found at a single viscosity jump. They are outside Girin's
theory and outside the profile path, which is inviscid. What is not known is whether they can spray: their wave travels
at a depth inside the bridge, and the solver gives growth, not the surface amplitude at which crests would break.

**At speeds scaled to a kink Weber number of 30.**

![Viscous growth of the sampled flight columns](kinked_dispersion_orr_sommerfeld/flight_growth.png)

The left panel shows the own-speed growth just described against the floor Reynolds number (the seven columns without a
kink are not drawn). The right panel is the test of the floor law on real profiles: each blue point is a scaled flight
column's viscous growth over its inviscid growth on the whole profile, over the synthetic columns (grey) and the fitted
curve. The flight points follow the curve: the mean difference is 0.14 of the inviscid growth, and on the verdict ("keeps
at least half its growth") the curve agrees with the viscous solver in all 18 columns. The points well below the curve near
a floor Reynolds number of 2 are mostly the d050 columns whose fastest wave lay at an end of the four wavenumbers tried.
Eight of the 18 did, seven at the top of the set (1.41 times the inviscid fastest wavenumber) and one at the bottom. A
search beyond the end (`checks.json`, at 2 and 2.83 times, or 0.35 and 0.25 times) raised the largest growth found by 3
to 78 percent in five of them, and in two of those the growth was still rising at 2.83 times; in the other three it
fell beyond the end. No verdict changes with the higher growths, and the mean errors in the table below move by at most
0.013 (today's model from 0.26 to 0.27, the smooth factor from 0.14 to 0.13).

All 18 scaled columns are fluid-like under the defaults (floor Reynolds numbers from 1.6 to 394, median 3.2). In the
six whose layer above the wall is not dominated by slurry (the five d100 columns and one d050 parts-melting column) the
gate passes the mode and the profile path's growth is the inviscid one, against a viscous 75 to 94 percent of it. In
the other twelve, all d050 columns, the gate's damping is 2 to 29 times the growth, so the profile path predicts no
gradient-mode growth, while the viscous solver keeps 5 to 74 percent of the inviscid growth, between 113 and 340,000 per
second; five of the twelve keep less than a quarter. A gate computed over the layer above the first kink instead, which
I tried as the obvious repair, is worse: it passes the mode where the viscous growth is small and agrees with the viscous
verdict in 12 of 18 columns, with a mean error of 0.34 (`checks.json`).

The candidate constants were scored on these 18 columns in the same way as on the synthetic ones. The table below
compares the profile path as it is with two alternatives; it shows that raising `FLOOR_RE` to 3.5 makes the flight
verdicts worse under the present gate, that `KINK_DROP` of 0.5 changes nothing at the present `FLOOR_RE` with the gate,
and that the smooth factor is the best of all, but only with the first kink chosen at `KINK_DROP` = 0.2.

| Variant (with the present gate unless stated) | Verdicts agreeing with the viscous solver | Mean error over the inviscid growth | Floor-like columns |
|---|---|---|---|
| FLOOR_RE 1.0, KINK_DROP 0.2 (today) | 16 of 18 | 0.26 | 0 |
| FLOOR_RE 1.0, KINK_DROP 0.5 | 16 of 18 | 0.26 | 2 |
| FLOOR_RE 3.5, KINK_DROP 0.2 | 12 of 18 | 0.33 | 10 |
| FLOOR_RE 3.5, KINK_DROP 0.5, without the gate | 13 of 18 | 0.32 | 10 |
| Smooth factor Re_f / (Re_f + 3.5) on the whole profile's growth, no gate, no wall, first kink at KINK_DROP 0.2 | 18 of 18 | 0.14 | none needed |
| The same, first kink at KINK_DROP 0.1 | 14 of 18 | 0.28 | none needed |
| The same, first kink at KINK_DROP 0.5 | 13 of 18 | 0.23 | none needed |

Why does a higher threshold do worse here? At 3.5 ten scaled d050 columns become floor-like. In eight of them the wall
moves from the rigid front, 750 to 7,400 micrometres down, to the bridge's base, 15 to 185 micrometres down; in the other
two it moves to 400 and 835 micrometres. The bridge above that wall is itself curved, so
the inviscid growth of the profile above it is still 94 to 97 percent of the whole profile's in five of them; and the gate,
now computed over the thin bridge alone, lets four of them through at that growth, where the viscous solver says 25 to 42
percent. The floor test was designed for a straight runny layer over a kink; in the flight columns the runny layer is
curved, and the curvature alone carries most of the instability.

**Density jumps at faces.** Each cell keeps its own density in the equation, but at a face the solver ignores a jump in
density below 1 percent (without gravity a density jump enters only through inertia, and very small jumps made
eigenvalues near infinity appear). Across the 69 columns only 8 faces exceed 1 percent, the largest by 1.8 percent, while
a column's density varies by about 7 percent from top to bottom in steps of 0.1 to 0.3 percent. Keeping every jump
changes the growth of five flight cases (three at their own speeds and two scaled) by 0.007 to 0.04 percent
(`checks.json`), so this simplification costs well under a tenth of a percent.

### 3.5 Resolution and run time

Resolution was checked by doubling every region's Chebyshev degree and solving again at the same wavenumber. Every
validation gave the same value to at least six digits. In the sweeps, four synthetic columns (one each from Re_1 of 3,
30, 100 and 10,000) changed their growth by at most three parts in a hundred billion, and seven flight columns (six at
their own speeds, spread through the list, and the smallest scaled column) by at most six parts in a million, the
largest change being in the slowest-growing column checked (0.37 per second); the wave speed moved by at most
1.4 × 10⁻⁸ of the surface speed. In the first version of this report those seven were the only flight columns checked,
which left 17 of the 18 scaled columns, the ones that carry the comparisons of section 3.4, unchecked. `floor_sweep
checks` has since re-solved 15 of the 18 scaled columns at doubled resolution with ARPACK, including the four whose
verdicts decide the comparison at `FLOOR_RE` 3.5 (d050 fixture columns 12, 20 and 32, and spread-60 column 56): their
growth moved by at most 1.2 × 10⁻⁸ of itself and their wave speed by at most 2.3 × 10⁻¹⁰ of the top speed. The other
three (fixture columns 4 and 24 and parts-melting column 44) remain unchecked: at doubled resolution ARPACK did not reach
far enough from its shift, and the dense solve that would then take over is the kind that ran for over 25 minutes on the
first attempt, so the check declines it. Those three still rest on the eigenfunction-tail test alone. A unit-tier test
also re-solves a fixture column at doubled resolution.

The extraction of the 69 flight columns, loading each saved state once, took 41 seconds (timed by the shell, not
recorded). The synthetic sweep's 80 columns took 222 seconds on four processes, and 50 seconds when regenerated on a quiet
machine (`timing_synthetic.json`). The flight sweep's 87 cases (69 at their own speeds and 18 scaled) took 27 minutes on
four processes, 6,452 process-seconds in all, on a machine heavily shared with other jobs; one own-speed column took from
half a minute to six minutes, and one scaled solve up to three minutes, against one to three seconds on the quiet
machine later. Its timing file, `timing_flight.json`, was written by hand from the run's log and the files' times, as its
`provenance` key says, because the run was stopped during its resolution step. The resolution re-solves took 1 and 15
seconds (`timing_synthetic.json`, `timing_resolution_flight.json`), and all of `floor_sweep checks` 54 seconds.

## 4. What it means for the research

The floor test's physics is right but its form is not. The floor Reynolds number does decide how much of Girin's growth
survives a viscous slurry, on idealized columns and on real ones, so the spec's scaling argument is confirmed as the
right variable. But the slurry stops moving with the wave gradually, over two decades of that number, so a switch at any
value misjudges the growth by about a fifth of the inviscid value on average, and at 1.0 the switch sits where a quarter
of the growth still survives.

For the flights so far this changes nothing that has been sprayed: at the slip-regime speeds of d050 at 178.5 s and d100
at 32.5 s surface tension holds the gradient mode back in every sampled column, with or without the floor. It matters
for the phase you expect to spray, under hard deceleration. There, in d050-like columns, the mode first appears while
the slurry is half-way through the transition (a floor Reynolds number of about 0.7 to 2.5 at the Weber threshold), so
the floor will cut the growth to between a sixth and two fifths, not to zero; in d100-like columns it hardly matters. And
in d050-like columns the profile path's viscous gate, not the floor test, is what currently decides whether the mode
sprays: it removes the mode altogether in all twelve scaled d050 columns, where the viscous solver keeps 5 to 74 percent
of it (five of the twelve keep less than a quarter).

The own-speed modes inside the bridge are a separate question that this check opens rather than answers. If they could
break the surface they would be a spraying mechanism the model does not have, active already in the slip regime, at
growth rates of hundreds per second in the d100 columns; if they only stir the bridge they matter for heat and mixing,
not for droplets.

## 5. What you need to decide

**`FLOOR_RE`.** There are three options. The first is to keep 1.0. Nothing changes in any run; columns with a floor
Reynolds number between 1 and about 3.5 keep being treated as fluid (their full inviscid growth) although the viscous
growth is a fifth to a half of it, and columns below 1 as floors (no gradient growth) although a tenth to a quarter
survives. The second is to raise it to 3.5, the half-way point of the idealized columns (2.9 to 5.7 across all series).
At the flights' own speeds this makes five more of the 69 sampled columns floor-like (60 instead of 55) and changes no
spraying, since the mode is held back by surface tension; under hard deceleration it moves the wall to the bridge's base
in about half the d050 columns, and with the present gate it would let four of eighteen scaled columns spray at nearly
the full inviscid growth where the viscous answer is a third, so the flight verdicts get worse (12 of 18 agreeing instead
of 16). The third is to replace the switch, for the gradient mode's growth, by the factor Re_f / (Re_f + 3.5) applied to
the inviscid growth on the whole profile, with Re_f at the first kink found with `KINK_DROP` = 0.2, and to drop the
viscous gate for that mode. That
is a change to the model, not a constant: it needs a spec amendment and code in `dispersion`, `melt_profile` and `spray`,
the wall would still be needed for the release depth and for the gas-driven mode's film depth, and the 3.5 rests on these
columns alone. It matches the viscous solver best, with a mean error of 0.05 on the idealized columns (against 0.20 to
0.28 for the thresholds tried, 0.1 to 30) and 0.14 on the flight columns with all 18 verdicts right (against 0.26 and
16 today); both flight figures rest on four wavenumbers per column, and with the search extended beyond its end they are
0.13 and 0.27.

I recommend the third option, and keeping `FLOOR_RE` at 1.0 until it is built, because raising the threshold on its own
makes the flight columns' verdicts worse under the present gate.

**`KINK_DROP`.** I recommend keeping 0.2, and the evidence for it is now direct. With the present gate it changes no
verdict at `FLOOR_RE` = 1, on the idealized columns or on the scaled flight columns. Without the gate it does: on the
idealized columns at `FLOOR_RE` = 1 the binary test agrees with the viscous verdict on 43, 45, 45, 46 and 46 of 57
columns for `KINK_DROP` of 0.05, 0.1, 0.2, 0.5 and 1, and at `FLOOR_RE` = 10 even with the gate on 52 against 51. More
important, under the third option `KINK_DROP` picks the first kink whose floor Reynolds number the smooth factor uses:
at 0.5 that kink changes in 58 of the 69 sampled columns at their own speeds and in 14 of the 18 scaled ones, and the
factor's agreement with the viscous solver falls from 18 to 13 of 18 (mean error from 0.14 to 0.23); at 0.1 it falls to
14 of 18 (0.28). A limitation: the factor's constant 3.5 was fitted on columns whose kink is unambiguous, and the flight
comparison used 0.2 from the start, so 0.2 has an advantage over the other values tried; still, none of them did better.
Its other effect is on where the walls land: at the flights' own speeds, 0.5 would move 46 of the 69 sampled walls
closer to the surface (median depth from 74 to 37 micrometres) and 0.1 would make 23 floor-like columns fluid-like.
That moves the release depth of the gradient mode and the film depth of the gas-driven mode, but no viscous result here
says which walls are right, because at those speeds the gradient mode does not exist.

**The viscous gate.** This was not one of the two constants, but it decides the outcome in the d050 columns under hard
deceleration. The options are to keep it as it is (it removes the mode in 12 of 18 scaled columns where the viscous
solver keeps 5 to 74 percent), to compute its effective viscosity over the layer above the first kink (which I tested
and which is worse, 12 of 18 verdicts right and a mean error of 0.34), or to drop it for the gradient mode as part of the
third option above. I recommend the last.

**The modes inside the bridge.** You may want a follow-up that computes the surface amplitude of these modes relative to
their internal amplitude, and their sensitivity to the slurry's shear-thinning and to the deceleration's normal
acceleration, before deciding whether they need a place in the spray model. Leaving them as a noted open item costs
nothing now.

## 6. What was verified, what was assumed, and what is unknown

Verified: the solver against Lamb's capillary wave, Girin's profile at two Reynolds numbers (with equal densities and
with a density jump at the kink), Hooper and Boyd's short-wave growth and Charru and Hinch's sign rule with Yih's wave
speed, each well inside its tolerance; the dense solve against scipy's QZ algorithm, and ARPACK against the dense solve
on a fixture column where ARPACK really runs; every validation, four synthetic cases, seven flight cases at their own
speeds and 15 of the 18 scaled flight cases unchanged at doubled resolution (section 3.5; the other three scaled cases
are unchecked); and the collapse of the growth against the floor Reynolds number across six top-layer Reynolds numbers,
two slurry depths and two Weber numbers.

Assumed: a Newtonian disturbance viscosity at each cell's base shear rate (the slurry law is shear-thinning, which would
make the slurry yield more to the wave and the floor weaker); no gravity or deceleration in the disturbance (deceleration
pressing the melt onto the body would steady long waves, not the short ones); no disturbance of the gas; no surface
tension at internal faces; density jumps below 1 percent ignored at faces (which costs under 0.05 percent of the growth
on the five flight cases tried); the surface-condition closure of section 3.1;
the base profile built from cell-centre speeds by the skins' stress-continuity rule; and, for the scaled flight columns,
viscosities held fixed while the speeds were multiplied by 14 to 930, which makes those columns what-ifs, not flight
states. The verdict threshold "keeps at least half its inviscid growth" is my choice; the mean errors do not depend on
it.

Measured once or in one setting: each sweep was run once; the flight sample covers two flight times of two flights;
the curve's constant 3.5 comes from Couette columns with one sharp viscosity step, and the flight columns agree with it
only to about 0.14; the scaled flight columns were solved at four wavenumbers each; the cell-scale broken-line check ran
on five columns. Not recorded in any committed script or data: the extraction's 41 seconds, the two development
observations of section 3.1, and the printed values of the validation table (reproducible from the tests' functions);
every other number in this report comes from the committed sweep outputs, `summary.json` or `checks.json`. Unknown:
whether the modes inside the
bridge can spray; what Albert and Charru (2000) would add on surface tension in two-layer Couette flow; and whether the
straight-film mode found at high Reynolds numbers (section 3.3) is known in the literature.

## Files

- `analysis/kinked_dispersion/orr_sommerfeld.py`: the solver.
- `analysis/kinked_dispersion/floor_sweep.py`: the sampling of flight columns (`extract`), the sweeps (`run`), the
  flight cases' resolution re-solves (`resolution`), the further checks (`checks`) and this report's figures and summary
  (`plot`).
- `tests/test_kinked_dispersion_orr_sommerfeld.py`: the validations; 17 run in the unit tier (about 20 seconds), and
  Girin's fastest wave over all wavenumbers is marked `reference` (about six minutes on a loaded machine).
- `reports/kinked_dispersion_orr_sommerfeld/`: `flight_columns.npz` (the 69 sampled columns, 0.27 MB), `synthetic.csv`,
  `flight.csv`, `resolution_synthetic.csv`, `resolution_flight.csv`, `checks.json`, `timing_synthetic.json`,
  `timing_resolution_flight.json`, the hand-written `timing_flight.json`, `summary.json` and the two figures.
