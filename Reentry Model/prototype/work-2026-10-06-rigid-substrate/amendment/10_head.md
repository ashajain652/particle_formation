## Amendment of 2026-10-06 — the rigid substrate's four history columns and its frame field

> Part of the rigid-substrate amendment (sub-plan 09's amendment of this date holds the design; facts 78–87 in
> `00-shared-context.md` the measurements). The code below is the tested code, as a diff against the throwaway copy
> described there.

`MELT_COLUMNS` gains `nonrigid_depth_mean_mm`, `slurry_thick_fraction`, `rigid_thin_fraction` and `slurry_held_mass_kg`
(defined in sub-plan 09's amendment; nan with `--rigid-substrate off`, so a run with the switch off has every column it
had before, bit for bit, and four more that are nan). `write_vtk_frame` writes `nonrigid_depth` per patch on the surface
frame: the spray step's own value, carried across the step's deaths by `on_current_surface` like `delta_m`, nan where
the step did not evaluate it. The coupled run's frame test checks the field exists, and the conjugate-depth frame test
that it is the step's own, carried.

