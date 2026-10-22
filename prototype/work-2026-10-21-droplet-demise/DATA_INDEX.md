# Droplet demise: index of outputs

What every output of this folder is, which are current, how to regenerate them, and what each column means.
Interpretation lives in `REPORT.md` (demise model, drag law, mass loss) and `VERTICAL_WINDS.md` (updrafts,
settling to the ground); sources in `../../citations/`.

**Outputs are not in git.** The repository ignores `prototype/**/runs/` and `prototype/**/*.png`, so the folders
below exist only on the machine that ran them (Asha's Mac, 2026-10-21 to 2026-10-22). On a fresh clone, regenerate
them with the commands below; all scripts run with
`PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python` from this folder.

## Model and settings (all current outputs)

Molten AA7075 droplets released with the parent's velocity (DRAMA `Main Body` descent,
`Generic_Satellite Reentry/output/dmf_output.json`) at 90-50 km in 2 km steps (21 altitudes), US76 atmosphere,
908 K release temperature, emissivity 0.1, accommodation 0.8, evaporation coefficient 1, Al/Zn/Mg evaporation,
Singh et al. (2021) drag (the default from 2026-10-22), regime-based heating. A run ends at the final state:
speed within 1 % of the local terminal velocity and temperature within 1 K of equilibrium. Wind runs then follow
the droplet to the ground (`settle`).

## Folders under `runs/`

| folder | status | contents | regenerate (time) |
|---|---|---|---|
| `demise_sweep/` | current | 15 diameters 1-200 um log-spaced x 21 altitudes, no wind, Singh drag | `$PY run_demise_sweep.py` (about 3 min) |
| `demise_sweep_nano/` | current | 15 diameters 1-500 nm log-spaced x 21 altitudes, no wind, Singh drag | `$PY run_demise_sweep.py --nano` (about 30 s) |
| `demise_sweep_w+5cms_s+0.3mms/` | current | 1-200 um, upward bounding wind (+5 cm/s mesosphere, +0.3 mm/s stratosphere), followed to the ground | `$PY run_demise_sweep.py --wind 5 --wind-strat 0.3` |
| `demise_sweep_w+0cms_s+0mms/` | current | 1-200 um, no wind, followed to the ground | `--wind 0 --wind-strat 0` |
| `demise_sweep_w-5cms_s-0.3mms/` | current | 1-200 um, downward bounding wind, followed to the ground | `--wind -5 --wind-strat -0.3` |
| `demise_sweep_nano_w...mms/` (3) | current | the same three wind cases for 1-500 nm | add `--nano` |
| `*_meso-only/` (6) | superseded | first wind runs: mesospheric wind only, full strength 55-90 km, no stratospheric wind | not regenerable with current code |

Inside the two no-wind folders, files with a `_blend` suffix are the same runs with the earlier composed drag law
(`--drag blend`), kept for comparison; `final_states_mach_blend.csv` is the heating hand-over sensitivity (three
diameters). `sweep*.log` are the run logs. Each folder has `histories/`, one CSV per case,
`demise_d<diameter in um>um_h<altitude>km[_suffix].csv`.

## `final_states.csv`: one row per case

Units are SI unless the name says otherwise.

- **Case and release:** `d0_um` initial diameter [um]; `h0_km` release altitude; `V0_ms` release speed relative to
  the air; `gamma0_deg` release flight-path angle (negative down); `We_release` release Weber number; `skipped`
  True if the Weber pre-check (We >= 12) skipped it; `heating_mode`, `drag_law` settings used; `wall_s`, `nfev`
  run cost.
- **Final state:** `end_reason` (`steady`, `ground`, `vaporised`, `t_max`); `t_final_s` time from release;
  `h_final_km` altitude; `V_final_ms` speed relative to the air; `v_terminal_final_ms` local terminal velocity;
  `T_final_K` temperature; `T_eq_final_K` local equilibrium temperature; `d_final_um` diameter.
- **Mass:** `mass0_kg` initial; `mass_final_kg`; `mass_lost_kg`; `mass_lost_frac` (share of initial mass lost to
  vapour); `lost_<Al|Zn|Mg>_kg` and `lost_<el>_frac_of_element` (share of that element's initial mass lost).
- **Extremes along the flight:** `T_peak_K`; `We_max` (Weber number, maximum along the flight); `peak_decel_g`
  (multiples of 9.80665 m/s2); `boiled` True if the alloy's vapour pressure ever exceeded the ambient pressure.
- **Wind runs only:** `w_meso_ms`, `w_strat_ms` wind case [m/s, + up]; `settle_start_km` altitude where settling
  starts; `v_sink_start_ms` sinking speed over the ground there (terminal velocity minus updraft);
  `updraft_at_start_ms`; `min_net_sink_ms` smallest sinking speed below the start; `end_settle` (`ground` or
  `suspended`); `h_suspended_km` altitude where a suspended droplet hovers (NaN if it lands); `t_fall_s` settling
  time; `t_ground_s` time from release to the ground (inf if suspended); `t_tropopause_s` time from release to
  15 km (inf if held above it).

## `histories/*.csv`: one row per output time of one case

`t_s` time from release; `h_km`, `lat_deg`, `lon_deg` position; `V_ms` speed relative to the air; `T_K`; `d_um`;
`mass_kg`, `mass_frac`; `m_<el>_kg` and `mdot_<el>_kgs` element mass and evaporation rate; `q_conv_W`
(convective heating, + in), `q_rad_W` (radiation out), `q_evap_W` (latent heat out); `w_physics_heating` weight of
the validated PhysicsHeating in the heat input; `vapour_over_ambient` summed vapour pressure over ambient pressure
(above 1 = boiling); `decel_g`; `C_D`; `Kn` (Knudsen number, aerosol convention 2 lambda/D); `Ma`; `speed_ratio`;
`blend_w`; `We`; `rho_air`, `p_air`, `T_air`; `v_terminal_ms`, `T_eq_K` (from the second row on). Histories cover
the flight to the final state only; the settling phase of the wind runs is in `final_states.csv`.

## Figures (`runs/demise_sweep/figures/`, and the nanometre ones in `runs/demise_sweep_nano/figures/`)

| figure | content | script |
|---|---|---|
| 1-3, 5, 6 | final-state maps, mass lost, histories, terminal velocity / altitude / final mass / time, mass-loss percentage | `plot_demise_sweep.py [runs/demise_sweep_nano]` |
| 4_heating_choice_blend | heating hand-over sensitivity (earlier drag law) | from the `_blend` runs |
| 7 | Singh drag coefficient across regimes | `plot_drag_regimes.py` |
| 8_time_to_{ground,tropopause}_w...mms | per wind case: time to the ground / tropopause, held-aloft cells coloured by hover altitude (15-100 km scale) | `plot_wind_settling.py` |
| 9 | time to ground and tropopause against diameter, all sizes, three wind cases | `plot_wind_settling.py` |
| 10 | maximum Weber number of every case against the break-up threshold | `plot_weber.py` |

## Checks

`$PY verify_droplet_demise.py` (model, settling; about 2 min) and `$PY check_singh_drag.py` (drag law limits) must
print ALL PASS.
