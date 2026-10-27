"""Builds a Box by hand for the box tests, so they run neither the droplet model nor the Stage A Monte Carlo.

`make_box` makes a Trail of length L holding N_Al aluminium atoms (no emission samples) and a hot start that puts the
whole inventory into one species (`start`) at radius r0, with zero Stage A diagnostics. The ambient is the WACCM
column of the given time of day at altitude z. Every other keyword is a Settings field (for example
reactions=("R6",), t_end_s=1e3, eddy_scale=0.0, freeze_radius=True, chemistry=False, particles=False).
"""
from reentry_model.wake_column.ambient import AmbientColumn
from reentry_model.wake_column.box import Box
from reentry_model.wake_column.settings import Settings
from reentry_model.wake_column.slowing import HotStartResult
from reentry_model.wake_column.species import load_species
from reentry_model.wake_column.trail import Trail

S = load_species()


def make_box(*, z=80.0, time_of_day="night", N_Al=1e12, start="Al", r0=0.05, L=1e9, initial_bins=None,
             **settings_fields):
    trail = Trail(d0_um=1.0, h0_km=z, L_m=L, N_Al=N_Al, m_Al_kg=0.0, m_Zn_kg=0.0, m_Mg_kg=0.0, z_mean_km=z,
                  samples=())
    hot = HotStartResult(fractions={start: 1.0}, r0_m=r0, mean_collisions=0.0, mean_stop_time_s=0.0,
                         share_above_ionisation=0.0, share_above_AlO_bond=0.0, per_sample=())
    amb = AmbientColumn.load("waccm", time_of_day).at(z)
    return Box(trail, hot, amb, time_of_day, Settings(**settings_fields), species=S, initial_bins=initial_bins)
