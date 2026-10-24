#!/bin/bash
# 2026-10-08 runs: four seeds to 120 s (fixed proto3), the Scheil frames flight, its export control, the energy decomposition
P=/Users/ashajain/MIT/particle_formation/prototype; O=/Users/ashajain/MIT/particle_formation/reentry_model_output
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
common="run --diameter 100 --altitude 77.500133 --velocity 7.5 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics --removal girin"
mkdir -p $O/proto3_verification/to120_fixed $O/proto3_verification/scheil_export_control $O/diag_energy
for s in 12345 1 2 3; do
  (cd $P/proto3 && nohup $PY ../rebuild/diag/capped_count.py $O/proto3_verification/to120_fixed/capped_seed$s.json -- $common --material AA7075_range --t-max 120 --seed $s --outdir $O/proto3_verification/to120_fixed > $O/proto3_verification/to120_fixed/V6_seed$s.log 2>&1 &)
done
(cd $P/proto3 && nohup $PY -m reentry_model $common --material AA7075_scheil --deep-runoff off --molten-cascade on --frames-every 1 --outdir $O > $O/proto3_verification/S2_d100_scheil_frames.log 2>&1 &)
(cd $P/proto3_export_control && nohup $PY -m reentry_model $common --material AA7075_scheil --deep-runoff off --molten-cascade on --outdir $O/proto3_verification/scheil_export_control > $O/proto3_verification/S2_d100_scheil_export_control.log 2>&1 &)
for m in scheil range; do
  M=AA7075_$m
  (cd $P/proto3_pre_20261008 && nohup $PY ../rebuild/diag/energy_trace.py $O/diag_energy/${m}_decomp.npz --decomp -- $common --material $M --deep-runoff off --molten-cascade on --outdir $O/diag_energy/${m}_decomp > $O/diag_energy/${m}_decomp.log 2>&1 &)
done
