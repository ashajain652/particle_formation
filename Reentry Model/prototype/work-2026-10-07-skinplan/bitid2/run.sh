#!/bin/bash
S=/private/tmp/claude-501/-Users-ashajain-Documents-University-Documents--MIT-Graduate-Work-Research-Space-Sustainability-Particle-Wake-Evolution/930f32b6-e08a-42f0-a4ce-595da1b5cc42/scratchpad/skinplan
BASE="/Users/ashajain/Documents/University Documents /MIT Graduate Work/Research/Space Sustainability/Particle Wake Evolution/prototype/work-2026-10-07-dt-continuum/code"
PY=/Users/ashajain/miniforge3/envs/drama_env/bin/python
ARGS=(run --diameter 100 --velocity 7.5 --altitude 77.500133 --flight-path-angle -0.959331 --atmosphere us76 --thermal fem --melt on --heating physics --t-max 55)
echo "power at start: $(pmset -g batt | head -1)"
(cd "$BASE" && caffeinate -i $PY -m reentry_model "${ARGS[@]}" --outdir $S/bitid2/orig > $S/bitid2/orig.log 2>&1; echo "exit $?" >> $S/bitid2/orig.log) &
(cd $S/rebase/code && caffeinate -i $PY -m reentry_model "${ARGS[@]}" --outdir $S/bitid2/new > $S/bitid2/new.log 2>&1; echo "exit $?" >> $S/bitid2/new.log) &
wait
$PY $S/bitid/compare.py $S/bitid2/orig $S/bitid2/new
echo "power at end: $(pmset -g batt | head -1)"
