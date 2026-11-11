#!/bin/zsh
# launch.sh <outdir> <logdir> <label> [reentry_model run arguments...]: one sub-plan 22 Task 5 flight through fly.py,
# run from the planning copy's code/ with no bytecode written, under caffeinate (no idle sleep while it runs) and
# /usr/bin/time -l (its wall time and peak memory at the end of the log); the start, the end, the exit code and the
# power source logged in <logdir>/<label>.log. $PY chooses the interpreter (default: the laptop's drama_env).
here="${0:A:h}"
out="$1"; logs="$2"; label="$3"; shift 3
mkdir -p "$out" "$logs"
log="$logs/$label.log"
py="${PY:-/Users/ashajain/miniforge3/envs/drama_env/bin/python}"
echo "start $(date '+%F %T') $(pmset -g batt | head -1)" > "$log"
cd "$here/../code" && PYTHONDONTWRITEBYTECODE=1 caffeinate -i /usr/bin/time -l "$py" -B "$here/fly.py" "$out" "$label" "$@" >> "$log" 2>&1
echo "exit $?" >> "$log"
echo "end $(date '+%F %T') $(pmset -g batt | head -1)" >> "$log"
