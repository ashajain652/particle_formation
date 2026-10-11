"""Task 9 measurement campaign: run the queued flights at most SLOTS at a time, logging each run's start, end, exit code
and power source. A run whose directory already holds a JSON result is skipped, so the queue can be restarted.

    python run_queue.py QUEUE_FILE SLOTS

QUEUE_FILE lines: label | working directory | model arguments (after `-m reentry_model`). Blank lines and # comments are
ignored. Each run writes into ../runs/<label>/ and logs to ../runs/<label>.log.
"""
import glob
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(os.path.dirname(HERE), "runs")
PY = "/Users/ashajain/miniforge3/envs/drama_env/bin/python"


def power():
    out = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True).stdout.splitlines()
    return out[0] if out else "unknown"


def stamp(msg):
    line = "{}  {}".format(time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)


queue = []
for raw in open(sys.argv[1]):
    raw = raw.strip()
    if not raw or raw.startswith("#"):
        continue
    label, cwd, args = (part.strip() for part in raw.split("|", 2))
    queue.append((label, cwd, args.split()))
slots = int(sys.argv[2])
os.makedirs(RUNS, exist_ok=True)

running = {}
while queue or running:
    for label, (proc, fh) in list(running.items()):
        if proc.poll() is not None:
            fh.write("exit {}\npower at end: {}\n".format(proc.returncode, power()))
            fh.close()
            stamp("end   {} exit {}  ({})".format(label, proc.returncode, power()))
            del running[label]
    while queue and len(running) < slots:
        label, cwd, args = queue.pop(0)
        out = os.path.join(RUNS, label)
        if glob.glob(os.path.join(out, "*.json")):
            stamp("skip  {} (result present)".format(label))
            continue
        os.makedirs(out, exist_ok=True)
        fh = open(os.path.join(RUNS, label + ".log"), "w")
        fh.write("power at start: {}\ncwd: {}\nargs: {}\n".format(power(), cwd, " ".join(args)))
        fh.flush()
        proc = subprocess.Popen(["caffeinate", "-i", PY, "-m", "reentry_model"] + args + ["--outdir", out],
                                cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
        running[label] = (proc, fh)
        stamp("start {}  ({})".format(label, power()))
    time.sleep(20)
stamp("queue empty")
