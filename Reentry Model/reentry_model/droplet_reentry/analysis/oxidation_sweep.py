"""Bounding sweep of the oxidation model: oxygen sticking efficiency x share of the oxide left as a surface skin,
over the 21 release altitudes (90-50 km) of the demise sweep, for a 10 um droplet; plus a size check (1 and 100 um)
at the two extremes of the skin share with full sticking.

sticking 0 is the unoxidised model (the skin share is then irrelevant and run once). f_skin 0 is the Guillen &
Williams (2005) picture, every oxide entrained into the melt and the surface kept fresh; f_skin 1 is a skin that
closes the surface within a few nanometres of growth.

  output/droplet_reentry_output/oxidation_sweep/final_states.csv
"""
import argparse
import csv
import os
import time
from multiprocessing import Pool

from reentry_model.droplet_reentry import OUTPUT
from reentry_model.droplet_reentry import model as dd
from reentry_model.droplet_reentry.analysis.demise_sweep import ALTITUDES

STICKING = (0.0, 0.01, 0.1, 1.0)
SKIN = (0.0, 0.01, 0.1, 1.0)
D_MAIN = 10.0
SIZE_CHECK = dict(d=(1.0, 100.0), h=(90, 78, 66, 54), sticking=1.0, skin=(0.0, 1.0))
OUT = os.path.join(OUTPUT, "oxidation_sweep")


def jobs():
    out = [(D_MAIN, h, 0.0, 0.0) for h in ALTITUDES]
    out += [(D_MAIN, h, g, f) for g in STICKING[1:] for f in SKIN for h in ALTITUDES]
    out += [(d, h, SIZE_CHECK["sticking"], f) for d in SIZE_CHECK["d"] for f in SIZE_CHECK["skin"]
            for h in SIZE_CHECK["h"]]
    return out


def worker(job):
    d, h, g, f = job
    t0 = time.time()
    try:
        s, _ = dd.run_case(d, h, sticking=g, f_skin=f)
    except Exception as exc:
        return dict(d0_um=d, h0_km=h, sticking=g, f_skin=f, end_reason="error: {!r}".format(exc))
    s.update(sticking=g, f_skin=f, wall_s=time.time() - t0)
    if g == 0.0:                       # the unoxidised model: everything vaporised, nothing oxidised
        s.update(vaporised_kg=s["mass_lost_kg"], vaporised_frac=s["mass_lost_frac"], oxide_mass_frac=0.0,
                 Al_oxidised_frac_of_element=0.0, m_ox_skin_kg=0.0, m_ox_entrained_kg=0.0, E_ox_over_E_conv=0.0)
    return s


def main():
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args()
    os.makedirs(OUT, exist_ok=True)
    js = sorted(jobs(), key=lambda j: -j[0])
    t0 = time.time()
    with Pool() as pool:
        out = pool.map(worker, js, chunksize=1)
    print("{} cases in {:.0f} s; end reasons {}".format(len(out), time.time() - t0,
                                                       sorted({str(o.get("end_reason")) for o in out})))
    out.sort(key=lambda r: (r["d0_um"], r["sticking"], r["f_skin"], -r["h0_km"]))
    keys = []
    for r in out:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(OUT, "final_states.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main()
