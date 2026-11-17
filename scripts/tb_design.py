"""Sampling design for the Taylor-bubble sweep (UNUSED, the sweep was dropped 27 Oct 2026), plus the case set
S, T1, T2, T3 in VV_CASES. No solver is run.

The sweep is a scrambled Sobol sequence of 128 points (a power of two) over log Eo, log Nf and
log rho_r. Any leading part of a Sobol sequence of power-of-two length is itself well spread, so if the
sweep runs over budget it can be stopped early at 32 or 64 points and remain a valid design. A Latin
hypercube stopped early leaves gaps.

Each point carries the indicators to look at before committing:
  film      estimated developed film thickness, and cells across it on the medium mesh
  cost      estimated minutes for the case running alone, by the model in mesh_design.py
  u_cap     capillary velocity sigma / mu_l, in units of sqrt(g D). Spurious currents from the
            surface-tension discretisation scale with it, so it is the index of contamination risk.
            It is largest where surface tension is strong AND viscosity is low: low Eo with high Nf.
            The static-bubble test (V1) measures the actual contamination and sets the threshold.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from config import WORK
from mesh_design import NR, case_minutes, cells_across, radial_cells
from tb_params import fluid, morton
from tb_references import FR_DUMITRESCU, film_thickness

RANGES = {"Eo": (10.0, 500.0), "Nf": (10.0, 2000.0), "rho_r": (100.0, 1000.0)}
N_POINTS, SEED = 128, 11
THROUGHPUT_4 = 2.07
OUT = WORK / "design.json"
FACES = radial_cells(NR, 0.25)[0]

# Case set after the re-plan of 27 Oct 2026. V3 to V6 run at T1, the hard point (Eo 300, Nf 2000) is dropped.
VV_CASES = [
    dict(id="S", purpose="static bubble (V1a), worst spurious-current corner", Eo=10.0, Nf=2000.0, rho_r=1000.0),
    dict(id="T1", purpose="inertial regime, V2 to V6 and the cross-code comparison, Fr ref 0.338",
         Eo=300.0, Nf=1000.0, rho_r=1000.0),
    dict(id="T2", purpose="viscous regime, thick film, Fr ref 0.250", Eo=100.0, Nf=50.0, rho_r=1000.0),
    dict(id="T3", purpose="surface tension matters, Eo below 40, Fr ref 0.277 (provisional)",
         Eo=20.0, Nf=500.0, rho_r=1000.0),
]


def point(Eo, Nf, rho_r):
    f = fluid(Eo, Nf, rho_r)
    d = film_thickness(FR_DUMITRESCU, Nf)
    return dict(Eo=Eo, Nf=Nf, rho_r=rho_r, sigma=f["sigma"], mu_l=f["mu_l"], morton=morton(f),
                u_cap=f["sigma"] / f["mu_l"], film=d, cells_in_film=cells_across(FACES, d),
                minutes_plan=case_minutes(Nf),
                minutes_levers=case_minutes(Nf, nx=336, length=14.0, travel=8.0))


# UNUSED since the re-plan of 27 Oct 2026: the sweep is dropped. Kept for reference.
def design():
    u = qmc.Sobol(d=3, scramble=True, seed=SEED).random_base2(int(math.log2(N_POINTS)))
    pts = []
    for i, (a, b, c) in enumerate(u):
        lo = {k: math.log10(v[0]) for k, v in RANGES.items()}
        hi = {k: math.log10(v[1]) for k, v in RANGES.items()}
        Eo = 10 ** (lo["Eo"] + a * (hi["Eo"] - lo["Eo"]))
        Nf = 10 ** (lo["Nf"] + b * (hi["Nf"] - lo["Nf"]))
        rr = 10 ** (lo["rho_r"] + c * (hi["rho_r"] - lo["rho_r"]))
        pts.append(dict(id=i, kind="sobol", **point(Eo, Nf, rr)))
    vv = [dict(v, **point(v["Eo"], v["Nf"], v["rho_r"])) for v in VV_CASES]
    return pts, vv


if __name__ == "__main__":
    # The Sobol design and the sweep cost printout below are unused since the sweep was dropped (27 Oct 2026).
    pts, vv = design()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(dict(ranges=RANGES, seed=SEED, sweep=pts, verification=vv), indent=1))
    u = np.array([p["u_cap"] for p in pts])
    cells = np.array([p["cells_in_film"] for p in pts])
    print(f"{len(pts)} Sobol points written to {OUT}")
    print(f"\ncells across the film: min {cells.min()}, median {int(np.median(cells))} (target 8 to 10)")
    print("\ncapillary velocity sigma/mu_l, the spurious-current risk index (units of sqrt(gD)):")
    print(f"  min {u.min():.2f}   median {np.median(u):.1f}   max {u.max():.0f}")
    for t in (1, 10, 50, 100):
        print(f"  points with u_cap > {t:>3}: {int((u > t).sum()):3d} of {len(pts)}")
    for key, label in (("minutes_plan", "planned settings (pipe 20 D, travel 12 D)"),
                       ("minutes_levers", "with the three levers (pipe 14 D, travel 8 D)")):
        m = np.array([p[key] for p in pts])
        print(f"\nestimated sweep, {label}, 1 outer corrector:")
        print(f"  per case {m.min():.0f} to {m.max():.0f} min alone; all {len(pts)} cases, 4 at a time: "
              f"about {m.sum() / THROUGHPUT_4 / 60:.0f} h of wall time")
        for k in (32, 64):
            print(f"  first {k} points only: about {m[:k].sum() / THROUGHPUT_4 / 60:.0f} h")
    print("\nverification cases (separate from the sweep):")
    for v in vv:
        print(f"  {v['id']:<6} Eo={v['Eo']:5g} Nf={v['Nf']:5g}  u_cap={v['u_cap']:7.1f}  "
              f"film={v['film']:.3f} ({v['cells_in_film']} cells)  {v['purpose']}")
