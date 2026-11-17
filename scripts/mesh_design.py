"""Mesh design and cost check for the Taylor bubble, by arithmetic only. No solver is run.

Two questions the project plan leaves open, answered before any CPU is spent:
  1. Does the planned radial mesh put enough cells across the liquid film? (risk 1 in the plan)
  2. What will the sweep cost? (the plan budgets 15 core-hours and 4 to 6 minutes per case)

The cost model is deliberately simple and its assumptions are printed. Its purpose is to tell a factor
of two from a factor of twenty, which decides whether the sweep plan must change. The real number
needs one timed run, which the plan already requires before the sweep.
"""
import math

import numpy as np

from tb_references import FR_DUMITRESCU, film_reynolds, film_thickness

R = 0.5
NR, NX, L = 96, 480, 20.0          # planned baseline mesh: radial cells, axial cells, pipe length in D
ALPHA_CO = 0.2                     # interFoam interface Courant limit from the case plan
TRAVEL = 12.0                      # pipe diameters the nose travels, for transient plus a clean fit
NF_RANGE = (10, 50, 200, 500, 1000, 2000)

# per-step cost, anchored to a measurement: 0.153 s per SIMPLE iteration on 24 000 cells, one case
# alone on this laptop (Project 1 benchmark). Scaled linearly with cells; one PIMPLE outer corrector
# is taken as one such iteration; interface advection adds about 20%.
S_PER_IT_24K = 0.153
THROUGHPUT_4 = 2.07                # measured: four parallel cases give 2.07x the throughput of one


def radial_cells(nr, grading):
    """Cell faces from the axis to the wall. grading = last cell / first cell, below 1 is wall-fine."""
    if abs(grading - 1.0) < 1e-12:
        s = np.full(nr, R / nr)
    else:
        q = grading ** (1.0 / (nr - 1))
        s = (R * (q - 1.0) / (q ** nr - 1.0)) * q ** np.arange(nr)
    return np.concatenate([[0.0], np.cumsum(s)]), s


def cells_across(faces, delta):
    centres = 0.5 * (faces[1:] + faces[:-1])
    return int(np.sum(centres > R - delta))


def case_minutes(Nf, nr=NR, nx=NX, length=L, travel=TRAVEL, outer=1):
    """Estimated wall minutes for one case running alone, by the cost model in this file's docstring."""
    d = film_thickness(FR_DUMITRESCU, Nf)
    u_peak = max(1.5 * Nf * d ** 2 / 3.0, FR_DUMITRESCU)
    dt = ALPHA_CO * (length / nx) / u_peak
    steps = (travel / FR_DUMITRESCU) / dt
    return steps * S_PER_IT_24K * (nr * nx / 24000) * 1.2 * outer / 60


if __name__ == "__main__":
    print("FILM THICKNESS, laminar Nusselt film with volume conservation, at the inertial-limit Fr")
    print(f"{'Nf':>6} {'delta/D':>8} {'Re_film':>8}  regime")
    films = {}
    for Nf in NF_RANGE:
        d = film_thickness(FR_DUMITRESCU, Nf)
        films[Nf] = d
        re = film_reynolds(d, Nf)
        regime = "laminar" if re < 20 else ("wavy laminar" if re < 400 else "TRANSITIONAL, estimate unreliable")
        print(f"{Nf:6g} {d:8.4f} {re:8.1f}  {regime}")
    dmin = films[max(NF_RANGE)]

    print(f"\nCELLS ACROSS THE THINNEST FILM (Nf = {max(NF_RANGE)}, delta = {dmin:.4f} D), {NR} radial cells")
    print(f"{'grading':>8} {'cells in film':>14} {'wall cell':>10} {'axis cell':>10} {'axis/wall':>10}")
    for g in (1.0, 0.5, 0.25, 0.1):
        faces, s = radial_cells(NR, g)
        print(f"{g:8g} {cells_across(faces, dmin):14d} {s[-1]:10.5f} {s[0]:10.5f} {s[0] / s[-1]:10.1f}")
    print("target from the case plan: 8 to 10 cells across the thinnest film")

    dx = L / NX
    cells = NR * NX
    print(f"\nCOST, {NR} x {NX} = {cells} cells, axial cell {dx:.4f} D, alphaCo {ALPHA_CO}")
    print(f"{'Nf':>6} {'U_film':>7} {'dt':>9} {'steps':>8} | per case, minutes, alone:  1 outer  3 outer")
    worst = 0.0
    for Nf in NF_RANGE:
        d = films[Nf]
        u_peak = 1.5 * (Nf * d ** 2 / 3.0)                 # Nusselt film peak = 1.5 x mean
        u_peak = max(u_peak, FR_DUMITRESCU)                 # never below the bubble speed itself
        dt = ALPHA_CO * dx / u_peak
        steps = (TRAVEL / FR_DUMITRESCU) / dt
        c1 = steps * S_PER_IT_24K * (cells / 24000) * 1.2 / 60
        worst = max(worst, c1)
        print(f"{Nf:6g} {u_peak:7.3f} {dt:9.2e} {steps:8.0f} |                            {c1:7.1f}  {3 * c1:7.1f}")

    n_traj = 120
    for outer, label in ((1, "1 outer corrector"), (3, "3 outer correctors, as in the case plan")):
        per = worst * outer
        sweep_h = n_traj * per / 60 / THROUGHPUT_4
        print(f"\nsweep of {n_traj} at the worst case, {label}: {per:.0f} min per case, "
              f"about {sweep_h:.0f} h wall, {sweep_h * 2:.0f} core-hours")
    print("plan budget: 15 core-hours, 4 to 6 minutes per case")
    print("\nASSUMPTIONS: axial Courant limit in the film sets dt; per-step cost scaled from a measured")
    print("steady iteration; the nose region's radial Courant limit is ignored and may make dt smaller.")
