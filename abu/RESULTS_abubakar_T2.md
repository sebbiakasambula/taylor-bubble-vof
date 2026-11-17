# Abubakar and Matar (2021), Fig. 6a, as a second reference at T2 (13 Nov 2026)

Source: Habib Abubakar and Omar Matar, "Taylor bubble motion in stagnant and flowing liquids in vertical pipes.
Part I: Steady-states", arXiv 2109.09091, Fig. 6a, p. 13 (`docs/2109.09091v1.pdf`). BibTeX key abubakar2021taylor.

## What they computed

- Sharp-interface FEM (FreeFem++), steady state in a frame moving with the nose, axisymmetric.
- Gas density and viscosity neglected, constant bubble pressure (PDF p. 4), fixed bubble volume (p. 5). This study runs a
  density ratio of 1000 and a viscosity ratio of 100.
- Same groups as here: Ub = ub / sqrt(gD), Nf = rho sqrt(gD^3) / mu, Eo = rho g D^2 / sigma.
- Validated against Bugg and Saad (2002) at Eo 98.33, Nf 88.95: Ub 0.2928 against 0.303 measured, -3.4% (p. 8).
- Their eq. (3.1) is the correlation of Viana et al. (2003), plotted as the lines. The markers are their simulations.

## Digitising

The figure is a raster of 585 by 440 px. Axis calibration from the tick marks (x residual below 0.26 in Eo, y exact).
Markers found by colour after an opening that removes the 1.5 px lines (`digit.py`, output `digit_out.txt`).
Check (`linechk.py`): the plotted Viana lines read off the same way agree with `tb_references.fr_viana` within
+0.0002 to +0.0005 at Eo 80, 120 and 200 for Nf 40, 60, 80. So a marker is good to about +-0.0005.

Markers at Eo = 100 (a marker column exists there, so no interpolation in Eo):

| Nf | Ub (A&M) | Viana | ratio |
|---|---|---|---|
| 40 | 0.2336 | 0.2291 | 1.0194 |
| 60 | 0.2692 | 0.2659 | 1.0125 |
| 80 | 0.2879 | 0.2864 | 1.0054 |

Interpolation to Nf 50: linear 0.2514, quadratic through 40, 60, 80 0.2535, ratio to Viana interpolated and
multiplied by Viana at Nf 50 0.2544. Taken as **0.2531, range 0.2514 to 0.2544**. Viana at T2: 0.2504.

## T2 runs against it (5 to 8.4 s, medium mesh unless stated)

| run | Fr | against A&M 0.2531 (range) | against Viana |
|---|---|---|---|
| OpenFOAM baseline | 0.2184 | -13.7% (-14.2 to -13.1) | -12.8% |
| OpenFOAM + predictor | 0.2298 | -9.2% | -8.2% |
| OpenFOAM A3, predictor + 3 outer loops | 0.2331 | -7.9% | -6.9% |
| OpenFOAM A1, isoAlpha + predictor | 0.2554 | +0.9% | +2.0% |
| OpenFOAM A2, plicRDF + predictor | 0.2557 | +1.0% (+0.5 to +1.7) | +2.1% |
| Fluent | 0.2632 | +4.0% (+3.4 to +4.7) | +5.1% |
| OpenFOAM baseline, fine mesh | 0.2116 | -16.4% | -15.5% |
| Fluent, fine mesh | 0.2592 | +2.4% | +3.5% |
| Fluent, extrapolated (matched rows) | 0.2527 | -0.2% | +0.9% |

## Limits

- A&M is another simulation, not a measurement. Its own deviation from experiment was -3.4% at Eo 98, Nf 89.
- Passive gas against a resolved gas phase here. Steady state against a transient run.
- Interpolated in Nf between 40 and 60 (range 0.0030 wide).
