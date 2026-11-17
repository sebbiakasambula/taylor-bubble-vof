# Set 18. Validation uncertainty at T2 (ASME V&V 20 form), 17 Nov

All in per cent of the quantity. 5 to 8.4 s, medium mesh unless stated.

- base three meshes, matched: Fr 0.2266, 0.2184, 0.2116, order 0.6224, R 0.8168, GCI medium 20.4871%, extrapolated 0.1826; film 0.1361, 0.1340, 0.1328, order 1.7206, R 0.5660, GCI medium 2.5048%, extrapolated 0.1313
- fl three meshes, matched: Fr 0.2696, 0.2632, 0.2592, order 1.4422, R 0.6211, GCI medium 4.9573%, extrapolated 0.2527; film 0.1426, 0.1420, 0.1409, order 1.6661, R 1.7599, GCI medium 0.7614%, extrapolated 0.1394
- Fluent film: changes coarse to medium -0.4499%, medium to fine -0.7954%; Roache two-grid from the medium-fine pair: p 2 2.5009%, p 1 5.9968%
- pred two meshes, all rows: Fr coarse 0.2398, medium 0.2298 (+4.3574%), GCI Fs 3 p 2 13.8035%, p 1 33.0642%; film coarse 0.1375, medium 0.1366 (+0.6919%), GCI p 2 2.1920%, p 1 5.2506%
- geo two meshes, all rows: Fr coarse 0.2597, medium 0.2557 (+1.5405%), GCI Fs 3 p 2 4.8801%, p 1 11.6896%; film coarse 0.1414, medium 0.1409 (+0.3811%), GCI p 2 1.2072%, p 1 2.8917%
- 10 degree wedge against the baseline: Fr 0.4681%, film 0.2520% (used as u_wedge for OpenFOAM)
- Fluent Courant 0.25 run (the step stayed at its cap, so a repeat), matched: 0.2617% (not added; inside u_fit)
- references, Fr: Abubakar and Matar 0.2531 (u_D 0.62%); White and Beardmore 0.2671 (u_D 0.90%)
- references, film: Abubakar and Matar 0.1403 (u_D 0.20%); Llewellin et al. 0.1394 (u_D nan%)

## Fr

| configuration | Fr | u_num | u_fit | u_wedge | E vs A&M | u_val vs A&M | E vs W&B | u_val vs W&B |
|---|---|---|---|---|---|---|---|---|
| OpenFOAM baseline | 0.2184 | 20.49 | 0.69 | 0.47 | -13.73 | 20.51 | -18.25 | 20.52 |
| OpenFOAM predictor (p 2) | 0.2298 | 13.80 | 0.21 | 0.47 | -9.20 | 13.83 | -13.96 | 13.84 |
| OpenFOAM predictor (p 1) | 0.2298 | 33.06 | 0.21 | 0.47 | -9.20 | 33.07 | -13.96 | 33.08 |
| OpenFOAM geometric and predictor (p 2) | 0.2557 | 4.88 | 0.23 | 0.47 | +1.05 | 4.95 | -4.25 | 4.99 |
| OpenFOAM geometric and predictor (p 1) | 0.2557 | 11.69 | 0.23 | 0.47 | +1.05 | 11.72 | -4.25 | 11.74 |
| Fluent | 0.2632 | 4.96 | 0.43 | 0.00 | +3.98 | 5.01 | -1.47 | 5.06 |

## Film

| configuration | film | u_num | u_wedge | E vs A&M | u_val vs A&M | E vs Llewellin |
|---|---|---|---|---|---|---|
| OpenFOAM baseline | 0.1340 | 2.50 | 0.25 | -4.48 | 2.53 | -3.86 |
| OpenFOAM predictor (p 2) | 0.1366 | 2.19 | 0.25 | -2.64 | 2.22 | -2.01 |
| OpenFOAM predictor (p 1) | 0.1366 | 5.25 | 0.25 | -2.64 | 5.26 | -2.01 |
| OpenFOAM geometric and predictor (p 2) | 0.1409 | 1.21 | 0.25 | +0.43 | 1.25 | +1.08 |
| OpenFOAM geometric and predictor (p 1) | 0.1409 | 2.89 | 0.25 | +0.43 | 2.91 | +1.08 |
| Fluent (p 2) | 0.1420 | 2.50 | 0.00 | +1.21 | 2.51 | +1.87 |
| Fluent (p 1) | 0.1420 | 6.00 | 0.00 | +1.21 | 6.00 | +1.87 |

Extrapolated (three meshes): OpenFOAM baseline Fr 0.1826 (-27.87% against A&M), film 0.1313 (-6.39%); Fluent Fr 0.2527 (-0.15%)
