## Set 7. OpenFOAM configuration ladder at T2, recomputed 14 Nov

Baseline T2_prod over 5 to 8.4 s: Fr 0.2184 from its 4 snapshots (5.8238, 6.6557, 7.4877, 8.3197), film 0.1340.

| run | change | Fr all (n) | Fr matched | gap to Fluent all (%) | phi all | vs Abubakar all (%) | vs White and Beardmore all (%) | vs Viana all (%) | gap matched (%) | vs Abubakar matched (%) | film all | largest gas drift (%) | starting gas share vs T2_prod (%) | wall time (h) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T2_prod | baseline | 0.2184 (4) | 0.2184 | 20.5197 | 0.0000 | -13.7271 | -18.2490 | -12.8003 | 20.5197 | -13.7271 | 0.1340 | 0.0071 | 0.0000 | 10.6367 |
| T2_mompred | momentum predictor | 0.2298 (18) | 0.2310 | 14.5144 | 0.2556 | -9.2028 | -13.9619 | -8.2274 | 13.9450 | -8.7491 | 0.1366 | 0.0019 | 0.0000 | 2.2847 |
| T2_A3_outer3 | predictor + 3 outer loops | 0.2331 (18) | 0.2338 | 12.9149 | 0.3282 | -7.9166 | -12.7431 | -6.9274 | 12.5764 | -7.6397 | 0.1370 | 0.0006 | 0.0000 | 6.3342 |
| T2_A1_iso | interIsoFoam isoAlpha + predictor | 0.2554 (18) | 0.2557 | 3.0409 | 0.8267 | 0.9074 | -4.3816 | 1.9914 | 2.9235 | 1.0225 | 0.1408 | 0.0000 | 0.0000 | 2.0378 |
| T2_A2_plic | interIsoFoam plicRDF + predictor | 0.2557 (18) | 0.2558 | 2.9004 | 0.8345 | 1.0452 | -4.2511 | 2.1306 | 2.8769 | 1.0683 | 0.1409 | 0.0000 | 0.0000 | 2.9467 |
| T2_A2c_coarse | A2 on the coarse mesh | 0.2597 (18) | 0.2571 | 1.3393 | 0.9224 | 2.6018 | -2.7761 | 3.7040 | 2.3570 | 1.5817 | 0.1414 | 0.0000 | 1.7414 | 0.7417 |
| T2_B1_wedge10 | baseline on a 10 degree wedge | 0.2194 (18) | 0.2196 | 19.9582 | 0.0228 | -13.3232 | -17.8664 | -12.3921 | 19.8384 | -13.2366 | 0.1337 | 0.0015 | 2.4120 | 1.6719 |
| T2_A9_halfdt | baseline at half the time step | 0.2273 (18) | 0.2267 | 15.7963 | 0.1988 | -10.2079 | -14.9144 | -9.2434 | 16.1087 | -10.4495 | 0.1355 | 0.0005 | 0.0000 | 3.2444 |

### Derived

- A3 against the predictor run: 1.4166%; A1 against A2 (all rows): -0.1363%; A2 against the predictor run: phi 0.7776
- Fluent against Abubakar 3.9759%, against White and Beardmore -1.4740%; Abubakar band 0.2514 to 0.2544
- B1 against the baseline: all 0.4681%, matched 0.5685%; A9 against the baseline: all 4.0790%, matched 3.7990%
- A9 against T2_prod over the same 18 snapshot times is not possible (T2_prod wrote 4 in the window); A9 maxDeltaT 0.0004262375 against 0.000852475
- A3 log: 3996 steps with 3 outer loops in 5 to 8.4 s
  p_rgh median initial residual, loops 1 to 3: 0.005426, 0.000617, 0.0001482 (first over third 36.6193)
  Uy median initial residual, loops 1 to 3: 0.001596, 0.000137, 5.949e-05 (first over third 26.8254)
  Ux median initial residual, loops 1 to 3: 0.0002125, 0.0002293, 0.0003446 (first over third 0.6168)
  alpha.water median initial residual, loops 1 to 3: 4.397e-05, 3.557e-05, 3.554e-05 (first over third 1.2372)
