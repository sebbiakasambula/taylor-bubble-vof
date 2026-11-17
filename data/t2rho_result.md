# T2_rhoST reading (R2.1 to R2.8), window 5 to 8.4 s

Fixed inputs: baseline Fr 0.218357, Fluent Fr 0.263163, baseline film 0.134017, Fluent film 0.142004, T2_mompred Fr 0.229808. Window rule: 4.999999 <= t <= 8.401. Standard error = standard error of the least-squares slope.

## Self-checks

- baseline Fr from m_T2_prod.csv: 0.218357, n 4, s.e. 0.001500 (base_5to84.md: 0.218357, 4, 0.001500)
- Fluent Fr from m_F_T2.csv: 0.263163, n 5, s.e. 0.001132 (base_5to84.md: 0.263163, 5, 0.001132)
- T2_mompred Fr from T2_mompred_5to84.npy: 0.229808 (expected 0.229808), phi 0.2556 (expected 0.2556)
- t2rho.py on T2_rhoST (T2_rhoST_5to84.npy, n 18): Fr 0.153969456

## R2.1 Fr_test

- Fr_test = 0.153969 (least-squares nose slope, m_T2_rhoST.csv), n = 18, standard error 0.000529
- t2rho.py value 0.153969456; difference 0.000e+00, required <= 1e-9: PASS
- rows: t 5.0 to 8.4

## R2.2 to R2.4

- R2.2 dFr = 100 (Fr_test / 0.218357 - 1) = -29.487 %
- R2.3 phi = (Fr_test - 0.218357) / (0.263163 - 0.218357) = -1.4370 (self-test with T2_mompred Fr 0.229808: 0.2556, expected 0.2556)
- R2.4 remaining gap = 100 (0.263163 / Fr_test - 1) = 70.92 %, against 20.52 % before

## R2.5 film (corrected, wall 0.5 cos 2.5 deg), mean over the window

- film = 0.120992 D, valid rows 18 of 18
- change against the baseline 0.134017: -0.013025 D (-9.72 %)
- change against Fluent 0.142004: -0.021012 D (-14.80 %)
- t2rho.py film (uncorrected, wall 0.5): 0.121468; film_raw mean here 0.121468

## R2.6 refit on the 4 rows nearest 5.82, 6.66, 7.49, 8.32 s

- rows used: t = 5.8, 6.6, 7.4, 8.4
- Fr (4 rows) = 0.155598, n = 4, standard error 0.000616
- compare: Fr (18 rows) = 0.153969, standard error 0.000529; baseline (4 rows) 0.218357, standard error 0.001500
- dFr (4 rows) = -28.742 %, phi (4 rows) = -1.4007

## R2.7 gas drift (m_T2_rhoST_all.csv)

- rows 42, t 0.2 to 8.4; G first 1.756143114e-02, G last 1.756160535e-02
- drift = 100 (G_last / G_first - 1) = 0.0010 %; largest |G / G_first - 1| over all rows 0.0010 %

## R2.8 field quantities, 5 to 8.4 s rows

| quantity | T2_rhoST (n = 18) | T2_prod (n = 4) | Fluent T2 (n = 5) |
|---|---|---|---|
| mean L_b (D) | 2.9439 | 3.2215 | 3.3672 |
| mean out_b (%) | 0.0387 | 0.0390 | 0.0747 |
| mean wake_b (%) | 0.0387 | 0.0390 | 0.0747 |
| mean smear (%) | 0.6913 | 0.5844 | 0.6207 |
| median pk_all (m/s) | 0.6690 | 0.8703 | 0.9630 |
| median pk_liq (m/s) | 0.3281 | 0.3941 | 0.4376 |
| median pk_gas (m/s) | 0.6690 | 0.8703 | 0.9630 |

Per-row values, T2_rhoST:

| t | nose | tail_ax | L_b | film | smear % | pk_all |
|---|---|---|---|---|---|---|
| 5.0 | 5.5916 | 2.6702 | 2.9041 | 0.1182 | 0.623 | 0.6677 |
| 5.2 | 5.6253 | 2.6975 | 2.8962 | 0.1199 | 0.628 | 0.6691 |
| 5.4 | 5.6590 | 2.7138 | 2.9298 | 0.1210 | 0.731 | 0.6710 |
| 5.6 | 5.6866 | 2.7384 | 2.9158 | 0.1221 | 0.731 | 0.6747 |
| 5.8 | 5.7150 | 2.7539 | 2.9441 | 0.1225 | 0.795 | 0.6799 |
| 6.0 | 5.7490 | 2.7779 | 2.9365 | 0.1225 | 0.667 | 0.6851 |
| 6.2 | 5.7822 | 2.7956 | 2.9697 | 0.1222 | 0.719 | 0.6890 |
| 6.4 | 5.8082 | 2.8235 | 2.9540 | 0.1216 | 0.668 | 0.6900 |
| 6.6 | 5.8385 | 2.8456 | 2.9843 | 0.1212 | 0.641 | 0.6870 |
| 6.8 | 5.8726 | 2.8789 | 2.9768 | 0.1208 | 0.751 | 0.6796 |
| 7.0 | 5.9051 | 2.9162 | 2.9676 | 0.1206 | 0.742 | 0.6690 |
| 7.2 | 5.9302 | 2.9558 | 2.9511 | 0.1205 | 0.712 | 0.6623 |
| 7.4 | 5.9620 | 2.9950 | 2.9412 | 0.1206 | 0.660 | 0.6625 |
| 7.6 | 5.9962 | 3.0299 | 2.9337 | 0.1207 | 0.698 | 0.6626 |
| 7.8 | 6.0277 | 3.0565 | 2.9652 | 0.1208 | 0.706 | 0.6627 |
| 8.0 | 6.0527 | 3.0876 | 2.9485 | 0.1208 | 0.711 | 0.6628 |
| 8.2 | 6.0856 | 3.1206 | 2.9398 | 0.1209 | 0.630 | 0.6629 |
| 8.4 | 6.1197 | 3.1493 | 2.9322 | 0.1209 | 0.627 | 0.6631 |

## Decision rule

|dFr| = 29.487 %; threshold 2 %: |dFr| > 2 % is TRUE (shift of Fr_test from the baseline is -29.487 %).
Combined with the standard errors: Fr_test s.e. 0.000529 (0.34 % of Fr_test); baseline s.e. 0.001500.

