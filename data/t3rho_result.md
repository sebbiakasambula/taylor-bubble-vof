# T3_rhoST reading (R3.1 to R3.9), window 5 to 8.4 s

Files: m_T3_rhoST_5to84.csv, m_T3_rhoST_all.csv. Baselines from base_5to84.csv: T3 OpenFOAM Fr 0.320527, Fluent Fr 0.282005, films 0.077342 and 0.073169. Window rule: 4.999999 <= t <= 8.401. Standard error = standard error of the least-squares slope.

## R3.1 Fr_test

- Fr_test = nan (least-squares nose slope), n = 18, standard error nan (nan % of Fr_test); rows t 5.0000 to 8.4000
- baseline: Fr 0.320527, n 5, standard error 0.010495; Fluent: Fr 0.282005, n 5, standard error 0.001299

## R3.2 to R3.3

- R3.2 dFr = 100 (Fr_test / 0.320527 - 1) = nan %; |dFr| > 2 %: FALSE
- R3.3 phi = (Fr_test - 0.320527) / (0.282005 - 0.320527) = nan
- remaining gap = 100 (0.282005 / Fr_test - 1) = nan %, against -12.02 % before

## R3.4 film (corrected, wall 0.5 cos 2.5 deg), mean over the window

- film = nan D, valid rows 0 of 18
- change against OpenFOAM 0.077342: +nan D (+nan %)
- change against Fluent 0.073169: +nan D (+nan %)

## R3.5 and R3.8 field quantities, window rows

| quantity | T3_rhoST test (n = 18) | T3_prod baseline (n = 5) | Fluent T3 (n = 5) |
|---|---|---|---|
| mean out_b (%) | 64.2827 | 3.9306 | 0.0058 |
| mean wake_b (%) | 21.8377 | 2.3187 | 0.0058 |
| mean smear (%) | 10.8434 | 4.3398 | 0.5115 |
| mean L_b (D) | 0.6674 | 3.2964 | 3.5833 |
| min L_b (D) | 0.0411 | 3.0174 | 3.5423 |
| max L_b (D) | 2.4844 | 3.3928 | 3.6439 |
| median pk_all (m/s) | 7.5674 | 4.8755 | 1.7315 |
| median pk_liq (m/s) | 3.3737 | 1.4508 | 1.1970 |
| median pk_gas (m/s) | 7.3103 | 4.7100 | 1.7315 |

base_5to84.csv, OpenFOAM T3: out_b 3.9306, wake_b 2.3187, smear 4.3398; Fluent T3: out_b 0.0058, wake_b 0.0058, smear 0.5115

## R3.6 gas drift (m_T3_rhoST_all.csv)

- rows 42, t 0.2000 to 8.4000; G first 2.580445381e-02, G last 2.580533910e-02
- drift = 100 (G_last / G_first - 1) = 0.0034 %; largest |G / G_first - 1| over all rows 0.0034 %

## R3.7 refit on the rows nearest to the T3_prod snapshot times in the window

- T3_prod times: 5.2643, 6.0163, 6.7683, 7.5204, 8.2724
- rows used: t = 5.2000, 6.0000, 6.8000, 7.6000, 8.2000 (5 distinct)
- Fr = 0.366008, n = 5, standard error 0.326572; dFr = 14.189 %, phi = -1.1806

## R3.9 wake at the last window row

| run | t | wake_b (%) | wk_fastgas_b (%) | fast cells (U > 0.3) |
|---|---|---|---|---|
| T3_rhoST test | 8.4000 | 17.4931 | 66.969 | 3636 |
| T3_prod baseline | 8.2724 | 1.7081 | 59.669 | 2601 |
| Fluent T3 | 8.1068 | 0.0058 | 0.000 | 173 |

Per-row values:

| t | nose | tail_b | L_b | film | smear % | out_b % | wake_b % | pk_all |
|---|---|---|---|---|---|---|---|---|
| 5.0000 | 6.3160 | 6.2708 | 0.0452 | nan | 5.505 | 99.992 | 35.278 | 5.2483 |
| 5.2000 | 6.0296 | 5.8542 | 0.1754 | nan | 5.794 | 99.460 | 51.025 | 9.8327 |
| 5.4000 | 6.3848 | 6.1458 | 0.2390 | nan | 5.838 | 97.204 | 38.794 | 7.8845 |
| 5.6000 | 6.3201 | 6.2708 | 0.0492 | nan | 5.781 | 99.778 | 37.033 | 8.8249 |
| 5.8000 | 3.7612 | 3.2292 | 0.5320 | nan | 5.597 | 26.671 | 24.978 | 6.1707 |
| 6.0000 | 6.5637 | 6.4792 | 0.0845 | nan | 5.725 | 99.975 | 33.656 | 7.7876 |
| 6.2000 | nan | nan | nan | nan | 6.497 | nan | nan | 5.2432 |
| 6.4000 | 6.6733 | 6.3125 | 0.3608 | nan | 8.107 | 94.235 | 34.661 | 7.3472 |
| 6.6000 | 6.6451 | 5.5208 | 1.1243 | nan | 9.324 | 34.368 | 1.203 | 7.2733 |
| 6.8000 | 6.6316 | 5.2292 | 1.4024 | nan | 10.111 | 29.724 | 3.053 | 11.9225 |
| 7.0000 | 6.6302 | 4.1458 | 2.4844 | nan | 11.490 | 26.565 | 18.899 | 25.9813 |
| 7.2000 | 6.7890 | 6.3125 | 0.4765 | nan | 13.631 | 86.346 | 29.306 | 23.0084 |
| 7.4000 | 5.9224 | 4.8125 | 1.1099 | nan | 15.200 | 30.108 | 9.476 | 9.2946 |
| 7.6000 | 5.8096 | 4.6875 | 1.1221 | nan | 16.197 | 24.185 | 5.267 | 4.9731 |
| 7.8000 | 5.4548 | 4.5625 | 0.8923 | nan | 16.839 | 21.250 | 5.063 | 5.2857 |
| 8.0000 | 5.0856 | 4.1875 | 0.8981 | nan | 17.607 | 23.638 | 13.105 | 6.2570 |
| 8.2000 | 7.8952 | 7.8542 | 0.0411 | nan | 17.923 | 99.637 | 12.950 | 5.6650 |
| 8.4000 | 8.0795 | 7.7708 | 0.3087 | nan | 18.016 | 99.671 | 17.493 | 9.5502 |
