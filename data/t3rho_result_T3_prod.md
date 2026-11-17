# T3_prod reading (R3.1 to R3.9), window 5 to 8.4 s

Files: m_T3_prod_5to84.csv, m_T3_prod_all_to8.4.csv. Baselines from base_5to84.csv: T3 OpenFOAM Fr 0.320527, Fluent Fr 0.282005, films 0.077342 and 0.073169. Window rule: 4.999999 <= t <= 8.401. Standard error = standard error of the least-squares slope.

## R3.1 Fr_test

- Fr_test = 0.320527 (least-squares nose slope), n = 5, standard error 0.010495 (3.27 % of Fr_test); rows t 5.2643 to 8.2724
- baseline: Fr 0.320527, n 5, standard error 0.010495; Fluent: Fr 0.282005, n 5, standard error 0.001299

## R3.2 to R3.3

- R3.2 dFr = 100 (Fr_test / 0.320527 - 1) = 0.000 %; |dFr| > 2 %: FALSE
- R3.3 phi = (Fr_test - 0.320527) / (0.282005 - 0.320527) = 0.0000
- remaining gap = 100 (0.282005 / Fr_test - 1) = -12.02 %, against -12.02 % before

## R3.4 film (corrected, wall 0.5 cos 2.5 deg), mean over the window

- film = 0.077342 D, valid rows 5 of 5
- change against OpenFOAM 0.077342: +0.000000 D (+0.00 %)
- change against Fluent 0.073169: +0.004173 D (+5.70 %)

## R3.5 and R3.8 field quantities, window rows

| quantity | T3_prod test (n = 5) | T3_prod baseline (n = 5) | Fluent T3 (n = 5) |
|---|---|---|---|
| mean out_b (%) | 3.9306 | 3.9306 | 0.0058 |
| mean wake_b (%) | 2.3187 | 2.3187 | 0.0058 |
| mean smear (%) | 4.3398 | 4.3398 | 0.5115 |
| mean L_b (D) | 3.2964 | 3.2964 | 3.5833 |
| min L_b (D) | 3.0174 | 3.0174 | 3.5423 |
| max L_b (D) | 3.3928 | 3.3928 | 3.6439 |
| median pk_all (m/s) | 4.8755 | 4.8755 | 1.7315 |
| median pk_liq (m/s) | 1.4508 | 1.4508 | 1.1970 |
| median pk_gas (m/s) | 4.7100 | 4.7100 | 1.7315 |

base_5to84.csv, OpenFOAM T3: out_b 3.9306, wake_b 2.3187, smear 4.3398; Fluent T3: out_b 0.0058, wake_b 0.0058, smear 0.5115

## R3.6 gas drift (m_T3_prod_all_to8.4.csv)

- rows 11, t 0.7520 to 8.2724; G first 2.580459525e-02, G last 2.580731531e-02
- drift = 100 (G_last / G_first - 1) = 0.0105 %; largest |G / G_first - 1| over all rows 0.0105 %

## R3.7 refit on the rows nearest to the T3_prod snapshot times in the window

- T3_prod times: 5.2643, 6.0163, 6.7683, 7.5204, 8.2724
- rows used: t = 5.2643, 6.0163, 6.7683, 7.5204, 8.2724 (5 distinct)
- Fr = 0.320527, n = 5, standard error 0.010495; dFr = 0.000 %, phi = 0.0000

## R3.9 wake at the last window row

| run | t | wake_b (%) | wk_fastgas_b (%) | fast cells (U > 0.3) |
|---|---|---|---|---|
| T3_prod test | 8.2724 | 1.7081 | 59.669 | 2601 |
| T3_prod baseline | 8.2724 | 1.7081 | 59.669 | 2601 |
| Fluent T3 | 8.1068 | 0.0058 | 0.000 | 173 |

Per-row values:

| t | nose | tail_b | L_b | film | smear % | out_b % | wake_b % | pk_all |
|---|---|---|---|---|---|---|---|---|
| 5.2643 | 6.5473 | 3.1875 | 3.3598 | 0.0835 | 4.077 | 3.667 | 2.685 | 7.3527 |
| 6.0163 | 6.8299 | 3.8125 | 3.0174 | 0.0713 | 4.333 | 3.934 | 2.425 | 6.7459 |
| 6.7683 | 7.0376 | 3.6875 | 3.3501 | 0.0785 | 4.472 | 4.045 | 2.606 | 4.8755 |
| 7.5204 | 7.2578 | 3.8958 | 3.3619 | 0.0824 | 4.317 | 4.032 | 2.170 | 4.6848 |
| 8.2724 | 7.5386 | 4.1458 | 3.3928 | 0.0710 | 4.501 | 3.975 | 1.708 | 4.7100 |
