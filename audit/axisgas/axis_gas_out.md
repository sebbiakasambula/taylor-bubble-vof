## A. Largest cell Courant number from the logs (computed at the start of each step)

| run | first step | max to 0.01 s | max to 0.05 s | max to 0.2 s | max to 1 s | first t with Co > 1 |
|---|---|---|---|---|---|---|
| T2_C2_fine | 0.0000 (t 1.20e-04) | 0.010 | 0.012 | 0.020 | 0.031 | none |
| T2_M2f_mompred_fine | 0.0000 (t 1.20e-04) | 0.013 | 0.015 | 0.351 | 0.712 | none |
| T2_A2f_fine | 0.0000 (t 1.20e-04) | 0.012 | 0.016 | 0.370 | 0.592 | none |
| T2_mompred | 0.0000 (t 1.20e-04) | 0.023 | 0.023 | 0.029 | 0.039 | none |
| T2_prod | 0.0000 (t 1.20e-04) | 0.023 | 0.023 | 0.025 | 0.036 | none |
| T2_C1_coarse | 0.0000 (t 1.20e-04) | 0.011 | 0.012 | 0.024 | 0.039 | none |
| T2_M1c_mompred_coarse | 0.0000 (t 1.20e-04) | 0.011 | 0.015 | 0.030 | 0.040 | none |

## B. Fastest gas cells (alpha < 0.5) in the reconstructed fields

Row 0 is the cell row next to the axis. x - nose in D (negative = behind the nose). Ux axial (up is +), Ur radial.

### T2_C2_fine (134 rows from axis to wall, first row centre r 0.0046 D, its height 0.0092 D)
| t | nose x | max gas |U| | row | x - nose | Ux | Ur | cells with gas |U| > 2 m/s, by row (0,1,2,3,4+) | max liquid |U| |
|---|---|---|---|---|---|---|---|---|
| 0.2 | 4.866 | 0.25 | 65 | -3.333 | +0.23 | -0.11 | [0, 0, 0, 0, 0] | 0.12 |
| 1 | 4.985 | 0.83 | 0 | -2.768 | +0.83 | -0.00 | [0, 0, 0, 0, 0] | 0.37 |
| 2 | 5.223 | 1.02 | 0 | -2.768 | +1.02 | -0.00 | [0, 0, 0, 0, 0] | 0.44 |
| 2.6 | 5.342 | 1.04 | 0 | -2.768 | +1.04 | -0.00 | [0, 0, 0, 0, 0] | 0.45 |

### T2_M2f_mompred_fine (134 rows from axis to wall, first row centre r 0.0046 D, its height 0.0092 D)
| t | nose x | max gas |U| | row | x - nose | Ux | Ur | cells with gas |U| > 2 m/s, by row (0,1,2,3,4+) | max liquid |U| |
|---|---|---|---|---|---|---|---|---|
| 0.2 | 4.866 | 7.72 | 4 | -0.804 | -7.72 | +0.01 | [49, 53, 63, 67, 2335] | 0.25 |
| 1 | 5.253 | 9.54 | 0 | -1.667 | -9.54 | +0.00 | [77, 75, 73, 71, 1680] | 0.67 |
| 2 | 5.789 | 25.01 | 3 | -1.875 | +25.01 | -0.00 | [121, 122, 122, 123, 5418] | 2.65 |
| 2.6 | 7.574 | 122.30 | 0 | -1.369 | +122.29 | +0.25 | [105, 105, 105, 105, 2019] | 4.08 |

### T2_A2f_fine (134 rows from axis to wall, first row centre r 0.0046 D, its height 0.0092 D)
| t | nose x | max gas |U| | row | x - nose | Ux | Ur | cells with gas |U| > 2 m/s, by row (0,1,2,3,4+) | max liquid |U| |
|---|---|---|---|---|---|---|---|---|
| 0.2 | 4.866 | 4.41 | 0 | -1.577 | +4.41 | -0.00 | [53, 44, 41, 39, 469] | 0.28 |
| 1 | 5.402 | 7.35 | 4 | -1.250 | -7.35 | +0.00 | [82, 85, 88, 89, 3058] | 0.65 |
| 2 | 6.027 | 10.45 | 3 | -1.994 | +10.45 | -0.01 | [127, 127, 118, 116, 4004] | 1.72 |
| 2.6 | 6.949 | 8.77 | 1 | -1.488 | +8.77 | -0.00 | [152, 152, 151, 151, 2651] | 1.21 |

### T2_mompred (96 rows from axis to wall, first row centre r 0.0064 D, its height 0.0128 D)
| t | nose x | max gas |U| | row | x - nose | Ux | Ur | cells with gas |U| > 2 m/s, by row (0,1,2,3,4+) | max liquid |U| |
|---|---|---|---|---|---|---|---|---|
| 0.2 | 4.854 | 0.40 | 2 | -2.500 | +0.40 | -0.00 | [0, 0, 0, 0, 0] | 0.18 |
| 1 | 5.021 | 1.01 | 0 | -2.708 | +1.01 | +0.00 | [0, 0, 0, 0, 0] | 0.43 |
| 2 | 5.229 | 1.07 | 0 | -2.708 | +1.07 | +0.00 | [0, 0, 0, 0, 0] | 0.46 |
| 2.6 | 5.354 | 1.01 | 0 | -2.750 | +1.01 | -0.00 | [0, 0, 0, 0, 0] | 0.44 |


## C. Onset: largest Courant number in windows of the first 0.3 s

| run | 0.000-0.025 | 0.025-0.050 | 0.050-0.075 | 0.075-0.100 | 0.100-0.125 | 0.125-0.150 | 0.150-0.175 | 0.175-0.200 | 0.200-0.225 | 0.225-0.250 | 0.250-0.275 | 0.275-0.300 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T2_C2_fine | 0.011 | 0.012 | 0.014 | 0.016 | 0.017 | 0.020 | 0.020 | 0.020 | 0.020 | 0.020 | 0.021 | 0.021 |
   - T2_C2_fine: first step with Co max > 0.025 at t = never
   - T2_C2_fine: first step with Co max > 0.050 at t = never
   - T2_C2_fine: first step with Co max > 0.100 at t = never
| T2_M2f_mompred_fine | 0.013 | 0.015 | 0.017 | 0.020 | 0.023 | 0.088 | 0.282 | 0.351 | 0.334 | 0.373 | 0.379 | 0.223 |
   - T2_M2f_mompred_fine: first step with Co max > 0.025 at t = 0.1361 s
   - T2_M2f_mompred_fine: first step with Co max > 0.050 at t = 0.1459 s
   - T2_M2f_mompred_fine: first step with Co max > 0.100 at t = 0.1520 s
| T2_A2f_fine | 0.012 | 0.016 | 0.018 | 0.020 | 0.023 | 0.029 | 0.084 | 0.370 | 0.436 | 0.491 | 0.400 | 0.478 |
   - T2_A2f_fine: first step with Co max > 0.025 at t = 0.1479 s
   - T2_A2f_fine: first step with Co max > 0.050 at t = 0.1587 s
   - T2_A2f_fine: first step with Co max > 0.100 at t = 0.1758 s

## D. Structure of the gas velocity in the main bubble

Sign changes = share of neighbouring cell pairs whose axial velocity has opposite sign (along the axis row inside the bubble, and across the bubble in the column of the fastest cell). Jump = mean |difference| between neighbours over the mean |Ux|. A smooth flow gives few sign changes and a small jump; cell-to-cell noise gives many.

| run | t | gas cells in bubble | share > 2 m/s | median gas |U| | 95th pct | axis row: sign changes, jump | across: sign changes, jump | Ux along the axis row near the peak (10 cells) |
|---|---|---|---|---|---|---|---|---|
| T2_C2_fine | 0.2 | 8043 | 0.0% | 0.13 | 0.23 | 0.01, 0.02 | 0.03, 0.20 | +0.2 +0.2 +0.2 +0.2 +0.2 +0.2 +0.2 +0.2 +0.2 +0.2 |
| T2_C2_fine | 1 | 8058 | 0.0% | 0.36 | 0.81 | 0.00, 0.01 | 0.01, 0.04 | +0.8 +0.8 +0.8 +0.8 +0.8 +0.8 +0.8 +0.8 +0.8 +0.8 |
| T2_C2_fine | 2 | 8021 | 0.0% | 0.41 | 0.96 | 0.00, 0.02 | 0.01, 0.04 | +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 |
| T2_M2f_mompred_fine | 0.2 | 8055 | 31.9% | 1.37 | 5.12 | 0.02, 0.06 | 0.03, 0.07 | -5.7 -6.0 -6.2 -6.3 -6.4 -6.4 -6.4 -6.4 -6.3 -6.2 |
| T2_M2f_mompred_fine | 1 | 8197 | 24.1% | 1.09 | 4.71 | 0.02, 0.08 | 0.03, 0.07 | -7.6 -8.3 -8.9 -9.2 -9.5 -9.5 -9.5 -9.4 -9.1 -8.7 |
| T2_M2f_mompred_fine | 2 | 8376 | 70.5% | 4.43 | 18.67 | 0.00, 0.02 | 0.01, 0.07 | +21.0 +21.4 +21.8 +22.0 +22.2 +22.2 +22.2 +22.1 +21.8 +21.6 |
| T2_A2f_fine | 0.2 | 8052 | 8.0% | 0.62 | 2.33 | 0.05, 0.06 | 0.04, 0.12 | +4.2 +4.2 +4.3 +4.4 +4.4 +4.4 +4.4 +4.3 +4.3 +4.2 |
| T2_A2f_fine | 1 | 8280 | 41.1% | 1.58 | 5.41 | 0.03, 0.07 | 0.03, 0.07 | -6.2 -6.4 -6.6 -6.7 -6.8 -6.8 -6.8 -6.6 -6.4 -6.0 |
| T2_A2f_fine | 2 | 8505 | 52.8% | 2.17 | 8.29 | 0.01, 0.05 | 0.02, 0.07 | +8.7 +8.9 +9.2 +9.5 +9.7 +9.7 +9.6 +9.5 +9.4 +9.2 |

## E. Largest p_rgh initial residual in the step / largest continuity error (sum local), from the logs (14 Nov 22:53)

| run | t=0.0005 | 0.002 | 0.01 | 0.03 | 0.06 | 0.1 | 0.13 | 1.0 |
|---|---|---|---|---|---|---|---|---|
| T2_C1_coarse | 1.0e-03 / 1e-08 | 2.2e-03 / 1e-07 | 3.0e-03 / 5e-07 | 2.0e-03 / 9e-08 | 1.4e-03 / 9e-08 | 1.1e-03 / 1e-07 | 1.2e-03 / 1e-07 | 4.6e-03 / 4e-07 |
| T2_prod | 1.4e-03 / 1e-08 | 2.3e-03 / 1e-07 | 1.4e-03 / 4e-08 | 9.6e-04 / 3e-08 | 8.5e-04 / 3e-08 | 9.4e-04 / 4e-08 | 1.0e-03 / 4e-08 | 4.0e-03 / 1e-07 |
| T2_C2_fine | 2.5e-03 / 1e-08 | 3.2e-03 / 1e-07 | 1.2e-03 / 1e-08 | 8.7e-04 / 1e-08 | 6.0e-04 / 1e-08 | 6.5e-04 / 1e-08 | 7.1e-04 / 1e-08 | 3.3e-03 / 5e-08 |
| T2_M1c_mompred_coarse | 1.2e-03 / 2e-08 | 2.4e-03 / 7e-08 | 3.7e-03 / 6e-07 | 2.5e-03 / 4e-07 | 1.8e-03 / 3e-07 | 1.6e-03 / 2e-07 | 2.0e-03 / 3e-07 | 6.8e-03 / 9e-07 |
| T2_mompred | 1.8e-03 / 3e-08 | 2.7e-03 / 1e-07 | 1.8e-03 / 1e-07 | 1.3e-03 / 1e-07 | 1.4e-03 / 1e-07 | 1.6e-03 / 1e-07 | 1.7e-03 / 1e-07 | 6.0e-03 / 3e-07 |
| T2_M2f_mompred_fine | 3.3e-03 / 4e-08 | 4.1e-03 / 1e-07 | 1.7e-03 / 6e-08 | 1.3e-03 / 4e-08 | 2.1e-03 / 7e-08 | 3.2e-02 / 1e-06 | 2.3e-01 / 1e-05 | 6.4e-01 / 7e-05 |
| T2_A2c_coarse | 1.2e-03 / 2e-08 | 2.3e-03 / 8e-08 | 3.6e-03 / 6e-07 | 2.3e-03 / 4e-07 | 1.6e-03 / 3e-07 | 1.9e-03 / 3e-07 | 2.0e-03 / 3e-07 | 7.0e-03 / 1e-06 |
| T2_A2_plic | 1.8e-03 / 3e-08 | 2.7e-03 / 1e-07 | 1.8e-03 / 1e-07 | 1.3e-03 / 1e-07 | 1.3e-03 / 9e-08 | 1.6e-03 / 1e-07 | 1.9e-03 / 1e-07 | 6.3e-03 / 4e-07 |
| T2_A2f_fine | 3.2e-03 / 4e-08 | 4.0e-03 / 1e-07 | 1.6e-03 / 5e-08 | 1.1e-03 / 3e-08 | 1.6e-03 / 5e-08 | 1.2e-02 / 3e-07 | 9.8e-02 / 4e-06 | 6.7e-01 / 1e-04 |
