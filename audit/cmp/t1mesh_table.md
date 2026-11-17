# t1mesh_table

## M1 to M10

| Quantity | OF stretched | OF square | F stretched | F square |
|---|---|---|---|---|
| M1 Fr | 0.3293 | 0.3411 | 0.3376 | 0.3391 |
| M1 Fr, square vs stretched (%) | – | +3.6 | – | +0.4 |
| M2 film (D), valid rows in brackets | 0.0656 (24) | 0.0662 (15) | 0.0658 (28) | 0.0654 (11) |
| M3 L_b mean (D) | 3.03 | 2.96 | 3.35 | 3.23 |
| M3 L_b min (D) | 2.71 | 2.60 | 3.08 | 2.99 |
| M3 L_b max (D) | 3.33 | 3.36 | 3.89 | 3.40 |
| M3 L_b mean, square minus stretched (D) | – | -0.07 | – | -0.12 |
| M3 L_b mean, square vs stretched (%) | – | -2.2 | – | -3.5 |
| M4 out_b mean (% of gas) | 14.2 | 15.0 | 6.2 | 7.2 |
| M4 out_b max (% of gas) | 20.6 | 22.3 | 8.7 | 10.9 |
| M4 wake_b mean (% of gas) | 10.8 | 12.2 | 5.3 | 6.6 |
| M4 wake_b max (% of gas) | 14.4 | 16.9 | 6.4 | 9.8 |
| M5 median pk_all (m/s) | 2.36 | 2.33 | 2.02 | 2.03 |
| M5 median pk_liq (m/s) | 1.62 | 1.60 | 1.66 | 1.66 |
| M5 median pk_gas (m/s) | 2.30 | 2.24 | 2.02 | 2.03 |
| M5 pk_all, square vs stretched (%) | – | -1.0 | – | +0.6 |
| M5 pk_liq, square vs stretched (%) | – | -1.5 | – | 0.0 |
| M5 pk_gas, square vs stretched (%) | – | -2.8 | – | +0.6 |
| M6 median u_core (m/s) | 0.81 (24) | 0.88 (24) | 0.82 (28) | 0.82 (11) |
| M7 median fast_core (% of core gas) | 0.4 (24) | 1.2 (24) | 0.0 (28) | 0.0 (11) |
| M8 pk_region counts | bubble 23, wake 1 | bubble 17, wake 7 | bubble 23, wake 5 | bubble 9, wake 2 |
| M9 mean gap05 (D) | 0.1082 | 0.1112 | 0.1100 | 0.1098 |
| M10 mean skirt_b (D) | 0.60 | 0.54 | 0.38 | 0.49 |
| M10 share of rows with skirt_b > 0.1 (%, rows) | 91.7 (22 of 24) | 95.8 (23 of 24) | 96.4 (27 of 28) | 100.0 (11 of 11) |

## M11 figure snapshots, value (percentile rank)

| Quantity | OF stretched | OF square | F stretched | F square |
|---|---|---|---|---|
| time (s) | 19.7004 | 19.7004 | 19.6629 | 19.1868 |
| code snapshot | t = 19.700352 s | t = 19.700352 s | step 21645 | step 8190 |
| L_b (D) | 2.76 (8) | 2.95 (42) | 3.36 (54) | 3.16 (36) |
| out_b (% of gas) | 20.6 (100) | 20.6 (88) | 8.7 (100) | 9.4 (91) |
| wake_b (% of gas) | 11.8 (63) | 13.3 (58) | 6.4 (100) | 8.3 (91) |
| pk_all (m/s) | 1.80 (8) | 3.06 (79) | 1.94 (4) | 2.65 (82) |
| pk_liq (m/s) | 1.60 (21) | 1.70 (75) | 1.90 (96) | 2.33 (100) |
| pk_gas (m/s) | 1.80 (8) | 3.06 (79) | 1.94 (4) | 2.65 (82) |
| u_core (m/s) | 0.75 (12) | 0.84 (4) | 0.81 (14) | 0.81 (27) |
| fast_core (%) | 0.0 (15) | 1.2 (38) | 0.0 (21) | 0.9 (91) |
| skirt_b (D) | 1.08 (92) | 0.54 (58) | 0.21 (18) | 0.53 (55) |

Window 5 to 20 s. OF stretched = V4_medium, OF square = VA_square, F stretched = Fluent T1_Case, F square = Fluent T1 square (with completion). Fr = least-squares slope of the nose against t. Square vs stretched = 100 (square / stretched - 1) from unrounded values. u_core and fast_core are medians over the rows where they are defined (count in brackets). M11: value of the snapshot used in the figures, with its percentile rank (scipy percentileofscore, kind 'rank', 0 to 100) among the window rows of the same run, in brackets.
