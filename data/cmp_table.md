# cmp_table

| Quantity | T1 OF | T1 F | T2 OF | T2 F | T3 OF | T3 F |
|---|---|---|---|---|---|---|
| C1 Fr | 0.3293 | 0.3376 | 0.2182 | 0.2637 | 0.3170 | 0.2829 |
| C2 Fr vs Viana (%) | -2.7 | -0.2 | -12.9 | +5.3 | +14.4 | +2.1 |
| C3 Fr, F vs OF (%) | – | +2.5 | – | +20.8 | – | -10.7 |
| C4 film (D), valid rows in brackets | 0.0656 (24) | 0.0658 (28) | 0.1340 (18) | 0.1421 (20) | 0.0778 (19) | 0.0732 (22) |
| C5 film vs Llewellin (%) | +6.7 | +6.9 | -3.9 | +2.0 | n/a | n/a |
| C6 film, F vs OF (%) | – | +0.2 | – | +6.1 | – | -6.0 |
| C7 L_b mean (D) | 3.03 | 3.35 | 3.22 | 3.37 | 3.00 | 3.58 |
| C7 L_b min (D) | 2.71 | 3.08 | 3.20 | 3.35 | 2.04 | 3.54 |
| C7 L_b max (D) | 3.33 | 3.89 | 3.24 | 3.40 | 3.39 | 3.64 |
| C9 gas outside the body, mean (% of gas) | 14.2 | 6.2 | 0.1 | 0.1 | 6.9 | 0.0 |
| C10 gas in the 2 D wake, mean (% of gas) | 10.8 | 5.3 | 0.0 | 0.1 | 2.3 | 0.0 |
| C11 smeared gas, mean (% of gas) | 8.8 | 1.8 | 0.6 | 0.6 | 7.2 | 0.5 |
| C12 median peak speed, all cells (m/s) | 2.36 | 2.02 | 0.87 | 0.96 | 5.17 | 1.73 |
| C13 median peak speed, liquid cells (m/s) | 1.62 | 1.66 | 0.39 | 0.44 | 2.03 | 1.20 |
| C14 median peak speed, gas cells (m/s) | 2.30 | 2.02 | 0.87 | 0.96 | 5.17 | 1.73 |
| C15 window rows | 24 | 28 | 18 | 20 | 20 | 22 |

Window 5 to 20 s (5 - 1e-6 <= t <= 20.001). Medium mesh: OF = V4_medium, T2_prod, T3_prod; F = Fluent T1, T2, T3. Fr = least-squares slope of the nose against t. Viana Fr 0.338404, 0.250410, 0.277025; Llewellin film 0.061508, 0.139409 (none for T3). Percentages from unrounded values, F vs OF = 100 (F / OF - 1).
