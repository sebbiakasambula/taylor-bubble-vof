# wake_t3

## W1 snapshot

| Quantity | OF 20.305 s | F step 78000 |
|---|---|---|
| wake_ax (% of gas) | 3.98 | 0.00 |
| wake_b (% of gas) | 1.99 | 0.00 |
| out_ax (% of gas) | 10.93 | 0.01 |
| out_b (% of gas) | 8.97 | 0.01 |
| wk_nfast_ax (cells) | 2510 | 315 |
| wk_fastgas_ax (% of gas) | 76.57 | 0.00 |
| wk_nfast_b (cells) | 2507 | 315 |
| wk_fastgas_b (% of gas) | 76.55 | 0.00 |

## W2 mean wake_b over the window

| Quantity | OF | F |
|---|---|---|
| mean wake_b, 5 to 20 s (% of gas) | 2.289 (20 rows) | 0.001 (22 rows) |

## W3 fast cells in the wake containing gas

| Quantity | OF | F |
|---|---|---|
| median wk_fastgas_b over rows with wk_nfast_b >= 10 (%) | 74.0 (20 rows) | 0.0 (22 rows) |

## W4 OpenFOAM T3 over the long run (m_T3_prod_all.csv)

| Row nearest | t (s) | L_b (D) | smear (%) |
|---|---|---|---|
| t = 5 s | 5.264 | 3.36 | 4.1 |
| t = 20 s | 20.305 | 3.17 | 9.5 |
| t = 40 s | 39.858 | 2.52 | 18.1 |

| Quantity | 5 to 20 s | 5 to 40 s |
|---|---|---|
| L_b min (D) | 2.04 | 2.04 |
| L_b max (D) | 3.39 | 3.39 |
| smear min (%) | – | 4.1 |
| smear max (%) | – | 18.2 |
| rows | 20 | 47 |

OF = m_T3_prod.csv (W1 to W3) and m_T3_prod_all.csv (W4); F = Fluent T3. W1: OF at t = 20.305026 s, F at step 78000 (t = 20.2748 s). Windows: 5 - 1e-6 <= t <= 20.001 (W2, W3, and 5 to 20 in W4), <= 40.001 (5 to 40 in W4); the OF row at 20.305 s is outside the 5 to 20 s window.
