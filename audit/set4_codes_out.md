## Set 4. OpenFOAM against Fluent on the medium mesh, recomputed 14 Nov

cmp_table.csv and wake_t3.csv regenerated from the recomputed per-snapshot tables (identical to the stored ones, see set 3). Window 5 to 20 s unless stated.

### cmp_table (regenerated), 4 decimal places

| id | case | code | statistic | value | n |
|---|---|---|---|---|---|
| C1 | T1 | OF | Fr | 0.3293 | 24 |
| C2 | T1 | OF | Fr_vs_Viana_pct | -2.6856 | 24 |
| C4 | T1 | OF | film_mean | 0.0656 | 24 |
| C5 | T1 | OF | film_vs_Llewellin_pct | 6.7073 | 24 |
| C7 | T1 | OF | L_b_mean | 3.0320 | 24 |
| C7 | T1 | OF | L_b_min | 2.7051 | 24 |
| C7 | T1 | OF | L_b_max | 3.3304 | 24 |
| C9 | T1 | OF | out_b_mean | 14.2113 | 24 |
| C10 | T1 | OF | wake_b_mean | 10.8397 | 24 |
| C11 | T1 | OF | smear_mean | 8.7578 | 24 |
| C12 | T1 | OF | pk_all_median | 2.3560 | 24 |
| C13 | T1 | OF | pk_liq_median | 1.6192 | 24 |
| C14 | T1 | OF | pk_gas_median | 2.3030 | 24 |
| C15 | T1 | OF | n_window_rows | 24.0000 | 24 |
| C1 | T1 | F | Fr | 0.3376 | 28 |
| C2 | T1 | F | Fr_vs_Viana_pct | -0.2279 | 28 |
| C4 | T1 | F | film_mean | 0.0658 | 28 |
| C5 | T1 | F | film_vs_Llewellin_pct | 6.9055 | 28 |
| C7 | T1 | F | L_b_mean | 3.3496 | 28 |
| C7 | T1 | F | L_b_min | 3.0778 | 28 |
| C7 | T1 | F | L_b_max | 3.8934 | 28 |
| C9 | T1 | F | out_b_mean | 6.1992 | 28 |
| C10 | T1 | F | wake_b_mean | 5.2924 | 28 |
| C11 | T1 | F | smear_mean | 1.8171 | 28 |
| C12 | T1 | F | pk_all_median | 2.0207 | 28 |
| C13 | T1 | F | pk_liq_median | 1.6612 | 28 |
| C14 | T1 | F | pk_gas_median | 2.0207 | 28 |
| C15 | T1 | F | n_window_rows | 28.0000 | 28 |
| C3 | T1 | F_vs_OF | Fr_F_vs_OF_pct | 2.5256 | 28 |
| C6 | T1 | F_vs_OF | film_F_vs_OF_pct | 0.1858 | 28 |
| C1 | T2 | OF | Fr | 0.2182 | 18 |
| C2 | T2 | OF | Fr_vs_Viana_pct | -12.8601 | 18 |
| C4 | T2 | OF | film_mean | 0.1340 | 18 |
| C5 | T2 | OF | film_vs_Llewellin_pct | -3.8687 | 18 |
| C7 | T2 | OF | L_b_mean | 3.2201 | 18 |
| C7 | T2 | OF | L_b_min | 3.1971 | 18 |
| C7 | T2 | OF | L_b_max | 3.2362 | 18 |
| C9 | T2 | OF | out_b_mean | 0.0512 | 18 |
| C10 | T2 | OF | wake_b_mean | 0.0436 | 18 |
| C11 | T2 | OF | smear_mean | 0.5970 | 18 |
| C12 | T2 | OF | pk_all_median | 0.8701 | 18 |
| C13 | T2 | OF | pk_liq_median | 0.3940 | 18 |
| C14 | T2 | OF | pk_gas_median | 0.8701 | 18 |
| C15 | T2 | OF | n_window_rows | 18.0000 | 18 |
| C1 | T2 | F | Fr | 0.2637 | 20 |
| C2 | T2 | F | Fr_vs_Viana_pct | 5.2900 | 20 |
| C4 | T2 | F | film_mean | 0.1421 | 20 |
| C5 | T2 | F | film_vs_Llewellin_pct | 1.9518 | 20 |
| C7 | T2 | F | L_b_mean | 3.3695 | 20 |
| C7 | T2 | F | L_b_min | 3.3482 | 20 |
| C7 | T2 | F | L_b_max | 3.3967 | 20 |
| C9 | T2 | F | out_b_mean | 0.0740 | 20 |
| C10 | T2 | F | wake_b_mean | 0.0740 | 20 |
| C11 | T2 | F | smear_mean | 0.6251 | 20 |
| C12 | T2 | F | pk_all_median | 0.9629 | 20 |
| C13 | T2 | F | pk_liq_median | 0.4376 | 20 |
| C14 | T2 | F | pk_gas_median | 0.9629 | 20 |
| C15 | T2 | F | n_window_rows | 20.0000 | 20 |
| C3 | T2 | F_vs_OF | Fr_F_vs_OF_pct | 20.8288 | 20 |
| C6 | T2 | F_vs_OF | film_F_vs_OF_pct | 6.0548 | 20 |
| C1 | T3 | OF | Fr | 0.3170 | 20 |
| C2 | T3 | OF | Fr_vs_Viana_pct | 14.4244 | 20 |
| C4 | T3 | OF | film_mean | 0.0778 | 19 |
| C7 | T3 | OF | L_b_mean | 3.0003 | 20 |
| C7 | T3 | OF | L_b_min | 2.0399 | 20 |
| C7 | T3 | OF | L_b_max | 3.3928 | 20 |
| C9 | T3 | OF | out_b_mean | 6.8698 | 20 |
| C10 | T3 | OF | wake_b_mean | 2.2894 | 20 |
| C11 | T3 | OF | smear_mean | 7.1542 | 20 |
| C12 | T3 | OF | pk_all_median | 5.1696 | 20 |
| C13 | T3 | OF | pk_liq_median | 2.0256 | 20 |
| C14 | T3 | OF | pk_gas_median | 5.1696 | 20 |
| C15 | T3 | OF | n_window_rows | 20.0000 | 20 |
| C1 | T3 | F | Fr | 0.2829 | 22 |
| C2 | T3 | F | Fr_vs_Viana_pct | 2.1317 | 22 |
| C4 | T3 | F | film_mean | 0.0732 | 22 |
| C7 | T3 | F | L_b_mean | 3.5790 | 22 |
| C7 | T3 | F | L_b_min | 3.5423 | 22 |
| C7 | T3 | F | L_b_max | 3.6439 | 22 |
| C9 | T3 | F | out_b_mean | 0.0058 | 22 |
| C10 | T3 | F | wake_b_mean | 0.0013 | 22 |
| C11 | T3 | F | smear_mean | 0.4793 | 22 |
| C12 | T3 | F | pk_all_median | 1.7316 | 22 |
| C13 | T3 | F | pk_liq_median | 1.1969 | 22 |
| C14 | T3 | F | pk_gas_median | 1.7316 | 22 |
| C15 | T3 | F | n_window_rows | 22.0000 | 22 |
| C3 | T3 | F_vs_OF | Fr_F_vs_OF_pct | -10.7431 | 22 |
| C6 | T3 | F_vs_OF | film_F_vs_OF_pct | -5.9591 | 22 |

### wake_t3 (regenerated), 4 decimal places

| id | item | code | statistic | value | n |
|---|---|---|---|---|---|
| W1 | snapshot | OF | t | 20.3050 | 1 |
| W1 | snapshot | OF | wake_ax | 3.9825 | 1 |
| W1 | snapshot | OF | wake_b | 1.9946 | 1 |
| W1 | snapshot | OF | out_ax | 10.9284 | 1 |
| W1 | snapshot | OF | out_b | 8.9674 | 1 |
| W1 | snapshot | OF | wk_nfast_ax | 2510.0000 | 1 |
| W1 | snapshot | OF | wk_fastgas_ax | 76.5737 | 1 |
| W1 | snapshot | OF | wk_nfast_b | 2507.0000 | 1 |
| W1 | snapshot | OF | wk_fastgas_b | 76.5457 | 1 |
| W1 | snapshot | F | t | 20.2748 | 1 |
| W1 | snapshot | F | wake_ax | 0.0000 | 1 |
| W1 | snapshot | F | wake_b | 0.0000 | 1 |
| W1 | snapshot | F | out_ax | 0.0058 | 1 |
| W1 | snapshot | F | out_b | 0.0058 | 1 |
| W1 | snapshot | F | wk_nfast_ax | 315.0000 | 1 |
| W1 | snapshot | F | wk_fastgas_ax | 0.0000 | 1 |
| W1 | snapshot | F | wk_nfast_b | 315.0000 | 1 |
| W1 | snapshot | F | wk_fastgas_b | 0.0000 | 1 |
| W2 | window_5to20 | OF | wake_b_mean | 2.2894 | 20 |
| W3 | window_5to20 | OF | wk_fastgas_b_median_nfast_ge_10 | 73.9783 | 20 |
| W2 | window_5to20 | F | wake_b_mean | 0.0013 | 22 |
| W3 | window_5to20 | F | wk_fastgas_b_median_nfast_ge_10 | 0.0000 | 22 |
| W4 | nearest_t_5 | OF | t | 5.2643 | 1 |
| W4 | nearest_t_5 | OF | L_b | 3.3598 | 1 |
| W4 | nearest_t_5 | OF | smear | 4.0769 | 1 |
| W4 | nearest_t_20 | OF | t | 20.3050 | 1 |
| W4 | nearest_t_20 | OF | L_b | 3.1685 | 1 |
| W4 | nearest_t_20 | OF | smear | 9.5000 | 1 |
| W4 | nearest_t_40 | OF | t | 39.8580 | 1 |
| W4 | nearest_t_40 | OF | L_b | 2.5187 | 1 |
| W4 | nearest_t_40 | OF | smear | 18.1235 | 1 |
| W4 | window_5to20 | OF | L_b_min | 2.0399 | 20 |
| W4 | window_5to20 | OF | L_b_max | 3.3928 | 20 |
| W4 | window_5to40 | OF | L_b_min | 2.0399 | 47 |
| W4 | window_5to40 | OF | L_b_max | 3.3928 | 47 |
| W4 | window_5to40 | OF | smear_min | 4.0769 | 47 |
| W4 | window_5to40 | OF | smear_max | 18.1647 | 47 |

### Derived here

- T1: Fr OF 0.3293, F 0.3376, F against OF 2.5256% (1.0960 times the medium GCI 2.3043%); against Viana 0.3384: OF -2.6856% (|dFr| 0.0091), F -0.2279% (|dFr| 0.0008)
  film OF 0.0656, F 0.0658, F against OF 0.1858%; against Llewellin 0.0615: OF 6.7073%, F 6.9055%
  bubble length OF 3.0320, F 3.3496, F against OF 10.4755%
- T2: Fr OF 0.2182, F 0.2637, F against OF 20.8288% (9.0391 times the medium GCI 2.3043%); against Viana 0.2504: OF -12.8601% (|dFr| 0.0322), F 5.2900% (|dFr| 0.0132)
  film OF 0.1340, F 0.1421, F against OF 6.0548%; against Llewellin 0.1394: OF -3.8687%, F 1.9518%
  bubble length OF 3.2201, F 3.3695, F against OF 4.6393%
- T3: Fr OF 0.3170, F 0.2829, F against OF -10.7431% (4.6622 times the medium GCI 2.3043%); against Viana 0.2770: OF 14.4244% (|dFr| 0.0400), F 2.1317% (|dFr| 0.0059)
  film OF 0.0778, F 0.0732, F against OF -5.9591%
  bubble length OF 3.0003, F 3.5790, F against OF 19.2867%
- Mesh shape at T1 against the medium GCI: 3.5784 / 2.3043 = 1.5529; codes at T1: 1.0960
- T1 fine mesh: Fr 5 s to end 0.3267, against Viana -3.4666%; film 5 to 20 s 0.0655, against Llewellin 6.4092%
- T2_prod nose 5 s to end (t last 53.2460): Fr 0.2180 (n 58)
- T3_prod nose 5 to 40 s: Fr 0.3150 (n 47), against Viana 13.7198%; nose at t = 39.8580: 17.4562 D
- T3_prod at t = 5.2643: L_b 3.3598 D, smear 4.0769%
- T3_prod at t = 39.8580: L_b 2.5187 D, smear 18.1235%
- T3_prod at t = 5.2643: smear 4.0769%
- T3_prod at t = 39.8580: smear 18.1235%
- T3_prod L_b over 5 to 20 s: 2.0399 to 3.3928 D (n 20)
- T2_prod gas centroid 5 s to end: Fr 0.2179 (n 58)

### Runs table

| case | OF end (s) | OF maxDeltaT | OF wall time (h) | OF largest gas drift (%) | Fluent last snapshot (s) | Fluent final data (s) | Fluent step cap | cap against OF (%) | Fluent largest gas drift (%) |
|---|---|---|---|---|---|---|---|---|---|
| T1 | 39.4007 | 0.000912233 | 5.7300 | 0.0029 | 19.6629 | 20.0042 | 0.00091 | -0.2448 | 1e-05 |
| T2 | 53.2460 | 0.000852475 | 10.6367 | 0.0071 | 20.9160 | 16.9720 | 0.00085 | -0.2903 | 8.61e-05 |
| T3 | 48.1304 | 0.00025709 | 25.0861 | 0.0385 | 22.9788 | 23.1543 | 0.00026 | 1.1319 | 2.8e-06 |

### Fluent continuity residual at the end of each step (scaled, from the residual history in the .dat.h5)

- T1, 20 iterations: 21 step ends, median 0.004117, mean 0.004451, min 0.00261, max 0.00818 (last iteration counter 440400)
- T1, 40 iterations (step 24000): 3 step ends, median 1.375e-05, mean 1.452e-05, min 8.85e-06, max 2.1e-05 (last iteration counter 519600)
- T1 square: 15 step ends, median 0.01998, mean 0.02096, min 0.0168, max 0.0253 (last iteration counter 176780)
- T2: 26 step ends, median 0.000283, mean 0.0002913, min 0.000188, max 0.000449 (last iteration counter 400000)
- T3: 26 step ends, median 3.391e-06, mean 6.538e-06, min 2.86e-06, max 2.13e-05 (last iteration counter 1700000)

### Fluent T1 with 40 iterations per step

- 20 snapshots, t 20.0861 to 21.9880; Fr 0.3379 (s.e. 0.001008, 0.2983% of Fr); 20 iterations over 5 to 20 s: 0.3376; difference 0.0804%
- film: 40 iterations 0.0660 (n 20), 20 iterations 5 to 20 s 0.0658; change 0.3826%
