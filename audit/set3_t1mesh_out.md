## Set 3. T1 mesh shape (OpenFOAM and Fluent), recomputed 14 Nov

Per-snapshot tables recomputed from the raw files (audit/m), identical to the stored ones. `cmp_table.py` rerun on them (`code/analysis/audit/cmp/`) gives t1mesh_table.csv, cmp_table.csv, base_5to84.csv and wake_t3.csv identical to the stored tables (largest difference 0). Window 5 to 20 s unless stated.

### t1mesh_table (regenerated), 4 decimal places

| id | statistic | OF stretched | OF square | F stretched | F square |
|---|---|---|---|---|---|
| M1 | Fr | 0.3293 (n 24) | 0.3411 (n 24) | 0.3376 (n 28) | 0.3391 (n 11) |
| M2 | film_mean | 0.0656 (n 24) | 0.0662 (n 15) | 0.0658 (n 28) | 0.0654 (n 11) |
| M3 | L_b_mean | 3.0320 (n 24) | 2.9650 (n 24) | 3.3496 (n 28) | 3.2335 (n 11) |
| M3 | L_b_min | 2.7051 (n 24) | 2.5991 (n 24) | 3.0778 (n 28) | 2.9894 (n 11) |
| M3 | L_b_max | 3.3304 (n 24) | 3.3611 (n 24) | 3.8934 (n 28) | 3.4011 (n 11) |
| M4 | out_b_mean | 14.2113 (n 24) | 15.0328 (n 24) | 6.1992 (n 28) | 7.1925 (n 11) |
| M4 | out_b_max | 20.6390 (n 24) | 22.3332 (n 24) | 8.6558 (n 28) | 10.8649 (n 11) |
| M4 | wake_b_mean | 10.8397 (n 24) | 12.2006 (n 24) | 5.2924 (n 28) | 6.5966 (n 11) |
| M4 | wake_b_max | 14.3855 (n 24) | 16.8889 (n 24) | 6.4312 (n 28) | 9.7931 (n 11) |
| M5 | pk_all_median | 2.3560 (n 24) | 2.3330 (n 24) | 2.0207 (n 28) | 2.0329 (n 11) |
| M5 | pk_liq_median | 1.6192 (n 24) | 1.5955 (n 24) | 1.6612 (n 28) | 1.6619 (n 11) |
| M5 | pk_gas_median | 2.3030 (n 24) | 2.2384 (n 24) | 2.0207 (n 28) | 2.0329 (n 11) |
| M6 | u_core_median | 0.8065 (n 24) | 0.8757 (n 24) | 0.8164 (n 28) | 0.8177 (n 11) |
| M7 | fast_core_median | 0.3936 (n 24) | 1.1851 (n 24) | 0.0392 (n 28) | 0.0095 (n 11) |
| M8 | pk_region_bubble | 23.0000 (n 24) | 17.0000 (n 24) | 23.0000 (n 28) | 9.0000 (n 11) |
| M8 | pk_region_wake | 1.0000 (n 24) | 7.0000 (n 24) | 5.0000 (n 28) | 2.0000 (n 11) |
| M8 | pk_region_ahead | 0.0000 (n 24) | 0.0000 (n 24) | 0.0000 (n 28) | 0.0000 (n 11) |
| M8 | pk_region_beside | 0.0000 (n 24) | 0.0000 (n 24) | 0.0000 (n 28) | 0.0000 (n 11) |
| M8 | pk_region_none | 0.0000 (n 24) | 0.0000 (n 24) | 0.0000 (n 28) | 0.0000 (n 11) |
| M9 | gap05_mean | 0.1082 (n 24) | 0.1112 (n 24) | 0.1100 (n 28) | 0.1098 (n 11) |
| M10 | skirt_b_mean | 0.5972 (n 24) | 0.5447 (n 24) | 0.3780 (n 28) | 0.4915 (n 11) |
| M10 | skirt_b_gt_0.1_share_pct | 91.6667 (n 24) | 95.8333 (n 24) | 96.4286 (n 28) | 100.0000 (n 11) |
| M10 | skirt_b_gt_0.1_count | 22.0000 (n 24) | 23.0000 (n 24) | 27.0000 (n 28) | 11.0000 (n 11) |
| M1 | Fr_sq_vs_str_pct |  | 3.5784 (n 24) |  | 0.4387 (n 11) |
| M3 | L_b_mean_sq_minus_str_D |  | -0.0670 (n 24) |  | -0.1161 (n 11) |
| M3 | L_b_mean_sq_vs_str_pct |  | -2.2110 (n 24) |  | -3.4661 (n 11) |
| M5 | pk_all_sq_vs_str_pct |  | -0.9789 (n 24) |  | 0.6021 (n 11) |
| M5 | pk_liq_sq_vs_str_pct |  | -1.4611 (n 24) |  | 0.0441 (n 11) |
| M5 | pk_gas_sq_vs_str_pct |  | -2.8033 (n 24) |  | 0.6021 (n 11) |
| M11 | t | 19.7004 (n 24) | 19.7004 (n 24) | 19.6629 (n 28) | 19.1868 (n 11) |
| M11 | L_b | 2.7574 (n 24) | 2.9512 (n 24) | 3.3553 (n 28) | 3.1629 (n 11) |
| M11 | L_b_pctrank | 8.3333 (n 24) | 41.6667 (n 24) | 53.5714 (n 28) | 36.3636 (n 11) |
| M11 | out_b | 20.6390 (n 24) | 20.5747 (n 24) | 8.6558 (n 28) | 9.4150 (n 11) |
| M11 | out_b_pctrank | 100.0000 (n 24) | 87.5000 (n 24) | 100.0000 (n 28) | 90.9091 (n 11) |
| M11 | wake_b | 11.7943 (n 24) | 13.3007 (n 24) | 6.4312 (n 28) | 8.3432 (n 11) |
| M11 | wake_b_pctrank | 62.5000 (n 24) | 58.3333 (n 24) | 100.0000 (n 28) | 90.9091 (n 11) |
| M11 | pk_all | 1.7964 (n 24) | 3.0552 (n 24) | 1.9385 (n 28) | 2.6518 (n 11) |
| M11 | pk_all_pctrank | 8.3333 (n 24) | 79.1667 (n 24) | 3.5714 (n 28) | 81.8182 (n 11) |
| M11 | pk_liq | 1.6043 (n 24) | 1.6965 (n 24) | 1.9029 (n 28) | 2.3253 (n 11) |
| M11 | pk_liq_pctrank | 20.8333 (n 24) | 75.0000 (n 24) | 96.4286 (n 28) | 100.0000 (n 11) |
| M11 | pk_gas | 1.7964 (n 24) | 3.0552 (n 24) | 1.9385 (n 28) | 2.6518 (n 11) |
| M11 | pk_gas_pctrank | 8.3333 (n 24) | 79.1667 (n 24) | 3.5714 (n 28) | 81.8182 (n 11) |
| M11 | u_core | 0.7471 (n 24) | 0.8433 (n 24) | 0.8102 (n 28) | 0.8095 (n 11) |
| M11 | u_core_pctrank | 12.5000 (n 24) | 4.1667 (n 24) | 14.2857 (n 28) | 27.2727 (n 11) |
| M11 | fast_core | 0.0000 (n 24) | 1.1547 (n 24) | 0.0000 (n 28) | 0.9427 (n 11) |
| M11 | fast_core_pctrank | 14.5833 (n 24) | 37.5000 (n 24) | 21.4286 (n 28) | 90.9091 (n 11) |
| M11 | skirt_b | 1.0833 (n 24) | 0.5417 (n 24) | 0.2083 (n 28) | 0.5313 (n 11) |
| M11 | skirt_b_pctrank | 91.6667 (n 24) | 58.3333 (n 24) | 17.8571 (n 28) | 54.5455 (n 11) |

### Square mesh in OpenFOAM, whole run

- VA_square Fr 5 s to end 0.3390 (n 56); 10 s to end 0.3385 (n 48)
- V4_fine Fr 5 s to end 0.3267; 10 s to end 0.3269
- Square against fine, 5 s to end: 3.7802%, which is 3.5739 times the fine-grid GCI of 1.0577%
- Square against fine, 10 s to end: 3.5344%
- VA_square mesh: 48 radial by 1920 axial, 92160 cells, dx/dr at wall 1.0000, cells in film 6. Axial cell of V4_medium / VA_square = 4.0000; cells in film V4_medium 22
- V4_medium at 19.700352 s: nose 11.3199, window centre 8.8199, 3 columns, film 0.0646 D (min 0.0645, max 0.0648), cells across 22.8498
- VA_square at 19.700352 s: nose 11.5606, window centre 9.0606, 12 columns, film 0.0644 D (min 0.0642, max 0.0646), cells across 6.1866

### Fluent square mesh

- Gas volume, first snapshot: stretched 0.318675437, square 0.318670277, difference -0.0016%
- T1_square/maximum-velocity-rfile.out: step size over 6533 consecutive rows after step 200: median 0.00239605 s, 10th to 90th percentile 0.00133614 to 0.00311704 s
- T1_square_w_completion/maximum-velocity-rfile_6001.out: step size over 3148 consecutive rows after step 200: median 0.00227858 s, 10th to 90th percentile 0.00104445 to 0.00301699 s
- Fluent against OpenFOAM, 5 to 20 s: stretched 2.5256%, square -0.5822%
- Square against stretched, 5 to 20 s: OpenFOAM 3.5784%, Fluent 0.4387%
- Fluent bubble against OpenFOAM: stretched 10.4755%, square 9.0577%
- gap05 mean, square against stretched: OpenFOAM 2.7628%, Fluent -0.1670%
- out_b mean OpenFOAM / Fluent: stretched 2.2924, square 2.0900; square minus stretched: OpenFOAM 0.8215, Fluent 0.9933 points
