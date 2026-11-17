## Set 15. T2_M1c_mompred_coarse, the momentum predictor on the coarse T2 mesh (hint 2)

| run | what | Fr all (s.e., n) | Fr matched (max offset s) | film | smear (%) | L_b | out_b (%) | pk_all median | pk_liq median |
|---|---|---|---|---|---|---|---|---|---|
| T2_C1_coarse | baseline, coarse | 0.2272 (0.0007, 18) | 0.2266 (0.0877) | 0.1360 | 1.0358 | 3.2949 | 0.2014 | 0.9208 | 0.4028 |
| T2_prod | baseline, medium | 0.2184 (0.0015, 4) | 0.2184 (0.0000) | 0.1340 | 0.5844 | 3.2215 | 0.0390 | 0.8703 | 0.3941 |
| T2_M1c_mompred_coarse | predictor, coarse | 0.2398 (0.0007, 18) | 0.2427 (0.0877) | 0.1375 | 0.7747 | 3.3281 | 0.1991 | 0.9407 | 0.4138 |
| T2_mompred | predictor, medium | 0.2298 (0.0005, 18) | 0.2310 (0.0877) | 0.1366 | 0.9743 | 3.2657 | 0.0399 | 0.8956 | 0.4069 |
| T2_A2c_coarse | A2, coarse | 0.2597 (0.0008, 18) | 0.2571 (0.0877) | 0.1414 | 1.0473 | 3.3676 | 0.1588 | 0.9693 | 0.4320 |
| T2_A2_plic | A2, medium | 0.2557 (0.0006, 18) | 0.2558 (0.0877) | 0.1409 | 0.9160 | 3.3166 | 0.0414 | 0.9702 | 0.4309 |

- baseline, medium against coarse: Fr all -3.9105%, matched -3.6389%; film all -1.4364%
- predictor, medium against coarse: Fr all -4.1754%, matched -4.8333%; film all -0.6872%
- A2, medium against coarse: Fr all -1.5171%, matched -0.5054%; film all -0.3796%
- predictor against baseline on the coarse mesh: Fr all 5.5351%, matched 7.0975%; film 1.1553%
- predictor against baseline on the medium mesh: Fr all 5.2442%, matched 5.7700%; film 1.9243%
- M1c gas drift, last against t = 0: 0.0034%; starting gas against T2_C1_coarse: 0.0000%
- M1c end time 8.4000 s, 43 snapshots, wall time 0.5500 h
