## Set 8. Fluent tests at T2, recomputed 14 Nov

Per-snapshot tables rebuilt from the raw exports in audit/flcopy, identical to fl2/m_F_T2_*.csv. Baseline Fluent T2 over 5 to 8.4 s: Fr 0.2632 (n 5), out_b 0.0747%, smear 0.6207%, film 0.1420.

| run | change | Fr all (n) | Fr matched | largest offset (s) | shift all (%) | shift matched (%) | phi_F all | phi_F matched | vs Abubakar matched (%) | film | out_b (%) | smear (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A6 | Green-Gauss node based gradient | 0.2659 (17) | 0.2655 | 0.0892 | 1.0334 | 0.9034 | -0.0607 | -0.0531 | 4.9152 | 0.1424 | 0.0784 | 0.5179 |
| A7 | Modified HRIC | 0.2584 (17) | 0.2597 | 0.0892 | -1.7911 | -1.3000 | 0.1052 | 0.0764 | 2.6242 | 0.1404 | 0.9873 | 1.8084 |
| A8 | global Courant 0.25 | 0.2637 (17) | 0.2639 | 0.0892 | 0.1925 | 0.2616 | -0.0113 | -0.0154 | 4.2479 | 0.1421 | 0.0797 | 0.6253 |

- A7 matched against the OpenFOAM predictor run (all rows 0.229808): 13.0256%; against it matched (0.230956): 12.4638%
- Fluent baseline against A2 (all rows 0.255745): 2.9005%; Fluent baseline against Abubakar 3.9759%
