## Set 16. T2_M2f_mompred_fine, the momentum predictor on the fine T2 mesh (hint 2)

Ran 14 Nov 11:56 to 21:59, exit 0, to 8.4 s, 43 snapshots. Case files identical to T2_C2_fine except momentumPredictor yes (checked 14 Nov 19:14). Metrics: of_metrics.py (clear 1.5); regions: intact_check.py.

Starting gas G: baseline, fine 0.017740, predictor, fine 0.017740, A2, fine 0.017740. M2f against the fine baseline +0.0000%. Gas drift 0 to 8.4 s: baseline, fine +0.0008%, predictor, fine +0.0030%, A2, fine -0.0004%.

| t | baseline, fine: regions, largest %, axis nose, L_b | predictor, fine: regions, largest %, axis nose, L_b | A2, fine: regions, largest %, axis nose, L_b |
|---|---|---|---|
| 0.2 | 1, 99.84, 4.870, 3.367 | 1, 99.82, 4.880, 3.377 | 1, 99.83, 4.878, 3.375 |
| 1.0 | 1, 99.57, 5.012, 3.420 | 1, 99.50, 5.267, 3.645 | 1, 99.61, 5.406, 3.784 |
| 2.0 | 1, 99.55, 5.224, 3.334 | 1, 99.38, 5.816, 3.867 | 1, 99.52, 6.040, 4.061 |
| 2.4 | 1, 99.54, 5.306, 3.297 | 1, 99.22, 6.857, 4.818 | 1, 99.49, 6.720, 4.652 |
| 2.6 | 1, 99.52, 5.353, 3.284 | 1, 99.24, 7.589, 5.490 | 1, 99.47, 6.968, 4.840 |
| 2.8 | 1, 99.46, 5.391, 3.263 | 2, 80.86, 8.303, 0.699 | 1, 99.40, 7.145, 4.987 |
| 3.0 | 1, 99.53, 5.437, 3.250 | 2, 80.81, 8.934, 0.378 | 1, 99.45, 7.254, 5.037 |
| 3.6 | 1, 99.54, 5.563, 3.197 | 2, 80.83, 6.636, 4.389 | 1, 99.48, 7.464, 5.039 |
| 3.8 | 1, 99.47, 5.602, 3.206 | 3, 65.50, 6.670, 1.000 | 1, 99.49, 7.523, 4.948 |
| 4.0 | 1, 99.49, 5.649, 3.193 | 3, 65.53, 6.704, 0.885 | 1, 99.56, 7.576, 4.764 |
| 5.0 | 1, 99.52, 5.860, 3.285 | 4, 65.52, 6.917, 0.771 | 1, 99.49, 8.470, 3.395 |
| 5.2 | 1, 99.54, 5.898, 3.294 | 4, 65.47, 6.963, 0.787 | 2, 85.54, 9.239, 0.742 |
| 6.0 | 1, 99.47, 6.070, 3.228 | 4, 65.57, 10.567, 0.016 | 2, 85.62, 8.499, 1.460 |
| 7.0 | 1, 99.51, 6.281, 3.260 | 4, 65.47, 11.406, 0.767 | 2, 85.57, 8.806, 3.850 |
| 8.4 | 1, 99.51, 6.577, 3.259 | 4, 65.37, 12.608, 1.462 | 3, 85.64, 8.710, 2.267 |

Axis nose speed before the first split (2.0 to 2.6 s): predictor 2.954 D/s, fine baseline 0.214 D/s.
Film before the split (D), predictor against baseline: t 1.0: 0.1455 / 0.1459, t 1.6: 0.1425 / 0.1457, t 2.0: 0.1440 / 0.1433.
Largest gas speed (m/s), predictor against baseline: t 1.0: 9.542 / 0.835, t 2.0: 25.013 / 1.022, t 2.4: 15.294 / 1.040, t 2.6: 122.295 / 1.041.

Reading (facts only): on the fine mesh the predictor run splits between 2.6 and 2.8 s and again between 3.6 and 3.8 s; the fine baseline stays one bubble; A2 splits once between 5.0 and 5.2 s. No fine Fr exists for the predictor configuration, so its grid study has two points (coarse, medium). Why is not known.
