## Set 9. Tests at T3 (not in the paper yet), recomputed 14 Nov

### OpenFOAM A5 (A2 settings at T3), body nose

- Fr all 0.3178 (s.e. 0.005689, n 18), -0.8485%, phi 0.0706; matched (nearest the 5 T3_prod times) 0.3151 (n 5), -1.6921%, phi 0.1408
- T3_prod body nose 0.3205 (identical to its axis nose 0.320527); baseline s.e. 0.010495 (3.27% of Fr)
- Fluent against A5: all -11.2654%, matched -10.5040%; against the baseline -12.0183%
- gas outside the main body: A5 mean 0.9933%, T3_prod 4.3775%
- smeared gas (metrics smear, independent of the nose): A5 0.9480%, T3_prod 4.3398%
- median largest gas speed: A5 8.5794, T3_prod 4.7100 m/s
- A5 largest gas drift 0.0000%

### Fluent D1 and D2

Fluent T3 baseline over 5 to 8.4 s: Fr 0.2820 (n 5), film 0.0732, out_b 0.0058%, smear 0.5115%.

| run | Fr all (n) | Fr matched | shift all (%) | shift matched (%) | phi_F all | paired slope of the difference (%) | film | out_b (%) | smear (%) |
|---|---|---|---|---|---|---|---|---|---|
| D1 Green-Gauss node based gradient | 0.2842 (20) | 0.2847 | 0.7882 | 0.9526 | 0.0577 | 0.9526 +- 0.5332 | 0.0734 | 0.0144 | 0.4976 |
  window to 7.79 s: all 0.1589%, matched 0.4102%
  window to 7.11 s: all 1.1104%, matched 0.7596%
| D2 Modified HRIC | 0.2753 (20) | 0.2755 | -2.3669 | -2.3178 | -0.1733 | -2.3178 +- 0.7156 | 0.0722 | 1.1618 | 2.1312 |
  window to 7.79 s: all -3.0249%, matched -3.1703%
  window to 7.11 s: all -2.0177%, matched -1.7320%
