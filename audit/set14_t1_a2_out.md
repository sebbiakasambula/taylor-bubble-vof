## Set 14. T1med_A2plic, the A2 configuration at T1 (computed 14 Nov)

Set-up against T1med_mompred: only `application interIsoFoam` and the fvSolution alpha block of T2_A2_plic (plicRDF). Mesh points and initial alpha identical; endTime 10, writeInterval 0.25 as T1med_mompred.

| Fr, 5 to 10 s | Fr | s.e. | n |
|---|---|---|---|
| V4_medium (baseline) | 0.3286 | 0.0008 | 8 |
| T1med_mompred | 0.3295 | 0.0003 | 21 |
| T1med_A2plic, all rows | 0.3161 | 0.0005 | 21 |
| F_T1 (Fluent) | 0.3375 | 0.0006 | 9 |
| T1med_A2plic, matched to the 8 baseline times (max offset 0.1189 s) | 0.3154 | | 8 |

- A2 against the baseline: -3.7910% all rows, -4.0131% matched
- A2 against the predictor run: -4.0584%
- A2 against Fluent: -6.3178% (baseline -2.6263%, predictor -2.3549%)
- gap to Fluent as the paper gives it (Fluent / OF - 1): baseline 2.6972%, predictor 2.4117%, A2 6.7438%
- phi (A2 - OF)/(F - OF): -1.4055 all rows, -1.4879 matched

| Quantity, 5 to 10 s | V4_medium | T1med_mompred | T1med_A2plic | A2 vs baseline | F_T1 |
|---|---|---|---|---|---|
| Film (mean) | 0.0656 | 0.0653 | 0.0638 | -2.6808% | 0.0657 |
| Smeared gas (mean, %) | 4.5871 | 4.6668 | 1.3859 | -69.7872% | 1.4169 |
| L_b (mean) | 3.0474 | 3.1697 | 3.2675 | 7.2240% | 3.4077 |
| out_b (mean, %) | 9.2529 | 8.7833 | 1.5213 | -83.5589% | 4.7856 |
| Largest speed (median) | 2.1349 | 2.5259 | 2.3042 | 7.9293% | 2.0262 |
| Largest liquid speed (median) | 1.6261 | 1.6259 | 1.6263 | 0.0141% | 1.6648 |

- gas drift T1med_A2plic, last against t = 0: -0.0001%
- end time 10.0000 s, 41 snapshots, wall time 3.1481 h
- regions: max 17, max second region 0.4662%, max gas outside the largest after 1 s 2.5928%
