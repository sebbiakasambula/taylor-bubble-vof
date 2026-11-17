## Set 5. T2 grid study, OpenFOAM and Fluent, recomputed 14 Nov

r21 = 1.3979, r32 = 1.3954 (square root of the cell-count ratio).

### OpenFOAM

| mesh | rows in window | Fr all | s.e. | Fr matched | largest offset (s) | film all | film matched | gas at first row | gas against medium (%) |
|---|---|---|---|---|---|---|---|---|---|
| coarse | 18 | 0.2272 | 0.000742 | 0.2266 | 0.0877 | 0.1360 | 0.1361 | 0.017867241 | 1.7414 |
| medium | 4 | 0.2184 | 0.001500 | 0.2184 | 0.0000 | 0.1340 | 0.1340 | 0.017561431 | 0.0000 |
| fine | 18 | 0.2114 | 0.000367 | 0.2116 | 0.0877 | 0.1329 | 0.1328 | 0.017739786 | 1.0156 |
- Fr matched: coarse to medium -3.6389%, medium to fine -3.0843%, monotone, order 0.6224, GCI fine 17.1610%, GCI medium 20.4871%, extrapolated 0.1826
- Fr all rows: coarse to medium -3.9105%, medium to fine -3.1937%, monotone, order 0.7421, GCI fine 14.6132%, GCI medium 18.1386%, extrapolated 0.1867
- film matched: 0.1361, 0.1340, 0.1328; changes -1.5273% and -0.8778%
- film all rows: 0.1360, 0.1340, 0.1329; changes -1.4364% and -0.8665%
- with the medium Fr from 5 s to the end (0.2180, n 58): order 0.8872, GCI fine 10.9374%, GCI medium 14.2901%

### Fluent

| mesh | rows in window | Fr all | s.e. | Fr matched | largest offset (s) | film all | film matched | gas at first row | gas against medium (%) |
|---|---|---|---|---|---|---|---|---|---|
| coarse | 17 | 0.2690 | 0.000860 | 0.2696 | 0.0884 | 0.1428 | 0.1426 | 0.203172337 | 0.0141 |
| medium | 5 | 0.2632 | 0.001132 | 0.2632 | 0.0000 | 0.1420 | 0.1420 | 0.203143794 | 0.0000 |
| fine | 17 | 0.2587 | 0.000470 | 0.2592 | 0.0884 | 0.1410 | 0.1409 | 0.203148493 | 0.0023 |
- Fr matched: coarse to medium -2.3878%, medium to fine -1.5194%, monotone, order 1.4422, GCI fine 3.1052%, GCI medium 4.9573%, extrapolated 0.2527
- Fr all rows: coarse to medium -2.1790%, medium to fine -1.6823%, monotone, order 0.8569, GCI fine 6.4332%, GCI medium 8.4279%, extrapolated 0.2454
- film matched: 0.1426, 0.1420, 0.1409; changes -0.4499% and -0.7954%
- film all rows: 0.1428, 0.1420, 0.1410; changes -0.5487% and -0.7419%

- OpenFOAM GCI medium / Fluent GCI medium: matched 4.1327, with the 5 s to end medium fit 2.8826
- OpenFOAM starting gas against Fluent medium: coarse -91.2046%, medium -91.3552%, fine -91.2674%
- Fluent starting gas against Fluent medium: coarse 0.0141%, fine 0.0023%
- T2_C1_coarse maxDeltaT 0.00140711
- T2_prod maxDeltaT 0.000852475
- T2_C2_fine maxDeltaT 0.000515677
- Fluent caps (fl2/run_new.py RUNS): coarse 0.001407, medium 0.00085, fine 0.000516
