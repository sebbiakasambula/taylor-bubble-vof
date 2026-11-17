## Set 6. Momentum predictor and density weighting, recomputed 14 Nov

### Momentum predictor at T2 (5 to 8.4 s)

- Fr without 0.2184 (n 4), with 0.2298 (n 18): 5.2442%; gap to Fluent (0.263163) 20.5197% to 14.5144%; phi 0.2556 (closed share of the gap 0.2927)
- film without 0.1340, with 0.1366
- wall time T2_mompred 2.2847 h

### Momentum predictor at T1

- medium 5 to 10 s: without 0.3286, with 0.3295, 0.2787%; Fluent 5 to 10 s 0.3375; gap 2.6972% to 2.4117%
- medium film 5 to 10 s: without 0.0656, with 0.0653
- square 10 to 15 s: without 0.3386, with 0.3387, 0.0272%
- square 5 to 10 s: without 0.3489, with 0.3392, -2.7873%
- square against medium over 5 to 10 s: with the predictor 2.9323%, without 6.1787%
- wall time T1med_mompred 1.9953 h, T1sq_mompred 5.9439 h; end times 10.0000 s and 15.0000 s

### Density weighting at T2 (5 to 8.4 s)

- Fr without 0.2184, with 0.1540, -29.4872%; gap to Fluent 20.5197% to 70.9190%; against Viana -38.5131%
- film 0.1340 to 0.1210 (-9.7192%); bubble length 3.2215 to 2.9439 (-8.6177%)
- largest gas drift 0.0010%; out_b mean 0.0387%

### Density weighting at T3, quantities from the per-snapshot tables

- smeared gas over 5 to 8.4 s: without 4.0769 to 4.5007%, with 5.5055 to 18.0156%
- median largest speed 5 to 8.4 s: without 4.8755, with 7.5674 m/s
- largest gas drift with the weight 0.0034%
- with the weight at t = 1.0000: largest speed 1.4306, largest liquid speed 0.9081, nose 4.9874
- with the weight at t = 2.2000: largest speed 9.2972, largest liquid speed 2.1559, nose 5.0790
- with the weight at t = 2.4000: largest speed 7.6872, largest liquid speed 1.4834, nose 5.0007
- with the weight at t = 5.0000: largest speed 5.2483, largest liquid speed 3.7723, nose 6.3160
- with the weight at t = 8.4000: largest speed 9.5502, largest liquid speed 3.7396, nose 8.0795
- without the weight at t = 2.2561: largest speed 3.2272, largest liquid speed 1.2572, nose 5.5559
- without the weight at t = 5.2643: largest speed 7.3527, largest liquid speed 1.8514, nose 6.5473
- nose at t = 0: without 4.9167, with 4.9167

### Build check of the rebuilt solver

- T2_rhoST_check logs ['log.checkMesh', 'log.decomposePar', 'log.interFoam']; steps read 64 and 64; steps with identical deltaT and phase-volume lines: 64
