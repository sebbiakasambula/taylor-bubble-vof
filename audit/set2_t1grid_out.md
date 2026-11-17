## Set 2. T1 grid study and time step (OpenFOAM), recomputed 14 Nov

Script `code/analysis/audit/set2_t1grid.py`, per-snapshot tables `code/analysis/audit/m/<run>.csv`.

### Meshes

| mesh | radial | axial | cells | dx | dx/dr at wall | cells in film (r = 0.5) | same, wall at 0.4995 | domain length |
|---|---|---|---|---|---|---|---|---|
| V4_coarse | 69 | 343 | 23667 | 0.0583 | 17.4531 | 16 | 16 | 20.0000 |
| V4_medium | 96 | 480 | 46080 | 0.0417 | 17.3409 | 22 | 22 | 20.0000 |
| V4_fine | 134 | 672 | 90048 | 0.0298 | 17.2813 | 30 | 31 | 20.0000 |

### Fr, film and gas per run

| run | snapshots from 5 s | t last | Fr 5 s to end | s.e. | film 5 to 20 s (D) | film rows 5 to 20 s | snapshots after 20 s | with a film | gas drift t = 0 to end (%) | largest drift (%) | nose at end | outlet minus nose (D) | maxDeltaT | wall time (h) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V4_coarse | 56 | 39.4007 | 0.3372 | 0.000049 | 0.0667 | 24 | 32 | 6 | 0.0032 | 0.0032 | 18.2377 | 1.7623 | 0.00149416 | 1.7603 |
| V4_medium | 56 | 39.4007 | 0.3300 | 0.000048 | 0.0656 | 24 | 32 | 13 | 0.0029 | 0.0029 | 17.8279 | 2.1721 | 0.000912233 | 5.7300 |
| V4_fine | 56 | 39.4007 | 0.3267 | 0.000064 | 0.0655 | 24 | 32 | 5 | 0.0032 | 0.0032 | 17.7528 | 2.2472 | 0.000550102 | 24.2314 |
| V3_halfdt | 56 | 39.4007 | 0.3306 | 0.000040 | 0.0657 | 24 | 32 | 12 | 0.0082 | 0.0082 | 17.8468 | 2.1532 | 0.000456117 | 10.9094 |

### Grid convergence (Celik et al. 2008, Fs 1.25)

- r21 = 1.3979, r32 = 1.3954
- Fr: coarse 0.3372, medium 0.3300, fine 0.3267; coarse to medium -2.1489%, medium to fine -1.0058%
- Fr: monotone, order 2.3548, GCI fine 1.0577%, GCI medium from e32 2.3043%, GCI medium from the fine-medium pair 2.3043%, extrapolated 0.3239
- Film: coarse 0.0667, medium 0.0656, fine 0.0655; coarse to medium -1.6010%, medium to fine -0.2793%, order 5.2954

### Fr over 5 s windows (start every 2.5 s, at least 3 snapshots)

- V4_coarse: 7.5000: 0.3358, 10.0000: 0.3366, 12.5000: 0.3371, 15.0000: 0.3368, 17.5000: 0.3369, 20.0000: 0.3377, 22.5000: 0.3382, 25.0000: 0.3381, 27.5000: 0.3373, 30.0000: 0.3367, 32.5000: 0.3373, 35.0000: 0.3367, 37.5000: 0.3375, 40.0000: 0.3383. Largest deviation from the 5 s to end Fr: -0.4303%; last window (centre 40.0000, 4 snapshots) 0.3128%
- V4_medium: 7.5000: 0.3286, 10.0000: 0.3291, 12.5000: 0.3300, 15.0000: 0.3286, 17.5000: 0.3297, 20.0000: 0.3295, 22.5000: 0.3320, 25.0000: 0.3294, 27.5000: 0.3308, 30.0000: 0.3301, 32.5000: 0.3302, 35.0000: 0.3294, 37.5000: 0.3300, 40.0000: 0.3327. Largest deviation from the 5 s to end Fr: 0.8337%; last window (centre 40.0000, 4 snapshots) 0.8337%
- V4_fine: 7.5000: 0.3254, 10.0000: 0.3253, 12.5000: 0.3253, 15.0000: 0.3252, 17.5000: 0.3253, 20.0000: 0.3272, 22.5000: 0.3273, 25.0000: 0.3263, 27.5000: 0.3265, 30.0000: 0.3276, 32.5000: 0.3266, 35.0000: 0.3261, 37.5000: 0.3276, 40.0000: 0.3299. Largest deviation from the 5 s to end Fr: 0.9763%; last window (centre 40.0000, 4 snapshots) 0.9763%

### Time step

- V4_medium maxDeltaT 0.000912233 s, V3_halfdt maxDeltaT 0.000456117 s (ratio 0.5000)
- V4_medium: 43774 steps, 87.3943% at the cap. At the capped steps, 0.2 / largest interface Courant number: median 1.4571, 10th percentile 1.1387, smallest 0.9633
- Wall time on 2 processes (ClockTime): V4_medium 5.7300 h, V3_halfdt 10.9094 h
- Fr V3_halfdt 0.3306 against V4_medium 0.3300: 0.1760%
