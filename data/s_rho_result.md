# B11: V1a_static (standard interFoam) against S_rhoST (density-weighted surface tension)

Numbers only. All values are computed by `scripts/b11_srho.py` (per-snapshot table `b11_snap.csv`, report generator `scripts/b11_report.py`). Cell values come from the written fields of `<P2_RUNS>/<run>/<t>/{alpha.water,U,p}` and `0/{C,V}`, read with `scal`/`vec` and `tb_post` exactly as `scripts/snap.py` does (liquid = alpha > 0.99, gas = alpha < 0.01, interface = the rest; `tb_post.effective_radius`, `tb_post.laplace_jump(Cx, Cy, alpha, p, 2.0, 0.25)` with the default margin 0.3; sigma = 0.0999 from transportProperties in both runs). maxSpeed monitor = column 'max' of `postProcessing/maxSpeed/0/fieldMinMax.dat`; liquidVolume = `postProcessing/liquidVolume/0/volFieldValue.dat`.

Case check (diff of the two case folders, whitespace and line endings ignored): the only differences are `application interFoamRhoST` in system/controlDict and `rhoWeightedSurfaceTension yes;` in the PIMPLE dict of system/fvSolution; constant/ is identical, 0/C and 0/alpha.water are byte-identical. Both ran 23 061 steps to 3.95482 s with the same time step (monitor spacing 1.7142e-3 to 1.7151e-3 s per 10 steps), so the arithmetic mean of monitor samples equals the time mean (the trapezoid mean agrees to 3 digits, listed below). Mesh: 96 x 96 = 9216 cells, x uniform (dx = 0.04167), r stretched (dr = 0.0061 at r = 0.25), cell aspect ratio at the interface about 6.8. Boundaries: x = 0 wall, x = 4 open (totalPressure 0, pressureInletOutletVelocity), tube wall no-slip.

## 1. Reproduction of the paper row for V1a_static (to t = 1 s)

| Quantity | Paper (Table, stretched row) | My value, V1a_static | Source |
|---|---|---|---|
| Mean maximum speed, 0.3 to 1 s | 7.27 m/s | **7.268 m/s** (n = 408 samples; trapezoid time mean 7.269) | maxSpeed monitor, samples from t = 0.3017 to 0.9997 s |
| Liquid peak, snapshots to 1 s | 0.41 to 0.60 m/s | **0.408 to 0.605 m/s** (5 snapshots, t = 0.198 to 0.989 s) | max of abs(U) over cells with alpha > 0.99 |
| Pressure jump, snapshots to 1 s | 0.72 to 0.78 Pa (Laplace 0.799) | **0.715 to 0.782 Pa** (Laplace 2 sigma / Reff = 0.7987 Pa) | tb_post.laplace_jump |

Match: 7.268 rounds to 7.27; 0.408 to 0.605 rounds to 0.41 to 0.60; 0.7150 to 0.7822 rounds to 0.72 to 0.78 (0.71501 is just above the 0.715 rounding edge, so this end of the range is sensitive to the last digit). The method reproduces the paper row, so the same code was used on S_rhoST below.

## 2. Same three quantities, same window (to t = 1 s, mean over 0.3 to 1 s)

| Quantity | V1a_static | S_rhoST | S_rhoST / V1a_static |
|---|---|---|---|
| Mean maximum speed, 0.3 to 1 s (m/s) | 7.268 | 2.454 | 0.338 |
| Max of the maximum speed in that window (m/s) | 10.195 | 8.025 | 0.787 |
| Liquid peak, 5 snapshots (m/s) | 0.408 to 0.605 (mean 0.473) | 0.355 to 0.888 (mean 0.618) | mean 1.307; min 0.871; max 1.469 |
| Gas peak, 5 snapshots (m/s) | 5.747 to 8.098 (mean 6.806) | 1.970 to 2.932 (mean 2.536) | mean 0.373 |
| Pressure jump, 5 snapshots (Pa) | 0.715 to 0.782 (mean 0.733, -8.3% of Laplace) | 0.570 to 0.775 (mean 0.653, -18.3% of Laplace) | mean 0.891; min 0.797; max 0.991 |

## 3. Whole run

maxSpeed monitor: t = 0.3 to 3.9546 s (2131 samples each). Snapshots: the 19 written snapshots with t >= 0.3 s (0.395 to 3.955 s) are the main rows, and the 20 snapshots from 0.198 s (the set behind the paper's 0.728 Pa) are given in a second block. The pressure jump of S_rhoST is undefined (the function returns None: no gas cell inside 0.7 R of (x = 2, r = 0)) in the 7 snapshots from t = 2.768 s on, because the bubble has moved away from x = 2 (see item 5 and the Expert note). Its jump statistics therefore use 12 snapshots (t >= 0.3) or 13 (all 20), against 19 or 20 for V1a. A like-for-like row on the same 12 snapshots is added.

| Quantity | V1a_static | S_rhoST | S_rhoST / V1a_static |
|---|---|---|---|
| maxSpeed, time mean 0.3 to 3.95 s (m/s) | 7.558 (trapezoid 7.558) | 2.864 (trapezoid 2.864) | 0.379 |
| maxSpeed, max over 0.3 to 3.95 s (m/s) | 13.313 | 8.025 | 0.603 |
| maxSpeed, max over the whole run incl. start-up (m/s, at t) | 45.515 (t = 0.0325) | 8.025 (t = 0.3788) | 0.176 |
| maxSpeed, time mean per window (m/s) | 0.3-1 s: 7.268; 1-2 s: 7.687; 2-3 s: 7.459; 3-3.96 s: 7.740 | 0.3-1 s: 2.454; 1-2 s: 2.943; 2-3 s: 2.990; 3-3.96 s: 2.948 | 0.338, 0.383, 0.401, 0.381 |
| Liquid peak, mean (min to max), snapshots t >= 0.3 s, n = 19 (m/s) | 0.545 (0.383 to 0.793) | 0.854 (0.355 to 1.175) | mean 1.568; min 0.928; max 1.481 |
| Gas peak, mean (min to max), snapshots t >= 0.3 s, n = 19 (m/s) | 7.329 (4.809 to 9.532) | 3.028 (1.889 to 5.273) | mean 0.413; min 0.393; max 0.553 |
| Pressure jump, mean (min to max), fixed centre x = 2, snapshots t >= 0.3 s, n = 19 (Pa) | 0.729 (0.667 to 0.822), n = 19, -8.8% of Laplace | 0.586 (0.476 to 0.775), n = 12, -26.6% of Laplace | mean 0.804; min 0.713; max 0.942 |
| Liquid peak, mean (min to max), all 20 snapshots from 0.198 s (m/s) | 0.538 (0.383 to 0.793) | 0.846 (0.355 to 1.175) | mean 1.573; min 0.928; max 1.481 |
| Gas peak, mean (min to max), all 20 snapshots from 0.198 s (m/s) | 7.257 (4.809 to 9.532) | 2.975 (1.889 to 5.273) | mean 0.410; min 0.393; max 0.553 |
| Pressure jump, mean (min to max), fixed centre x = 2, all 20 snapshots from 0.198 s (Pa) | 0.728 (0.667 to 0.822), n = 20, -8.8% of Laplace | 0.585 (0.476 to 0.775), n = 13, -26.8% of Laplace | mean 0.803; min 0.713; max 0.942 |
| Pressure jump on the same 12 snapshots (t = 0.395 to 2.571 s) where S_rhoST is defined (Pa) | 0.724 (0.667 to 0.819), sd 0.047 | 0.586 (0.476 to 0.775), sd 0.084 | mean 0.809; min 0.713; max 0.947 |
| Supplement, NOT the paper definition: same function with the centre moved to the gas centroid, t >= 0.3 s, n = 19 (Pa) | 0.729 (0.666 to 0.823), -8.8% | 0.495 (0.220 to 0.775), -38.1% | mean 0.679 |
| Laplace value 2 sigma / Reff, mean over snapshots t >= 0.3 s (Pa) | 0.7987 | 0.7987 | |

## 4. Where the largest speed sits in the last snapshot (t = 3.95482 s)

| Quantity | V1a_static | S_rhoST |
|---|---|---|
| Largest speed (m/s) | 8.097 | 1.889 |
| Class of that cell | gas (alpha<0.01), alpha = 2.18e-05 | gas (alpha<0.01), alpha = 4.99e-06 |
| Position x, r | 2.0208, 0.1855 | 1.1458, 0.2191 |
| Velocity components (Ux, Ur) (m/s) | (-8.001, 1.245) | (-1.833, 0.455) |
| Largest speed in gas / interface / liquid (m/s) | 8.097 / 3.869 / 0.628 | 1.889 / 1.707 / 1.043 |
| Mean abs(U) over interface cells, unweighted (m/s) (n cells) | 0.468 (n = 120) | 0.480 (n = 190) |
| Same, volume-weighted (m/s) | 0.496 | 0.566 |
| Mean abs(U) over interface cells, mean of the 20 snapshots (m/s); range of cell count | 0.379; n = 102 to 124 | 0.464; n = 93 to 202 |
| Monitor at its last sample (t = 3.9546 s): max speed at (x, r) | 8.142 at (2.0208, 0.1855) | 1.884 at (1.1458, 0.2191) |

The maximum is in the gas, next to the interface (r about 0.19 to 0.22 against a bubble radius of 0.25), in both runs. The monitor values differ slightly from the snapshot values only because the last monitor sample is 1.7e-4 s earlier than the last write.

## 5. Shape

| Quantity | V1a_static | S_rhoST |
|---|---|---|
| Reff at t = 0 (from 0/alpha.water) (m) | 0.250166 | 0.250166 |
| Reff at first written snapshot, t = 0.1977 s (m) | 0.250166 | 0.250166 |
| Reff at end, t = 3.9548 s (m) | 0.250174 | 0.250171 |
| Reff, min to max over the 20 snapshots (m) | 0.250166 to 0.250174 | 0.250166 to 0.250171 |
| Smeared gas share at t = 0 (%) | 0.000 | 0.000 |
| Smeared gas share at t = 0.1977 s (%) | 3.656 | 3.300 |
| Smeared gas share at end, t = 3.9548 s (%) | 3.526 | 7.387 |
| Smeared gas share, mean (min to max) over 20 snapshots (%) | 3.541 (2.682 to 4.487) | 4.405 (3.010 to 7.387) |
| Smeared share, mean of the last 5 snapshots, t = 3.17 to 3.95 s (%) | 3.592 | 6.087 |
| Gas centroid x at end (m); start = 2.0000 | 1.9372 | 1.1106 |
| Gas centroid x, min to max over the run (m) | 1.9372 to 2.0226 | 1.1106 to 2.0073 |
| Centroid velocity, linear fit over snapshots t >= 2.3 s (m/s) | -0.042 | -0.407 |
| Extent of cells with alpha < 0.5 in x, min to max over snapshots (m) | 1.646 to 2.354 | 0.938 to 2.479 |
| Largest r of cells with alpha < 0.5, min to max over snapshots (m) | 0.213 to 0.256 | 0.164 to 0.285 |

Reff is computed from the gas volume (tb_post.effective_radius), so it only tracks volume (its snapshot-to-snapshot change is the gas-volume drift of item 6) and says nothing about shape. The smeared share starts at 0 because setFields marks whole cells. Centroid and extent rows come from `b11_srho.py` (the extent of alpha < 0.5 uses cell centres, so values are quantised to the cell size, 0.0417 in x and about 0.006 in r); `b11_drift.py` and `b11_ux.py` (logs in this folder if kept) were used for the cross-checks quoted in the Expert note.

## 6. Gas volume drift (wedge volume, m^3)

| Quantity | V1a_static | S_rhoST |
|---|---|---|
| Gas volume at t = 0, snapshot 0/alpha.water | 9.11417733e-04 | 9.11417733e-04 |
| Gas volume at 3.95482 s, snapshot | 9.11506424e-04 | 9.11473723e-04 |
| Drift, snapshot last vs t = 0 (%) | +0.00973 | +0.00614 |
| Drift max abs over the 20 snapshots (%) | 0.00973 | 0.00614 |
| liquidVolume monitor (gas = 4.35778720e-02 minus volIntegrate(alpha.water)): drift of last (3.9546 s) vs first (0.0033 s) sample, n = 1153 (%) | +0.00973 | +0.00612 |
| Same, max abs over all samples (%) | 0.00973 | 0.00612 |
| Monitor last vs t = 0 snapshot (%) | +0.00973 | +0.00614 |

Both runs gain a trace of gas (monotonic, 0.010% and 0.006% in 4 s); neither is a mass-conservation problem. The two sources agree to 4 digits.

## Per-snapshot numbers (both runs)

| t (s) | V1a liq | S liq | V1a gas | S gas | V1a jump | S jump (x=2) | S jump (centroid) | V1a xc | S xc | V1a smeared % | S smeared % | V1a mean itf | S mean itf |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1977 | 0.411 | 0.696 | 5.893 | 1.970 | 0.720 | 0.570 | 0.570 | 2.0000 | 2.0000 | 3.656 | 3.300 | 0.486 | 0.508 |
| 0.3955 | 0.408 | 0.662 | 5.747 | 2.912 | 0.782 | 0.602 | 0.602 | 1.9997 | 2.0001 | 3.746 | 4.971 | 0.366 | 0.436 |
| 0.5932 | 0.507 | 0.355 | 7.278 | 2.175 | 0.715 | 0.703 | 0.703 | 1.9996 | 2.0007 | 2.682 | 3.882 | 0.434 | 0.323 |
| 0.7910 | 0.434 | 0.489 | 7.015 | 2.692 | 0.717 | 0.775 | 0.775 | 2.0022 | 2.0020 | 4.113 | 3.295 | 0.373 | 0.246 |
| 0.9887 | 0.605 | 0.888 | 8.098 | 2.932 | 0.728 | 0.614 | 0.614 | 2.0088 | 2.0042 | 3.840 | 3.053 | 0.396 | 0.619 |
| 1.1864 | 0.536 | 0.754 | 6.044 | 2.402 | 0.781 | 0.476 | 0.476 | 2.0186 | 2.0061 | 3.606 | 3.010 | 0.283 | 0.463 |
| 1.3842 | 0.545 | 0.664 | 8.327 | 3.167 | 0.667 | 0.504 | 0.504 | 2.0217 | 2.0073 | 2.886 | 4.474 | 0.426 | 0.523 |
| 1.5819 | 0.478 | 0.844 | 5.876 | 4.654 | 0.671 | 0.515 | 0.515 | 2.0226 | 2.0011 | 3.834 | 4.991 | 0.300 | 0.392 |
| 1.7797 | 0.383 | 0.918 | 4.809 | 3.151 | 0.703 | 0.591 | 0.591 | 2.0223 | 1.9772 | 4.487 | 3.167 | 0.363 | 0.317 |
| 1.9774 | 0.522 | 0.773 | 8.930 | 3.242 | 0.727 | 0.552 | 0.545 | 2.0184 | 1.9274 | 3.224 | 3.852 | 0.398 | 0.384 |
| 2.1752 | 0.632 | 0.660 | 5.787 | 3.109 | 0.819 | 0.584 | 0.568 | 2.0101 | 1.8536 | 3.133 | 3.781 | 0.304 | 0.425 |
| 2.3729 | 0.572 | 0.845 | 6.749 | 2.609 | 0.680 | 0.579 | 0.571 | 2.0020 | 1.7591 | 3.369 | 3.791 | 0.406 | 0.473 |
| 2.5706 | 0.451 | 1.058 | 6.421 | 2.242 | 0.699 | 0.540 | 0.536 | 1.9947 | 1.6542 | 3.085 | 4.619 | 0.374 | 0.473 |
| 2.7684 | 0.468 | 1.175 | 8.154 | 5.273 | 0.725 | None | 0.485 | 1.9868 | 1.5513 | 3.848 | 3.261 | 0.391 | 0.529 |
| 2.9661 | 0.645 | 1.106 | 7.040 | 2.395 | 0.822 | None | 0.357 | 1.9744 | 1.4579 | 3.360 | 4.229 | 0.370 | 0.558 |
| 3.1639 | 0.628 | 1.124 | 9.532 | 3.322 | 0.685 | None | 0.388 | 1.9649 | 1.3782 | 3.337 | 3.910 | 0.365 | 0.506 |
| 3.3616 | 0.515 | 0.958 | 9.517 | 2.698 | 0.747 | None | 0.392 | 1.9576 | 1.3058 | 3.302 | 5.766 | 0.333 | 0.586 |
| 3.5593 | 0.596 | 0.922 | 9.111 | 4.122 | 0.722 | None | 0.305 | 1.9506 | 1.2337 | 3.776 | 6.980 | 0.400 | 0.499 |
| 3.7571 | 0.793 | 0.986 | 6.715 | 2.543 | 0.749 | None | 0.253 | 1.9448 | 1.1723 | 4.018 | 6.391 | 0.338 | 0.548 |
| 3.9548 | 0.628 | 1.043 | 8.097 | 1.889 | 0.705 | None | 0.220 | 1.9372 | 1.1106 | 3.526 | 7.387 | 0.468 | 0.480 |

## Expert note

Facts first (each from the tables above), then hypotheses, kept apart. Nothing here is verified beyond the numbers.

Facts.

1. Maximum speed: S_rhoST is lower, by a factor 0.34 (0.3 to 1 s, 2.454 against 7.268 m/s) and 0.38 over the whole run (2.864 against 7.558 m/s). The reduction sits in the gas: gas peak mean 3.028 against 7.329 m/s (ratio 0.41, snapshots t >= 0.3 s). Neither run decays: S_rhoST rises from 2.45 (0.3 to 1 s) to 2.94 to 2.99 m/s afterwards (windows 1 to 2, 2 to 3, 3 to 3.96 s); V1a stays at 7.3 to 7.7 m/s. The start-up peak is 45.5 m/s in V1a (t = 0.0325 s) and does not occur in S_rhoST (largest value of the whole run 8.025 m/s at t = 0.379 s).
2. The liquid is not quieter, it is faster: liquid peak mean 0.854 against 0.545 m/s (ratio 1.57, t >= 0.3 s), largest single snapshot 1.175 against 0.793 m/s. To 1 s the mean is 0.618 against 0.473 m/s (1.31). The mean abs(U) over interface cells is also not lower: 0.464 against 0.379 m/s (mean of 20 snapshots), 0.480 against 0.468 m/s in the last snapshot. So the 3x reduction in max speed is a reduction of the gas maximum only; by the paper's liquid-peak measure S_rhoST is worse, not better.
3. S_rhoST has a net translation of the bubble. The gas centroid goes from x = 2.0000 to 1.1106 (-0.889 m, 3.6 bubble radii) at about -0.41 m/s from t = 2.3 s on (linear fit over 2.3 to 3.95 s); the gas-volume-weighted mean Ux in the snapshots is -0.35 to -0.54 m/s between 1.98 and 3.56 s and -0.243 m/s at the end (`b11_ux.log`). V1a moves 0.063 m in total (centroid range 1.937 to 2.023 m; gas-weighted Ux -0.01 to -0.07 m/s after 1.8 s), about 14 times less. There is no gravity and no imposed flow, so the physical expectation is zero. The drift starts after t = 1.6 s (centroid 2.001 at 1.58 s, 1.977 at 1.78 s, 1.927 at 1.98 s). The mesh is uniform in x, so the drift does not move the bubble into coarser cells.
4. The pressure-jump test is not usable for S_rhoST as it stands. With the fixed centre x = 2 and R = 0.25 the function returns None from t = 2.768 s on (no gas cell within 0.7 R of the centre; bubble centroid 1.55 and moving). On the bubble-at-rest part of the record (8 snapshots to 1.58 s, centroid within 0.007 of x = 2) the jump is already lower and much more scattered: S_rhoST mean 0.595 (sd 0.103, range 0.476 to 0.775) against V1a mean 0.723 (sd 0.043, range 0.667 to 0.782) Pa, Laplace 0.799. With the centre moved to the centroid the S_rhoST jump falls to 0.220 Pa in the last snapshot (supplement row; not the paper definition, because the function assumes a sphere and a liquid ring up to 2 R that the tube wall at r = 0.5 truncates once the bubble is deformed and moving). I would not quote any S_rhoST jump after about 1.6 s.
5. Shape: Reff is blind to this (0.250166 at the start, 0.250171 at the end in S_rhoST, 0.250174 in V1a; it follows the gas volume only). The smeared gas share, which is the shape-sensitive number, ends at 7.387% in S_rhoST (mean of the last 5 snapshots 6.087%) against 3.526% (3.592%) in V1a, from 3.300% and 3.656% at the first write. S_rhoST also deforms more: at t = 0.791 s the cells with alpha < 0.5 span x = 1.521 to 2.479 and reach r = 0.164 only (an elongated bubble, a capillary oscillation), against x = 1.771 to 2.229 and r = 0.238 in V1a at the same time (`b11_drift.log`).
6. Mass is not an issue in either run: gas volume +0.0097% (V1a) and +0.0061% (S_rhoST) over 3.95 s, from both the snapshots and the liquidVolume monitor.

Hypotheses (not tested).

- The density weighting multiplies the face force by about 2 rho_f / (rho_l + rho_g), which is about 0.002 on the gas side of the interface and about 2 on the liquid side. Its integral across the interface is unchanged (integral of 2 (alpha rho_l + (1-alpha) rho_g) / (rho_l + rho_g) over alpha from 0 to 1 is exactly 1), so the Laplace jump should be kept in equilibrium; what changes is where the force acts. Removing most of the force from the very light gas cells is consistent with the gas maximum dropping 2.4 to 3 times while the liquid does not improve. That is a reading of the numbers, not a demonstrated mechanism.
- A steady drift with no gravity and no imposed flow needs a net force on the fluid. If the weighted face force is no longer exactly balanced by the pressure gradient on the same faces, a net spurious force appears; the drift of 0.4 m/s with a rest period to 1.6 s, then growth, fits a slow build-up of such an imbalance together with the larger shape oscillation, but I have not tested it (a check would be to sum the discrete face force and the pressure-gradient force over the gas and see whether the net x-component is nonzero, or to repeat with the outlet closed).
- Because the drift is large and one-sided, the 3x lower max speed of S_rhoST probably should not be read as lower spurious currents. A fair statement from these two runs: the density-weighted force removes the gas-side speed peak, but in this single realisation it introduces a net bubble translation, larger liquid speeds, a larger smeared share and a worse pressure jump.

Caveats. One run of each, one mesh, so no estimate of run-to-run scatter; snapshot statistics are over 19 or 20 samples with a liquid-peak standard deviation of 0.22 m/s (S_rhoST) and 0.10 m/s (V1a), so the liquid-peak ratio of 1.57 is robust only in sign and rough size. The maxSpeed ratios rest on 2131 monitor samples. The V1a numbers in section 1 reproduce the paper row; no paper number was changed.

## Files written (all in `p3/`)

`b11_srho.py` (main computation), `b11_report.py` (builds this file and `b11_snap.csv` from `b11_res.pkl`), `b11_drift.py` and `b11_drift.log` (centroid, extent and sampling cell counts per snapshot), `b11_ux.py` and `b11_ux.log` (gas-weighted Ux), `b11_run.log`, `b11_res.pkl`, `b11_snap.csv`, `s_rho_result.md`. No existing file was edited.
